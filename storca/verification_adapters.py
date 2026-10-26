"""Default ORCA/Arkane adapters for the family-neutral verification engine."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Callable

from .arkane_runner import VerifiedArkaneRouteSpec, run_verified_arkane_route
from .decomposition_visuals import render_xyz_trajectory_gif, write_visual_manifest
from .methods import arkane_model_chemistry
from .method_applicability import evaluate_method_applicability
from .reaction_candidates import ReactionCandidate
from .reaction_path import run_generic_reaction_path
from .rmg_evidence import normalize_rmg_label
from .route_verify import route_spec_from_candidate_route, route_spec_from_rmg_route
from .verification_policy import GENERIC_ORCA_PATH_VERIFIER_POLICY_ID


_SUPPORTED_BATH_LABELS = {
    "nitrogen": "nitrogen",
    "n2": "nitrogen",
    "oxygen": "oxygen",
    "o2": "oxygen",
    "water": "water",
    "h2o": "water",
    "argon": "argon",
    "ar": "argon",
    "helium": "helium",
    "he": "helium",
}



def _normalized_species_labels(route: dict, side: str, count: int) -> tuple[str, ...]:
    records = (route.get("resolved_endpoints") or {}).get(side) or []
    if len(records) != count or any(not item.get("label") for item in records):
        raise ValueError(f"Resolved RMG {side} labels do not match the verified stationary points")
    labels = tuple(normalize_rmg_label(str(item["label"])) for item in records)
    # The same normalized label may repeat only for a stoichiometric copy of
    # the same structure.  Arkane defines that species once and retains the
    # coefficient through repeated reaction-list entries.
    structures: dict[str, str] = {}
    for label, record in zip(labels, records):
        smiles = str(record.get("smiles") or "")
        if label in structures and structures[label] != smiles:
            raise ValueError(f"RMG label normalization aliases different structures as {label!r}")
        structures[label] = smiles
    return labels


def _validated_stationary_labels(route: dict, side: str, stationary: list[dict]) -> tuple[str, ...]:
    """Bind ORCA stationary artifacts to the exact ordered RMG endpoints."""
    records = (route.get("resolved_endpoints") or {}).get(side) or []
    labels = _normalized_species_labels(route, side, len(stationary))
    try:
        from rdkit import Chem
    except ImportError as error:
        raise RuntimeError("RDKit is required to validate the ORCA-to-Arkane species handoff") from error

    def canonical(smiles: object) -> str:
        molecule = Chem.MolFromSmiles(str(smiles or ""))
        if molecule is None:
            raise ValueError("ORCA stationary-point handoff contains invalid SMILES")
        return Chem.MolToSmiles(molecule, canonical=True)

    for label, endpoint, artifact in zip(labels, records, stationary):
        artifact_label = artifact.get("original_label") or artifact.get("label")
        if artifact_label and normalize_rmg_label(str(artifact_label)) != label:
            raise ValueError(f"ORCA stationary-point order/label mismatch for {side} species {label!r}")
        if canonical(endpoint.get("smiles")) != canonical(artifact.get("smiles")):
            raise ValueError(f"ORCA stationary-point structure mismatch for {side} species {label!r}")
        if artifact.get("status") != "validated_minimum":
            raise ValueError(f"Arkane endpoint {label!r} is not a validated ORCA local minimum")
        if endpoint.get("multiplicity") is not None \
                and int(endpoint["multiplicity"]) != int(artifact.get("multiplicity", 0)):
            raise ValueError(f"ORCA stationary-point multiplicity mismatch for {side} species {label!r}")
    return labels


def _validate_reached_endpoint_artifact_contract(verification: dict, inputs: dict) -> dict | None:
    """Reopen the exact reached-minimum artifacts, not just matching SMILES.

    Legacy path records have no contract. New reached-endpoint paths must bind
    every optimized coordinate and frequency output before Arkane can read
    them; copying a different conformer's output under the same label fails.
    """
    if "reached_endpoint_artifact_contract" not in verification:
        return None
    contract = verification["reached_endpoint_artifact_contract"]
    if (not isinstance(contract, dict) or contract.get("schema_version") != 1
            or contract.get("kind") != "reached_endpoint_artifact_contract"):
        raise ValueError("Invalid reached-endpoint artifact contract")

    def validate_artifact(expected: object, actual_path: object, label: str) -> None:
        if not isinstance(expected, dict) or not expected.get("path") or not expected.get("sha256"):
            raise ValueError(f"Reached-endpoint contract lacks hashed {label}")
        path = Path(str(expected["path"]))
        digest = str(expected["sha256"])
        if (not path.is_absolute() or path.is_symlink() or not path.is_file()
                or len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest)):
            raise ValueError(f"Reached-endpoint contract has invalid {label} artifact")
        if not actual_path or Path(str(actual_path)).resolve() != path.resolve():
            raise ValueError(f"Arkane reached-endpoint {label} artifact differs from retained contract")
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Arkane reached-endpoint {label} artifact hash changed")

    for side in ("reactants", "products"):
        expected = contract.get(side)
        actual = inputs.get(side)
        if (not isinstance(expected, list) or not isinstance(actual, list)
                or not expected or len(expected) != len(actual)):
            raise ValueError(f"Reached-endpoint contract {side} cardinality differs")
        for index, (record, stationary) in enumerate(zip(expected, actual)):
            if not isinstance(record, dict) or not isinstance(stationary, dict):
                raise ValueError(f"Invalid reached-endpoint {side} record")
            for name in ("optimized_xyz", "frequency_output"):
                validate_artifact(record.get(name), stationary.get(name), f"{side}[{index}].{name}")
            if "optimization_output" in record:
                validate_artifact(record["optimization_output"], stationary.get("optimization_output"),
                                  f"{side}[{index}].optimization_output")
            if "frequency_hessian" in record:
                output = Path(str(stationary["frequency_output"]))
                validate_artifact(record["frequency_hessian"], output.with_suffix(".hess"),
                                  f"{side}[{index}].frequency_hessian")
    expected_ts = contract.get("transition_state")
    actual_ts = inputs.get("transition_state") or {}
    if not isinstance(expected_ts, dict) or not isinstance(actual_ts, dict):
        raise ValueError("Reached-endpoint contract lacks transition-state frequency binding")
    validate_artifact(expected_ts.get("frequency_output"), actual_ts.get("frequency_output"),
                      "transition_state.frequency_output")
    if "frequency_hessian" in expected_ts:
        output = Path(str(actual_ts["frequency_output"]))
        validate_artifact(expected_ts["frequency_hessian"], output.with_suffix(".hess"),
                          "transition_state.frequency_hessian")
    return contract


def _retain_reached_endpoint_artifact_admission(contract: dict, workdir: Path) -> dict:
    """Retain the admission without rewriting any old calculation or manifest."""
    path = Path(workdir) / "reached-endpoint-artifact-admission.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    record = {"schema_version": 1, "kind": "reached_endpoint_artifact_admission",
              "status": "passed", "contract": contract}
    if path.is_symlink():
        raise ValueError("Reached-endpoint artifact admission must be a retained regular file")
    serialized = json.dumps(record, sort_keys=True, indent=2, allow_nan=False) + "\n"
    try:
        with path.open("x") as stream:
            stream.write(serialized)
    except FileExistsError:
        if path.read_text() != serialized:
            raise ValueError("Retained reached-endpoint artifact admission differs")
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _bath_gas(condition_contract: dict) -> tuple[tuple[str, float], ...]:
    composition = condition_contract.get("composition") or {}
    target = str(condition_contract.get("target_label") or "stability")
    retained: dict[str, float] = {}
    unsupported = []
    for raw_label, raw_fraction in composition.items():
        fraction = float(raw_fraction)
        if not math.isfinite(fraction) or fraction < 0.0:
            raise ValueError(f"Invalid bath-gas mole fraction for {raw_label!r}")
        if str(raw_label) == target or fraction == 0.0:
            continue
        normalized = normalize_rmg_label(str(raw_label)).casefold()
        label = _SUPPORTED_BATH_LABELS.get(normalized)
        if label is None:
            unsupported.append(str(raw_label))
            continue
        retained[label] = retained.get(label, 0.0) + fraction
    if unsupported:
        raise ValueError(
            "Arkane has no explicit bath-gas structure for condition species: "
            + ", ".join(sorted(unsupported))
        )
    total = sum(retained.values())
    if total <= 0.0:
        raise ValueError("The immutable condition contract has no non-target Arkane bath gas")
    return tuple(sorted((label, fraction / total) for label, fraction in retained.items()))


def _grid_with_condition(values: tuple[float, ...], condition: float) -> tuple[float, ...]:
    grid = sorted({float(value) for value in values} | {float(condition)})
    return tuple(grid)


def _retain_method_applicability(decision: dict, workdir: Path) -> dict:
    """Create one immutable decision, or reuse only its exact scientific match."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    path = workdir / "method-applicability.json"
    if path.is_symlink():
        raise ValueError("Method applicability evidence must be a retained regular file")
    try:
        with path.open("x") as stream:
            stream.write(json.dumps(decision, indent=2, sort_keys=True, allow_nan=False) + "\n")
    except FileExistsError:
        try:
            retained = json.loads(path.read_text())
        except (OSError, ValueError) as error:
            raise ValueError("Retained method applicability evidence is unreadable or corrupt") from error
        if retained != decision:
            raise ValueError("Retained method applicability decision differs from the requested context")
    return {
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "decision_scientific_hash": decision["scientific_hash"],
        "policy_id": decision["policy_id"], "policy_hash": decision["policy_hash"],
    }


def _condition_collision_bound(route: dict, *, molecularity: int,
                               temperature_K: float, pressure_bar: float) -> tuple[float | None, str | None]:
    """Retain only a bimolecular collision ceiling evaluated at this contract."""
    if molecularity != 2:
        return None, None
    evidence = route.get("retained_rate_evidence") or {}
    units = str(evidence.get("collision_limit_units") or "")
    if units != "m^3/(mol*s)":
        raise ValueError("Bimolecular Arkane verification requires an SI molar collision-limit value")
    try:
        value = float(evidence["collision_limit_si"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("Bimolecular Arkane verification lacks a numeric collision limit") from error
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("Bimolecular Arkane collision limit must be positive and finite")
    source = str(evidence.get("collision_limit_source") or "").strip()
    if not source:
        raise ValueError("Bimolecular Arkane collision limit lacks provenance")
    evaluated = evidence.get("evaluated_at") or {}
    try:
        retained_temperature = float(evaluated["temperature_K"])
        retained_pressure = float(evaluated["pressure_bar"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("Bimolecular collision evidence lacks its evaluated condition") from error
    if not math.isclose(retained_temperature, temperature_K, rel_tol=1e-9, abs_tol=1e-9) \
            or not math.isclose(retained_pressure, pressure_bar, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError("Bimolecular collision evidence does not match the immutable condition")
    return value, source


def _validated_trajectory_visual(path_result: dict, workdir: Path) -> dict | None:
    verification = path_result.get("route_verification") or {}
    trajectory = verification.get("trajectory")
    path_classification = str(verification.get("path_classification") or "")
    allowed = {
        "verified_barriered_path",
        "barrierless_capture_candidate",
        "barrierless_dissociation_candidate",
        "barrierless_path_candidate",
    }
    if not trajectory or path_classification not in allowed:
        return None
    trajectory_path = Path(trajectory)
    if not trajectory_path.is_file():
        return None
    visual_dir = Path(workdir) / "visuals"
    label = (
        "ORCA-VERIFIED TS/IRC PATH"
        if path_classification == "verified_barriered_path"
        else "ORCA PATH CLASSIFICATION — RATE NOT YET VERIFIED"
    )
    animation = render_xyz_trajectory_gif(
        trajectory_path,
        visual_dir / "reaction-path.gif",
        title=f"{verification.get('route_id')}: {path_classification.replace('_', ' ')}",
        evidence_label=label,
    )
    manifest = write_visual_manifest(
        visual_dir / "visuals.json",
        route_id=str(verification.get("route_id")),
        evidence_level=path_classification,
        animation=animation,
        trajectory=trajectory_path,
    )
    return {
        "animation": str(animation),
        "trajectory": str(trajectory_path),
        "manifest": str(manifest),
        "evidence_label": label,
    }


def make_orca_route_verifier(
    parent_smiles: str,
    *,
    ncores: int,
    method_profile: dict,
    condition_contract: dict | None = None,
    timeout_seconds: float | None = 14400.0,
    orientations: int = 3,
    neb_images: int = 8,
    rmg_env: str | None = None,
) -> Callable[[dict, Path], dict]:
    """Return the default route callback used by :mod:`verification_engine`."""
    if ncores < 1:
        raise ValueError("ORCA route verification needs at least one core")

    def verify(route: dict, workdir: Path) -> dict:
        spec = (
            route_spec_from_candidate_route(parent_smiles, route)
            if route.get("native_candidate_id")
            else route_spec_from_rmg_route(parent_smiles, route)
        )
        recipe_mapping = None
        if route.get("rmg_mapping_source") and not spec.atom_mapping:
            from .rmg_mapping import apply_rmg_recipe_mapping
            spec, recipe_mapping = apply_rmg_recipe_mapping(
                spec, route["rmg_mapping_source"], Path(workdir) / "mapping-evidence",
                rmg_env=rmg_env, timeout_seconds=min(float(timeout_seconds or 120.0), 120.0),
            )
        print(
            f"[STORCA] ORCA verification started for {spec.route_id}: "
            f"{route.get('reaction_equation')}",
            flush=True,
        )
        result = run_generic_reaction_path(
            spec,
            workdir,
            ncores=ncores,
            method_keywords=method_profile.get("orca_keywords"),
            timeout_seconds=timeout_seconds,
            orientations=orientations,
            neb_images=neb_images,
            condition_contract=condition_contract,
            **({"retained_unexpected_proposal": route["unexpected_discovery_evidence"]["proposal"]}
               if route.get("unexpected_discovery_evidence") else {}),
        )
        visuals = _validated_trajectory_visual(result, workdir)
        if recipe_mapping is not None:
            result["recipe_atom_mapping"] = recipe_mapping
        if visuals:
            result["validated_visuals"] = visuals
        verification = result.get("route_verification") or {}
        print(
            f"[STORCA] ORCA verification finished for {spec.route_id}: "
            f"{verification.get('status')} / {verification.get('path_classification')}",
            flush=True,
        )
        return result

    return verify


def make_arkane_rate_builder(
    condition_contract: dict,
    *,
    method_profile: dict,
    rmg_env: str | None,
) -> Callable[[dict, dict, Path], dict]:
    """Return a condition-bound Arkane callback for a verified ORCA path."""
    temperature = float(condition_contract["temperature_K"])
    pressure = float(condition_contract["pressure_bar"])
    model_chemistry = arkane_model_chemistry(method_profile)

    def build(route: dict, path_result: dict, workdir: Path) -> dict:
        verification = path_result.get("route_verification") or {}
        if verification.get("status") == "verified_multistep_path":
            multistep = verification.get("multistep_path") or {}
            segments = list(multistep.get("segments") or [])
            if not segments or any(
                item.get("status") != "verified_transition_state_path"
                or (item.get("verification") or {}).get("status") != "verified_transition_state_path"
                or (item.get("verification") or {}).get("path_classification") != "verified_barriered_path"
                or item.get("rate_input_status") != "ready"
                or not item.get("arkane_inputs")
                for item in segments
            ):
                raise ValueError(
                    "Multistep Arkane generation requires every elementary segment to have "
                    "verified TS/frequency/IRC evidence and retained endpoint frequencies"
                )
            built_segments: list[dict] = []
            for index, segment in enumerate(segments, start=1):
                candidate = ReactionCandidate.from_dict(segment["native_candidate"])
                segment_route = candidate.as_route_dict(
                    parent_smiles=candidate.reactants[0],
                    target_label=(candidate.reactant_labels[0] if candidate.reactant_labels else "segment"),
                )
                segment_verification = dict(segment.get("verification") or {})
                segment_verification.update(
                    arkane_inputs=segment["arkane_inputs"],
                    native_candidate_id=candidate.candidate_id,
                    parent_multistep_route_id=route.get("route_id"),
                )
                segment_verification.setdefault("selected_multiplicity", candidate.multiplicity)
                segment_verification.setdefault("mapping", segment.get("segment_mapping") or {})
                segment_path_result = {
                    **{key: path_result[key] for key in (
                        "environment_applicability", "surface_resolution_scientific_context",
                    ) if key in path_result},
                    "route": segment.get("segment_route") or segment_route,
                    "route_verification": segment_verification,
                }
                segment_rate = build(
                    segment_route,
                    segment_path_result,
                    Path(workdir) / f"segment-{index:02d}-{candidate.candidate_id}",
                )
                built_segments.append({
                    "index": index,
                    "label": segment.get("label"),
                    "route": segment_route,
                    "path_result": segment_path_result,
                    "rate_result": segment_rate,
                    "native_candidate_id": candidate.candidate_id,
                })
            return {
                "schema_version": 1,
                "kind": "storca_verified_multistep_rate_bundle",
                "status": "completed_multistep",
                "parent_route_id": route.get("route_id"),
                "condition_contract": condition_contract,
                "segments": built_segments,
                "effective_single_barrier_used": False,
                "network_insertion_required": True,
                "rate_claim_allowed": True,
                "lifetime_claim_allowed": False,
                "interpretation": (
                    "Every elementary segment has its own condition-bound Arkane result. "
                    "No macroscopic observable is allowed until all segment reactions and the "
                    "intermediate are inserted and the repaired network is repropagated."
                ),
            }
        if verification.get("status") != "verified_transition_state_path":
            raise ValueError(
                "Arkane TST requires a verified TS/frequency/IRC path; "
                f"received {verification.get('status')!r}"
            )
        inputs = verification.get("arkane_inputs") or {}
        reactants = list(inputs.get("reactants") or [])
        products = list(inputs.get("products") or [])
        transition_state = inputs.get("transition_state") or {}
        if not reactants or not products or not transition_state.get("frequency_output"):
            raise ValueError("Verified path lacks separated stationary-point artifacts required by Arkane")
        reached_contract = _validate_reached_endpoint_artifact_contract(verification, inputs)
        if "conformer_rate_basis" in verification or reached_contract is not None:
            basis = verification.get("conformer_rate_basis")
            if (not isinstance(basis, dict)
                    or basis.get("kind") != "reached_endpoint_conformer_rate_basis"
                    or not isinstance(basis.get("network_reinjection_allowed"), bool)):
                raise ValueError("Invalid or missing reached-endpoint conformer rate basis")
        reached_evidence = (
            _retain_reached_endpoint_artifact_admission(reached_contract, workdir)
            if reached_contract is not None else None
        )
        reactant_labels = _validated_stationary_labels(route, "reactants", reactants)
        product_labels = _validated_stationary_labels(route, "products", products)
        collision_limit, collision_source = _condition_collision_bound(
            route, molecularity=len(reactants), temperature_K=temperature, pressure_bar=pressure,
        )
        selected_multiplicity = verification.get("selected_multiplicity")
        if selected_multiplicity is None:
            selected_multiplicity = route.get("multiplicity")
        if selected_multiplicity is None:
            raise ValueError("Verified ORCA path does not retain its selected total-spin multiplicity")
        temperatures = _grid_with_condition(
            (250.0, 273.15, 298.15, 323.15, 400.0, 500.0, 750.0, 1000.0), temperature,
        )
        pressures = _grid_with_condition((1.0e-6, 1.0e-4, 1.0e-2, 1.0, 10.0), pressure)
        # Follow VerifiedArkaneRouteSpec's existing well/channel branch. HPL
        # TST uses no bath or collision-energy-transfer model; requiring a
        # non-target gas here would reject a pure-gas isomerization before
        # any rate was built. This does not establish HPL applicability.
        use_pdep = {len(reactants), len(products)} == {1, 2}
        bath_gas = _bath_gas(condition_contract) if use_pdep else ()
        arkane_spec = VerifiedArkaneRouteSpec(
            label=str(route["route_id"]),
            reactant_labels=reactant_labels,
            reactant_smiles=tuple(str(item["smiles"]) for item in reactants),
            reactant_orca_outputs=tuple(Path(item["frequency_output"]) for item in reactants),
            reactant_multiplicities=tuple(int(item["multiplicity"]) for item in reactants),
            product_labels=product_labels,
            product_smiles=tuple(str(item["smiles"]) for item in products),
            product_orca_outputs=tuple(Path(item["frequency_output"]) for item in products),
            product_multiplicities=tuple(int(item["multiplicity"]) for item in products),
            transition_state_orca_output=Path(transition_state["frequency_output"]),
            transition_state_multiplicity=int(selected_multiplicity),
            model_chemistry=model_chemistry,
            reaction_equation=str(
                route.get("reaction_equation") or route.get("printed_reaction_equation") or ""
            ) or None,
            temperatures_K=temperatures,
            pressures_bar=pressures,
            bath_gas=bath_gas,
            frequency_scale_factor=float(method_profile.get("harmonic_scale_factor", 1.0)),
            condition_temperature_K=temperature,
            condition_pressure_bar=pressure,
            collision_limit_m3_mol_s=(float(collision_limit) if collision_limit is not None else None),
            collision_limit_source=(str(collision_source) if collision_source else None),
            explicit_pressure_dependence_model=(
                route.get("pressure_dependence_model")
                or verification.get("pressure_dependence_model")
            ),
            conformer_rate_basis=(
                verification["conformer_rate_basis"]
                if "conformer_rate_basis" in verification else None
            ),
        )
        applicability = evaluate_method_applicability(
            route, condition_contract=condition_contract, method_profile=method_profile,
            path_result=path_result, rate_theory="harmonic_gas_tst",
            tunneling_model=arkane_spec.tunneling,
            pressure_dependent=arkane_spec.use_pressure_dependence,
        )
        applicability_evidence = _retain_method_applicability(applicability, workdir)
        if not applicability["central_rate_allowed"]:
            raise ValueError(
                "Ordinary Arkane TST method applicability failed: "
                + ", ".join(applicability["blocking_reasons"])
            )
        print(f"[STORCA] Arkane kinetics started for {route['route_id']}", flush=True)
        result = run_verified_arkane_route(arkane_spec, workdir, rmg_env=rmg_env)
        if reached_contract is not None:
            _validate_reached_endpoint_artifact_contract(verification, inputs)
        # Derived callback evidence accompanies the returned repair record.
        # Existing checksummed generated-library artifacts are not rewritten.
        result = {
            **result, "method_applicability": applicability,
            "method_applicability_evidence": applicability_evidence,
            **({"reached_endpoint_artifact_evidence": reached_evidence}
               if reached_evidence is not None else {}),
            **({"conformer_rate_basis": dict(arkane_spec.conformer_rate_basis)}
               if arkane_spec.conformer_rate_basis is not None else {}),
        }
        if isinstance(result.get("generated_kinetics_library"), dict):
            result["generated_kinetics_library"] = {
                **result["generated_kinetics_library"],
                "method_applicability": applicability,
                "method_applicability_evidence": applicability_evidence,
                **({"reached_endpoint_artifact_evidence": reached_evidence}
                   if reached_evidence is not None else {}),
                **({"conformer_rate_basis": dict(arkane_spec.conformer_rate_basis)}
                   if arkane_spec.conformer_rate_basis is not None else {}),
            }
        print(f"[STORCA] Arkane kinetics finished for {route['route_id']}: {result.get('status')}", flush=True)
        if result.get("status") != "completed":
            nested = result.get("generated_kinetics_library") or {}
            raise RuntimeError(
                result.get("failure_reason") or nested.get("failure_reason")
                or "Arkane verification did not complete"
            )
        return result

    return build


__all__ = ["make_arkane_rate_builder", "make_orca_route_verifier"]
