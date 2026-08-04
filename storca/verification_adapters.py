"""Default ORCA/Arkane adapters for the family-neutral verification engine."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Callable

from .arkane_runner import VerifiedArkaneRouteSpec, run_verified_arkane_route
from .decomposition_visuals import render_xyz_trajectory_gif, write_visual_manifest
from .methods import arkane_model_chemistry
from .reaction_path import run_generic_reaction_path
from .rmg_evidence import normalize_rmg_label
from .route_verify import route_spec_from_rmg_route


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
) -> Callable[[dict, Path], dict]:
    """Return the default route callback used by :mod:`verification_engine`."""
    if ncores < 1:
        raise ValueError("ORCA route verification needs at least one core")

    def verify(route: dict, workdir: Path) -> dict:
        spec = route_spec_from_rmg_route(parent_smiles, route)
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
        )
        visuals = _validated_trajectory_visual(result, workdir)
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
    bath_gas = _bath_gas(condition_contract)
    model_chemistry = arkane_model_chemistry(method_profile)

    def build(route: dict, path_result: dict, workdir: Path) -> dict:
        verification = path_result.get("route_verification") or {}
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
        )
        print(f"[STORCA] Arkane kinetics started for {route['route_id']}", flush=True)
        result = run_verified_arkane_route(arkane_spec, workdir, rmg_env=rmg_env)
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
