"""Bounded graph-edit discovery with deliberately non-kinetic xTB screening.

This backend proposes connectivity hypotheses independently of RMG.  It does
not estimate a transition-state barrier, a rate coefficient, or a lifetime.
Only the existing ORCA path verifier plus Arkane/repropagation may promote a
surviving candidate into kinetic evidence.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import math
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Callable, Mapping

from .reaction_candidates import ReactionCandidate, canonical_smiles, endpoint_state
from .reaction_discovery import DiscoveryResult


HARTREE_TO_KCAL_MOL = 627.509474
GRAPH_XTB_DISCOVERY_SCHEMA_VERSION = 3

# Conservative names supported by both the documented GFN2-ALPB model and
# the installed 6.7 series. No guessed solvent parameters or gas fallback.
ALPB_SOLVENTS = frozenset({
    "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene", "ch2cl2",
    "chcl3", "cs2", "dioxane", "dmf", "dmso", "ether", "ethylacetate", "furane",
    "hexane", "methanol", "nitromethane", "octanol", "woctanol", "phenol",
    "toluene", "thf", "water",
})


def _screening_environment(value: Mapping[str, Any] | None) -> dict[str, Any]:
    """Physical screen model, separate from the requested ORCA rate model."""
    value = dict(value or {})
    phase = str(value.get("phase") or "").casefold()
    aliases = [str(value.get(key) or "").strip().casefold() for key in ("phase_model", "model", "environment_kind")]
    for alias in aliases:
        if not alias:
            continue
        if any(term in alias for term in ("surface", "solid", "interface", "electrochem", "explicit", "bulk", "neat")):
            raise ValueError("Outside implicit gas/solution graph/xTB screening environment")
        declared = "gas" if "gas" in alias else "solution" if "solution" in alias or "solvent" in alias else None
        if declared is None:
            raise ValueError("Unresolved graph/xTB screening environment model")
        if phase and phase != declared:
            raise ValueError("Conflicting graph/xTB screening environment aliases")
        phase = declared
    solvent = str(value.get("solvent") or "").strip().casefold()
    model = str(value.get("solvent_model") or "").upper()
    if solvent or model:
        if phase not in {"", "solution"} or "gas" in str(value.get("phase_model") or "").casefold():
            raise ValueError("Conflicting graph/xTB screening environment phase and solvent")
        if not solvent or model not in {"SMD", "CPCM", "ALPB"}:
            raise ValueError("Solution graph/xTB discovery needs a named solvent and declared implicit model")
        if solvent not in ALPB_SOLVENTS:
            raise ValueError(f"No supported GFN2-ALPB screening parameters for solvent: {solvent}")
        return {"phase": "solution", "solvent": solvent, "solvent_model": "ALPB", "reference_state": "gsolv"}
    if phase not in {"", "gas"}:
        raise ValueError("Graph/xTB discovery cannot silently use gas screening for this environment")
    return {"phase": "gas"}


@dataclass(frozen=True)
class GraphXtbBudget:
    """Explicitly finite edit and endpoint-screening budget."""

    max_generated_graph_candidates: int = 128
    max_structurally_valid_candidates: int = 48
    max_products_per_reactant_state: int = 48
    max_xtb_optimizations: int = 48
    max_xtb_failures: int = 8
    max_candidates_for_verification: int = 4
    max_reactive_partners: int = 4
    max_bimolecular_generated_candidates: int = 128
    max_bimolecular_structurally_valid_candidates: int = 48
    energy_window_kcal_mol: float = 80.0
    timeout_seconds: float = 300.0
    overall_timeout_seconds: float | None = None
    allow_dissociation: bool = True
    allow_bond_formation: bool = True
    allow_bond_rearrangement: bool = True
    allow_hydrogen_transfer: bool = True
    allow_electrocyclic_bond_rearrangement: bool = True
    allow_bimolecular_radical_association: bool = True
    allow_bimolecular_hydrogen_abstraction: bool = True

    def __post_init__(self) -> None:
        for name in (
            "max_generated_graph_candidates", "max_structurally_valid_candidates",
            "max_products_per_reactant_state", "max_xtb_optimizations",
            "max_xtb_failures", "max_candidates_for_verification",
            "max_reactive_partners", "max_bimolecular_generated_candidates",
            "max_bimolecular_structurally_valid_candidates",
        ):
            if int(getattr(self, name)) < 0:
                raise ValueError(f"{name} must be non-negative")
        if not math.isfinite(self.energy_window_kcal_mol):
            raise ValueError("energy_window_kcal_mol must be finite")
        if not math.isfinite(self.timeout_seconds) or self.timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be positive and finite")
        if self.overall_timeout_seconds is not None and (
            isinstance(self.overall_timeout_seconds, bool)
            or not math.isfinite(self.overall_timeout_seconds) or self.overall_timeout_seconds <= 0
        ):
            raise ValueError("overall_timeout_seconds must be positive and finite when declared")

    @classmethod
    def coerce(cls, value: "GraphXtbBudget | Mapping[str, Any] | None") -> "GraphXtbBudget":
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            allowed = {name: value[name] for name in cls.__dataclass_fields__ if name in value}
            return cls(**allowed)
        raise TypeError("Graph xTB discovery budget must be GraphXtbBudget, mapping, or None")

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": GRAPH_XTB_DISCOVERY_SCHEMA_VERSION,
            "kind": "graph_xtb_discovery_budget",
            **self.__dict__,
        }


XtbEndpointRunner = Callable[..., dict[str, Any]]


def _sha256_file(path: Path) -> str | None:
    path = Path(path)
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, value: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    return path


def _extract_xtb_energy(output: str, json_path: Path | None = None) -> float:
    """Parse a finite xTB total electronic energy in Hartree."""
    matches = re.findall(
        r"TOTAL\s+ENERGY\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
        output,
        flags=re.IGNORECASE,
    )
    if not matches:
        raise ValueError("xTB output contains no finite TOTAL ENERGY")
    energy = float(matches[-1])
    if not math.isfinite(energy):
        raise ValueError("xTB output contains a non-finite total energy")
    if json_path and json_path.is_file():
        payload = json.loads(json_path.read_text())
        value = payload.get("total energy") if isinstance(payload, dict) else None
        if (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
                or not math.isclose(value, energy, rel_tol=0.0, abs_tol=5.1e-9)):
            # Native 6.7 JSON prints eight decimals, stdout twelve. This is
            # source-rounding agreement, not a scientific uncertainty bound.
            raise ValueError("xTB JSON/stdout total energy is absent, invalid or inconsistent")
    return energy


def _observe_xtb_screen(output: str, expected: Mapping[str, Any]) -> dict[str, Any]:
    if not re.search(r"\bHamiltonian\s*[:=]?\s*GFN2-xTB\b", output, re.IGNORECASE):
        raise ValueError("Native xTB output does not establish GFN2 Hamiltonian")
    if expected["phase"] == "gas":
        if re.search(r"\bSolvation model\s*[:=]?\s*(ALPB|GBSA|COSMO)", output, re.IGNORECASE):
            raise ValueError("Native gas screen unexpectedly used a solvation model")
    else:
        if not re.search(r"\bSolvation model\s*[:=]?\s*ALPB\b", output, re.IGNORECASE):
            raise ValueError("Native solution screen does not establish ALPB model")
        solvent = re.search(r"\bSolvent\s*[:=]?\s+(\S+)", output, re.IGNORECASE)
        if not solvent or solvent.group(1).casefold() != expected["solvent"]:
            raise ValueError("Native ALPB solvent differs from the screen contract")
        if not re.search(r"\bReference state\s*[:=]?\s*gsolv\s*\[\s*1\s*M\s+gas/solution\s*\]", output, re.IGNORECASE):
            raise ValueError("Native ALPB gsolv reference state was not observed")
    return dict(expected)


def _validate_xtb_json_invocation(path: Path, command: list[str], version: str) -> None:
    if path.is_file():
        value = json.loads(path.read_text())
        if (value.get("method") != "GFN2-xTB" or str(value.get("xtb version") or "").split()[0:1] != [version]
                or value.get("program call") != " ".join(command)):
            raise ValueError("Native xTB JSON method/version/invocation differs from its contract")


def _validate_xtb_endpoint_state(output: str, json_path: Path, charge: int, unpaired: int) -> None:
    """Check observed native state, not merely intended --chrg/--uhf flags.

    This binds the charge/unpaired-electron calculation setting, not spin
    purity, ground-state ordering or a single-reference applicability claim.
    Repeated native summaries must all agree; missing diagnostics refuse.
    """
    for label, expected in ((r"net\s+charge", charge), (r"unpaired\s+electrons", unpaired)):
        values = re.findall(r"\b" + label + r"\s*[:=]?\s*([+-]?\d+)(?=\s|$)", output, flags=re.IGNORECASE)
        if not values or any(int(value) != expected for value in values):
            raise ValueError("Native xTB observed charge/unpaired-electron setting differs from its contract")
    if json_path.is_file():
        value = json.loads(json_path.read_text()).get("number of unpaired electrons")
        if type(value) is not int or value != unpaired:
            raise ValueError("Native xTB JSON unpaired-electron setting differs from its contract")


def _xtb_version(executable: str, timeout_seconds: float = 30.0) -> str | None:
    try:
        result = subprocess.run(
            [executable, "--version"], capture_output=True, text=True, timeout=timeout_seconds,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    text = f"{result.stdout}\n{result.stderr}"
    match = re.search(r"xtb\s+version\s+([^\s]+)", text, flags=re.IGNORECASE)
    return match.group(1) if result.returncode == 0 and match else None


def _endpoint_geometry(smiles: str, path: Path) -> Path:
    """Generate a retained RDKit seed geometry only for xTB endpoint screening."""
    from src.molecule_tools import smiles_to_xyz

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file():
        return path
    return smiles_to_xyz(smiles, path, n_confs=10)


def run_xtb_endpoint_screen(
    smiles: str,
    workdir: Path,
    *,
    charge: int,
    multiplicity: int,
    timeout_seconds: float,
    executable: str | None = None,
    ncores: int = 1,
    environment: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Run or reuse one GFN2-xTB endpoint optimization for screening only."""
    if multiplicity < 1:
        raise ValueError("xTB endpoint multiplicity must be positive")
    workdir = Path(workdir)
    seed: Path | None = None
    try:
        screening_environment = _screening_environment(environment)
        deadline = time.monotonic() + float(timeout_seconds)
        if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
            raise ValueError("xTB screening timeout must be finite and positive")
        workdir.mkdir(parents=True, exist_ok=True)
        seed = _endpoint_geometry(smiles, workdir / "initial.xyz")
        from src.orca_runner import find_xtb

        xtb = str(Path(executable or find_xtb()).resolve())
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise RuntimeError("xTB endpoint budget exhausted during seed acquisition")
        version = _xtb_version(xtb, timeout_seconds=min(30.0, remaining))
        if version is None:
            raise RuntimeError("Could not determine xTB executable version")
        contract = {
            "schema_version": GRAPH_XTB_DISCOVERY_SCHEMA_VERSION,
            "method": "GFN2-xTB",
            "optimization": "loose",
            "smiles": smiles,
            "charge": int(charge),
            "unpaired_electrons": int(multiplicity - 1),
            "timeout_seconds": float(timeout_seconds),
            "ncores": int(ncores),
            "input_xyz_sha256": _sha256_file(seed),
            "xtb_executable": xtb,
            "xtb_version": version,
            "screening_environment": screening_environment,
        }
        contract_path = workdir / "screen-contract.json"
        record_path = workdir / "screen-result.json"
        if contract_path.is_file() and record_path.is_file():
            retained_contract = json.loads(contract_path.read_text())
            retained_record = json.loads(record_path.read_text())
            if retained_contract == contract:
                if retained_record.get("status") == "completed":
                    artifacts = retained_record.get("artifacts") or {}
                    for key, expected in (("input_xyz", seed), ("output", workdir / "xtb.out"),
                                          ("optimized_xyz", workdir / "optimized.xyz")):
                        if (artifacts.get(key) != str(expected) or not artifacts.get(key + "_sha256")
                                or _sha256_file(expected) != artifacts[key + "_sha256"]):
                            raise RuntimeError("Retained xTB screening artifact hash/identity changed")
                    json_output = workdir / "xtbout.json"
                    if json_output.is_file() != bool(artifacts.get("json_output")) or (
                        json_output.is_file() and _sha256_file(json_output) != artifacts.get("json_output_sha256")
                    ):
                        raise RuntimeError("Retained xTB JSON energy source changed")
                    observed = _extract_xtb_energy((workdir / "xtb.out").read_text(), json_output)
                    if observed != retained_record.get("energy_hartree"):
                        raise RuntimeError("Retained xTB screening energy differs from its native source")
                    if retained_record.get("screening_environment") != screening_environment:
                        raise RuntimeError("Retained xTB screening environment changed")
                    _observe_xtb_screen((workdir / "xtb.out").read_text(), screening_environment)
                    _validate_xtb_json_invocation(json_output, retained_record["command"], version)
                    _validate_xtb_endpoint_state((workdir / "xtb.out").read_text(), json_output, charge, multiplicity - 1)
                    from .route_geometry import validate_declared_connectivity
                    if not validate_declared_connectivity(workdir / "optimized.xyz", (smiles,))["valid"]:
                        raise RuntimeError("Retained xTB endpoint no longer represents its declared molecular identity")
                return retained_record
            raise RuntimeError("Refusing to reuse an xTB endpoint directory with a different calculation contract")
        if contract_path.exists() or record_path.exists():
            raise RuntimeError("Incomplete retained xTB endpoint contract; use a new owned directory")
        _write_json(contract_path, contract)
        output = workdir / "xtb.out"
        command = [
            xtb, seed.name, "--gfn", "2", "--opt", "loose", "--chrg", str(charge),
            "--uhf", str(multiplicity - 1), "--json", "--norestart",
        ]
        if screening_environment["phase"] == "solution":
            command.extend(["--alpb", screening_environment["solvent"], "gsolv"])
        thread_environment = {
            "OMP_NUM_THREADS": str(max(1, int(ncores))),
            "MKL_NUM_THREADS": str(max(1, int(ncores))),
        }
        try:
            result = subprocess.run(
                command, cwd=workdir, capture_output=True, text=True,
                timeout=max(0.001, deadline - time.monotonic()),
                env={**__import__("os").environ, **thread_environment},
            )
            output.write_text((result.stdout or "") + ("\n" if result.stdout and result.stderr else "") + (result.stderr or ""))
            if result.returncode != 0:
                raise RuntimeError(f"xTB exited with return code {result.returncode}")
            text = output.read_text(errors="replace")
            if "GEOMETRY OPTIMIZATION CONVERGED" not in text.upper():
                raise RuntimeError("xTB endpoint optimizer convergence was not observed")
            energy = _extract_xtb_energy(text, workdir / "xtbout.json")
            observed_environment = _observe_xtb_screen(text, screening_environment)
            _validate_xtb_json_invocation(workdir / "xtbout.json", command, version)
            _validate_xtb_endpoint_state(text, workdir / "xtbout.json", charge, multiplicity - 1)
            optimized = workdir / "xtbopt.xyz"
            retained_optimized = workdir / "optimized.xyz"
            if not optimized.is_file():
                raise RuntimeError("xTB endpoint lacks its optimized geometry")
            shutil.copy2(optimized, retained_optimized)
            from .route_geometry import validate_declared_connectivity
            geometry_validation = validate_declared_connectivity(retained_optimized, (smiles,))
            if geometry_validation.get("valid") is not True:
                raise RuntimeError("xTB optimized endpoint changed or failed declared molecular connectivity")
            record = {
                "status": "completed",
                "energy_hartree": energy,
                "method": "GFN2-xTB",
                "command": command,
                "xtb_version": version,
                "screening_environment": observed_environment,
                "geometry_validation": geometry_validation,
                "artifacts": {
                    "input_xyz": str(seed),
                    "output": str(output),
                    "optimized_xyz": str(retained_optimized) if retained_optimized.is_file() else None,
                    "input_xyz_sha256": _sha256_file(seed),
                    "output_sha256": _sha256_file(output),
                    "optimized_xyz_sha256": _sha256_file(retained_optimized),
                    "json_output": str(workdir / "xtbout.json") if (workdir / "xtbout.json").is_file() else None,
                    "json_output_sha256": _sha256_file(workdir / "xtbout.json"),
                },
                "warnings": [
                    "GFN2-xTB endpoint energy is screening evidence only; it is not a barrier or kinetic rate.",
                    "ALPB total energies include solvation terms but no thermal barrier; ALPB is not SMD/CPCM rate evidence.",
                ],
            }
        except Exception as error:
            record = {
                "status": "xtb_failed",
                "failure_reason": f"{type(error).__name__}: {error}",
                "method": "GFN2-xTB",
                "command": command,
                "xtb_version": version,
                "screening_environment": screening_environment,
                "artifacts": {
                    "input_xyz": str(seed),
                    "output": str(output),
                    "input_xyz_sha256": _sha256_file(seed),
                    "output_sha256": _sha256_file(output),
                },
            }
        _write_json(record_path, record)
        return record
    except Exception as error:
        return {
            "status": "xtb_failed",
            "failure_reason": f"{type(error).__name__}: {error}",
            "method": "GFN2-xTB",
            "command": None,
            "xtb_version": None,
            "artifacts": {
                "input_xyz": str(seed) if seed is not None else None,
                "input_xyz_sha256": _sha256_file(seed) if seed is not None else None,
            },
        }


def validate_graph_xtb_candidate_sources(candidate: ReactionCandidate) -> dict[str, Any]:
    """Reopen a selected graph screen without jobs, version reads or writes.

    A candidate's copied status/environment is not evidence that its native
    endpoint sources still exist or agree. This admission check is deliberately
    separate from discovery generation: synthetic/test backends may generate
    hypotheses, but cannot acquire native screening authority from metadata.
    It establishes screening provenance only, never a saddle, rate or coverage.
    """
    if "graph_xtb" not in candidate.discovery_sources:
        raise ValueError("Source reopening requires a graph/xTB discovery assertion")
    if (candidate.screening_status != "candidate_unverified"
            or candidate.verification_status != "candidate_unverified"):
        raise ValueError("Only a surviving unverified graph/xTB candidate is eligible")
    expected_environment = _screening_environment(candidate.environment)
    screen = candidate.provenance.get("xtb_screen") or {}
    if (screen.get("status") != "completed" or screen.get("method") != "GFN2-xTB"
            or screen.get("screening_environment") != expected_environment):
        raise ValueError("Candidate lacks completed environment-matched native xTB screen evidence")
    endpoints = [(side, index, smiles) for side, values in (
        ("reactant", candidate.reactants), ("product", candidate.products))
        for index, smiles in enumerate(values)]
    records = screen.get("endpoint_records")
    if not isinstance(records, list) or len(records) != len(endpoints):
        raise ValueError("Native xTB screen needs exactly one record per stoichiometric endpoint")
    source_hashes: dict[str, str] = {}
    totals = {"reactant": 0.0, "product": 0.0}
    try:
        for (side, index, smiles), record in zip(endpoints, records):
            if (record.get("side") != side or type(record.get("endpoint_index")) is not int
                    or record["endpoint_index"] != index or record.get("smiles") != smiles
                    or record.get("status") != "completed" or record.get("method") != "GFN2-xTB"
                    or record.get("screening_environment") != expected_environment):
                raise ValueError("Native xTB endpoint identity/status/environment differs from candidate")
            artifacts = record["artifacts"]
            output = Path(artifacts["output"])
            root = output.parent
            contract_path, result_path = root / "screen-contract.json", root / "screen-result.json"
            contract = json.loads(contract_path.read_text())
            retained = json.loads(result_path.read_text())
            native_record = {key: value for key, value in record.items()
                             if key not in {"side", "endpoint_index", "smiles"}}
            if retained != native_record:
                raise ValueError("Copied xTB endpoint record differs from retained native result")
            states = endpoint_state(smiles)
            supplied_charge = candidate.reactant_charges if side == "reactant" else candidate.product_charges
            supplied_spin = candidate.reactant_multiplicities if side == "reactant" else candidate.product_multiplicities
            side_count = len(candidate.reactants if side == "reactant" else candidate.products)
            charge = supplied_charge[index] if len(supplied_charge) == side_count else states["charge"]
            spin = supplied_spin[index] if len(supplied_spin) == side_count else states["multiplicity"]
            if (charge is None or spin is None or states["charge"] != charge
                    or (states["multiplicity"] is not None and states["multiplicity"] != spin)
                    or contract.get("schema_version") != GRAPH_XTB_DISCOVERY_SCHEMA_VERSION
                    or contract.get("method") != "GFN2-xTB" or contract.get("optimization") != "loose"
                    or contract.get("smiles") != smiles or type(contract.get("charge")) is not int
                    or contract.get("charge") != charge or type(contract.get("unpaired_electrons")) is not int
                    or contract.get("unpaired_electrons") != spin - 1
                    or contract.get("screening_environment") != expected_environment
                    or contract.get("xtb_version") != record.get("xtb_version")):
                raise ValueError("Native xTB calculation contract differs from candidate state/environment")
            command = [contract["xtb_executable"], "initial.xyz", "--gfn", "2", "--opt", "loose",
                       "--chrg", str(charge), "--uhf", str(spin - 1), "--json", "--norestart"]
            if expected_environment["phase"] == "solution":
                command.extend(["--alpb", expected_environment["solvent"], "gsolv"])
            if record.get("command") != command:
                raise ValueError("Native xTB invocation differs from its physical contract")
            for key, expected in (("input_xyz", root / "initial.xyz"), ("output", root / "xtb.out"),
                                  ("optimized_xyz", root / "optimized.xyz")):
                digest = artifacts.get(key + "_sha256")
                if artifacts.get(key) != str(expected) or not digest or _sha256_file(expected) != digest:
                    raise ValueError("Native xTB screening source hash/identity changed: " + key)
                source_hashes[str(expected)] = digest
            if contract.get("input_xyz_sha256") != artifacts["input_xyz_sha256"]:
                raise ValueError("Native xTB seed differs from the calculation contract")
            json_output = root / "xtbout.json"
            if bool(artifacts.get("json_output")) != json_output.is_file():
                raise ValueError("Native xTB JSON source presence changed")
            if json_output.is_file():
                if (artifacts["json_output"] != str(json_output) or not artifacts.get("json_output_sha256")
                        or _sha256_file(json_output) != artifacts["json_output_sha256"]):
                    raise ValueError("Native xTB JSON source hash/identity changed")
                source_hashes[str(json_output)] = artifacts["json_output_sha256"]
            text = output.read_text()
            if "GEOMETRY OPTIMIZATION CONVERGED" not in text.upper():
                raise ValueError("Native xTB endpoint optimizer convergence is absent")
            _observe_xtb_screen(text, expected_environment)
            _validate_xtb_json_invocation(json_output, command, record["xtb_version"])
            _validate_xtb_endpoint_state(text, json_output, charge, spin - 1)
            energy = _extract_xtb_energy(text, json_output)
            if (isinstance(record.get("energy_hartree"), bool)
                    or record.get("energy_hartree") != energy):
                raise ValueError("Copied xTB endpoint energy differs from its native source")
            from .route_geometry import validate_declared_connectivity
            if not validate_declared_connectivity(root / "optimized.xyz", (smiles,))["valid"]:
                raise ValueError("Native xTB endpoint no longer represents candidate connectivity")
            for path in (contract_path, result_path):
                source_hashes[str(path)] = _sha256_file(path)
            totals[side] += energy
        for side, total in totals.items():
            field = (side + "_screen_total_energy_hartree" if expected_environment["phase"] == "solution"
                     else side + "_electronic_energy_hartree")
            if isinstance(screen.get(field), bool) or screen.get(field) != total:
                raise ValueError("Copied xTB side-energy sum differs from native endpoint evidence")
        delta = (totals["product"] - totals["reactant"]) * HARTREE_TO_KCAL_MOL
        if (isinstance(screen.get("estimated_reaction_energy_kcal_mol"), bool)
                or screen.get("estimated_reaction_energy_kcal_mol") != delta
                or candidate.estimated_reaction_energy_kcal_mol != delta):
            raise ValueError("Copied xTB reaction-screen energy differs from native evidence")
    except (OSError, KeyError, TypeError, AttributeError, json.JSONDecodeError) as error:
        raise ValueError("Missing or invalid native xTB candidate source: " + str(error)) from error
    return {"kind": "storca_graph_xtb_source_reopening", "status": "screen_sources_reopened",
            "screening_environment": expected_environment, "endpoint_count": len(endpoints),
            "source_hashes": source_hashes, "native_jobs_executed": 0,
            "interpretation": "Source-consistent endpoint screening only; no saddle, rate or network-coverage authority."}


def _inventory(molecules: list[Any]) -> tuple[dict[tuple[int, int], int], int]:
    """Return an explicit-H isotope inventory and total formal charge."""
    from rdkit import Chem

    atoms: dict[tuple[int, int], int] = {}
    charge = 0
    for molecule in molecules:
        explicit = Chem.AddHs(molecule)
        for atom in explicit.GetAtoms():
            key = (atom.GetAtomicNum(), int(atom.GetIsotope()))
            atoms[key] = atoms.get(key, 0) + 1
        charge += int(Chem.GetFormalCharge(molecule))
    return atoms, charge


def _mapped_fragments(molecule: Any) -> tuple[str, ...]:
    """Return atom-mapped fragments in canonical graph form."""
    from rdkit import Chem

    fragments = Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True)
    return tuple(sorted(
        Chem.MolToSmiles(fragment, canonical=True, isomericSmiles=True)
        for fragment in fragments
    ))


def _heavy_atom_map_pairs(molecule: Any) -> tuple[tuple[int, int], ...]:
    return tuple(
        (int(atom.GetAtomMapNum()), int(atom.GetAtomMapNum()))
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() != 1 and atom.GetAtomMapNum()
    )


def _coupled_multiplicities(multiplicities: tuple[int, ...]) -> tuple[int, ...]:
    """Return spin-coupled total multiplicities without choosing a surface."""
    if not multiplicities or any(value < 1 for value in multiplicities):
        return ()
    # Store 2S as integers. Coupling S_a and S_b admits every S between
    # |S_a-S_b| and S_a+S_b in unit steps.
    totals = {int(multiplicities[0]) - 1}
    for multiplicity in multiplicities[1:]:
        partner_spin2 = int(multiplicity) - 1
        totals = {
            coupled
            for current in totals
            for coupled in range(
                abs(current - partner_spin2), current + partner_spin2 + 1, 2,
            )
        }
    return tuple(sorted(value + 1 for value in totals))


class GraphXtbDiscoveryBackend:
    """Bounded target and declared-reactive-partner graph/xTB discovery.

    The unimolecular edit space is retained.  When the chemical-system
    contract supplies reactive partners, the backend additionally proposes
    radical association and H-atom abstraction for each declared pair.  It
    does not infer encounter rates, solvent catalysis, arbitrary multi-bond
    edits, or reaction kinetics.
    """

    name = "graph_xtb"

    def __init__(
        self,
        *,
        xtb_runner: XtbEndpointRunner | None = None,
        executable: str | None = None,
        ncores: int = 1,
        workdir: Path | None = None,
    ) -> None:
        if ncores < 1:
            raise ValueError("Graph xTB discovery needs at least one core")
        self.xtb_runner = xtb_runner or run_xtb_endpoint_screen
        self.executable = executable
        self.ncores = ncores
        self.workdir = Path(workdir) if workdir is not None else None

    @staticmethod
    def _system_value(system: str | Mapping[str, Any], key: str, default: Any) -> Any:
        return system.get(key, default) if isinstance(system, Mapping) else default

    def _enumerate(
        self,
        smiles: str,
        *,
        charge: int,
        multiplicity: int,
        target_label: str,
        environment: Mapping[str, Any] | None,
        budget: GraphXtbBudget,
    ) -> tuple[list[ReactionCandidate], list[dict[str, Any]], list[str]]:
        from rdkit import Chem

        source = Chem.MolFromSmiles(smiles)
        if source is None:
            raise ValueError(f"Invalid graph-discovery SMILES: {smiles}")
        source = Chem.Mol(source)
        for atom in source.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        source_fragments = _mapped_fragments(source)
        source_inventory, source_charge = _inventory([source])
        if source_charge != charge:
            raise ValueError("Declared graph-discovery charge does not match target SMILES")

        candidates: list[ReactionCandidate] = []
        records: list[dict[str, Any]] = []
        seen_products: set[tuple[str, ...]] = set()
        truncation: list[str] = []
        raw_count = 0

        def emit(edit: dict[str, Any], product: Any) -> bool:
            nonlocal raw_count
            if raw_count >= budget.max_generated_graph_candidates:
                if "maximum_generated_graph_candidates" not in truncation:
                    truncation.append("maximum_generated_graph_candidates")
                return False
            raw_count += 1
            record = {"raw_index": raw_count, "edit": dict(edit)}
            try:
                # Invalid graph edits are expected screening outcomes, not
                # user-facing RDKit errors.  Keep their reason in the record
                # while suppressing one warning per rejected edit.
                from rdkit import rdBase

                with rdBase.BlockLogs():
                    sanitization = Chem.SanitizeMol(product, catchErrors=True)
                if sanitization != Chem.SanitizeFlags.SANITIZE_NONE:
                    raise ValueError(f"rdkit_sanitize:{sanitization}")
                fragments = _mapped_fragments(product)
                product_molecules = list(Chem.GetMolFrags(product, asMols=True, sanitizeFrags=True))
                product_inventory, product_charge = _inventory(product_molecules)
                unmasked_fragments = tuple(sorted(canonical_smiles(item) for item in fragments))
                if source_inventory != product_inventory:
                    record.update(status="generated_but_pruned", reason="element_or_hydrogen_inventory_mismatch")
                elif product_charge != source_charge:
                    record.update(status="generated_but_pruned", reason="total_charge_mismatch")
                elif tuple(sorted(canonical_smiles(item) for item in source_fragments)) == unmasked_fragments:
                    record.update(status="generated_but_pruned", reason="products_identical_to_reactants")
                elif len(fragments) > 1 and edit["kind"] != "single_bond_break":
                    record.update(status="generated_but_pruned", reason="unexpected_fragmentation_for_edit_class")
                elif fragments in seen_products:
                    record.update(status="generated_but_pruned", reason="duplicate_graph_isomorph")
                elif len(candidates) >= budget.max_structurally_valid_candidates or len(candidates) >= budget.max_products_per_reactant_state:
                    if "maximum_structurally_valid_candidates" not in truncation:
                        truncation.append("maximum_structurally_valid_candidates")
                    record.update(status="search_truncated", reason="structural_candidate_budget_exhausted")
                else:
                    seen_products.add(fragments)
                    product_states = [endpoint_state(value) for value in fragments]
                    product_charges = tuple(
                        int(state["charge"]) for state in product_states
                        if state.get("charge") is not None
                    )
                    product_multiplicities = tuple(
                        int(state["multiplicity"]) for state in product_states
                        if state.get("multiplicity") is not None
                    )
                    candidate = ReactionCandidate.create(
                        source_fragments,
                        fragments,
                        total_charge=source_charge,
                        multiplicity=multiplicity,
                        atom_mapping=_heavy_atom_map_pairs(source),
                        discovery_assertions=({
                            "source": "graph_xtb",
                            "edit": dict(edit),
                            "structural_pruning": "passed",
                            "search_scope": "single_target_unimolecular_graph_edit",
                        },),
                        screening_status="structurally_valid_unscreened",
                        provenance={
                            "atom_mapping_representation": "smiles_atom_map_numbers",
                            "graph_enumeration": {
                                "edit": dict(edit),
                                "target_smiles": canonical_smiles(smiles),
                                "search_scope": "single_target_unimolecular_graph_edit",
                            },
                        },
                        warnings=(
                            "Connectivity proposal only; no transition-state barrier or rate was inferred from graph editing.",
                        ),
                        verification_status="candidate_unverified",
                        kinetic_relevance="unknown",
                        reactant_labels=(target_label,),
                        reactant_charges=(int(charge),),
                        product_charges=(
                            product_charges if len(product_charges) == len(fragments) else ()
                        ),
                        reactant_multiplicities=(int(multiplicity),),
                        product_multiplicities=(
                            product_multiplicities
                            if len(product_multiplicities) == len(fragments) else ()
                        ),
                        environment=dict(environment or {}),
                    )
                    candidates.append(candidate)
                    record.update(
                        status="structurally_valid",
                        candidate_id=candidate.candidate_id,
                        reactants=list(source_fragments),
                        products=list(fragments),
                    )
            except Exception as error:
                record.update(status="generated_but_pruned", reason=f"structural_validation_failed:{type(error).__name__}:{error}")
            records.append(record)
            return not truncation

        heavy = [atom.GetIdx() for atom in source.GetAtoms() if atom.GetAtomicNum() != 1]

        def ordinary_four_carbon_sites(indices: tuple[int, ...]) -> bool:
            if len(set(indices)) != 4:
                return False
            for index in indices:
                atom = source.GetAtomWithIdx(index)
                if (atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.GetFormalCharge()
                        or atom.GetIsotope() or atom.GetNumRadicalElectrons()
                        or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED):
                    return False
            # Do not erase specified reacting stereochemistry to invent a path.
            for bond in source.GetBonds():
                if bond.GetBeginAtomIdx() in indices and bond.GetEndAtomIdx() in indices:
                    if bond.GetStereo() != Chem.BondStereo.STEREONONE:
                        return False
            return True

        if budget.allow_electrocyclic_bond_rearrangement:
            # A connectivity hypothesis only: neither orbital allowedness nor a
            # stereochemical branch, saddle, barrier or rate is inferred here.
            for ring in source.GetRingInfo().AtomRings():
                if len(ring) != 4 or not ordinary_four_carbon_sites(tuple(ring)):
                    continue
                if any(source.GetRingInfo().NumAtomRings(index) != 1 for index in ring):
                    continue
                bonds = [source.GetBondBetweenAtoms(ring[i], ring[(i + 1) % 4]) for i in range(4)]
                doubles = [i for i, bond in enumerate(bonds) if bond.GetBondType() == Chem.BondType.DOUBLE]
                if len(doubles) != 1 or any(bond.GetIsAromatic() or bond.GetBondType() not in {
                    Chem.BondType.SINGLE, Chem.BondType.DOUBLE} for bond in bonds):
                    continue
                offset = doubles[0]
                a, b, c, d = (ring[(offset + i) % 4] for i in range(4))
                if (any(source.GetAtomWithIdx(index).GetHybridization() != Chem.HybridizationType.SP2 for index in (a, b))
                        or any(source.GetAtomWithIdx(index).GetHybridization() != Chem.HybridizationType.SP3 for index in (c, d))):
                    continue
                edited = Chem.RWMol(source)
                edited.RemoveBond(c, d)
                for left, right, kind in ((a, b, Chem.BondType.SINGLE),
                                          (b, c, Chem.BondType.DOUBLE), (d, a, Chem.BondType.DOUBLE)):
                    edited.GetBondBetweenAtoms(left, right).SetBondType(kind)
                if not emit({"kind": "four_pi_electrocyclic_ring_opening", "reactive_atoms": [a, b, c, d],
                             "broken_bond": [c, d], "bond_order_changes": [[a, b, 2, 1], [b, c, 1, 2], [d, a, 1, 2]]},
                            edited.GetMol()):
                    return candidates, records, truncation

            for central in source.GetBonds():
                if central.GetBondType() != Chem.BondType.SINGLE or central.GetIsAromatic():
                    continue
                b, c = central.GetBeginAtomIdx(), central.GetEndAtomIdx()
                left = [atom.GetIdx() for atom in source.GetAtomWithIdx(b).GetNeighbors()
                        if atom.GetIdx() != c and source.GetBondBetweenAtoms(atom.GetIdx(), b).GetBondType() == Chem.BondType.DOUBLE]
                right = [atom.GetIdx() for atom in source.GetAtomWithIdx(c).GetNeighbors()
                         if atom.GetIdx() != b and source.GetBondBetweenAtoms(c, atom.GetIdx()).GetBondType() == Chem.BondType.DOUBLE]
                for a in left:
                    for d in right:
                        sites = (a, b, c, d)
                        if (not ordinary_four_carbon_sites(sites)
                                or any(source.GetAtomWithIdx(index).IsInRing() for index in sites)
                                or any(source.GetAtomWithIdx(index).GetHybridization() != Chem.HybridizationType.SP2 for index in sites)
                                or source.GetBondBetweenAtoms(a, d) is not None):
                            continue
                        edited = Chem.RWMol(source)
                        for i, j, kind in ((a, b, Chem.BondType.SINGLE),
                                           (b, c, Chem.BondType.DOUBLE), (c, d, Chem.BondType.SINGLE)):
                            edited.GetBondBetweenAtoms(i, j).SetBondType(kind)
                        edited.AddBond(d, a, Chem.BondType.SINGLE)
                        if not emit({"kind": "four_pi_electrocyclic_ring_closure", "reactive_atoms": list(sites),
                                     "formed_bond": [d, a], "bond_order_changes": [[a, b, 2, 1], [b, c, 1, 2], [c, d, 2, 1]]},
                                    edited.GetMol()):
                            return candidates, records, truncation

        def homolytic_edit(
            broken_bonds: tuple[tuple[int, int], ...],
            formed_bonds: tuple[tuple[int, int], ...] = (),
        ) -> Any:
            """Apply single-bond edits without silently creating/deleting H.

            RDKit otherwise fills a valence exposed by a bond deletion with an
            implicit hydrogen (for example, ethane -> two methane molecules).
            A standalone single-bond break is represented as homolysis in this
            deliberately restricted first backend, so the exposed atoms retain
            their hydrogen inventory and acquire the corresponding radical
            electron.  Net bond-rearrangement changes retain that treatment
            only on atoms with an uncompensated broken bond.
            """
            edited = Chem.RWMol(source)
            for left, right in broken_bonds:
                edited.RemoveBond(left, right)
            for left, right in formed_bonds:
                edited.AddBond(left, right, Chem.BondType.SINGLE)
            net_bond_loss: dict[int, int] = {}
            for left, right in broken_bonds:
                net_bond_loss[left] = net_bond_loss.get(left, 0) + 1
                net_bond_loss[right] = net_bond_loss.get(right, 0) + 1
            for left, right in formed_bonds:
                net_bond_loss[left] = net_bond_loss.get(left, 0) - 1
                net_bond_loss[right] = net_bond_loss.get(right, 0) - 1
            for index, loss in net_bond_loss.items():
                if loss > 0:
                    atom = edited.GetAtomWithIdx(index)
                    implicit_hydrogens = atom.GetNumImplicitHs()
                    atom.SetNoImplicit(True)
                    # Turning off implicit H without retaining the current
                    # count would transform CH3 into bare C.  Retain those
                    # nuclei explicitly before marking the homolytic radical.
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() + implicit_hydrogens)
                    atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons() + loss)
            return edited.GetMol()

        if budget.allow_dissociation:
            for bond in sorted(source.GetBonds(), key=lambda item: (item.GetBeginAtomIdx(), item.GetEndAtomIdx())):
                left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if bond.GetIsAromatic() or bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                if left not in heavy or right not in heavy:
                    continue
                if not emit(
                    {"kind": "single_bond_break", "broken_bond": [left, right]},
                    homolytic_edit(((left, right),)),
                ):
                    return candidates, records, truncation

        if budget.allow_bond_formation:
            for position, left in enumerate(heavy):
                for right in heavy[position + 1:]:
                    if source.GetBondBetweenAtoms(left, right) is not None:
                        continue
                    if source.GetAtomWithIdx(left).GetIsAromatic() or source.GetAtomWithIdx(right).GetIsAromatic():
                        continue
                    edited = Chem.RWMol(source)
                    edited.AddBond(left, right, Chem.BondType.SINGLE)
                    if not emit({"kind": "single_bond_formation", "formed_bond": [left, right]}, edited.GetMol()):
                        return candidates, records, truncation

        if budget.allow_bond_rearrangement:
            bonds = [
                (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                for bond in source.GetBonds()
                if not bond.GetIsAromatic() and bond.GetBondType() == Chem.BondType.SINGLE
                and bond.GetBeginAtomIdx() in heavy and bond.GetEndAtomIdx() in heavy
            ]
            nonbonds = [
                (left, right)
                for position, left in enumerate(heavy)
                for right in heavy[position + 1:]
                if source.GetBondBetweenAtoms(left, right) is None
                and not source.GetAtomWithIdx(left).GetIsAromatic()
                and not source.GetAtomWithIdx(right).GetIsAromatic()
            ]
            for broken_left, broken_right in bonds:
                for formed_left, formed_right in nonbonds:
                    if {broken_left, broken_right} == {formed_left, formed_right}:
                        continue
                    if not emit({
                        "kind": "single_bond_break_and_formation",
                        "broken_bond": [broken_left, broken_right],
                        "formed_bond": [formed_left, formed_right],
                    }, homolytic_edit(
                        ((broken_left, broken_right),), ((formed_left, formed_right),),
                    )):
                        return candidates, records, truncation

        if budget.allow_hydrogen_transfer:
            explicit = Chem.AddHs(source)
            for atom in explicit.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    atom.SetAtomMapNum(10000 + atom.GetIdx())
            heavy_explicit = [atom.GetIdx() for atom in explicit.GetAtoms() if atom.GetAtomicNum() != 1]
            for hydrogen in [atom for atom in explicit.GetAtoms() if atom.GetAtomicNum() == 1]:
                neighbours = [item.GetIdx() for item in hydrogen.GetNeighbors() if item.GetAtomicNum() != 1]
                if len(neighbours) != 1:
                    continue
                donor = neighbours[0]
                for acceptor in heavy_explicit:
                    if acceptor == donor:
                        continue
                    edited = Chem.RWMol(explicit)
                    edited.RemoveBond(donor, hydrogen.GetIdx())
                    edited.AddBond(acceptor, hydrogen.GetIdx(), Chem.BondType.SINGLE)
                    donor_atom = edited.GetAtomWithIdx(donor)
                    donor_atom.SetNoImplicit(True)
                    donor_atom.SetNumRadicalElectrons(donor_atom.GetNumRadicalElectrons() + 1)
                    acceptor_atom = edited.GetAtomWithIdx(acceptor)
                    if acceptor_atom.GetNumRadicalElectrons() > 0:
                        acceptor_atom.SetNumRadicalElectrons(acceptor_atom.GetNumRadicalElectrons() - 1)
                    # Keep the transferred H explicit.  Removing it would let
                    # RDKit recreate implicit hydrogens according to an
                    # arbitrary valence completion and erase the graph edit.
                    product = edited.GetMol()
                    if not emit({
                        "kind": "hydrogen_transfer",
                        "hydrogen_atom": hydrogen.GetIdx(),
                        "donor_atom": donor,
                        "acceptor_atom": acceptor,
                    }, product):
                        return candidates, records, truncation
        return candidates, records, truncation

    def _enumerate_bimolecular(
        self,
        smiles: str,
        *,
        charge: int,
        multiplicity: int,
        target_label: str,
        reactive_partners: tuple[Mapping[str, Any], ...],
        environment: Mapping[str, Any] | None,
        budget: GraphXtbBudget,
    ) -> tuple[
        list[ReactionCandidate], list[dict[str, Any]], list[str], list[dict[str, Any]],
    ]:
        """Enumerate a deliberately small condition-specific pairwise edit space."""
        from rdkit import Chem, rdBase

        candidates: list[ReactionCandidate] = []
        records: list[dict[str, Any]] = []
        truncation: list[str] = []
        partner_assessments: list[dict[str, Any]] = []
        seen: set[tuple[str, tuple[str, ...], int]] = set()
        raw_count = 0

        if len(reactive_partners) > budget.max_reactive_partners:
            truncation.append("maximum_reactive_partners")
        selected_partners = reactive_partners[:budget.max_reactive_partners]

        target = Chem.MolFromSmiles(smiles)
        if target is None:
            raise ValueError(f"Invalid graph-discovery target SMILES: {smiles}")
        target_state = endpoint_state(smiles)
        if target_state.get("charge") != charge:
            raise ValueError("Declared graph-discovery charge does not match target SMILES")
        target_explicit = Chem.AddHs(Chem.Mol(target))
        for atom in target_explicit.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)

        for partner_index, partner_record in enumerate(selected_partners):
            label = str(partner_record.get("label") or f"reactive-partner-{partner_index + 1}")
            partner_smiles = str(partner_record.get("smiles") or "")
            assessment = {
                "label": label,
                "smiles": partner_smiles,
                "status": "unresolved",
                "attempted_reaction_classes": [],
                "raw_generation_attempts": 0,
                "structurally_valid_candidates": 0,
            }
            partner = Chem.MolFromSmiles(partner_smiles) if partner_smiles else None
            if partner is None:
                assessment.update(
                    status="invalid_partner_smiles",
                    reason="Reactive partner does not have a valid SMILES graph.",
                )
                partner_assessments.append(assessment)
                continue
            inferred = endpoint_state(partner_smiles)
            partner_charge = partner_record.get("charge", inferred.get("charge"))
            partner_multiplicity = partner_record.get(
                "multiplicity", inferred.get("multiplicity"),
            )
            if partner_charge is None or partner_multiplicity is None:
                assessment.update(
                    status="electronic_state_unresolved",
                    reason=(
                        "Reactive-partner charge or multiplicity is not uniquely encoded; "
                        "declare both in the condition contract."
                    ),
                )
                partner_assessments.append(assessment)
                continue
            partner_charge = int(partner_charge)
            partner_multiplicity = int(partner_multiplicity)
            if partner_multiplicity < 1:
                assessment.update(
                    status="electronic_state_invalid",
                    reason="Declared reactive-partner multiplicity must be positive.",
                )
                partner_assessments.append(assessment)
                continue
            if inferred.get("charge") is not None and int(inferred["charge"]) != partner_charge:
                assessment.update(
                    status="electronic_state_mismatch",
                    reason="Declared reactive-partner charge disagrees with its SMILES graph.",
                )
                partner_assessments.append(assessment)
                continue

            partner_explicit = Chem.AddHs(Chem.Mol(partner))
            map_offset = 1000 * (partner_index + 1)
            for atom in partner_explicit.GetAtoms():
                atom.SetAtomMapNum(map_offset + atom.GetIdx() + 1)
            combined = Chem.CombineMols(target_explicit, partner_explicit)
            target_atom_count = target_explicit.GetNumAtoms()
            target_indices = range(target_atom_count)
            partner_indices = range(target_atom_count, combined.GetNumAtoms())
            target_heavy = [
                index for index in target_indices
                if combined.GetAtomWithIdx(index).GetAtomicNum() != 1
            ]
            partner_heavy = [
                index for index in partner_indices
                if combined.GetAtomWithIdx(index).GetAtomicNum() != 1
            ]
            target_hydrogens = [
                index for index in target_indices
                if combined.GetAtomWithIdx(index).GetAtomicNum() == 1
            ]
            partner_hydrogens = [
                index for index in partner_indices
                if combined.GetAtomWithIdx(index).GetAtomicNum() == 1
            ]
            reactants = (
                Chem.MolToSmiles(target_explicit, canonical=True, isomericSmiles=True),
                Chem.MolToSmiles(partner_explicit, canonical=True, isomericSmiles=True),
            )
            source_inventory, source_charge = _inventory([target_explicit, partner_explicit])
            reactant_multiplicities = (multiplicity, partner_multiplicity)
            reactant_total_multiplicities = _coupled_multiplicities(
                reactant_multiplicities,
            )
            partner_candidate_count = 0
            partner_raw_start = raw_count
            pair_has_radical_site = any(
                combined.GetAtomWithIdx(index).GetNumRadicalElectrons() > 0
                for index in (*target_heavy, *partner_heavy)
            )

            def emit(edit: dict[str, Any], product: Any) -> bool:
                nonlocal raw_count, partner_candidate_count
                if raw_count >= budget.max_bimolecular_generated_candidates:
                    if "maximum_bimolecular_generated_candidates" not in truncation:
                        truncation.append("maximum_bimolecular_generated_candidates")
                    return False
                raw_count += 1
                record = {
                    "raw_index": raw_count,
                    "scope": "condition_specific_bimolecular_graph_edit",
                    "reactive_partner": label,
                    "edit": dict(edit),
                }
                try:
                    with rdBase.BlockLogs():
                        sanitization = Chem.SanitizeMol(product, catchErrors=True)
                    if sanitization != Chem.SanitizeFlags.SANITIZE_NONE:
                        raise ValueError(f"rdkit_sanitize:{sanitization}")
                    fragments = _mapped_fragments(product)
                    product_molecules = list(
                        Chem.GetMolFrags(product, asMols=True, sanitizeFrags=True)
                    )
                    product_inventory, product_charge = _inventory(product_molecules)
                    unmasked_reactants = tuple(sorted(canonical_smiles(item) for item in reactants))
                    unmasked_products = tuple(sorted(canonical_smiles(item) for item in fragments))
                    if source_inventory != product_inventory:
                        raise ValueError("element_or_hydrogen_inventory_mismatch")
                    if product_charge != source_charge:
                        raise ValueError("total_charge_mismatch")
                    if unmasked_reactants == unmasked_products:
                        raise ValueError("products_identical_to_reactants")
                    product_states = tuple(endpoint_state(value) for value in fragments)
                    if any(
                        state.get("charge") is None or state.get("multiplicity") is None
                        for state in product_states
                    ):
                        raise ValueError("product_electronic_state_unresolved")
                    product_charges = tuple(int(state["charge"]) for state in product_states)
                    product_multiplicities = tuple(
                        int(state["multiplicity"]) for state in product_states
                    )
                    compatible_total_multiplicities = tuple(sorted(
                        set(reactant_total_multiplicities)
                        & set(_coupled_multiplicities(product_multiplicities))
                    ))
                    if not compatible_total_multiplicities:
                        raise ValueError("reactant_product_spin_surfaces_do_not_intersect")
                    if (
                        partner_candidate_count >= budget.max_bimolecular_structurally_valid_candidates
                    ):
                        if "maximum_bimolecular_structurally_valid_candidates" not in truncation:
                            truncation.append(
                                "maximum_bimolecular_structurally_valid_candidates"
                            )
                        record.update(
                            status="search_truncated",
                            reason="bimolecular_structural_candidate_budget_exhausted",
                        )
                        records.append(record)
                        return False
                    emitted_ids = []
                    for total_multiplicity in compatible_total_multiplicities:
                        identity = (label, unmasked_products, total_multiplicity)
                        if identity in seen:
                            continue
                        seen.add(identity)
                        candidate = ReactionCandidate.create(
                            reactants,
                            fragments,
                            total_charge=source_charge,
                            multiplicity=total_multiplicity,
                            atom_mapping=_heavy_atom_map_pairs(combined),
                            discovery_assertions=({
                                "source": "graph_xtb",
                                "edit": dict(edit),
                                "structural_pruning": "passed",
                                "search_scope": (
                                    "condition_specific_bimolecular_radical_association_"
                                    "and_hydrogen_abstraction"
                                ),
                                "reactive_partner": label,
                            },),
                            screening_status="structurally_valid_unscreened",
                            provenance={
                                "atom_mapping_representation": "smiles_atom_map_numbers",
                                "graph_enumeration": {
                                    "edit": dict(edit),
                                    "target_smiles": canonical_smiles(smiles),
                                    "reactive_partner": {
                                        "label": label,
                                        "smiles": canonical_smiles(partner_smiles),
                                        "charge": partner_charge,
                                        "multiplicity": partner_multiplicity,
                                    },
                                    "search_scope": (
                                        "condition_specific_bimolecular_radical_"
                                        "association_and_hydrogen_abstraction"
                                    ),
                                    "compatible_total_multiplicity": total_multiplicity,
                                },
                            },
                            warnings=(
                                "Pairwise connectivity proposal only; no encounter, capture, barrier, or rate was inferred.",
                            ),
                            verification_status="candidate_unverified",
                            kinetic_relevance="unknown",
                            reactant_labels=(target_label, label),
                            reactant_charges=(charge, partner_charge),
                            product_charges=product_charges,
                            reactant_multiplicities=reactant_multiplicities,
                            product_multiplicities=product_multiplicities,
                            environment=dict(environment or {}),
                        )
                        candidates.append(candidate)
                        emitted_ids.append(candidate.candidate_id)
                        partner_candidate_count += 1
                    if emitted_ids:
                        record.update(
                            status="structurally_valid",
                            candidate_ids=emitted_ids,
                            reactants=list(reactants),
                            products=list(fragments),
                            compatible_total_multiplicities=list(
                                compatible_total_multiplicities
                            ),
                        )
                    else:
                        record.update(
                            status="generated_but_pruned",
                            reason="duplicate_graph_and_spin_state",
                        )
                except Exception as error:
                    record.update(
                        status="generated_but_pruned",
                        reason=(
                            f"structural_validation_failed:{type(error).__name__}:{error}"
                        ),
                    )
                records.append(record)
                return not truncation

            if budget.allow_bimolecular_radical_association:
                assessment["attempted_reaction_classes"].append(
                    "bimolecular_radical_association"
                )
                for left in target_heavy:
                    for right in partner_heavy:
                        left_atom = combined.GetAtomWithIdx(left)
                        right_atom = combined.GetAtomWithIdx(right)
                        if (
                            left_atom.GetNumRadicalElectrons() < 1
                            or right_atom.GetNumRadicalElectrons() < 1
                        ):
                            continue
                        edited = Chem.RWMol(combined)
                        edited.AddBond(left, right, Chem.BondType.SINGLE)
                        for index in (left, right):
                            atom = edited.GetAtomWithIdx(index)
                            atom.SetNumRadicalElectrons(
                                atom.GetNumRadicalElectrons() - 1
                            )
                        if not emit({
                            "kind": "bimolecular_radical_association",
                            "formed_bond": [left, right],
                        }, edited.GetMol()):
                            break
                    if truncation:
                        break

            if budget.allow_bimolecular_hydrogen_abstraction and not truncation:
                assessment["attempted_reaction_classes"].append(
                    "bimolecular_hydrogen_abstraction"
                )
                directions = (
                    (target_hydrogens, partner_heavy, "target_to_partner"),
                    (partner_hydrogens, target_heavy, "partner_to_target"),
                )
                for hydrogens, acceptors, direction in directions:
                    for hydrogen_index in hydrogens:
                        hydrogen = combined.GetAtomWithIdx(hydrogen_index)
                        neighbours = [atom.GetIdx() for atom in hydrogen.GetNeighbors()]
                        if len(neighbours) != 1:
                            continue
                        donor = neighbours[0]
                        for acceptor in acceptors:
                            acceptor_atom = combined.GetAtomWithIdx(acceptor)
                            if acceptor_atom.GetNumRadicalElectrons() < 1:
                                continue
                            edited = Chem.RWMol(combined)
                            edited.RemoveBond(donor, hydrogen_index)
                            edited.AddBond(acceptor, hydrogen_index, Chem.BondType.SINGLE)
                            donor_atom = edited.GetAtomWithIdx(donor)
                            donor_atom.SetNumRadicalElectrons(
                                donor_atom.GetNumRadicalElectrons() + 1
                            )
                            edited_acceptor = edited.GetAtomWithIdx(acceptor)
                            edited_acceptor.SetNumRadicalElectrons(
                                edited_acceptor.GetNumRadicalElectrons() - 1
                            )
                            if not emit({
                                "kind": "bimolecular_hydrogen_abstraction",
                                "direction": direction,
                                "hydrogen_atom": hydrogen_index,
                                "donor_atom": donor,
                                "acceptor_atom": acceptor,
                            }, edited.GetMol()):
                                break
                        if truncation:
                            break
                    if truncation:
                        break

            partner_attempt_count = raw_count - partner_raw_start
            if truncation:
                partner_status = "search_truncated"
                partner_reason = None
            elif not pair_has_radical_site:
                partner_status = "outside_supported_pairwise_reaction_classes"
                partner_reason = (
                    "Both reactants are closed-shell, while this backend currently "
                    "supports only radical association and radical H abstraction."
                )
            else:
                partner_status = "completed_declared_search"
                partner_reason = None
            assessment.update(
                status=partner_status,
                raw_generation_attempts=partner_attempt_count,
                structurally_valid_candidates=partner_candidate_count,
            )
            if partner_reason is not None:
                assessment["reason"] = partner_reason
            partner_assessments.append(assessment)
            if truncation:
                break

        return candidates, records, truncation, partner_assessments

    def _screen_candidate(
        self,
        candidate: ReactionCandidate,
        *,
        candidate_dir: Path,
        budget: GraphXtbBudget,
        optimization_counter: list[int],
        failure_counter: list[int],
        truncation: list[str],
        deadline: float | None = None,
    ) -> ReactionCandidate:
        endpoints = [("reactant", index, smiles) for index, smiles in enumerate(candidate.reactants)] + [
            ("product", index, smiles) for index, smiles in enumerate(candidate.products)
        ]
        if optimization_counter[0] + len(endpoints) > budget.max_xtb_optimizations:
            if "maximum_xtb_optimizations" not in truncation:
                truncation.append("maximum_xtb_optimizations")
            return replace(
                candidate,
                screening_status="search_truncated",
                warnings=tuple(sorted(set((*candidate.warnings, "xTB endpoint budget exhausted before screening.")))),
            )
        if failure_counter[0] >= budget.max_xtb_failures:
            if "maximum_xtb_failures" not in truncation:
                truncation.append("maximum_xtb_failures")
            return replace(
                candidate,
                screening_status="search_truncated",
                warnings=tuple(sorted(set((*candidate.warnings, "xTB failure budget exhausted before screening.")))),
            )
        records: list[dict[str, Any]] = []
        screening_environment = _screening_environment(candidate.environment)
        reactant_energy = 0.0
        product_energy = 0.0
        for side, index, smiles in endpoints:
            if deadline is not None and time.monotonic() >= deadline:
                if "overall_discovery_timeout" not in truncation:
                    truncation.append("overall_discovery_timeout")
                return replace(candidate, screening_status="search_truncated",
                               provenance={**candidate.provenance, "xtb_screen": {"status": "search_truncated", "endpoint_records": records}},
                               warnings=tuple(sorted(set((*candidate.warnings, "Aggregate discovery deadline exhausted; completed endpoint evidence retained without rate credit.")))))
            state = endpoint_state(smiles)
            supplied_charges = (
                candidate.reactant_charges
                if side == "reactant" else candidate.product_charges
            )
            supplied_multiplicities = (
                candidate.reactant_multiplicities
                if side == "reactant" else candidate.product_multiplicities
            )
            endpoint_charge = (
                supplied_charges[index]
                if len(supplied_charges) == len(candidate.reactants if side == "reactant" else candidate.products)
                else state.get("charge")
            )
            endpoint_multiplicity = (
                supplied_multiplicities[index]
                if len(supplied_multiplicities) == len(candidate.reactants if side == "reactant" else candidate.products)
                else state.get("multiplicity")
            )
            if endpoint_charge is None or endpoint_multiplicity is None:
                result = {
                    "status": "xtb_failed",
                    "failure_reason": "endpoint_charge_or_multiplicity_unresolved_from_smiles",
                    "smiles": smiles,
                }
            else:
                optimization_counter[0] += 1
                solvent_kwargs = {"environment": candidate.environment} if screening_environment["phase"] == "solution" else {}
                result = self.xtb_runner(
                    smiles,
                    candidate_dir / f"{side}-{index + 1:02d}",
                    charge=int(endpoint_charge),
                    multiplicity=int(endpoint_multiplicity),
                    timeout_seconds=min(budget.timeout_seconds, max(0.001, deadline - time.monotonic()))
                    if deadline is not None else budget.timeout_seconds,
                    executable=self.executable,
                    ncores=self.ncores,
                    **solvent_kwargs,
                )
            result = dict(result)
            if screening_environment["phase"] == "solution" and result.get("screening_environment") != screening_environment:
                result = {**result, "status": "xtb_failed", "failure_reason": "Missing or mismatched observed solvent-screen evidence"}
            result.update({"side": side, "endpoint_index": index, "smiles": smiles})
            records.append(result)
            energy = result.get("energy_hartree")
            if (result.get("status") != "completed" or isinstance(energy, bool)
                    or not isinstance(energy, (int, float)) or not math.isfinite(energy)):
                failure_counter[0] += 1
                provenance = {
                    **candidate.provenance,
                    "xtb_screen": {
                        "status": "xtb_failed",
                        "method": "GFN2-xTB",
                        "endpoint_records": records,
                    },
                }
                return replace(
                    candidate,
                    screening_status="xtb_failed",
                    provenance=provenance,
                    warnings=tuple(sorted(set((*candidate.warnings, "xTB endpoint screening failed; candidate is not eligible for automatic verification.")))),
                )
            if side == "reactant":
                reactant_energy += float(result["energy_hartree"])
            else:
                product_energy += float(result["energy_hartree"])
        delta = (product_energy - reactant_energy) * HARTREE_TO_KCAL_MOL
        screen = {
            "status": "completed",
            "method": "GFN2-xTB",
            "endpoint_records": records,
            "screening_environment": screening_environment,
            "reactant_electronic_energy_hartree": reactant_energy,
            "product_electronic_energy_hartree": product_energy,
            "estimated_reaction_energy_kcal_mol": delta,
            "energy_window_kcal_mol": budget.energy_window_kcal_mol,
            "interpretation": (
                "Endpoint electronic-energy screen only. It does not establish a transition-state barrier, "
                "a rate coefficient, or kinetic relevance."
                + (" ALPB totals include solvation terms in the explicit gsolv reference state; they are not pure electronic energies, thermal free energies or SMD/CPCM rate evidence."
                   if screening_environment["phase"] == "solution" else "")
                + (
                    " For this bimolecular proposal, independently optimized fragment "
                    "energies also omit encounter, standard-state, capture, and entropy effects."
                    if len(candidate.reactants) > 1 else ""
                )
            ),
        }
        if screening_environment["phase"] == "solution":
            screen["reactant_screen_total_energy_hartree"] = screen.pop("reactant_electronic_energy_hartree")
            screen["product_screen_total_energy_hartree"] = screen.pop("product_electronic_energy_hartree")
        if delta > budget.energy_window_kcal_mol:
            return replace(
                candidate,
                estimated_reaction_energy_kcal_mol=delta,
                screening_status="xtb_screened_out",
                provenance={**candidate.provenance, "xtb_screen": screen},
                warnings=tuple(sorted(set((*candidate.warnings, "Candidate exceeds configured xTB endergonic screening window.")))),
            )
        return replace(
            candidate,
            estimated_reaction_energy_kcal_mol=delta,
            screening_status="candidate_unverified",
            provenance={**candidate.provenance, "xtb_screen": screen},
            warnings=tuple(sorted(set((*candidate.warnings, "xTB endpoint energy retained for screening only; kinetic relevance remains unknown.")))),
        )

    def discover(
        self,
        system: str | Mapping[str, Any],
        conditions: Mapping[str, Any] | None,
        budget: GraphXtbBudget | Mapping[str, Any] | None = None,
        *,
        workdir: Path | None = None,
    ) -> DiscoveryResult:
        search_budget = GraphXtbBudget.coerce(budget)
        deadline = time.monotonic() + search_budget.overall_timeout_seconds if search_budget.overall_timeout_seconds is not None else None
        smiles = str(system if isinstance(system, str) else system.get("smiles") or system.get("target_smiles") or "")
        if not smiles:
            raise ValueError("Graph xTB discovery needs target SMILES")
        inferred_state = endpoint_state(smiles)
        declared_charge = self._system_value(system, "charge", inferred_state.get("charge"))
        declared_multiplicity = self._system_value(
            system, "multiplicity", inferred_state.get("multiplicity"),
        )
        if declared_charge is None:
            raise ValueError(
                "Graph xTB discovery could not infer target charge; declare it in the chemical-system contract"
            )
        if declared_multiplicity is None:
            raise ValueError(
                "Graph xTB discovery could not infer a unique target multiplicity; declare it in the chemical-system contract"
            )
        charge = int(declared_charge)
        multiplicity = int(declared_multiplicity)
        if multiplicity < 1:
            raise ValueError("Graph xTB target multiplicity must be positive")
        target_label = str(self._system_value(system, "target_label", "stability"))
        environment = self._system_value(system, "environment", {})
        if not isinstance(environment, Mapping):
            raise ValueError("Graph xTB discovery environment metadata must be a mapping")
        screening_environment = _screening_environment(environment)
        raw_partners = self._system_value(system, "reactive_partners", ())
        if not isinstance(raw_partners, (list, tuple)) or any(
            not isinstance(item, Mapping) for item in raw_partners
        ):
            raise ValueError("Graph xTB reactive_partners must be a sequence of mappings")
        reactive_partners = tuple(dict(item) for item in raw_partners)
        root = Path(workdir or self.workdir or Path.cwd() / "graph-xtb-discovery")
        root.mkdir(parents=True, exist_ok=True)
        try:
            candidates, records, truncation = self._enumerate(
                smiles, charge=charge, multiplicity=multiplicity,
                target_label=target_label, environment=environment, budget=search_budget,
            )
            (
                bimolecular_candidates,
                bimolecular_records,
                bimolecular_truncation,
                partner_assessments,
            ) = self._enumerate_bimolecular(
                smiles,
                charge=charge,
                multiplicity=multiplicity,
                target_label=target_label,
                reactive_partners=reactive_partners,
                environment=environment,
                budget=search_budget,
            )
            # Condition-specific candidates are screened first so a large
            # intramolecular edit space cannot silently consume the entire xTB
            # budget before declared environmental chemistry is considered.
            candidates = [*bimolecular_candidates, *candidates]
            records = [*bimolecular_records, *records]
            truncation = list(dict.fromkeys([
                *bimolecular_truncation, *truncation,
            ]))
        except ImportError as error:
            result = DiscoveryResult(
                backend=self.name,
                status="backend_unavailable",
                provenance={
                    "reason": f"RDKit is required for graph enumeration: {error}",
                    "budget": search_budget.as_dict(),
                    "search_scope": (
                        "single_target_unimolecular_and_condition_specific_"
                        "bimolecular_graph_edit"
                        if reactive_partners else "single_target_unimolecular_graph_edit"
                    ),
                },
                warnings=("No graph candidates were generated because the RDKit backend is unavailable.",),
            )
            _write_json(root / "graph-xtb-discovery.json", result.as_dict())
            return result
        optimization_counter = [0]
        failure_counter = [0]
        screened = [
            self._screen_candidate(
                candidate,
                candidate_dir=root / candidate.candidate_id,
                budget=search_budget,
                optimization_counter=optimization_counter,
                failure_counter=failure_counter,
                truncation=truncation,
                deadline=deadline,
            )
            for candidate in candidates
        ]
        status = "discovery_truncated" if truncation else "completed"
        result = DiscoveryResult(
            backend=self.name,
            status=status,
            candidates=tuple(screened),
            records=tuple(records),
            provenance={
                "schema_version": GRAPH_XTB_DISCOVERY_SCHEMA_VERSION,
                "target_smiles": canonical_smiles(smiles),
                "target_charge": charge,
                "target_multiplicity": multiplicity,
                "condition_contract_id": (conditions or {}).get("contract_id"),
                "screening_environment": screening_environment,
                "budget": search_budget.as_dict(),
                "search_scope": (
                    "single_target_unimolecular_and_condition_specific_"
                    "bimolecular_graph_edit"
                    if reactive_partners else "single_target_unimolecular_graph_edit"
                ),
                "reaction_classes": [
                    "single_bond_break", "single_bond_formation",
                    "single_bond_break_and_formation", "hydrogen_transfer",
                    *(["four_pi_electrocyclic_ring_opening", "four_pi_electrocyclic_ring_closure"]
                      if search_budget.allow_electrocyclic_bond_rearrangement else []),
                    *(
                        [
                            "bimolecular_radical_association",
                            "bimolecular_hydrogen_abstraction",
                        ] if reactive_partners else []
                    ),
                ],
                "reactive_partner_assessments": partner_assessments,
                "raw_generation_attempts": len(records),
                "structurally_valid_candidates": len(candidates),
                "xtb_endpoint_optimizations": optimization_counter[0],
                "xtb_failures": failure_counter[0],
                "truncation_reasons": truncation,
                "no_candidate_interpretation": (
                    "No additional candidate reaction was found within the configured graph-enumeration and xTB screening bounds."
                ),
            },
            warnings=(
                "Graph enumeration is bounded and incomplete; it does not establish reaction-network completeness.",
                "xTB supplies endpoint screening only and is not a reaction predictor or kinetic model.",
                "Four-pi edits cover only isolated ordinary carbon cyclobutenes and acyclic carbon conjugated dienes without specified reacting stereochemistry; orbital/TS validity and broader pericyclic coverage are unresolved.",
            ),
        )
        _write_json(root / "graph-xtb-discovery.json", result.as_dict())
        return result


__all__ = [
    "GRAPH_XTB_DISCOVERY_SCHEMA_VERSION", "GraphXtbBudget", "GraphXtbDiscoveryBackend",
    "run_xtb_endpoint_screen", "validate_graph_xtb_candidate_sources",
]
