"""Narrow, hash-bound recovery of a demonstrated ORCA 6.1 input dialect bug."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re


NEB_PREOPT_REPAIR_POLICY = "orca-6.1-neb-preopt-parser-repair-v1"
NEB_RESTART_PREOPT_REPAIR_POLICY = "orca-6.1-neb-restart-preopt-parser-repair-v1"
NEB_MANUAL = "https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/neb.html"


def _artifact(path):
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def inspect_neb_preopt_parser_failure(input_path: Path) -> dict | None:
    """No calculation: require exact syntax failure, immutable inputs, dead lease.

    This is not a generic retry of ORCA errors. A completed electronic or NEB
    calculation, changed dependency or indeterminate process is ineligible.
    """
    inp = Path(input_path)
    output, contract_path, lease_path = (inp.with_suffix(suffix) for suffix in
                                         (".out", ".contract.json", ".process.json"))
    if not all(p.is_file() for p in (inp, output, contract_path, lease_path)):
        return None
    text = inp.read_text()
    out = output.read_text(errors="replace")
    if ("NEB-TS" not in text.upper()
            or len(re.findall(r"(?im)^\s*PreOptEnds\s+true\s*$", text)) != 1
            or not re.search(r"Program Version\s+6\.1\.\d+", out)
            or not re.search(r"Unknown identifier in NEB block[^\n]*\n\s*Last token\s*:\s*PREOPTENDS", out, re.I)
            or re.search(r"ORCA TERMINATED NORMALLY|SCF CONVERGED|FINAL SINGLE POINT ENERGY|NEB ITERATION", out, re.I)):
        return None
    contract = json.loads(contract_path.read_text())
    lease = json.loads(lease_path.read_text())
    if lease.get("status") != "terminated" or type(lease.get("returncode")) is not int \
            or lease["returncode"] == 0:
        return None
    for key, group in (("pid", False), ("process_group_id", True)):
        identifier = lease.get(key)
        if type(identifier) is not int or identifier <= 0:
            return None
        try:
            os.kill(-identifier if group else identifier, 0)
        except ProcessLookupError:
            pass
        except (PermissionError, OSError):
            return None
        else:
            return None
    names = set(re.findall(r'Product\s+"([^"\n]+)"', text, re.I))
    xyz = re.search(r"\*\s+xyzfile\s+[-+]?\d+\s+\d+\s+(\S+)", text, re.I)
    if not xyz or len(names) != 1:
        return None
    names.add(xyz.group(1).strip('"'))
    digest = hashlib.sha256(inp.read_bytes())
    dependencies = []
    for name in sorted(names):
        path = inp.parent / name
        if not path.is_file():
            return None
        digest.update(name.encode())
        digest.update(path.read_bytes())
        dependencies.append(_artifact(path))
    if contract.get("sha256") != digest.hexdigest():
        return None
    return {"policy_id": NEB_PREOPT_REPAIR_POLICY, "source_input": _artifact(inp),
            "source_output": _artifact(output), "source_job_contract": _artifact(contract_path),
            "source_process_lease": _artifact(lease_path), "dependencies": dependencies,
            "electronic_structure_started": False, "scientific_model_changed": False,
            "manual": NEB_MANUAL}


def prepare_neb_preopt_input_repair(input_path: Path, expected_text: str) -> Path:
    inp = Path(input_path)
    evidence = inspect_neb_preopt_parser_failure(inp)
    normalized, count = re.subn(r"(?im)^(\s*)PreOptEnds(\s+true\s*)$", r"\1PreOpt\2", inp.read_text())
    if evidence is None or count != 1 or normalized != expected_text:
        raise RuntimeError(f"Conflicting retained ORCA input cannot be repaired: {inp}")
    repaired = inp.with_name(f"{inp.stem}-syntax-repair-001.inp")
    record_path = repaired.with_suffix(".syntax-repair.json")
    record = {"kind": "orca_neb_input_syntax_repair", "schema_version": 1, **evidence,
              "repaired_input": {"path": str(repaired.resolve()),
                                 "sha256": hashlib.sha256(expected_text.encode()).hexdigest()},
              "correction": "PreOptEnds true -> PreOpt true; all other input bytes unchanged"}
    if repaired.exists() or record_path.exists():
        if (not repaired.is_file() or not record_path.is_file()
                or repaired.read_text() != expected_text or json.loads(record_path.read_text()) != record):
            raise RuntimeError(f"Retained ORCA syntax repair evidence differs: {repaired}")
        return repaired
    with repaired.open("x") as stream:
        stream.write(expected_text)
    with record_path.open("x") as stream:
        json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return repaired


def inspect_neb_restart_preopt_parser_failure(input_path: Path) -> dict | None:
    """Inspect the exact ORCA 6.1 restart/preoptimization conflict, read-only.

    Admission proves an execution-input defect, not a path or a rate. The
    retained native band, both endpoints and their executed input contract must
    still match; an active or indeterminate worker always blocks recovery.
    """
    inp = Path(input_path)
    output, contract_path, lease_path = (inp.with_suffix(suffix) for suffix in
                                         (".out", ".contract.json", ".process.json"))
    if not all(path.is_file() for path in (inp, output, contract_path, lease_path)):
        return None
    try:
        raw = inp.read_bytes()
        text = raw.decode("utf-8")  # Preserve line endings in the eventual repair.
        out = output.read_text(errors="replace")
        restart = re.findall(
            r'(?im)^[ \t]*Restart_ALLXYZFile[ \t]+"([^"\r\n]+)"[ \t]*\r?$', text,
        )
        product = re.findall(r'(?im)^[ \t]*Product[ \t]+"([^"\r\n]+)"[ \t]*\r?$', text)
        images = re.findall(r"(?im)^[ \t]*NImages[ \t]+(\d+)[ \t]*\r?$", text)
        xyz = re.findall(
            r'(?im)^[ \t]*\*[ \t]+xyzfile[ \t]+[-+]?\d+[ \t]+\d+[ \t]+'
            r'("[^"\r\n]+"|[^\s"]+)[ \t]*\r?$', text,
        )
        preopt = re.findall(r"(?im)^[ \t]*PreOpt\b[^\r\n]*", text)
        if ("NEB-TS" not in text.upper()
                or len(re.findall(r"(?im)^[ \t]*%neb[ \t]*\r?$", text)) != 1
                or len(restart) != 1 or len(product) != 1 or len(images) != 1 or len(xyz) != 1
                or int(images[0]) < 3 or len(preopt) != 1
                or not re.fullmatch(r"[ \t]*PreOpt[ \t]+true[ \t]*", preopt[0], re.I)
                or not re.search(r"Program Version\s+6\.1\.\d+\b", out)
                or not re.search(
                    r"(?m)^Preoptimization of end points requested\. "
                    r"This option is not possible with restart\.[ \t]*$", out,
                )
                or not re.search(
                    r"(?m)^\[file orca_main/main_driver_opt_neb\.cpp, line \d+\]: "
                    r"Combination not possible[ \t]*$", out,
                )
                or re.search(
                    r"ORCA TERMINATED NORMALLY|SCF CONVERGED|FINAL SINGLE POINT ENERGY|"
                    r"NEB ITERATION|THE NEB OPTIMIZATION HAS CONVERGED|"
                    r"(?m:^[ \t]*(?:Warning:.*?[ \t]+)?LBFGS[ \t]+\d+[ \t]+\d+)",
                    out, re.I,
                )):
            return None
        contract = json.loads(contract_path.read_text())
        lease = json.loads(lease_path.read_text())
        if not isinstance(contract, dict) or not isinstance(lease, dict):
            return None
        if (contract.get("kind") != "orca_job_input_contract"
                or Path(contract.get("input", "")).resolve() != inp.resolve()
                or lease.get("kind") != "storca_external_process_lease" or lease.get("engine") != "orca"
                or lease.get("status") != "terminated" or type(lease.get("returncode")) is not int
                or lease["returncode"] == 0
                or lease.get("input_sha256") != hashlib.sha256(raw).hexdigest()
                or Path(lease.get("input", "")).resolve() != inp.resolve()
                or Path(lease.get("output", "")).resolve() != output.resolve()):
            return None
        started, ended = float(lease["started_at_epoch"]), float(lease["ended_at_epoch"])
        if not all(math.isfinite(value) for value in (started, ended)) \
                or started <= 0 or ended < started or output.stat().st_mtime > ended + 1:
            return None
        for key, group in (("pid", False), ("process_group_id", True)):
            identifier = lease.get(key)
            if type(identifier) is not int or identifier <= 0:
                return None
            try:
                os.kill(-identifier if group else identifier, 0)
            except ProcessLookupError:
                pass
            else:
                return None
        names = set(re.findall(
            r'(?:Product|Hess_Filename|Restart_ALLXYZFile)\s+"([^"\r\n]+)"', text, re.I,
        ))
        names.add(xyz[0].strip('"'))
        digest = hashlib.sha256(raw)
        dependencies = []
        paths = []
        for name in sorted(names):
            path = inp.parent / name
            if Path(name).is_absolute() or not path.is_file() or path.parent.resolve() != inp.parent.resolve():
                return None
            digest.update(name.encode())
            digest.update(path.read_bytes())
            dependencies.append(_artifact(path))
            paths.append(str(path.resolve()))
        references = contract.get("referenced_files")
        if (not isinstance(references, list) or not all(isinstance(path, str) for path in references)
                or sorted(str(Path(path).resolve()) for path in references) != sorted(paths)
                or contract.get("sha256") != digest.hexdigest()):
            return None
        # Reuse the central finite/complete/ordered native-band validator,
        # rather than creating a second geometry acceptance policy here.
        from storca.reaction_path import _validate_orca_neb_restart_band
        from storca.route_geometry import read_xyz

        band = inp.parent / restart[0]
        validation = _validate_orca_neb_restart_band(band, expected_images=int(images[0]) + 2)
        reactant_elements, reactant_coordinates, _ = read_xyz(inp.parent / xyz[0].strip('"'))
        product_elements, product_coordinates, _ = read_xyz(inp.parent / product[0])
        if (not validation["valid"] or validation["elements"] != reactant_elements
                or product_elements != reactant_elements
                or not all(math.isfinite(float(value)) for coordinates in
                           (reactant_coordinates, product_coordinates) for value in coordinates.ravel())):
            return None
        return {
            "policy_id": NEB_RESTART_PREOPT_REPAIR_POLICY,
            "source_input": _artifact(inp), "source_output": _artifact(output),
            "source_job_contract": _artifact(contract_path), "source_process_lease": _artifact(lease_path),
            "source_band": _artifact(band), "dependencies": dependencies, "band_validation": validation,
            "electronic_structure_started": False, "scientific_model_changed": False,
            "rate_claim_allowed": False, "lifetime_claim_allowed": False,
            "renews_compute_budget": False, "manual": NEB_MANUAL,
            "reason": "ORCA 6.1 rejects endpoint preoptimization when loading a native restart band",
        }
    except (OSError, TypeError, ValueError, KeyError):
        return None


def prepare_neb_restart_preopt_input_repair(input_path: Path) -> Path:
    """Append one exact PreOpt true -> false correction; never overwrite."""
    inp = Path(input_path)
    evidence = inspect_neb_restart_preopt_parser_failure(inp)
    if evidence is None:
        raise RuntimeError(f"Retained ORCA restart/preopt conflict cannot be repaired: {inp}")
    original = inp.read_bytes()
    corrected, count = re.subn(
        rb"(?im)^([ \t]*PreOpt[ \t]+)true([ \t]*\r?$)", rb"\1false\2", original,
    )
    if count != 1 or hashlib.sha256(original).hexdigest() != evidence["source_input"]["sha256"]:
        raise RuntimeError(f"Retained ORCA restart input changed during repair: {inp}")
    repaired = inp.with_name(f"{inp.stem}-restart-preopt-repair-001.inp")
    record_path = repaired.with_suffix(".syntax-repair.json")
    record = {
        "kind": "orca_neb_restart_preopt_input_repair", "schema_version": 1, **evidence,
        "repaired_input": {"path": str(repaired.resolve()), "sha256": hashlib.sha256(corrected).hexdigest()},
        "correction": "PreOpt true -> PreOpt false for native band restart; all other input bytes unchanged",
    }
    if repaired.exists() or record_path.exists():
        try:
            matches = (repaired.is_file() and record_path.is_file() and repaired.read_bytes() == corrected
                       and json.loads(record_path.read_text()) == record)
        except (OSError, TypeError, ValueError):
            matches = False
        if not matches:
            raise RuntimeError(f"Retained ORCA restart/preopt repair evidence differs: {repaired}")
        return repaired
    with repaired.open("xb") as stream:
        stream.write(corrected)
    with record_path.open("x") as stream:
        json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return repaired
