"""Bound, explicit runtime selection; never evidence of loaded native workers.

Only the declared MPI inventory and selected ORCA helpers are reopened. No
binary is executed, no loader tracing is enabled, and no calculation authorized.
"""
from __future__ import annotations

from contextlib import contextmanager
import os
from pathlib import Path

from . import rate_transfer_preparation as prep
from .configured_mpi_runtime import _prefix_path, _selection_environment


HELPERS = ("orca", "orca_scf", "orca_leanscf", "orca_leanscf_mpi",
           "orca_scfgrad_mpi", "orca_mdci", "orca_mdci_mpi", "orca_numfreq", "orca_neb")
REQUIRED = {"orca", "orca_scfgrad_mpi", "orca_mdci_mpi"}
KIND = "storca_rate_transfer_runtime_selection"
FLAGS = {"native_jobs_launched": 0, "authorized": False,
         "loaded_runtime_proved_by_configuration": False,
         "loader_trace_required_for_acceptance": False}


def _environment(prefix, binary):
    environment = {
        "PATH": ":".join((str(prefix / "bin"), str(binary.parent),
                          "/opt/homebrew/anaconda3/condabin", "/opt/homebrew/bin",
                          "/usr/local/bin", "/usr/bin", "/bin", "/usr/sbin", "/sbin")),
        "DYLD_LIBRARY_PATH": str(prefix / "lib"),
        "DYLD_FALLBACK_LIBRARY_PATH": str(prefix / "lib"),
        **{name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                  "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
    }
    environment["OMPI_MCA_mca_base_env_list"] = ";".join(
        name + "=" + value for name, value in environment.items())
    _selection_environment(environment, str(prefix))
    return environment


def _source(source_binding):
    source = prep._json(prep._reopen(source_binding))
    binary = prep._reopen(source.get("orca"))
    runtime = source.get("runtime")
    if not isinstance(runtime, dict) or not {"prefix", "inventory"}.issubset(runtime):
        raise ValueError("Source contract lacks a declared MPI runtime inventory")
    prefix = Path(_prefix_path(runtime["prefix"]))
    inventory = runtime["inventory"]
    if not isinstance(inventory, dict) or not inventory:
        raise ValueError("Nonempty source MPI inventory required")
    for key, record in inventory.items():
        if (not isinstance(record, dict) or set(record) != {"path", "sha256", "mtime_ns"}
                or key != record["path"] or type(record["mtime_ns"]) is not int
                or record["mtime_ns"] < 0):
            raise ValueError("Malformed source MPI inventory record")
        path = prep._reopen({name: record[name] for name in ("path", "sha256")})
        if (not path.is_relative_to(prefix) or path.stat().st_mtime_ns != record["mtime_ns"]):
            raise ValueError("Source MPI inventory path or modification time changed")
    # Resolve declared symlinks to a bound inventory member, not an unbound PATH result.
    for path in (prefix / "bin/mpirun", prefix / "lib/libmpi.40.dylib"):
        if str(path.resolve(strict=True)) not in inventory:
            raise ValueError("Selected MPI launcher/library is absent from bound inventory")
    if binary.name != "orca":
        raise ValueError("Source ORCA binary must be the declared orca driver")
    helpers = {name: prep.binding(binary.parent / name) for name in HELPERS
               if (binary.parent / name).is_file()}
    if not REQUIRED.issubset(helpers):
        raise ValueError("Selected ORCA helper inventory lacks required parallel workers")
    return {"schema_version": 1, "kind": KIND, "source_contract": source_binding,
            "orca_binary": source["orca"], "prefix": str(prefix),
            "mpi_inventory": inventory, "orca_helpers": helpers,
            "selected_environment": _environment(prefix, binary), **FLAGS,
            "scope": "Bound runtime selection only; ownership, native completion and physical gates remain separate."}


def prepare_runtime(source_contract: Path, destination: Path) -> dict:
    """Publish an immutable selection from retained hashes, with zero native work."""
    value = _source(prep.binding(source_contract))
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    prep._publish(destination, value)
    return prep.binding(destination)


def reopen_runtime(binding, *, orca_binary, extra_environment):
    """Reopen the exact selection, source hashes, helper set and selected environment."""
    value = prep._json(prep._reopen(binding))
    expected = _source(value.get("source_contract"))
    if (prep.canonical_hash(value) != prep.canonical_hash(expected)
            or value["orca_binary"] != orca_binary
            or value["selected_environment"] != extra_environment):
        raise ValueError("Runtime selection schema, inventory or execution environment differs")
    return value


@contextmanager
def runtime_environment(binding):
    """Select only the bound driver; the owned dispatch receives other variables."""
    selected = prep._json(prep._reopen(binding))
    value = reopen_runtime(binding, orca_binary=selected.get("orca_binary"),
                           extra_environment=selected.get("selected_environment"))
    binary = value["orca_binary"]["path"]
    prior = os.environ.get("STORCA_ORCA_BIN")
    if prior is not None and prior != binary:
        raise ValueError("Explicit STORCA_ORCA_BIN contradicts the bound runtime selection")
    try:
        if prior is None:
            os.environ["STORCA_ORCA_BIN"] = binary
        yield value
    finally:
        if prior is None:
            os.environ.pop("STORCA_ORCA_BIN", None)
        else:
            os.environ["STORCA_ORCA_BIN"] = prior
