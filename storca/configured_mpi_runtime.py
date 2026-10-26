"""Pure receipt for explicit MPI selection; loader traces are diagnostic only.

Selecting an environment does not prove which library a native worker loaded.
The caller must independently bind its runtime inventory, process ownership,
native completion, and scientific checks. This module does not inspect files,
launch processes, reconstruct parallel loader output, or admit chemical claims.
"""
from __future__ import annotations

import os
from pathlib import Path, PurePosixPath
import re
from typing import Callable


_ENVIRONMENT_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")
_LIBRARY_VARIABLES = ("DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH")
_RELAY_VARIABLE = "OMPI_MCA_mca_base_env_list"


def _prefix_path(prefix: str | os.PathLike[str]) -> str:
    """Compare configured paths lexically without filesystem resolution."""
    raw = os.fspath(prefix)
    if (not isinstance(raw, str) or not raw.startswith("/")
            or raw.startswith("//")
            or any(part in (".", "..") for part in raw.split("/"))
            or any(ord(char) < 32 for char in raw)
            or any(char in raw for char in (":", ";"))):
        raise ValueError("Unambiguous absolute configured MPI prefix required")
    return str(PurePosixPath(raw))


def _selection_environment(environment: dict, prefix: str) -> None:
    if not isinstance(environment, dict):
        raise TypeError("Selected MPI environment must be a dictionary")
    for name, value in environment.items():
        if (not isinstance(name, str) or not _ENVIRONMENT_NAME.fullmatch(name)
                or not isinstance(value, str)
                or any(ord(char) < 32 for char in value)):
            raise ValueError("Selected MPI environment contains a malformed entry")

    expected = {"PATH": str(PurePosixPath(prefix) / "bin")}
    expected.update({name: str(PurePosixPath(prefix) / "lib")
                     for name in _LIBRARY_VARIABLES})
    for name, first_path in expected.items():
        value = environment.get(name)
        if value is None or value.split(":", 1)[0] != first_path:
            raise ValueError(f"Configured MPI {name} must start with {first_path}")

    # Open MPI's relay list is optional here, not evidence of rank execution.
    # When supplied, its exact assignments cannot contradict the selected
    # environment or omit the explicitly selected library paths.
    if _RELAY_VARIABLE not in environment:
        return
    relay = {}
    for assignment in environment[_RELAY_VARIABLE].split(";"):
        name, separator, value = assignment.partition("=")
        if (not separator or not _ENVIRONMENT_NAME.fullmatch(name)
                or name in relay):
            raise ValueError("Configured MPI worker environment relay is malformed or duplicate")
        if name not in environment or environment[name] != value:
            raise ValueError("Configured MPI worker environment relay contradicts selection")
        relay[name] = value
    required = set(_LIBRARY_VARIABLES)
    if "DYLD_PRINT_LIBRARIES" in environment:
        required.add("DYLD_PRINT_LIBRARIES")
    if not required.issubset(relay):
        raise ValueError("Configured MPI worker environment relay omits selected values")


def configured_mpi_receipt(
    native: str,
    prefix: str | Path,
    selected_environment: dict,
    *,
    trace_reader: Callable[[str, Path], dict],
) -> dict:
    """Check explicit selection and retain a non-blocking loader diagnostic.

    Configuration inconsistencies refuse before reading diagnostic text.
    A historical reader's ``ValueError`` means its trace evidence is unresolved,
    not that an otherwise owned, normally completed calculation failed. Other
    exceptions propagate so implementation defects are not silently hidden.
    No loader text, successful or otherwise, replaces the caller's completion,
    inventory, geometry, gradient, Hessian, or reaction-validation gates.
    """
    if not isinstance(native, str):
        raise TypeError("Native loader output must be text")
    selected_prefix = _prefix_path(prefix)
    _selection_environment(selected_environment, selected_prefix)
    try:
        evidence = trace_reader(native, Path(selected_prefix))
    except ValueError as error:
        loader_trace = {"status": "unresolved_diagnostic", "reason": str(error)}
    else:
        loader_trace = {"status": "observed", "evidence": evidence}
    return {
        "kind": "configured_mpi_runtime_diagnostic",
        "schema_version": 1,
        "selected_prefix": selected_prefix,
        "configured_launcher": str(PurePosixPath(selected_prefix) / "bin/mpirun"),
        "configured_library": str(PurePosixPath(selected_prefix) / "lib/libmpi.40.dylib"),
        "selection_consistent": True,
        "loader_trace_required_for_acceptance": False,
        "loaded_runtime_proved_by_configuration": False,
        "loader_trace": loader_trace,
        "scope": ("Explicit environment selection only; caller separately binds runtime "
                  "inventory, process ownership, native completion, and physical checks. "
                  "Configuration does not prove actual loaded runtime or scientific validity."),
    }
