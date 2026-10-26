"""Pure, bounded evidence for the MPI library observed at ORCA gradient workers.

This module neither launches processes nor reads files. The caller must bind
the native output, ORCA executable, resolved private MPI prefix and runtime
inventory, and check process ownership and completion independently. A loader
receipt does not verify a saddle, rate, complete runtime, or crash cause.
"""
from __future__ import annotations

import os
from pathlib import PurePosixPath
import re
from typing import NamedTuple


class DyldToken(NamedTuple):
    """An intact header and evidence-relevant image path, observed in text."""

    pid: str
    image_uuid: str
    path: str


_UUID = (r"[0-9A-Fa-f]{8}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-"
         r"[0-9A-Fa-f]{4}-[0-9A-Fa-f]{12}")
_TOKEN = re.compile(
    rf"dyld\[(\d+)\]:\s*<({_UUID})>\s+"
    r"(/[^\s<>]+?)(?=dyld\[|\s|$)"
)
_GRADIENT_IMAGE = "orca_scfgrad_mpi"
_MPI_IMAGE = "libmpi.40.dylib"
_RELEVANT_IMAGES = {_GRADIENT_IMAGE, _MPI_IMAGE}
_REFUSAL = "Private MPI loading at actual gradient workers not proved"


def parse_dyld_tokens(native: str) -> tuple[DyldToken, ...]:
    """Return intact tokens naming the gradient executable or ABI-40 MPI image.

    Parallel stderr can glue another ``dyld[`` header directly to an image path.
    Stop the regex *before* that header: trimming a greedy capture afterwards
    would discard its PID and prevent the next match from seeing the record.
    Do not reconstruct damaged headers, UUIDs or paths. Unrelated image paths
    are omitted because an interleaved framework fragment cannot establish its
    full endpoint; a relevant image must have its complete audited basename.
    PIDs are canonicalized so repeated/zero-padded IDs cannot inflate coverage.
    """
    if not isinstance(native, str):
        raise TypeError("Native loader output must be text")
    tokens = []
    for match in _TOKEN.finditer(native):
        pid, image_uuid, path = match.groups()
        if int(pid) <= 0 or PurePosixPath(path).name not in _RELEVANT_IMAGES:
            continue
        tokens.append(DyldToken(str(int(pid)), image_uuid, path))
    return tuple(tokens)


def _absolute_lexical_path(value: str | os.PathLike[str]) -> str:
    """Normalize redundant slashes, refusing ambiguous parent/dot aliases."""
    raw = os.fspath(value)
    if (not isinstance(raw, str) or not raw.startswith("/") or raw.startswith("//")
            or any(part in (".", "..") for part in raw.split("/"))
            or any(ord(char) < 32 for char in raw)):
        raise ValueError("Unambiguous absolute private MPI image path required")
    return str(PurePosixPath(raw))


def private_mpi_runtime_witness(
    native: str, mpi_prefix: str | os.PathLike[str]
) -> dict:
    """Require >=4 observed gradient PIDs, each bound to the private ABI-40 image.

    Every complete observed ABI-40 library token, including those at other
    processes, must identify the expected private image. This retains the
    previous witness's minimum-worker and mixed-runtime refusal gates. Unlike
    that filesystem-resolving reader, this pure reader compares absolute paths
    lexically and does not admit symlink aliases or parent/dot paths. The caller
    must supply the resolved, independently sealed prefix and reopen its files.

    Missing or garbled required evidence remains unresolved. Success covers
    identified workers only, not unobserved ranks or every loaded library.
    """
    prefix = _absolute_lexical_path(mpi_prefix)
    expected = str(PurePosixPath(prefix) / "lib" / _MPI_IMAGE)
    tokens = parse_dyld_tokens(native)
    workers = {token.pid for token in tokens
               if PurePosixPath(token.path).name == _GRADIENT_IMAGE}
    libraries = [token for token in tokens
                 if PurePosixPath(token.path).name == _MPI_IMAGE]
    witnessed = {token.pid for token in libraries}
    if len(workers) < 4 or not libraries or not workers.issubset(witnessed):
        raise ValueError(_REFUSAL)
    try:
        paths_match = all(_absolute_lexical_path(token.path) == expected
                          for token in libraries)
    except (TypeError, ValueError) as error:
        raise ValueError(_REFUSAL) from error
    if not paths_match:
        raise ValueError(_REFUSAL)
    return {
        "gradient_worker_PIDs": sorted(workers, key=int),
        "MPI_loading_witnesses": len(libraries),
        "loaded_MPI_library": expected,
        "all_identified_gradient_workers_witnessed": True,
        "path_comparison": "absolute lexical; caller binds resolved runtime files",
        "scope": ("Complete native loader tokens at identified gradient workers only; "
                  "not a proof of all unobserved workers, every loaded library, "
                  "or original crash causation"),
    }
