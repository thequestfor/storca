"""Opt-in, bounded native dispatch for the known-reaction transfer benchmark.

The immutable journal counts attempts, including failed executions. It sums
owned execution time rather than time between workflow invocations. It does
not grant authorization, classify chemistry, or validate a physical rate.
"""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import stat
import threading
import time

from src.orca_job_settings import NativeJobSettings, explicit_native_input_scope, native_settings_scope


CATEGORIES = ("optimization_path", "frequency", "irc", "single_point")
_ACTIVE = ContextVar("storca_explicit_rate_transfer_control", default=None)


class TransferControlError(RuntimeError):
    """Budget, ownership, source integrity, or authorization is unresolved."""


def _bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode()


def _hash(data):
    return hashlib.sha256(data).hexdigest()


def _read(path):
    if path.is_symlink() or not path.is_file():
        raise TransferControlError(f"Missing regular evidence file: {path}")
    data = path.read_bytes()
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise TransferControlError("Duplicate journal JSON key")
            result[key] = value
        return result
    try:
        value = json.loads(data, object_pairs_hook=pairs)
        if _bytes(value) != data:
            raise TransferControlError("Noncanonical or incomplete journal record")
        return value
    except (ValueError, UnicodeError) as error:
        raise TransferControlError("Invalid retained control record") from error


def _write_new(path, value):
    with path.open("xb") as stream:
        stream.write(_bytes(value))
        stream.flush()
        os.fsync(stream.fileno())
    descriptor = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _positive(value, field):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise TransferControlError(f"{field} must be finite and positive")


def _contract(value):
    # Round-trip first so caller mutation cannot change active policy.
    value = json.loads(_bytes(value))
    expected = {"schema_version", "kind", "campaign_id", "authorized", "resources", "limits", "extra_environment"}
    if set(value) != expected or value["schema_version"] != 1 or value["kind"] != "storca_rate_transfer_execution_control":
        raise TransferControlError("Unsupported transfer execution control contract")
    if not isinstance(value["campaign_id"], str) or not value["campaign_id"] or type(value["authorized"]) is not bool:
        raise TransferControlError("Explicit campaign identity and authorization are required")
    resources = value["resources"]
    if not isinstance(resources, dict) or set(resources) != {"ncores", "maxcore_mb_per_rank", "reserve_disk_bytes", "reserve_memory_bytes"}:
        raise TransferControlError("Complete native resource safeguards are required")
    for field, item in resources.items():
        if type(item) is not int or item < (1 if field in {"ncores", "maxcore_mb_per_rank"} else 0):
            raise TransferControlError(f"Invalid resource {field}")
    limits = value["limits"]
    if not isinstance(limits, dict) or set(limits) != {*CATEGORIES, "per_job_seconds", "total_wall_seconds", "total_core_seconds"}:
        raise TransferControlError("Complete categorized compute limits are required")
    for category in CATEGORIES:
        if type(limits[category]) is not int or limits[category] < 0:
            raise TransferControlError("Attempt caps must be nonnegative integers")
    for field in ("per_job_seconds", "total_wall_seconds", "total_core_seconds"):
        _positive(limits[field], field)
    if not isinstance(value["extra_environment"], dict) or any(
        not isinstance(key, str) or not isinstance(item, str)
        for key, item in value["extra_environment"].items()
    ):
        raise TransferControlError("Native environment must be explicit string pairs")
    return value


def classify_native_input(text):
    """Count one engine dispatch, not its images or displaced workers."""
    words = re.findall(r"(?im)^\s*!([^\n]*)", text)
    tokens = {item.casefold() for line in words for item in line.split()}
    if not words:
        raise TransferControlError("Cannot classify a native job without keywords")
    groups = []
    if tokens & {"freq", "anfreq", "numfreq"}:
        groups.append("frequency")
    if tokens & {"opt", "optts", "neb-ts", "neb", "neb-ci", "neb-ci-ts", "copt"}:
        groups.append("optimization_path")
    if "irc" in tokens:
        groups.append("irc")
    if len(groups) > 1:
        raise TransferControlError("Combined native job categories must be split explicitly")
    if groups:
        return groups[0]
    # These jobs have a different computational meaning and cannot masquerade
    # as an inexpensive SP under this narrowly scoped contract.
    if tokens & {"md", "engrad", "goat", "numgrad", "compound", "eprnmr", "scan"}:
        raise TransferControlError("Unsupported native dispatch category")
    if any(re.search(r"(?im)^\s*%" + block + r"\b", text) for block in ("md", "compound", "goat")):
        raise TransferControlError("Unsupported native dispatch block")
    if not tokens & {"wb97x-d3", "b3lyp", "hf", "rhf", "uhf", "ccsd(t)", "ccsd", "mp2"}:
        raise TransferControlError("Single-point input lacks a supported explicit electronic method")
    return "single_point"


def validate_control_contract(value):
    """Validate and copy a proposed contract without authorizing execution."""
    return _contract(value)


def active_execution_control():
    return _ACTIVE.get()


def _is_h2_ccsd_input(text):
    """Recognize only the inline, neutral singlet H2 canonical-CC input."""
    tokens = {t.casefold() for line in re.findall(r"(?im)^\s*!([^\n]*)", text) for t in line.split()}
    match = re.search(r"(?im)^\s*\*\s+xyz\s+0\s+1\s*\n(.*?)^\s*\*\s*$", text, re.DOTALL)
    rows = [line.split() for line in match.group(1).splitlines() if line.strip()] if match else []
    return (bool({"rhf", "ccsd(t)"}.issubset(tokens)) and len(rows) == 2
            and all(len(row) == 4 and row[0] == "H" for row in rows))


def _serial_h2_retry_text(text):
    if not _is_h2_ccsd_input(text) or text.count("%pal\n  nprocs 4\nend\n") != 1:
        raise ValueError("Serial exception requires the exact four-rank H2 CCSD(T) input")
    return text.replace("%pal\n  nprocs 4\nend\n", "%pal\n  nprocs 1\nend\n", 1)


class NativeExecutionControl:
    def __init__(self, root, contract, *, execution_amendment=None):
        self.root = Path(root).resolve()
        self._policy = _contract(contract)
        self.contract_sha256 = _hash(_bytes(self._policy))
        self._lock = None
        self._owner = None
        self._execution_amendment = json.loads(_bytes(execution_amendment))

    def validate_execution_amendment(self, binding=None):
        """One explicitly authorized PMIx workaround, on the original ledger.

        This does not rewrite the original contract or refund a failed job.
        Authority is limited to frozen single-point inputs and one extra SP.
        Historical readers validate the same authority from each reservation.
        """
        from . import rate_transfer_preparation as prep
        binding = self._execution_amendment if binding is None else binding
        try:
            value = prep._json(prep._reopen(binding))
            fields = {"schema_version", "kind", "authorized", "original_control_contract",
                      "original_single_point_cap", "revised_single_point_cap", "journal_prefix",
                      "prefix_spending", "extra_environment", "runtime_source", "allowed_inputs",
                      "failed_input", "failed_output"}
            serial = value.get("kind") == "storca_h2_serial_single_point_amendment"
            memory = value.get("kind") == "storca_ts_memory_single_point_amendment"
            cap = self.contract["limits"]["single_point"]
            if serial or memory:
                fields |= {"previous_amendment", "serial_input"}
                if memory:
                    fields |= {"memory_input", "maxcore_mb_per_rank"}
                parent = prep._json(prep._reopen(value["previous_amendment"]))
                parent_kind = "storca_h2_serial_single_point_amendment" if memory else "storca_pmix_hash_single_point_amendment"
                if parent.get("kind") != parent_kind:
                    raise ValueError("Exact predecessor execution amendment required")
                parent = self.validate_execution_amendment(value["previous_amendment"])
                cap = parent["revised_single_point_cap"]
            if (set(value) != fields or value["schema_version"] != 1
                    or value["kind"] != ("storca_ts_memory_single_point_amendment" if memory else
                                         "storca_h2_serial_single_point_amendment" if serial else "storca_pmix_hash_single_point_amendment")
                    or value["authorized"] is not True
                    or type(value["original_single_point_cap"]) is not int
                    or type(value["revised_single_point_cap"]) is not int
                    or value["original_single_point_cap"] != cap
                    or value["revised_single_point_cap"] != cap + 1
                    or value["extra_environment"] != {"PMIX_MCA_gds": "hash"}
                    or value["original_control_contract"] != prep.binding(self.root / "contract.json")):
                raise ValueError("Amendment changes more than the authorized PMIx/SP scope")
            plugin = Path(self.contract["extra_environment"]["DYLD_LIBRARY_PATH"]) / "pmix/mca_gds_hash.so"
            if prep._reopen(value["runtime_source"]) != plugin:
                raise ValueError("PMIx hash component is not from the selected runtime")
            events = self._events()
            prefix = value["journal_prefix"]
            if not prefix or len(prefix) > len(events):
                raise ValueError("Missing retained amendment checkpoint")
            paths = sorted((self.root / "events").iterdir())[:len(prefix)]
            if [prep.binding(path) for path in paths] != prefix:
                raise ValueError("Amendment journal prefix differs")
            reservations = {e["attempt_id"]: e for e in events[:len(prefix)] if e["event"] == "reserved"}
            terminals = {e["attempt_id"]: e for e in events[:len(prefix)] if e["event"] == "terminal"}
            spending = {"counts": {c: sum(r["category"] == c for r in reservations.values()) for c in CATEGORIES},
                        "wall_seconds": sum(t["elapsed_seconds"] for t in terminals.values()),
                        "core_seconds": sum(t["elapsed_seconds"] * reservations[k]["ncores"] for k, t in terminals.items())}
            if (len(reservations) != len(terminals) or set(reservations) != set(terminals)
                    or len(prefix) != 2 * len(reservations) or spending != value["prefix_spending"]):
                raise ValueError("Amendment requires exact closed spending, without refunds")
            inputs = value["allowed_inputs"]
            if len(inputs) != 5 or len({row["path"] for row in inputs}) != 5:
                raise ValueError("Exactly five distinct frozen energy inputs required")
            for row in inputs:
                path = prep._reopen(row)
                if classify_native_input(path.read_text()) != "single_point":
                    raise ValueError("Only single-point inputs may use this amendment")
            failed = prep._reopen(value["failed_input"])
            failed_output = prep._reopen(value["failed_output"])
            matches = [r for r in reservations.values() if r["input"] == str(failed)]
            retry_matches = any(row["sha256"] == value["failed_input"]["sha256"]
                                and row["path"] != str(failed) for row in inputs)
            if serial:
                retry = prep._reopen(value["serial_input"])
                prior_inputs = parent["allowed_inputs"]
                if (value["serial_input"] not in inputs or retry == failed
                        or retry.read_text() != _serial_h2_retry_text(failed.read_text())
                        or value["failed_input"] not in prior_inputs
                        or len(prefix) <= len(parent["journal_prefix"])
                        or [r for r in inputs if r != value["serial_input"]]
                           != [r for r in prior_inputs if r != value["failed_input"]]
                        or "Number of processes (4) in parallel calculation exceeds number of pairs (1)" not in failed_output.read_text()):
                    raise ValueError("Serial H2 retry changed chemistry, unrelated inputs or failure evidence")
                retry_matches = True
            if memory:
                retry = prep._reopen(value["memory_input"])
                text = failed.read_text()
                if (type(value["maxcore_mb_per_rank"]) is not int or value["maxcore_mb_per_rank"] != 1500
                        or value["serial_input"] != parent["serial_input"]
                        or value["memory_input"] not in inputs or retry == failed
                        or value["failed_input"] not in parent["allowed_inputs"]
                        or text.count("%maxcore 1000\n") != 1
                        or retry.read_text() != text.replace("%maxcore 1000\n", "%maxcore 1500\n", 1)
                        or "%pal\n  nprocs 4\nend\n" not in text
                        or len(prefix) <= len(parent["journal_prefix"])
                        or [r for r in inputs if r != value["memory_input"]]
                           != [r for r in parent["allowed_inputs"] if r != value["failed_input"]]
                        or not all(marker in failed_output.read_text() for marker in
                                   ("ERROR - OUT OF MEMORY", "MINIMUM REQUIRED:  1359.7 MB", "MAXCORE:           1000.0 MB"))):
                    raise ValueError("TS memory retry changed chemistry, unrelated inputs or failure evidence")
                retry_matches = True
            if (len(matches) != 1 or failed_output != failed.with_suffix(".out")
                    or matches[0]["category"] != "single_point"
                    or terminals[matches[0]["attempt_id"]]["outcome"] != "failed"
                    or value["failed_output"] not in terminals[matches[0]["attempt_id"]]["sources"]
                    or not retry_matches):
                raise ValueError("Exact charged failed input/output and identical fresh retry required")
            return value
        except (ValueError, KeyError, TypeError, OSError) as error:
            raise TransferControlError("Invalid execution amendment: " + str(error)) from error

    @property
    def contract(self):
        return json.loads(_bytes(self._policy))

    def initialize(self):
        self.root.mkdir(parents=True, exist_ok=True)
        with self._exclusive():
            identity = self.root / "contract.json"
            events = self.root / "events"
            if identity.exists():
                if _read(identity) != self.contract:
                    raise TransferControlError("Immutable transfer control contract changed")
                if events.is_symlink() or not events.is_dir():
                    raise TransferControlError("Retained event directory missing")
            else:
                if events.exists() and (not events.is_dir() or any(events.iterdir())):
                    raise TransferControlError("Missing contract cannot reset retained attempts")
                events.mkdir(exist_ok=True)
                _write_new(identity, self.contract)
            self._events()
        return {"path": str(self.root / "contract.json"), "sha256": self.contract_sha256}

    @contextmanager
    def _exclusive(self):
        if self._lock is not None:
            raise TransferControlError("Transfer controller already active")
        descriptor = os.open(self.root / "controller.lock", os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0), 0o600)
        try:
            if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                raise TransferControlError("Native controller lock is not a regular file")
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise TransferControlError("Another transfer controller owns this campaign") from error
            self._lock = descriptor
            self._owner = (os.getpid(), threading.get_ident())
            yield
        finally:
            self._lock = None
            self._owner = None
            os.close(descriptor)

    def _events(self):
        if _read(self.root / "contract.json") != self.contract:
            raise TransferControlError("Retained native contract mutated")
        paths = sorted((self.root / "events").iterdir())
        previous = self.contract_sha256
        events = []
        for number, path in enumerate(paths, 1):
            event = _read(path)
            digest = _hash(_bytes(event))
            if path.name != f"{number:06d}-{digest}.json" or event.get("previous_sha256") != previous:
                raise TransferControlError("Transfer attempt chain has a mutation or gap")
            if event.get("event") not in {"reserved", "terminal"}:
                raise TransferControlError("Unknown transfer journal event")
            events.append(event)
            previous = digest
        return events

    def _append(self, fields):
        if self._lock is None or self._owner != (os.getpid(), threading.get_ident()):
            raise TransferControlError("Exclusive controller ownership required")
        events = self._events()
        previous = _hash(_bytes(events[-1])) if events else self.contract_sha256
        event = {**fields, "previous_sha256": previous}
        digest = _hash(_bytes(event))
        _write_new(self.root / "events" / f"{len(events)+1:06d}-{digest}.json", event)
        return event

    def snapshot(self):
        events = self._events()
        reservations = {}
        amendments = {}
        counts = {c: 0 for c in CATEGORIES}
        for event in events:
            if event["event"] == "reserved":
                attempt = event.get("attempt_id")
                if type(attempt) is not int or attempt != len(reservations) + 1:
                    raise TransferControlError("Reservation attempt sequence changed")
                if event.get("category") not in CATEGORIES or type(event.get("ncores")) is not int:
                    raise TransferControlError("Invalid reserved native category or resources")
                category = event["category"]
                counts[category] += 1
                cap = self.contract["limits"][category]
                binding = event.get("execution_amendment")
                expected_cores = self.contract["resources"]["ncores"]
                if binding is not None:
                    key = _bytes(binding)
                    if key not in amendments:
                        amendments[key] = self.validate_execution_amendment(binding)
                    amended = amendments[key]
                    if (category != "single_point"
                            or {"path": event["input"], "sha256": event["native_input_file_sha256"]} not in amended["allowed_inputs"]
                            or (amended.get("memory_input") is not None
                                and amended["memory_input"] != {"path": event["input"], "sha256": event["native_input_file_sha256"]})
                            or attempt <= len(amended["journal_prefix"]) // 2):
                        raise TransferControlError("Reservation exceeds amendment scope")
                    cap = amended["revised_single_point_cap"]
                    if amended.get("serial_input") == {"path": event["input"], "sha256": event["native_input_file_sha256"]}:
                        expected_cores = 1
                if event["ncores"] != expected_cores:
                    raise TransferControlError("Reserved core count exceeds its exact input authority")
                if counts[category] > cap:
                    raise TransferControlError("Journal exceeds authorized native attempt cap")
                reservations[attempt] = event
        terminals = {}
        for event in events:
            if event["event"] == "terminal":
                attempt = event["attempt_id"]
                if attempt not in reservations or attempt in terminals:
                    raise TransferControlError("Duplicate or orphan native terminal receipt")
                if not math.isfinite(event["elapsed_seconds"]) or event["elapsed_seconds"] < 0:
                    raise TransferControlError("Invalid native elapsed time")
                terminals[attempt] = event
        return {"counts": {category: sum(event["category"] == category for event in reservations.values()) for category in CATEGORIES},
                "wall_seconds": sum(event["elapsed_seconds"] for event in terminals.values()),
                "core_seconds": sum(event["elapsed_seconds"] * reservations[key]["ncores"] for key, event in terminals.items()),
                "pending": [event for key, event in reservations.items() if key not in terminals],
                "reservations": reservations, "terminals": terminals}

    @contextmanager
    def activate(self, native_settings):
        if not isinstance(native_settings, NativeJobSettings):
            raise TypeError("Explicit NativeJobSettings required")
        if native_settings.maxcore_mb_per_rank != self.contract["resources"]["maxcore_mb_per_rank"]:
            raise TransferControlError("Native settings and memory contract disagree")
        if _ACTIVE.get() is not None:
            raise TransferControlError("Nested transfer controllers are not supported")
        self.initialize()
        with self._exclusive(), native_settings_scope(native_settings):
            settings_path = self.root / "native-settings.json"
            settings_identity = {"sha256": native_settings.sha256, "settings": native_settings.as_dict()}
            if settings_path.exists():
                if _read(settings_path) != settings_identity:
                    raise TransferControlError("Frozen native settings changed between continuations")
            else:
                _write_new(settings_path, settings_identity)
            if self.snapshot()["pending"]:
                raise TransferControlError("An unclosed native reservation remains; retained execution must be reconciled explicitly")
            # Some path-level caches return without calling the dispatcher.
            # Reopen all previously owned sources before allowing those caches.
            from .reaction_path import _orca_job_input_evidence
            for reservation in self.snapshot()["reservations"].values():
                path = Path(reservation["input"])
                try:
                    digest = _orca_job_input_evidence(path)[0].hexdigest()
                except OSError as error:
                    raise TransferControlError("Retained native input evidence is missing") from error
                self.validate_reuse(path, digest)
            token = _ACTIVE.set(self)
            try:
                yield self
            finally:
                _ACTIVE.reset(token)

    def validate_reuse(self, input_path, input_sha256):
        """Never silently adopt unknown native work or changed paid evidence."""
        path = Path(input_path).resolve()
        snapshot = self.snapshot()
        matches = [item for item in snapshot["reservations"].values() if item["input"] == str(path)]
        if not matches:
            if path.with_suffix(".out").exists():
                raise TransferControlError("Retained native output has no owned terminal receipt")
            return
        for reservation in matches:
            if reservation["input_sha256"] != input_sha256:
                raise TransferControlError("Reserved native input changed")
            terminal = snapshot["terminals"].get(reservation["attempt_id"])
            if terminal is None:
                raise TransferControlError("Native attempt ownership is unresolved")
            for binding in terminal["sources"]:
                retained = Path(binding["path"])
                if not retained.is_file() or _hash(retained.read_bytes()) != binding["sha256"]:
                    raise TransferControlError("Owned native output or process receipt changed")

    def _memory_available(self):
        # psutil is optional in production; no estimate is invented if absent.
        try:
            import psutil
            return int(psutil.virtual_memory().available)
        except ImportError as error:
            raise TransferControlError("Native memory safeguard unavailable") from error

    def reserve(self, input_path, input_sha256, *, timeout_seconds, ncores):
        if self._lock is None or self._owner != (os.getpid(), threading.get_ident()) or _ACTIVE.get() is not self:
            raise TransferControlError("Explicit active controller required")
        if not self.contract["authorized"]:
            raise TransferControlError("Native transfer execution has not been authorized")
        snapshot = self.snapshot()
        if snapshot["pending"]:
            raise TransferControlError("Previous native reservation is unresolved")
        resource, limits = self.contract["resources"], self.contract["limits"]
        input_text = Path(input_path).read_text()
        maxcore = re.findall(r"(?im)^\s*%maxcore\s+(\d+)\s*$", input_text)
        category = classify_native_input(input_text)
        amendment = None
        if self._execution_amendment is not None:
            amendment = self.validate_execution_amendment()
            from .rate_transfer_preparation import binding
            if category != "single_point" or binding(input_path) not in amendment["allowed_inputs"]:
                raise TransferControlError("Native input is outside the energy-only amendment")
            if amendment.get("memory_input") is not None and binding(input_path) != amendment["memory_input"]:
                raise TransferControlError("Only the exact TS memory retry may launch under this amendment")
            limits["single_point"] = amendment["revised_single_point_cap"]
        expected_cores = (1 if amendment and amendment.get("serial_input") == binding(input_path)
                          else resource["ncores"])
        expected_memory = (amendment["maxcore_mb_per_rank"]
                           if amendment and amendment.get("memory_input") == binding(input_path)
                           else resource["maxcore_mb_per_rank"])
        if maxcore != [str(expected_memory)]:
            raise TransferControlError("Native MaxCore differs from the frozen contract or exact TS exception")
        if type(ncores) is not int or ncores != expected_cores:
            raise TransferControlError("Native core count differs from the frozen contract or exact H2 exception")
        if any(item["input"] == str(Path(input_path).resolve()) for item in snapshot["reservations"].values()):
            raise TransferControlError("An attempted native basename cannot be overwritten; use an append-only retry")
        if snapshot["counts"][category] >= limits[category]:
            raise TransferControlError(f"Native {category} attempt cap reached")
        if shutil.disk_usage(input_path.parent).free < resource["reserve_disk_bytes"]:
            raise TransferControlError("Native disk reserve is unavailable")
        # Retain the conservative original campaign guard even for serial H2.
        required_memory = resource["reserve_memory_bytes"] + resource["ncores"] * expected_memory * 1024 * 1024
        if self._memory_available() < required_memory:
            raise TransferControlError("Native memory reserve is unavailable")
        timeout = min(limits["per_job_seconds"], limits["total_wall_seconds"] - snapshot["wall_seconds"],
                      (limits["total_core_seconds"] - snapshot["core_seconds"]) / ncores)
        if timeout_seconds is not None:
            _positive(timeout_seconds, "requested timeout")
            timeout = min(timeout, timeout_seconds)
        _positive(timeout, "remaining native timeout")
        record = self._append({"event": "reserved", "attempt_id": len(snapshot["reservations"]) + 1,
                               "category": category, "input": str(Path(input_path).resolve()),
                               "input_sha256": input_sha256, "ncores": ncores,
                               "native_input_file_sha256": _hash(Path(input_path).read_bytes()),
                               "timeout_seconds": timeout, "started_at_epoch": time.time(),
                               **({"execution_amendment": self._execution_amendment} if amendment else {})})
        return record

    def finish(self, reservation, *, elapsed_seconds, outcome):
        sources = []
        path = Path(reservation["input"])
        if not path.is_file() or _hash(path.read_bytes()) != reservation["native_input_file_sha256"]:
            raise TransferControlError("Native input changed during owned execution")
        lease_path = path.with_suffix(".process.json")
        if lease_path.exists():
            try:
                lease = json.loads(lease_path.read_text())
                if (lease.get("status") != "terminated" or lease.get("returncode") is None
                        or lease.get("input_sha256") != reservation["native_input_file_sha256"]):
                    raise TransferControlError("Native process lease remains live, ambiguous, or mismatched")
            except (OSError, ValueError, TypeError) as error:
                raise TransferControlError("Native process lease cannot be closed defensibly") from error
        retained_paths = {path.with_suffix(suffix) for suffix in (".inp", ".out", ".process.json", ".xyz", ".hess")}
        for pattern in (f"{path.stem}*_MEP*.xyz", f"{path.stem}*_MEP*.allxyz", f"{path.stem}*_IRC*.xyz"):
            retained_paths.update(path.parent.glob(pattern))
        for retained in sorted(retained_paths):
            if retained.is_file():
                sources.append({"path": str(retained), "sha256": _hash(retained.read_bytes())})
        self._append({"event": "terminal", "attempt_id": reservation["attempt_id"],
                      "elapsed_seconds": elapsed_seconds, "outcome": outcome,
                      "sources": sources})

    def reconcile_retained_attempt(self, attempt_id):
        """Close a controller-crash reservation only from stopped native proof.

        This never launches, retries, refunds an attempt, or estimates an
        unrecorded duration. Missing/ambiguous leases remain unresolved.
        """
        self.initialize()
        with self._exclusive():
            snapshot = self.snapshot()
            pending = [item for item in snapshot["pending"] if item["attempt_id"] == attempt_id]
            if len(pending) != 1:
                raise TransferControlError("Exactly one pending reservation is required")
            reservation = pending[0]
            path = Path(reservation["input"])
            if not path.is_file() or _hash(path.read_bytes()) != reservation["native_input_file_sha256"]:
                raise TransferControlError("Reserved native input changed")
            from .reaction_path import _orca_job_input_evidence
            if _orca_job_input_evidence(path)[0].hexdigest() != reservation["input_sha256"]:
                raise TransferControlError("Referenced native input evidence changed")
            lease_path = path.with_suffix(".process.json")
            try:
                lease = json.loads(lease_path.read_text())
                pid = lease["pid"]
                group = lease["process_group_id"]
                start, end = lease["started_at_epoch"], lease["ended_at_epoch"]
                if (lease.get("status") != "terminated" or type(pid) is not int or pid <= 0
                        or type(group) is not int or group <= 0 or lease.get("returncode") is None
                        or lease.get("input_sha256") != reservation["native_input_file_sha256"]
                        or not math.isfinite(start) or not math.isfinite(end) or end < start
                        or start < reservation["started_at_epoch"] - 1):
                    raise TransferControlError("Retained process lease does not establish a stopped owned job")
                for identifier, probe in ((pid, os.kill), (group, os.killpg)):
                    try:
                        probe(identifier, 0)
                    except ProcessLookupError:
                        continue
                    except PermissionError as error:
                        raise TransferControlError("Cannot disprove native process ownership") from error
                    raise TransferControlError("Owned native process or worker group is still live")
            except (OSError, KeyError, TypeError, ValueError) as error:
                raise TransferControlError("Native terminal process proof is missing or invalid") from error
            elapsed = end - reservation["started_at_epoch"]
            self.finish(reservation, elapsed_seconds=elapsed,
                        outcome="reconciled_stopped_native_execution_not_physically_validated")
        return self.snapshot()

    def dispatch(self, input_path, input_sha256, *, timeout_seconds, ncores, operation):
        """Charge before dispatch and close execution, including thrown errors."""
        reservation = self.reserve(input_path, input_sha256, timeout_seconds=timeout_seconds,
                                   ncores=ncores)
        started = time.monotonic()
        outcome = "failed"
        try:
            environment = dict(self.contract["extra_environment"])
            if self._execution_amendment is not None:
                environment.update(self.validate_execution_amendment()["extra_environment"])
                forwarding = environment.get("OMPI_MCA_mca_base_env_list", "")
                environment["OMPI_MCA_mca_base_env_list"] = ";".join(
                    item for item in (forwarding, "PMIX_MCA_gds=hash") if item)
            result = operation(reservation["timeout_seconds"], environment)
            outcome = "returned_native_work_not_physically_validated"
            return result
        finally:
            self.finish(reservation, elapsed_seconds=time.monotonic() - started, outcome=outcome)


def controlled_run_or_resume(input_path, *, timeout_seconds=None, explicit_single_point_settings=False):
    """Public controlled entry point, sharing the path executor and its reuse."""
    if active_execution_control() is None:
        raise TransferControlError("Use an explicit NativeExecutionControl.activate scope")
    from .reaction_path import _run_or_resume
    if explicit_single_point_settings:
        if classify_native_input(Path(input_path).read_text()) != "single_point":
            raise TransferControlError("Explicit input settings mode is restricted to single points")
        with explicit_native_input_scope():
            return _run_or_resume(Path(input_path), timeout_seconds=timeout_seconds)
    return _run_or_resume(Path(input_path), timeout_seconds=timeout_seconds)
