"""Durable attempt accounting for the bounded OH/methane angular experiment.

Reservations, terminal execution receipts and derived analyses are separate,
immutable hash-chained events. A reservation spends budget even without a
launched job. Incomplete/corrupt writes block further spending, never disappear
from accounting. This is not a process supervisor or scientific validator.
Callers own native leases, evidence hashes and permission to retry/reanalyse.

Use a dedicated ledger directory and retain returned record references in the
campaign's sealed evidence. Local hash chains detect mutation and gaps, but no
local ledger can prove that an entire trailing history was not rolled back;
externally bound references are the anchor for that case.
"""
from __future__ import annotations

from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import threading


_KINDS = {"saddles", "NumFreq"}
_HASH = re.compile(r"[0-9a-f]{64}\Z")
_EVENT_FILE = re.compile(r"(\d{6})-([0-9a-f]{64})\.json\Z")
_TERMINAL = {"normal", "completed", "failed"}


class AngularAttemptError(ValueError):
    """Accounting, ownership or retained evidence cannot safely be established."""


def _encode(value):
    try:
        return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                           allow_nan=False) + "\n").encode("utf-8")
    except (TypeError, ValueError) as error:
        raise AngularAttemptError("Strict finite JSON record required") from error


def _digest(content):
    return hashlib.sha256(content).hexdigest()


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise AngularAttemptError("Duplicate JSON key in ledger")
        result[key] = value
    return result


def _decode(content):
    try:
        value = json.loads(content, object_pairs_hook=_unique_object,
                           parse_constant=lambda _: (_ for _ in ()).throw(
                               AngularAttemptError("Nonfinite JSON in ledger")))
    except (UnicodeError, ValueError) as error:
        raise AngularAttemptError("Incomplete or invalid ledger record") from error
    if _encode(value) != content:
        raise AngularAttemptError("Ledger record bytes are not canonical")
    return value


def _metadata(value):
    if not isinstance(value, dict):
        raise AngularAttemptError("Receipt/metadata must be a JSON object")

    def keys(item):
        if isinstance(item, dict):
            if any(not isinstance(key, str) for key in item):
                raise AngularAttemptError("JSON object keys must be strings")
            for child in item.values():
                keys(child)
        elif isinstance(item, (list, tuple)):
            for child in item:
                keys(child)

    keys(value)
    return _decode(_encode(value))


def _fsync_directory(path):
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _exclusive_write(path, content):
    # A partial final-name record is deliberate: interruption must block resume,
    # not leave an ignored temporary file and silently refund an attempt.
    with path.open("xb") as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())
    _fsync_directory(path.parent)


def _read(path):
    if path.is_symlink() or not path.is_file():
        raise AngularAttemptError("Ledger record missing or not a regular file")
    return path.read_bytes()


class AttemptLedger:
    """One immutable campaign identity and an append-only event journal.

    All mutations require ``with ledger.controller():``. The stable lock inode
    is never unlinked; flock exclusion therefore survives competing opens and
    is released by the OS if the controller dies. No PID-based stale recovery,
    process termination or automatic job retry is performed.
    """

    def __init__(self, root, *, prior, caps, contract_sha256):
        for name, values in (("prior", prior), ("caps", caps)):
            if (not isinstance(values, dict) or set(values) != _KINDS
                    or any(type(value) is not int or value < 0 for value in values.values())):
                raise AngularAttemptError(f"Exact nonnegative integer {name} required")
        if any(prior[key] > caps[key] for key in _KINDS):
            raise AngularAttemptError("Prior spending exceeds campaign cap")
        if not isinstance(contract_sha256, str) or not _HASH.fullmatch(contract_sha256):
            raise AngularAttemptError("Exact campaign contract SHA256 required")
        self.root = Path(root).resolve()
        self.events = self.root / "events"
        self.identity_path = self.root / "identity.json"
        self.lock_path = self.root / "controller.lock"
        self._identity = {"schema_version": 1, "kind": "OH_CH4_angular_attempt_ledger",
                          "contract_sha256": contract_sha256,
                          "prior": dict(prior), "caps": dict(caps)}
        self._identity_bytes = _encode(self._identity)
        self._identity_hash = _digest(self._identity_bytes)
        self._owner = None
        self._lock_fd = None

    def initialize(self):
        """Create or reopen the exact identity; never reset an existing ledger."""
        self.root.mkdir(parents=True, exist_ok=True)
        _fsync_directory(self.root.parent)
        with self._controller_lock():
            if self.identity_path.exists():
                self._check_identity()
                if self.events.is_symlink() or not self.events.is_dir():
                    raise AngularAttemptError("Existing ledger events directory missing")
            else:
                if self.events.exists() and (self.events.is_symlink()
                                             or not self.events.is_dir()
                                             or any(self.events.iterdir())):
                    raise AngularAttemptError("Missing identity cannot reset retained events")
                self.events.mkdir(exist_ok=True)
                _fsync_directory(self.root)
                _exclusive_write(self.identity_path, self._identity_bytes)
            self._scan()
        return {"path": str(self.identity_path), "sha256": self._identity_hash}

    @contextmanager
    def _controller_lock(self):
        if self._owner is not None:
            raise AngularAttemptError("Controller already held by this ledger object")
        if not self.root.is_dir():
            raise AngularAttemptError("Initialize ledger before taking controller")
        flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0)
        descriptor = os.open(self.lock_path, flags, 0o600)
        acquired = False
        try:
            if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                raise AngularAttemptError("Controller lock is not a regular file")
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise AngularAttemptError("Another campaign controller is active") from error
            _fsync_directory(self.root)
            self._owner = (os.getpid(), threading.get_ident())
            self._lock_fd = descriptor
            acquired = True
            yield
        finally:
            if acquired:
                self._owner = None
                self._lock_fd = None
            os.close(descriptor)  # Stable lockfile stays; OS releases this flock.

    @contextmanager
    def controller(self):
        """Take exclusive controller ownership, without stale/PID recovery."""
        with self._controller_lock():
            self._scan()
            yield self

    def _require_controller(self):
        if self._owner != (os.getpid(), threading.get_ident()) or self._lock_fd is None:
            raise AngularAttemptError("Mutation requires this thread's controller context")
        current = self.lock_path.stat()
        owned = os.fstat(self._lock_fd)
        if (current.st_dev, current.st_ino) != (owned.st_dev, owned.st_ino):
            raise AngularAttemptError("Controller lock inode changed")

    def _check_identity(self):
        if _read(self.identity_path) != self._identity_bytes:
            raise AngularAttemptError("Immutable ledger identity differs")

    def _scan(self):
        self._check_identity()
        if self.events.is_symlink() or not self.events.is_dir():
            raise AngularAttemptError("Ledger events directory missing")
        names = sorted(path.name for path in self.events.iterdir())
        records, reservations, executions, analyses = [], {}, {}, {}
        counts = dict(self._identity["prior"])
        previous = self._identity_hash
        for number, name in enumerate(names, 1):
            match = _EVENT_FILE.fullmatch(name)
            if match is None or int(match[1]) != number:
                raise AngularAttemptError("Unexpected, missing or incomplete ledger event")
            path = self.events / name
            content = _read(path)
            if _digest(content) != match[2]:
                raise AngularAttemptError("Ledger event changed or interrupted write")
            record = _decode(content)
            base = {"schema_version", "ledger_sha256", "sequence", "previous_sha256", "event"}
            if (not isinstance(record, dict) or record.get("schema_version") != 1
                    or record.get("ledger_sha256") != self._identity_hash
                    or type(record.get("sequence")) is not int or record["sequence"] != number
                    or record.get("previous_sha256") != previous):
                raise AngularAttemptError("Ledger event identity/chain differs")
            reference = {"path": str(path), "sha256": match[2],
                         "sequence": number, "id": path.stem}
            event = record["event"]
            if event == "reserve":
                if (set(record) != base | {"point_id", "kind", "metadata"}
                        or not isinstance(record["kind"], str) or record["kind"] not in _KINDS):
                    raise AngularAttemptError("Malformed reservation")
                self._point(record["point_id"])
                _metadata(record["metadata"])
                if any((row["point_id"], row["kind"]) == (record["point_id"], record["kind"])
                       for row in reservations.values()):
                    raise AngularAttemptError("Duplicate point/kind reservation")
                reservations[path.stem] = {**reference, "point_id": record["point_id"],
                                            "kind": record["kind"], "metadata": record["metadata"]}
                counts[record["kind"]] += 1
            elif isinstance(event, str) and event in {"execution", "analysis"}:
                if (set(record) != base | {"attempt_id", "payload"}
                        or not isinstance(record["attempt_id"], str)
                        or record["attempt_id"] not in reservations):
                    raise AngularAttemptError("Receipt has no earlier reservation")
                payload = _metadata(record["payload"])
                target = executions if event == "execution" else analyses
                if record["attempt_id"] in target:
                    raise AngularAttemptError("Immutable receipt recorded twice")
                if event == "analysis" and record["attempt_id"] not in executions:
                    raise AngularAttemptError("Analysis lacks terminal execution receipt")
                if (event == "execution" and (not isinstance(payload.get("status"), str)
                                              or payload["status"] not in _TERMINAL)):
                    raise AngularAttemptError("Explicit terminal execution status required")
                target[record["attempt_id"]] = {**reference, "payload": payload}
            else:
                raise AngularAttemptError("Unknown ledger event")
            records.append({**reference, "record": record})
            previous = match[2]
        if any(counts[key] > self._identity["caps"][key] for key in _KINDS):
            raise AngularAttemptError("Retained reservations exceed budget")
        if names != sorted(path.name for path in self.events.iterdir()):
            raise AngularAttemptError("Ledger changed while reading; reopen it")
        return records, reservations, executions, analyses, counts

    @staticmethod
    def _point(point_id):
        if (not isinstance(point_id, str) or not point_id or len(point_id) > 160
                or point_id != point_id.strip() or point_id in {".", ".."}
                or any(char in point_id for char in "/\\")
                or any(ord(char) < 32 for char in point_id)):
            raise AngularAttemptError("Explicit safe point identifier required")

    @staticmethod
    def _id(attempt):
        value = attempt.get("id") if isinstance(attempt, dict) else attempt
        if not isinstance(value, str) or not _EVENT_FILE.fullmatch(value + ".json"):
            raise AngularAttemptError("Retained attempt identifier required")
        return value

    def _append(self, event, **fields):
        self._require_controller()
        records, *_ = self._scan()
        number = len(records) + 1
        if number > 999999:
            raise AngularAttemptError("Experiment ledger event bound exceeded")
        record = {"schema_version": 1, "ledger_sha256": self._identity_hash,
                  "sequence": number,
                  "previous_sha256": records[-1]["sha256"] if records else self._identity_hash,
                  "event": event, **fields}
        content = _encode(record)
        digest = _digest(content)
        path = self.events / f"{number:06d}-{digest}.json"
        _exclusive_write(path, content)
        return {"path": str(path), "sha256": digest, "sequence": number, "id": path.stem}

    def reservations(self):
        return list(self._scan()[1].values())

    def counts(self):
        """Whole-increment spending, including unlaunched/reserved attempts."""
        return self._scan()[4]

    def evidence(self):
        """Bind the exact immutable history externally, excluding the lockfile."""
        records = self._scan()[0]
        return [{"path": str(self.identity_path), "sha256": self._identity_hash},
                *({"path": row["path"], "sha256": row["sha256"]} for row in records)]

    def reserve(self, point_id, kind, metadata):
        self._require_controller()
        self._point(point_id)
        if kind not in _KINDS:
            raise AngularAttemptError("Known angular attempt kind required")
        copied = _metadata(metadata)  # Validation precedes the durable spend.
        _, reservations, executions, _, counts = self._scan()
        if any((row["point_id"], row["kind"]) == (point_id, kind)
               for row in reservations.values()):
            raise AngularAttemptError("Point/kind already reserved; resume without redispatch")
        if set(reservations) - set(executions):
            raise AngularAttemptError("Unresolved reserved execution must be reconciled first")
        if counts[kind] >= self._identity["caps"][kind]:
            raise AngularAttemptError("Cumulative angular attempt cap reached")
        reference = self._append("reserve", point_id=point_id, kind=kind, metadata=copied)
        return {**reference, "point_id": point_id, "kind": kind, "metadata": copied}

    def _receipt(self, attempt, event):
        attempt_id = self._id(attempt)
        _, reservations, executions, analyses, _ = self._scan()
        if attempt_id not in reservations:
            raise AngularAttemptError("Attempt is not retained in this ledger")
        receipt = (executions if event == "execution" else analyses).get(attempt_id)
        return None if receipt is None else receipt["payload"]

    def execution(self, attempt):
        return self._receipt(attempt, "execution")

    def analysis(self, attempt):
        return self._receipt(attempt, "analysis")

    def _seal(self, attempt, payload, event):
        self._require_controller()
        attempt_id, copied = self._id(attempt), _metadata(payload)
        _, reservations, executions, analyses, _ = self._scan()
        if attempt_id not in reservations:
            raise AngularAttemptError("Attempt is not retained in this ledger")
        if event == "analysis" and attempt_id not in executions:
            raise AngularAttemptError("Seal terminal execution before derived analysis")
        if (event == "execution" and (not isinstance(copied.get("status"), str)
                                      or copied["status"] not in _TERMINAL)):
            raise AngularAttemptError("Explicit terminal execution status required; ambiguous stays reserved")
        existing = (executions if event == "execution" else analyses).get(attempt_id)
        if existing is not None:
            if _encode(existing["payload"]) != _encode(copied):
                raise AngularAttemptError("Immutable receipt differs; cannot replace it")
            return {key: value for key, value in existing.items() if key != "payload"}
        return self._append(event, attempt_id=attempt_id, payload=copied)

    def seal_execution(self, attempt, payload):
        """Seal caller-validated status normal/completed/failed, never ambiguity."""
        return self._seal(attempt, payload, "execution")

    def seal_analysis(self, attempt, payload):
        """Seal one derived receipt; parser failure can leave this absent for replay."""
        return self._seal(attempt, payload, "analysis")

    def reopen(self, reference):
        """Reopen an externally retained exact event reference, detecting rollback."""
        if (not isinstance(reference, dict)
                or type(reference.get("sequence")) is not int):
            raise AngularAttemptError("Exact event reference required")
        for row in self._scan()[0]:
            if all(reference.get(key) == row[key] for key in ("path", "sha256", "sequence", "id")):
                return row["record"]
        raise AngularAttemptError("Retained event reference missing or changed")
