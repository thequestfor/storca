"""Explicit, scoped native controls for the bounded rate-transfer experiment.

This opt-in object carries controls that cannot be represented by a method
keyword list. Outside ``native_settings_scope`` the old builders are unchanged.
It is not a method selector or an applicability certificate.
"""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import re


@dataclass(frozen=True)
class NativeJobSettings:
    maxcore_mb_per_rank: int = 1000
    frequency_method: str = "NumFreq"
    central_diff: bool = True
    increment_bohr: float = 0.001
    frequency_restart: bool = False
    frequency_scale: float = 1.0
    project_translation_rotation: bool = True
    quasi_rrho: bool = False
    cnvsoscf: bool = False
    conv_check_mode: int = 0
    conv_forced: int = 1
    scf_max_iterations: int = 500
    optimization_convergence: str = "VeryTightOpt"

    def __post_init__(self):
        for field in ("maxcore_mb_per_rank", "scf_max_iterations"):
            if type(getattr(self, field)) is not int or getattr(self, field) < 1:
                raise ValueError(f"{field} must be a positive integer")
        for field in ("central_diff", "frequency_restart", "project_translation_rotation",
                      "quasi_rrho", "cnvsoscf"):
            if type(getattr(self, field)) is not bool:
                raise ValueError(f"{field} must be boolean")
        if self.frequency_method != "NumFreq":
            raise ValueError("This transfer profile requires explicit NumFreq")
        if self.optimization_convergence != "VeryTightOpt":
            raise ValueError("The frozen transfer profile requires VeryTightOpt")
        for field in ("increment_bohr", "frequency_scale"):
            value = getattr(self, field)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{field} must be finite and positive")
        if type(self.conv_check_mode) is not int or self.conv_check_mode != 0:
            raise ValueError("The frozen transfer SCF profile requires ConvCheckMode 0")
        if type(self.conv_forced) is not int or self.conv_forced != 1:
            raise ValueError("The frozen transfer SCF profile requires ConvForced 1")

    def as_dict(self):
        return {"schema_version": 1, "kind": "storca_native_job_settings", **asdict(self)}

    @classmethod
    def from_dict(cls, value):
        value = dict(value)
        if value.pop("schema_version", None) != 1 or value.pop("kind", None) != "storca_native_job_settings":
            raise ValueError("Unsupported native settings contract")
        return cls(**value)

    @property
    def sha256(self):
        return hashlib.sha256(json.dumps(self.as_dict(), sort_keys=True,
                                        separators=(",", ":"), allow_nan=False).encode()).hexdigest()


_SETTINGS = ContextVar("storca_explicit_native_job_settings", default=None)


def active_native_settings():
    return _SETTINGS.get()


@contextmanager
def native_settings_scope(settings):
    if not isinstance(settings, NativeJobSettings):
        raise TypeError("An explicit NativeJobSettings object is required")
    if _SETTINGS.get() is not None:
        raise RuntimeError("Nested native settings scopes are not supported")
    token = _SETTINGS.set(settings)
    try:
        yield settings
    finally:
        _SETTINGS.reset(token)


@contextmanager
def explicit_native_input_scope():
    """An already complete SP input, not a replacement implicit method profile."""
    token = _SETTINGS.set(None)
    try:
        yield
    finally:
        _SETTINGS.reset(token)


def apply_native_job_settings(text, native_settings=None):
    """Render the complete profile, refusing conflicting retained controls.

    Simple frequency keywords are normalized to exactly one NumFreq. Existing
    temperature-only %freq blocks are merged; the matching fully rendered
    form is idempotent. Other conflicting SCF/frequency blocks are refused.
    """
    settings = native_settings or active_native_settings()
    if settings is None:
        return text
    if not isinstance(settings, NativeJobSettings):
        raise TypeError("native_settings must be NativeJobSettings")
    lines = text.splitlines()
    keyword_indices = [index for index, line in enumerate(lines) if line.lstrip().startswith("!")]
    if not keyword_indices:
        raise ValueError("Native settings require an ORCA simple-keyword line")
    frequency_tokens = []
    for index in keyword_indices:
        words = lines[index].split()
        for word in words[1:]:
            if word.casefold() in {"freq", "anfreq", "numfreq"}:
                frequency_tokens.append(word)
        lines[index] = " ".join(words)
    if len(frequency_tokens) > 1:
        raise ValueError("Conflicting duplicate frequency job keywords")
    frequency = bool(frequency_tokens)
    all_tokens = {word.casefold() for index in keyword_indices for word in lines[index].split()[1:]}
    # Legacy path builders append TightSCF to a supplied method. It must not
    # override that method's explicitly tighter SCF convergence keyword.
    if "extremescf" in all_tokens:
        for index in keyword_indices:
            lines[index] = re.sub(r"\b(?:TightSCF|VeryTightSCF)\b", "", lines[index], flags=re.I)
            lines[index] = " ".join(lines[index].split())
    if all_tokens & {"opt", "optts", "neb-ts", "neb", "neb-ci", "neb-ci-ts", "copt"}:
        for index in keyword_indices:
            lines[index] = re.sub(r"\b(?:LooseOpt|NormalOpt|TightOpt|VeryTightOpt)\b", "", lines[index], flags=re.I)
            lines[index] = " ".join(lines[index].split())
        lines[keyword_indices[0]] += " " + settings.optimization_convergence
    if frequency:
        index = next(index for index in keyword_indices if re.search(
            r"\b(?:Freq|AnFreq|NumFreq)\b", lines[index], re.I))
        lines[index] = re.sub(r"\b(?:Freq|AnFreq|NumFreq)\b", settings.frequency_method,
                              lines[index], flags=re.I)
    text = "\n".join(lines)
    maxcore = f"%maxcore {settings.maxcore_mb_per_rank}"
    declarations = re.findall(r"(?im)^\s*%maxcore\s+([^\n]+)", text)
    if declarations and (len(declarations) != 1 or declarations[0].strip() != str(settings.maxcore_mb_per_rank)):
        raise ValueError("Native MaxCore conflicts with the frozen settings")
    extras = [] if declarations else [maxcore]
    boolean = lambda value: "true" if value else "false"
    scf = (f"%scf\n  CNVSOSCF {boolean(settings.cnvsoscf)}\n"
           f"  ConvCheckMode {settings.conv_check_mode}\n  ConvForced {settings.conv_forced}\n"
           f"  MaxIter {settings.scf_max_iterations}\nend")
    blocks = re.findall(r"%scf\b.*?\bend\b", text, re.I | re.S)
    if blocks:
        if len(blocks) != 1 or " ".join(blocks[0].casefold().split()) != " ".join(scf.casefold().split()):
            raise ValueError("Native SCF controls conflict with the frozen settings")
    else:
        extras.append(scf)
    if frequency:
        controls = (f"  CentralDiff {boolean(settings.central_diff)}\n"
                    f"  Increment {settings.increment_bohr:.15g}\n"
                    f"  Restart {boolean(settings.frequency_restart)}\n"
                    f"  ScalFreq {settings.frequency_scale:.15g}\n"
                    f"  ProjectTR {boolean(settings.project_translation_rotation)}\n"
                    f"  QuasiRRHO {boolean(settings.quasi_rrho)}")
        blocks = re.findall(r"%freq\b(.*?)\bend\b", text, re.I | re.S)
        if len(blocks) > 1:
            raise ValueError("Duplicate native frequency blocks")
        if blocks:
            body = blocks[0]
            temperatures = re.findall(r"(?im)^\s*Temp\s+([-+0-9.Ee]+)\s*$", body)
            stripped = re.sub(r"(?im)^\s*Temp\s+[-+0-9.Ee]+\s*$", "", body).strip()
            if stripped and " ".join(stripped.casefold().split()) != " ".join(controls.casefold().split()):
                raise ValueError("Native frequency controls conflict with the frozen settings")
            if len(temperatures) > 1:
                raise ValueError("Multiple native thermochemistry temperatures")
            temp = f"\n  Temp {temperatures[0]}" if temperatures else ""
            text = re.sub(r"%freq\b.*?\bend\b", "%freq\n" + controls + temp + "\nend",
                          text, flags=re.I | re.S)
        else:
            extras.append("%freq\n" + controls + "\nend")
    # Put controls before coordinates, rather than after the geometry block.
    if extras:
        coordinate = re.search(r"(?m)^\s*\*\s+xyz(?:file)?\b", text)
        position = coordinate.start() if coordinate else len(text)
        text = text[:position].rstrip() + "\n" + "\n".join(extras) + "\n" + text[position:]
    return text.rstrip() + "\n"
