#!/usr/bin/env python
"""Launch RMG with process-local execution and evidence-preservation adapters.

RMG initializes its database in the parent process.  Python's macOS ``spawn``
default starts clean workers which therefore lack that database; Unix ``fork``
preserves the initialized state and is RMG's normal Linux behaviour.
External core seeds retain loader-established identities across restart loads;
gas Chebyshev and ordinary Arrhenius exports preserve the in-memory model's
numerical precision.
No installed RMG files or physical rate models are changed.
"""

from __future__ import annotations

import multiprocessing as multiprocessing
import functools
import logging
import math
import re
import sys


def chebyshev_export_adapter(original, chebyshev_type):
    """Preserve binary64 Chebyshev data in RMG's gas-phase Chemkin export.

    RMG's fixed three-decimal pressure format can turn a positive low-pressure
    bound into zero. Re-render only the polynomial block from the in-memory
    model, retaining RMG's reaction identity, comments and duplicate handling.
    This is serialization, not a fit or a change to the physical rate model.
    """
    @functools.wraps(original)
    def write(reaction, species_list, verbose=True, java_library=False, commented=False):
        text = original(reaction, species_list, verbose=verbose,
                        java_library=java_library, commented=commented)
        kinetics = reaction.kinetics
        if not isinstance(kinetics, chebyshev_type):
            return text
        bounds = [float(getattr(kinetics, name).value_si)
                  for name in ("Tmin", "Tmax", "Pmin", "Pmax")]
        if (any(not math.isfinite(value) or value <= 0.0 for value in bounds)
                or bounds[0] >= bounds[1] or bounds[2] >= bounds[3]):
            raise ValueError("Chebyshev export requires finite positive ordered T/P bounds")
        rows, columns = int(kinetics.degreeT), int(kinetics.degreeP)
        coeffs = kinetics.coeffs.value_si
        if rows < 1 or columns < 1 or tuple(coeffs.shape) != (rows, columns):
            raise ValueError("Chebyshev export coefficient shape is inconsistent")
        values = [float(coeffs[i, j]) for i in range(rows) for j in range(columns)]
        if any(not math.isfinite(value) for value in values):
            raise ValueError("Chebyshev export coefficients must be finite")
        molecularity = len(reaction.reactants)
        if molecularity not in (1, 2, 3):
            raise ValueError("Chebyshev gas export molecularity is unsupported")
        # RMG stores log10(k_SI); Chemkin's gas MOLES convention uses cm, mol, s.
        # Copy values rather than mutating the kinetic object on repeated saves.
        values[0] += 6 * (molecularity - 1)
        prefix = "!" if commented else ""
        pattern = re.compile(r"^\s*" + re.escape(prefix) + r"\s*(TCHEB|PCHEB|CHEB)\s*/")
        lines = text.splitlines(keepends=True)
        indices = [index for index, line in enumerate(lines) if pattern.match(line)]
        kinds = [pattern.match(lines[index]).group(1) for index in indices]
        if (len(indices) < 4 or indices != list(range(indices[0], indices[-1] + 1))
                or kinds[:3] != ["TCHEB", "PCHEB", "CHEB"]
                or any(kind != "CHEB" for kind in kinds[3:])):
            raise ValueError("Unrecognized RMG Chebyshev export layout; refusing lossy export")
        block = [
            f"{prefix}    TCHEB/ {bounds[0]:.17g} {bounds[1]:.17g} /\n",
            f"{prefix}    PCHEB/ {bounds[2] / 101325.0:.17g} {bounds[3] / 101325.0:.17g} /\n",
            f"{prefix}    CHEB/ {rows} {columns} /\n",
        ]
        for start in range(0, len(values), 3):
            block.append(prefix + "    CHEB/ " + " ".join(
                format(value, ".17g") for value in values[start:start + 3]
            ) + " /\n")
        return "".join(lines[:indices[0]] + block + lines[indices[-1] + 1:])
    return write


def arrhenius_export_adapter(original, arrhenius_type):
    """Preserve ordinary gas Arrhenius parameters in process-local exports.

    RMG normally writes n and Ea with only three decimals.  This can move a
    retained verified rate beyond its reinjection tolerance despite preserving
    reaction identity.  Replace only that numeric triplet using the same T0,
    SI-to-cm/mol/s and J-to-kcal/mol conversions as RMG; do not refit, mutate the
    model, patch installed packages, or relax a scientific acceptance check.
    Specialized/subclass kinetics and Java-library exports are unchanged.
    """
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
    parameter_line = re.compile(
        r"^(?P<head>.*(?:<=>|=>|=).*?[ \t]+)" + number + r"[ \t]+"
        + number + r"[ \t]+" + number
        + r"(?P<tail>[ \t]*(?:![^\r\n]*)?)(?P<ending>\r?\n|)$"
    )

    @functools.wraps(original)
    def write(reaction, species_list, verbose=True, java_library=False, commented=False):
        text = original(reaction, species_list, verbose=verbose,
                        java_library=java_library, commented=commented)
        kinetics = reaction.kinetics
        if type(kinetics) is not arrhenius_type or java_library:
            return text
        surface_check = getattr(reaction, "is_surface_reaction", None)
        if callable(surface_check) and surface_check():
            return text
        molecularity = len(reaction.reactants)
        if molecularity not in (1, 2, 3):
            raise ValueError("Ordinary Arrhenius gas export molecularity is unsupported")
        try:
            preexponential = float(kinetics.A.value_si)
            exponent = float(kinetics.n.value_si)
            activation_energy = float(kinetics.Ea.value_si)
            reference_temperature = float(kinetics.T0.value_si)
            conversion = float(kinetics.A.get_conversion_factor_from_si_to_cm_mol_s())
            if (not all(math.isfinite(value) for value in
                        (preexponential, exponent, activation_energy, reference_temperature, conversion))
                    or reference_temperature <= 0.0 or conversion <= 0.0):
                raise ValueError("parameters must be finite with positive T0 and conversion")
            if not math.isclose(conversion, 1e6 ** (molecularity - 1), rel_tol=1e-12):
                raise ValueError("rate units are inconsistent with gas molecularity")
            values = (preexponential / reference_temperature ** exponent * conversion,
                      exponent, activation_energy / 4184.0)
            if any(not math.isfinite(value) for value in values):
                raise ValueError("converted parameters must be finite")
        except (AttributeError, TypeError, ValueError, OverflowError, ZeroDivisionError) as error:
            raise ValueError("Invalid ordinary Arrhenius gas export: %s" % error) from error
        lines = text.splitlines(keepends=True)
        matched = []
        for index, line in enumerate(lines):
            if commented:
                if not line.startswith("!"):
                    continue
                body = line[1:]
            else:
                if line.lstrip().startswith("!"):
                    continue
                body = line
            if body.lstrip().startswith("!"):
                continue
            match = parameter_line.match(body)
            if match is not None:
                matched.append((index, match))
        if len(matched) != 1:
            raise ValueError("Unrecognized RMG ordinary Arrhenius export layout; refusing lossy export")
        index, match = matched[0]
        lines[index] = (("!" if commented else "") + match.group("head")
                        + " ".join(format(value, ".17g") for value in values)
                        + match.group("tail") + match.group("ending"))
        return "".join(lines)
    return write


def library_load_adapter(original):
    """Retain loader-established aliases across RMG's restart-library load.

    RMG resets external_library_labels on every load_libraries call, even
    while keeping previously loaded libraries. Bind retained aliases to the
    actual library object so an overwritten short name cannot be misused.
    """
    @functools.wraps(original)
    def load(kinetics, *args, **kwargs):
        result = original(kinetics, *args, **kwargs)
        bindings = dict(getattr(kinetics, "_storca_external_seed_bindings", {}))
        for source, key in getattr(kinetics, "external_library_labels", {}).items():
            if key not in kinetics.libraries:
                raise ValueError("RMG external library alias points to an unloaded library")
            bindings[source] = (key, kinetics.libraries[key])
        kinetics._storca_external_seed_bindings = bindings
        return result
    return load


def seed_core_adapter(original, kinetics_provider):
    """Resolve external seeds through RMG's own loader registry, never basename guessing."""
    @functools.wraps(original)
    def add_seed(model, seed_mechanism, *args, **kwargs):
        kinetics = kinetics_provider()
        key = seed_mechanism
        if key not in kinetics.libraries:
            aliases = getattr(kinetics, "external_library_labels", {})
            bindings = getattr(kinetics, "_storca_external_seed_bindings", {})
            if key in aliases or key in bindings:
                key = aliases[key] if key in aliases else bindings[key][0]
                if key not in kinetics.libraries:
                    raise ValueError("RMG external seed alias points to an unloaded library")
                retained = bindings.get(seed_mechanism)
                if retained is not None and (retained[0] != key or retained[1] is not kinetics.libraries[key]):
                    raise ValueError("RMG external seed library identity changed after loading")
                logging.info("STORCA external seed registry resolution: %s -> %s", seed_mechanism, key)
        return original(model, key, *args, **kwargs)
    return add_seed


def main() -> None:
    requested = "auto"
    if "--storca-start-method" in sys.argv:
        index = sys.argv.index("--storca-start-method")
        try:
            requested = sys.argv[index + 1]
        except IndexError as error:
            raise SystemExit("--storca-start-method requires a value") from error
        del sys.argv[index:index + 2]
    if requested not in {"auto", "fork", "spawn"}:
        raise SystemExit("STORCA RMG start method must be auto, fork, or spawn")
    method = "fork" if requested == "auto" and sys.platform == "darwin" else requested
    if method != "auto":
        multiprocessing.set_start_method(method)
    # Import only after selecting the process start method.
    from rmgpy.__main__ import main as rmg_main
    import rmgpy.data.rmg
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.data.kinetics.database import KineticsDatabase
    import rmgpy.chemkin
    from rmgpy.kinetics import Arrhenius, Chebyshev
    rmgpy.chemkin.write_kinetics_entry = chebyshev_export_adapter(
        rmgpy.chemkin.write_kinetics_entry, Chebyshev,
    )
    logging.info("STORCA gas Chebyshev export policy: binary64-precision-v1")
    rmgpy.chemkin.write_kinetics_entry = arrhenius_export_adapter(
        rmgpy.chemkin.write_kinetics_entry, Arrhenius,
    )
    logging.info("STORCA ordinary gas Arrhenius export policy: binary64-precision-v1")
    KineticsDatabase.load_libraries = library_load_adapter(KineticsDatabase.load_libraries)
    CoreEdgeReactionModel.add_seed_mechanism_to_core = seed_core_adapter(
        CoreEdgeReactionModel.add_seed_mechanism_to_core,
        lambda: rmgpy.data.rmg.database.kinetics,
    )
    rmg_main()


if __name__ == "__main__":
    main()
