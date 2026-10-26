#!/usr/bin/env python
"""RMG-environment side of STORCA's JSON artifact bridge.

This file intentionally has no STORCA imports: it is executed by RMG 3/Python
3.9 and RMG 4/Python 3.11 as well as the app's Python 3.12 environment.
"""

from __future__ import print_function

import json
import hashlib
import math
import sys
from pathlib import Path


def _species_identifier(species, get_species_identifier):
    try:
        return get_species_identifier(species)
    except Exception:
        return species.label


def _finite_nonnegative(value):
    """Optional numerical evidence is unavailable, not physical, if invalid."""
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) and number >= 0.0 else None


def _equilibrium_constant_evidence(reaction, temperature):
    """Evaluate native thermochemistry directly, never a fitted reverse rate.

    RMG's Kc concentration convention can be dimensional when stoichiometric
    molecule counts differ.  This record does not claim rate calibration,
    pressure applicability, or conformer-equilibrated species thermochemistry.
    """
    result = {
        "evaluated_equilibrium_constant_Kc": None,
        "equilibrium_constant_unavailable_reason": None,
        "equilibrium_constant_evaluated_at": None,
        "equilibrium_constant_source": "RMG Reaction.get_equilibrium_constant(T, type='Kc')",
        "equilibrium_constant_basis": "native_species_thermochemistry",
        "equilibrium_constant_convention": (
            "RMG concentration equilibrium Kc; dimensionality depends on stoichiometry and environment"
        ),
    }
    if temperature is None:
        result["equilibrium_constant_unavailable_reason"] = "temperature_not_requested"
        return result
    try:
        if isinstance(temperature, bool):
            raise ValueError("Temperature cannot be boolean")
        evaluated_temperature = float(temperature)
        if not math.isfinite(evaluated_temperature) or evaluated_temperature <= 0.0:
            raise ValueError("Temperature must be positive and finite")
    except (TypeError, ValueError, OverflowError):
        result["equilibrium_constant_unavailable_reason"] = "equilibrium_temperature_nonpositive_nonfinite_or_unavailable"
        return result
    result["equilibrium_constant_evaluated_at"] = {"temperature_K": evaluated_temperature}
    try:
        raw = reaction.get_equilibrium_constant(evaluated_temperature, type="Kc")
        if isinstance(raw, (bool, str, bytes)):
            raise ValueError("Equilibrium constant must be numerical")
        number = _finite_nonnegative(raw)
        if number is None or number <= 0.0:
            result["equilibrium_constant_unavailable_reason"] = "equilibrium_constant_nonpositive_nonfinite_or_unavailable"
        else:
            result["evaluated_equilibrium_constant_Kc"] = number
    except Exception as error:
        result["equilibrium_constant_unavailable_reason"] = "equilibrium_constant_evaluation_failed: %s" % type(error).__name__
    return result


def _inspection_model_artifacts(payload):
    """Optional source snapshot; absent mock/legacy paths never abort inspect."""
    result = {}
    try:
        for key in ("chemkin", "dictionary", "transport"):
            value = payload.get(key)
            if value is None:
                continue
            path = Path(value).resolve(strict=True)
            if not path.is_file():
                return None
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
            result[key] = {"path": str(path), "sha256": digest.hexdigest()}
    except (OSError, TypeError, ValueError):
        return None
    return result if "chemkin" in result else None


def inspect(payload):
    from rmgpy.chemkin import load_chemkin_file, get_species_identifier
    equilibrium_model_artifacts_before = _inspection_model_artifacts(payload)
    identity_payload = None
    if payload.get("replacement_libraries"):
        try:
            from .rmg_identity_bridge import build_replacement_correspondences, mechanism_artifacts
        except ImportError:  # Executed directly inside the owning RMG environment.
            from rmg_identity_bridge import build_replacement_correspondences, mechanism_artifacts
        identity_payload = dict(payload)
        identity_payload["_identity_mechanism_artifacts_before"] = mechanism_artifacts(payload)
    species, reactions = load_chemkin_file(
        payload["chemkin"], dictionary_path=payload.get("dictionary"),
        transport_path=payload.get("transport"),
    )
    records = []
    for item in species:
        molecule = item.molecule[0] if item.molecule else None
        net_charge = (sum(int(getattr(atom, "charge", 0)) for atom in molecule.atoms)
                      if molecule is not None else None)
        records.append({"label": item.label, "chemkin_identifier": _species_identifier(item, get_species_identifier),
                        "multiplicity": item.multiplicity, "reactive": bool(item.reactive),
                        "smiles": molecule.to_smiles() if molecule is not None else None,
                        "net_charge": net_charge})
    temperature = payload.get("temperature_K")
    pressure_bar = payload.get("pressure_bar")
    reaction_records = []
    for reaction_index, reaction in enumerate(reactions):
        kinetics = getattr(reaction, "kinetics", None)
        reactant_labels = [_species_identifier(item, get_species_identifier) for item in reaction.reactants]
        product_labels = [_species_identifier(item, get_species_identifier) for item in reaction.products]
        arrow = "<=>" if reaction.reversible else "=>"
        equation = "+".join(reactant_labels) + arrow + "+".join(product_labels)
        molecularity = len(reaction.reactants)
        reverse_molecularity = len(reaction.products)
        rate = None
        reverse_rate = None
        collision_limit = None
        reverse_collision_limit = None
        reverse_rate_unavailable_reason = None
        if temperature is not None:
            try:
                try:
                    rate = reaction.get_rate_coefficient(float(temperature), float(pressure_bar) * 1e5)
                except TypeError:
                    rate = reaction.get_rate_coefficient(float(temperature))
                rate = _finite_nonnegative(rate)
                if rate is None:
                    raise ValueError("forward rate coefficient is non-finite, negative, or unavailable")
            except Exception as error:
                raise ValueError(
                    "Reaction index %d (%s), kinetics %s: %s" %
                    (reaction_index, equation, type(kinetics).__name__, error)
                ) from error
            if bool(reaction.reversible):
                try:
                    reverse_kinetics = reaction.generate_reverse_rate_coefficient()
                    try:
                        reverse_rate = reverse_kinetics.get_rate_coefficient(
                            float(temperature), float(pressure_bar) * 1e5,
                        )
                    except TypeError:
                        reverse_rate = reverse_kinetics.get_rate_coefficient(float(temperature))
                    reverse_rate = _finite_nonnegative(reverse_rate)
                    if reverse_rate is None:
                        reverse_rate_unavailable_reason = "reverse_rate_nonfinite_negative_or_unavailable"
                except Exception as error:
                    reverse_rate = None
                    reverse_rate_unavailable_reason = "reverse_rate_evaluation_failed: %s" % type(error).__name__
            if molecularity == 2:
                try:
                    collision_limit = reaction.calculate_coll_limit(temp=float(temperature), reverse=False)
                # Loading a retained Chemkin file outside the full RMG driver
                # may leave the transport database uninitialized.  The rate
                # inspection remains valid; collision validation then stays
                # explicitly unavailable instead of aborting the bridge.
                except Exception:
                    collision_limit = None
            if reverse_molecularity == 2:
                try:
                    reverse_collision_limit = reaction.calculate_coll_limit(
                        temp=float(temperature), reverse=True,
                    )
                except Exception:
                    reverse_collision_limit = None
        collision_limit = _finite_nonnegative(collision_limit)
        reverse_collision_limit = _finite_nonnegative(reverse_collision_limit)
        # A zero collision ceiling is not a usable positive physical bound.
        collision_limit = collision_limit if collision_limit else None
        reverse_collision_limit = reverse_collision_limit if reverse_collision_limit else None
        units = {1: "s^-1", 2: "m^3/(mol*s)", 3: "m^6/(mol^2*s)"}.get(molecularity)
        reverse_units = {1: "s^-1", 2: "m^3/(mol*s)", 3: "m^6/(mol^2*s)"}.get(reverse_molecularity)
        def endpoint(item):
            molecule = item.molecule[0] if item.molecule else None
            return {
                "label": _species_identifier(item, get_species_identifier),
                "smiles": molecule.to_smiles() if molecule is not None else None,
                "multiplicity": item.multiplicity,
                "net_charge": (sum(int(getattr(atom, "charge", 0)) for atom in molecule.atoms)
                               if molecule is not None else None),
            }

        reaction_records.append({
            "reaction_index": reaction_index,
            "route_id": "network-%d" % (reaction_index + 1),
            "equation": equation,
            "reaction_equation": equation,
            "reactant_labels": reactant_labels,
            "product_labels": product_labels,
            "resolved_endpoints": {
                "reactants": [endpoint(item) for item in reaction.reactants],
                "products": [endpoint(item) for item in reaction.products],
            },
            "reversible": reaction.reversible,
            "molecularity": molecularity,
            "reverse_molecularity": reverse_molecularity,
            "evaluated_forward_rate_coefficient_si": rate,
            "evaluated_rate_units": units,
            "evaluated_reverse_rate_coefficient_si": reverse_rate,
            "reverse_rate_unavailable_reason": reverse_rate_unavailable_reason,
            "evaluated_reverse_rate_units": reverse_units,
            "forward_collision_limit_si": collision_limit,
            "forward_collision_limit_units": "m^3/(mol*s)" if collision_limit is not None else None,
            "collision_limit_source": (
                "RMG Reaction.calculate_coll_limit using retained species transport parameters"
                if collision_limit is not None or reverse_collision_limit is not None else None
            ),
            "reverse_collision_limit_si": reverse_collision_limit,
            "reverse_collision_limit_units": "m^3/(mol*s)" if reverse_collision_limit is not None else None,
            "evaluated_at": ({"temperature_K": float(temperature), "pressure_bar": float(pressure_bar)}
                             if temperature is not None and pressure_bar is not None else None),
            "source_library": str(getattr(reaction, "library", "") or "") or None,
            "kinetics_comment": str(getattr(kinetics, "comment", "") or "") or None,
            **_equilibrium_constant_evidence(reaction, temperature),
        })
    result = {"status": "completed", "operation": "inspect", "species": records,
              "reactions": reaction_records}
    if equilibrium_model_artifacts_before is not None:
        artifacts_after = _inspection_model_artifacts(payload)
        unchanged = artifacts_after == equilibrium_model_artifacts_before
        result["equilibrium_model_evidence"] = {
            "schema_version": 1, "kind": "native_rmg_equilibrium_model_evidence",
            "status": "completed" if unchanged else "unavailable",
            "source_artifacts": equilibrium_model_artifacts_before,
            "source_artifacts_unchanged": unchanged,
            "unavailable_reason": None if unchanged else "model_source_artifacts_changed_during_inspection",
            "evaluation_method": "RMG Reaction.get_equilibrium_constant(T, type='Kc')",
            "equilibrium_constant_evaluation_conditions": next(
                (record["equilibrium_constant_evaluated_at"] for record in reaction_records
                 if record["equilibrium_constant_evaluated_at"] is not None), None,
            ),
            "equilibrium_constants_available": sum(
                record["evaluated_equilibrium_constant_Kc"] is not None for record in reaction_records
            ),
            "equilibrium_constants_unavailable": sum(
                record["evaluated_equilibrium_constant_Kc"] is None for record in reaction_records
            ),
            "thermochemistry_basis": "species thermochemistry loaded from the retained Chemkin model",
            "condition_contract_id": payload.get("condition_contract_id"),
            "condition_binding_status": "caller_supplied_reference_only",
            "physical_rate_uncertainty_resolved": False,
            "pressure_applicability_verified": False,
            "conformer_completeness_verified": False,
        }
    if identity_payload is not None:
        result["replacement_identity"] = build_replacement_correspondences(
            identity_payload, species, reactions, reaction_records,
        )
    return result


def _cantera_solution(species, reactions):
    import cantera as ct
    cantera_species = [item.to_cantera(use_chemkin_identifier=True) for item in species]
    # RMG returns one Cantera reaction for most entries, but falloff and a few
    # pressure-dependent reaction types return a sequence.  Passing that
    # nested sequence to Cantera makes an otherwise valid retained mechanism
    # fail to load (notably the dry-air NNO mechanism).
    cantera_reactions = []
    for item in reactions:
        converted = item.to_cantera(species_list=species, use_chemkin_identifier=True)
        if isinstance(converted, (list, tuple)):
            cantera_reactions.extend(converted)
        else:
            cantera_reactions.append(converted)
    try:  # Cantera 3
        return ct.Solution(thermo="ideal-gas", kinetics="gas", species=cantera_species, reactions=cantera_reactions)
    except Exception:  # Cantera 2.6
        return ct.Solution(thermo="IdealGas", kinetics="GasKinetics", species=cantera_species, reactions=cantera_reactions)


def _stoich_coefficient(gas, kind, species_index, reaction_index):
    """Read a Cantera stoichiometric coefficient across supported versions."""
    accessor = getattr(gas, "%s_stoich_coeff" % kind, None)
    if callable(accessor):
        return float(accessor(species_index, reaction_index))
    matrix = getattr(gas, "%s_stoich_coeffs" % kind)
    return float(matrix[species_index, reaction_index])


def _reaction_equations(gas):
    equations = gas.reaction_equations
    return list(equations() if callable(equations) else equations)


def _target_flux_snapshot(gas, reactor, target_index):
    """Return instantaneous gross target destruction/production by reaction.

    Cantera reports rates of progress in kmol m-3 s-1.  Multiplication by the
    retained reactor volume converts each value to kmol s-1.  Forward and
    reverse rates are kept separate so a reversible cycle is not mistaken for
    net target loss.
    """
    forward = list(gas.forward_rates_of_progress)
    reverse = list(gas.reverse_rates_of_progress)
    net = list(gas.net_rates_of_progress)
    volume = float(reactor.volume)
    records = []
    for reaction_index, (forward_rate, reverse_rate, net_rate) in enumerate(
        zip(forward, reverse, net)
    ):
        reactant = _stoich_coefficient(gas, "reactant", target_index, reaction_index)
        product = _stoich_coefficient(gas, "product", target_index, reaction_index)
        change = product - reactant
        forward_extent = max(0.0, float(forward_rate)) * volume
        reverse_extent = max(0.0, float(reverse_rate)) * volume
        destruction = max(0.0, -change) * forward_extent + max(0.0, change) * reverse_extent
        production = max(0.0, change) * forward_extent + max(0.0, -change) * reverse_extent
        # Integrate Cantera's net rate directly.  Computing this later as the
        # difference of two large gross integrals is catastrophically
        # ill-conditioned for fast near-equilibrium reactions.
        net_target_loss = -change * float(net_rate) * volume
        records.append({
            "destruction_kmol_s": destruction,
            "production_kmol_s": production,
            "forward_rate_of_progress_kmol_m3_s": float(forward_rate),
            "reverse_rate_of_progress_kmol_m3_s": float(reverse_rate),
            "net_rate_of_progress_kmol_m3_s": float(net_rate),
            "net_target_loss_kmol_s": net_target_loss,
            "forward_extent_kmol_s": forward_extent,
            "reverse_extent_kmol_s": reverse_extent,
            "target_stoichiometric_change": change,
        })
    return records


def _logarithmic_times(horizon, points, earliest_time):
    """Build a bounded log grid with an explicit fast-chemistry floor."""
    horizon = float(horizon)
    points = max(3, int(points))
    earliest_time = min(horizon, max(float(earliest_time), horizon * 1e-30, 1e-300))
    if earliest_time == horizon:
        return [0.0, horizon]
    start = math.log10(earliest_time)
    stop = math.log10(horizon)
    times = [0.0] + [10 ** (start + (stop - start) * i / float(points - 2))
                     for i in range(points - 1)]
    times[-1] = horizon
    return sorted(set(times))


def _integrate_flux_snapshots(profile, flux_snapshots, reaction_count):
    """Trapezoid-integrate retained rates on one propagation pass."""
    keys = (
        "destruction_kmol_s", "production_kmol_s", "forward_extent_kmol_s",
        "reverse_extent_kmol_s", "net_target_loss_kmol_s",
    )
    totals = {key: [0.0] * reaction_count for key in keys}
    for point_index in range(1, len(profile)):
        delta_time = profile[point_index]["time_seconds"] - profile[point_index - 1]["time_seconds"]
        previous = flux_snapshots[point_index - 1]
        current = flux_snapshots[point_index]
        for reaction_index in range(reaction_count):
            for key in keys:
                totals[key][reaction_index] += 0.5 * delta_time * (
                    previous[reaction_index][key] + current[reaction_index][key]
                )
    return totals


def _flux_closure(
    initial_amount,
    final_amount,
    integrated_destruction,
    integrated_production,
    integrated_net_target_loss=None,
):
    amount_loss = initial_amount - final_amount
    integrated_net_loss = (
        sum(integrated_net_target_loss)
        if integrated_net_target_loss is not None
        else sum(integrated_destruction) - sum(integrated_production)
    )
    absolute_error = abs(integrated_net_loss - amount_loss)
    scale = max(abs(amount_loss), abs(integrated_net_loss), abs(initial_amount) * 1e-12, 1e-300)
    return amount_loss, integrated_net_loss, absolute_error, absolute_error / scale


def propagate(payload):
    """Reintegrate a retained mechanism and retain target-loss attribution."""
    import cantera as ct
    from rmgpy.chemkin import load_chemkin_file, get_species_identifier
    species, reactions = load_chemkin_file(payload["chemkin"], dictionary_path=payload.get("dictionary"))
    multipliers = payload.get("reaction_multipliers") or {}
    for raw_index, multiplier in multipliers.items():
        reaction_index = int(raw_index)
        value = float(multiplier)
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("Reaction multipliers must be finite and nonnegative")
    identifiers = {item.label: _species_identifier(item, get_species_identifier) for item in species}
    requested = payload["initial_mole_fractions"]
    composition = {}
    for label, fraction in requested.items():
        identifier = identifiers.get(label, label)
        composition[identifier] = float(fraction)
    if abs(sum(composition.values()) - 1.0) > 1e-8:
        raise ValueError("Initial mole fractions must sum to one")
    target = identifiers.get(payload["target_label"], payload["target_label"])
    horizon = float(payload["target_duration_seconds"])
    if not math.isfinite(horizon) or horizon <= 0.0:
        raise ValueError("target_duration_seconds must be finite and positive")
    points = max(3, min(int(payload.get("points", 160)), 500))
    max_points = max(points, min(int(payload.get("flux_max_points", 2000)), 5000))
    max_refinements = max(0, min(int(payload.get("flux_refinement_attempts", 3)), 5))
    max_internal_steps = max(
        1000, min(int(payload.get("flux_max_internal_steps", 500000)), 2000000),
    )
    closure_target = float(payload.get("flux_refinement_closure_target", 0.05))
    if not math.isfinite(closure_target) or not 0.0 < closure_target < 1.0:
        raise ValueError("flux_refinement_closure_target must be between zero and one")
    retention = float(payload.get("retention_fraction", 0.95))
    temperature_K = float(payload["temperature_K"])
    pressure_bar = float(payload["pressure_bar"])
    if not all(math.isfinite(value) and value > 0.0 for value in (temperature_K, pressure_bar)):
        raise ValueError("Propagation temperature and pressure must be finite and positive")
    base_solver_rtol = float(payload.get("flux_solver_rtol", 1e-9))
    base_solver_atol = float(payload.get("flux_solver_atol", 1e-15))
    if (not math.isfinite(base_solver_rtol) or base_solver_rtol <= 0.0
            or not math.isfinite(base_solver_atol) or base_solver_atol <= 0.0):
        raise ValueError("Flux solver tolerances must be positive and finite")

    def run_pass(earliest_time, pass_points, solver_rtol, solver_atol):
        gas = _cantera_solution(species, reactions)
        for raw_index, multiplier in multipliers.items():
            gas.set_multiplier(float(multiplier), int(raw_index))
        for label in composition:
            if label not in gas.species_names:
                raise ValueError("Initial species not found in retained mechanism: %s" % label)
        gas.TPX = temperature_K, pressure_bar * 1e5, composition
        try:
            reactor = ct.IdealGasConstPressureReactor(gas, energy="off")
        except TypeError:
            reactor = ct.IdealGasConstPressureReactor(contents=gas, energy="off")
        # Read the actual reactor phase (also correct if a newer Cantera
        # constructor clones the supplied Solution). The stability condition
        # specifies fixed T/P, not adiabatic constant-pressure chemistry.
        gas = reactor.phase if hasattr(reactor, "phase") else reactor.thermo
        observed_temperatures, observed_pressures = [temperature_K], [pressure_bar]

        def check_conditions():
            actual_T, actual_P = float(reactor.T), float(gas.P) / 1e5
            if (bool(reactor.energy_enabled)
                    or not math.isfinite(actual_T) or not math.isfinite(actual_P)
                    or abs(actual_T - temperature_K) > max(1e-8, temperature_K * 1e-10)
                    or abs(actual_P - pressure_bar) > max(1e-12, pressure_bar * 1e-9)):
                raise ValueError("Cantera violated the isothermal constant-pressure condition contract")
            observed_temperatures[:] = [min(observed_temperatures[0], actual_T),
                                        max(observed_temperatures[-1], actual_T)]
            observed_pressures[:] = [min(observed_pressures[0], actual_P),
                                     max(observed_pressures[-1], actual_P)]

        check_conditions()
        network = ct.ReactorNet([reactor])
        network.rtol = solver_rtol
        network.atol = solver_atol
        if target not in gas.species_names:
            raise ValueError("Target species not found in retained mechanism: %s" % payload["target_label"])
        index = gas.species_index(target)

        def amount_kmol():
            return reactor.mass * gas.Y[index] / gas.molecular_weights[index]

        initial = amount_kmol()
        requested_times = _logarithmic_times(horizon, pass_points, earliest_time)
        initial_record = {
            "time_seconds": 0.0,
            "target_amount_kmol": initial,
            "target_fraction_remaining": 1.0 if initial else None,
            "target_mole_fraction": float(gas.X[index]),
        }
        integration_profile = [initial_record]
        pass_snapshots = [_target_flux_snapshot(gas, reactor, index)]
        internal_step_count = 0
        for requested_time in requested_times[1:]:
            while float(network.time) < requested_time:
                remaining = requested_time - float(network.time)
                try:
                    network.max_time_step = remaining
                except Exception:
                    # Cantera versions predating the writable property still
                    # expose adaptive ``step``; exact requested profile points
                    # then become best-effort while the final horizon remains
                    # checked below.
                    pass
                prior_time = float(network.time)
                advanced_time = float(network.step())
                check_conditions()
                if not math.isfinite(advanced_time) or advanced_time <= prior_time:
                    raise RuntimeError("Cantera adaptive step made no finite forward progress")
                internal_step_count += 1
                if internal_step_count > max_internal_steps:
                    raise RuntimeError(
                        "Cantera internal-step flux integration budget exhausted"
                    )
                amount = amount_kmol()
                integration_profile.append({
                    "time_seconds": advanced_time,
                    "target_amount_kmol": amount,
                    "target_fraction_remaining": amount / initial if initial else None,
                    "target_mole_fraction": float(gas.X[index]),
                })
                pass_snapshots.append(_target_flux_snapshot(gas, reactor, index))
        # ``ReactorNet.step`` may advance slightly beyond its requested output
        # boundary on some Cantera versions.  Truncate the final quadrature
        # segment to the immutable horizon instead of silently propagating a
        # longer condition contract.
        if integration_profile[-1]["time_seconds"] > horizon:
            left_record, right_record = integration_profile[-2:]
            left_time = left_record["time_seconds"]
            right_time = right_record["time_seconds"]
            fraction = (horizon - left_time) / (right_time - left_time)

            def interpolate(left, right):
                return left + fraction * (right - left)

            integration_profile[-1] = {
                key: (
                    horizon if key == "time_seconds"
                    else interpolate(float(left_record[key]), float(right_record[key]))
                    if isinstance(left_record.get(key), (int, float))
                    and isinstance(right_record.get(key), (int, float))
                    else right_record.get(key)
                )
                for key in set(left_record) | set(right_record)
            }
            truncated_snapshot = []
            for left_reaction, right_reaction in zip(
                pass_snapshots[-2], pass_snapshots[-1]
            ):
                truncated_snapshot.append({
                    key: interpolate(
                        float(left_reaction[key]), float(right_reaction[key])
                    )
                    for key in left_reaction
                })
            pass_snapshots[-1] = truncated_snapshot

        # Report the requested logarithmic profile exactly.  Internal adaptive
        # steps remain the integration basis; linear interpolation here is
        # presentation/threshold-crossing evidence and never flux quadrature.
        pass_profile = []
        bracket = 1
        for requested_time in requested_times:
            while (
                bracket < len(integration_profile) - 1
                and integration_profile[bracket]["time_seconds"] < requested_time
            ):
                bracket += 1
            right = integration_profile[bracket] if requested_time else integration_profile[0]
            left = integration_profile[max(0, bracket - 1)] if requested_time else right
            if right["time_seconds"] == left["time_seconds"]:
                fraction = 0.0
            else:
                fraction = (
                    (requested_time - left["time_seconds"])
                    / (right["time_seconds"] - left["time_seconds"])
                )
            pass_profile.append({
                "time_seconds": requested_time,
                "target_amount_kmol": (
                    left["target_amount_kmol"]
                    + fraction * (right["target_amount_kmol"] - left["target_amount_kmol"])
                ),
                "target_fraction_remaining": (
                    left["target_fraction_remaining"]
                    + fraction * (
                        right["target_fraction_remaining"]
                        - left["target_fraction_remaining"]
                    )
                ),
                "target_mole_fraction": (
                    left["target_mole_fraction"]
                    + fraction * (
                        right["target_mole_fraction"] - left["target_mole_fraction"]
                    )
                ),
            })
        pass_crossing = None
        for prior, current in zip(pass_profile, pass_profile[1:]):
            if prior["target_fraction_remaining"] >= retention >= current["target_fraction_remaining"]:
                denominator = prior["target_fraction_remaining"] - current["target_fraction_remaining"]
                pass_crossing = current["time_seconds"] if denominator == 0 else prior["time_seconds"] + (
                    (prior["target_fraction_remaining"] - retention) / denominator
                    * (current["time_seconds"] - prior["time_seconds"])
                )
                break
        totals = _integrate_flux_snapshots(
            integration_profile, pass_snapshots, gas.n_reactions,
        )
        closure = _flux_closure(
            initial, integration_profile[-1]["target_amount_kmol"],
            totals["destruction_kmol_s"], totals["production_kmol_s"],
            totals["net_target_loss_kmol_s"],
        )
        return (
            gas, initial, pass_profile, pass_snapshots, pass_crossing, totals,
            closure, internal_step_count, integration_profile[-1]["time_seconds"],
            {"policy_id": "isothermal-ideal-gas-constant-pressure-v1",
             "temperature_K": temperature_K, "pressure_bar": pressure_bar,
             "energy_enabled": False, "condition_validation": "passed",
             "observed_temperature_range_K": observed_temperatures,
             "observed_pressure_range_bar": observed_pressures},
        )

    earliest_time = horizon * 1e-12
    refinement_history = []
    for refinement in range(max_refinements + 1):
        solver_rtol = base_solver_rtol * (0.1 ** refinement)
        solver_atol = base_solver_atol * (0.1 ** refinement)
        (
            gas, initial_amount, profile, flux_snapshots, crossing, totals,
            closure, internal_step_count, simulated_time, reactor_model,
        ) = run_pass(earliest_time, points, solver_rtol, solver_atol)
        amount_loss, integrated_net_loss, closure_absolute_error, closure_error = closure
        crossing_time_resolved = bool(
            crossing is None or profile[1]["time_seconds"] <= crossing * 1e-2
        )
        refinement_history.append({
            "attempt": refinement + 1,
            "sample_count": len(profile),
            "earliest_nonzero_time_seconds": profile[1]["time_seconds"],
            "estimated_time_to_retention_seconds": crossing,
            "numerical_closure_relative_error": closure_error,
            "crossing_time_resolved": crossing_time_resolved,
            "internal_step_count": internal_step_count,
            "solver_rtol": solver_rtol,
            "solver_atol": solver_atol,
            "simulated_time_seconds": simulated_time,
        })
        if (closure_error <= closure_target and crossing_time_resolved) or refinement == max_refinements:
            break
        # A crossing from the coarse pass identifies the missed fast-chemistry
        # scale.  Start four decades before it; repeated failures push the
        # floor three more decades earlier while increasing quadrature density.
        anchor = crossing or profile[1]["time_seconds"]
        earliest_time = min(earliest_time * 1e-3, anchor * 1e-4)
        points = min(max_points, max(points * 2, 320))

    equations = _reaction_equations(gas)
    integrated_destruction = totals["destruction_kmol_s"]
    integrated_production = totals["production_kmol_s"]
    integrated_forward_extent = totals["forward_extent_kmol_s"]
    integrated_reverse_extent = totals["reverse_extent_kmol_s"]
    integrated_net_target_loss = totals["net_target_loss_kmol_s"]
    total_destruction = sum(integrated_destruction)
    total_production = sum(integrated_production)
    positive_net_target_loss = [max(0.0, value) for value in integrated_net_target_loss]
    total_positive_net_target_loss = sum(positive_net_target_loss)
    reaction_flux = []
    for reaction_index, equation in enumerate(equations):
        destruction = integrated_destruction[reaction_index]
        production = integrated_production[reaction_index]
        if (destruction <= 0.0 and production <= 0.0 and
                integrated_forward_extent[reaction_index] <= 0.0 and
                integrated_reverse_extent[reaction_index] <= 0.0):
            continue
        reaction_flux.append({
            "reaction_index": reaction_index,
            "reaction_equation": equation,
            "target_stoichiometric_change": flux_snapshots[-1][reaction_index]["target_stoichiometric_change"],
            "integrated_forward_extent_kmol": integrated_forward_extent[reaction_index],
            "integrated_reverse_extent_kmol": integrated_reverse_extent[reaction_index],
            "integrated_net_extent_kmol": (integrated_forward_extent[reaction_index] -
                                            integrated_reverse_extent[reaction_index]),
            "integrated_target_destruction_kmol": destruction,
            "integrated_target_production_kmol": production,
            "integrated_net_target_loss_kmol": integrated_net_target_loss[reaction_index],
            "positive_net_target_loss_fraction": (
                positive_net_target_loss[reaction_index] / total_positive_net_target_loss
                if total_positive_net_target_loss > 0.0 else None
            ),
            "sampled_rates_of_progress_kmol_m3_s": {
                "final_forward": flux_snapshots[-1][reaction_index]["forward_rate_of_progress_kmol_m3_s"],
                "final_reverse": flux_snapshots[-1][reaction_index]["reverse_rate_of_progress_kmol_m3_s"],
                "final_net": flux_snapshots[-1][reaction_index]["net_rate_of_progress_kmol_m3_s"],
                "maximum_forward": max(item[reaction_index]["forward_rate_of_progress_kmol_m3_s"]
                                       for item in flux_snapshots),
                "maximum_reverse": max(item[reaction_index]["reverse_rate_of_progress_kmol_m3_s"]
                                       for item in flux_snapshots),
                "maximum_absolute_net": max(abs(item[reaction_index]["net_rate_of_progress_kmol_m3_s"])
                                            for item in flux_snapshots),
            },
            "gross_target_destruction_fraction": (
                destruction / total_destruction if total_destruction > 0.0 else None
            ),
        })
    reaction_flux.sort(key=lambda item: item["integrated_target_destruction_kmol"], reverse=True)
    finite_flux = all(math.isfinite(value) for value in (
        total_destruction, total_production, initial_amount, amount_loss, integrated_net_loss,
        closure_absolute_error, closure_error,
    ))
    refinement_converged = closure_error <= closure_target and crossing_time_resolved
    flux_attribution = {
        "status": "completed" if finite_flux and refinement_converged else "incomplete_adaptive_refinement",
        "basis": "cantera_adaptive_internal_step_direct_net_target_rate_integration",
        "rate_basis": "kmol_s-1_after_reactor_volume",
        "sample_count": len(profile),
        "internal_step_count": internal_step_count,
        "adaptive_refinement": {
            "status": (
                "converged"
                if refinement_converged
                else "exhausted_without_closure_or_crossing_resolution"
            ),
            "closure_target": closure_target,
            "attempts": refinement_history,
        },
        "total_integrated_target_destruction_kmol": total_destruction,
        "total_integrated_target_production_kmol": total_production,
        "total_integrated_positive_net_target_loss_kmol": total_positive_net_target_loss,
        "initial_target_amount_kmol": initial_amount,
        "integrated_net_target_loss_kmol": integrated_net_loss,
        "amount_profile_target_loss_kmol": amount_loss,
        "numerical_closure_absolute_error_kmol": closure_absolute_error,
        "numerical_closure_relative_error": closure_error,
        "reactions": reaction_flux,
    }
    return {"status": "completed", "operation": "propagate", "target_chemkin_identifier": target,
            "reactor_model": reactor_model,
            "profile_basis": "target_amount_kmol", "target_profile": profile,
            "target_loss_fraction": max(0.0, 1.0 - profile[-1]["target_fraction_remaining"]),
            "estimated_time_to_retention_seconds": crossing,
            "coverage": {"requested_seconds": horizon, "simulated_seconds": profile[-1]["time_seconds"], "complete": True},
            "reaction_flux_attribution": flux_attribution}


def main():
    try:
        if len(sys.argv) == 3 and sys.argv[1] == "--payload":
            with open(sys.argv[2]) as handle:
                payload = json.load(handle)
        else:
            payload = json.load(sys.stdin)
        operation = payload.get("operation")
        if operation == "inspect":
            result = inspect(payload)
        elif operation == "propagate":
            result = propagate(payload)
        elif operation in ("linear_model_snapshot", "linear_model_diagnostics"):
            try:
                from . import native_linear_diagnostics
            except ImportError:
                import native_linear_diagnostics
            result = getattr(native_linear_diagnostics, operation)(payload)
        elif operation == "trace_family_recipe":
            try:
                from .rmg_recipe_bridge import trace_family_recipe
            except ImportError:
                from rmg_recipe_bridge import trace_family_recipe
            result = trace_family_recipe(payload)
        else:
            raise ValueError("Unsupported bridge operation")
        serialized = json.dumps(result, sort_keys=True, default=str, allow_nan=False)
    except Exception as error:
        result = {"status": "failed", "failure_reason": "%s: %s" % (type(error).__name__, error)}
        serialized = json.dumps(result, sort_keys=True, allow_nan=False)
    sys.stdout.write(serialized)


if __name__ == "__main__":
    main()
