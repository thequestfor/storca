#!/usr/bin/env python
"""RMG-environment side of STORCA's JSON artifact bridge.

This file intentionally has no STORCA imports: it is executed by RMG 3/Python
3.9 and RMG 4/Python 3.11 as well as the app's Python 3.12 environment.
"""

from __future__ import print_function

import json
import math
import sys


def _species_identifier(species, get_species_identifier):
    try:
        return get_species_identifier(species)
    except Exception:
        return species.label


def inspect(payload):
    from rmgpy.chemkin import load_chemkin_file, get_species_identifier
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
        molecularity = len(reaction.reactants)
        reverse_molecularity = len(reaction.products)
        rate = None
        reverse_rate = None
        collision_limit = None
        reverse_collision_limit = None
        if temperature is not None:
            try:
                rate = reaction.get_rate_coefficient(float(temperature), float(pressure_bar) * 1e5)
            except TypeError:
                rate = reaction.get_rate_coefficient(float(temperature))
            if bool(reaction.reversible):
                try:
                    reverse_kinetics = reaction.generate_reverse_rate_coefficient()
                    try:
                        reverse_rate = reverse_kinetics.get_rate_coefficient(
                            float(temperature), float(pressure_bar) * 1e5,
                        )
                    except TypeError:
                        reverse_rate = reverse_kinetics.get_rate_coefficient(float(temperature))
                except Exception:
                    reverse_rate = None
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
        units = {1: "s^-1", 2: "m^3/(mol*s)", 3: "m^6/(mol^2*s)"}.get(molecularity)
        reverse_units = {1: "s^-1", 2: "m^3/(mol*s)", 3: "m^6/(mol^2*s)"}.get(reverse_molecularity)
        kinetics = getattr(reaction, "kinetics", None)
        reactant_labels = [_species_identifier(item, get_species_identifier) for item in reaction.reactants]
        product_labels = [_species_identifier(item, get_species_identifier) for item in reaction.products]
        arrow = "<=>" if reaction.reversible else "=>"
        equation = "+".join(reactant_labels) + arrow + "+".join(product_labels)

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
        })
    return {"status": "completed", "operation": "inspect", "species": records,
            "reactions": reaction_records}


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
    volume = float(reactor.volume)
    records = []
    for reaction_index, (forward_rate, reverse_rate) in enumerate(zip(forward, reverse)):
        reactant = _stoich_coefficient(gas, "reactant", target_index, reaction_index)
        product = _stoich_coefficient(gas, "product", target_index, reaction_index)
        change = product - reactant
        forward_extent = max(0.0, float(forward_rate)) * volume
        reverse_extent = max(0.0, float(reverse_rate)) * volume
        destruction = max(0.0, -change) * forward_extent + max(0.0, change) * reverse_extent
        production = max(0.0, change) * forward_extent + max(0.0, -change) * reverse_extent
        records.append({
            "destruction_kmol_s": destruction,
            "production_kmol_s": production,
            "forward_rate_of_progress_kmol_m3_s": float(forward_rate),
            "reverse_rate_of_progress_kmol_m3_s": float(reverse_rate),
            "net_rate_of_progress_kmol_m3_s": float(forward_rate) - float(reverse_rate),
            "forward_extent_kmol_s": forward_extent,
            "reverse_extent_kmol_s": reverse_extent,
            "target_stoichiometric_change": change,
        })
    return records


def propagate(payload):
    """Reintegrate a retained mechanism and retain target-loss attribution."""
    import cantera as ct
    from rmgpy.chemkin import load_chemkin_file, get_species_identifier
    species, reactions = load_chemkin_file(payload["chemkin"], dictionary_path=payload.get("dictionary"))
    gas = _cantera_solution(species, reactions)
    for raw_index, multiplier in (payload.get("reaction_multipliers") or {}).items():
        reaction_index = int(raw_index)
        value = float(multiplier)
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("Reaction multipliers must be finite and nonnegative")
        gas.set_multiplier(value, reaction_index)
    identifiers = {item.label: _species_identifier(item, get_species_identifier) for item in species}
    requested = payload["initial_mole_fractions"]
    composition = {}
    for label, fraction in requested.items():
        identifier = identifiers.get(label, label)
        if identifier not in gas.species_names:
            raise ValueError("Initial species not found in retained mechanism: %s" % label)
        composition[identifier] = float(fraction)
    if abs(sum(composition.values()) - 1.0) > 1e-8:
        raise ValueError("Initial mole fractions must sum to one")
    gas.TPX = float(payload["temperature_K"]), float(payload["pressure_bar"]) * 1e5, composition
    try:
        reactor = ct.IdealGasConstPressureReactor(gas)
    except TypeError:
        reactor = ct.IdealGasConstPressureReactor(contents=gas)
    network = ct.ReactorNet([reactor])
    target = identifiers.get(payload["target_label"], payload["target_label"])
    if target not in gas.species_names:
        raise ValueError("Target species not found in retained mechanism: %s" % payload["target_label"])
    index = gas.species_index(target)
    horizon = float(payload["target_duration_seconds"])
    points = max(2, min(int(payload.get("points", 160)), 500))
    times = [0.0] + [horizon * (10 ** (math.log10(1e-12) + (12.0 * i / float(points - 2)))) for i in range(points - 1)]
    # Ensure exact horizon and remove duplicate underflow values.
    times[-1] = horizon
    times = sorted(set(times))

    def amount_kmol():
        return reactor.mass * gas.Y[index] / gas.molecular_weights[index]

    initial_amount = amount_kmol()
    profile = []
    flux_snapshots = []
    for time in times:
        if time:
            network.advance(time)
        amount = amount_kmol()
        profile.append({"time_seconds": time, "target_amount_kmol": amount,
                        "target_fraction_remaining": amount / initial_amount if initial_amount else None,
                        "target_mole_fraction": float(gas.X[index])})
        flux_snapshots.append(_target_flux_snapshot(gas, reactor, index))
    retention = float(payload.get("retention_fraction", 0.95))
    crossing = None
    for prior, current in zip(profile, profile[1:]):
        if prior["target_fraction_remaining"] >= retention >= current["target_fraction_remaining"]:
            denominator = prior["target_fraction_remaining"] - current["target_fraction_remaining"]
            crossing = current["time_seconds"] if denominator == 0 else prior["time_seconds"] + (
                (prior["target_fraction_remaining"] - retention) / denominator * (current["time_seconds"] - prior["time_seconds"])
            )
            break
    equations = _reaction_equations(gas)
    integrated_destruction = [0.0] * gas.n_reactions
    integrated_production = [0.0] * gas.n_reactions
    integrated_forward_extent = [0.0] * gas.n_reactions
    integrated_reverse_extent = [0.0] * gas.n_reactions
    for point_index in range(1, len(profile)):
        delta_time = profile[point_index]["time_seconds"] - profile[point_index - 1]["time_seconds"]
        previous = flux_snapshots[point_index - 1]
        current = flux_snapshots[point_index]
        for reaction_index in range(gas.n_reactions):
            integrated_destruction[reaction_index] += 0.5 * delta_time * (
                previous[reaction_index]["destruction_kmol_s"] + current[reaction_index]["destruction_kmol_s"]
            )
            integrated_production[reaction_index] += 0.5 * delta_time * (
                previous[reaction_index]["production_kmol_s"] + current[reaction_index]["production_kmol_s"]
            )
            integrated_forward_extent[reaction_index] += 0.5 * delta_time * (
                previous[reaction_index]["forward_extent_kmol_s"] + current[reaction_index]["forward_extent_kmol_s"]
            )
            integrated_reverse_extent[reaction_index] += 0.5 * delta_time * (
                previous[reaction_index]["reverse_extent_kmol_s"] + current[reaction_index]["reverse_extent_kmol_s"]
            )
    total_destruction = sum(integrated_destruction)
    total_production = sum(integrated_production)
    positive_net_target_loss = [max(0.0, destruction - production)
                                for destruction, production in zip(integrated_destruction, integrated_production)]
    total_positive_net_target_loss = sum(positive_net_target_loss)
    amount_loss = initial_amount - profile[-1]["target_amount_kmol"]
    integrated_net_loss = total_destruction - total_production
    # A relative closure metric alone is misleading when a reversible system
    # has essentially zero net conversion: two tiny, nearly cancelling values
    # can have a large relative difference.  Retain both metrics so the policy
    # layer can distinguish numerical dust from material target loss without
    # discarding the underlying gross turnover evidence.
    closure_absolute_error = abs(integrated_net_loss - amount_loss)
    closure_scale = max(abs(amount_loss), abs(integrated_net_loss), abs(initial_amount) * 1e-12, 1e-300)
    closure_error = closure_absolute_error / closure_scale
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
            "integrated_net_target_loss_kmol": destruction - production,
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
    flux_attribution = {
        "status": "completed" if all(math.isfinite(value) for value in (
            total_destruction, total_production, initial_amount, amount_loss, integrated_net_loss,
            closure_absolute_error, closure_error
        )) else "invalid",
        "basis": "time_integrated_cantera_rates_of_progress",
        "rate_basis": "kmol_s-1_after_reactor_volume",
        "sample_count": len(profile),
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
        else:
            raise ValueError("Unsupported bridge operation")
    except Exception as error:
        result = {"status": "failed", "failure_reason": "%s: %s" % (type(error).__name__, error)}
    sys.stdout.write(json.dumps(result, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
