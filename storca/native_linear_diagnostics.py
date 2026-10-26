"""Read-only fixed first-order model diagnostics in the owning RMG environment.

Imports are deliberately lazy: neither RMG, Cantera, NumPy nor SciPy becomes an
application dependency.  These are declared-model comparisons, not physical
uncertainty bounds, mechanism-completeness evidence, or lifetime certificates.
"""

import hashlib
import math
from collections import Counter
from pathlib import Path


POLICY_ID = "native-retained-linear-diagnostics-v1"
MAX_SPECIES = 32
MAX_CHANNELS = 128
MAX_SCENARIOS = 16
PARENT_AGREEMENT = 1e-6


def _number(value, name, positive=False):
    if isinstance(value, (bool, str, bytes)):
        raise ValueError("%s must be numerical" % name)
    number = float(value)
    if not math.isfinite(number) or number < 0.0 or (positive and number == 0.0):
        raise ValueError("%s must be finite and %s" % (name, "positive" if positive else "nonnegative"))
    return number


def _digest(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _exports(payload):
    result = []
    for role in ("chemkin", "dictionary", "transport"):
        value = payload.get(role)
        if value is None and role == "transport":
            continue
        if not value:
            raise ValueError("native_linear_%s_source_missing" % role)
        path = Path(value).resolve(strict=True)
        if not path.is_file():
            raise ValueError("native_linear_source_not_file")
        result.append({"path": str(path), "sha256": _digest(path)})
    return result


def _validate_sources(records):
    if not isinstance(records, list) or not records:
        raise ValueError("native_linear_source_artifacts_missing")
    seen = set()
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("native_linear_source_artifact_invalid")
        path = Path(record.get("path") or "")
        if not path.is_absolute() or str(path) in seen or not path.is_file():
            raise ValueError("native_linear_source_artifact_path_invalid")
        seen.add(str(path))
        if _digest(path) != record.get("sha256"):
            raise ValueError("native_linear_source_artifact_hash_mismatch")


def _bridge():
    try:
        from . import rmg_bridge
    except ImportError:
        import rmg_bridge
    return rmg_bridge


def _state_record(species, identifier):
    if not species.molecule or not species.molecule[0].atoms:
        raise ValueError("native_linear_species_graph_missing")
    molecule = species.molecule[0]
    multiplicity = species.multiplicity
    if (isinstance(multiplicity, bool) or not isinstance(multiplicity, int)
            or multiplicity < 1 or molecule.multiplicity != multiplicity):
        raise ValueError("native_linear_species_electronic_state_invalid")
    composition, isotopes = Counter(), Counter()
    charge = 0
    for atom in molecule.atoms:
        symbol, isotope = atom.element.symbol, atom.element.isotope
        if (symbol in ("X", "R", "L", "e") or isinstance(isotope, bool)
                or not isinstance(isotope, int)):
            raise ValueError("native_linear_species_state_outside_policy")
        if isinstance(atom.charge, bool) or int(atom.charge) != atom.charge:
            raise ValueError("native_linear_species_charge_invalid")
        composition[symbol] += 1
        isotopes["%s:%d" % (symbol, isotope)] += 1
        charge += int(atom.charge)
    return {
        "chemkin_identifier": identifier, "label": species.label,
        "multiplicity": multiplicity, "net_charge": charge,
        "adjacency_list": molecule.to_adjacency_list(),
        "atom_composition": dict(sorted(composition.items())),
        "isotopic_composition": dict(sorted(isotopes.items())),
    }


def linear_model_snapshot(payload):
    """Retain operative Cantera constants and exact native graph/state IDs."""
    from rmgpy.chemkin import load_chemkin_file, get_species_identifier
    bridge = _bridge()
    before = _exports(payload)
    temperature = _number(payload.get("temperature_K"), "temperature_K", True)
    pressure = _number(payload.get("pressure_bar"), "pressure_bar", True)
    species, reactions = load_chemkin_file(
        payload["chemkin"], dictionary_path=payload["dictionary"],
        transport_path=payload.get("transport"),
    )
    if not 1 <= len(species) <= MAX_SPECIES or len(reactions) > MAX_CHANNELS:
        raise ValueError("native_linear_model_size_outside_policy")
    ids = [bridge._species_identifier(item, get_species_identifier) for item in species]
    if any(not isinstance(value, str) or not value for value in ids) or len(set(ids)) != len(ids):
        raise ValueError("native_linear_species_identifier_invalid_or_duplicate")
    records = [_state_record(item, identifier) for item, identifier in zip(species, ids)]
    record_by_id = {item["chemkin_identifier"]: item for item in records}
    initial = payload.get("initial_fractions", payload.get("initial_mole_fractions"))
    if not isinstance(initial, dict) or not initial or set(initial) - set(ids):
        raise ValueError("native_linear_initial_inventory_invalid")
    initial = {identifier: _number(initial.get(identifier, 0.0), "initial_fraction") for identifier in ids}
    if abs(sum(initial.values()) - 1.0) > 1e-12:
        raise ValueError("native_linear_initial_inventory_not_normalized")
    parent = payload.get("parent_id", payload.get("target_label"))
    if parent not in initial or initial[parent] <= 0.0:
        raise ValueError("native_linear_parent_identifier_invalid")
    # Flattening is correct for general propagation but is deliberately not
    # admitted here: native reaction indices must map one-to-one to Cantera.
    for reaction in reactions:
        converted = reaction.to_cantera(species_list=species, use_chemkin_identifier=True)
        if isinstance(converted, (tuple, list)) and len(converted) != 1:
            raise ValueError("native_linear_native_cantera_reaction_mapping_not_one_to_one")
    gas = bridge._cantera_solution(species, reactions)
    if gas.n_reactions != len(reactions) or list(gas.species_names) != ids:
        raise ValueError("native_linear_native_cantera_inventory_mapping_invalid")
    gas.TPX = temperature, pressure * 1e5, initial
    operative_forward = [float(value) for value in gas.forward_rate_constants]
    operative_reverse = [float(value) for value in gas.reverse_rate_constants]
    # A one-to-one graph reaction is not necessarily a linear kinetic law:
    # species-dependent collider efficiencies can make the operative constant
    # composition dependent even at fixed T/P.  Challenge every pure-species
    # composition before admitting a constant generator.  No approximation of
    # a composition-dependent law as a fixed rate is made here.
    for identifier in ids:
        gas.TPX = temperature, pressure * 1e5, {identifier: 1.0}
        for retained, tested in zip(operative_forward + operative_reverse,
                                   list(gas.forward_rate_constants) + list(gas.reverse_rate_constants)):
            tested = _number(float(tested), "operative_composition_challenge_rate")
            if not math.isclose(retained, tested, rel_tol=1e-12, abs_tol=0.0):
                raise ValueError("native_linear_composition_dependent_rate_law")
    gas.TPX = temperature, pressure * 1e5, initial
    channels = []
    for index, reaction in enumerate(reactions):
        if len(reaction.reactants) != 1 or len(reaction.products) != 1:
            raise ValueError("native_linear_non_first_order_stoichiometry")
        left = bridge._species_identifier(reaction.reactants[0], get_species_identifier)
        right = bridge._species_identifier(reaction.products[0], get_species_identifier)
        if left == right or left not in record_by_id or right not in record_by_id:
            raise ValueError("native_linear_channel_endpoint_invalid")
        for key in ("atom_composition", "isotopic_composition", "net_charge"):
            if record_by_id[left][key] != record_by_id[right][key]:
                raise ValueError("native_linear_channel_atom_isotope_or_charge_unbalanced")
        actual = gas.reaction(index)
        if (dict(actual.reactants) != {left: 1.0} or dict(actual.products) != {right: 1.0}
                or bool(actual.reversible) != bool(reaction.reversible)
                or getattr(actual, "orders", {})):
            raise ValueError("native_linear_native_cantera_channel_mapping_invalid")
        # Implicit third bodies are absent from the displayed stoichiometry;
        # forward_rate_constants can exclude their concentration factor.  Do
        # not mistake that coefficient for an operative first-order rate.
        reaction_type = str(getattr(actual, "reaction_type", "")).lower()
        if (getattr(actual, "third_body", None) is not None
                or any(marker in reaction_type for marker in ("three-body", "falloff", "chemically-activated"))):
            raise ValueError("native_linear_implicit_collider_rate_law_unsupported")
        cross_check = bridge._equilibrium_constant_evidence(reaction, temperature)
        forward = _number(operative_forward[index], "operative_forward_rate")
        reverse = _number(operative_reverse[index], "operative_reverse_rate")
        try:
            native_forward = _number(float(reaction.get_rate_coefficient(temperature, pressure * 1e5)),
                                     "native_forward_cross_check")
        except (AttributeError, ValueError, TypeError, OverflowError):
            native_forward = None
        cross_check["native_forward_rate_s^-1"] = native_forward
        cross_check["cantera_operative_forward_rate_s^-1"] = forward
        cross_check["operative_native_forward_relative_difference"] = (
            abs(forward - native_forward) / native_forward if native_forward else None
        )
        cross_check["conversion_difference_in_reverse_omission_bound"] = False
        cross_check["cantera_operative_forward_reverse_ratio"] = (
            forward / reverse if reverse > 0.0 and math.isfinite(forward / reverse) else None
        )
        directions = [("forward", left, right, forward)]
        if reaction.reversible:
            directions.append(("reverse", right, left, reverse))
        elif reverse != 0.0:
            raise ValueError("native_linear_irreversible_channel_has_operative_reverse")
        for direction, reactant, product, rate in directions:
            channels.append({
                "channel_id": "native-%d-%s" % (index, direction),
                "reaction_index": index, "native_reaction_index": index,
                "direction": direction, "reactants": {reactant: 1}, "products": {product: 1},
                "molecularity": 1, "rate_s^-1": rate, "rate_units": "s^-1",
                "reversible": bool(reaction.reversible),
                "native_equilibrium_cross_check": cross_check,
            })
    if len(channels) > MAX_CHANNELS:
        raise ValueError("native_linear_channel_count_outside_policy")
    if _exports(payload) != before:
        raise ValueError("native_linear_sources_changed_during_snapshot")
    result = {
        "schema_version": 1, "status": "completed",
        "kind": "native_fixed_first_order_rate_snapshot", "policy_id": POLICY_ID,
        "condition": {"temperature_K": temperature, "pressure_bar": pressure},
        "condition_contract_id": payload.get("condition_contract_id"),
        "species": records, "channels": channels, "source_artifacts": before,
        "initial_fractions": initial, "parent_id": parent,
        "rate_basis": "cantera_operative_rate_constants",
        "fixed_condition_rate_composition_invariance_checked": True,
        "source_artifacts_unchanged": True,
        "physical_rate_bounds_established": False, "uncertainty_calibrated": False,
        "terminal_robustness_eligible": False, "lifetime_claim_allowed": False,
    }
    if payload.get("include_rate_source_snapshot") is True:
        try:
            from .native_rate_sources import collect_rate_source_snapshot
        except ImportError:
            from native_rate_sources import collect_rate_source_snapshot
        result["rate_source_snapshot"] = collect_rate_source_snapshot(payload, species, reactions, records)
    return result


def _validated_inventory(inventory, parent):
    ids = inventory.get("species_ids")
    if not isinstance(ids, list) or not 1 <= len(ids) <= MAX_SPECIES or len(set(ids)) != len(ids):
        raise ValueError("native_linear_species_ids_invalid")
    if any(not isinstance(identifier, str) or not identifier for identifier in ids):
        raise ValueError("native_linear_species_ids_invalid")
    initial = inventory.get("initial_fractions")
    if not isinstance(initial, dict) or set(initial) - set(ids):
        raise ValueError("native_linear_initial_fractions_invalid")
    initial = [_number(initial.get(identifier, 0.0), "initial_fraction") for identifier in ids]
    if abs(sum(initial) - 1.0) > 1e-12 or parent not in ids or initial[ids.index(parent)] <= 0.0:
        raise ValueError("native_linear_initial_fractions_or_parent_invalid")
    channels = inventory.get("channels")
    if not isinstance(channels, list) or len(channels) > MAX_CHANNELS:
        raise ValueError("native_linear_channels_invalid")
    seen, directed = set(), []
    for channel in channels:
        identifier = channel.get("channel_id")
        left, right = channel.get("reactants"), channel.get("products")
        if (not isinstance(identifier, str) or not identifier or identifier in seen
                or not isinstance(left, dict) or not isinstance(right, dict)
                or len(left) != 1 or len(right) != 1 or channel.get("molecularity") != 1
                or isinstance(channel.get("molecularity"), bool)):
            raise ValueError("native_linear_channel_invalid")
        source, destination = next(iter(left)), next(iter(right))
        if (source not in ids or destination not in ids or source == destination
                or left[source] != 1 or right[destination] != 1
                or isinstance(left[source], bool) or isinstance(right[destination], bool)):
            raise ValueError("native_linear_channel_not_conservative_first_order")
        directed.append((identifier, ids.index(source), ids.index(destination),
                         _number(channel.get("rate_s^-1"), "channel_rate")))
        seen.add(identifier)
    return ids, initial, directed


def _times(horizon, channels, requested=None):
    if requested is not None:
        if not isinstance(requested, list) or not 2 <= len(requested) <= 500:
            raise ValueError("native_linear_retained_grid_invalid")
        values = [_number(value, "retained_grid_time") for value in requested]
        if values[0] != 0.0 or values[-1] != horizon or any(right <= left for left, right in zip(values, values[1:])):
            raise ValueError("native_linear_retained_grid_must_include_exact_horizon")
        return values
    fastest = max([channel[3] for channel in channels] or [0.0])
    earliest = min(horizon * 1e-12, 1e-3 / fastest if fastest else horizon * 1e-12)
    return _bridge()._logarithmic_times(horizon, 120, max(earliest, horizon * 1e-30))


def _solve_pass(ids, initial, channels, multipliers, times, method, refined):
    import numpy as np
    from scipy.integrate import solve_ivp
    count, channel_count = len(ids), len(channels)
    jacobian = np.zeros((count + channel_count, count + channel_count))
    incidence = np.zeros((count, channel_count))
    for position, (identifier, source, destination, rate) in enumerate(channels):
        rate *= multipliers.get(identifier, 1.0)
        if not math.isfinite(rate):
            raise ValueError("native_linear_scenario_rate_overflow")
        jacobian[destination, source] += rate
        jacobian[source, source] -= rate
        jacobian[count + position, source] = rate
        incidence[source, position], incidence[destination, position] = -1.0, 1.0
    if not np.isfinite(jacobian).all():
        raise ValueError("native_linear_scenario_generator_nonfinite")
    y0 = np.array(initial + [0.0] * channel_count)
    atol = np.array([1e-26 if refined else 1e-24] * count + [1e-18 if refined else 1e-16] * channel_count)
    result = solve_ivp(
        lambda time, populations: jacobian.dot(populations), (0.0, times[-1]), y0,
        method=method, jac=jacobian, t_eval=times,
        rtol=1e-11 if refined else 1e-9, atol=atol,
    )
    if not result.success or result.y.shape != (len(y0), len(times)) or result.nfev > 100000:
        raise RuntimeError("native_linear_solver_failed_or_budget_exhausted")
    if not np.isfinite(result.y).all():
        raise RuntimeError("native_linear_solver_nonfinite")
    populations, flux = result.y[:count], result.y[count:]
    conservation = float(np.max(np.abs(populations.sum(axis=0) - 1.0)))
    balance = float(np.max(np.abs(populations - np.array(initial)[:, None] - incidence.dot(flux))))
    minimum = float(np.min(populations))
    minimum_flux = float(np.min(flux)) if channel_count else 0.0
    if conservation > 1e-8 or balance > 1e-8 or minimum < -1e-22 or minimum_flux < -1e-16:
        raise RuntimeError("native_linear_solver_conservation_positivity_or_flux_balance_failed")
    return populations, flux, {
        "method": method, "rtol": 1e-11 if refined else 1e-9,
        "population_atol": 1e-26 if refined else 1e-24,
        "integrated_flux_atol": 1e-18 if refined else 1e-16,
        "function_evaluations": int(result.nfev), "jacobian_evaluations": int(result.njev),
        "minimum_population": minimum, "maximum_conservation_error": conservation,
        "maximum_integrated_flux_balance_error": balance,
        "exact_horizon_reached": float(result.t[-1]) == times[-1],
        "clipping_renormalization_or_qss_used": False,
    }


def solve_fixed_linear_models(inventory, parent, horizon, scenarios, retained_grid_times=None):
    """Small deterministic solver, also usable for analytical unit tests.

    Source validation and native baseline agreement belong to the bridge
    operation; calling this pure solver alone cannot certify a native model.
    """
    import numpy as np
    ids, initial, channels = _validated_inventory(inventory, parent)
    horizon = _number(horizon, "target_duration_seconds", True)
    times = _times(horizon, channels, retained_grid_times)
    if not isinstance(scenarios, list) or not 1 <= len(scenarios) <= MAX_SCENARIOS:
        raise ValueError("native_linear_scenario_count_invalid")
    results, seen = [], set()
    known = {item[0] for item in channels}
    for scenario in scenarios:
        identifier = scenario.get("scenario_id")
        multipliers = scenario.get("multipliers", {})
        if (not isinstance(identifier, str) or not identifier or identifier in seen
                or not isinstance(multipliers, dict) or set(multipliers) - known):
            raise ValueError("native_linear_scenario_invalid")
        multipliers = {key: _number(value, "directed_channel_multiplier") for key, value in multipliers.items()}
        seen.add(identifier)
        coarse_x, coarse_flux, coarse_checks = _solve_pass(ids, initial, channels, multipliers, times, "BDF", False)
        fine_x, fine_flux, fine_checks = _solve_pass(ids, initial, channels, multipliers, times, "BDF", True)
        independent_x, independent_flux, independent_checks = _solve_pass(ids, initial, channels, multipliers, times, "Radau", True)
        population_refinement = float(np.max(np.abs(fine_x - coarse_x)))
        population_independent = float(np.max(np.abs(fine_x - independent_x)))
        flux_scale = np.maximum(np.abs(fine_flux), np.abs(independent_flux))
        flux_refinement = np.abs(fine_flux - coarse_flux)
        flux_independent = np.abs(fine_flux - independent_flux)
        tiny_mask = (np.abs(fine_x) > 1e-24) & (np.abs(fine_x) <= 1e-10)
        tiny_scale = np.maximum(np.abs(fine_x), 1e-24)
        tiny_refinement = np.abs(fine_x - coarse_x)
        tiny_independent = np.abs(fine_x - independent_x)
        if (population_refinement > PARENT_AGREEMENT or population_independent > PARENT_AGREEMENT
                or np.any(flux_refinement > 1e-9 + 1e-6 * flux_scale)
                or np.any(flux_independent > 1e-9 + 1e-6 * flux_scale)
                or np.any(tiny_refinement[tiny_mask] > (1e-24 + 1e-5 * tiny_scale)[tiny_mask])
                or np.any(tiny_independent[tiny_mask] > (1e-24 + 1e-5 * tiny_scale)[tiny_mask])):
            raise RuntimeError("native_linear_solver_tolerance_or_independent_agreement_failed")
        parent_index = ids.index(parent)
        profile = [{
            "time_seconds": time,
            "target_fraction_remaining": float(fine_x[parent_index, index] / initial[parent_index]),
            "species_fractions": {species_id: float(fine_x[position, index]) for position, species_id in enumerate(ids)},
        } for index, time in enumerate(times)]
        results.append({
            "scenario_id": identifier, "multipliers": multipliers,
            "horizon_parent_fraction": profile[-1]["target_fraction_remaining"],
            "target_profile": profile,
            "integrated_channel_flux": {channel[0]: float(fine_flux[index, -1]) for index, channel in enumerate(channels)},
            "numerical_checks": {
                "status": "passed", "passes": [coarse_checks, fine_checks, independent_checks],
                "maximum_population_refinement_difference": population_refinement,
                "maximum_population_independent_difference": population_independent,
                "maximum_integrated_flux_refinement_difference": float(np.max(flux_refinement)) if channels else 0.0,
                "maximum_integrated_flux_independent_difference": float(np.max(flux_independent)) if channels else 0.0,
                "population_absolute_agreement_tolerance": PARENT_AGREEMENT,
                "flux_absolute_agreement_tolerance": 1e-9,
                "flux_relative_agreement_tolerance": 1e-6,
                "tiny_population_resolution_checked": True,
                "tiny_population_absolute_agreement_tolerance": 1e-24,
                "tiny_population_relative_agreement_tolerance": 1e-5,
                "maximum_tiny_population_refinement_relative_difference": (
                    float(np.max((tiny_refinement / tiny_scale)[tiny_mask])) if np.any(tiny_mask) else 0.0
                ),
                "maximum_tiny_population_independent_relative_difference": (
                    float(np.max((tiny_independent / tiny_scale)[tiny_mask])) if np.any(tiny_mask) else 0.0
                ),
                "native_baseline_agreement_checked": False,
            },
        })
    return results


def _cantera_baseline(payload, times):
    import cantera as ct
    from rmgpy.chemkin import load_chemkin_file
    species, reactions = load_chemkin_file(
        payload["chemkin"], dictionary_path=payload["dictionary"], transport_path=payload.get("transport"),
    )
    profiles, checks = [], []
    for rtol, atol in ((1e-10, 1e-24), (1e-12, 1e-26)):
        gas = _bridge()._cantera_solution(species, reactions)
        gas.TPX = payload["temperature_K"], payload["pressure_bar"] * 1e5, payload["initial_fractions"]
        try:
            reactor = ct.IdealGasConstPressureReactor(gas, energy="off")
        except TypeError:
            reactor = ct.IdealGasConstPressureReactor(contents=gas, energy="off")
        gas = reactor.phase if hasattr(reactor, "phase") else reactor.thermo
        network = ct.ReactorNet([reactor])
        network.rtol, network.atol = rtol, atol
        parent = gas.species_index(payload["parent_id"])
        initial_amount = reactor.mass * gas.Y[parent] / gas.molecular_weights[parent]
        profile = []
        for time in times:
            if time:
                network.advance(time)
            if (reactor.energy_enabled or abs(float(reactor.T) - payload["temperature_K"]) > 1e-7
                    or abs(float(gas.P) / 1e5 - payload["pressure_bar"]) > 1e-10):
                raise RuntimeError("native_linear_cantera_baseline_condition_contract_violated")
            value = reactor.mass * gas.Y[parent] / gas.molecular_weights[parent] / initial_amount
            if not math.isfinite(value) or value < -1e-22 or abs(float(network.time) - time) > max(1e-12, time * 1e-12):
                raise RuntimeError("native_linear_cantera_baseline_nonfinite_negative_or_horizon_invalid")
            profile.append(value)
        profiles.append(profile)
        checks.append({"rtol": rtol, "atol": atol, "energy_enabled": False, "exact_requested_times": True})
    difference = max(abs(left - right) for left, right in zip(*profiles))
    if difference > PARENT_AGREEMENT:
        raise RuntimeError("native_linear_cantera_baseline_refinement_failed")
    return profiles[-1], {"passes": checks, "maximum_refinement_difference": difference}


def linear_model_diagnostics(payload):
    """Reopen every input and compare a fixed native baseline and scenarios."""
    result = {
        "schema_version": 1, "kind": "native_linear_model_diagnostics", "policy_id": POLICY_ID,
        "status": "numerical_unresolved", "scenarios": [],
        "physical_rate_bounds_established": False, "uncertainty_calibrated": False,
        "terminal_robustness_eligible": False, "lifetime_claim_allowed": False,
    }
    try:
        sources = payload.get("source_artifacts")
        _validate_sources(sources)
        condition = payload.get("condition") or {}
        native = dict(payload)
        native["temperature_K"] = payload.get("temperature_K", condition.get("temperature_K"))
        native["pressure_bar"] = payload.get("pressure_bar", condition.get("pressure_bar"))
        if condition and any(condition.get(key) != native[key] for key in ("temperature_K", "pressure_bar")):
            raise ValueError("native_linear_exact_condition_mismatch")
        inventory = payload.get("inventory")
        if not isinstance(inventory, dict):
            raise ValueError("native_linear_inventory_missing")
        native["initial_fractions"] = inventory.get("initial_fractions")
        native["parent_id"] = payload.get("parent_id", inventory.get("parent_id"))
        snapshot = linear_model_snapshot(native)
        if inventory.get("condition") is not None and inventory["condition"] != snapshot["condition"]:
            raise ValueError("native_linear_inventory_condition_mismatch")
        records_by_path = {record["path"]: record["sha256"] for record in sources}
        if any(records_by_path.get(record["path"]) != record["sha256"] for record in snapshot["source_artifacts"]):
            raise ValueError("native_linear_snapshot_exports_not_bound")
        ids = [item["chemkin_identifier"] for item in snapshot["species"]]
        if inventory.get("species_ids") != ids or inventory.get("initial_fractions") != snapshot["initial_fractions"]:
            raise ValueError("native_linear_snapshot_inventory_mismatch")
        if inventory.get("species") is not None:
            state_fields = ("chemkin_identifier", "multiplicity", "net_charge", "adjacency_list",
                            "atom_composition", "isotopic_composition")
            if ([{key: record.get(key) for key in state_fields} for record in inventory["species"]]
                    != [{key: record.get(key) for key in state_fields} for record in snapshot["species"]]):
                raise ValueError("native_linear_snapshot_species_graph_or_state_mismatch")
        fields = ("channel_id", "reaction_index", "direction", "reactants", "products", "molecularity", "rate_s^-1")
        supplied = inventory.get("channels", [])
        if (len(supplied) != len(snapshot["channels"])
                or [{key: record.get(key) for key in fields} for record in supplied]
                != [{key: record.get(key) for key in fields} for record in snapshot["channels"]]):
            raise ValueError("native_linear_snapshot_operative_channel_mismatch")
        scenarios = payload.get("scenarios")
        if not isinstance(scenarios, list) or not scenarios or scenarios[0].get("multipliers", {}) != {}:
            raise ValueError("native_linear_first_scenario_must_be_unmodified_baseline")
        horizon = _number(payload.get("target_duration_seconds"), "target_duration_seconds", True)
        solved = solve_fixed_linear_models(
            inventory, native["parent_id"], horizon, scenarios, payload.get("retained_grid_times"),
        )
        times = [item["time_seconds"] for item in solved[0]["target_profile"]]
        cantera_profile, cantera_checks = _cantera_baseline(native, times)
        native_difference = max(abs(item["target_fraction_remaining"] - value)
                                for item, value in zip(solved[0]["target_profile"], cantera_profile))
        if native_difference > PARENT_AGREEMENT:
            raise RuntimeError("native_linear_cantera_and_linear_baseline_agreement_failed")
        _validate_sources(sources)
        if _exports(native) != snapshot["source_artifacts"]:
            raise ValueError("native_linear_sources_changed_during_diagnostics")
        solved[0]["numerical_checks"].update(
            native_baseline_agreement_checked=True,
            maximum_native_cantera_parent_difference=native_difference,
            native_cantera_checks=cantera_checks,
        )
        result.update(status="completed", scenarios=solved,
                      condition=snapshot["condition"], source_artifacts=sources,
                      numerical_checks={"status": "passed", "source_artifacts_unchanged": True,
                                        "native_baseline_checked": True,
                                        "maximum_native_cantera_parent_difference": native_difference})
    except ImportError as error:
        result.update(status="unavailable", failure_reason="native_linear_solver_dependency_unavailable: %s" % error)
    except (ValueError, RuntimeError, OSError, TypeError, KeyError, OverflowError) as error:
        result["failure_reason"] = "%s: %s" % (type(error).__name__, error)
    return result
