"""Score exact native rates against the frozen OH/methane measurement slice.

Descriptive development comparison only: no fitting to measurements, promotion,
calibrated confidence interval, pressure applicability or lifetime claim.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
REFERENCE = ROOT / 'benchmarks/stability/reference_sources/2026-10-07-berg2024-oh-methane/source-audit.json'
AVOGADRO_EXACT = 6.02214076e23


def artifact(path):
    path = Path(path).resolve(strict=True)
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def reopen(record):
    observed = artifact(record['path'])
    if observed != {key: record[key] for key in ('path', 'sha256')}:
        raise ValueError('Closeout source hash changed: ' + record['path'])
    if 'mtime_ns' in record and Path(record['path']).stat().st_mtime_ns != record['mtime_ns']:
        raise ValueError('Closeout source modification time changed: ' + record['path'])


def closeout(work):
    """Reopen the evidence graph; stored success booleans alone are insufficient."""
    from storca.arkane_rotor import validate_rotor_receipt
    contract = json.loads((work / 'contract.json').read_text())
    for key in ('driver', 'revision_driver', 'route_evidence'):
        reopen(contract[key])
    spec = importlib.util.spec_from_file_location('oh_rate_completion', Path(__file__).with_name('complete_oh_methane_rate_audit.py'))
    completion = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(completion)
    revision = completion.validate_revision_receipt(work)
    for path, values in contract['prior_inventory'].items():
        reopen({'path': path, 'sha256': values[0], 'mtime_ns': values[1]})
    energy_contract = json.loads((work / 'energy-contract.json').read_text())
    for key in ('driver', 'source_contract', 'baseline_admission', 'orca'):
        reopen(energy_contract[key])
    energy_result = json.loads((work / 'energy-result.json').read_text())
    if energy_result['status'] != 'completed' or not energy_result['source_unchanged']:
        raise ValueError('Correlated-energy stage not closed')
    for point in energy_result['points'].values():
        for key in ('input', 'output', 'lease'):
            reopen(point[key])
    stage_rows = {}
    for stage in ('baseline', 'composite'):
        result_path = work / (stage + '-result.json')
        stage_result = json.loads(result_path.read_text())
        if stage_result['status'] != 'completed' or not stage_result['source_unchanged']:
            raise ValueError('Rate stage not closed: ' + stage)
        for key in ('statmech_audit', 'execution', 'recipe', 'energy_stage', 'builder_revision_admission', 'completion_driver'):
            if key in stage_result:
                reopen(stage_result[key])
        if stage == 'composite' and stage_result['builder_revision_admission'] != revision:
            raise ValueError('Composite stage used a different source revision')
        folder = work / stage
        validate_rotor_receipt(folder)
        provenance = json.loads((folder / 'provenance.json').read_text())
        reopen(provenance['composite_recipe'])
        recipe = json.loads((folder / 'composite-rate-recipe.json').read_text())
        for point in recipe['points']:
            for key in ('frequency_input', 'frequency_output', 'hessian', 'energy_input', 'energy_output'):
                reopen(point[key])
        for bound in recipe['recipe']['route_evidence']:
            reopen(bound)
        audit = json.loads((folder / 'native-statmech-audit.json').read_text())
        reopen(audit['input'])
        reopen(audit['contract'])
        reopen(audit['helper'])
        for bound in audit['prior_source_inventory'].values():
            reopen(bound)
        for key in ('interpreter', 'symmetry_executable'):
            reopen(audit['native'][key])
        for bound in audit['native']['modules'].values():
            reopen(bound)
        stage_rows[stage] = artifact(result_path)
    return {'status': 'passed', 'stage_results': stage_rows,
            'builder_revision_admission': revision,
            'energy_contract': artifact(work / 'energy-contract.json'),
            'energy_result': artifact(work / 'energy-result.json'),
            'external_energy_sources_reopened': True, 'stored_success_flags_only': False}


def score(audit, observations):
    if audit.get('status') != 'passed' or not audit.get('all_admitted_and_previous_artifacts_unchanged'):
        raise ValueError('A passed, preserved native rate diagnostic is required')
    if len(observations) != 14 or len({point['point_id'] for point in observations}) != 14:
        raise ValueError('All fourteen frozen measurement rows are required')
    result = {}
    for name in ('classical', 'Wigner', 'Eckart'):
        model = audit['rates'][name]
        by_T = {row['temperature_K']: row for row in model['grid']}
        if len(by_T) != len(model['grid']):
            raise ValueError('Duplicate native evaluation temperature')
        rows = []
        for point in observations:
            if point['rate_units'] != 'cm^3 molecule^-1 s^-1':
                raise ValueError('Unexpected frozen experimental units')
            T = point['temperature_K']
            if T not in by_T:
                raise ValueError('Direct rate absent at exact measured temperature')
            predicted = by_T[T]['forward_direct_m3_mol_s']
            measured = point['rate'] * AVOGADRO_EXACT * 1e-6
            if any(not math.isfinite(v) or v <= 0 for v in (predicted, measured)):
                raise ValueError('Invalid finite positive prediction/reference')
            ratio = predicted / measured
            rows.append({'point_id': point['point_id'], 'temperature_K': T,
                         'direct_predicted_m3_mol_s': predicted, 'measured_m3_mol_s': measured,
                         'prediction_over_measurement': ratio,
                         'absolute_log_error': abs(math.log(ratio)),
                         'fit_over_direct': by_T[T]['fit_over_direct'],
                         'tunneling_factor': by_T[T]['tunneling_factor']})
        result[name] = {'rows': rows, 'mean_absolute_log_error': sum(r['absolute_log_error'] for r in rows) / 14,
                        'ratio_min': min(r['prediction_over_measurement'] for r in rows),
                        'ratio_max': max(r['prediction_over_measurement'] for r in rows),
                        'frozen_primary_point': next(r for r in rows if r['point_id'] == 'S1-A-01'),
                        'maximum_absolute_log_fit_error_original_grid': model['maximum_absolute_log_fit_error_original_grid']}
    return result


def oh_soc_classical_control(audit):
    """Declared asymptotic sensitivity, not a computed spin-orbit correction.

    Spin-free OH energy is treated as the centroid of equal-degeneracy levels
    split by140cm-1, with quenched TS SOC. Shift the reactant energy AND use
    both populations. This formula applies only to bare TST here: changing
    barriers also changes Eckart transmission, which is not inferred.
    """
    split_J_molecule = 6.62607015e-34 * 299792458.0 * 100.0 * 140.0
    rows = []
    for row in audit['rates']['classical']['grid']:
        T = row['temperature_K']
        exponent = split_J_molecule / (1.380649e-23 * T)
        population_factor = 1.0 / (1.0 + math.exp(-exponent))
        combined = math.exp(-exponent / 2.0) * population_factor
        rows.append({'temperature_K': T, 'population_only_factor': population_factor,
                     'energy_and_population_factor': combined,
                     'classical_control_m3_mol_s': row['forward_direct_m3_mol_s'] * combined})
    return {'kind': 'asymptotic_OH_SOC_bare_TST_sensitivity_not_computed_SOC',
            'splitting_cm_inverse': 140.0, 'OH_ground_lowering_J_mol': split_J_molecule * AVOGADRO_EXACT / 2.0,
            'TS_SOC_quenched_assumed': True, 'applies_to_Eckart': False,
            'source': 'https://www2.sci.u-szeged.hu/czako/papers/PCCP_22_14560_2020.pdf',
            'rows': rows, 'accuracy_claim_allowed': False}


def report(work):
    closure = closeout(work)
    frozen = json.loads(REFERENCE.read_text())['audited_observation_sets'][0]
    if not frozen['transcription_frozen']:
        raise ValueError('Measured source is not frozen')
    stages = {}
    for stage in ('baseline', 'composite'):
        path = work / stage / 'native-statmech-audit.json'
        audit = json.loads(path.read_text())
        if audit['extra_evaluation_temperature_source'] != artifact(REFERENCE):
            raise ValueError('Frozen reference source changed between evaluation and comparison')
        recipe_path = work / stage / 'composite-rate-recipe.json'
        recipe = json.loads(recipe_path.read_text())
        for direction in ('forward', 'reverse'):
            native = audit['barriers_J_mol'][direction + '_E0']
            expected = recipe['diagnostic_barriers'][direction]['native_arkane_expected_E0_barrier_J_mol']
            if abs(native - expected) > 0.1:
                raise ValueError('Native E0 does not reproduce bound recipe within printed-frequency rounding allowance')
        stages[stage] = {'native_audit': artifact(path), 'comparison': score(audit, frozen['observations']),
                         'composite_recipe': artifact(recipe_path),
                         'native_E0_recipe_rounding_allowance_J_mol': 0.1,
                         'native_E0_barriers_J_mol': audit['barriers_J_mol'],
                         'counting_warnings': audit['warnings'],
                         'OH_SOC_classical_sensitivity': oh_soc_classical_control(audit),
                         'optical_count_sensitivity': {'declared_TS_optical_count': audit['stationary_points'][audit['transition_state']]['imposed_optical_isomers'],
                                                       'alternative_TS_count': 2,
                                                       'rate_multiplier_if_independent_mirror_paths_required': 2,
                                                       'alternative_not_applied': True,
                                                       'meaning': 'A counting sensitivity, not evidence that double counting is appropriate.'}}
    document = {'schema_version': 1, 'kind': 'oh_methane_native_rate_vs_frozen_reference',
                'status': 'completed_conditional_development_comparison',
                'evidence_closeout': closure,
                'reference': artifact(REFERENCE), 'report_driver': artifact(__file__), 'stages': stages,
                'Avogadro_constant_conversion': AVOGADRO_EXACT,
                'independent_experimental_studies': 1, 'measurement_rows': 14,
                'experimental_condition': 'Nominal1000hPa humid synthetic air; total methane-dependent OH loss',
                'calculated_condition': 'Single adiabatic surface, isolated ideal-gas RRHO/TST direct abstraction',
                'environment_applicability': 'conditional_not_certified',
                'reference_used_for_fitting_or_method_selection': False,
                'physical_rate_accuracy_claim_allowed': False, 'calibrated_uncertainty_claim_allowed': False,
                'network_reinjection_allowed': False, 'lifetime_claim_allowed': False,
                'limitations': ['This is a development study, not a prospective holdout or fourteen independent validations.',
                                'Electronic-state, soft-mode/RRHO, recrossing and tunneling approximations remain unresolved.',
                                'OH SOC and optical-count controls are assumptions, not confidence bounds or silently corrected rates.',
                                'Finite QZ basis and correlated single points at DFT geometries are not ground truth.',
                                'One connected reaction does not establish total pathway coverage or a molecular lifetime.']}
    with (work / 'rate-comparison.json').open('x') as stream:
        json.dump(document, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    for stage, value in stages.items():
        print(stage, {name: row['frozen_primary_point']['prediction_over_measurement']
                      for name, row in value['comparison'].items()})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdir', type=Path)
    report(parser.parse_args().workdir.resolve(strict=True))
