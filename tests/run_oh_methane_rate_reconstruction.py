"""Explicit-method revision of the append-only OH/methane rate experiment.

The first failed string-only method declaration remains unchanged. This driver
uses the source-separated recipe API and never rewrites ORCA native evidence.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location('bounded_oh_composite_driver', Path(__file__).with_name('run_oh_methane_composite_audit.py'))
bounded = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bounded)
from src.orca_runner import run_arkane
from storca.arkane_rotor import validate_rotor_receipt
from storca.composite_rate_recipe import CompositeRateRecipe, bind_artifact, bind_stationary_point, write_isolated_arkane_recipe

FIRST_ATTEMPT = ROOT / 'runs/oh-methane-composite-rate-20261011-5ivCKk'
REFERENCE = ROOT / 'benchmarks/stability/reference_sources/2026-10-07-berg2024-oh-methane/source-audit.json'
SYMMETRIES = {'OH': 1, 'CH4': 12, 'H2O': 2, 'CH3': 6, 'TS': 1}
SMILES = {'OH': '[OH]', 'CH4': 'C', 'H2O': 'O', 'CH3': '[CH3]', 'TS': None}


def recipe(points, energies, method):
    bindings = tuple(bind_stationary_point(
        name, 'transition_state' if name == 'TS' else 'species', Path(points[name]['output']['path']),
        Path(energies[name]['output']['path']), energy_hartree=energies[name]['energy_Eh'],
        charge=points[name]['charge'], multiplicity=points[name]['multiplicity'],
        external_symmetry=SYMMETRIES[name], optical_isomers=1, smiles=SMILES[name])
        for name in bounded.POINTS)
    return CompositeRateRecipe(
        label='oh-methane-refined-diagnostic', points=bindings, reactant_labels=('OH', 'CH4'),
        product_labels=('H2O', 'CH3'), transition_state_label='TS', electronic_method=method,
        vibrational_method='WB97X-D3/def2-TZVPP', frequency_scale_factor=1.0,
        temperatures_K=bounded.GRID, tunneling='Eckart', environment='gas',
        surface_model='single_adiabatic_surface',
        electronic_state_assumption='Neutral doublet reaction; RHF singlets/UHF doublets. Spin-only electronic partition baseline, no computed SOC.',
        channel_counting_assumption='Declared rotational symmetries OH1/CH4_12/H2O2/CH3_6/TS1, optical1; no extra methane times-four multiplier. TS chirality/counting requires native diagnostic review.',
        soft_mode_assumption='All positive physical modes retained as harmonic; weak relative OH angular mode is not deleted, floored or inferred to be a covalent rotor.',
        route_evidence=(bind_artifact(bounded.SOURCE / 'numerical-continuation-result.json'),))


def run_rate(work, stage, points, energy_points, method, rmg_env):
    folder = work / stage
    staged = write_isolated_arkane_recipe(recipe(points, energy_points, method), folder)
    try:
        result = run_arkane(folder / 'input.py', rmg_env=rmg_env, timeout_seconds=1800,
                            statmech_contract=folder / 'arkane-rotor-contract.json')
        bounded.write_new(folder / 'execution-result.json', {'status': 'completed', **result, 'staging': staged})
        validate_rotor_receipt(folder)
        audit = bounded.inspect_statmech(folder, rmg_env)
        return audit
    except Exception as error:
        if not (folder / 'execution-result.json').exists():
            bounded.write_new(folder / 'execution-result.json', {
                'status': 'failed', 'failure_reason': type(error).__name__ + ': ' + str(error),
                'stdout': getattr(error, 'stdout', None), 'stderr': getattr(error, 'stderr', None),
                'staging': staged})
        raise


def baseline(work, rmg_env):
    contract = json.loads((FIRST_ATTEMPT / 'contract.json').read_text())
    bounded.verify(contract['driver'])
    points = bounded.retained_points()
    if points != contract['points'] or bounded.inventory(bounded.SOURCE) != contract['prior_inventory']:
        raise ValueError('Preserved first-stage source bindings changed')
    contract.update(revision_driver=bounded.artifact(__file__),
                    composite_builder=bounded.artifact(ROOT / 'storca/composite_rate_recipe.py'),
                    failed_method_declaration_attempt={
                        'contract': bounded.artifact(FIRST_ATTEMPT / 'contract.json'),
                        'result': bounded.artifact(FIRST_ATTEMPT / 'baseline/execution-result.json')},
                    rmg_env=rmg_env, method_declaration_revision='Explicit native LevelOfTheory with software=orca; physical method unchanged')
    bounded.write_new(work / 'contract.json', contract)
    run_rate(work, 'baseline', points, points, 'WB97X-D3/def2-TZVPP', rmg_env)
    unchanged = bounded.inventory(bounded.SOURCE) == contract['prior_inventory']
    bounded.write_new(work / 'baseline-result.json', {
        'status': 'completed' if unchanged else 'failed_source_mutation', 'source_unchanged': unchanged,
        'execution': bounded.artifact(work / 'baseline/execution-result.json'),
        'statmech_audit': bounded.artifact(work / 'baseline/native-statmech-audit.json'),
        'electronic_structure_jobs_launched': 0, 'production_defaults_changed': False,
        'physical_rate_accuracy_claim_allowed': False, 'lifetime_claim_allowed': False})
    if not unchanged:
        raise ValueError('Original scientific evidence changed')
    print('BASELINE: explicit-method native import and rate reconstruction complete; zero ORCA jobs', flush=True)


def composite(work, rmg_env):
    contract = json.loads((work / 'contract.json').read_text())
    for key in ('driver', 'revision_driver', 'composite_builder', 'route_evidence'):
        bounded.verify(contract[key])
    energy = json.loads((work / 'energy-result.json').read_text())
    if energy['status'] != 'completed' or not energy['source_unchanged']:
        raise ValueError('Matched five-point energy stage incomplete')
    for name, result in energy['points'].items():
        bounded.verify(result['input'])
        bounded.verify(result['output'])
        parsed = bounded.read_cc_result(Path(result['output']['path']).parent, name)
        if parsed['energy_Eh'] != result['energy_Eh']:
            raise ValueError('Bound energy differs from native result')
    points = bounded.retained_points()
    if points != contract['points']:
        raise ValueError('Frequency geometry/state basis changed')
    run_rate(work, 'composite', points, energy['points'], 'CCSD(T)/def2-QZVPP', rmg_env)
    unchanged = bounded.inventory(bounded.SOURCE) == contract['prior_inventory']
    bounded.write_new(work / 'composite-result.json', {
        'status': 'completed' if unchanged else 'failed_source_mutation', 'source_unchanged': unchanged,
        'statmech_audit': bounded.artifact(work / 'composite/native-statmech-audit.json'),
        'energy_stage': bounded.artifact(work / 'energy-result.json'),
        'recipe': bounded.artifact(work / 'composite/composite-rate-recipe.json'),
        'production_defaults_changed': False, 'network_reinjection_allowed': False,
        'physical_rate_accuracy_claim_allowed': False, 'lifetime_claim_allowed': False})
    if not unchanged:
        raise ValueError('Original scientific evidence changed')
    print('COMPOSITE: matched five-species native rate reconstruction complete', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('baseline', 'energies', 'composite'))
    parser.add_argument('workdir', type=Path)
    parser.add_argument('--rmg-env', default='rmg4_env')
    args = parser.parse_args()
    work = args.workdir.resolve()
    if ROOT / 'runs' not in work.parents:
        raise ValueError('Owned repository runs child required')
    cache = work / 'matplotlib-cache'
    cache.mkdir(parents=True, exist_ok=True)
    os.environ['MPLCONFIGDIR'] = str(cache)
    if args.stage == 'baseline':
        baseline(work, args.rmg_env)
    elif args.stage == 'energies':
        bounded.energies(work)
    else:
        composite(work, args.rmg_env)


if __name__ == '__main__':
    main()
