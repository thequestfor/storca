"""Append-only, five-point OH/methane rate audit; never changes production defaults.

Stages deliberately separate a zero-ORCA-job native import check from the five
new correlated single points. A failed stage is retained, not overwritten or
automatically retried. Existing frequency/IRC evidence remains immutable.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import sys
import threading
import time
import shutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.orca_runner import find_orca, run_orca, run_bounded_kinetics_command
from storca.arkane_runner import VerifiedArkaneRouteSpec, run_verified_arkane_route

SOURCE = ROOT / 'runs/oh-methane-geometry-frequency-20261010-lX9ABY'
POINTS = ('OH', 'CH4', 'H2O', 'CH3', 'TS')
STATES = {'OH': (0, 2), 'CH4': (0, 1), 'H2O': (0, 1), 'CH3': (0, 2), 'TS': (0, 2)}
GRID = (250.0, 273.15, 298.15, 323.15, 340.3, 400.0, 500.0, 750.0, 1000.0)
EH_KJ_MOL = 2625.4996394799


def artifact(path):
    path = Path(path).resolve(strict=True)
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(bound):
    if artifact(bound['path']) != bound:
        raise ValueError('Bound source changed: ' + bound['path'])


def write_new(path, document):
    with Path(path).open('x') as stream:
        json.dump(document, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def inventory(folder):
    return {str(p.resolve()): [artifact(p)['sha256'], p.stat().st_mtime_ns]
            for p in Path(folder).rglob('*') if p.is_file()}


def geometry(path):
    blocks = re.findall(r'(?ims)^\s*\*\s+xyz\s+([-+]?\d+)\s+(\d+)\s*\n(.*?)^\s*\*\s*$',
                        Path(path).read_text())
    if len(blocks) != 1:
        raise ValueError('One inline geometry required')
    charge, multiplicity, coordinates = blocks[0]
    rows = [line.split() for line in coordinates.splitlines() if line.strip()]
    if not rows or any(len(row) != 4 or not all(math.isfinite(float(v)) for v in row[1:]) for row in rows):
        raise ValueError('Invalid Cartesian geometry')
    return int(charge), int(multiplicity), coordinates.rstrip() + '\n'


def retained_points():
    report = json.loads((SOURCE / 'independent-result.json').read_text())
    route = json.loads((SOURCE / 'numerical-continuation-result.json').read_text())
    if (not route.get('route_verified') or not route.get('numerical_local_derivative_admission')
            or not route.get('prior_scientific_artifacts_unchanged')):
        raise ValueError('Numerical route admission absent')
    points = {}
    for name in POINTS:
        row = report['jobs'][name]['num001_frequency' if name == 'TS' else 'analytic_frequency']
        for field in ('input', 'output', 'hessian'):
            verify(row[field])
        charge, multiplicity, coordinates = geometry(row['input']['path'])
        if (charge, multiplicity) != STATES[name]:
            raise ValueError('Retained electronic state mismatch: ' + name)
        if len(row['negative_physical_mode_indices']) != int(name == 'TS'):
            raise ValueError('Retained physical mode inventory mismatch: ' + name)
        points[name] = {**row, 'charge': charge, 'multiplicity': multiplicity, 'coordinates': coordinates}
    return points


def inspect_statmech(folder, rmg_env):
    # The selected env is a native dependency, not an implicit global default.
    interpreter = Path('/opt/homebrew/anaconda3/envs') / rmg_env / 'bin/python'
    if not interpreter.is_file():
        raise ValueError('Specified RMG interpreter absent: ' + str(interpreter))
    script = ROOT / 'tests/oh_methane_statmech_inspect.py'
    completed = run_bounded_kinetics_command(
        [str(interpreter), '-B', str(script), str(folder / 'input.py'),
         str(folder / 'arkane-rotor-contract.json'), str(folder / 'native-statmech-audit.json')],
        cwd=folder, timeout_seconds=1800)
    write_new(folder / 'statmech-replay-execution.json', {
        'returncode': completed.returncode, 'stdout': completed.stdout, 'stderr': completed.stderr,
        'interpreter': artifact(interpreter), 'script': artifact(script)})
    if completed.returncode != 0:
        raise ValueError('Native statmech replay refused; retained its execution streams')
    result = json.loads((folder / 'native-statmech-audit.json').read_text())
    if result.get('status') != 'passed':
        raise ValueError('Native statmech audit incomplete')
    return result


def baseline(work, rmg_env):
    work.mkdir(parents=True, exist_ok=True)
    points = retained_points()
    before = inventory(SOURCE)
    write_new(work / 'contract.json', {
        'kind': 'oh_methane_consistent_rate_audit', 'schema_version': 1,
        'driver': artifact(__file__), 'route_evidence': artifact(SOURCE / 'numerical-continuation-result.json'),
        'points': points, 'rmg_env': rmg_env, 'baseline_frequency_scale_factor': 1.0,
        'prior_inventory': before,
        'planned_ORCA_jobs': ['ccsd-t-qzvpp/' + name for name in POINTS],
        'resources': {'cores': 2, 'maxcore_MB_per_worker': 2500, 'per_job_seconds': 21600,
                      'overall_seconds': 86400, 'starting_free_bytes': 40 * 1024**3,
                      'runtime_reserve_bytes': 4 * 1024**3},
        'claim_scope': 'isolated_single_surface_gas_RRHO_rate_diagnostic_not_calibration_or_lifetime',
        'old_geometry_correlated_energies_reused': False,
        'production_defaults_changed': False})
    output = work / 'baseline'
    output.mkdir()
    spec = VerifiedArkaneRouteSpec(
        label='oh-methane-refined-diagnostic', reactant_labels=('OH', 'CH4'),
        reactant_smiles=('[OH]', 'C'), reactant_orca_outputs=tuple(Path(points[n]['output']['path']) for n in ('OH', 'CH4')),
        reactant_multiplicities=(2, 1), product_labels=('H2O', 'CH3'), product_smiles=('O', '[CH3]'),
        product_orca_outputs=tuple(Path(points[n]['output']['path']) for n in ('H2O', 'CH3')),
        product_multiplicities=(1, 2), transition_state_orca_output=Path(points['TS']['output']['path']),
        transition_state_multiplicity=2, model_chemistry='WB97X-D3/def2-TZVPP',
        temperatures_K=GRID, pressures_bar=(1.0,), bath_gas=(), pressure_dependent=False,
        interpolation_model=('PDepArrhenius',),
        tunneling='Eckart', frequency_scale_factor=1.0, condition_temperature_K=340.3,
        collision_limit_m3_mol_s=190571421.80130282,
        collision_limit_source='Retained OH/methane RMG transport collision ceiling; diagnostic only')
    result = run_verified_arkane_route(spec, output, rmg_env=rmg_env)
    write_new(output / 'execution-result.json', result)
    if result['status'] != 'completed':
        raise ValueError('Baseline Arkane stage refused: ' + str(result.get('failure_reason') or
                         result.get('generated_kinetics_library', {}).get('failure_reason')))
    native = inspect_statmech(output, rmg_env)
    unchanged = inventory(SOURCE) == before
    write_new(work / 'baseline-result.json', {
        'status': 'completed' if unchanged else 'failed_source_mutation', 'source_unchanged': unchanged,
        'execution': artifact(output / 'execution-result.json'),
        'statmech_audit': artifact(output / 'native-statmech-audit.json'),
        'electronic_structure_jobs_launched': 0, 'production_defaults_changed': False,
        'physical_rate_accuracy_claim_allowed': False, 'lifetime_claim_allowed': False})
    if not unchanged:
        raise ValueError('Original scientific evidence changed')
    print('BASELINE NATIVE IMPORT/RATE AUDIT PASSED; no ORCA dispatch', flush=True)


def read_cc_result(folder, name):
    text = (folder / 'energy.out').read_text(errors='replace')
    if (text.count('ORCA TERMINATED NORMALLY') != 1 or 'SCF CONVERGED AFTER' not in text
            or '--- The Coupled-Cluster iterations have converged ---' not in text):
        raise ValueError('Missing native normal SCF/CC completion')
    energy = re.findall(r'^\s*FINAL SINGLE POINT ENERGY\s+([-+0-9.Ee]+)\s*$', text, re.M)
    totals = re.findall(r'E\(CCSD\(T\)\)\s+\.\.\.\s+([-+0-9.Ee]+)', text)
    if (len(energy) != 1 or len(totals) != 1 or not math.isfinite(float(energy[0]))
            or abs(float(energy[0]) - float(totals[0])) > 1e-8):
        raise ValueError('CCSD(T) final energy missing or inconsistent')
    actual = text.split('****END OF INPUT****', 1)[-1]
    ref = re.findall(r'Reference Wavefunction\s+\.\.\.\s+(\w+)', actual)
    expected = 'RHF' if STATES[name][1] == 1 else 'UHF'
    frozen = re.findall(r'Frozen core treatment\s+\.\.\.\s+chemical core \((\d+) el\)', actual)
    if ref != [expected] or frozen != [str(4 if name == 'TS' else 2)] or 'QR-MO GENERATION' in actual:
        raise ValueError('Reference or frozen-core convention mismatch')
    charge = re.findall(r'Total Charge\s+Charge\s+\.{2,}\s+([-+]?\d+)', actual)
    multiplicity = re.findall(r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)', actual)
    basis = re.findall(r'Number of basis functions\s+\.\.\.\s+(\d+)', actual)
    correlated = re.findall(r'Number of correlated electrons\s+\.\.\.\s+(\d+)', actual)
    expected_basis = {'OH': 87, 'CH4': 177, 'H2O': 117, 'CH3': 147, 'TS': 264}[name]
    expected_electrons = {'OH': 7, 'CH4': 8, 'H2O': 8, 'CH3': 7, 'TS': 15}[name]
    if (charge != ['0'] or multiplicity != [str(STATES[name][1])]
            or not basis or set(basis) != {str(expected_basis)}
            or correlated != [str(expected_electrons)]):
        raise ValueError('Native charge/spin/QZ basis/correlated-electron inventory mismatch')
    table = re.search(r'Iter\s+E\(tot\).*?\n(.*?)--- The Coupled-Cluster iterations have converged', text, re.S)
    rows = [] if table is None else re.findall(
        r'^\s*\d+\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([\d.]+)', table[1], re.M)
    if not rows or float(rows[-1][3]) > 1e-6:
        raise ValueError('CC residual missing or >1e-6')
    return {'energy_Eh': float(energy[0]), 'reference': expected, 'frozen_core_electrons': int(frozen[0]),
            'basis_function_count': expected_basis, 'correlated_electrons': expected_electrons,
            'CC_residual': float(rows[-1][3]), 'input': artifact(folder / 'energy.inp'),
            'output': artifact(folder / 'energy.out'),
            'T1': [float(v) for v in re.findall(r'T1 diagnostic\s+\.\.\.\s+([\d.]+)', actual)],
            'preceding_SCF_S_squared': [float(v) for v in re.findall(r'Expectation value of <S\*\*2>\s*:\s*([\d.]+)', actual)]}


def energies(work):
    contract = json.loads((work / 'contract.json').read_text())
    verify(contract['driver'])
    verify(contract['route_evidence'])
    baseline_result = json.loads((work / 'baseline-result.json').read_text())
    if baseline_result['status'] != 'completed' or not baseline_result['source_unchanged']:
        raise ValueError('Baseline native rate check must pass before correlated dispatch')
    points = retained_points()
    if points != contract['points'] or inventory(SOURCE) != contract['prior_inventory']:
        raise ValueError('Retained source identity changed')
    if (work / 'energy-result.json').exists() or (work / 'energy-contract.json').exists():
        raise ValueError('Fresh energy stage required; failed attempts are not overwritten')
    resources = contract['resources']
    if shutil.disk_usage(work).free < resources['starting_free_bytes']:
        raise ValueError('Five-point sequential QZ audit requires 40GiB starting free space')
    write_new(work / 'energy-contract.json', {
        'driver': artifact(__file__), 'source_contract': artifact(work / 'contract.json'),
        'baseline_admission': artifact(work / 'baseline-result.json'), 'orca': artifact(find_orca()),
        'recipe': 'canonical CCSD(T)/def2-QZVPP VeryTightSCF FrozenCore NoAutoStart; UseQROs false',
        'fixed_frequency_input_geometries': True, 'resources': resources,
        'initial_free_bytes': shutil.disk_usage(work).free,
        'memory_policy': 'MaxCore is per worker, not a hard total RSS limit',
        'limitations': ['Finite basis; no CBS or post-CCSD(T) correction.',
                        'Single points on validated DFT geometries; no high-level saddle claim.',
                        'No computed spin-orbit correction or reference stability proof.']})
    started = time.monotonic()
    results, failure = {}, None
    try:
        for name in POINTS:
            remaining = resources['overall_seconds'] - (time.monotonic() - started)
            if remaining <= 0 or shutil.disk_usage(work).free < resources['starting_free_bytes']:
                raise ValueError('Overall time or per-job 40GiB starting-space requirement reached')
            point = points[name]
            folder = work / 'ccsd-t-qzvpp' / name
            folder.mkdir(parents=True, exist_ok=False)
            inp = folder / 'energy.inp'
            text = ('! CCSD(T) def2-QZVPP NoAutoStart VeryTightSCF FrozenCore\n'
                    '%pal\n nprocs 2\nend\n%maxcore 2500\n'
                    '%mdci\n DoSDiagnostic true\n UseQROs false\nend\n'
                    f"* xyz {point['charge']} {point['multiplicity']}\n" + point['coordinates'] + '*\n')
            with inp.open('x') as stream:
                stream.write(text)
            if geometry(inp) != geometry(point['input']['path']):
                raise ValueError('New input differs from exact frequency-input coordinates/state')
            done, stopped = threading.Event(), []
            input_hash = artifact(inp)['sha256']
            def watch():
                while not done.wait(2):
                    if shutil.disk_usage(work).free >= resources['runtime_reserve_bytes']:
                        continue
                    stopped.append('Disk reserve reached during owned native job')
                    try:
                        lease = json.loads(inp.with_suffix('.process.json').read_text())
                        if (lease['status'] == 'running' and lease['input_sha256'] == input_hash
                                and Path(lease['input']) == inp.resolve()
                                and os.getpgid(lease['pid']) == lease['process_group_id']):
                            os.killpg(lease['process_group_id'], signal.SIGINT)
                    except (OSError, ValueError, KeyError):
                        pass
                    return
            def progress(line):
                if re.search(r'SCF CONVERGED|Coupled.Cluster.*converged|T1 diagnostic|FINAL SINGLE POINT|'
                             r'ORCA TERMINATED|Memory needed for triples|Final InCoreLevel|error termination', line):
                    print(line.rstrip(), flush=True)
                return bool(stopped)
            job_start = time.monotonic()
            watcher = threading.Thread(target=watch, daemon=True)
            watcher.start()
            print(f'START {name}/QZ: 2 MPI workers,2500MB/worker, up to6h; exact refined coordinates', flush=True)
            try:
                run_orca(inp, timeout_seconds=min(remaining, resources['per_job_seconds']), stream_output=progress,
                         extra_env={key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                                       'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')})
                if stopped:
                    raise ValueError(stopped[0])
                lease = json.loads(inp.with_suffix('.process.json').read_text())
                if lease['returncode'] != 0 or lease['input_sha256'] != input_hash:
                    raise ValueError('Owned native execution failed or source changed')
                result = read_cc_result(folder, name)
                result.update(elapsed_seconds=time.monotonic() - job_start,
                              lease=artifact(inp.with_suffix('.process.json')))
                write_new(folder / 'result.json', result)
                results[name] = result
                print(f"DONE {name}/QZ: {result['energy_Eh']:.12f}Eh", flush=True)
            finally:
                done.set()
                watcher.join(timeout=3)
    except Exception as error:
        failure = type(error).__name__ + ': ' + str(error)
    unchanged = inventory(SOURCE) == contract['prior_inventory']
    result = {'status': 'completed' if not failure and unchanged else 'failed', 'failure_reason': failure,
              'points': results, 'source_unchanged': unchanged, 'elapsed_seconds': time.monotonic() - started,
              'physical_rate_accuracy_claim_allowed': False, 'lifetime_claim_allowed': False}
    if len(results) == 5:
        result['forward_electronic_barrier_kJ_mol'] = (results['TS']['energy_Eh'] - results['OH']['energy_Eh'] - results['CH4']['energy_Eh']) * EH_KJ_MOL
        result['reverse_electronic_barrier_kJ_mol'] = (results['TS']['energy_Eh'] - results['H2O']['energy_Eh'] - results['CH3']['energy_Eh']) * EH_KJ_MOL
    write_new(work / 'energy-result.json', result)
    if result['status'] != 'completed':
        raise ValueError(failure or 'Prior scientific evidence changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('baseline', 'energies'))
    parser.add_argument('workdir', type=Path)
    parser.add_argument('--rmg-env', default='rmg4_env')
    args = parser.parse_args()
    work = args.workdir.resolve()
    if ROOT / 'runs' not in work.parents:
        raise ValueError('Audit workdir must be an owned child of repository runs/')
    # Local writable cache avoids modifying global matplotlib state.
    cache = work / 'matplotlib-cache'
    cache.mkdir(parents=True, exist_ok=True)
    os.environ['MPLCONFIGDIR'] = str(cache)
    if args.stage == 'baseline':
        baseline(work, args.rmg_env)
    else:
        energies(work)


if __name__ == '__main__':
    main()
