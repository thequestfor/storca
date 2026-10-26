"""Append-only completion after the precisely bounded CC energy-parser repair.

Never rewrite the original contract or repeat electronic-structure jobs. This
experiment-specific receipt admits exactly one f-string quantifier correction,
not arbitrary changed source code or a changed scientific model.
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
spec = importlib.util.spec_from_file_location('oh_rate_reconstruction', Path(__file__).with_name('run_oh_methane_rate_reconstruction.py'))
reconstruction = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reconstruction)
bounded = reconstruction.bounded
from storca.composite_rate_recipe import validate_composite_rate_recipe

NEW_LINE = '        totals = re.findall(rf"E\\(CCSD\\(T\\)\\)\\s+\\.{{2,}}\\s+({NUMBER})", native_energy)\n'
OLD_LINE = NEW_LINE.replace(r'\.{{2,}}', r'\.{2,}')


def prove_quantifier_revision(current_source, expected_old_hash):
    if current_source.count(NEW_LINE.encode()) != 1 or OLD_LINE.encode() in current_source:
        raise ValueError('Expected exactly one corrected CC parser line')
    reconstructed = current_source.replace(NEW_LINE.encode(), OLD_LINE.encode(), 1)
    reverse_hash = hashlib.sha256(reconstructed).hexdigest()
    if reverse_hash != expected_old_hash:
        raise ValueError('Source revision is not exclusively the declared quantifier repair')
    return {'old_line': OLD_LINE, 'new_line': NEW_LINE,
            'reverse_reconstructed_original_sha256': reverse_hash,
            'current_sha256': hashlib.sha256(current_source).hexdigest()}


def revision_evidence(work):
    contract = json.loads((work / 'contract.json').read_text())
    for key in ('driver', 'revision_driver', 'route_evidence'):
        bounded.verify(contract[key])
    if bounded.inventory(bounded.SOURCE) != contract['prior_inventory']:
        raise ValueError('Original scientific evidence changed')
    source = Path(contract['composite_builder']['path'])
    delta = prove_quantifier_revision(source.read_bytes(), contract['composite_builder']['sha256'])
    points = bounded.retained_points()
    if points != contract['points']:
        raise ValueError('Frequency geometry/state basis changed')
    baseline = json.loads((work / 'baseline-result.json').read_text())
    if baseline['status'] != 'completed' or not baseline['source_unchanged']:
        raise ValueError('Original baseline not closed')
    for key in ('execution', 'statmech_audit'):
        bounded.verify(baseline[key])
    assessment = validate_composite_rate_recipe(reconstruction.recipe(points, points, 'WB97X-D3/def2-TZVPP'))
    normalized = json.loads(json.dumps(assessment, allow_nan=False))
    manifest = work / 'baseline/composite-rate-recipe.json'
    if normalized != json.loads(manifest.read_text()):
        raise ValueError('Repaired builder does not reproduce the retained DFT recipe exactly')
    return {'schema_version': 1, 'kind': 'bounded_source_revision_admission',
            'original_contract': bounded.artifact(work / 'contract.json'),
            'original_builder': contract['composite_builder'],
            'admitted_builder': bounded.artifact(source), 'exact_delta': delta,
            'admission_driver': bounded.artifact(__file__),
            'baseline_manifest': bounded.artifact(manifest),
            'baseline_result': bounded.artifact(work / 'baseline-result.json'),
            'baseline_assessment_exactly_reproduced': True,
            'DFT_branch_unaffected': True, 'physical_model_changed': False,
            'reason': 'Escape the regex quantifier in the CCSD(T)-only f-string; no other source delta admitted.',
            'electronic_structure_jobs_launched': 0, 'original_contract_rewritten': False}


def validate_revision_receipt(work):
    expected = revision_evidence(work)
    path = work / 'builder-revision-admission.json'
    if json.loads(path.read_text()) != expected:
        raise ValueError('Source revision receipt does not reproduce its bounded evidence')
    return bounded.artifact(path)


def match_native_energy(parsed, point):
    if (set(point) != set(parsed) | {'elapsed_seconds', 'lease'}
            or not all(point[key] == value for key, value in parsed.items())
            or not isinstance(point['elapsed_seconds'], (int, float))
            or not math.isfinite(point['elapsed_seconds']) or point['elapsed_seconds'] < 0):
        raise ValueError('Bound correlated result differs from native evidence')


def complete(work, rmg_env):
    receipt = work / 'builder-revision-admission.json'
    if not receipt.exists():
        bounded.write_new(receipt, revision_evidence(work))
    revision = validate_revision_receipt(work)
    energy = json.loads((work / 'energy-result.json').read_text())
    if energy['status'] != 'completed' or not energy['source_unchanged']:
        raise ValueError('Matched five-point energy stage incomplete')
    if set(energy['points']) != set(bounded.POINTS):
        raise ValueError('All five matched energies required')
    for name, point in energy['points'].items():
        for key in ('input', 'output', 'lease'):
            bounded.verify(point[key])
        parsed = bounded.read_cc_result(Path(point['output']['path']).parent, name)
        match_native_energy(parsed, point)
    reconstruction.run_rate(work, 'composite', bounded.retained_points(), energy['points'],
                            'CCSD(T)/def2-QZVPP', rmg_env)
    validate_revision_receipt(work)
    bounded.write_new(work / 'composite-result.json', {
        'status': 'completed', 'source_unchanged': True,
        'statmech_audit': bounded.artifact(work / 'composite/native-statmech-audit.json'),
        'execution': bounded.artifact(work / 'composite/execution-result.json'),
        'energy_stage': bounded.artifact(work / 'energy-result.json'),
        'recipe': bounded.artifact(work / 'composite/composite-rate-recipe.json'),
        'builder_revision_admission': revision, 'completion_driver': bounded.artifact(__file__),
        'production_defaults_changed': False, 'network_reinjection_allowed': False,
        'physical_rate_accuracy_claim_allowed': False, 'lifetime_claim_allowed': False})
    print('COMPOSITE: matched five-point rate reconstructed; source revision explicitly admitted', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdir', type=Path)
    parser.add_argument('--rmg-env', default='rmg4_env')
    args = parser.parse_args()
    work = args.workdir.resolve(strict=True)
    if ROOT / 'runs' not in work.parents:
        raise ValueError('Owned repository runs child required')
    complete(work, args.rmg_env)
