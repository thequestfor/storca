"""Strict replacement correspondence, executed in the owning RMG environment.

No application imports or dependency on the application's Python version.  This
proves endpoint identity and source attribution, never rate accuracy, atom
mapping along a path, stereochemistry, or kinetic-model applicability.
"""

import hashlib
import json
from collections import Counter
from pathlib import Path


POLICY_ID = "rmg-resonance-strict-v1"
SOURCE_FILES = ("reactions.py", "dictionary.txt", "generated-kinetics.json")


class IdentityUnresolved(ValueError):
    """A specific, auditable reason not to assert correspondence."""


def _digest(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def mechanism_artifacts(payload):
    """Snapshot exact inputs before the mechanism is interpreted by RMG."""
    result = {}
    for key in ("chemkin", "dictionary", "transport"):
        value = payload.get(key)
        if value:
            path = Path(value).resolve(strict=True)
            result[key] = {"path": str(path), "sha256": _digest(path)}
    if not all(key in result for key in ("chemkin", "dictionary")):
        raise IdentityUnresolved("mechanism_structure_artifacts_missing")
    return result


def _check_sources(request):
    library = Path(request.get("library") or "")
    if not library.is_absolute() or not library.is_dir():
        raise IdentityUnresolved("replacement_library_path_invalid")
    library = library.resolve(strict=True)
    artifacts = request.get("source_artifacts")
    if not isinstance(artifacts, dict) or set(artifacts) != set(SOURCE_FILES):
        raise IdentityUnresolved("replacement_source_artifacts_incomplete")
    for name in SOURCE_FILES:
        record = artifacts[name]
        if not isinstance(record, dict):
            raise IdentityUnresolved("replacement_source_artifact_invalid")
        path = Path(record.get("path") or "")
        if (not path.is_absolute() or path.resolve(strict=True) != library / name
                or not path.is_file()):
            raise IdentityUnresolved("replacement_source_artifact_path_mismatch")
        if _digest(path) != record.get("sha256"):
            raise IdentityUnresolved("replacement_source_artifact_hash_mismatch")
    manifest = json.loads((library / "generated-kinetics.json").read_text())
    if (manifest.get("kind") != "storca_generated_kinetics"
            or not str(manifest.get("verification_status") or "").startswith("verified")
            or not (manifest.get("purpose") == "kinetic_repair" or manifest.get("replaces_reaction"))):
        raise IdentityUnresolved("replacement_manifest_not_verified_repair")
    if not request.get("route_id") or manifest.get("route_id") != request["route_id"]:
        raise IdentityUnresolved("replacement_manifest_route_id_mismatch")
    if manifest.get("library") and Path(manifest["library"]).resolve() != library:
        raise IdentityUnresolved("replacement_manifest_library_path_mismatch")
    checksums = manifest.get("checksums") or {}
    if not all(name in checksums for name in ("reactions.py", "dictionary.txt")):
        raise IdentityUnresolved("replacement_manifest_checksums_missing")
    for name, digest in checksums.items():
        if Path(name).name != name or name in (".", ".."):
            raise IdentityUnresolved("replacement_manifest_checksum_path_invalid")
        if _digest(library / name) != digest:
            raise IdentityUnresolved("replacement_manifest_checksum_mismatch")
    return library, manifest


def _check_declared_stereo(value, key=""):
    from rdkit import Chem
    if isinstance(value, dict):
        for child_key, child in value.items():
            _check_declared_stereo(child, str(child_key).lower())
    elif isinstance(value, list):
        for child in value:
            _check_declared_stereo(child, key)
    elif isinstance(value, str) and "smiles" in key and value:
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            raise IdentityUnresolved("declared_structure_stereo_diagnosis_failed")
        if Chem.FindPotentialStereo(molecule):
            raise IdentityUnresolved("stereochemistry_not_preserved_by_rmg_artifacts")
    elif "stereo" in key and value not in (None, False, "", "none", "unspecified"):
        raise IdentityUnresolved("declared_stereochemistry_not_preserved_by_rmg_artifacts")


def _state(molecule):
    multiplicity = molecule.multiplicity
    if not isinstance(multiplicity, int) or multiplicity < 1 or not molecule.atoms:
        raise IdentityUnresolved("species_electronic_state_unavailable")
    composition = Counter()
    charge = 0
    for atom in molecule.atoms:
        symbol, isotope = atom.element.symbol, atom.element.isotope
        if symbol in ("X", "R", "L", "e") or not isinstance(isotope, int):
            raise IdentityUnresolved("species_graph_state_outside_identity_policy")
        composition[(symbol, isotope)] += 1
        atom_charge = atom.charge
        if int(atom_charge) != atom_charge:
            raise IdentityUnresolved("species_charge_unavailable")
        charge += int(atom_charge)
    return (multiplicity, charge, tuple(sorted(composition.items())))


def _prepared_species(species, cache):
    key = id(species)
    if key in cache:
        return cache[key][1:]
    if not getattr(species, "molecule", None):
        raise IdentityUnresolved("species_structure_unavailable")
    # Do not mutate the already-inspected mechanism or its original graph.
    prepared = species.copy(deep=True)
    initial_state = _state(prepared.molecule[0])
    if species.multiplicity != initial_state[0]:
        raise IdentityUnresolved("species_multiplicity_inconsistent")
    prepared.generate_resonance_structures(keep_isomorphic=True)
    if not prepared.molecule:
        raise IdentityUnresolved("species_resonance_generation_failed")
    for molecule in prepared.molecule:
        if _state(molecule) != initial_state:
            raise IdentityUnresolved("species_resonance_state_changed")
    # Keep the original object alive as well: source libraries are loaded one
    # at a time and a recycled object id must never reuse an earlier identity.
    cache[key] = (species, prepared, initial_state)
    return cache[key][1:]


def _check_species_stereo(prepared):
    from rdkit import Chem
    # Any supported resonance form revealing stereo sensitivity is enough to
    # withhold a claim; picking a stereo-blind resonance form would lose proof.
    for molecule in prepared.molecule:
        try:
            stereo = Chem.FindPotentialStereo(molecule.to_rdkit_mol())
        except Exception as error:
            raise IdentityUnresolved("species_stereo_diagnosis_failed") from error
        if stereo:
            raise IdentityUnresolved("stereochemistry_not_preserved_by_rmg_artifacts")


def _equivalent(source, candidate, cache):
    left, left_state = _prepared_species(source, cache)
    if not getattr(candidate, "molecule", None):
        raise IdentityUnresolved("species_structure_unavailable")
    if left_state != _state(candidate.molecule[0]):
        return False
    right, right_state = _prepared_species(candidate, cache)
    return left_state == right_state and left.is_isomorphic(right, strict=True)


def _match_side(source, model, cache):
    """Occurrence-level bijection, including repeated species; no degeneracy claim."""
    if len(source) != len(model):
        return None
    possibilities = [[index for index, candidate in enumerate(model)
                      if _equivalent(item, candidate, cache)] for item in source]
    if any(not candidates for candidates in possibilities):
        return None

    def assign(position, used, mapping):
        if position == len(source):
            return mapping
        for index in possibilities[position]:
            if index not in used:
                result = assign(position + 1, used | {index}, mapping + [index])
                if result is not None:
                    return result
        return None

    return assign(0, set(), [])


def _species_evidence(species, cache):
    prepared, state = _prepared_species(species, cache)
    return {
        "label": species.label,
        "multiplicity": state[0],
        "net_charge": state[1],
        "isotopic_composition": [
            {"element": element, "isotope": isotope, "count": count}
            for (element, isotope), count in state[2]
        ],
        "adjacency_list": species.molecule[0].to_adjacency_list(),
        "resonance_adjacency_lists": sorted(set(
            molecule.to_adjacency_list() for molecule in prepared.molecule
        )),
    }


def _occurrence_evidence(source, model, indexes, cache):
    from rmgpy.chemkin import get_species_identifier
    return [{"source_position": index, "model_position": model_index,
             "source_label": source[index].label,
             "model_label": get_species_identifier(model[model_index]),
             "source": _species_evidence(source[index], cache),
             "model": _species_evidence(model[model_index], cache)}
            for index, model_index in enumerate(indexes)]


def _load_source(library):
    try:
        from rmgpy.data.kinetics.database import KineticsDatabase
        database = KineticsDatabase()
        database.load_libraries(str(library.parent), libraries=[str(library)])
    except Exception as error:
        # A library that RMG could not interpret has not reached structural
        # matching. Preserve that distinction without repairing trusted inputs
        # or relaxing identity checks in response to a loader error.
        raise IdentityUnresolved("replacement_source_library_load_failed") from error
    if len(database.libraries) != 1:
        raise IdentityUnresolved("replacement_source_library_registration_ambiguous")
    registry_name, loaded = next(iter(database.libraries.items()))
    loaded_reactions = loaded.get_library_reactions()
    if len(loaded_reactions) != 1:
        raise IdentityUnresolved("replacement_source_reaction_count_not_one")
    return str(registry_name), loaded_reactions[0]


def _one_correspondence(request, reactions, records, cache):
    result = {"route_id": request.get("route_id"), "status": "unresolved",
              "blocking_reasons": [], "library": request.get("library"),
              "source_artifacts": request.get("source_artifacts"),
              "source_library": None, "source_equation": None}
    try:
        library, manifest = _check_sources(request)
        stereo_blockers = []
        try:
            _check_declared_stereo(manifest)
        except IdentityUnresolved as error:
            stereo_blockers.append(str(error))
        registry_name, source = _load_source(library)
        result.update(source_library=registry_name, source_equation=str(source))
        for item in source.reactants + source.products:
            try:
                _check_species_stereo(_prepared_species(item, cache)[0])
            except IdentityUnresolved as error:
                stereo_blockers.append(str(error))
        matches = []
        for index, candidate in enumerate(reactions):
            for direction, reactants, products in (
                    ("forward", candidate.reactants, candidate.products),
                    ("reverse", candidate.products, candidate.reactants)):
                left = _match_side(source.reactants, reactants, cache)
                if left is None:
                    continue
                right = _match_side(source.products, products, cache)
                if right is not None:
                    matches.append((index, direction, left, right))
        result["structural_candidates"] = [
            {"reaction_index": index, "applied_direction": direction,
             "equation": records[index].get("equation")}
            for index, direction, _, _ in matches
        ]
        if not matches:
            raise IdentityUnresolved("replacement_structural_match_not_found")
        if len(matches) != 1:
            raise IdentityUnresolved("replacement_structural_match_ambiguous")
        index, direction, left, right = matches[0]
        candidate, record = reactions[index], records[index]
        if bool(candidate.reversible) != bool(source.reversible):
            raise IdentityUnresolved("replacement_reversibility_mismatch")
        if direction == "reverse" and not bool(source.reversible):
            raise IdentityUnresolved("replacement_irreversible_orientation_mismatch")
        if str(getattr(candidate, "library", "") or "") != registry_name:
            raise IdentityUnresolved("replacement_source_library_provenance_mismatch")
        if record.get("source_library") != registry_name or record.get("reaction_index") != index:
            raise IdentityUnresolved("replacement_inspection_record_mismatch")
        model_reactants = candidate.reactants if direction == "forward" else candidate.products
        model_products = candidate.products if direction == "forward" else candidate.reactants
        for item in model_reactants + model_products:
            try:
                _check_species_stereo(_prepared_species(item, cache)[0])
            except IdentityUnresolved as error:
                stereo_blockers.append(str(error))
        if stereo_blockers:
            result["blocking_reasons"] = sorted(set(stereo_blockers))
            return result
        left_labels = record["reactant_labels"] if direction == "forward" else record["product_labels"]
        right_labels = record["product_labels"] if direction == "forward" else record["reactant_labels"]
        equation = "+".join(left_labels) + ("<=>" if source.reversible else "=>") + "+".join(right_labels)
        # Re-read all three trusted source files after RMG interpreted them.
        _check_sources(request)
        result.update(status="matched", matched_reaction_index=index, applied_direction=direction,
                      source_oriented_model_equation=equation,
                      species_correspondence={
                          "reactants": _occurrence_evidence(source.reactants, model_reactants, left, cache),
                          "products": _occurrence_evidence(source.products, model_products, right, cache),
                      })
    except IdentityUnresolved as error:
        result["blocking_reasons"] = [str(error)]
        if error.__cause__ is not None:
            cause = error.__cause__
            result["diagnostic_error"] = "%s: %s" % (type(cause).__name__, cause)
    except Exception as error:
        result["blocking_reasons"] = ["replacement_identity_diagnosis_failed"]
        result["diagnostic_error"] = "%s: %s" % (type(error).__name__, error)
    return result


def build_replacement_correspondences(payload, species, reactions, reaction_records):
    """Return source-bound identity evidence; unsuccessful diagnosis fails closed."""
    requests = payload.get("replacement_libraries") or []
    before = payload.get("_identity_mechanism_artifacts_before")
    # Direct callers must explicitly provide a pre-load snapshot.  Computing
    # both hashes after parsing would not prove what RMG actually read.
    report = {
        "schema_version": 1, "kind": "storca_replacement_correspondences",
        "policy_id": POLICY_ID, "mechanism_artifacts": before or {},
        "inspection_sha256": hashlib.sha256(json.dumps(
            reaction_records, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode("utf-8")).hexdigest(),
        "correspondences": [],
    }
    cache = {}
    report["correspondences"] = [
        _one_correspondence(request, reactions, reaction_records, cache) for request in requests
    ]
    try:
        if not before or mechanism_artifacts(payload) != before:
            raise IdentityUnresolved("mechanism_artifacts_changed_during_identity_inspection")
        # Also recheck sources after the entire operation, not merely after an
        # individual library.  Later loads cannot silently invalidate earlier claims.
        for request in requests:
            _check_sources(request)
    except Exception as error:
        reason = str(error) if isinstance(error, IdentityUnresolved) else "identity_artifact_recheck_failed"
        for entry in report["correspondences"]:
            entry["status"] = "unresolved"
            entry["blocking_reasons"] = sorted(set(entry["blocking_reasons"] + [reason]))
            for key in ("matched_reaction_index", "applied_direction", "source_oriented_model_equation",
                        "species_correspondence"):
                entry.pop(key, None)
    return report
