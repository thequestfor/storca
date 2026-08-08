#!/usr/bin/env python
"""Extract direct target/environment reactions from installed RMG libraries.

This module runs inside the selected RMG Conda environment and deliberately
has no STORCA imports.  Loading a complete combustion library into the model
edge can admit thousands of unrelated species or violate the screen's bounded
species constraints.  The extractor retains only library reactions whose
entire reactant side (or reversible product side) is composed of the declared
initial reactive species and includes the target.
"""

from __future__ import print_function

import argparse
import copy
import hashlib
import json
import os
from collections import OrderedDict


def _checksum(path):
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def _matches(species, molecule):
    return any(candidate.is_isomorphic(molecule) for candidate in species.molecule)


def _side_is_declared(side, target, declared):
    return bool(side) and any(_matches(species, target) for species in side) and all(
        any(_matches(species, molecule) for molecule in declared) for species in side
    )


def _species_key(species):
    smiles = [molecule.to_smiles() for molecule in species.molecule]
    return min(smiles) if smiles else species.label


def _reaction_key(reaction):
    left = tuple(sorted(_species_key(species) for species in reaction.reactants))
    right = tuple(sorted(_species_key(species) for species in reaction.products))
    if reaction.reversible and right < left:
        left, right = right, left
    return left, right, bool(reaction.reversible)


def _species_record(species, target):
    molecule = species.molecule[0]
    return {
        "label": "stability" if _matches(species, target) else species.label,
        "source_label": species.label,
        "smiles": molecule.to_smiles(),
        "multiplicity": int(molecule.multiplicity),
        "charge": int(sum(atom.charge for atom in molecule.atoms)),
    }


def _direction_record(reactants, products, target, source_library, source_label):
    reactant_records = [_species_record(species, target) for species in reactants]
    product_records = [_species_record(species, target) for species in products]

    def equation_side(records):
        counts = OrderedDict()
        for record in records:
            label = record["label"]
            counts[label] = counts.get(label, 0) + 1
        return "+".join(
            ((str(count) + " ") if count > 1 else "") + label
            for label, count in counts.items()
        )

    return {
        "reaction_equation": equation_side(reactant_records) + "=>" + equation_side(product_records),
        "reactants": reactant_records,
        "products": product_records,
        "source_library": source_library,
        "source_entry_label": source_label,
    }


def extract_subset(target_smiles, declared_smiles, library_names, output_dir):
    from rmgpy import settings
    from rmgpy.data.kinetics.database import KineticsDatabase
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy.molecule import Molecule

    target = Molecule().from_smiles(target_smiles)
    declared = [Molecule().from_smiles(smiles) for smiles in declared_smiles]
    if not any(target.is_isomorphic(molecule) for molecule in declared):
        declared.insert(0, target)

    database = KineticsDatabase()
    library_root = os.path.join(settings["database.directory"], "kinetics", "libraries")
    database.load_libraries(library_root, libraries=library_names)

    selected = []
    seen = set()
    for library_name in library_names:
        library = database.libraries[library_name]
        for entry in library.entries.values():
            reaction = entry.item
            direct_forward = _side_is_declared(reaction.reactants, target, declared)
            direct_reverse = reaction.reversible and _side_is_declared(reaction.products, target, declared)
            if not (direct_forward or direct_reverse):
                continue
            key = _reaction_key(reaction)
            if key in seen:
                continue
            seen.add(key)
            retained = copy.deepcopy(entry)
            provenance = "Source RMG database library: {0}; source entry index: {1}".format(
                library_name, entry.index,
            )
            retained.long_desc = (retained.long_desc.strip() + "\n" + provenance).strip()
            selected.append((library_name, retained))

    os.makedirs(output_dir, exist_ok=True)
    reaction_file = os.path.join(output_dir, "reactions.py")
    dictionary_file = os.path.join(output_dir, "dictionary.txt")
    subset = KineticsLibrary(
        label="storca_reference_subset",
        name="STORCA direct reference subset",
        short_desc="Direct reactions extracted from installed RMG reference libraries",
        long_desc="Only reactions directly reachable from the declared initial species are retained.",
    )
    subset.entries = OrderedDict((index + 1, entry) for index, (_, entry) in enumerate(selected))
    subset.save(reaction_file)
    # ``save`` also writes dictionary.txt, but retain an explicit compatibility
    # call for RMG releases whose library saver does not.
    if not os.path.isfile(dictionary_file):
        subset.save_dictionary(dictionary_file)

    manifest = {
        "schema_version": 1,
        "kind": "storca_rmg_reference_subset",
        "target_smiles": target_smiles,
        "declared_reactive_smiles": declared_smiles,
        "source_libraries": library_names,
        "selected_reaction_count": len(selected),
        "selected_reactions": [
            {
                "source_library": library_name,
                "label": entry.label,
                "candidate_directions": [
                    *(
                        [_direction_record(entry.item.reactants, entry.item.products, target, library_name, entry.label)]
                        if _side_is_declared(entry.item.reactants, target, declared) else []
                    ),
                    *(
                        [_direction_record(entry.item.products, entry.item.reactants, target, library_name, entry.label)]
                        if entry.item.reversible and _side_is_declared(entry.item.products, target, declared) else []
                    ),
                ],
            }
            for library_name, entry in selected
        ],
        "library": os.path.abspath(output_dir),
        "checksums": {
            "reactions.py": _checksum(reaction_file),
            "dictionary.txt": _checksum(dictionary_file),
        },
        "interpretation": (
            "Curated RMG database kinetics retained as candidate evidence; "
            "not an ORCA-verified rate or stability conclusion."
        ),
    }
    manifest_file = os.path.join(output_dir, "reference-subset.json")
    with open(manifest_file, "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--declared-smiles", action="append", default=[])
    parser.add_argument("--library", action="append", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    manifest = extract_subset(
        args.target_smiles,
        args.declared_smiles or [args.target_smiles],
        args.library,
        args.output,
    )
    print("STORCA_RMG_REFERENCE_SUBSET=" + json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
