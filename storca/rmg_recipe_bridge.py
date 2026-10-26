"""Source-bound family atom traces, executed in the owning RMG environment.

No rate is estimated and no path is verified here. Unsupported families,
stereochemistry, or non-unique recipe outcomes remain discovery limitations.
"""

from collections import Counter
from pathlib import Path


POLICY_ID = "rmg-family-atom-trace-v1"


def trace_family_recipe(payload):
    from rdkit import Chem
    from rmgpy import settings
    from rmgpy.chemkin import get_species_identifier, load_chemkin_file
    from rmgpy.data.kinetics.database import KineticsDatabase
    from rmgpy.molecule import Molecule
    from rmgpy.molecule.converter import to_rdkit_mol
    import rmgpy
    import rmgpy.chemkin
    import rmgpy.data.kinetics.family
    import rmgpy.molecule.converter
    try:
        from .rmg_identity_bridge import _digest, mechanism_artifacts
    except ImportError:
        from rmg_identity_bridge import _digest, mechanism_artifacts

    request = payload["request"]
    before = mechanism_artifacts(payload)
    if before != payload.get("source_artifacts"):
        raise ValueError("recipe_mechanism_artifact_hash_mismatch")
    index = request.get("reaction_index")
    if type(index) is not int or index < 0:
        raise ValueError("recipe_reaction_index_invalid")
    _, reactions = load_chemkin_file(payload["chemkin"], dictionary_path=payload["dictionary"])
    if index >= len(reactions):
        raise ValueError("recipe_reaction_index_out_of_range")
    reaction = reactions[index]
    direction = request.get("direction")
    if direction not in ("forward", "reverse"):
        raise ValueError("recipe_reaction_direction_invalid")
    if direction == "reverse" and not reaction.reversible:
        raise ValueError("recipe_irreversible_reaction_cannot_be_reversed")
    native_sides = (reaction.reactants, reaction.products)
    if direction == "reverse":
        native_sides = native_sides[::-1]
    endpoints = request["endpoints"]
    molecules = []
    for side, native_species in zip(("reactants", "products"), native_sides):
        declared = endpoints[side]
        if Counter(item["label"] for item in declared) != Counter(
                get_species_identifier(item) for item in native_species):
            raise ValueError("recipe_endpoint_label_or_stoichiometry_mismatch")
        remaining = list(native_species)
        side_molecules = []
        for item in declared:
            parsed = Chem.MolFromSmiles(item["smiles"])
            if parsed is None or Chem.FindPotentialStereo(parsed):
                raise ValueError("recipe_endpoint_invalid_or_stereochemistry_unsupported")
            molecule = Molecule().from_smiles(item["smiles"])
            if (type(item["charge"]) is not int or type(item["multiplicity"]) is not int
                    or molecule.get_net_charge() != item["charge"]
                    or molecule.multiplicity != item["multiplicity"]):
                raise ValueError("recipe_endpoint_electronic_state_mismatch")
            matches = [i for i, species in enumerate(remaining)
                       if get_species_identifier(species) == item["label"]
                       and any(molecule.is_isomorphic(other, strict=True) for other in species.molecule)]
            if not matches:
                raise ValueError("recipe_endpoint_graph_or_isotope_mismatch")
            remaining.pop(matches[0])
            side_molecules.append(molecule)
        molecules.append(side_molecules)
    family_name = getattr(reaction, "family", None)
    base = {"schema_version": 1, "kind": "rmg_family_atom_trace", "policy_id": POLICY_ID,
            "request": request, "source_artifacts": before,
            "rate_claim_allowed": False, "path_claim_allowed": False}
    if not family_name:
        return dict(base, status="unresolved", reason="source_reaction_has_no_family_recipe")
    if request.get("reaction_family") and request["reaction_family"] != family_name:
        raise ValueError("recipe_source_reaction_family_mismatch")
    families_path = Path(settings["database.directory"]) / "kinetics" / "families"
    family_path = (families_path / family_name).resolve(strict=True)
    if family_path.parent != families_path.resolve(strict=True):
        raise ValueError("recipe_family_directory_invalid")
    # The loader reads groups/rules and optional forbidden structures. No
    # depositories or unrelated families are loaded. Retain all local .py
    # sources so the recipe and its matching groups cannot change unnoticed.
    family_sources = {str(path.resolve()): _digest(path) for path in sorted(family_path.glob("*.py"))}
    implementation_paths = [Path(__file__), Path(__file__).with_name("rmg_bridge.py"),
                            Path(__file__).with_name("rmg_identity_bridge.py"),
                            Path(rmgpy.chemkin.__file__), Path(rmgpy.data.kinetics.family.__file__),
                            Path(rmgpy.molecule.converter.__file__)]
    implementation_sources = {str(path.resolve()): _digest(path) for path in implementation_paths}
    database = KineticsDatabase()
    database.load_families(str(families_path), families=[family_name], depositories=[])
    reactants, products = molecules
    next_id = 1
    for molecule in reactants:
        for atom in molecule.atoms:
            atom.id = next_id
            next_id += 1
    generated = database.families[family_name].generate_reactions(
        reactants, products=products, prod_resonance=False, delete_labels=False,
    )
    if len(generated) > 64:
        raise ValueError("recipe_trace_count_budget_exhausted")

    def mapped(molecule):
        rd_molecule, indices = to_rdkit_mol(
            molecule, remove_h=False, return_mapping=True, save_order=True,
        )
        for atom, atom_index in indices.items():
            if atom.id < 1:
                raise ValueError("recipe_atom_identity_not_preserved")
            rd_molecule.GetAtomWithIdx(atom_index).SetAtomMapNum(atom.id)
        return {"smiles": Chem.MolToSmiles(rd_molecule, isomericSmiles=True),
                "charge": molecule.get_net_charge(), "multiplicity": molecule.multiplicity}

    traces = []
    for generated_reaction in generated:
        left, right = generated_reaction.reactants, generated_reaction.products
        if generated_reaction.is_forward is False:
            left, right = right, left
        elif generated_reaction.is_forward is not True:
            raise ValueError("recipe_family_direction_not_established")
        left_ids = [atom.id for molecule in left for atom in molecule.atoms]
        right_ids = [atom.id for molecule in right for atom in molecule.atoms]
        if (len(set(left_ids)) != len(left_ids) or set(left_ids) != set(right_ids)
                or len(right_ids) != len(left_ids) or set(left_ids) != set(range(1, next_id))):
            raise ValueError("recipe_all_nuclei_not_bijective")
        trace = {"reactants": [mapped(item) for item in left],
                 "products": [mapped(item) for item in right]}
        if trace not in traces:
            traces.append(trace)
    if mechanism_artifacts(payload) != before or any(_digest(path) != digest
            for path, digest in {**family_sources, **implementation_sources}.items()):
        raise ValueError("recipe_sources_changed_during_trace")
    return dict(base, status="completed" if traces else "unresolved",
                reason=None if traces else "source_recipe_does_not_reproduce_declared_endpoints",
                reaction_family=family_name, family_sources=family_sources,
                implementation_sources=implementation_sources, traces=traces,
                rmg_version=getattr(rmgpy, "__version__", None))
