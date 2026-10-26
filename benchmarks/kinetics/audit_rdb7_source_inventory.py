"""Read-only source audit; generates an unresolved reference inventory, not bounds."""
import ast
import csv
from decimal import Decimal
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
import tarfile
from rdkit import Chem, rdBase


def sha(data):
    return hashlib.sha256(data).hexdigest()


def require(condition, detail):
    if not condition:
        raise ValueError(f"Reference audit failed: {detail}")


def molecule_identity(smiles):
    parameters = Chem.SmilesParserParams()
    parameters.removeHs = False
    molecule = Chem.MolFromSmiles(smiles, parameters)
    if molecule is None:
        raise ValueError(f"Invalid source SMILES: {smiles}")
    mapping = sorted(atom.GetAtomMapNum() for atom in molecule.GetAtoms())
    if any(i <= 0 for i in mapping) or len(mapping) != len(set(mapping)):
        raise ValueError("Source mapping is not positive and bijective")
    atom_map = sorted((atom.GetAtomMapNum(), atom.GetAtomicNum()) for atom in molecule.GetAtoms())
    charge = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    radical_count = sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms())
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    canonical = Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
    fragments = sorted(canonical.split("."))
    return {"canonical_fragments": fragments, "charge": charge,
            "rdkit_radical_electron_count": radical_count, "atom_map": atom_map}


def literal_assignments(text):
    result = {}
    for node in ast.parse(text).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                result[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                pass
    return result


def conformer_record(text):
    records = []
    for node in ast.parse(text).body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name) and node.value.func.id == "conformer":
            data = {}
            for keyword in node.value.keywords:
                if keyword.arg in ("label", "E0", "spin_multiplicity", "optical_isomers", "frequency"):
                    data[keyword.arg] = ast.literal_eval(keyword.value)
                if keyword.arg == "modes":
                    data["mode_models"] = []
                    for mode in keyword.value.elts:
                        if not isinstance(mode, ast.Call) or not isinstance(mode.func, ast.Name):
                            raise ValueError("Unexpected cached mode syntax")
                        data["mode_models"].append(mode.func.id)
                        if mode.func.id in ("LinearRotor", "NonlinearRotor"):
                            require("external_rotor_symmetry" not in data, "one external rotor per conformer")
                            values = [key for key in mode.keywords if key.arg == "symmetry"]
                            require(len(values) == 1, "explicit cached external rotor symmetry")
                            data["external_rotor_symmetry"] = ast.literal_eval(values[0].value)
                        if mode.func.id == "HarmonicOscillator":
                            require("positive_frequencies" not in data, "one cached harmonic oscillator")
                            values = [key for key in mode.keywords if key.arg == "frequencies"]
                            require(len(values) == 1, "explicit cached harmonic frequencies")
                            data["positive_frequencies"] = ast.literal_eval(values[0].value)
            records.append(data)
    if len(records) != 1:
        raise ValueError("Expected one cached conformer")
    return records[0]


def forward_rate_rows(output, index):
    tables = re.findall(r"^#\s+Temp\.\s+k \(TST\)\s+Tunneling\s+k \(TST\+T\)\s+Units\n((?:#[^\n]*\n)+)", output, re.M)
    require(len(tables) == 1, (index, "one explicitly labeled forward table"))
    rows = re.findall(r"^#\s+([\d.]+) K\s+(\S+)\s+(\S+)\s+(\S+)\s+(s\^-1)\s*$", tables[0], re.M)
    require(len(rows) == 50, (index, "50 forward source rate-table rows"))
    return rows


def author_recipe_record(archive_data):
    """Audit a pinned author-code deposit without importing or executing it."""
    require(hashlib.md5(archive_data).hexdigest() == "3063db4589a1038f99af8f4b7f6340e8", "official author-code checksum")
    require(sha(archive_data) == "3b104924953bd7d741bf43ff3ccdaa1163d2e1ca85b95963edaa821f50c24c40", "reviewed author-code SHA256")
    member = "reactants_products_ts_refined-1.0.0/scripts/run_arkane.py"
    with tarfile.open(fileobj=io.BytesIO(archive_data), mode="r:gz") as archive:
        data = archive.extractfile(member).read()
    require(sha(data) == "cb694ad72c480911d956b3f4f65554b4e2c0728d5de5ae02caffc7d8b2a0e857", "reviewed author rate script SHA256")
    tree = ast.parse(data.decode())
    assignments = [node for node in ast.walk(tree) if isinstance(node, ast.Assign)]
    three_parameter_values = [ast.literal_eval(node.value) for node in assignments
                             if any(isinstance(target, ast.Name) and target.id == "three_params" for target in node.targets)]
    require(three_parameter_values == [False], "author two-parameter fit declaration")
    scale_sources = [ast.unparse(node.value) for node in assignments
                     if any(isinstance(target, ast.Name) and target.id == "freq_scale_factor" for target in node.targets)]
    require(scale_sources == ["assign_frequency_scale_factor(LOT)"], "author dynamic scale source")
    defaults = {}
    methods = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Attribute) and node.func.attr == "add_argument" and node.args:
            flag = ast.literal_eval(node.args[0])
            if flag in ("--T_min", "--T_max", "--T_count", "--method", "--basis", "--software"):
                default = [keyword.value for keyword in node.keywords if keyword.arg == "default"]
                require(len(default) == 1, (flag, "author argument default"))
                defaults[flag] = ast.literal_eval(default[0])
        if isinstance(node.func, ast.Attribute) and node.func.attr == "compute_high_p_rate_coefficient":
            methods.append(node.func.attr)
        if isinstance(node.func, ast.Name) and node.func.id == "ArkaneAdapter":
            bac = [keyword.value for keyword in node.keywords if keyword.arg == "bac_type"]
            require(len(bac) == 1 and ast.literal_eval(bac[0]) is None, "author BAC disabled for rate/thermo adapter")
    require(methods == ["compute_high_p_rate_coefficient"], "author high-pressure rate producer")
    require({flag: defaults[flag] for flag in ("--T_min", "--T_max", "--T_count")} == {"--T_min": 300, "--T_max": 2000, "--T_count": 50}, "author temperature defaults")
    return {
        "source_url": "https://zenodo.org/api/records/5652085/files/reactants_products_ts_refined-1.0.0.tar.gz/content",
        "source_DOI": "10.5281/zenodo.5652085", "source_version": "v1.0.0",
        "archive_sha256": sha(archive_data), "member": member, "member_sha256": sha(data),
        "read_not_executed": True, "three_params": False, "bac_type": None,
        "rate_method": methods[0], "argument_defaults": defaults,
        "frequency_scale_source": scale_sources[0], "numeric_frequency_scale": None,
        "exact_dependency_commits_verified": False, "historical_scale_application_resolved": False,
    }


def main():
    source_dir, destination = map(Path, sys.argv[1:])
    csv_data = (source_dir / "rdb7-rigid-v1.0.1.csv").read_bytes()
    archive_data = (source_dir / "rdb7-rigid-v1.0.1.tar.gz").read_bytes()
    temperatures_data = (source_dir / "rdb7-temperatures-v1.0.1.csv").read_bytes()
    author_archive = source_dir / "rdb7-authors-code-v1.0.0.tar.gz"
    author_recipe = author_recipe_record(author_archive.read_bytes()) if author_archive.is_file() else None
    require(sha(csv_data) == "7934d83ff7b9217142ee6b53499de281c416929b0916c75ccf625c5a804b9b99", "raw CSV hash")
    require(hashlib.md5(archive_data).hexdigest() == "c9578110caa2c07e6c8db2ee70fdb34e", "official archive checksum")
    require(sha(archive_data) == "35c088276a27f70ac914f4c66dab7b5020982f4c390987432ca79d5e5c2cbfb8", "reviewed archive SHA256")
    require(hashlib.md5(temperatures_data).hexdigest() == "56753d9e0e4f1c4fca908da91cf590fa", "official temperature checksum")
    require(sha(temperatures_data) == "11517f9165e0ef31d1544f01d1c4f5e29e6dff11399ffec7bde95bbdfe6e4f4a", "reviewed temperature SHA256")
    temp_rows = list(csv.reader(io.StringIO(temperatures_data.decode())))
    temperatures = {row[0]: float(row[1]) for row in temp_rows[1:]}
    source_rows = list(csv.DictReader(io.StringIO(csv_data.decode())))
    require(len(source_rows) == 105 and len({row["idx"] for row in source_rows}) == 105, "105 unique source indices")
    require(set(temperatures) == {f"T{i}" for i in range(50)} and all(math.isfinite(value) and value > 0 for value in temperatures.values()), "50 positive finite temperature labels")
    p1_pair = sorted([molecule_identity("[CH:1]1=[CH:2][CH2:3][CH2:4]1")["canonical_fragments"],
                      molecule_identity("[CH2:1]=[CH:2][CH:3]=[CH2:4]")["canonical_fragments"]])
    entries = []
    versions = set()
    max_rounding = 0.0
    archive = tarfile.open(fileobj=io.BytesIO(archive_data), mode="r:gz")
    members = {member.name: member for member in archive.getmembers()
               if member.isfile() and not Path(member.name).name.startswith("._")}
    def read(name):
        return archive.extractfile(members[name]).read()
    for source in source_rows:
        index = int(source["idx"])
        prefix = f"arkane_logs_rigid/rxn{index:06d}"
        reactant = molecule_identity(source["rsmi"])
        products = molecule_identity(source["psmi"])
        require(reactant["atom_map"] == products["atom_map"], (index, "atom-map conservation"))
        require(reactant["charge"] == products["charge"], (index, "charge conservation"))
        require(len(reactant["canonical_fragments"]) == 1, (index, "one reactant"))
        pair = sorted([reactant["canonical_fragments"], products["canonical_fragments"]])
        output_name = f"{prefix}/rxns/ts{index:06d}/arkane/output.py"
        output_data = read(output_name)
        output = output_data.decode()
        ts_conformer = conformer_record(output)
        positive, units = ts_conformer["positive_frequencies"]
        require(units == "cm^-1" and positive and all(math.isfinite(f) and f > 0 for f in positive), (index, "positive finite TS modes"))
        imaginary, imaginary_units = ts_conformer["frequency"]
        require(imaginary_units == "cm^-1" and math.isfinite(imaginary) and imaginary < 0, (index, "negative finite TS imaginary mode"))
        rate_rows = forward_rate_rows(output, index)
        for i, row in enumerate(rate_rows):
            temperature, plain_rate, tunneling, corrected_rate = map(float, row[:4])
            expected_temperature = temperatures[f"T{i}"]
            require(math.isclose(temperature, expected_temperature, rel_tol=5e-6), (index, i, "temperature matches retained precision"))
            rate = float(source[f"k(T{i})"])
            require(math.isfinite(rate) and rate > 0 and math.isfinite(corrected_rate) and corrected_rate > 0, (index, i, "positive finite rates"))
            # Source CSV has three/four significant figures; compare display precision,
            # not a physical rate-error tolerance. Normalizing integer trailing zeroes
            # gives a conservative last-displayed-nonzero-place interpretation.
            difference = abs(rate - corrected_rate) / corrected_rate
            csv_half_step = Decimal(10) ** Decimal(source[f"k(T{i})"]).normalize().as_tuple().exponent / 2
            log_half_step = Decimal(10) ** Decimal(row[3]).as_tuple().exponent / 2
            require(abs(Decimal(source[f"k(T{i})"]) - Decimal(row[3])) <= csv_half_step + log_half_step, (index, i, rate, corrected_rate, "rate display-rounding intervals do not overlap"))
            max_rounding = max(max_rounding, difference)
        inputs = []
        cached_states = []
        cached_statmech = []
        for name in sorted(members):
            if not name.startswith(prefix + "/"):
                continue
            if name.endswith("_arkane_input_no_BAC.py"):
                data = read(name)
                inputs.append({"member": name, "sha256": sha(data), "declared_settings": literal_assignments(data.decode())})
            elif name.endswith("/arkane/output.py"):
                data = read(name)
                conformer = conformer_record(data.decode())
                frequencies, units = conformer.get("positive_frequencies", ([], "cm^-1"))
                require(units == "cm^-1" and all(math.isfinite(f) and f > 0 for f in frequencies), (index, name, "cached finite positive harmonic modes"))
                cached_statmech.append({
                    "member": name, "sha256": sha(data),
                    "kind": "transition_state" if "/rxns/" in name else "species",
                    "label": conformer["label"], "E0": conformer["E0"],
                    "spin_multiplicity": conformer["spin_multiplicity"],
                    "optical_isomers": conformer["optical_isomers"],
                    "mode_models": conformer["mode_models"],
                    "external_rotor_symmetry": conformer.get("external_rotor_symmetry"),
                    "positive_harmonic_mode_count": len(frequencies),
                    "minimum_positive_frequency_cm_inverse": min(frequencies) if frequencies else None,
                    "imaginary_frequency": conformer.get("frequency"),
                    "frequency_scale_factor_applied": None,
                    "ZPE_scale_factor_applied": None,
                    "scale_application_status": "cached_values_available_unscaled_QM_logs_and_global_recipe_not_joined",
                })
            elif name.endswith(".yml"):
                data = read(name)
                body = data.decode()
                versions.add(re.search(r"^RMG_version: (.+)$", body, re.M)[1])
                charge = int(re.search(r"^charge: (-?\d+)$", body, re.M)[1])
                graph_multiplicity = int(re.search(r"^multiplicity: (\d+)$", body, re.M)[1])
                statmech_multiplicity = int(re.search(r"^  spin_multiplicity: (\d+)$", body, re.M)[1])
                cached_states.append({"member": name, "sha256": sha(data), "charge": charge,
                                      "graph_metadata_multiplicity": graph_multiplicity,
                                      "statmech_conformer_multiplicity": statmech_multiplicity,
                                      "source_graph_statmech_multiplicity_consistent": graph_multiplicity == statmech_multiplicity,
                                      "ground_state_verified": False})
        entries.append({
            "source_reaction_index": index,
            "source_reactant_atom_mapped_smiles": source["rsmi"],
            "source_product_atom_mapped_smiles": source["psmi"],
            "canonical_reactant_fragments": reactant["canonical_fragments"],
            "canonical_product_fragments": products["canonical_fragments"],
            "undirected_graph_pair_sha256": sha(json.dumps(pair, separators=(",", ":")).encode()),
            "molecularity_forward": 1,
            "product_count": len(products["canonical_fragments"]),
            "net_charge": reactant["charge"],
            "rdkit_endpoint_radical_electron_counts": [reactant["rdkit_radical_electron_count"], products["rdkit_radical_electron_count"]],
            "p1_cyclobutene_butadiene_pair_overlap": pair == p1_pair,
            "source_arkane_rate_member": output_name,
            "source_arkane_rate_member_sha256": sha(output_data),
            "source_forward_rate_units": "s^-1",
            "raw_rate_column_count_verified": 50,
            "source_TS_declared_multiplicity": ts_conformer["spin_multiplicity"],
            "source_TS_imaginary_frequency_cm_inverse": imaginary,
            "source_TS_positive_mode_count": len(positive),
            "source_TS_minimum_positive_frequency_cm_inverse": min(positive),
            "source_species_cached_states": cached_states,
            "source_cached_statmech_records": cached_statmech,
            "source_species_TS_input_records": inputs,
            "current_producer_HPL_branch_compatible": len(products["canonical_fragments"]) == 1,
            "current_producer_optical_isomer_policy_mismatch": any(record["optical_isomers"] != 1 for record in cached_statmech),
            "ordinary_closed_shell_graph_candidate_before_reference_and_method_audit": (
                len(products["canonical_fragments"]) == 1
                and not reactant["rdkit_radical_electron_count"]
                and not products["rdkit_radical_electron_count"]
                and all(state["source_graph_statmech_multiplicity_consistent"] for state in cached_states)
            ),
            "orbital_reference_classification": {"forward_barrier": None, "reverse_barrier": None, "reaction_energy": None},
            "orbital_reference_join_status": "unresolved_refined_spreadsheet_unavailable",
            "ground_state_and_mechanism_audit_status": "pending",
            "panel_overlap_audit_status": "only_P1_pair_checked_other_evaluation_identities_and_scaffolds_pending",
            "calibration_admitted": False,
        })
    pair_counts = {}
    for entry in entries:
        pair_counts.setdefault(entry["undirected_graph_pair_sha256"], []).append(entry["source_reaction_index"])
    inventory = {
        "schema_version": 2, "kind": "storca_unresolved_rdb7_rigid_reference_inventory",
        "audit_date": "2026-10-02", "status": "raw_source_inventory_verified_orbital_join_unresolved",
        "support_claim_allowed": False, "calibration_profile_admitted": False,
        "attribution": "Spiekermann, Pattanaik and Green, Scientific Data 9, 417 (2022), DOI10.1038/s41597-022-01529-6; Zenodo v1.0.1 DOI10.5281/zenodo.6618262, CC BY4.0. Inventory adds source hashes, structural identity keys, and unresolved audit dispositions; rates are not replaced or fitted.",
        "source_urls": {
            "raw_rates": "https://zenodo.org/api/records/6618262/files/ccsdtf12_dz_rigid.csv/content",
            "raw_arkane_archive": "https://zenodo.org/api/records/6618262/files/ccsdtf12_dz_rigid.tar.gz/content",
            "temperatures": "https://zenodo.org/api/records/6618262/files/arkane_temperatures.csv/content",
            "orbital_audit": "https://doi.org/10.1039/D5CP01181G"},
        "source_sha256": {"raw_rates": sha(csv_data), "raw_arkane_archive": sha(archive_data), "temperatures": sha(temperatures_data)},
        "software": {"identity_audit_RDKit": rdBase.rdkitVersion, "reference_cached_RMG_versions": sorted(versions)},
        "author_code_recipe": author_recipe,
        "temperatures_K": temperatures,
        "summary": {"raw_reaction_count": len(entries), "unresolved_orbital_classification_count": len(entries),
                    "admitted_calibration_reaction_count": 0,
                    "P1_cyclobutene_butadiene_pair_overlap_count": sum(entry["p1_cyclobutene_butadiene_pair_overlap"] for entry in entries),
                    "unique_undirected_graph_pair_count": len(pair_counts),
                    "duplicate_undirected_graph_pairs": [ids for ids in pair_counts.values() if len(ids) > 1],
                    "endpoint_graph_radical_electron_warning_indices": [entry["source_reaction_index"] for entry in entries if any(entry["rdkit_endpoint_radical_electron_counts"])],
                    "source_graph_statmech_multiplicity_inconsistency_indices": [entry["source_reaction_index"] for entry in entries if any(not state["source_graph_statmech_multiplicity_consistent"] for state in entry["source_species_cached_states"])],
                    "product_count_distribution": {str(count): sum(entry["product_count"] == count for entry in entries) for count in sorted({entry["product_count"] for entry in entries})},
                    "current_producer_HPL_branch_compatible_count_before_other_audits": sum(entry["current_producer_HPL_branch_compatible"] for entry in entries),
                    "ordinary_closed_shell_graph_candidate_count_before_reference_and_method_audit": sum(entry["ordinary_closed_shell_graph_candidate_before_reference_and_method_audit"] for entry in entries),
                    "reference_TS_optical_isomer_count_distribution": {str(count): sum(record["optical_isomers"] == count for entry in entries for record in entry["source_cached_statmech_records"] if record["kind"] == "transition_state") for count in sorted({record["optical_isomers"] for entry in entries for record in entry["source_cached_statmech_records"] if record["kind"] == "transition_state"})},
                    "current_producer_optical_isomer_policy_mismatch_indices": [entry["source_reaction_index"] for entry in entries if entry["current_producer_optical_isomer_policy_mismatch"]],
                    "maximum_CSV_vs_printed_Arkane_relative_rounding_difference": max_rounding},
        "limits": ["RDB7_refined.xlsx was not obtained; no row inferred easy from favorable rates or barriers.",
                   "Graph identity, declared multiplicity, positive cached modes, and rate-table consistency are source audits, not verification of ground state, TS connection, conformer completeness, or method accuracy.",
                   "Global Arkane calculation input is absent from the small archive. Frequency scaling application, path degeneracy, numerical electronic settings and reference electronic logs remain pending.",
                   "Cached rotor symmetry, optical-isomer count, mode models and E0 are retained, but cached scale-application fields remain null until historical source logs and the global dependency recipe are joined.",
                   "Current STORCA shape-compatible HPL comparison is limited to1->1 routes;1->2 uses pressure dependence and1->3 exceeds current endpoint cardinality. Graph candidates are not calibrated or admitted rows.",
                   "At row3474 the graph metadata multiplicity is3 but input/conformer multiplicity is1; do not treat graph multiplicity as the electronic state used to compute the reference rate.",
                   "P1 pair exclusion checks graph identities in both directions; independent scaffold grouping and other evaluation-identity exclusions remain pending.",
                   "Original forward rate model is gas high-pressure RRHO/Eckart. No measured physical uncertainty or B3LYP/def2-SVP residuals exist in this inventory."],
        "reactions": entries,
    }
    with destination.open("x") as stream:
        json.dump(inventory, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(inventory["summary"], indent=2))


if __name__ == "__main__":
    main()
