"""Route normalization, selection, and verification preparation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
import re
from typing import Any


@dataclass(frozen=True)
class RouteSpec:
    """Engine-neutral description of one decomposition route.

    ``atom_mapping`` maps compact, concatenated heavy-atom indices on the
    reactant side to the corresponding compact heavy-atom indices on the
    product side.  Every hydrogen nucleus—including a standalone ``[H]``,
    ``[H+]``, or ``[H-]`` fragment—is mapped separately by the explicit-H
    layer.  The heavy-atom mapping may be omitted: the generic path engine
    attempts a bounded graph-edit mapping and fails closed when the best
    correspondence is ambiguous.  Fragment charges and multiplicities retain
    the electronic-state contract used to select a common total-spin surface;
    they are never silently replaced by neutral singlets.
    """

    route_id: str
    source: str
    parent_smiles: str
    reactant_smiles: tuple[str, ...]
    product_smiles: tuple[str, ...]
    reactant_labels: tuple[str, ...] = ()
    product_labels: tuple[str, ...] = ()
    reaction_equation: str | None = None
    charge: int = 0
    multiplicity: int = 1
    reactant_charges: tuple[int, ...] = ()
    product_charges: tuple[int, ...] = ()
    reactant_multiplicities: tuple[int, ...] = ()
    product_multiplicities: tuple[int, ...] = ()
    surface_multiplicities: tuple[int, ...] = ()
    broken_bonds: tuple[tuple[int, int], ...] = ()
    formed_bonds: tuple[tuple[int, int], ...] = ()
    atom_mapping: tuple[tuple[int, int], ...] = ()
    source_stage: str | None = None
    kinetic_relevance: str | None = None
    estimated_t95_seconds: float | None = None
    target_loss_fraction: float | None = None
    score: float = 0.0
    evidence_level: str = "candidate"
    protocol: str = "unsupported_atom_mapping"
    limitation: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _radical_multiplicity(smiles: str) -> int:
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid route species SMILES: {smiles}")
    return max(1, sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()) + 1)


def _homolysis_spec(parent_smiles: str, route: dict, *, source: str, stage: str | None,
                    t95: float | None, target_loss: float | None, index: int) -> RouteSpec:
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(parent_smiles)
    if molecule is None:
        raise ValueError(f"Invalid parent SMILES: {parent_smiles}")
    products = tuple(route.get("fragment_radical_smiles") or route.get("product_smiles") or ())
    bond = tuple(int(value) for value in route.get("bond_atom_indices", ()))
    if len(products) < 2 or len(bond) != 2:
        raise ValueError("A homolysis route needs radical products and one retained bond index pair")
    mapping = tuple(
        (compact_index, compact_index)
        for compact_index, _ in enumerate(
            atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() != 1
        )
    )
    accessibility = float(route.get("accessibility") or 0.0)
    relevance = route.get("kinetic_relevance")
    evidence = "kinetically_relevant" if relevance == "kinetically_relevant_candidate" or t95 is not None else "candidate"
    score = (1e12 / max(t95, 1e-12) if t95 is not None else 0.0) + accessibility
    base_route_id = str(route.get("route_id") or f"{source}-route-{index + 1}")
    return RouteSpec(
        route_id=f"{stage}:{base_route_id}" if stage else base_route_id,
        source=source,
        parent_smiles=parent_smiles,
        reactant_smiles=(parent_smiles,),
        product_smiles=products,
        reactant_labels=("parent",),
        product_labels=tuple(f"product-{index + 1}" for index in range(len(products))),
        product_multiplicities=tuple(_radical_multiplicity(value) for value in products),
        broken_bonds=((bond[0], bond[1]),),
        atom_mapping=mapping,
        source_stage=stage,
        kinetic_relevance=relevance,
        estimated_t95_seconds=t95,
        target_loss_fraction=target_loss,
        score=score,
        evidence_level=evidence,
        protocol="relaxed_dissociation_scan",
        limitation=(
            "The radical pair is followed on one declared total-spin surface. "
            "A relaxed ground-state scan is not excited-state photodynamics."
        ),
    )


def _rmg_spec(parent_smiles: str, route: dict, *, stage: str | None,
              t95: float | None, target_loss: float | None, index: int) -> RouteSpec:
    """Normalize an RMG route without inventing missing molecular mappings."""
    relevance = route.get("kinetic_relevance")
    evidence = "kinetically_relevant" if relevance == "kinetically_relevant_candidate" else "candidate"
    # Newer evidence producers may retain these explicitly.  Legacy reports
    # contain Chemkin labels only and are deliberately non-executable.
    resolved = route.get("resolved_endpoints") or {}
    reactant_records = resolved.get("reactants") or []
    product_records = resolved.get("products") or []
    reactants = tuple(route.get("reactant_smiles") or [item.get("smiles") for item in reactant_records if item.get("smiles")])
    products = tuple(route.get("product_smiles") or [item.get("smiles") for item in product_records if item.get("smiles")])
    mapping = tuple(tuple(pair) for pair in route.get("atom_mapping", ()))
    broken = tuple(tuple(pair) for pair in route.get("broken_bonds", ()))
    formed = tuple(tuple(pair) for pair in route.get("formed_bonds", ()))
    reactant_charges = tuple(int(value) for value in (
        route.get("reactant_charges") or [item.get("charge", 0) for item in reactant_records]
    ))
    product_charges = tuple(int(value) for value in (
        route.get("product_charges") or [item.get("charge", 0) for item in product_records]
    ))
    reactant_multiplicities = tuple(int(value) for value in (
        route.get("reactant_multiplicities") or
        [item["multiplicity"] for item in reactant_records if item.get("multiplicity") is not None]
    ))
    product_multiplicities = tuple(int(value) for value in (
        route.get("product_multiplicities") or
        [item["multiplicity"] for item in product_records if item.get("multiplicity") is not None]
    ))
    reactant_labels = tuple(
        str(value) for value in (
            route.get("reactant_labels") or [item.get("label") for item in reactant_records if item.get("label")]
        )
    )
    product_labels = tuple(
        str(value) for value in (
            route.get("product_labels") or [item.get("label") for item in product_records if item.get("label")]
        )
    )
    declared_multiplicity = route.get("multiplicity", route.get("total_multiplicity"))
    # Zero means no total surface was declared.  The geometry engine will fail
    # closed instead of quietly turning a radical route into a singlet.
    multiplicity = int(declared_multiplicity) if declared_multiplicity is not None else 0
    if "charge" in route:
        total_charge = int(route["charge"])
    elif reactant_charges and product_charges and sum(reactant_charges) == sum(product_charges):
        total_charge = sum(reactant_charges)
    else:
        total_charge = 0
    executable = bool(reactants and products)
    if executable:
        protocol = "generic_endpoint_path"
        limitation = None
    else:
        protocol = "unsupported_atom_mapping"
        limitation = (
            "The retained RMG route does not contain a verified atom mapping and resolved endpoint SMILES; "
            "STORCA will not guess them from Chemkin labels."
        )
    base_route_id = str(route.get("route_id") or f"rmg-route-{index + 1}")
    # A system-level target retention crossing cannot be assigned to every
    # candidate reaction.  Use one only when future flux analysis resolves it.
    route_t95 = _float_or_none(route.get("route_t95_seconds"))
    return RouteSpec(
        route_id=f"{stage}:{base_route_id}" if stage else base_route_id,
        source="rmg",
        parent_smiles=parent_smiles,
        reactant_smiles=reactants or (parent_smiles,),
        product_smiles=products,
        reactant_labels=reactant_labels,
        product_labels=product_labels,
        reaction_equation=route.get("reaction_equation"),
        charge=total_charge,
        multiplicity=multiplicity,
        reactant_charges=reactant_charges,
        product_charges=product_charges,
        reactant_multiplicities=reactant_multiplicities,
        product_multiplicities=product_multiplicities,
        surface_multiplicities=tuple(int(value) for value in route.get("surface_multiplicities", ())),
        broken_bonds=broken,
        formed_bonds=formed,
        atom_mapping=mapping,
        source_stage=stage,
        kinetic_relevance=relevance,
        estimated_t95_seconds=route_t95,
        target_loss_fraction=target_loss,
        score=(1e12 / max(route_t95, 1e-12) if route_t95 is not None else 0.0)
        + (1.0 if relevance == "kinetically_relevant_candidate" else 0.0),
        evidence_level=evidence,
        protocol=protocol,
        limitation=limitation,
    )


def route_spec_from_rmg_route(parent_smiles: str, route: dict, *, stage: str | None = None) -> RouteSpec:
    """Public normalization entry point for one resolved RMG route record."""
    return _rmg_spec(
        parent_smiles,
        route,
        stage=stage,
        t95=_float_or_none(route.get("route_t95_seconds")),
        target_loss=_float_or_none(route.get("target_loss_fraction")),
        index=int(route.get("route_index", 0)),
    )


def _float_or_none(value: Any) -> float | None:
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def _routes_from_stability(report: dict, *, stage: str | None = None,
                           t95: float | None = None, parent_smiles: str | None = None) -> list[RouteSpec]:
    parent = parent_smiles or report.get("smiles") or report.get("provenance", {}).get("smiles")
    if not parent:
        return []
    rmg = report.get("rmg_evidence", {})
    target_loss = _float_or_none((rmg.get("solver_profile") or {}).get("target_loss_fraction"))
    return [
        _rmg_spec(parent, route, stage=stage, t95=t95, target_loss=target_loss, index=index)
        for index, route in enumerate(rmg.get("candidate_routes", []))
    ]


def discover_routes(report_path: Path) -> tuple[dict, list[RouteSpec]]:
    """Read a stability, ladder, or computational-light result."""
    report_path = Path(report_path)
    if report_path.is_dir():
        candidates = [
            report_path / "stability-ladder.json",
            report_path / "stability.json",
            report_path / "computational-light.json",
        ]
        report_path = next((path for path in candidates if path.is_file()), report_path)
    if not report_path.is_file():
        raise FileNotFoundError(f"No supported result JSON found at {report_path}")
    report = json.loads(report_path.read_text())
    kind = report.get("kind")
    routes: list[RouteSpec] = []
    if kind == "condition_ladder_stability_screen":
        parent = report["smiles"]
        first_t95 = report.get("first_t95_stage")
        for stage in report.get("stages", []):
            stage_id = stage.get("id")
            t95 = _float_or_none(
                (stage.get("kinetic_lifetime") or {}).get("estimated_time_to_retention_seconds")
                or stage.get("estimated_time_to_retention_seconds")
            )
            for index, route in enumerate((stage.get("intrinsic_initiation") or {}).get("candidates", [])):
                routes.append(_homolysis_spec(
                    parent, route, source="intrinsic_initiation", stage=stage_id,
                    t95=t95, target_loss=None, index=index,
                ))
            routes.extend(_routes_from_stability(stage, stage=stage_id, t95=t95, parent_smiles=parent))
            light = stage.get("computational_light") or {}
            central = ((stage.get("computational_light_kinetics") or {}).get("profiles") or {}).get("central", {})
            light_t95 = _float_or_none(central.get("estimated_time_to_retention_seconds")) or t95
            light_loss = _float_or_none(central.get("target_loss_fraction"))
            for index, route in enumerate(light.get("candidate_routes", [])):
                routes.append(_homolysis_spec(parent, route, source="computational_light",
                                               stage=stage_id, t95=light_t95,
                                               target_loss=light_loss, index=index))
        # A route in the actual first-t95 stage wins ties deterministically.
        routes = [
            RouteSpec(**{**route.as_dict(), "score": route.score + (1e6 if route.source_stage == first_t95 else 0.0)})
            for route in routes
        ]
    elif kind == "computational_light_model":
        for index, route in enumerate(report.get("candidate_routes", [])):
            routes.append(_homolysis_spec(report["smiles"], route, source="computational_light",
                                           stage=None, t95=None, target_loss=None, index=index))
    else:
        metadata_path = report_path.parent / "metadata.json"
        metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
        routes.extend(_routes_from_stability(
            report,
            parent_smiles=report.get("smiles") or metadata.get("smiles"),
        ))
    return report, routes


def select_explanatory_route(routes: list[RouteSpec], route_id: str | None = None) -> RouteSpec:
    if not routes:
        raise ValueError("The retained report contains no resolved decomposition candidates")
    if route_id is not None:
        matches = [route for route in routes if route.route_id == route_id]
        if not matches and re.fullmatch(r"\d+", route_id):
            index = int(route_id)
            if 0 <= index < len(routes):
                return routes[index]
        if len(matches) != 1:
            raise ValueError(f"Route {route_id!r} was not found or is ambiguous")
        return matches[0]
    executable = [route for route in routes if route.protocol != "unsupported_atom_mapping"]
    pool = executable or routes
    return max(pool, key=lambda route: (route.score, route.route_id))


def prepare_explanation(report_path: Path, route_id: str | None = None) -> dict:
    report, routes = discover_routes(report_path)
    selected = select_explanatory_route(routes, route_id)
    return {
        "schema_version": 2,
        "kind": "decomposition_explanation",
        "status": "prepared" if selected.protocol != "unsupported_atom_mapping" else "unsupported_route_physics",
        "selected_route": selected.as_dict(),
        "candidate_routes": [route.as_dict() for route in routes],
        "condition_contract": report.get("condition_contract"),
        "limitations": [selected.limitation] if selected.limitation else [],
    }


def prepare_route_verification(run_dir: Path, route_index: int) -> dict:
    """Retain the original explicit RMG verification dossier interface."""
    run_dir = Path(run_dir)
    screen = json.loads((run_dir / "stability.json").read_text())
    routes = screen.get("rmg_evidence", {}).get("candidate_routes", [])
    if route_index < 0 or route_index >= len(routes):
        raise ValueError(f"Route index must be between 0 and {len(routes) - 1}")
    route = routes[route_index]
    if route.get("kinetic_relevance") != "kinetically_relevant_candidate":
        raise ValueError(
            "Route is not kinetically relevant under the declared conditions; "
            "ORCA TS verification is not justified yet"
        )
    dossier = {
        "schema_version": 1,
        "kind": "orca_route_verification_preparation",
        "status": (route.get("orca_verification") or {}).get("status", "requires_endpoint_and_ts_guess"),
        "selected_route_index": route_index,
        "route": route,
        "resolved_endpoints": route.get("resolved_endpoints"),
        "verification_state": route.get("orca_verification"),
        "condition_contract": screen.get("condition_contract"),
        "required_steps": [
            "Generate and optimize reactant/product endpoint structures with ORCA.",
            "Provide or generate a chemically plausible transition-state guess.",
            "Optimize the TS and require exactly one imaginary frequency.",
            "Run IRC calculations toward both expected endpoints.",
            "Calculate the condition-specific rate and t95 from the verified pathway.",
        ],
    }
    path = run_dir / "route-verification.json"
    path.write_text(json.dumps(dossier, indent=2, sort_keys=True) + "\n")
    return {**dossier, "path": path}
