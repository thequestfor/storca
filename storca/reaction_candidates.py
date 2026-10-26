"""Small source-neutral reaction-candidate records.

The stability workflow still uses RMG routes as its default kinetic-model
representation.  This module deliberately does less: it preserves candidate
identity, discovery assertions, and screening evidence before a route has a
rate or belongs to an RMG model.  A candidate is never evidence that a
reaction exists at an observable kinetic rate.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field, replace
import hashlib
import json
import re
from typing import Any, Iterable


REACTION_CANDIDATE_SCHEMA_VERSION = 1


def _jsonable(value: Any) -> Any:
    """Return a deterministic JSON-compatible copy for retained evidence."""
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in sorted(value.items(), key=lambda item: str(item[0]))}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    return value


def canonical_smiles(smiles: str) -> str:
    """Canonicalize a graph identity without treating atom-map labels as chemistry.

    RDKit is intentionally imported lazily.  Candidate artifacts remain
    readable in a minimal STORCA environment even when a graph backend cannot
    run there.
    """
    value = str(smiles).strip()
    try:
        from rdkit import Chem

        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            return value
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    except Exception:
        # Retained artifacts should still deduplicate mapped and unmapped
        # copies when RDKit is unavailable.  Atom-map syntax is terminal in a
        # bracket atom, so this conservative fallback cannot alter charge or
        # isotope annotations.
        return re.sub(r":\d+(?=\])", "", value)


def _side_identity(smiles: Iterable[str]) -> tuple[tuple[str, int], ...]:
    counts = Counter(canonical_smiles(value) for value in smiles)
    return tuple(sorted(counts.items()))


def reaction_net_identity(reactants: Iterable[str], products: Iterable[str]) -> tuple:
    """Return a direction-independent net graph/stoichiometry identity.

    The direction used for verification remains on :class:`ReactionCandidate`.
    This identity only decides whether independent discovery backends proposed
    the same elementary connectivity change.
    """
    sides = (_side_identity(reactants), _side_identity(products))
    return tuple(sorted(sides))


def candidate_id_for(
    reactants: Iterable[str], products: Iterable[str], *, total_charge: int | None,
    multiplicity: int | None = None,
) -> str:
    payload = {
        "net_identity": reaction_net_identity(reactants, products),
        "total_charge": total_charge,
        "multiplicity": multiplicity,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return "reaction-" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def endpoint_state(smiles: str) -> dict[str, int | None]:
    """Read the minimal state contract available directly from a SMILES graph."""
    try:
        from rdkit import Chem

        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Invalid endpoint SMILES: {smiles}")
        radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms())
        # A species with several radical electrons can have multiple physical
        # couplings.  Leave it unresolved rather than inventing one state.
        multiplicity = radical_electrons + 1 if radical_electrons <= 1 else None
        return {
            "charge": int(Chem.GetFormalCharge(molecule)),
            "multiplicity": multiplicity,
        }
    except Exception:
        return {"charge": None, "multiplicity": None}


def _state_vector(
    smiles: tuple[str, ...], supplied: tuple[int, ...], key: str,
) -> tuple[int, ...]:
    if len(supplied) == len(smiles):
        return tuple(int(value) for value in supplied)
    values = [endpoint_state(value).get(key) for value in smiles]
    if any(value is None for value in values):
        return ()
    return tuple(int(value) for value in values if value is not None)


def _assertion_key(assertion: dict[str, Any]) -> str:
    return json.dumps(_jsonable(assertion), sort_keys=True, separators=(",", ":"), default=str)


@dataclass(frozen=True)
class ReactionCandidate:
    """A minimal, source-neutral elementary-reaction hypothesis.

    ``atom_mapping`` is an atom-map *identifier* relation when present.  The
    representation is retained in provenance because the existing route
    geometry layer expects compact heavy-atom indices, whereas a graph backend
    normally supplies mapped SMILES.  This prevents a map-number convention
    from being silently misused as a compact index convention.
    """

    candidate_id: str
    reactants: tuple[str, ...]
    products: tuple[str, ...]
    total_charge: int | None = None
    multiplicity: int | None = None
    atom_mapping: tuple[tuple[int, int], ...] = ()
    discovery_assertions: tuple[dict[str, Any], ...] = ()
    estimated_reaction_energy_kcal_mol: float | None = None
    estimated_barrier_kcal_mol: float | None = None
    screening_status: str = "candidate_unverified"
    provenance: dict[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    verification_status: str = "candidate_unverified"
    kinetic_relevance: str = "unknown"
    reactant_labels: tuple[str, ...] = ()
    product_labels: tuple[str, ...] = ()
    reactant_charges: tuple[int, ...] = ()
    product_charges: tuple[int, ...] = ()
    reactant_multiplicities: tuple[int, ...] = ()
    product_multiplicities: tuple[int, ...] = ()
    # Descriptive environment metadata belongs to the candidate evidence even
    # when the first discovery backend has no solvent/surface physics.  It is
    # never interpreted as a rate model without an explicit downstream lane.
    environment: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.candidate_id:
            raise ValueError("ReactionCandidate needs a stable candidate_id")
        if not self.reactants or not self.products:
            raise ValueError("ReactionCandidate needs at least one reactant and product")
        if any(not str(value).strip() for value in (*self.reactants, *self.products)):
            raise ValueError("ReactionCandidate endpoint SMILES must be nonempty")
        if self.multiplicity is not None and self.multiplicity < 1:
            raise ValueError("ReactionCandidate multiplicity must be positive when known")

    @classmethod
    def create(
        cls,
        reactants: Iterable[str],
        products: Iterable[str],
        *,
        total_charge: int | None = None,
        multiplicity: int | None = None,
        **kwargs: Any,
    ) -> "ReactionCandidate":
        reactant_tuple = tuple(str(value) for value in reactants)
        product_tuple = tuple(str(value) for value in products)
        return cls(
            candidate_id=candidate_id_for(
                reactant_tuple, product_tuple, total_charge=total_charge,
                multiplicity=multiplicity,
            ),
            reactants=reactant_tuple,
            products=product_tuple,
            total_charge=total_charge,
            multiplicity=multiplicity,
            **kwargs,
        )

    @property
    def net_identity(self) -> tuple:
        return reaction_net_identity(self.reactants, self.products)

    @property
    def discovery_sources(self) -> tuple[str, ...]:
        values = {
            str(item.get("source"))
            for item in self.discovery_assertions
            if isinstance(item, dict) and item.get("source")
        }
        return tuple(sorted(values))

    @property
    def reactant_stoichiometry(self) -> dict[str, int]:
        return dict(_side_identity(self.reactants))

    @property
    def product_stoichiometry(self) -> dict[str, int]:
        return dict(_side_identity(self.products))

    def with_assertions(self, assertions: Iterable[dict[str, Any]]) -> "ReactionCandidate":
        """Return a candidate retaining every distinct backend assertion."""
        records: list[dict[str, Any]] = []
        seen: set[str] = set()
        for item in (*self.discovery_assertions, *tuple(assertions)):
            if not isinstance(item, dict):
                continue
            retained = _jsonable(item)
            key = _assertion_key(retained)
            if key not in seen:
                records.append(retained)
                seen.add(key)
        return replace(self, discovery_assertions=tuple(records))

    def as_dict(self) -> dict[str, Any]:
        """Serialize without losing the distinction between a candidate and a rate."""
        return {
            "schema_version": REACTION_CANDIDATE_SCHEMA_VERSION,
            "kind": "storca_reaction_candidate",
            **asdict(self),
            "discovery_sources": list(self.discovery_sources),
            "reactant_stoichiometry": self.reactant_stoichiometry,
            "product_stoichiometry": self.product_stoichiometry,
            "net_identity": _jsonable(self.net_identity),
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "ReactionCandidate":
        """Restore a retained candidate artifact with its scientific identity intact.

        This is deliberately strict about the deterministic candidate identity.
        A shared discovery artifact may be reused by several condition lanes,
        but it must never be silently rebound to a different reaction graph,
        charge, or multiplicity contract.
        """
        if not isinstance(value, dict):
            raise TypeError("ReactionCandidate artifact must be a mapping")
        kind = value.get("kind")
        if kind not in {None, "storca_reaction_candidate"}:
            raise ValueError(f"Unsupported reaction-candidate artifact kind: {kind}")
        reactants = tuple(str(item) for item in value.get("reactants") or ())
        products = tuple(str(item) for item in value.get("products") or ())
        total_charge = value.get("total_charge")
        multiplicity = value.get("multiplicity")
        if total_charge is not None:
            total_charge = int(total_charge)
        if multiplicity is not None:
            multiplicity = int(multiplicity)
        expected_id = candidate_id_for(
            reactants, products, total_charge=total_charge, multiplicity=multiplicity,
        )
        retained_id = str(value.get("candidate_id") or expected_id)
        if retained_id != expected_id:
            raise ValueError(
                "Reaction-candidate artifact identity does not match its endpoint and state contract"
            )
        return cls(
            candidate_id=retained_id,
            reactants=reactants,
            products=products,
            total_charge=total_charge,
            multiplicity=multiplicity,
            atom_mapping=tuple(
                tuple(int(entry) for entry in pair)
                for pair in value.get("atom_mapping") or ()
            ),
            discovery_assertions=tuple(
                _jsonable(item) for item in value.get("discovery_assertions") or ()
                if isinstance(item, dict)
            ),
            estimated_reaction_energy_kcal_mol=(
                float(value["estimated_reaction_energy_kcal_mol"])
                if value.get("estimated_reaction_energy_kcal_mol") is not None else None
            ),
            estimated_barrier_kcal_mol=(
                float(value["estimated_barrier_kcal_mol"])
                if value.get("estimated_barrier_kcal_mol") is not None else None
            ),
            screening_status=str(value.get("screening_status") or "candidate_unverified"),
            provenance=_jsonable(dict(value.get("provenance") or {})),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            verification_status=str(value.get("verification_status") or "candidate_unverified"),
            kinetic_relevance=str(value.get("kinetic_relevance") or "unknown"),
            reactant_labels=tuple(str(item) for item in value.get("reactant_labels") or ()),
            product_labels=tuple(str(item) for item in value.get("product_labels") or ()),
            reactant_charges=tuple(int(item) for item in value.get("reactant_charges") or ()),
            product_charges=tuple(int(item) for item in value.get("product_charges") or ()),
            reactant_multiplicities=tuple(
                int(item) for item in value.get("reactant_multiplicities") or ()
            ),
            product_multiplicities=tuple(
                int(item) for item in value.get("product_multiplicities") or ()
            ),
            environment=_jsonable(dict(value.get("environment") or {})),
        )

    def as_route_dict(self, *, parent_smiles: str, target_label: str = "stability") -> dict[str, Any]:
        """Adapt one candidate to the existing generic ORCA route verifier.

        This is intentionally an adapter, not a second path-verification
        engine.  Graph-only candidates retain ``kinetic_relevance=unknown``;
        they cannot support an observable lifetime until a verified rate is
        injected and a full RMG rerun accepts it.
        """
        reactant_charges = _state_vector(self.reactants, self.reactant_charges, "charge")
        product_charges = _state_vector(self.products, self.product_charges, "charge")
        reactant_multiplicities = _state_vector(
            self.reactants, self.reactant_multiplicities, "multiplicity",
        )
        product_multiplicities = _state_vector(
            self.products, self.product_multiplicities, "multiplicity",
        )
        total_charge = self.total_charge
        if total_charge is None and reactant_charges and sum(reactant_charges) == sum(product_charges):
            total_charge = sum(reactant_charges)
        labels = self.reactant_labels or tuple(
            target_label if index == 0 else f"{self.candidate_id}-r{index + 1}"
            for index in range(len(self.reactants))
        )
        product_labels = self.product_labels or tuple(
            f"{self.candidate_id}-p{index + 1}" for index in range(len(self.products))
        )
        equation = " + ".join(labels) + " <=> " + " + ".join(product_labels)
        atom_mapping_representation = str((self.provenance or {}).get("atom_mapping_representation") or "")
        compact_mapping = (
            [list(pair) for pair in self.atom_mapping]
            if atom_mapping_representation == "compact_heavy_atom_indices" else []
        )
        resolved_reactants = [
            {
                "label": label,
                "smiles": smiles,
                "charge": reactant_charges[index] if len(reactant_charges) == len(self.reactants) else None,
                "multiplicity": (
                    reactant_multiplicities[index]
                    if len(reactant_multiplicities) == len(self.reactants) else None
                ),
            }
            for index, (label, smiles) in enumerate(zip(labels, self.reactants))
        ]
        resolved_products = [
            {
                "label": label,
                "smiles": smiles,
                "charge": product_charges[index] if len(product_charges) == len(self.products) else None,
                "multiplicity": (
                    product_multiplicities[index]
                    if len(product_multiplicities) == len(self.products) else None
                ),
            }
            for index, (label, smiles) in enumerate(zip(product_labels, self.products))
        ]
        return {
            "schema_version": REACTION_CANDIDATE_SCHEMA_VERSION,
            "kind": "storca_candidate_route_adapter",
            "native_candidate_id": self.candidate_id,
            "route_id": self.candidate_id,
            "source": "+".join(self.discovery_sources) or "unknown",
            "discovery_sources": list(self.discovery_sources),
            "discovery_assertions": _jsonable(self.discovery_assertions),
            "parent_smiles": parent_smiles,
            "reaction_equation": equation,
            "reactant_smiles": list(self.reactants),
            "product_smiles": list(self.products),
            "reactant_labels": list(labels),
            "product_labels": list(product_labels),
            "reactant_charges": list(reactant_charges),
            "product_charges": list(product_charges),
            "reactant_multiplicities": list(reactant_multiplicities),
            "product_multiplicities": list(product_multiplicities),
            "charge": total_charge,
            # Zero deliberately asks route_geometry to select a unique common
            # surface, or fail closed if several surfaces remain possible.
            "multiplicity": self.multiplicity or 0,
            "atom_mapping": compact_mapping,
            "resolved_endpoints": {
                "reactants": resolved_reactants,
                "products": resolved_products,
            },
            "atom_balance": {"status": "passed", "source": "graph_candidate_inventory_gate"},
            "orca_verification": {
                "status": "ready_for_bounded_generic_mapping",
                "can_execute": True,
                "atom_mapping_resolved": bool(self.atom_mapping),
                "reason": "Candidate is a reaction hypothesis; ORCA path, TS, frequency, and IRC evidence remain required.",
            },
            "estimated_reaction_energy_kcal_mol": self.estimated_reaction_energy_kcal_mol,
            "estimated_barrier_kcal_mol": self.estimated_barrier_kcal_mol,
            "screening_status": self.screening_status,
            "verification_status": self.verification_status,
            "kinetic_relevance": self.kinetic_relevance,
            "environment": _jsonable(self.environment),
            "warnings": list(self.warnings),
            "provenance": _jsonable(self.provenance),
        }


def candidate_from_rmg_route(parent_smiles: str, route: dict) -> ReactionCandidate | None:
    """Adapt one resolved RMG route without changing RMG's route semantics."""
    resolved = route.get("resolved_endpoints") or {}
    reactants = tuple(route.get("reactant_smiles") or [
        item.get("smiles") for item in resolved.get("reactants") or [] if item.get("smiles")
    ])
    products = tuple(route.get("product_smiles") or [
        item.get("smiles") for item in resolved.get("products") or [] if item.get("smiles")
    ])
    if not reactants or not products:
        return None
    reactant_charges = tuple(int(value) for value in route.get("reactant_charges") or [
        item.get("charge") for item in resolved.get("reactants") or [] if item.get("charge") is not None
    ])
    product_charges = tuple(int(value) for value in route.get("product_charges") or [
        item.get("charge") for item in resolved.get("products") or [] if item.get("charge") is not None
    ])
    total_charge = route.get("charge")
    if total_charge is None and len(reactant_charges) == len(reactants) and sum(reactant_charges) == sum(product_charges):
        total_charge = sum(reactant_charges)
    multiplicity = route.get("multiplicity")
    if multiplicity is not None and int(multiplicity) < 1:
        multiplicity = None
    assertion = {
        "source": "rmg",
        "route_id": route.get("route_id"),
        "stable_route_identity": route.get("stable_route_identity"),
        "reaction_equation": route.get("reaction_equation"),
        "kinetic_relevance": route.get("kinetic_relevance"),
        "retained_rate_evidence": _jsonable(route.get("retained_rate_evidence") or {}),
    }
    return ReactionCandidate.create(
        reactants,
        products,
        total_charge=int(total_charge) if total_charge is not None else None,
        multiplicity=int(multiplicity) if multiplicity is not None else None,
        atom_mapping=tuple(tuple(int(value) for value in pair) for pair in route.get("atom_mapping") or ()),
        discovery_assertions=(assertion,),
        screening_status="candidate_unverified",
        provenance={
            "rmg_route_id": route.get("route_id"),
            "rmg_reaction_equation": route.get("reaction_equation"),
            "atom_mapping_representation": "compact_heavy_atom_indices",
            "parent_smiles": parent_smiles,
        },
        verification_status=str((route.get("orca_verification") or {}).get("status") or "candidate_unverified"),
        kinetic_relevance=str(route.get("kinetic_relevance") or "unknown"),
        reactant_labels=tuple(str(value) for value in route.get("reactant_labels") or [
            item.get("label") for item in resolved.get("reactants") or [] if item.get("label")
        ]),
        product_labels=tuple(str(value) for value in route.get("product_labels") or [
            item.get("label") for item in resolved.get("products") or [] if item.get("label")
        ]),
        reactant_charges=reactant_charges,
        product_charges=product_charges,
        reactant_multiplicities=tuple(int(value) for value in route.get("reactant_multiplicities") or [
            item.get("multiplicity") for item in resolved.get("reactants") or [] if item.get("multiplicity") is not None
        ]),
        product_multiplicities=tuple(int(value) for value in route.get("product_multiplicities") or [
            item.get("multiplicity") for item in resolved.get("products") or [] if item.get("multiplicity") is not None
        ]),
    )


def chemically_equivalent(left: ReactionCandidate, right: ReactionCandidate) -> bool:
    """Conservative cross-backend duplicate test.

    Atom maps improve later path verification but RMG normally does not retain
    one.  Thus map absence cannot prevent a graph/RMG merge when graph,
    stoichiometry, and known state contracts agree.
    """
    if left.net_identity != right.net_identity:
        return False
    if left.total_charge is not None and right.total_charge is not None and left.total_charge != right.total_charge:
        return False
    if left.multiplicity is not None and right.multiplicity is not None and left.multiplicity != right.multiplicity:
        return False
    return True


def _screening_priority(status: str) -> int:
    return {
        "candidate_unverified": 5,
        "xtb_screened_out": 4,
        "xtb_failed": 3,
        "generated_but_pruned": 2,
        "not_generated": 1,
    }.get(status, 0)


def merge_reaction_candidates(candidates: Iterable[ReactionCandidate]) -> list[ReactionCandidate]:
    """Merge independent discovery assertions without discarding provenance."""
    merged: list[ReactionCandidate] = []
    for incoming in candidates:
        index = next((position for position, current in enumerate(merged)
                      if chemically_equivalent(current, incoming)), None)
        if index is None:
            merged.append(incoming)
            continue
        current = merged[index]
        primary = incoming if incoming.atom_mapping and not current.atom_mapping else current
        assertions = (*current.discovery_assertions, *incoming.discovery_assertions)
        warnings = tuple(sorted(set((*current.warnings, *incoming.warnings))))
        merged_provenance = {
            **_jsonable(primary.provenance),
            "merged_candidate_ids": sorted({current.candidate_id, incoming.candidate_id}),
            "merged_discovery_sources": sorted(set((*current.discovery_sources, *incoming.discovery_sources))),
        }
        merged[index] = replace(
            primary,
            candidate_id=candidate_id_for(
                primary.reactants, primary.products,
                total_charge=(primary.total_charge if primary.total_charge is not None else
                              (incoming.total_charge if incoming.total_charge is not None else current.total_charge)),
                multiplicity=(primary.multiplicity if primary.multiplicity is not None else incoming.multiplicity),
            ),
            total_charge=(primary.total_charge if primary.total_charge is not None else incoming.total_charge),
            multiplicity=(primary.multiplicity if primary.multiplicity is not None else incoming.multiplicity),
            discovery_assertions=(),
            provenance=merged_provenance,
            warnings=warnings,
            screening_status=(
                primary.screening_status
                if _screening_priority(primary.screening_status) >= _screening_priority(incoming.screening_status)
                else incoming.screening_status
            ),
        ).with_assertions(assertions)
    return merged


__all__ = [
    "REACTION_CANDIDATE_SCHEMA_VERSION", "ReactionCandidate", "candidate_from_rmg_route",
    "candidate_id_for", "canonical_smiles", "chemically_equivalent", "endpoint_state",
    "merge_reaction_candidates", "reaction_net_identity",
]
