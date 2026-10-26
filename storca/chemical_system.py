"""Shared molecular identity for multi-lane stability campaigns.

Amounts and thermodynamic conditions deliberately do not live here: they are
properties of a lane environment.  This contract only establishes which
chemical components are shared by the independently interpreted lanes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import re
from typing import Any, Mapping


CHEMICAL_SYSTEM_SCHEMA_VERSION = 1
COMPONENT_ROLES = {
    "target", "product", "environment", "reactant", "catalyst", "impurity", "additive",
}
_IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {
            str(key): _jsonable(item)
            for key, item in sorted(value.items(), key=lambda entry: str(entry[0]))
        }
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    return value


def _contract_id(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        _jsonable(payload), sort_keys=True, separators=(",", ":"), allow_nan=False,
    )
    return "sha256:" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _required_text(value: Any, *, name: str) -> str:
    text = str(value or "").strip()
    if not text:
        raise ValueError(f"{name} is required")
    return text


def _component_identifier(value: Any, *, name: str) -> str:
    identifier = _required_text(value, name=name)
    if not _IDENTIFIER.fullmatch(identifier):
        raise ValueError(f"{name} must contain only letters, numbers, '.', '_', or '-'")
    return identifier


def _canonical_molecular_identity(smiles: Any) -> tuple[str, int]:
    value = _required_text(smiles, name="smiles")
    try:
        from rdkit import Chem
    except ImportError as error:
        raise RuntimeError("RDKit is required to validate campaign chemical identity") from error
    molecule = Chem.MolFromSmiles(value)
    if molecule is None:
        raise ValueError(f"Invalid component SMILES: {value}")
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return (
        Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True),
        int(Chem.GetFormalCharge(molecule)),
    )


@dataclass(frozen=True)
class ChemicalComponent:
    """One state-specific molecular identity shared by campaign lanes."""

    component_id: str
    smiles: str
    charge: int
    multiplicity: int
    role: str
    identity_provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        component_id = _component_identifier(self.component_id, name="component_id")
        canonical, derived_charge = _canonical_molecular_identity(self.smiles)
        charge = int(self.charge)
        multiplicity = int(self.multiplicity)
        role = _required_text(self.role, name=f"{component_id}.role").lower()
        if charge != derived_charge:
            raise ValueError(
                f"Component {component_id} charge does not match its SMILES: "
                f"declared={charge}, derived={derived_charge}"
            )
        if multiplicity < 1:
            raise ValueError(f"Component {component_id} multiplicity must be positive")
        if role not in COMPONENT_ROLES:
            raise ValueError(
                f"Component {component_id} role must be one of {sorted(COMPONENT_ROLES)}"
            )
        object.__setattr__(self, "component_id", component_id)
        object.__setattr__(self, "smiles", canonical)
        object.__setattr__(self, "charge", charge)
        object.__setattr__(self, "multiplicity", multiplicity)
        object.__setattr__(self, "role", role)
        object.__setattr__(
            self, "identity_provenance", _jsonable(dict(self.identity_provenance)),
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ChemicalComponent":
        if not isinstance(value, Mapping):
            raise TypeError("chemical component must be a mapping")
        component = cls(
            component_id=value.get("component_id"),
            smiles=value.get("smiles") or value.get("canonical_smiles"),
            charge=value.get("charge", 0),
            multiplicity=value.get("multiplicity", 1),
            role=value.get("role"),
            identity_provenance=value.get("identity_provenance") or {},
        )
        retained_canonical = value.get("canonical_smiles")
        if retained_canonical is not None:
            canonical, _ = _canonical_molecular_identity(retained_canonical)
            if canonical != component.smiles:
                raise ValueError(
                    f"Component {component.component_id} retained canonical_smiles does not "
                    "match its SMILES"
                )
        retained_contract = value.get("contract_id")
        if retained_contract is not None and retained_contract != component.as_dict()["contract_id"]:
            raise ValueError(f"Component {component.component_id} contract_id does not match its content")
        return component

    def _payload(self) -> dict[str, Any]:
        return {
            "schema_version": CHEMICAL_SYSTEM_SCHEMA_VERSION,
            "kind": "storca_chemical_component",
            "component_id": self.component_id,
            "smiles": self.smiles,
            "canonical_smiles": self.smiles,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "role": self.role,
            "identity_provenance": _jsonable(self.identity_provenance),
        }

    def as_dict(self) -> dict[str, Any]:
        payload = self._payload()
        return {**payload, "contract_id": _contract_id(payload)}


@dataclass(frozen=True)
class ChemicalSystem:
    """Immutable component set common to all requested environments."""

    target_component_id: str
    components: tuple[ChemicalComponent, ...]
    material_context: Mapping[str, Any] = field(default_factory=dict)
    user_declared_applicability_flags: Mapping[str, bool] = field(default_factory=dict)

    def __post_init__(self) -> None:
        target = _component_identifier(
            self.target_component_id, name="target_component_id",
        )
        components = tuple(sorted(self.components, key=lambda item: item.component_id))
        if not components:
            raise ValueError("ChemicalSystem needs at least one component")
        identities = [item.component_id for item in components]
        if len(identities) != len(set(identities)):
            raise ValueError("ChemicalSystem component_id values must be unique")
        if target not in identities:
            raise ValueError("target_component_id does not identify a retained component")
        target_roles = [item.component_id for item in components if item.role == "target"]
        if target_roles != [target]:
            raise ValueError(
                "ChemicalSystem needs exactly one component with role='target', matching "
                "target_component_id"
            )
        state_identities = [
            (item.smiles, item.charge, item.multiplicity) for item in components
        ]
        if len(state_identities) != len(set(state_identities)):
            raise ValueError(
                "ChemicalSystem contains duplicate molecular state identities; use lane "
                "composition to represent additional amounts"
            )
        flags = {
            str(key): bool(value)
            for key, value in dict(self.user_declared_applicability_flags).items()
        }
        if any(not str(key).strip() for key in flags):
            raise ValueError("Applicability flag names must be nonempty")
        object.__setattr__(self, "target_component_id", target)
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "material_context", _jsonable(dict(self.material_context)))
        object.__setattr__(
            self, "user_declared_applicability_flags", dict(sorted(flags.items())),
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ChemicalSystem":
        if not isinstance(value, Mapping):
            raise TypeError("chemical_system must be a mapping")
        system = cls(
            target_component_id=value.get("target_component_id"),
            components=tuple(
                ChemicalComponent.from_dict(item) for item in value.get("components") or []
            ),
            material_context=value.get("material_context") or {},
            user_declared_applicability_flags=(
                value.get("user_declared_applicability_flags") or {}
            ),
        )
        retained_contract = value.get("contract_id")
        if retained_contract is not None and retained_contract != system.as_dict()["contract_id"]:
            raise ValueError("ChemicalSystem contract_id does not match its content")
        return system

    @property
    def target(self) -> ChemicalComponent:
        return next(
            item for item in self.components if item.component_id == self.target_component_id
        )

    def _payload(self) -> dict[str, Any]:
        return {
            "schema_version": CHEMICAL_SYSTEM_SCHEMA_VERSION,
            "kind": "storca_chemical_system",
            "target_component_id": self.target_component_id,
            "components": [item.as_dict() for item in self.components],
            "material_context": _jsonable(self.material_context),
            "user_declared_applicability_flags": _jsonable(
                self.user_declared_applicability_flags,
            ),
        }

    def as_dict(self) -> dict[str, Any]:
        payload = self._payload()
        return {**payload, "contract_id": _contract_id(payload)}


__all__ = [
    "CHEMICAL_SYSTEM_SCHEMA_VERSION", "COMPONENT_ROLES", "ChemicalComponent",
    "ChemicalSystem",
]
