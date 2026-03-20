"""Extract 2D coordinates and bond data from RDKit molecules."""

from __future__ import annotations

from typing import NamedTuple

from rdkit import Chem


class AtomRecord(NamedTuple):
    """One atom after 2D embedding (coordinates in Å).

    ``bond_order_sum`` is the sum of incident bond orders (single=1, aromatic=1.5,
    double=2, triple=3). ``degree`` is the number of bonded neighbors.
    """

    symbol: str
    atom_idx: int
    atomic_num: int
    bond_order_sum: float
    degree: int
    coord: tuple[float, float]


class BondRecord(NamedTuple):
    """One bond with 2D endpoints."""

    order: float
    start_xy: tuple[float, float]
    end_xy: tuple[float, float]
    z_begin: int
    z_end: int
    begin_atom_idx: int
    end_atom_idx: int


def extract_coordinates(mol: Chem.Mol) -> list[AtomRecord]:
    """Build per-atom records from the molecule's first conformer."""
    atom_data: list[AtomRecord] = []

    for atom in mol.GetAtoms():
        atom_type = atom.GetSymbol()
        atom_number = atom.GetAtomicNum()
        atom_id = atom.GetIdx()
        degree = atom.GetDegree()
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        coord = (pos.x, pos.y)
        n_bonds = 0.0
        for neighbor in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                n_bonds += 1
            elif bond.GetBondType() == Chem.rdchem.BondType.AROMATIC:
                n_bonds += 1.5
            elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                n_bonds += 2
            elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
                n_bonds += 3

        atom_data.append(AtomRecord(atom_type, atom_id, atom_number, n_bonds, degree, coord))

    return atom_data


def extract_bonds(mol: Chem.Mol) -> list[BondRecord]:
    """Build per-bond records from the molecule's first conformer."""
    bond_data: list[BondRecord] = []

    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        begin_id = begin_atom.GetIdx()
        end_id = end_atom.GetIdx()
        order = bond.GetBondTypeAsDouble()
        begin_pos = mol.GetConformer().GetAtomPosition(begin_atom.GetIdx())
        end_pos = mol.GetConformer().GetAtomPosition(end_atom.GetIdx())
        start_xy = (begin_pos.x, begin_pos.y)
        end_xy = (end_pos.x, end_pos.y)
        z1 = begin_atom.GetAtomicNum()
        z2 = end_atom.GetAtomicNum()
        bond_data.append(BondRecord(order, start_xy, end_xy, z1, z2, begin_id, end_id))

    return bond_data
