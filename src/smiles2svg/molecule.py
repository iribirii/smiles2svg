"""Extract 2D coordinates and bond data from RDKit molecules."""

from __future__ import annotations

from rdkit import Chem


def extract_coordinates(mol: Chem.Mol) -> list[tuple]:
    """Per-atom tuples: symbol, idx, atomic number, bond order sum, degree, (x, y)."""
    atom_data: list[tuple] = []

    for atom in mol.GetAtoms():
        atom_type = atom.GetSymbol()
        atom_number = atom.GetAtomicNum()
        atom_id = atom.GetIdx()
        atom_n_neighbors = atom.GetDegree()
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        atom_coords = (pos.x, pos.y)
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

        atom_data.append(
            (atom_type, atom_id, atom_number, n_bonds, atom_n_neighbors, atom_coords)
        )

    return atom_data


def extract_bonds(mol: Chem.Mol) -> list[tuple]:
    """Per-bond tuples: order, ((x1,y1),(x2,y2)), Z1, Z2, begin_idx, end_idx."""
    bond_data: list[tuple] = []

    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        begin_id = begin_atom.GetIdx()
        end_id = end_atom.GetIdx()
        bond_type = bond.GetBondTypeAsDouble()
        begin_pos = mol.GetConformer().GetAtomPosition(begin_atom.GetIdx())
        end_pos = mol.GetConformer().GetAtomPosition(end_atom.GetIdx())
        bond_coords = ((begin_pos.x, begin_pos.y), (end_pos.x, end_pos.y))
        atom1 = begin_atom.GetAtomicNum()
        atom2 = end_atom.GetAtomicNum()
        bond_data.append((bond_type, bond_coords, atom1, atom2, begin_id, end_id))

    return bond_data
