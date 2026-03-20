"""High-level SMILES → 2D layout → SVG pipeline."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from smiles2svg.drawing import draw_molecule
from smiles2svg.molecule import extract_bonds, extract_coordinates


def render_smiles_to_svg(
    smiles: str,
    *,
    style: str = "plain",
    color: str = "default",
    bond_color: str = "#717171",
    font: str = "Calibri, Helvetica, Arial, sans-serif",
    add_hydrogens: bool = False,
    canvas_size: int = 800,
    padding: float = 20.0,
    output_path: str | os.PathLike[str] | None = None,
) -> str:
    """
    Build a 2D structure from ``smiles`` and return the SVG document as a string.

    If ``output_path`` is given, the SVG is also written to that path.
    If omitted, a temporary file is used internally and removed after reading.

    Coordinates are rendered in **pixel space** (xyzrender-style): the molecule is
    fitted to ``canvas_size`` with ``padding``. With ``scale`` = px/Å from that fit,
    bond stroke and label ``font-size`` are fixed fractions of Å (× scale), so layout
    stays proportional when the canvas or molecule size changes.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    try:
        Chem.SanitizeMol(mol)
    except Chem.MolSanitizeException as e:
        raise ValueError(f"Invalid or unsanitizable SMILES: {smiles!r}") from e
    AllChem.Compute2DCoords(mol)
    atom_data = extract_coordinates(mol)
    bond_data = extract_bonds(mol)

    out = Path(output_path) if output_path is not None else None
    if out is None:
        fd, tmp = tempfile.mkstemp(suffix=".svg")
        os.close(fd)
        path = Path(tmp)
        try:
            return draw_molecule(
                atom_data,
                bond_data,
                style,
                color,
                font,
                bond_color,
                add_hydrogens,
                output_path=path,
                canvas_size=canvas_size,
                padding=padding,
            )
        finally:
            path.unlink(missing_ok=True)

    return draw_molecule(
        atom_data,
        bond_data,
        style,
        color,
        font,
        bond_color,
        add_hydrogens,
        output_path=out,
        canvas_size=canvas_size,
        padding=padding,
    )
