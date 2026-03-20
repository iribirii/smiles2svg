#!/usr/bin/env python3
"""Regenerate SVG figures under docs/_static/images/.

Run from the repository root (with the package importable), for example::

    pip install -e .
    python docs/scripts/generate_doc_figures.py

Or::

    PYTHONPATH=src python docs/scripts/generate_doc_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    src = root / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))

    from smiles2svg import render_smiles_to_svg

    out_dir = root / "docs" / "_static" / "images"
    out_dir.mkdir(parents=True, exist_ok=True)

    figures: list[tuple[str, str, str, dict]] = [
        ("caffeine_plain", "Cn1cnc2n(C)c(=O)n(C)c(=O)c12", "plain", {}),
        ("caffeine_hetero", "Cn1cnc2n(C)c(=O)n(C)c(=O)c12", "names_hetero", {}),
        ("caffeine_all", "Cn1cnc2n(C)c(=O)n(C)c(=O)c12", "names_all", {}),
        ("caffeine_stroke", "Cn1cnc2n(C)c(=O)n(C)c(=O)c12", "stroke", {}),
        ("caffeine_h", "Cn1cnc2n(C)c(=O)n(C)c(=O)c12", "plain", {"add_hydrogens": True}),
        ("benzene", "c1ccccc1", "plain", {}),
        ("benzene_labels", "c1ccccc1", "names_all", {}),
        ("ethanol", "CCO", "names_hetero", {}),
        ("water", "O", "names_hetero", {}),
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "names_hetero", {}),
        ("toluene", "Cc1ccccc1", "plain", {}),
        ("naphthalene", "c1ccc2ccccc2c1", "plain", {}),
        ("glycine", "NCC(=O)O", "names_hetero", {}),
    ]

    for stem, smiles, style, extra in figures:
        path = out_dir / f"{stem}.svg"
        render_smiles_to_svg(
            smiles,
            style=style,
            output_path=path,
            canvas_size=640,
            padding=16,
            **extra,
        )
        print("wrote", path)


if __name__ == "__main__":
    main()
