#!/usr/bin/env python3
"""Regenerate SVG figures under docs/_static/images/ (run from repo root: python docs/scripts/generate_doc_figures.py)."""

from __future__ import annotations

import os
import sys

# Repo root: .../smiles2svg
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

from smiles2svg import render_smiles_to_svg  # noqa: E402

OUT = os.path.join(ROOT, "docs", "_static", "images")

FIGURES: list[tuple[str, str, str, dict]] = [
    # (filename stem, SMILES, style, extra kwargs)
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


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    for stem, smiles, style, extra in FIGURES:
        path = os.path.join(OUT, f"{stem}.svg")
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
