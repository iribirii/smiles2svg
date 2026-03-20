"""Render 2D molecular structures from SMILES strings as SVG."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from smiles2svg.molecule import AtomRecord, BondRecord
from smiles2svg.render import render_smiles_to_svg

try:
    __version__ = version("smiles2svg")
except PackageNotFoundError:
    __version__ = "0.2.0"

__all__ = ["__version__", "AtomRecord", "BondRecord", "render_smiles_to_svg"]
