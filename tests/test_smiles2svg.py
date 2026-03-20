"""Smoke tests for smiles2svg."""

import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def test_render_smiles_to_svg_returns_svg():
    from smiles2svg import render_smiles_to_svg

    out = render_smiles_to_svg("CCO", style="plain")
    assert "<svg" in out.lower()


def test_invalid_smiles():
    from smiles2svg import render_smiles_to_svg

    with pytest.raises(ValueError, match="Invalid SMILES"):
        render_smiles_to_svg("not-a-valid-smiles-string!!!")


def test_cli_parser_smiles():
    from smiles2svg.cli import build_parser

    args = build_parser().parse_args(["-s", "CCO", "--style", "stroke"])
    assert args.smiles == "CCO"
    assert args.style == "stroke"
