"""Smoke tests for smiles2svg."""

import pytest


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
