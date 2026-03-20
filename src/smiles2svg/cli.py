"""Command-line interface for smiles2svg."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from cairosvg import svg2png

from smiles2svg.render import render_smiles_to_svg


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate 2D molecular structure diagrams (SVG, optional PNG) from SMILES.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--smiles",
        type=str,
        default=None,
        help="SMILES string for one structure",
    )
    parser.add_argument(
        "-f",
        "--smiles-file",
        type=str,
        default=None,
        help="Text file with one SMILES per line",
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        default="SMILES",
        help="Base name for output files (single SMILES only; ignored for batch file input)",
    )
    parser.add_argument(
        "--style",
        type=str,
        choices=["plain", "names_hetero", "names_all", "stroke"],
        default="plain",
        help="Atom labeling style",
    )
    parser.add_argument(
        "--color",
        type=str,
        default="default",
        help="Fill color for all atoms, or 'default' for per-element colors",
    )
    parser.add_argument(
        "--png",
        action="store_true",
        help="Also write a PNG (requires Cairo)",
    )
    parser.add_argument(
        "--png-width",
        type=int,
        default=900,
        metavar="PX",
        help="PNG width in pixels",
    )
    parser.add_argument(
        "--font",
        type=str,
        default="Calibri, Helvetica, Arial, sans-serif",
        help="Font stack for atom labels (names_* styles); comma-separated like CSS font-family",
    )
    parser.add_argument(
        "--bond-color",
        type=str,
        default="#717171",
        help='Bond stroke color, e.g. "#RRGGBB"',
    )
    parser.add_argument(
        "--add-hydrogens",
        action="store_true",
        help="Draw implicit H on typical H-bond donors (O, N, P)",
    )
    parser.add_argument(
        "--canvas-size",
        type=int,
        default=800,
        metavar="PX",
        help="Target canvas size (px; longer molecule axis fits in size − 2×padding)",
    )
    parser.add_argument(
        "--padding",
        type=float,
        default=20.0,
        metavar="PX",
        help="Padding inside the SVG canvas (px)",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.smiles is None and args.smiles_file is None:
        parser.error("Provide either --smiles or --smiles-file")
    if args.smiles is not None and args.smiles_file is not None:
        parser.error("Use only one of --smiles or --smiles-file")

    if args.smiles is not None:
        smiles_list = [args.smiles]
    else:
        path = Path(args.smiles_file)
        if not path.is_file():
            print(f"Error: file not found: {path}", file=sys.stderr)
            return 1
        with path.open(encoding="utf-8") as f:
            smiles_list = [line.strip() for line in f if line.strip()]

    use_custom_name = args.name != "SMILES" and len(smiles_list) == 1

    for smile in smiles_list:
        if use_custom_name:
            svg_name = f"{args.name}.svg"
            png_name = f"{args.name}.png"
        else:
            svg_name = f"{smile}.svg"
            png_name = f"{smile}.png"

        svg_document = render_smiles_to_svg(
            smile,
            style=args.style,
            color=args.color,
            bond_color=args.bond_color,
            font=args.font,
            add_hydrogens=args.add_hydrogens,
            canvas_size=args.canvas_size,
            padding=args.padding,
            output_path=svg_name,
        )

        if args.png:
            svg2png(bytestring=svg_document.encode("utf-8"), write_to=png_name, output_width=args.png_width)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
