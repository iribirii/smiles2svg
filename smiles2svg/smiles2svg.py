#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import AllChem
from cairosvg import svg2png

from smiles2svg.arg_parser import arg_parser
from smiles2svg.drawing_tools import draw_molecule 
from smiles2svg.mol_tools import extract_bonds, extract_coordinates


def main():
    args = arg_parser()

    if args.smiles != 'none':
        smiles = [args.smiles]
    elif args.smiles_file != 'none':
        with open(args.smiles_file, 'r') as f:
            smiles = [ line.rstrip('\n') for line in f.readlines()]

    for smile in smiles:
        if args.name == 'SMILES':
            svg_name = f'{smile}.svg' 
            png_name = f'{smile}.png' 
        else:
            svg_name = f'{args.name}.svg'
            png_name = f'{args.name}.png'

        # Load a mol from a SMILES string
        mol = Chem.MolFromSmiles(smile)
        draw_style = args.style
        draw_color = args.color
        bond_color = args.bond_color
        font = args.font

        Chem.Kekulize(mol)

        # Generate 2D coordinates using the ETKDG method
        AllChem.Compute2DCoords(mol)

        # Extract atom types, coordinates, and bond types
        atom_data = extract_coordinates(mol)
        bond_data = extract_bonds(mol)

        hydrogens = args.add_hydrogens

        svg = draw_molecule(svg_name, atom_data, bond_data, draw_style, draw_color, font, bond_color, hydrogens)

        if args.png:
            width = args.png_width

            svg2png(bytestring=svg, write_to=png_name, output_width=width)




