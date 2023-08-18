#!/usr/bin/env python3

import argparse
from ctypes import resize
from rdkit import Chem
from rdkit.Chem import AllChem
import svgwrite
import numpy as np
from cairosvg import svg2png
import math


def extract_coordinates(mol):
    atom_data = []

    for atom in mol.GetAtoms():
        atom_type = atom.GetSymbol()
        atom_number = atom.GetAtomicNum()
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        atom_coords = (pos.x, pos.y)
        # Count the number of bonds based on their types
        n_bonds = 0
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

        atom_data.append((atom_type, atom_number, n_bonds, atom_coords))

    return atom_data

def extract_bonds(mol):
    bond_data = []

    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        bond_type = bond.GetBondTypeAsDouble()  # Change to GetBondType() for string representation
        begin_pos = mol.GetConformer().GetAtomPosition(begin_atom.GetIdx())
        end_pos = mol.GetConformer().GetAtomPosition(end_atom.GetIdx())
        bond_coords = ((begin_pos.x, begin_pos.y), (end_pos.x, end_pos.y))
        atom1 = begin_atom.GetAtomicNum()
        atom2 = end_atom.GetAtomicNum()
        bond_data.append((bond_type, bond_coords, atom1, atom2))

    return bond_data

def atom_parameters(color):
    radii = [ 0.0,
      0.35, 0.20,
      1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 0.40,
      1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 0.90,
      2.20, 1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35,
      1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.00,
      2.35, 2.00, 1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40,
      1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 1.30,
      2.35, 1.98, 1.69, 1.65, 1.65, 1.64, 1.65, 1.66, 1.85,
                  1.61, 1.59, 1.59, 1.58, 1.57, 1.56, 1.70,
                  1.56, 1.44, 1.34, 1.30, 1.28, 1.26, 1.26, 1.29,
      1.34, 1.44, 1.55, 1.54, 1.52, 1.53, 1.53, 1.50,
      2.70, 2.23, 1.87, 1.78, 1.61, 1.40, 1.40, 1.40, 1.40,
                  1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40,
                  1.40, 1.40
           ] + [1.40]*23

    if color == 'default':
        colors = [ '#7f7f7f',
          '#d0d0d0',                                                                                                                                                                                 '#A0A0A0',
          '#A0A0A0', '#A0A0A0',                                                                                                               '#FF99FF', '#404040', '#2020FF', '#FF2020', '#00BB00', '#A0A0A0',
          '#880000', '#A0A0A0',                                                                                                               '#A0A0A0', '#090909', '#FF8800', '#F0F000', '#55FF55', '#A0A0A0',
          '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#AA3311', '#A0A0A0',
          '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#A0A0A0', '#802090', '#A0A0A0',
                ] + ['#A0A0A0']*64
    else:
        colors = [color]*128
    return radii, colors

def draw_molecule(filename, atoms, bonds, draw_style, draw_color, font, bond_color):

    # Get atomic radii and colors
    atom_radii, atom_colors = atom_parameters(draw_color)
    dwg = svgwrite.Drawing(filename)

    draw_bonds(dwg, bonds, bond_color)

    # Draw atoms
    if draw_style == 'plain':
        draw_atoms_plain(dwg, atoms, atom_radii, atom_colors)
    elif draw_style == 'names_hetero':
        draw_atoms_hetero_names(dwg, atoms, atom_radii, atom_colors, font)
    elif draw_style == 'names_all':
        draw_atoms_all_names(dwg, atoms, atom_radii, atom_colors, font)
    elif draw_style == 'stroke':
        draw_atoms_stroke(dwg, atoms, atom_radii, atom_colors, bond_color)

    # Set viewBox size
    coords = [ a[3] for a in atoms ]
    r_list = [ atom_radii[a[1]] for a in atoms ]
    set_viewbox(dwg, coords, r_list ) 

    dwg.save()

    return dwg.tostring()

def draw_bonds(dwg, bonds, color):
    for bond in bonds:
        bond_order = bond[0]
        bond_start, bond_end = bond[1]

        if bond_order == 1:
            line = dwg.line(start=bond_start, end=bond_end, stroke=color, stroke_width=0.1)
            dwg.add(line)
        elif bond_order == 2 or bond_order == 1.5:
            dist=0.2
            start_up, end_up = displace_perpendicular(bond_start, bond_end, dist)
            line_up = dwg.line(start=start_up, end=end_up, stroke=color, stroke_width=0.1)
            dwg.add(line_up)
            start_down, end_down = displace_perpendicular(bond_start, bond_end, -dist)
            line_down = dwg.line(start=start_down, end=end_down, stroke=color, stroke_width=0.1)
            dwg.add(line_down)
        elif bond_order == 3:
            dist=0.25
            line_mid = dwg.line(start=bond_start, end=bond_end, stroke=color, stroke_width=0.1)
            dwg.add(line_mid)
            start_up, end_up = displace_perpendicular(bond_start, bond_end, dist)
            line_up = dwg.line(start=start_up, end=end_up, stroke=color, stroke_width=0.1)
            dwg.add(line_up)
            start_down, end_down = displace_perpendicular(bond_start, bond_end, -dist)
            line_down = dwg.line(start=start_down, end=end_down, stroke=color, stroke_width=0.1)
            dwg.add(line_down)


def displace_perpendicular(start, end, dist):
    original_vector = np.array(end) - np.array(start)

    perpendicular_vector = np.array([-original_vector[1], original_vector[0]])
    perpendicular_vector = perpendicular_vector / np.linalg.norm(perpendicular_vector)

    new_start = np.array(start) + perpendicular_vector * dist
    new_end = np.array(end) + perpendicular_vector * dist

    return new_start.tolist(), new_end.tolist()

def calculate_angle(vec):
    angle_radians = math.atan2(vec[1], vec[0])
    angle_degrees = math.degrees(angle_radians)

    if angle_degrees < 0:
        angle_degrees = 360 + angle_degrees
    
    return angle_degrees

def draw_atoms_plain(dwg, atoms, atom_radii, atom_colors):
    for atom in atoms:
        a_name, a_number, n_bonds, coord = atom
        cx, cy = coord
        r = atom_radii[a_number]/1.5
        color = atom_colors[a_number]

        circle = dwg.circle(center=(cx,cy), r=r, fill=color, stroke='none')
        dwg.add(circle)

def draw_atoms_stroke(dwg, atoms, atom_radii, atom_colors, bond_color):
    for atom in atoms:
        a_name, a_number, n_bonds, coord = atom
        cx, cy = coord
        r = atom_radii[a_number]/1.5
        color = atom_colors[a_number]

        circle = dwg.circle(center=(cx,cy), r=r, fill=color, stroke=bond_color, stroke_width=0.1)
        dwg.add(circle)

def draw_atoms_hetero_names(dwg, atoms, atom_radii, atom_colors, font):
    for atom in atoms:
        a_name, a_number, n_bonds, coord = atom
        cx, cy = coord
        r = atom_radii[a_number]/1.5
        color = atom_colors[a_number]

        # if a_name in ['N','O','S','P']:
        if a_name != 'C': 
            r += 0.1
            circle = dwg.circle(center=(cx,cy), r=r, fill='#ffffff', stroke=color, stroke_width=0.1)
            text = dwg.text(a_name, insert=(cx-0.3,cy+0.25), font_size=0.8, font_family=font, fill=color)
        
            dwg.add(circle)
            dwg.add(text)
        else:
            circle = dwg.circle(center=(cx,cy), r=r, fill=color, stroke='none')
            dwg.add(circle)

def draw_atoms_all_names(dwg, atoms, atom_radii, atom_colors,font):
    for atom in atoms:
        a_name, a_number, n_bonds, coord = atom
        cx, cy = coord
        r = atom_radii[a_number]/1.5
        color = atom_colors[a_number]

        r += 0.1
        circle = dwg.circle(center=(cx,cy), r=r, fill='#ffffff', stroke=color, stroke_width=0.1)
        text = dwg.text(a_name, insert=(cx-0.3,cy+0.25), font_size=0.8, font_family=font, fill=color)
        
        dwg.add(circle)
        dwg.add(text)

def set_viewbox(dwg, coords, rs):
    r = max(rs)
    #Calculate the bounding box for all circle coordinates
    min_x = min(cx - 2*r for cx, cy in coords)
    min_y = min(cy - 2*r for cx, cy in coords)
    max_x = max(cx + 2*r for cx, cy in coords)
    max_y = max(cy + 2*r for cx, cy in coords)

    #Calculate viewBox and preserveAspectRatio values
    viewBox = (min_x, min_y, max_x - min_x, max_y - min_y)

    # Set the viewBox and preserveAspectRatio attributes
    dwg.viewbox(*viewBox)
    dwg.fit(scale='meet')

def main(args):
    if args.smiles != 'none':
        smiles = [args.smiles]
    elif args.smiles_file != 'none':
        with open(args.smiles_file, 'r') as f:
            smiles = [ line.rstrip('\n') for line in f.readlines()]

    for smile in smiles:
        svg_name = f'{smile}.svg' 

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

        svg = draw_molecule(svg_name, atom_data, bond_data, draw_style, draw_color, font, bond_color)

        if args.png:
            width = args.png_width

            svg2png(bytestring=svg, write_to=f'{smile}.png', output_width=900)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Generate 2D coordinates, draw, and save mol",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument(
            "-s","--smiles", 
            type=str, 
            default='none',
            help="SMILES string of the mol")
    parser.add_argument(
            "-f", "--smiles_file",
            type=str,
            default='none',
            help="Name of the file with all the SMILES codes")
    # parser.add_argument(
    #         "-n", "--name",
    #         type=str,
    #         default='none',
    #         help="Name for the SVG image")
    parser.add_argument(
            "--style", 
            type=str, 
            choices=['plain','names_hetero','names_all','stroke'], 
            default='plain', 
            help="Select the style for the atoms.")
    parser.add_argument(
            "--color", 
            type=str, 
            default='default', 
            help="Select a color for all the molecule. If 'default', atoms will have different colors depending on the element.")
    parser.add_argument(
            "--png",
            action='store_true',
            help="Saves the figure in png format as well.")
    parser.add_argument(
            "--png_width",
            type=int,
            default='900',
            help='Select the image width in pixels.')
    parser.add_argument(
            "--font",
            type=str,
            default='Calibri',
            help='Select the font for the atomic symbols.')
    parser.add_argument(
            "--bond_color",
            type=str,
            default='#717171',
            help='Select the color for the bonds format "#RRGGBB".')
    args = parser.parse_args()
    main(args)

