#!/usr/bin/env python3

import argparse

def arg_parser():
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
    parser.add_argument(
            "-n", "--name",
            type=str,
            default='SMILES',
            help="Name for the SVG (and PNG) image. Only available for 1 SMILES code, not for list of SMILES.")
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
    parser.add_argument(
            "--add_hydrogens",
            action='store_true',
            help="Adds Hydrogen atoms to the HB-Donors.")
    args = parser.parse_args()

    return args

