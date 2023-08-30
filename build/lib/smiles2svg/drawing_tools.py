#!/usr/bin/env python3

import svgwrite
import numpy as np
import math

from smiles2svg.utils import get_atom_radii, get_atom_colors

def draw_molecule(filename, atoms, bonds, draw_style, draw_color, font, bond_color, hydrogens):

    # Get atomic radii and colors
    atom_radii = get_atom_radii()
    atom_colors = get_atom_colors(draw_color)
    dwg = svgwrite.Drawing(filename)

    if draw_color != 'default':
        bond_color = atom_colors[0]

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

    if hydrogens:
        add_hydroges(dwg, atoms, bonds, atom_colors, atom_radii, draw_style, bond_color)

    # Set viewBox size
    coords = [ a[5] for a in atoms ]
    r_list = [ atom_radii[a[2]] for a in atoms ]
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
        a_name, a_id, a_number, n_neighbors, n_bonds, coord = atom
        cx, cy = coord
        r = atom_radii[a_number]/1.5
        color = atom_colors[a_number]

        circle = dwg.circle(center=(cx,cy), r=r, fill=color, stroke='none')
        dwg.add(circle)

def draw_atoms_stroke(dwg, atoms, atom_radii, atom_colors, bond_color):
    for atom in atoms:
        a_name, a_id, a_number, n_neighbors, n_bonds, coord = atom
        cx, cy = coord
        r = atom_radii[a_number]/1.5
        color = atom_colors[a_number]

        circle = dwg.circle(center=(cx,cy), r=r, fill=color, stroke=bond_color, stroke_width=0.1)
        dwg.add(circle)

def draw_atoms_hetero_names(dwg, atoms, atom_radii, atom_colors, font):
    for atom in atoms:
        a_name, a_id, a_number, n_neighbors, n_bonds, coord = atom
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
        a_name, a_id, a_number, n_neighbors, n_bonds, coord = atom
        cx, cy = coord
        r = atom_radii[a_number]/1.5
        color = atom_colors[a_number]

        r += 0.1
        circle = dwg.circle(center=(cx,cy), r=r, fill='#ffffff', stroke=color, stroke_width=0.1)
        text = dwg.text(a_name, insert=(cx-0.3,cy+0.25), font_size=0.8, font_family=font, fill=color)
        
        dwg.add(circle)
        dwg.add(text)

def add_hydroges(dwg, atoms, bonds, atom_colors, atom_radii, style, bond_color):
    if style == 'plain':
        r_extra = 0
        s_width = 0
    elif (style == 'names_hetero') or (style == 'names_all'):
        r_extra = 0.1
        s_width = 0
    elif (style == 'stroke'): 
        r_extra = 0
        s_width = 0.1

    h_r = atom_radii[1]/1.5
    h_color = atom_colors[1]

    for atom in atoms:
        a_name, a_id, a_number, n_bonds, n_neighbors, coord = atom

        if (a_name in ['O','S']) and (n_bonds == 1) and (n_neighbors == 1): 
            a_bond = [ x for x in bonds if (x[-2] == a_id) or (x[-1] == a_id)][0]
            a_r = atom_radii[a_number]/1.5

            if a_id == a_bond[-2]:
                start, end = a_bond[1]
                vec = np.array(end) - np.array(start)
            elif a_id == a_bond[-1]:
                end, start = a_bond[1]
                vec = np.array(end) - np.array(start)
            sign_x = np.sign(vec[0])
            sign_y = np.sign(vec[1])
            sign = sign_x * sign_y
            sign = 1 if abs(sign) == 0 else sign
            angle = sign * np.deg2rad(120)
            rot = np.array([[math.cos(angle), -math.sin(angle)],[math.sin(angle), math.cos(angle)]])
            rotated = np.dot(rot, vec)
            d = a_r + h_r + r_extra
            length = np.linalg.norm(rotated)
            final_vec = ( rotated / length ) * d
            h_center = start + final_vec
            circle = dwg.circle(center=(h_center), r=h_r, fill=h_color, stroke=bond_color, stroke_width=s_width)
            dwg.add(circle)

        if (a_name in ['N','P']) and (n_bonds == 2) and (n_neighbors == 1): 
            a_bond = [ x for x in bonds if (x[-2] == a_id) or (x[-1] == a_id)][0]
            a_r = atom_radii[a_number]/1.5
            h_color = atom_colors[1]

            if a_id == a_bond[-2]:
                start, end = a_bond[1]
                vec = np.array(end) - np.array(start)
            elif a_id == a_bond[-1]:
                end, start = a_bond[1]
                vec = np.array(end) - np.array(start)
            sign_x = np.sign(vec[0])
            sign_y = np.sign(vec[1])
            sign = sign_x * sign_y
            sign = 1 if abs(sign) == 0 else sign
            angle = sign * np.deg2rad(120)
            rot = np.array([[math.cos(angle), -math.sin(angle)],[math.sin(angle), math.cos(angle)]])
            rotated = np.dot(rot, vec)
            d = a_r + h_r + r_extra
            length = np.linalg.norm(rotated)
            final_vec = ( rotated / length ) * d
            h_center = start + final_vec
            circle = dwg.circle(center=(h_center), r=h_r, fill=h_color, stroke=bond_color, stroke_width=s_width)
            dwg.add(circle)

        if (a_name in ['N','P']) and (n_bonds == 1) and (n_neighbors == 1): 
            a_bond = [ x for x in bonds if (x[-2] == a_id) or (x[-1] == a_id)][0]
            a_r = atom_radii[a_number]/1.5
            h_color = atom_colors[1]
            d = a_r + h_r + r_extra

            if a_id == a_bond[-2]:
                start, end = a_bond[1]
                vec = np.array(end) - np.array(start)
            elif a_id == a_bond[-1]:
                end, start = a_bond[1]
                vec = np.array(end) - np.array(start)

            sign_x = np.sign(vec[0])
            sign_y = np.sign(vec[1])
            sign = sign_x * sign_y
            sign = 1 if abs(sign) == 0 else sign

            angle1 = sign * np.deg2rad(120)
            rot1 = np.array([[math.cos(angle1), -math.sin(angle1)],[math.sin(angle1), math.cos(angle1)]])
            rotated1 = np.dot(rot1, vec)
            length = np.linalg.norm(rotated1)
            final_vec1 = ( rotated1 / length ) * d
            h_center1 = start + final_vec1
            circle1 = dwg.circle(center=(h_center1), r=h_r, fill=h_color, stroke=bond_color, stroke_width=s_width)
            dwg.add(circle1)
            
            angle2 = -sign * np.deg2rad(120)
            rot2 = np.array([[math.cos(angle2), -math.sin(angle2)],[math.sin(angle2), math.cos(angle2)]])
            rotated2 = np.dot(rot2, vec)
            length = np.linalg.norm(rotated2)
            final_vec2 = ( rotated2 / length ) * d
            h_center2 = start + final_vec2
            circle2 = dwg.circle(center=(h_center2), r=h_r, fill=h_color, stroke=bond_color, stroke_width=s_width)
            dwg.add(circle2)

        if (a_name in ['N','P']) and (n_bonds == 2) and (n_neighbors == 2): 
            
            a_bond = [ x for x in bonds if (x[-2] == a_id) or (x[-1] == a_id)]
            b1, b2 = a_bond

            if a_id == b1[-2]:
                start, end = b1[1]
                vec1 = np.array(end) - np.array(start)
            elif a_id == b1[-1]:
                end, start = b1[1]
                vec1 = np.array(end) - np.array(start)

            if a_id == b2[-2]:
                start, end = b2[1]
                vec2 = np.array(end) - np.array(start)
            elif a_id == b2[-1]:
                end, start = b2[1]
                vec2 = np.array(end) - np.array(start)

            vec_bis = -1 * ( vec1 + vec2 ) / 2

            a_r = atom_radii[a_number]/1.5
            h_color = atom_colors[1]

            d = a_r + h_r + r_extra
            length = np.linalg.norm(vec_bis)

            final_vec = ( vec_bis / length ) * d
            h_center = start + final_vec
            circle = dwg.circle(center=(h_center), r=h_r, fill=h_color, stroke=bond_color, stroke_width=s_width)
            dwg.add(circle)


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

