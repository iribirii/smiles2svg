"""SVG drawing of 2D molecular graphs (pixel canvas, xyzrender-style scaling)."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Union

import numpy as np
import svgwrite

from smiles2svg.elements import get_atom_colors, get_atom_radii
from smiles2svg.layout import (
    BOND_DOUBLE_OFFSET,
    BOND_TRIPLE_OFFSET,
    PixelViewport,
    compute_scene_bounds,
    displace_perpendicular_mol,
    make_viewport,
)

PathLike = Union[str, Path]

# With ``font-size`` in px, many SVG stacks (incl. CairoSVG) ignore or mishandle
# ``dominant-baseline="central"``. Use the alphabetic baseline and shift Y so the
# visual center of typical Latin capitals sits on (cx, cy): cap height ~0.7em.
_LABEL_BASELINE_SHIFT_EM = 0.35


def _atom_symbol_text(
    dwg: svgwrite.Drawing,
    symbol: str,
    cx: float,
    cy: float,
    color: str,
    font: str,
    font_size_px: float,
) -> svgwrite.text.Text:
    """Horizontally center symbol on (cx, cy); vertical center via baseline offset (px fonts)."""
    baseline_y = cy + _LABEL_BASELINE_SHIFT_EM * font_size_px
    return dwg.text(
        symbol,
        insert=(cx, baseline_y),
        font_size=f"{font_size_px:.2f}px",
        font_family=font,
        font_weight="bold",
        fill=color,
        text_anchor="middle",
    )


def draw_molecule(
    atoms: list,
    bonds: list,
    draw_style: str,
    draw_color: str,
    font: str,
    bond_color: str,
    hydrogens: bool,
    *,
    output_path: PathLike,
    canvas_size: int = 800,
    padding: float = 20.0,
) -> str:
    """Write SVG to ``output_path`` and return the SVG as a string."""
    path = Path(output_path)
    atom_radii = get_atom_radii()
    atom_colors = get_atom_colors(draw_color)

    lo, hi = compute_scene_bounds(atoms, bonds, draw_style, atom_radii, hydrogens)
    vp = make_viewport(lo, hi, canvas_size, padding)

    dwg = svgwrite.Drawing(
        str(path),
        size=(f"{vp.canvas_w}px", f"{vp.canvas_h}px"),
        profile="full",
    )
    dwg.viewbox(0, 0, vp.canvas_w, vp.canvas_h)

    if draw_color != "default":
        bond_color = atom_colors[0]

    draw_bonds(dwg, bonds, bond_color, vp)
    if draw_style == "plain":
        draw_atoms_plain(dwg, atoms, atom_radii, atom_colors, vp)
    elif draw_style == "names_hetero":
        draw_atoms_hetero_names(dwg, atoms, atom_radii, atom_colors, font, vp)
    elif draw_style == "names_all":
        draw_atoms_all_names(dwg, atoms, atom_radii, atom_colors, font, vp)
    elif draw_style == "stroke":
        draw_atoms_stroke(dwg, atoms, atom_radii, atom_colors, bond_color, vp)

    if hydrogens:
        add_hydrogens_to_donors(dwg, atoms, bonds, atom_colors, atom_radii, draw_style, bond_color, vp)

    dwg.save()
    return dwg.tostring()


def draw_bonds(dwg: svgwrite.Drawing, bonds: list, color: str, vp: PixelViewport) -> None:
    bw = vp.bond_width_px
    for bond in bonds:
        bond_order = bond[0]
        bond_start, bond_end = bond[1]
        p1 = vp.to_px(bond_start[0], bond_start[1])
        p2 = vp.to_px(bond_end[0], bond_end[1])

        if bond_order == 1:
            line = dwg.line(start=p1, end=p2, stroke=color, stroke_width=bw)
            dwg.add(line)
        elif bond_order == 2 or bond_order == 1.5:
            dist = BOND_DOUBLE_OFFSET
            start_up, end_up = displace_perpendicular_mol(bond_start, bond_end, dist)
            start_down, end_down = displace_perpendicular_mol(bond_start, bond_end, -dist)
            line_up = dwg.line(
                start=vp.to_px(start_up[0], start_up[1]),
                end=vp.to_px(end_up[0], end_up[1]),
                stroke=color,
                stroke_width=bw,
            )
            dwg.add(line_up)
            line_down = dwg.line(
                start=vp.to_px(start_down[0], start_down[1]),
                end=vp.to_px(end_down[0], end_down[1]),
                stroke=color,
                stroke_width=bw,
            )
            dwg.add(line_down)
        elif bond_order == 3:
            dist = BOND_TRIPLE_OFFSET
            line_mid = dwg.line(start=p1, end=p2, stroke=color, stroke_width=bw)
            dwg.add(line_mid)
            start_up, end_up = displace_perpendicular_mol(bond_start, bond_end, dist)
            start_down, end_down = displace_perpendicular_mol(bond_start, bond_end, -dist)
            line_up = dwg.line(
                start=vp.to_px(start_up[0], start_up[1]),
                end=vp.to_px(end_up[0], end_up[1]),
                stroke=color,
                stroke_width=bw,
            )
            dwg.add(line_up)
            line_down = dwg.line(
                start=vp.to_px(start_down[0], start_down[1]),
                end=vp.to_px(end_down[0], end_down[1]),
                stroke=color,
                stroke_width=bw,
            )
            dwg.add(line_down)


def draw_atoms_plain(dwg, atoms, atom_radii, atom_colors, vp: PixelViewport):
    for atom in atoms:
        a_name, a_id, a_number, n_neighbors, n_bonds, coord = atom
        cx, cy = coord
        r = vp.radius_px(atom_radii[a_number] / 1.5)
        color = atom_colors[a_number]
        px, py = vp.to_px(cx, cy)
        circle = dwg.circle(center=(px, py), r=r, fill=color, stroke="none")
        dwg.add(circle)


def draw_atoms_stroke(dwg, atoms, atom_radii, atom_colors, bond_color, vp: PixelViewport):
    for atom in atoms:
        a_name, a_id, a_number, n_neighbors, n_bonds, coord = atom
        cx, cy = coord
        r = vp.radius_px(atom_radii[a_number] / 1.5)
        color = atom_colors[a_number]
        px, py = vp.to_px(cx, cy)
        circle = dwg.circle(center=(px, py), r=r, fill=color, stroke=bond_color, stroke_width=vp.bond_width_px)
        dwg.add(circle)


def draw_atoms_hetero_names(dwg, atoms, atom_radii, atom_colors, font, vp: PixelViewport):
    fs = vp.label_font_px
    for atom in atoms:
        a_name, a_id, a_number, n_neighbors, n_bonds, coord = atom
        cx, cy = coord
        r_mol = atom_radii[a_number] / 1.5
        color = atom_colors[a_number]

        if a_name != "C":
            r_mol += 0.1
            r = vp.radius_px(r_mol)
            px, py = vp.to_px(cx, cy)
            circle = dwg.circle(center=(px, py), r=r, fill="#ffffff", stroke=color, stroke_width=vp.bond_width_px)
            text = _atom_symbol_text(dwg, a_name, px, py, color, font, fs)
            dwg.add(circle)
            dwg.add(text)
        else:
            r = vp.radius_px(r_mol)
            px, py = vp.to_px(cx, cy)
            circle = dwg.circle(center=(px, py), r=r, fill=color, stroke="none")
            dwg.add(circle)


def draw_atoms_all_names(dwg, atoms, atom_radii, atom_colors, font, vp: PixelViewport):
    fs = vp.label_font_px
    for atom in atoms:
        a_name, a_id, a_number, n_neighbors, n_bonds, coord = atom
        cx, cy = coord
        r_mol = atom_radii[a_number] / 1.5 + 0.1
        color = atom_colors[a_number]
        r = vp.radius_px(r_mol)
        px, py = vp.to_px(cx, cy)
        circle = dwg.circle(center=(px, py), r=r, fill="#ffffff", stroke=color, stroke_width=vp.bond_width_px)
        text = _atom_symbol_text(dwg, a_name, px, py, color, font, fs)
        dwg.add(circle)
        dwg.add(text)


def add_hydrogens_to_donors(dwg, atoms, bonds, atom_colors, atom_radii, style, bond_color, vp: PixelViewport):
    if style == "plain":
        r_extra = 0
        s_width = 0
    elif style in ("names_hetero", "names_all"):
        r_extra = 0.1
        s_width = 0
    elif style == "stroke":
        r_extra = 0
        s_width = vp.bond_width_px
    else:
        r_extra = 0
        s_width = 0

    h_r_mol = atom_radii[1] / 1.5
    h_r = vp.radius_px(h_r_mol)
    h_color = atom_colors[1]

    for atom in atoms:
        a_name, a_id, a_number, n_bonds, n_neighbors, coord = atom

        if (a_name in ["O", "S"]) and (n_bonds == 1) and (n_neighbors == 1):
            a_bond = [x for x in bonds if (x[-2] == a_id) or (x[-1] == a_id)][0]
            a_r_mol = atom_radii[a_number] / 1.5

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
            rot = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
            rotated = np.dot(rot, vec)
            d = a_r_mol + h_r_mol + r_extra
            length = np.linalg.norm(rotated)
            final_vec = (rotated / length) * d
            h_center = np.array(start) + final_vec
            px, py = vp.to_px(float(h_center[0]), float(h_center[1]))
            circle = dwg.circle(center=(px, py), r=h_r, fill=h_color, stroke=bond_color, stroke_width=s_width)
            dwg.add(circle)

        if (a_name in ["N", "P"]) and (n_bonds == 2) and (n_neighbors == 1):
            a_bond = [x for x in bonds if (x[-2] == a_id) or (x[-1] == a_id)][0]
            a_r_mol = atom_radii[a_number] / 1.5

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
            rot = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
            rotated = np.dot(rot, vec)
            d = a_r_mol + h_r_mol + r_extra
            length = np.linalg.norm(rotated)
            final_vec = (rotated / length) * d
            h_center = np.array(start) + final_vec
            px, py = vp.to_px(float(h_center[0]), float(h_center[1]))
            circle = dwg.circle(center=(px, py), r=h_r, fill=h_color, stroke=bond_color, stroke_width=s_width)
            dwg.add(circle)

        if (a_name in ["N", "P"]) and (n_bonds == 1) and (n_neighbors == 1):
            a_bond = [x for x in bonds if (x[-2] == a_id) or (x[-1] == a_id)][0]
            a_r_mol = atom_radii[a_number] / 1.5
            d = a_r_mol + h_r_mol + r_extra

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
            rot1 = np.array([[math.cos(angle1), -math.sin(angle1)], [math.sin(angle1), math.cos(angle1)]])
            rotated1 = np.dot(rot1, vec)
            length = np.linalg.norm(rotated1)
            final_vec1 = (rotated1 / length) * d
            h_center1 = np.array(start) + final_vec1
            px1, py1 = vp.to_px(float(h_center1[0]), float(h_center1[1]))
            circle1 = dwg.circle(center=(px1, py1), r=h_r, fill=h_color, stroke=bond_color, stroke_width=s_width)
            dwg.add(circle1)

            angle2 = -sign * np.deg2rad(120)
            rot2 = np.array([[math.cos(angle2), -math.sin(angle2)], [math.sin(angle2), math.cos(angle2)]])
            rotated2 = np.dot(rot2, vec)
            length = np.linalg.norm(rotated2)
            final_vec2 = (rotated2 / length) * d
            h_center2 = np.array(start) + final_vec2
            px2, py2 = vp.to_px(float(h_center2[0]), float(h_center2[1]))
            circle2 = dwg.circle(center=(px2, py2), r=h_r, fill=h_color, stroke=bond_color, stroke_width=s_width)
            dwg.add(circle2)

        if (a_name in ["N", "P"]) and (n_bonds == 2) and (n_neighbors == 2):

            a_bond = [x for x in bonds if (x[-2] == a_id) or (x[-1] == a_id)]
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

            vec_bis = -1 * (vec1 + vec2) / 2

            a_r_mol = atom_radii[a_number] / 1.5
            d = a_r_mol + h_r_mol + r_extra
            length = np.linalg.norm(vec_bis)

            final_vec = (vec_bis / length) * d
            h_center = np.array(start) + final_vec
            px, py = vp.to_px(float(h_center[0]), float(h_center[1]))
            circle = dwg.circle(center=(px, py), r=h_r, fill=h_color, stroke=bond_color, stroke_width=s_width)
            dwg.add(circle)
