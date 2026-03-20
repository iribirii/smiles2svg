"""SVG drawing of 2D molecular graphs (pixel canvas, xyzrender-style scaling)."""

from __future__ import annotations

import math
from pathlib import Path

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
from smiles2svg.molecule import AtomRecord, BondRecord

PathLike = str | Path

# With ``font-size`` in px, many SVG stacks (incl. CairoSVG) ignore or mishandle
# ``dominant-baseline="central"``. Use the alphabetic baseline and shift Y so the
# visual center of typical Latin capitals sits on (cx, cy): cap height ~0.7em.
_LABEL_BASELINE_SHIFT_EM = 0.35


def _bonds_incident(atom_idx: int, bonds: list[BondRecord]) -> list[BondRecord]:
    return [b for b in bonds if atom_idx in (b.begin_atom_idx, b.end_atom_idx)]


def _start_and_vec_from_atom(
    bond: BondRecord, atom_idx: int
) -> tuple[tuple[float, float], np.ndarray]:
    if atom_idx == bond.begin_atom_idx:
        start_xy, neigh_xy = bond.start_xy, bond.end_xy
    else:
        start_xy, neigh_xy = bond.end_xy, bond.start_xy
    vec = np.array(neigh_xy, dtype=float) - np.array(start_xy, dtype=float)
    return start_xy, vec


def _geom_sign(vec: np.ndarray) -> float:
    g = float(np.sign(vec[0]) * np.sign(vec[1]))
    return 1.0 if abs(g) < 1e-15 else g


def _rot2(angle_rad: float) -> np.ndarray:
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[c, -s], [s, c]])


def _h_disk_from_start(
    start_xy: tuple[float, float], vec: np.ndarray, angle_rad: float, d: float
) -> tuple[float, float]:
    rotated = _rot2(angle_rad) @ vec
    ln = np.linalg.norm(rotated)
    final = (rotated / ln) * d
    h = np.array(start_xy, dtype=float) + final
    return float(h[0]), float(h[1])


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
    atoms: list[AtomRecord],
    bonds: list[BondRecord],
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
        add_hydrogens_to_donors(
            dwg, atoms, bonds, atom_colors, atom_radii, draw_style, bond_color, vp
        )

    dwg.save()
    return dwg.tostring()


def draw_bonds(
    dwg: svgwrite.Drawing, bonds: list[BondRecord], color: str, vp: PixelViewport
) -> None:
    bw = vp.bond_width_px
    for bond in bonds:
        bond_order = bond.order
        bond_start, bond_end = bond.start_xy, bond.end_xy
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


def draw_atoms_plain(
    dwg: svgwrite.Drawing,
    atoms: list[AtomRecord],
    atom_radii: list[float],
    atom_colors: list[str],
    vp: PixelViewport,
) -> None:
    for atom in atoms:
        cx, cy = atom.coord
        r = vp.radius_px(atom_radii[atom.atomic_num] / 1.5)
        color = atom_colors[atom.atomic_num]
        px, py = vp.to_px(cx, cy)
        circle = dwg.circle(center=(px, py), r=r, fill=color, stroke="none")
        dwg.add(circle)


def draw_atoms_stroke(
    dwg: svgwrite.Drawing,
    atoms: list[AtomRecord],
    atom_radii: list[float],
    atom_colors: list[str],
    bond_color: str,
    vp: PixelViewport,
) -> None:
    for atom in atoms:
        cx, cy = atom.coord
        r = vp.radius_px(atom_radii[atom.atomic_num] / 1.5)
        color = atom_colors[atom.atomic_num]
        px, py = vp.to_px(cx, cy)
        circle = dwg.circle(
            center=(px, py),
            r=r,
            fill=color,
            stroke=bond_color,
            stroke_width=vp.bond_width_px,
        )
        dwg.add(circle)


def draw_atoms_hetero_names(
    dwg: svgwrite.Drawing,
    atoms: list[AtomRecord],
    atom_radii: list[float],
    atom_colors: list[str],
    font: str,
    vp: PixelViewport,
) -> None:
    fs = vp.label_font_px
    for atom in atoms:
        cx, cy = atom.coord
        r_mol = atom_radii[atom.atomic_num] / 1.5
        color = atom_colors[atom.atomic_num]

        if atom.symbol != "C":
            r_mol += 0.1
            r = vp.radius_px(r_mol)
            px, py = vp.to_px(cx, cy)
            circle = dwg.circle(
                center=(px, py),
                r=r,
                fill="#ffffff",
                stroke=color,
                stroke_width=vp.bond_width_px,
            )
            text = _atom_symbol_text(dwg, atom.symbol, px, py, color, font, fs)
            dwg.add(circle)
            dwg.add(text)
        else:
            r = vp.radius_px(r_mol)
            px, py = vp.to_px(cx, cy)
            circle = dwg.circle(center=(px, py), r=r, fill=color, stroke="none")
            dwg.add(circle)


def draw_atoms_all_names(
    dwg: svgwrite.Drawing,
    atoms: list[AtomRecord],
    atom_radii: list[float],
    atom_colors: list[str],
    font: str,
    vp: PixelViewport,
) -> None:
    fs = vp.label_font_px
    for atom in atoms:
        cx, cy = atom.coord
        r_mol = atom_radii[atom.atomic_num] / 1.5 + 0.1
        color = atom_colors[atom.atomic_num]
        r = vp.radius_px(r_mol)
        px, py = vp.to_px(cx, cy)
        circle = dwg.circle(
            center=(px, py),
            r=r,
            fill="#ffffff",
            stroke=color,
            stroke_width=vp.bond_width_px,
        )
        text = _atom_symbol_text(dwg, atom.symbol, px, py, color, font, fs)
        dwg.add(circle)
        dwg.add(text)


def add_hydrogens_to_donors(
    dwg: svgwrite.Drawing,
    atoms: list[AtomRecord],
    bonds: list[BondRecord],
    atom_colors: list[str],
    atom_radii: list[float],
    style: str,
    bond_color: str,
    vp: PixelViewport,
) -> None:
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

    def put_h(mx: float, my: float) -> None:
        px, py = vp.to_px(mx, my)
        dwg.add(
            dwg.circle(
                center=(px, py),
                r=h_r,
                fill=h_color,
                stroke=bond_color,
                stroke_width=s_width,
            )
        )

    for atom in atoms:
        sym = atom.symbol
        aid = atom.atom_idx
        z = atom.atomic_num
        bond_order_sum = atom.bond_order_sum
        degree = atom.degree

        if (sym in ("O", "S")) and (bond_order_sum == 1) and (degree == 1):
            bond = _bonds_incident(aid, bonds)[0]
            a_r_mol = atom_radii[z] / 1.5
            start_xy, vec = _start_and_vec_from_atom(bond, aid)
            sign = _geom_sign(vec)
            d = a_r_mol + h_r_mol + r_extra
            hx, hy = _h_disk_from_start(start_xy, vec, sign * np.deg2rad(120), d)
            put_h(hx, hy)

        if (sym in ("N", "P")) and (bond_order_sum == 2) and (degree == 1):
            bond = _bonds_incident(aid, bonds)[0]
            a_r_mol = atom_radii[z] / 1.5
            start_xy, vec = _start_and_vec_from_atom(bond, aid)
            sign = _geom_sign(vec)
            d = a_r_mol + h_r_mol + r_extra
            hx, hy = _h_disk_from_start(start_xy, vec, sign * np.deg2rad(120), d)
            put_h(hx, hy)

        if (sym in ("N", "P")) and (bond_order_sum == 1) and (degree == 1):
            bond = _bonds_incident(aid, bonds)[0]
            a_r_mol = atom_radii[z] / 1.5
            d = a_r_mol + h_r_mol + r_extra
            start_xy, vec = _start_and_vec_from_atom(bond, aid)
            sign = _geom_sign(vec)
            hx1, hy1 = _h_disk_from_start(start_xy, vec, sign * np.deg2rad(120), d)
            hx2, hy2 = _h_disk_from_start(start_xy, vec, -sign * np.deg2rad(120), d)
            put_h(hx1, hy1)
            put_h(hx2, hy2)

        if (sym in ("N", "P")) and (bond_order_sum == 2) and (degree == 2):
            b1, b2 = _bonds_incident(aid, bonds)
            _, vec1 = _start_and_vec_from_atom(b1, aid)
            _, vec2 = _start_and_vec_from_atom(b2, aid)
            # ``start`` is the donor atom position (same for both incident bonds)
            if aid == b1.begin_atom_idx:
                start = np.array(b1.start_xy, dtype=float)
            else:
                start = np.array(b1.end_xy, dtype=float)
            vec_bis = -0.5 * (vec1 + vec2)
            a_r_mol = atom_radii[z] / 1.5
            d = a_r_mol + h_r_mol + r_extra
            ln = np.linalg.norm(vec_bis)
            final_vec = (vec_bis / ln) * d
            h_center = start + final_vec
            put_h(float(h_center[0]), float(h_center[1]))
