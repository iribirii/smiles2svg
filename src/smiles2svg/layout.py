"""Pixel canvas fitting: map RDKit 2D coordinates (Å) to SVG pixel space.

Fitting
  scale = (canvas_size − 2·padding) / max(span_x, span_y)   →  px per Å

Stroke and labels (canvas- and span-independent proportions)
  bond_width_px  = max(MIN, BOND_STROKE_ANGSTROM  × scale)
  label_font_px  = max(MIN, LABEL_FONT_ANGSTROM   × scale)

So doubling ``canvas_size`` doubles ``scale`` and doubles stroke/font in px — the drawing
scales uniformly. Changing only the molecule span changes ``scale`` so stroke and atom
radii stay in the same ratio to each other (unlike a fixed reference span of 6 Å).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from smiles2svg.molecule import AtomRecord, BondRecord

# Stroke width and label font size in Å (same units as RDKit 2D coordinates).
# Pixel size = factor × scale, where scale = (canvas_inner_px) / max_span_Å.
# This is *canvas-independent* and *span-independent*: bond/atom radius and text
# stay in fixed proportion to the geometry (like legacy unitless 0.1 / 0.8 in user space).
BOND_STROKE_ANGSTROM = 0.10
LABEL_FONT_ANGSTROM = 0.80

# Minimum readable sizes in px (when scale is small for very large molecules).
MIN_BOND_STROKE_PX = 0.8
MIN_LABEL_FONT_PX = 7.0

# Multi-bond perpendicular offset in molecule space (unchanged from legacy smiles2svg).
BOND_DOUBLE_OFFSET = 0.2
BOND_TRIPLE_OFFSET = 0.25

# Extra margin (Å) around scene when implicit H are drawn (conservative).
HYDROGEN_MARGIN = 1.8


def displace_perpendicular_mol(
    start: tuple[float, float], end: tuple[float, float], dist: float
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Offset bond line perpendicular by dist in molecule coordinates."""
    original_vector = np.array(end) - np.array(start)
    perpendicular_vector = np.array([-original_vector[1], original_vector[0]])
    perpendicular_vector = perpendicular_vector / np.linalg.norm(perpendicular_vector)
    new_start = np.array(start) + perpendicular_vector * dist
    new_end = np.array(end) + perpendicular_vector * dist
    ns = new_start.tolist()
    ne = new_end.tolist()
    return (float(ns[0]), float(ns[1])), (float(ne[0]), float(ne[1]))


def _effective_atom_radius(atom: AtomRecord, draw_style: str, atom_radii: list[float]) -> float:
    r = atom_radii[atom.atomic_num] / 1.5
    if draw_style == "names_all":
        r += 0.1
    elif draw_style == "names_hetero" and atom.symbol != "C":
        r += 0.1
    return r


def compute_scene_bounds(
    atoms: list[AtomRecord],
    bonds: list[BondRecord],
    draw_style: str,
    atom_radii: list[float],
    hydrogens: bool,
) -> tuple[np.ndarray, np.ndarray]:
    """Axis-aligned bounds (min_xy, max_xy) in molecule coordinates."""
    lo = np.array([math.inf, math.inf])
    hi = np.array([-math.inf, -math.inf])

    def expand_point(x: float, y: float) -> None:
        lo[0] = min(lo[0], x)
        lo[1] = min(lo[1], y)
        hi[0] = max(hi[0], x)
        hi[1] = max(hi[1], y)

    for atom in atoms:
        cx, cy = atom.coord
        r = _effective_atom_radius(atom, draw_style, atom_radii)
        expand_point(cx - r, cy - r)
        expand_point(cx + r, cy + r)

    for bond in bonds:
        bond_order = bond.order
        bond_start, bond_end = bond.start_xy, bond.end_xy
        expand_point(bond_start[0], bond_start[1])
        expand_point(bond_end[0], bond_end[1])
        if bond_order == 2 or bond_order == 1.5:
            su, eu = displace_perpendicular_mol(bond_start, bond_end, BOND_DOUBLE_OFFSET)
            sd, ed = displace_perpendicular_mol(bond_start, bond_end, -BOND_DOUBLE_OFFSET)
            expand_point(su[0], su[1])
            expand_point(eu[0], eu[1])
            expand_point(sd[0], sd[1])
            expand_point(ed[0], ed[1])
        elif bond_order == 3:
            su, eu = displace_perpendicular_mol(bond_start, bond_end, BOND_TRIPLE_OFFSET)
            sd, ed = displace_perpendicular_mol(bond_start, bond_end, -BOND_TRIPLE_OFFSET)
            expand_point(su[0], su[1])
            expand_point(eu[0], eu[1])
            expand_point(sd[0], sd[1])
            expand_point(ed[0], ed[1])

    if hydrogens:
        lo -= HYDROGEN_MARGIN
        hi += HYDROGEN_MARGIN

    if not np.isfinite(lo).all():
        lo[:] = 0.0
        hi[:] = 1.0
    span = max(hi[0] - lo[0], hi[1] - lo[1], 1e-9)
    if span < 1e-6:
        lo[:] = -0.5
        hi[:] = 0.5

    return lo, hi


def fit_canvas(
    lo: np.ndarray,
    hi: np.ndarray,
    canvas_size: int,
    padding: float,
) -> tuple[float, float, float, int, int]:
    """
    Return scale (px per Å), center (cx_mol, cy_mol), canvas_w, canvas_h.

    The larger molecule span (Å) is fitted to ``canvas_size - 2*padding`` px along that axis.
    """
    spans = hi - lo
    max_span = max(float(spans.max()), 1e-9)
    scale = (canvas_size - 2 * padding) / max_span
    cx_mol = float((lo[0] + hi[0]) / 2)
    cy_mol = float((lo[1] + hi[1]) / 2)
    canvas_w = max(1, int(np.ceil(spans[0] * scale + 2 * padding)))
    canvas_h = max(1, int(np.ceil(spans[1] * scale + 2 * padding)))
    return scale, cx_mol, cy_mol, canvas_w, canvas_h


@dataclass(frozen=True)
class PixelViewport:
    """Maps molecule coordinates (Å) to pixel SVG space (origin top-left, y down)."""

    scale: float  # px per Å (uniform x and y)
    cx_mol: float
    cy_mol: float
    canvas_w: int
    canvas_h: int
    padding: float
    bond_width_px: float
    label_font_px: float

    def to_px(self, x: float, y: float) -> tuple[float, float]:
        x_px = self.canvas_w / 2 + self.scale * (x - self.cx_mol)
        y_px = self.canvas_h / 2 + self.scale * (y - self.cy_mol)
        return x_px, y_px

    def radius_px(self, r_mol: float) -> float:
        return r_mol * self.scale


def make_viewport(
    lo: np.ndarray,
    hi: np.ndarray,
    canvas_size: int,
    padding: float,
) -> PixelViewport:
    scale, cx_mol, cy_mol, canvas_w, canvas_h = fit_canvas(lo, hi, canvas_size, padding)
    bond_w = BOND_STROKE_ANGSTROM * scale
    label = LABEL_FONT_ANGSTROM * scale
    return PixelViewport(
        scale=scale,
        cx_mol=cx_mol,
        cy_mol=cy_mol,
        canvas_w=canvas_w,
        canvas_h=canvas_h,
        padding=padding,
        bond_width_px=max(MIN_BOND_STROKE_PX, bond_w),
        label_font_px=max(MIN_LABEL_FONT_PX, label),
    )
