"""The same models in build123d, built a different way on purpose.

CadQuery sketches the L profile and extrudes it; here the bracket is the union of
two boxes (constructive solid geometry), then the same fillet and hole. If both
routes give the same solid, the parameters mean what I think they mean.

build123d has a builder mode (context managers that collect objects) and an
algebra mode (Part + Part). This file uses builder mode.
"""
from __future__ import annotations

from pathlib import Path

from build123d import (
    Align,
    Axis,
    Box,
    BuildPart,
    CenterOf,
    Cylinder,
    Locations,
    Mode,
    Part,
    export_step as _export_step,
    fillet,
)

from .params import Beam, Bracket

MIN, CEN = Align.MIN, Align.CENTER


def beam(p: Beam) -> Part:
    with BuildPart() as bp:
        Box(p.L, p.b, p.h, align=(MIN, CEN, CEN))
    return bp.part


def bracket(p: Bracket) -> Part:
    with BuildPart() as bp:
        Box(p.L, p.W, p.t, align=(MIN, CEN, MIN))  # arm
        Box(p.t, p.W, p.H, align=(MIN, CEN, MIN))  # upright, fused with the arm
        corner = bp.edges().filter_by(Axis.Y).sort_by_distance((p.t, 0.0, p.t))[0]
        fillet(corner, radius=p.r)
        with Locations((p.xh, 0.0, p.t / 2)):
            Cylinder(p.d / 2, p.t + 2.0, mode=Mode.SUBTRACT)
    return bp.part


def build(p: Beam | Bracket) -> Part:
    return bracket(p) if isinstance(p, Bracket) else beam(p)


def export_step(part: Part, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    _export_step(part, str(path))
    return path


def properties(part: Part) -> dict:
    """Mass properties through build123d's own API (unit density)."""
    com = part.center(CenterOf.MASS)
    bb = part.bounding_box(optimal=True)
    return {
        "volume": part.volume,
        "area": part.area,
        "com": [com.X, com.Y, com.Z],
        "bbox": [bb.min.X, bb.min.Y, bb.min.Z, bb.max.X, bb.max.Y, bb.max.Z],
        "inertia": [list(row) for row in part.matrix_of_inertia],
        "n_faces": len(part.faces()),
        "n_edges": len(part.edges()),
    }
