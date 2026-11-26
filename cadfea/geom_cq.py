"""The models in CadQuery: sketch a profile, extrude it, fillet an edge, cut a hole.

CadQuery is a fluent API: each call returns a new Workplane holding a stack of
objects, and selectors (strings like "|Y" or selector objects) pick edges and faces.
"""
from __future__ import annotations

from pathlib import Path

import cadquery as cq

from .params import Beam, Bracket


def beam(p: Beam) -> cq.Workplane:
    # The "YZ" workplane's normal is +X, so the extrusion runs from x=0 to x=L.
    return cq.Workplane("YZ").rect(p.b, p.h).extrude(p.L)


def bracket(p: Bracket) -> cq.Workplane:
    # The "XZ" workplane maps local (x, y) to global (X, Z); its normal is -Y.
    # both=True extrudes W/2 each way, so the part is centred on y = 0.
    profile = [(0, 0), (p.L, 0), (p.L, p.t), (p.t, p.t), (p.t, p.H), (0, p.H)]
    body = cq.Workplane("XZ").polyline(profile).close().extrude(p.W / 2, both=True)
    # The inner corner is the edge parallel to Y nearest the point (t, 0, t).
    body = (
        body.edges("|Y")
        .edges(cq.selectors.NearestToPointSelector((p.t, 0.0, p.t)))
        .fillet(p.r)
    )
    # The cutting cylinder overshoots both faces by 1 mm: coplanar faces in a
    # boolean are a classic source of slivers and failed cuts.
    cutter = (
        cq.Workplane("XY")
        .workplane(offset=-1.0)
        .center(p.xh, 0.0)
        .circle(p.d / 2)
        .extrude(p.t + 2.0)
    )
    return body.cut(cutter)


def build(p: Beam | Bracket) -> cq.Workplane:
    return bracket(p) if isinstance(p, Bracket) else beam(p)


def export_step(wp: cq.Workplane, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    cq.exporters.export(wp, str(path), exportType="STEP")
    return path


def properties(wp: cq.Workplane) -> dict:
    """Mass properties through CadQuery's own API (unit density)."""
    solid = wp.val()
    com = cq.Shape.centerOfMass(solid)
    bb = solid.BoundingBox()  # CadQuery uses BRepBndLib.AddOptimal here
    inertia = cq.Shape.matrixOfInertia(solid)
    return {
        "volume": solid.Volume(),
        "area": solid.Area(),
        "com": [com.x, com.y, com.z],
        "bbox": [bb.xmin, bb.ymin, bb.zmin, bb.xmax, bb.ymax, bb.zmax],
        "inertia": [list(row) for row in inertia],
        "n_faces": len(solid.Faces()),
        "n_edges": len(solid.Edges()),
    }
