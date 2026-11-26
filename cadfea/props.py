"""Mass properties three ways: read back from STEP, by closed form, and compared.

Reading the STEP files uses OCP (the Python binding of OpenCASCADE, the kernel
under both CadQuery and build123d) directly, so neither library's own code is
involved in the round trip. The closed-form values need no CAD at all.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .params import Beam, Bracket


def read_step(path: Path):
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Reader

    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_RetDone:
        raise RuntimeError(f"cannot read {path}")
    reader.TransferRoots()
    return reader.OneShape()


def step_properties(path: Path) -> dict:
    """Volume, area, centre of mass, inertia and bounding box of a STEP file."""
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer

    shape = read_step(path)
    vol = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, vol)
    surf = GProp_GProps()
    BRepGProp.SurfaceProperties_s(shape, surf)
    com = vol.CentreOfMass()
    m = vol.MatrixOfInertia()
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, box, False, False)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()

    def count(kind) -> int:
        exp, n = TopExp_Explorer(shape, kind), 0
        while exp.More():
            n += 1
            exp.Next()
        return n

    return {
        "volume": vol.Mass(),
        "area": surf.Mass(),
        "com": [com.X(), com.Y(), com.Z()],
        "bbox": [xmin, ymin, zmin, xmax, ymax, zmax],
        "inertia": [[m.Value(i, j) for j in (1, 2, 3)] for i in (1, 2, 3)],
        "n_solids": count(TopAbs_SOLID),
        "n_faces": count(TopAbs_FACE),
    }


def analytic(p: Beam | Bracket) -> dict:
    """Closed-form volume, surface area and centre of mass."""
    if isinstance(p, Beam):
        I_y = p.b * p.h**3 / 12.0
        return {
            "volume": p.L * p.b * p.h,
            "area": 2 * (p.L * p.b + p.L * p.h + p.b * p.h),
            "com": [p.L / 2, 0.0, 0.0],
            "bbox": [0.0, -p.b / 2, -p.h / 2, p.L, p.b / 2, p.h / 2],
            "I_section": I_y,
        }
    L, W, t, H, r, d = p.L, p.W, p.t, p.H, p.r, p.d
    # Profile in the x-z plane: arm + upright above the arm + the fillet's
    # spandrel (a square of side r minus a quarter disc of radius r).
    arm = (L * t, L / 2, t / 2)
    upright = (t * (H - t), t / 2, (t + H) / 2)
    s_area = r * r * (1 - math.pi / 4)
    s_off = r * (10 - 3 * math.pi) / (12 - 3 * math.pi)  # spandrel centroid from the corner
    spandrel = (s_area, t + s_off, t + s_off)
    parts = [arm, upright, spandrel]
    A = sum(a for a, _, _ in parts)
    hole_v = math.pi * (d / 2) ** 2 * t
    volume = A * W - hole_v
    xc = (sum(a * x for a, x, _ in parts) * W - hole_v * p.xh) / volume
    zc = (sum(a * z for a, _, z in parts) * W - hole_v * t / 2) / volume
    perimeter = L + t + (L - t - r) + math.pi * r / 2 + (H - t - r) + t + H
    area = 2 * A + perimeter * W - 2 * math.pi * (d / 2) ** 2 + math.pi * d * t
    return {
        "volume": volume,
        "area": area,
        "com": [xc, 0.0, zc],
        "bbox": [0.0, -W / 2, 0.0, L, W / 2, H],
    }


def rel_diff(a: dict, b: dict) -> dict:
    """Largest relative difference per property between two property dicts.

    Vectors and matrices are compared relative to their largest entry, so a
    centre of mass at y = 1e-15 against y = 0 does not count as a 100% error.
    """
    out = {}
    for key in ("volume", "area", "com", "bbox", "inertia"):
        if key not in a or key not in b:
            continue
        x, y = np.asarray(a[key], float), np.asarray(b[key], float)
        scale = max(np.max(np.abs(x)), np.max(np.abs(y)), 1e-300)
        out[key] = float(np.max(np.abs(x - y)) / scale)
    return out
