"""Build each output once and reuse it while its inputs are unchanged.

Scripts and tests both call these functions. Every output records a fingerprint
of what it was made from (the parameters for geometry and meshes, the mesh file
for solver results), so a stale file is rebuilt instead of trusted.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from . import paths
from .params import CASES, MESH_LEVELS, STEEL, Material, as_dict


def fingerprint(*parts) -> str:
    h = hashlib.sha1()
    for part in parts:
        if isinstance(part, Path):
            h.update(part.read_bytes())
        else:
            h.update(json.dumps(part, sort_keys=True).encode())
    return h.hexdigest()[:16]


def _read(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


# ---------------------------------------------------------------- geometry

def ensure_step(case: str, tool: str = "cadquery") -> Path:
    """STEP file for a case from one CAD tool (built if missing or stale)."""
    p = CASES[case]
    path = paths.step_path(case, tool)
    stamp = path.with_suffix(".json")
    fp = fingerprint(as_dict(p), tool)
    if path.exists() and (_read(stamp) or {}).get("fingerprint") == fp:
        return path
    if tool == "cadquery":
        from . import geom_cq as g
    else:
        from . import geom_b123d as g
    model = g.build(p)
    g.export_step(model, path)
    stamp.write_text(json.dumps({"fingerprint": fp, "properties": g.properties(model)}))
    return path


# ---------------------------------------------------------------- meshes

def mesh_fingerprint(case: str, level: str) -> str:
    return fingerprint(as_dict(CASES[case]), as_dict(MESH_LEVELS[case][level]), "mesh-v1")


def ensure_mesh(case: str, level: str) -> dict:
    """Gmsh mesh of a case at a level, from the CadQuery STEP file."""
    from .mesh import build_mesh, load_info

    fp = mesh_fingerprint(case, level)
    sidecar = paths.MESH_DIR / f"{case}_{level}.json"
    have = _read(sidecar)
    if have and have.get("fingerprint") == fp and (paths.MESH_DIR / f"{case}_{level}.msh").exists():
        return load_info(paths.MESH_DIR, case, level)
    info = build_mesh(case, level, ensure_step(case, "cadquery"), paths.MESH_DIR)
    info["fingerprint"] = fp
    sidecar.write_text(json.dumps(info))
    return load_info(paths.MESH_DIR, case, level)


def mesh_id(case: str, level: str) -> str:
    """Identifies the exact mesh a result was computed on."""
    return fingerprint(paths.MESH_DIR / f"{case}_{level}.msh")


# ---------------------------------------------------------------- solvers

def _mat_tag(mat: Material) -> str:
    return "" if mat == STEEL else f"_nu{mat.nu:g}".replace(".", "")


def _result_path(solver: str, case: str, level: str, element: str, mat: Material = STEEL) -> Path:
    default = {"calculix": "C3D10", "fenicsx": "P2"}[solver]
    tag = "" if element == default else f"_{element.lower()}"
    short = {"calculix": "ccx", "fenicsx": "fenicsx"}[solver]
    return paths.RESULTS_DIR / f"{short}_{case}_{level}{tag}{_mat_tag(mat)}.json"


def load_result(solver: str, case: str, level: str, element: str | None = None,
                mat: Material = STEEL) -> dict | None:
    element = element or {"calculix": "C3D10", "fenicsx": "P2"}[solver]
    return _read(_result_path(solver, case, level, element, mat))


def ensure_ccx(case: str, level: str, element: str = "C3D10", mat: Material = STEEL) -> dict:
    from . import ccx

    info = ensure_mesh(case, level)
    mid = mesh_id(case, level)
    path = _result_path("calculix", case, level, element, mat)
    have = _read(path)
    if have and have.get("mesh_id") == mid and have.get("material") == as_dict(mat):
        return have
    res = ccx.solve(info, paths.MESH_DIR, paths.CCX_DIR, element, mat)
    res.update(mesh_id=mid, material=as_dict(mat), size=info["size"],
               feature_size=info["feature_size"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1))
    return res


def ensure_fenicsx(case: str, level: str, degree: int = 2, solver: str = "lu") -> dict:
    from . import fenics

    info = ensure_mesh(case, level)
    mid = mesh_id(case, level)
    path = _result_path("fenicsx", case, level, f"P{degree}")
    have = _read(path)
    if have and have.get("mesh_id") == mid and have.get("material") == as_dict(STEEL):
        return have
    res = fenics.solve(paths.MESH_DIR / f"{case}_{level}.msh", info, degree=degree, solver=solver)
    res.update(mesh_id=mid, material=as_dict(STEEL), size=info["size"],
               feature_size=info["feature_size"], linear_solver=solver)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1))
    return res
