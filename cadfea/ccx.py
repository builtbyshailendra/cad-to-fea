"""CalculiX: write the input deck, run ccx, read the .dat and .frd results.

The mesh (nodes and C3D10 elements) comes from Gmsh's Abaqus export. This module
adds what a mesh file doesn't say: node sets, material, section, boundary
conditions, loads, the analysis step and which results to write.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

from .mesh import load_vectors
from .params import CASES, STEEL, Bracket, Material

ELEMENT_TYPES = ("C3D10", "C3D4")


def load_cases(case: str) -> list[str]:
    return ["shelf", "pull"] if isinstance(CASES[case], Bracket) else ["tip"]


def _nset(name: str, nodes: list[int]) -> str:
    lines = [f"*NSET, NSET={name}"]
    for i in range(0, len(nodes), 12):
        lines.append(", ".join(str(n) for n in nodes[i : i + 12]))
    return "\n".join(lines)


def _linear_mesh(mesh_inp: Path, out: Path) -> set[int]:
    """Rewrite a C3D10 mesh as C3D4 by dropping the six mid-edge nodes per element.

    Same corner nodes, same tets: only the element order changes. Returns the
    corner nodes (the only nodes the C3D4 elements use).
    """
    lines = mesh_inp.read_text().splitlines()
    res, in_elem, corners = [], False, set()
    for line in lines:
        if line.startswith("*"):
            in_elem = line.upper().startswith("*ELEMENT")
            res.append(line.replace("C3D10", "C3D4"))
            continue
        if in_elem and line.strip():
            fields = [f.strip() for f in line.split(",")]
            res.append(", ".join(fields[:5]))
            corners.update(int(f) for f in fields[1:5])
        else:
            res.append(line)
    out.write_text("\n".join(res) + "\n")
    return corners


def write_deck(info: dict, mesh_dir: Path, job_dir: Path, element: str = "C3D10",
               mat: Material = STEEL) -> Path:
    """Write <job>.inp next to a copy of the mesh, return its path."""
    case, level = info["case"], info["level"]
    job = f"{case}_{level}" + ("" if element == "C3D10" else f"_{element.lower()}")
    job_dir.mkdir(parents=True, exist_ok=True)
    src = mesh_dir / f"{case}_{level}_mesh.inp"
    mesh_file = job_dir / f"{job}_mesh.inp"
    nsets = info["nsets"]
    if element == "C3D10":
        shutil.copyfile(src, mesh_file)
    elif element == "C3D4":
        # Mid-edge nodes stay in the node list but belong to no element; keep
        # them out of the sets, or CalculiX would constrain or load free nodes.
        corners = _linear_mesh(src, mesh_file)
        nsets = {k: [n for n in v if n in corners] for k, v in nsets.items()}
    else:
        raise ValueError(f"element must be one of {ELEMENT_TYPES}")

    tip = info["tip_node"]
    parts = [
        f"** {case}, {level} mesh, {element}: written by cadfea/ccx.py",
        f"*INCLUDE, INPUT={mesh_file.name}",
    ]
    for name, nodes in nsets.items():
        parts.append(_nset(f"N{name}", nodes))
    parts.append(_nset("NTIP", [tip]))
    parts += [
        f"*MATERIAL, NAME={mat.name}",
        "*ELASTIC",
        f"{mat.E:.1f}, {mat.nu}",
        f"*SOLID SECTION, ELSET=SOLID, MATERIAL={mat.name}",
        # Every degree of freedom of every node on the root face is held: a
        # perfectly rigid wall. Outside a step, it applies to all steps.
        "*BOUNDARY",
        "NFIXED, 1, 3, 0.0",
    ]
    for i, lc in enumerate(load_cases(case)):
        loads = info["loads"][lc] if element == "C3D10" else _linear_loads(info, lc)
        parts += [
            f"** load case {lc}",
            "*STEP",
            "*STATIC",
            # OP=NEW replaces the previous step's loads instead of adding to them.
            "*CLOAD" + (", OP=NEW" if i else ""),
        ]
        for node, vec in loads.items():
            for dof, val in enumerate(vec, start=1):
                if abs(val) > 0.0:
                    parts.append(f"{node}, {dof}, {val:.10e}")
        parts += [
            "*NODE PRINT, NSET=NTIP",
            "U",
            "*NODE PRINT, NSET=NFIXED, TOTALS=ONLY",
            "RF",
            "*NODE FILE",
            "U",
            "*EL FILE",
            "S",
            "*END STEP",
        ]
    deck = job_dir / f"{job}.inp"
    deck.write_text("\n".join(parts) + "\n")
    return deck


def _linear_loads(info: dict, lc: str) -> dict[int, list[float]]:
    """Consistent loads for 3-node triangles: one third of each triangle's force per corner."""
    tris = np.asarray(info["load_tris"], dtype=np.int64)
    xyz = {int(k): np.array(v) for k, v in info["load_xyz"].items()}
    total = np.array(load_vectors(CASES[info["case"]])[lc])
    areas = np.array([0.5 * np.linalg.norm(np.cross(xyz[b] - xyz[a], xyz[c] - xyz[a]))
                      for a, b, c in tris[:, :3]])
    traction = total / areas.sum()
    out: dict[int, np.ndarray] = {}
    for tri, area in zip(tris, areas):
        for n in tri[:3]:
            out.setdefault(int(n), np.zeros(3))
            out[int(n)] += traction * area / 3.0
    return {n: v.tolist() for n, v in sorted(out.items())}


def run(deck: Path, threads: int | None = None) -> float:
    """Run ccx on a deck; return wall time in seconds. Raises if ccx fails."""
    exe = shutil.which("ccx") or shutil.which("ccx_2.21") or shutil.which("ccx_2.22")
    if exe is None:
        raise FileNotFoundError("ccx not found on PATH")
    env = dict(os.environ)
    n = str(threads or os.cpu_count() or 1)
    env.update(OMP_NUM_THREADS=n, CCX_NPROC_EQUATION_SOLVER=n, CCX_NPROC_STIFFNESS=n)
    t0 = time.perf_counter()
    proc = subprocess.run([exe, "-i", deck.stem], cwd=deck.parent, env=env,
                          capture_output=True, text=True)
    dt = time.perf_counter() - t0
    (deck.parent / f"{deck.stem}.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0 or "*ERROR" in proc.stdout:
        raise RuntimeError(f"ccx failed on {deck.name}:\n{proc.stdout[-3000:]}\n{proc.stderr[-2000:]}")
    return dt


# ---------------------------------------------------------------- results

_FLOAT = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def parse_dat(path: Path) -> list[dict]:
    """Read the *NODE PRINT tables: per step, the tip displacement and the
    total reaction force on the fixed face."""
    steps: list[dict] = []
    mode = None
    for line in path.read_text().splitlines():
        low = line.strip().lower()
        if low.startswith("displacements (vx,vy,vz) for set ntip"):
            steps.append({})
            mode = "u"
            continue
        if low.startswith("total force (fx,fy,fz) for set nfixed"):
            mode = "rf"
            continue
        if low.startswith("forces (fx,fy,fz) for set nfixed"):
            mode = None
            continue
        nums = re.findall(_FLOAT, line)
        if mode == "u" and len(nums) == 4:
            steps[-1]["tip_u"] = [float(x) for x in nums[1:]]
            mode = None
        elif mode == "rf" and len(nums) == 3:
            steps[-1]["reaction"] = [float(x) for x in nums]
            mode = None
    return steps


def parse_frd(path: Path) -> list[dict]:
    """Read nodal result blocks from a CalculiX .frd file (ASCII).

    Each block starts with a '-4  NAME' line and holds '-1' records with a node
    number (10 characters) and values (12 characters each, no separators:
    '-1.23456E-01-2.34567E-02' is two numbers). Returns one dict per block:
    {"name": "DISP" | "STRESS", "nodes": int array, "values": float array}.
    """
    blocks: list[dict] = []
    cur = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(" -4"):
                cur = {"name": line[5:13].strip(), "nodes": [], "values": []}
                blocks.append(cur)
            elif line.startswith(" -1") and cur is not None:
                cur["nodes"].append(int(line[3:13]))
                body = line[13:].rstrip("\n")
                cur["values"].append([float(body[i : i + 12]) for i in range(0, len(body), 12) if body[i : i + 12].strip()])
            elif line.startswith(" -3"):
                cur = None
    for b in blocks:
        b["nodes"] = np.asarray(b["nodes"], dtype=np.int64)
        b["values"] = np.asarray(b["values"], dtype=float)
    return blocks


def summarise(job_dir: Path, job: str, info: dict, seconds: float, element: str) -> dict:
    from . import post

    steps = parse_dat(job_dir / f"{job}.dat")
    blocks = parse_frd(job_dir / f"{job}.frd")
    disp = [b for b in blocks if b["name"] == "DISP"]
    stress = [b for b in blocks if b["name"] == "STRESS"]
    out = {"solver": "calculix", "element": element, "case": info["case"], "level": info["level"],
           "n_nodes": info["n_nodes"], "n_tets": info["n_tets"], "seconds": round(seconds, 2),
           "load_cases": {}}
    for i, lc in enumerate(load_cases(info["case"])):
        s = stress[i]
        nodal = {"nodes": s["nodes"], "stress": s["values"][:, :6]}
        out["load_cases"][lc] = {
            "tip_u": steps[i]["tip_u"],
            "reaction": steps[i]["reaction"],
            "max_abs_u": float(np.abs(disp[i]["values"][:, :3]).max()),
            "max_von_mises": float(post.von_mises(nodal["stress"]).max()),
            **post.feature_stresses(info, nodal),
        }
    return out


def solve(info: dict, mesh_dir: Path, job_root: Path, element: str = "C3D10") -> dict:
    deck = write_deck(info, mesh_dir, job_root, element)
    seconds = run(deck)
    return summarise(job_root, deck.stem, info, seconds, element)
