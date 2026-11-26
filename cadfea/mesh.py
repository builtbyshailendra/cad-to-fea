"""Gmsh: import the STEP file, name the faces, refine where stress peaks, mesh.

Writes, per case and level:
  <stem>.msh       Gmsh 4.1 binary mesh with physical groups (read by FEniCSx)
  <stem>_mesh.inp  Abaqus-format mesh written by Gmsh (read by CalculiX)
  <stem>.json      node sets, the consistent tip-load vector, and mesh statistics
                   (used to write the CalculiX deck and to check the mesh)

Physical groups (the names solvers see):
  3D  SOLID   the part
  2D  FIXED   the root face (x = 0)
  2D  LOAD    the tip face (x = L)
  2D  HOLE    the hole's cylindrical face (bracket only)
  2D  FILLET  the fillet face (bracket only)
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from .params import CASES, MESH_LEVELS, Beam, Bracket, MeshLevel

GROUP_TAGS = {"SOLID": 1, "FIXED": 11, "LOAD": 12, "HOLE": 13, "FILLET": 14}
TOL = 1e-6


def _face_info(gmsh) -> list[dict]:
    faces = []
    for _, tag in gmsh.model.getEntities(2):
        bb = gmsh.model.getBoundingBox(2, tag)
        faces.append({"tag": tag, "type": gmsh.model.getType(2, tag), "bbox": bb})
    return faces


def classify_faces(gmsh, p: Beam | Bracket) -> dict[str, list[int]]:
    """Find the named faces by where they are, not by their tag numbers.

    Face tags depend on how the CAD tool ordered the B-rep; positions don't.
    """
    groups: dict[str, list[int]] = {"FIXED": [], "LOAD": [], "HOLE": [], "FILLET": []}
    for f in _face_info(gmsh):
        xmin, ymin, zmin, xmax, ymax, zmax = f["bbox"]
        if f["type"] == "Plane" and xmax < TOL:
            groups["FIXED"].append(f["tag"])
        elif f["type"] == "Plane" and xmin > p.L - TOL:
            groups["LOAD"].append(f["tag"])
        elif isinstance(p, Bracket) and f["type"] == "Cylinder":
            near_hole = xmin > p.xh - p.d / 2 - 0.01 and xmax < p.xh + p.d / 2 + 0.01
            near_corner = xmin > p.t - 0.01 and xmax < p.t + p.r + 0.01 and zmax < p.t + p.r + 0.01
            if near_hole:
                groups["HOLE"].append(f["tag"])
            elif near_corner:
                groups["FILLET"].append(f["tag"])
    expected = ("FIXED", "LOAD", "HOLE", "FILLET") if isinstance(p, Bracket) else ("FIXED", "LOAD")
    for name in expected:
        if not groups[name]:
            raise RuntimeError(f"no face found for {name}")
    return {k: v for k, v in groups.items() if v}


def _size_field(gmsh, groups: dict[str, list[int]], level: MeshLevel, p: Beam | Bracket) -> None:
    """Small elements on the hole and fillet, growing to the global size further away.

    One Distance -> Threshold pair per feature, combined with a Min field: the
    element size anywhere is the smallest any feature asks for.
    """
    if not isinstance(p, Bracket):
        return
    f = gmsh.model.mesh.field
    fillet = level.feature if level.fillet is None else level.fillet
    thresholds = []
    for faces, size in ((groups["HOLE"], level.feature), (groups["FILLET"], fillet)):
        dist = f.add("Distance")
        f.setNumbers(dist, "SurfacesList", faces)
        f.setNumber(dist, "Sampling", 60)
        thr = f.add("Threshold")
        f.setNumber(thr, "InField", dist)
        f.setNumber(thr, "SizeMin", size)
        f.setNumber(thr, "SizeMax", level.size)
        # Full refinement within 1.5 mm of the surface, then a linear ramp to
        # the global size 10 mm away (about one hole diameter).
        f.setNumber(thr, "DistMin", 1.5)
        f.setNumber(thr, "DistMax", 10.0)
        thresholds.append(thr)
    both = f.add("Min")
    f.setNumbers(both, "FieldsList", thresholds)
    f.setAsBackgroundMesh(both)


def _nodes_of_group(gmsh, dim: int, tag: int) -> list[int]:
    nodes, _ = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
    return sorted(int(n) for n in nodes)


def load_triangles(gmsh, face_tags: list[int]) -> np.ndarray:
    """The 6-node triangles (node tags, Gmsh order: 3 corners, then mid-edges) of some faces."""
    tri_nodes = []
    for face in face_tags:
        types, _, node_tags = gmsh.model.mesh.getElements(2, face)
        for etype, nodes in zip(types, node_tags):
            if etype != 9:  # 9 = 6-node triangle
                raise RuntimeError(f"unexpected surface element type {etype}")
            tri_nodes.append(np.asarray(nodes, dtype=np.int64).reshape(-1, 6))
    return np.vstack(tri_nodes)


def consistent_load(tris: np.ndarray, xyz: np.ndarray, total: np.ndarray) -> tuple[dict[int, list[float]], float]:
    """Nodal forces equivalent to a uniform traction on 6-node triangles.

    For a quadratic triangle under uniform traction the consistent nodal loads are
    zero at the three corners and one third of the triangle's force at each
    mid-side node. Spreading the force equally over all nodes would be wrong.
    """
    a, b, c = xyz[tris[:, 0]], xyz[tris[:, 1]], xyz[tris[:, 2]]
    areas = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    traction = total / areas.sum()
    loads: dict[int, np.ndarray] = {}
    for tri, area in zip(tris, areas):
        for n in tri[3:]:
            loads.setdefault(int(n), np.zeros(3))
            loads[int(n)] += traction * area / 3.0
    return {n: v.tolist() for n, v in sorted(loads.items())}, float(areas.sum())


def _straight_volumes(xyz: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Volumes of the tets from their 4 corner nodes, as if every edge were straight."""
    a, b, c, d = (xyz[tets[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def _curved_volume(gmsh) -> tuple[float, float]:
    """Volume of the 10-node tets including their curved edges, and the smallest
    Jacobian determinant at any quadrature point (must be > 0: no inverted elements).

    The map from the reference tet is quadratic, so det J is cubic; Gauss4 is exact.
    """
    local, weights = gmsh.model.mesh.getIntegrationPoints(11, "Gauss4")
    _, dets, _ = gmsh.model.mesh.getJacobians(11, local)
    dets = np.asarray(dets).reshape(-1, len(weights))
    return float((dets * np.asarray(weights)).sum()), float(dets.min())


def build_mesh(case: str, level: str, step_file: Path, out_dir: Path) -> dict:
    import gmsh

    p = CASES[case]
    lv = MESH_LEVELS[case][level]
    stem = f"{case}_{level}"
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)  # deterministic meshes
        gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
        gmsh.model.add(stem)
        gmsh.model.occ.importShapes(str(step_file))
        gmsh.model.occ.synchronize()
        vols = gmsh.model.getEntities(3)
        if len(vols) != 1:
            raise RuntimeError(f"expected one volume, got {vols}")
        vol = vols[0][1]
        groups = classify_faces(gmsh, p)

        # A mesh node exactly where tip deflection is measured: embed a point
        # in the tip face, so both solvers report the same node.
        tx, ty, tz = p.tip_point
        tip = gmsh.model.occ.addPoint(tx, ty, tz)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.embed(0, [tip], 2, groups["LOAD"][0])

        gmsh.model.addPhysicalGroup(3, [vol], GROUP_TAGS["SOLID"], name="SOLID")
        for name, tags in groups.items():
            gmsh.model.addPhysicalGroup(2, tags, GROUP_TAGS[name], name=name)

        gmsh.option.setNumber("Mesh.MeshSizeMax", lv.size)
        gmsh.option.setNumber("Mesh.MeshSizeMin", min(lv.size, lv.feature))
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        _size_field(gmsh, groups, lv, p)
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay for surfaces
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay for the volume
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(2)  # 10-node tets; mid-edge nodes snap to the CAD surfaces
        t_mesh = time.perf_counter() - t0

        # Statistics and checks.
        tags, coords, _ = gmsh.model.mesh.getNodes()
        xyz = np.zeros((int(tags.max()) + 1, 3))
        xyz[tags.astype(np.int64)] = coords.reshape(-1, 3)
        etypes, etags, enodes = gmsh.model.mesh.getElements(3, vol)
        if list(etypes) != [11]:  # 11 = 10-node tetrahedron
            raise RuntimeError(f"unexpected volume element types {etypes}")
        tets = np.asarray(enodes[0], dtype=np.int64).reshape(-1, 10)
        vols_e = _straight_volumes(xyz, tets)
        curved_volume, min_det = _curved_volume(gmsh)
        q = gmsh.model.mesh.getElementQualities(np.asarray(etags[0]), "minSICN")

        tip_node = int(np.argmin(np.linalg.norm(xyz[tags.astype(np.int64)] - np.array(p.tip_point), axis=1)))
        tip_tag = int(tags[tip_node])
        tip_err = float(np.linalg.norm(xyz[tip_tag] - np.array(p.tip_point)))

        nsets = {name: _nodes_of_group(gmsh, 2, GROUP_TAGS[name]) for name in groups}
        tris = load_triangles(gmsh, groups["LOAD"])
        loads = {}
        for lc, vec in load_vectors(p).items():
            loads[lc], load_area = consistent_load(tris, xyz, np.array(vec))
        feature_nodes = sorted({n for g in ("HOLE", "FILLET") for n in nsets.get(g, [])} | {tip_tag})
        gmsh_version = gmsh.__version__

        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(str(out_dir / f"{stem}.msh"))
        # Abaqus format for CalculiX: nodes, C3D10 elements, an element set per
        # physical volume. Surface groups would come out as 2D elements that
        # CalculiX can't use, so only the volume group is kept for this file;
        # the node sets are written into the deck from the JSON instead.
        surface_groups = [(2, GROUP_TAGS[n]) for n in groups]
        gmsh.model.removePhysicalGroups(surface_groups)
        gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
        gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 0)
        gmsh.write(str(out_dir / f"{stem}_mesh.inp"))
    finally:
        gmsh.finalize()

    info = {
        "case": case,
        "level": level,
        "size": lv.size,
        "feature_size": lv.feature,
        "step_file": step_file.name,
        "n_nodes": int(len(tags)),
        "n_tets": int(len(tets)),
        "n_dofs": int(3 * len(tags)),
        "volume_straight": float(vols_e.sum()),
        "volume_curved": curved_volume,
        "min_jacobian_det": min_det,
        "min_tet_volume": float(vols_e.min()),
        "min_quality_sicn": float(np.min(q)),
        "mean_quality_sicn": float(np.mean(q)),
        "tip_node": tip_tag,
        "tip_node_offset": tip_err,
        "load_area": load_area,
        "gmsh": gmsh_version,
        "nsets": nsets,
        "loads": loads,
        "load_tris": tris.tolist(),
        "load_xyz": {int(n): xyz[n].tolist() for n in np.unique(tris)},
        "feature_xyz": {int(n): xyz[n].tolist() for n in feature_nodes},
        "seconds": round(t_mesh, 3),
    }
    (out_dir / f"{stem}.json").write_text(json.dumps(info))
    return info


def load_vectors(p: Beam | Bracket) -> dict[str, list[float]]:
    """Total force on the tip face for each load case (N)."""
    if isinstance(p, Bracket):
        return {"shelf": [0.0, 0.0, -p.P_shelf], "pull": [p.P_pull, 0.0, 0.0]}
    return {"tip": [0.0, 0.0, -p.P]}


def load_info(out_dir: Path, case: str, level: str) -> dict:
    """Read a mesh's JSON sidecar back, with integer node tags as keys."""
    data = json.loads((out_dir / f"{case}_{level}.json").read_text())
    data["loads"] = {lc: {int(k): v for k, v in d.items()} for lc, d in data["loads"].items()}
    data["load_xyz"] = {int(k): v for k, v in data["load_xyz"].items()}
    data["feature_xyz"] = {int(k): v for k, v in data["feature_xyz"].items()}
    return data


def read_inp(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Nodes (indexed by tag) and elements (by type) from an Abaqus-format mesh."""
    nodes: dict[int, list[float]] = {}
    elems: dict[str, list[list[int]]] = {}
    mode, etype = None, None
    with path.open() as fh:
        for line in fh:
            if line.startswith("*"):
                head = line.upper()
                if head.startswith("*NODE"):
                    mode = "node"
                elif head.startswith("*ELEMENT"):
                    mode, etype = "elem", head.split("TYPE=")[1].split(",")[0].strip()
                    elems.setdefault(etype, [])
                else:
                    mode = None
                continue
            if not line.strip():
                continue
            fields = line.replace(",", " ").split()
            if mode == "node":
                nodes[int(fields[0])] = [float(x) for x in fields[1:4]]
            elif mode == "elem":
                elems[etype].append([int(x) for x in fields])
    xyz = np.full((max(nodes) + 1, 3), np.nan)
    for tag, c in nodes.items():
        xyz[tag] = c
    return xyz, {k: np.asarray(v, dtype=np.int64) for k, v in elems.items()}


# Abaqus C3D10 node order: corners 1-4, then mid-edge nodes on edges
# 1-2, 2-3, 3-1, 1-4, 2-4, 3-4 (zero-based below).
ABAQUS_TET10_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
# Gmsh's own 10-node tet order has the last two edges the other way round.
GMSH_TET10_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (2, 3), (1, 3)]


def midnode_offsets(xyz: np.ndarray, tets: np.ndarray, edges) -> np.ndarray:
    """For each element, the largest distance of a mid-edge node from the midpoint
    of the edge it is supposed to sit on, relative to that edge's length."""
    worst = np.zeros(len(tets))
    for k, (i, j) in enumerate(edges):
        a, b, m = xyz[tets[:, i]], xyz[tets[:, j]], xyz[tets[:, 4 + k]]
        off = np.linalg.norm(m - (a + b) / 2, axis=1) / np.linalg.norm(b - a, axis=1)
        worst = np.maximum(worst, off)
    return worst
