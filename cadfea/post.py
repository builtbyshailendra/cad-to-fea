"""Stress numbers at the features, computed the same way for both solvers.

Input: nodal stresses (averaged over the elements around each node) for a set
of Gmsh node tags, columns xx, yy, zz, xy, yz, zx in MPa. Both solvers' results
are turned into this shape first, so the comparison is like for like.
"""
from __future__ import annotations

import numpy as np

from .params import CASES, Bracket


def tensor(s: np.ndarray) -> np.ndarray:
    xx, yy, zz, xy, yz, zx = s.T
    return np.stack([np.stack([xx, xy, zx], -1), np.stack([xy, yy, yz], -1),
                     np.stack([zx, yz, zz], -1)], -2)


def von_mises(s: np.ndarray) -> np.ndarray:
    xx, yy, zz, xy, yz, zx = s.T
    return np.sqrt(0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2)
                   + 3.0 * (xy**2 + yz**2 + zx**2))


def max_principal(s: np.ndarray) -> np.ndarray:
    return np.linalg.eigvalsh(tensor(s))[:, -1]


def feature_stresses(info: dict, nodal: dict) -> dict:
    """Peak stresses on the hole and fillet faces of the bracket.

    nodal: {"nodes": Gmsh node tags, "stress": (n, 6) array}.
    For the hole, sigma_xx is the hoop stress at the two points where it peaks
    under the pull load (y = +-d/2), so its maximum over the hole is what Kt uses.
    The through-thickness profile shows how that peak varies from the faces
    (z = 0 and z = t) to the mid-plane.
    """
    p = CASES[info["case"]]
    if not isinstance(p, Bracket):
        return {}
    row = {int(n): i for i, n in enumerate(nodal["nodes"])}
    out = {}
    for group in ("HOLE", "FILLET"):
        tags = [n for n in info["nsets"][group] if n in row]
        s = nodal["stress"][[row[n] for n in tags]]
        xyz = np.array([info["feature_xyz"][n] for n in tags])
        key = group.lower()
        vm = von_mises(s)
        s1 = max_principal(s)
        out[f"{key}_von_mises_max"] = float(vm.max())
        out[f"{key}_s1_max"] = float(s1.max())
        out[f"{key}_s1_max_at"] = xyz[int(np.argmax(s1))].round(3).tolist()
        if group == "HOLE":
            sxx = s[:, 0]
            i = int(np.argmax(sxx))
            out["hole_sxx_max"] = float(sxx[i])
            out["hole_sxx_max_at"] = xyz[i].round(3).tolist()
            on_face = (xyz[:, 2] < 1e-6) | (xyz[:, 2] > p.t - 1e-6)
            out["hole_sxx_face_max"] = float(sxx[on_face].max())
            mid = np.abs(xyz[:, 2] - p.t / 2) < p.t / 8
            out["hole_sxx_mid_max"] = float(sxx[mid].max())
            edges = np.linspace(0.0, p.t, 9)
            prof = []
            for lo, hi in zip(edges[:-1], edges[1:]):
                sel = (xyz[:, 2] >= lo - 1e-9) & (xyz[:, 2] <= hi + 1e-9)
                prof.append([round((lo + hi) / 2, 3), float(sxx[sel].max())])
            out["hole_sxx_profile"] = prof
    return out
