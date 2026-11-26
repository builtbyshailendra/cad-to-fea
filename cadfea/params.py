"""Every number that defines the models, in one place.

Units throughout: millimetres, newtons, megapascals (N/mm^2).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class Material:
    name: str = "STEEL"
    E: float = 210_000.0  # MPa
    nu: float = 0.3

    @property
    def G(self) -> float:
        return self.E / (2.0 * (1.0 + self.nu))


@dataclass(frozen=True)
class Beam:
    """Plain rectangular cantilever, root face at x=0, tip face at x=L.

    Cross-section b (along y) by h (along z), centred on the x axis.
    The tip load P acts in -z, spread uniformly over the tip face.
    """

    name: str = "beam"
    L: float = 200.0
    b: float = 20.0
    h: float = 20.0
    P: float = 1000.0

    @property
    def tip_point(self) -> tuple[float, float, float]:
        return (self.L, 0.0, 0.0)


@dataclass(frozen=True)
class Bracket:
    """L-shaped angle bracket: an upright bolted to a wall and an arm.

    The L profile lies in the x-z plane and is extruded W along y (centred).
    Upright: x in [0, t], z in [0, H]; its back face x=0 is fixed.
    Arm:     x in [0, L], z in [0, t]; its tip face x=L carries the load.
    Inner corner at (t, t) has a fillet of radius r.
    A through hole of diameter d (along z) sits in the arm at x = xh, y = 0.

    Two load cases on the tip face:
      shelf: P_shelf in -z (the arm bends, like a shelf bracket)
      pull:  P_pull in +x  (the arm is in tension, for the hole's Kt)
    """

    name: str = "bracket"
    L: float = 160.0
    W: float = 40.0
    t: float = 8.0
    H: float = 60.0
    r: float = 6.0
    d: float = 10.0
    xh: float = 100.0
    P_shelf: float = 200.0
    P_pull: float = 10_000.0

    @property
    def tip_point(self) -> tuple[float, float, float]:
        return (self.L, 0.0, self.t / 2.0)


STEEL = Material()
# Same steel with Poisson's ratio set to zero: an experiment, not a material.
# With nu = 0 the clamped root cannot restrain any sideways contraction.
NU0 = Material(name="NU0", nu=0.0)

# Two plain cantilevers: slender (L/h = 10), where Euler-Bernoulli and
# Timoshenko differ by 0.8%, and stocky (L/h = 4), where they differ by 4.8%.
BEAM = Beam()
STOCKY = Beam(name="stocky", L=80.0)
BRACKET = Bracket()
# The same bracket in 2 mm plate (t/d = 0.2 instead of 0.8): close to the
# plane-stress conditions Heywood's Kt assumes. Used for the Kt check only.
THIN = Bracket(name="thin", t=2.0)

CASES: dict[str, Beam | Bracket] = {c.name: c for c in (BEAM, STOCKY, BRACKET, THIN)}


@dataclass(frozen=True)
class MeshLevel:
    name: str
    size: float  # target element size away from features, mm
    feature: float  # target size on the hole surface (bracket only), mm
    fillet: float | None = None  # target size on the fillet surface; None: same as the hole


# Three levels per case for the convergence study. Element sizes in mm.
MESH_LEVELS: dict[str, dict[str, MeshLevel]] = {
    "beam": {
        "coarse": MeshLevel("coarse", 10.0, 10.0),
        "medium": MeshLevel("medium", 5.0, 5.0),
        "fine": MeshLevel("fine", 2.5, 2.5),
    },
    "stocky": {
        "coarse": MeshLevel("coarse", 10.0, 10.0),
        "medium": MeshLevel("medium", 5.0, 5.0),
        "fine": MeshLevel("fine", 2.5, 2.5),
    },
    "bracket": {
        "coarse": MeshLevel("coarse", 8.0, 2.0),
        "medium": MeshLevel("medium", 5.0, 1.0),
        "fine": MeshLevel("fine", 3.5, 0.6),
    },
    "thin": {
        "fine": MeshLevel("fine", 3.0, 0.4, fillet=1.0),
    },
}
LEVELS = ("coarse", "medium", "fine")


def levels(case: str) -> list[str]:
    return list(MESH_LEVELS[case])


def as_dict(obj) -> dict:
    return asdict(obj)
