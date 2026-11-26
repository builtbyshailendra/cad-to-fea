"""The known answers the solvers are checked against.

Cantilever tip deflection under an end load P (length L, section b x h):
  Euler-Bernoulli (bending only):   d_EB = P L^3 / (3 E I)
  Timoshenko (bending + shear):     d_T  = d_EB + P L / (k G A)
with I = b h^3 / 12, A = b h, G = E / (2 (1 + nu)) and k the shear correction
factor. For a rectangle, Cowper (1966) gives k = 10 (1 + nu) / (12 + 11 nu)
(0.850 at nu = 0.3); the older constant k = 5/6 (0.833) is also common.

Hole in a finite-width strip under tension (width W, hole diameter d):
  Heywood's formula, as given in Peterson's Stress Concentration Factors:
  Kt_net = 2 + (1 - d/W)^3, based on the net-section stress P / ((W - d) t).
It is a 2D (plane stress) result.
"""
from __future__ import annotations

from .params import STEEL, Beam, Bracket, Material


def cowper_k(nu: float) -> float:
    return 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu)


def euler_bernoulli_tip(p: Beam, mat: Material = STEEL) -> float:
    I = p.b * p.h**3 / 12.0
    return p.P * p.L**3 / (3.0 * mat.E * I)


def shear_tip(p: Beam, mat: Material = STEEL, k: float | None = None) -> float:
    k = cowper_k(mat.nu) if k is None else k
    return p.P * p.L / (k * mat.G * p.b * p.h)


def timoshenko_tip(p: Beam, mat: Material = STEEL, k: float | None = None) -> float:
    return euler_bernoulli_tip(p, mat) + shear_tip(p, mat, k)


def heywood_kt(p: Bracket) -> float:
    return 2.0 + (1.0 - p.d / p.W) ** 3


def net_stress(p: Bracket) -> float:
    """Nominal stress at the hole's net section under the pull load (MPa)."""
    return p.P_pull / ((p.W - p.d) * p.t)


def shelf_root_stress(p: Bracket) -> float:
    """Rough beam-theory bending stress in the arm where it meets the fillet (MPa).

    Treats the arm as a cantilever from the end of the fillet (x = t + r), which
    ignores the fillet's stiffening; a sanity check, not a target.
    """
    M = p.P_shelf * (p.L - (p.t + p.r))
    return 6.0 * M / (p.W * p.t**2)


def summary(mat: Material = STEEL) -> dict:
    from .params import BEAM, BRACKET, STOCKY

    out = {}
    for b in (BEAM, STOCKY):
        eb, ts = euler_bernoulli_tip(b, mat), timoshenko_tip(b, mat)
        out[b.name] = {
            "L_over_h": b.L / b.h,
            "euler_bernoulli": eb,
            "timoshenko_cowper": ts,
            "timoshenko_5_6": timoshenko_tip(b, mat, 5.0 / 6.0),
            "shear_share": (ts - eb) / ts,
        }
    out["bracket"] = {
        "heywood_kt": heywood_kt(BRACKET),
        "net_stress": net_stress(BRACKET),
        "shelf_root_stress": shelf_root_stress(BRACKET),
    }
    return out
