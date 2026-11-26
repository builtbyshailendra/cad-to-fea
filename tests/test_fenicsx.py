"""Step 4: FEniCSx on the same meshes (runs in the DOLFINx container in CI)."""
import numpy as np
import pytest

from cadfea import mesh, pipeline, theory
from cadfea.params import BEAM, CASES, LEVELS, STOCKY, levels

from test_ccx import TIMOSHENKO_TOL, check_kt

pytestmark = pytest.mark.needs("dolfinx")


@pytest.mark.parametrize("case", list(CASES))
def test_reactions_balance_the_load(case):
    for level in levels(case):
        r = pipeline.ensure_fenicsx(case, level)
        for lc, total in mesh.load_vectors(CASES[case]).items():
            np.testing.assert_allclose(r["load_cases"][lc]["reaction"], -np.array(total),
                                       rtol=1e-6, atol=1e-6 * np.abs(total).max())


@pytest.mark.parametrize("case", list(CASES))
def test_load_face_area(case):
    """FEniCSx integrates the traction over the tip face it sees: W t or b h exactly."""
    p = CASES[case]
    area = p.b * p.h if hasattr(p, "b") else p.W * p.t
    assert abs(pipeline.ensure_fenicsx(case, levels(case)[0])["load_area"] / area - 1) < 1e-12


@pytest.mark.parametrize("p", [BEAM, STOCKY], ids=["beam", "stocky"])
def test_cantilever_matches_timoshenko(p):
    """Same tolerances, and the same reason, as tests/test_ccx.py."""
    exact = theory.timoshenko_tip(p)
    d = -pipeline.ensure_fenicsx(p.name, "fine")["load_cases"]["tip"]["tip_u"][2]
    assert abs(d / exact - 1) < TIMOSHENKO_TOL[p.name], (d, exact)
    assert d < exact


@pytest.mark.parametrize("case", ["beam", "stocky"])
def test_p2_converges_and_p1_locks(case):
    q = [-pipeline.ensure_fenicsx(case, lv, 2)["load_cases"]["tip"]["tip_u"][2] for lv in LEVELS]
    lin = [-pipeline.ensure_fenicsx(case, lv, 1)["load_cases"]["tip"]["tip_u"][2] for lv in LEVELS]
    assert abs(q[2] - q[1]) <= abs(q[1] - q[0]) + 1e-12
    assert all(a < b for a, b in zip(lin, q))
    assert lin[0] / q[2] < 0.9
    assert lin[0] < lin[1] < lin[2]


def test_hole_kt_against_heywood():
    check_kt({(c, lv): pipeline.ensure_fenicsx(c, lv) for c, lv in [("thin", "fine"), ("bracket", "fine"), ("bracket", "medium")]})


def test_cg_with_multigrid_matches_lu():
    """Iterative CG + algebraic multigrid (with the rigid-body near-nullspace) gives the LU answer."""
    from cadfea import fenics, paths

    info = pipeline.ensure_mesh("beam", "medium")
    lu = pipeline.ensure_fenicsx("beam", "medium")
    cg = fenics.solve(paths.MESH_DIR / "beam_medium.msh", info, solver="cg")
    a, b = lu["load_cases"]["tip"]["tip_u"][2], cg["load_cases"]["tip"]["tip_u"][2]
    assert abs(b / a - 1) < 1e-8, (a, b)
    assert lu["load_cases"]["tip"]["ksp_iterations"] <= 1  # a direct solve
    assert 1 < cg["load_cases"]["tip"]["ksp_iterations"] < 200
