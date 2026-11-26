"""Step 4: FEniCSx on the same meshes (runs in the DOLFINx container in CI)."""
import numpy as np
import pytest

from cadfea import mesh, pipeline, theory
from cadfea.params import BEAM, BRACKET, CASES, LEVELS, STOCKY

pytestmark = pytest.mark.needs("dolfinx")


@pytest.mark.parametrize("case", list(CASES))
def test_reactions_balance_the_load(case):
    for level in LEVELS:
        r = pipeline.ensure_fenicsx(case, level)
        for lc, total in mesh.load_vectors(CASES[case]).items():
            np.testing.assert_allclose(r["load_cases"][lc]["reaction"], -np.array(total),
                                       rtol=1e-6, atol=1e-6 * np.abs(total).max())


@pytest.mark.parametrize("case", list(CASES))
def test_load_face_area(case):
    """FEniCSx integrates the traction over the tip face it sees: W t or b h exactly."""
    p = CASES[case]
    area = p.W * p.t if case == "bracket" else p.b * p.h
    assert abs(pipeline.ensure_fenicsx(case, "coarse")["load_area"] / area - 1) < 1e-12


@pytest.mark.parametrize("p", [BEAM, STOCKY], ids=["beam", "stocky"])
def test_cantilever_matches_timoshenko(p):
    exact = theory.timoshenko_tip(p)
    d = -pipeline.ensure_fenicsx(p.name, "fine")["load_cases"]["tip"]["tip_u"][2]
    assert abs(d / exact - 1) < 0.01, (d, exact)


@pytest.mark.parametrize("case", ["beam", "stocky"])
def test_p2_converges_and_p1_locks(case):
    q = [-pipeline.ensure_fenicsx(case, lv, 2)["load_cases"]["tip"]["tip_u"][2] for lv in LEVELS]
    lin = [-pipeline.ensure_fenicsx(case, lv, 1)["load_cases"]["tip"]["tip_u"][2] for lv in LEVELS]
    assert abs(q[2] - q[1]) <= abs(q[1] - q[0]) + 1e-12
    assert all(a < b for a, b in zip(lin, q))
    assert lin[0] / q[2] < 0.9
    assert lin[0] < lin[1] < lin[2]


def test_bracket_hole_kt():
    r = pipeline.ensure_fenicsx("bracket", "fine")["load_cases"]["pull"]
    kt = r["hole_sxx_max"] / theory.net_stress(BRACKET)
    assert abs(kt / theory.heywood_kt(BRACKET) - 1) < 0.05, kt
    assert r["hole_sxx_mid_max"] > r["hole_sxx_face_max"]
