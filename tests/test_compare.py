"""Step 5: CalculiX against FEniCSx on the same meshes.

Runs in the FEniCSx job, which downloads the CalculiX results from the first job.
"""
import numpy as np
import pytest

from cadfea import pipeline
from cadfea.params import CASES, LEVELS, levels

pytestmark = pytest.mark.needs("dolfinx")

PAIRS = [(c, lv) for c in CASES for lv in levels(c)]


def _both(case, level):
    a = pipeline.load_result("calculix", case, level)
    assert a is not None, "CalculiX results missing (they come from the first CI job)"
    b = pipeline.ensure_fenicsx(case, level)
    assert a["mesh_id"] == b["mesh_id"], "not the same mesh"
    return a, b


# Straight-sided 10-node tets: C3D10 and P2 are the same element, the loads are
# the same consistent vector, so the answers match to the 7 digits CalculiX
# prints (measured 2e-8 to 8e-8). Curved tets (hole, fillet) are integrated
# differently: CalculiX uses 4 points, FFCx picks a higher-degree rule. The
# gap shrinks as the mesh is refined (bracket: 1.8e-5 coarse, 1.4e-7 fine).
SAME_MESH_TOL = {"beam": 1e-6, "stocky": 1e-6, "bracket": 1e-4, "thin": 1e-4}


@pytest.mark.parametrize("case,level", PAIRS)
def test_same_mesh_same_tip_displacement(case, level):
    a, b = _both(case, level)
    for lc in a["load_cases"]:
        ua, ub = np.array(a["load_cases"][lc]["tip_u"]), np.array(b["load_cases"][lc]["tip_u"])
        assert np.abs(ua - ub).max() / np.abs(ub).max() < SAME_MESH_TOL[case], (lc, ua, ub)


@pytest.mark.parametrize("case,level", [(c, lv) for c, lv in PAIRS if c in ("bracket", "thin")])
def test_same_mesh_same_hole_stress(case, level):
    """Peak nodal stress at the hole, CalculiX vs FEniCSx.

    Peaks sit on curved elements, where the two codes differ most; the gap
    shrinks with the elements: 1.4% coarse, 0.8% medium, 0.03% fine (measured).
    """
    a, b = _both(case, level)
    tol = {"coarse": 0.02, "medium": 0.01, "fine": 0.005}[level]
    for lc in ("shelf", "pull"):
        sa, sb = a["load_cases"][lc]["hole_s1_max"], b["load_cases"][lc]["hole_s1_max"]
        assert abs(sa / sb - 1) < tol, (lc, sa, sb)


def test_linear_tets_agree_too():
    """C3D4 and P1 are also the same element: same answer to the printed digits."""
    for case in ("beam", "stocky"):
        for level in LEVELS:
            a = pipeline.load_result("calculix", case, level, "C3D4")
            b = pipeline.ensure_fenicsx(case, level, 1)
            assert a is not None
            ua, ub = a["load_cases"]["tip"]["tip_u"][2], b["load_cases"]["tip"]["tip_u"][2]
            assert abs(ua / ub - 1) < 1e-6
            assert a["n_dofs"] == b["n_dofs"]
