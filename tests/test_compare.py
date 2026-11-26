"""Step 5: CalculiX against FEniCSx on the same meshes.

Runs in the FEniCSx job, which downloads the CalculiX results from the first job.
"""
import numpy as np
import pytest

from cadfea import pipeline
from cadfea.params import CASES, LEVELS

pytestmark = pytest.mark.needs("dolfinx")

PAIRS = [(c, lv) for c in CASES for lv in LEVELS]


def _both(case, level):
    a = pipeline.load_result("calculix", case, level)
    assert a is not None, "CalculiX results missing (they come from the first CI job)"
    b = pipeline.ensure_fenicsx(case, level)
    assert a["mesh_id"] == b["mesh_id"], "not the same mesh"
    return a, b


@pytest.mark.parametrize("case,level", PAIRS)
def test_same_mesh_same_tip_displacement(case, level):
    """C3D10 and P2 on the same straight-sided tets are the same discrete problem."""
    a, b = _both(case, level)
    for lc in a["load_cases"]:
        ua, ub = np.array(a["load_cases"][lc]["tip_u"]), np.array(b["load_cases"][lc]["tip_u"])
        assert np.abs(ua - ub).max() / np.abs(ub).max() < 1e-3, (lc, ua, ub)


@pytest.mark.parametrize("level", LEVELS)
def test_same_mesh_same_hole_stress(level):
    a, b = _both("bracket", level)
    for lc in ("shelf", "pull"):
        sa, sb = a["load_cases"][lc]["hole_s1_max"], b["load_cases"][lc]["hole_s1_max"]
        assert abs(sa / sb - 1) < 0.02, (lc, sa, sb)
