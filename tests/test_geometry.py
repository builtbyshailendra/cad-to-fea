"""Step 1: the same parameters in CadQuery and build123d give the same solid."""
import dataclasses

import pytest

from cadfea import props
from cadfea.params import BRACKET, CASES

pytestmark = pytest.mark.needs("cadquery", "build123d")

# Both tools sit on the same kernel (OpenCASCADE 7.9) and the shapes are planes
# and cylinders, so agreement is at rounding level (about 1e-16 measured). The
# tolerance is far tighter than the 1e-6 the plan asked for, and still 1000x
# looser than what was measured.
TOL = 1e-12


@pytest.fixture(scope="module")
def models():
    from cadfea import geom_b123d, geom_cq

    return {name: (geom_cq.build(p), geom_b123d.build(p)) for name, p in CASES.items()}


@pytest.mark.parametrize("case", list(CASES))
def test_cadquery_and_build123d_agree(models, case):
    from cadfea import geom_b123d, geom_cq

    a = geom_cq.properties(models[case][0])
    b = geom_b123d.properties(models[case][1])
    diff = props.rel_diff(a, b)
    assert set(diff) == {"volume", "area", "com", "bbox", "inertia"}
    assert max(diff.values()) < TOL, diff
    assert (a["n_faces"], a["n_edges"]) == (b["n_faces"], b["n_edges"])


@pytest.mark.parametrize("case", list(CASES))
def test_matches_closed_form(models, case):
    from cadfea import geom_b123d, geom_cq

    exact = props.analytic(CASES[case])
    for got in (geom_cq.properties(models[case][0]), geom_b123d.properties(models[case][1])):
        diff = props.rel_diff(got, exact)
        assert max(diff.values()) < TOL, diff


@pytest.mark.parametrize("case", list(CASES))
def test_step_round_trip(models, case, tmp_path):
    """Write both STEP files, read them back with plain OpenCASCADE, compare."""
    from cadfea import geom_b123d, geom_cq

    sa = geom_cq.export_step(models[case][0], tmp_path / "cq.step")
    sb = geom_b123d.export_step(models[case][1], tmp_path / "bd.step")
    ra, rb = props.step_properties(sa), props.step_properties(sb)
    assert ra["n_solids"] == rb["n_solids"] == 1
    assert ra["n_faces"] == rb["n_faces"]
    assert max(props.rel_diff(ra, rb).values()) < TOL
    assert max(props.rel_diff(ra, geom_cq.properties(models[case][0])).values()) < 1e-12
    # STEP stores millimetres here; a unit mix-up would be off by 1000.
    assert "SI_UNIT(.MILLI.,.METRE.)" in sa.read_text()


def test_check_has_teeth():
    """A 1 micrometre change in the fillet radius must fail the comparison."""
    from cadfea import geom_b123d, geom_cq

    a = geom_cq.properties(geom_cq.build(BRACKET))
    b = geom_b123d.properties(geom_b123d.build(dataclasses.replace(BRACKET, r=BRACKET.r + 0.001)))
    diff = props.rel_diff(a, b)
    assert diff["volume"] > 1e-6  # 1.5e-6 measured: a million times the tolerance
    assert max(diff.values()) > TOL


def test_bracket_faces_are_what_the_mesher_expects(models):
    """10 faces: 8 planes, the fillet cylinder and the hole cylinder (one face with a seam)."""
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    for model in models["bracket"]:
        shape = model.val().wrapped if hasattr(model, "val") else model.wrapped
        kinds = []
        exp = TopExp_Explorer(shape, TopAbs_FACE)
        while exp.More():
            kinds.append(BRepAdaptor_Surface(TopoDS.Face_s(exp.Current())).GetType())
            exp.Next()
        assert kinds.count(GeomAbs_Plane) == 8
        assert kinds.count(GeomAbs_Cylinder) == 2
