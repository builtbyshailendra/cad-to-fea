"""Step 2: the Gmsh meshes are faithful to the CAD and ready for both solvers."""
import numpy as np
import pytest

from cadfea import mesh, paths, pipeline, props
from cadfea.params import CASES, LEVELS, BRACKET

pytestmark = pytest.mark.needs("cadquery", "gmsh")

ALL = [(c, lv) for c in CASES for lv in LEVELS]


@pytest.mark.parametrize("case,level", ALL)
def test_named_faces_and_counts(case, level):
    info = pipeline.ensure_mesh(case, level)
    expected = {"FIXED", "LOAD", "HOLE", "FILLET"} if case == "bracket" else {"FIXED", "LOAD"}
    assert set(info["nsets"]) == expected
    assert info["n_dofs"] == 3 * info["n_nodes"]
    assert info["min_jacobian_det"] > 0  # no inverted curved element
    assert info["tip_node_offset"] == 0.0  # the embedded point is a mesh node


@pytest.mark.parametrize("case,level", ALL)
def test_mesh_volume_matches_cad(case, level):
    """10-node tets follow the curved hole and fillet; straight tets would not."""
    info = pipeline.ensure_mesh(case, level)
    exact = props.analytic(CASES[case])["volume"]
    assert abs(info["volume_curved"] / exact - 1) < 1e-6
    if case == "bracket":
        # Straight-sided tets cut the hole's and fillet's arcs short: 4e-5 to 4e-4.
        assert abs(info["volume_straight"] / exact - 1) > 1e-5


@pytest.mark.parametrize("case,level", ALL)
def test_consistent_tip_load(case, level):
    info = pipeline.ensure_mesh(case, level)
    xyz, elems = mesh.read_inp(paths.MESH_DIR / f"{case}_{level}_mesh.inp")
    tris = np.asarray(info["load_tris"])
    corners, mids = set(tris[:, :3].ravel()), set(tris[:, 3:].ravel())
    for lc, total in mesh.load_vectors(CASES[case]).items():
        loads = info["loads"][lc]
        got = np.sum(list(loads.values()), axis=0)
        np.testing.assert_allclose(got, total, rtol=1e-12, atol=1e-9)
        # Quadratic triangles: the corners carry nothing, the mid-edge nodes everything.
        assert set(loads) == mids and not (set(loads) & corners)
        assert all(abs(xyz[n][0] - CASES[case].L) < 1e-9 for n in loads)  # all on the tip face


@pytest.mark.parametrize("case,level", [("bracket", "coarse"), ("beam", "medium")])
def test_inp_is_abaqus_ordered(case, level):
    """Gmsh writes C3D10 in Abaqus node order, which differs from its own.

    Each mid-edge node must sit (near) the middle of the edge Abaqus says it is
    on. Read with Gmsh's own edge order instead, the last two nodes land on
    the wrong edges, and the check fails by a mile.
    """
    info = pipeline.ensure_mesh(case, level)
    xyz, elems = mesh.read_inp(paths.MESH_DIR / f"{case}_{level}_mesh.inp")
    assert list(elems) == ["C3D10"]
    tets = elems["C3D10"][:, 1:]
    assert len(tets) == info["n_tets"]
    assert not np.isnan(xyz[1:]).any() and len(xyz) - 1 == info["n_nodes"]
    abaqus = mesh.midnode_offsets(xyz, tets, mesh.ABAQUS_TET10_EDGES)
    assert np.median(abaqus) < 1e-12  # straight edges: exactly in the middle
    assert abaqus.max() < 0.1  # curved edges near the hole: a little off
    gmsh_order = mesh.midnode_offsets(xyz, tets, mesh.GMSH_TET10_EDGES)
    assert np.median(gmsh_order) > 0.3


def test_msh_uses_gmsh_order():
    import gmsh

    pipeline.ensure_mesh("bracket", "coarse")
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.merge(str(paths.MESH_DIR / "bracket_coarse.msh"))
        tags, coords, _ = gmsh.model.mesh.getNodes()
        xyz = np.zeros((int(tags.max()) + 1, 3))
        xyz[tags.astype(np.int64)] = coords.reshape(-1, 3)
        _, _, nodes = gmsh.model.mesh.getElements(3)
        names = [gmsh.model.getPhysicalName(d, t) for d, t in gmsh.model.getPhysicalGroups()]
    finally:
        gmsh.finalize()
    tets = np.asarray(nodes[0], dtype=np.int64).reshape(-1, 10)
    assert np.median(mesh.midnode_offsets(xyz, tets, mesh.GMSH_TET10_EDGES)) < 1e-12
    assert sorted(names) == ["FILLET", "FIXED", "HOLE", "LOAD", "SOLID"]


def test_size_field_refines_the_hole():
    """Surface triangles on the hole are much smaller than on the far faces."""
    info = pipeline.ensure_mesh("bracket", "medium")
    assert len(info["nsets"]["HOLE"]) > 1000
    # Hole surface 251 mm2 vs root face 2400 mm2, yet more nodes on the hole.
    assert len(info["nsets"]["HOLE"]) > 2 * len(info["nsets"]["FIXED"])


def test_meshing_is_deterministic(tmp_path):
    step = pipeline.ensure_step("bracket", "cadquery")
    a = mesh.build_mesh("bracket", "coarse", step, tmp_path / "a")
    b = mesh.build_mesh("bracket", "coarse", step, tmp_path / "b")
    assert (a["n_nodes"], a["n_tets"]) == (b["n_nodes"], b["n_tets"])
    xa, _ = mesh.read_inp(tmp_path / "a" / "bracket_coarse_mesh.inp")
    xb, _ = mesh.read_inp(tmp_path / "b" / "bracket_coarse_mesh.inp")
    np.testing.assert_array_equal(xa[1:], xb[1:])


def test_build123d_step_meshes_the_same_part(tmp_path):
    """Meshing the other tool's STEP file: same volume, same named faces."""
    step = pipeline.ensure_step("bracket", "build123d")
    info = mesh.build_mesh("bracket", "coarse", step, tmp_path)
    exact = props.analytic(BRACKET)["volume"]
    assert abs(info["volume_curved"] / exact - 1) < 1e-6
    assert set(info["nsets"]) == {"FIXED", "LOAD", "HOLE", "FILLET"}
