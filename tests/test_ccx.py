"""Step 3: CalculiX. The deck and parsers are tested anywhere; the solves need ccx."""
import numpy as np
import pytest

from cadfea import ccx, mesh, paths, pipeline, theory
from cadfea.params import BEAM, BRACKET, CASES, LEVELS, NU0, STOCKY, THIN, levels

FRD_SAMPLE = """\
    1C
    2C                             2                                     1
 -1         1 0.00000E+00 0.00000E+00 0.00000E+00
 -1         2 1.00000E+00 0.00000E+00 0.00000E+00
 -3
  100CL  101 1.00000E+00         2                     0    1           1
 -4  DISP        4    1
 -5  D1          1    2    1    0
 -1         1 1.00000E-03-2.50000E-01 0.00000E+00 0.00000E+00
 -1         2-1.25000E-03-5.00000E-01 3.00000E-05 0.00000E+00
 -3
  100CL  102 1.00000E+00         2                     0    1           1
 -4  STRESS      6    1
 -5  SXX         1    4    1    1
 -1         1 1.00000E+02-2.00000E+01 0.00000E+00 5.00000E+00 0.00000E+00-1.00000E+00
 -1         2-1.00000E+02 2.00000E+01 0.00000E+00-5.00000E+00 0.00000E+00 1.00000E+00
 -3
 9999
"""

DAT_SAMPLE = """
                        S T E P       1


                                INCREMENT     1


 displacements (vx,vy,vz) for set NTIP and time  0.1000000E+01

        17  1.234567E-06  2.000000E-08 -9.876543E-01

 total force (fx,fy,fz) for set NFIXED and time  0.1000000E+01

        1.000000E-09 -2.000000E-10  1.000000E+03
"""


def test_frd_parser_reads_fixed_width_numbers(tmp_path):
    """'-2.50000E-01' right after '1.00000E-03' with no space: fixed columns, not split()."""
    f = tmp_path / "x.frd"
    f.write_text(FRD_SAMPLE)
    blocks = ccx.parse_frd(f)
    assert [b["name"] for b in blocks] == ["DISP", "STRESS"]
    np.testing.assert_array_equal(blocks[0]["nodes"], [1, 2])
    np.testing.assert_allclose(blocks[0]["values"][1], [-1.25e-3, -0.5, 3e-5, 0.0])
    np.testing.assert_allclose(blocks[1]["values"][0], [100, -20, 0, 5, 0, -1])


def test_dat_parser(tmp_path):
    f = tmp_path / "x.dat"
    f.write_text(DAT_SAMPLE + DAT_SAMPLE.replace("S T E P       1", "S T E P       2"))
    steps = ccx.parse_dat(f)
    assert len(steps) == 2
    assert steps[0]["tip_u"] == [1.234567e-06, 2e-08, -0.9876543]
    assert steps[1]["reaction"] == [1e-09, -2e-10, 1000.0]


@pytest.mark.needs("cadquery", "gmsh")
def test_deck_has_everything(tmp_path):
    info = pipeline.ensure_mesh("bracket", "coarse")
    deck, n_active = ccx.write_deck(info, paths.MESH_DIR, tmp_path)
    assert n_active == info["n_nodes"]
    text = deck.read_text()
    for card in ("*INCLUDE, INPUT=bracket_coarse_mesh.inp", "*NSET, NSET=NFIXED", "*NSET, NSET=NTIP",
                 "*ELASTIC", "*SOLID SECTION, ELSET=SOLID", "NFIXED, 1, 3, 0.0", "*CLOAD, OP=NEW"):
        assert card in text
    assert text.count("*STEP") == 2 and text.count("*END STEP") == 2
    # The CLOAD lines of each step add up to that load case's total force.
    for i, (lc, total) in enumerate(mesh.load_vectors(BRACKET).items()):
        block = text.split("*CLOAD")[i + 1].split("*NODE PRINT")[0].splitlines()[1:]
        f = np.zeros(3)
        for line in block:
            node, dof, val = line.split(",")
            f[int(dof) - 1] += float(val)
        np.testing.assert_allclose(f, total, rtol=1e-9, atol=1e-6)


@pytest.mark.needs("cadquery", "gmsh")
def test_c3d4_deck_uses_corner_nodes_only(tmp_path):
    info = pipeline.ensure_mesh("beam", "coarse")
    deck, n_active = ccx.write_deck(info, paths.MESH_DIR, tmp_path, element="C3D4")
    xyz, elems = mesh.read_inp(tmp_path / "beam_coarse_c3d4_mesh.inp")
    corners = set(elems["C3D4"][:, 1:].ravel())
    assert n_active == len(corners) < info["n_nodes"] / 3  # 10-node tets: most nodes are mid-edge
    assert elems["C3D4"].shape == (info["n_tets"], 5)
    text = deck.read_text()
    loaded = {int(line.split(",")[0]) for line in text.split("*CLOAD")[1].split("*NODE PRINT")[0].splitlines()[1:]}
    assert loaded <= corners


# ---------------------------------------------------------------- with ccx

pytestmark_solve = pytest.mark.needs("ccx", "cadquery", "gmsh")


@pytestmark_solve
@pytest.mark.parametrize("case", list(CASES))
def test_reactions_balance_the_load(case):
    for level in levels(case):
        r = pipeline.ensure_ccx(case, level)
        for lc, total in mesh.load_vectors(CASES[case]).items():
            np.testing.assert_allclose(r["load_cases"][lc]["reaction"], -np.array(total),
                                       rtol=1e-6, atol=1e-6 * np.abs(total).max())


# Tolerances against Timoshenko, and why they differ (see README, step 5):
# the fully clamped root face also stops the steel contracting sideways
# (Poisson's ratio), which beam theory's clamp does not. That makes the 3D
# beam stiffer by about 0.7% at L/h = 10 and 1.7% at L/h = 4 (measured). Run
# with nu = 0 and the gap closes to under 0.1%, so the solver and the mesh are
# fine; the difference is the boundary condition.
TIMOSHENKO_TOL = {"beam": 0.01, "stocky": 0.02}


@pytestmark_solve
@pytest.mark.parametrize("p", [BEAM, STOCKY], ids=["beam", "stocky"])
def test_cantilever_matches_timoshenko(p):
    exact = theory.timoshenko_tip(p)
    d = -pipeline.ensure_ccx(p.name, "fine")["load_cases"]["tip"]["tip_u"][2]
    assert abs(d / exact - 1) < TIMOSHENKO_TOL[p.name], (d, exact)
    assert d < exact  # the clamp makes it stiffer, never softer


@pytestmark_solve
@pytest.mark.parametrize("p", [BEAM, STOCKY], ids=["beam", "stocky"])
def test_poisson_zero_closes_the_gap(p):
    """With nu = 0 the 3D cantilever matches Timoshenko to 0.1% (measured 0.004% and 0.06%)."""
    exact = theory.timoshenko_tip(p, NU0)
    d = -pipeline.ensure_ccx(p.name, "fine", "C3D10", NU0)["load_cases"]["tip"]["tip_u"][2]
    assert abs(d / exact - 1) < 0.001, (d, exact)


@pytestmark_solve
def test_stocky_beam_needs_shear_deformation():
    """At L/h = 4 Euler-Bernoulli misses by 3.7-4.8%; the 3D answer sides with Timoshenko."""
    d0 = -pipeline.ensure_ccx("stocky", "fine", "C3D10", NU0)["load_cases"]["tip"]["tip_u"][2]
    eb0, ts0 = theory.euler_bernoulli_tip(STOCKY, NU0), theory.timoshenko_tip(STOCKY, NU0)
    assert abs(d0 - ts0) < 0.05 * abs(d0 - eb0)
    # With nu = 0.3 the clamp takes some of it back, but most of the shear
    # deflection is still there (measured 63%).
    d = -pipeline.ensure_ccx("stocky", "fine")["load_cases"]["tip"]["tip_u"][2]
    eb, ts = theory.euler_bernoulli_tip(STOCKY), theory.timoshenko_tip(STOCKY)
    assert d - eb > 0.5 * (ts - eb)


@pytestmark_solve
@pytest.mark.parametrize("case", ["beam", "stocky"])
def test_c3d10_converges_and_c3d4_locks(case):
    """Quadratic tets converge from the first mesh; linear tets are far too stiff."""
    q = [-pipeline.ensure_ccx(case, lv)["load_cases"]["tip"]["tip_u"][2] for lv in LEVELS]
    lin = [-pipeline.ensure_ccx(case, lv, "C3D4")["load_cases"]["tip"]["tip_u"][2] for lv in LEVELS]
    assert abs(q[2] - q[1]) <= abs(q[1] - q[0]) + 1e-12
    assert all(a < b for a, b in zip(lin, q))  # linear tets always stiffer
    assert lin[0] / q[2] < 0.9  # coarse C3D4: more than 10% too stiff
    assert lin[0] < lin[1] < lin[2]  # ...creeping up as the mesh is refined


def kt_of(r: dict, p) -> dict:
    """Hole Kt numbers from one solver result (pull load)."""
    pu, s = r["load_cases"]["pull"], theory.net_stress(p)
    prof = [v for _, v in pu["hole_sxx_profile"]]
    return {"peak": pu["hole_sxx_max"] / s, "face": pu["hole_sxx_face_max"] / s, "mean": sum(prof) / len(prof) / s}


def check_kt(results: dict) -> None:
    """Shared by the CalculiX and FEniCSx tests. results: {(case, level): result}.

    Heywood's Kt is a plane-stress (2D) number. In a plate of real thickness the
    hoop stress at the hole peaks at mid-thickness and drops where the hole
    meets the faces. So the tight check is the thin plate (t/d = 0.2, measured
    +1.7%); the 8 mm bracket (t/d = 0.8) must show the 3D effect: peak above
    Heywood (measured +5.4%), faces below (-8%), thickness average close (+3%).
    """
    ref = theory.heywood_kt(BRACKET)
    thin = kt_of(results[("thin", "fine")], THIN)
    thick = kt_of(results[("bracket", "fine")], BRACKET)
    thick_medium = kt_of(results[("bracket", "medium")], BRACKET)
    assert abs(thin["peak"] / ref - 1) < 0.025, thin
    assert abs(thin["mean"] / ref - 1) < 0.015, thin
    assert thick["face"] < ref < thick["mean"] < thick["peak"], thick
    assert thick["peak"] > thin["peak"] and thick["face"] < thin["face"]
    assert abs(thick["mean"] / ref - 1) < 0.04, thick
    assert abs(thick["peak"] / ref - 1) < 0.07, thick
    assert abs(thick["peak"] / thick_medium["peak"] - 1) < 0.01  # converged: medium vs fine


@pytestmark_solve
def test_hole_kt_against_heywood():
    check_kt({(c, lv): pipeline.ensure_ccx(c, lv) for c, lv in [("thin", "fine"), ("bracket", "fine"), ("bracket", "medium")]})
