"""The known answers, checked by hand before anything is compared with them."""
import math

from cadfea import theory
from cadfea.params import BEAM, BRACKET, STEEL, STOCKY


def test_cowper_shear_factor():
    assert math.isclose(theory.cowper_k(0.3), 13.0 / 15.3, rel_tol=1e-12)  # 0.8497
    assert math.isclose(theory.cowper_k(0.0), 10.0 / 12.0, rel_tol=1e-12)  # 5/6 at nu = 0


def test_beam_formulas_by_hand():
    # I = 20^4/12 = 13,333.3 mm^4; P L^3 / (3 E I) = 1000 * 200^3 / (3 * 210000 * 13333.3)
    assert math.isclose(theory.euler_bernoulli_tip(BEAM), 8e9 / 8.4e9, rel_tol=1e-12)
    G = 210000 / 2.6
    shear = 1000 * 200 / (13.0 / 15.3 * G * 400)
    assert math.isclose(theory.shear_tip(BEAM), shear, rel_tol=1e-12)
    assert math.isclose(theory.timoshenko_tip(BEAM), 8e9 / 8.4e9 + shear, rel_tol=1e-12)


def test_shear_matters_more_for_the_stocky_beam():
    s = theory.summary()
    assert 0.007 < s["beam"]["shear_share"] < 0.008  # 0.76% at L/h = 10
    assert 0.045 < s["stocky"]["shear_share"] < 0.047  # 4.6% at L/h = 4
    # The shear share scales with (h/L)^2, roughly.
    ratio = s["stocky"]["shear_share"] / s["beam"]["shear_share"]
    assert 5.5 < ratio < (10 / 4) ** 2 + 0.1


def test_heywood():
    assert math.isclose(theory.heywood_kt(BRACKET), 2 + 0.75**3, rel_tol=1e-12)  # 2.4219
    assert math.isclose(theory.net_stress(BRACKET), 10000 / (30 * 8), rel_tol=1e-12)


def test_material():
    assert math.isclose(STEEL.G, 80769.2307692, rel_tol=1e-9)
    assert STOCKY.L / STOCKY.h == 4.0 and BEAM.L / BEAM.h == 10.0
