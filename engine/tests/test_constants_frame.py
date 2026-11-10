"""E1.1: constants and Frame. Formulas cited from the compendium (docs/formulas-compendium/),
section "Cosmic Mass", paragraphs "The cosmic mass and its unit" and "Stoney mass"."""

import math

import mpmath
import pytest

from cosmicmass import constants
from cosmicmass.constants import SCALED, SI, Constants, scaled
from cosmicmass.frame import DEFAULT_PHI_BAR_DEG, Frame


def test_s_is_sqrt_K_over_G_to_six_digits():
    # "The cosmic mass and its unit": s = sqrt(K/G) = 1.16043e10 kg/C
    assert SI.s == math.sqrt(SI.K / SI.G_N)
    assert f"{SI.s:.5e}" == "1.16043e+10"


def test_stoney_mass():
    # "Stoney mass": m_S = e s = 1.85921e-9 kg
    assert SI.stoney_mass == SI.e * SI.s
    assert f"{SI.stoney_mass:.5e}" == "1.85921e-09"


def test_si_values_agree_with_the_verified_constants():
    # same numbers as verification/cosmic-mass/constants.py (CODATA 2022)
    assert SI.G_N == 6.67430e-11
    assert SI.c == 299792458.0
    assert SI.e == 1.602176634e-19
    assert SI.K == pytest.approx(8.9875517862e9, rel=1e-10)
    assert SI.units == "SI"


def test_scaled_units_have_s_equal_to_one():
    # decision D7: G_N = K = 1, so s = 1
    assert (SCALED.G_N, SCALED.K, SCALED.s, SCALED.units) == (1.0, 1.0, 1.0, "scaled")
    assert scaled(c=1e3).s == 1.0 and scaled(c=1e3).c == 1e3


def test_constants_refuse_bad_input():
    with pytest.raises(ValueError):
        Constants(G_N=1.0, K=1.0, c=1.0, e=1.0, units="cgs")
    with pytest.raises(ValueError):
        Constants(G_N=-1.0, K=1.0, c=1.0, e=1.0, units="SI")
    assert constants.SI is SI


def test_si_and_scaled_frames():
    f = Frame.si()
    assert f.s == SI.s and f.units == "SI"
    assert f.phi_bar_deg == pytest.approx(DEFAULT_PHI_BAR_DEG, rel=1e-15)
    g = Frame.scaled()
    assert g.s == 1.0 and g.units == "scaled" and g.phi_bar == math.pi / 2


def test_angle_is_stored_in_radians_and_read_and_written_in_degrees():
    f = Frame.from_degrees(2.0, 60.0)
    assert f.phi_bar == math.radians(60.0)
    assert f.phi_bar_deg == pytest.approx(60.0, rel=1e-15)
    d = f.to_dict()
    assert d == {"s": 2.0, "phi_bar_deg": pytest.approx(60.0, rel=1e-15), "units": "SI"}
    assert Frame.from_dict(d).phi_bar == pytest.approx(f.phi_bar, rel=1e-15)


def test_square_rule_coefficients():
    # "Rule of the Square and the product": G^2 = 2 s cos(phi_bar) G - s^2
    f = Frame.from_degrees(3.0, 60.0)
    p, q = f.square_rule()
    assert p == -9.0
    assert q == pytest.approx(3.0, rel=1e-15)  # 2 * 3 * cos(60 deg)
    # "Generalized complex numbers": q^2 + 4p < 0 (the elliptic case), at every angle
    for deg in range(10, 171, 5):
        p, q = Frame.from_degrees(SI.s, deg).square_rule()
        assert q * q + 4 * p < 0


@pytest.mark.parametrize("deg", [0.0, 180.0, -10.0, 200.0])
def test_frame_refuses_angles_outside_0_180(deg):
    with pytest.raises(ValueError):
        Frame.from_degrees(1.0, deg)


def test_frame_refuses_bad_length_and_units():
    with pytest.raises(ValueError):
        Frame(0.0, 1.0)
    with pytest.raises(ValueError):
        Frame(1.0, 1.0, "cgs")


def test_exact_frame_uses_mpmath():
    with mpmath.workdps(50):
        f = Frame.from_degrees(mpmath.mpf(2), 60)
        assert f.exact and not Frame.si().exact
        assert isinstance(f.cos_phi(), mpmath.mpf)
        assert abs(f.cos_phi() - mpmath.mpf(1) / 2) < mpmath.mpf(10) ** -49
        assert abs(f.phi_bar_deg - 60) < mpmath.mpf(10) ** -47
