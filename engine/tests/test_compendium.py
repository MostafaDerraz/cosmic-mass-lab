"""E1.6: the formulas of the compendium that the engine uses (docs/formulas-compendium/,
section "Cosmic Mass"). Each test names the paragraph of the formula it checks; the
compendium has no equation numbers. Formulas of the section "General Numbers" are tested
in test_gn.py, test_gnarray.py and test_readings_matter.py."""

import numpy as np
import pytest

from cosmicmass import readings as R
from cosmicmass.constants import SI
from cosmicmass.frame import Frame
from cosmicmass.gnarray import GNArray
from cosmicmass.matter import PRESETS, preset

SWEEP = [Frame.si(deg) for deg in range(10, 171, 5)]  # setting U1: sweep of the frame angle
PAIRS = [(a, b) for a in PRESETS.values() for b in PRESETS.values()]


@pytest.mark.parametrize("f", [Frame.si(), Frame.si(90.0), Frame.si(10.14)])
def test_full_product_of_two_cosmic_masses(f):
    # "Full product (both coordinates)":
    # mu_a mu_b = (m_a m_b - s^2 q_a q_b) + (m_a q_b + q_a m_b + 2 s cos(phi_bar) q_a q_b) G
    for a, b in PAIRS:
        prod = a.mu(f) * b.mu(f)
        assert R.Re(prod) == a.m * b.m - f.s * f.s * (a.q * b.q)
        coupling, share = a.m * b.q + a.q * b.m, 2 * f.s * f.cos_phi() * (a.q * b.q)
        assert R.Com(prod) == pytest.approx(coupling + share, rel=1e-14, abs=0)


def test_real_part_is_newton_plus_coulomb():
    # "Cosmic Interaction Law and its reading": F = G Re(mu_a mu_b)/r^2 = (G m_a m_b - K q_a q_b)/r^2
    f = Frame.si()
    for a, b in PAIRS:
        newton_coulomb = SI.G_N * a.m * b.m - SI.K * a.q * b.q
        assert SI.G_N * R.Re(a.mu(f) * b.mu(f)) == pytest.approx(newton_coulomb, rel=1e-14, abs=0)
    # two protons repel (Coulomb wins), a proton and an electron attract, two neutrons attract
    p, e, n = preset("proton").mu(f), preset("electron").mu(f), preset("neutron").mu(f)
    assert R.Re(p * p) < 0 < R.Re(p * e) and R.Re(n * n) > 0


def test_real_part_does_not_depend_on_the_frame_angle():
    # "... at every angle": Re[(a + bG)(c + dG)] = ac - s^2 bd
    for a, b in PAIRS:
        assert len({R.Re(a.mu(f) * b.mu(f)) for f in SWEEP}) == 1
    # while the Composer coordinate does depend on it, through the charge share
    p = preset("proton")
    assert len({R.Com(p.mu(f) * p.mu(f)) for f in SWEEP}) == len(SWEEP)


@pytest.mark.parametrize("f", [Frame.si(), Frame.si(90.0), Frame.si(120.0), Frame.scaled(79.86)])
def test_full_product_of_a_cosmic_mass_with_the_cosmic_field(f):
    # "Cosmic Field Law and law of motion", with E_cal = g + e G and Maxwell's E = -s^2 e;
    # docs/05, section 0: mu E_cal = [m g + q E] + [q g - m E/s^2 - (2 cos(phi_bar)/s) q E] G
    rng = np.random.default_rng(5)
    n, s = 500, f.s
    m, q = rng.uniform(0.1, 5, n), rng.uniform(-5, 5, n) / s
    g, e = rng.uniform(-5, 5, (n, 3)), rng.uniform(-5, 5, (n, 3)) / s
    E = -s * s * e
    prod = GNArray(m, q, f)[:, None] * GNArray(g, e, f)
    m_, q_ = m[:, None], q[:, None]
    assert np.allclose(R.Re(prod), m_ * g + q_ * E, rtol=1e-12, atol=1e-12)
    composer = q_ * g - m_ * E / s**2 - (2 * f.cos_phi() / s) * q_ * E
    assert np.allclose(R.Com(prod) * s, composer * s, rtol=1e-12, atol=1e-12)


def test_law_of_motion_real_part_for_a_body_at_rest():
    # "Cosmic Field Law and law of motion": Re[mu E_cal] = m g + q E, the Newton and Lorentz
    # forces on a body at rest; a neutral body feels g only, at every angle
    for f in SWEEP:
        field = GNArray([[0.0, -9.81, 0.0]], [[2e-15, 0.0, 0.0]], f)  # g and e
        neutron, proton = preset("neutron").mu(f), preset("proton").mu(f)
        assert np.array_equal(R.Re(neutron * field), preset("neutron").m * field.re)
        E = -f.s * f.s * field.com
        expected = preset("proton").m * field.re + preset("proton").q * E
        assert np.allclose(R.Re(proton * field), expected, rtol=1e-14, atol=0)


def test_two_stoney_masses_of_one_elementary_charge_have_a_zero_real_product():
    # "Stoney mass": m_S = e s; for mu = m_S + e G, Re(mu mu) = m_S^2 - s^2 e^2 = 0:
    # gravity and electricity balance exactly between two such bodies
    f = Frame.si()
    mu = GNArray([SI.stoney_mass], [SI.e], f)
    assert abs(R.Re(mu * mu)[0]) <= 1e-15 * SI.stoney_mass**2
