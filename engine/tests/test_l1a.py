"""E4.3 (L1a only): the reading L1a of the Composer coordinate, a hypothesis (D1).
docs/05, section 1: force = Re(mu E_cal) + lambda s Com(mu E_cal); compendium, "Full product"."""

import numpy as np
import pytest

from cosmicmass import readings as R
from cosmicmass.fields import FieldSample
from cosmicmass.frame import Frame
from cosmicmass.laws import LawSpec, available_laws, make_law
from cosmicmass.laws.l1 import L1a
from cosmicmass.state import State


def random_case(d, seed=21, n=30, s=None):
    rng = np.random.default_rng(seed)
    f = Frame.si() if s is None else Frame.from_degrees(s, 60.0)
    second = (n, 3) if d == 3 else (n,)
    fields = FieldSample(rng.normal(size=(n, d)), rng.normal(size=(n, d)) / f.s,
                         rng.normal(size=second), rng.normal(size=second) / f.s)
    state = State(m=rng.uniform(1, 2, n), q=rng.normal(size=n) / f.s, x=rng.normal(size=(n, d)),
                  v=rng.normal(size=(n, d)), fixed=np.zeros(n, bool), alive=np.ones(n, bool))
    return f, state, fields


def test_l1a_is_registered_as_a_hypothesis():
    assert {"L0", "L1a"} <= set(available_laws())
    law = make_law(LawSpec(reading="L1a"))
    assert isinstance(law, L1a) and law.hypothesis is True and law.id == "L1a"
    for option in ({"direction": "perpendicular"}, {"moving": "on"}, {"range": "screened"}):
        with pytest.raises(NotImplementedError):
            make_law(LawSpec(reading="L1a", **option))


@pytest.mark.parametrize("d", [2, 3])
def test_all_terms_equal_lambda_s_times_the_composer_coordinate_of_the_product(d):
    f, state, fields = random_case(d)
    pair = make_law(LawSpec(reading="L1a", strength=0.25)).forces(state, fields, f)
    product = state.mu(f)[:, None] * fields.E_cal(f)  # mu E_cal, both coordinates
    assert np.allclose(pair.composer, 0.25 * f.s * R.Com(product), rtol=1e-12, atol=1e-12)
    # the Real part is L0, untouched: the same array, bit for bit
    assert np.array_equal(pair.real, make_law("L0").forces(state, fields, f).real)
    assert np.array_equal(pair.total(), pair.real + pair.composer)


@pytest.mark.parametrize("d", [2, 3])
def test_reduces_to_L0_exactly_when_the_strength_is_zero(d):
    # the reduce-to-L0 test of ARCHITECTURE.md, section 7
    f, state, fields = random_case(d)
    reference = make_law("L0").forces(state, fields, f).total()
    for terms in ("all", "s q g", "-m E/s", "-2 cos q E"):
        for sign in "+-":
            law = make_law(LawSpec(reading="L1a", strength=0.0, terms=terms, sign=sign))
            assert np.array_equal(law.forces(state, fields, f).total(), reference)
            assert np.array_equal(law.inertia(state, f), state.m)


def test_the_three_terms_by_hand():
    # one body, m = 3, q = 5, in a frame s = 2, phi_bar = 60 deg, with g = (7, 0) and e = (0, 11)
    f = Frame.from_degrees(2.0, 60.0)
    state = State(m=[3.0], q=[5.0], x=[[0.0, 0.0]], v=[[0.0, 0.0]], fixed=[False], alive=[True])
    fields = FieldSample(np.array([[7.0, 0.0]]), np.array([[0.0, 11.0]]), np.zeros(1), np.zeros(1))
    force = lambda **kw: make_law(LawSpec(reading="L1a", **kw)).forces(state, fields, f).composer[0]  # noqa: E731
    assert force(terms="s q g") == pytest.approx([2 * 5 * 7, 0.0])  # s q g
    assert force(terms="-m E/s") == pytest.approx([0.0, 2 * 3 * 11])  # s m e = -m E / s with E = -s^2 e
    assert force(terms="-2 cos q E") == pytest.approx([0.0, 2 * 4 * 0.5 * 5 * 11])  # 2 s^2 cos q e
    assert force() == pytest.approx([70.0, 66.0 + 220.0])
    assert force(sign="-") == pytest.approx([-70.0, -286.0]) and force(strength=0.1) == pytest.approx([7.0, 28.6])
    assert force(strength="cos_phi_bar") == pytest.approx([35.0, 143.0])  # lambda = cos(60 deg)


def test_the_charge_share_vanishes_at_ninety_degrees_and_the_rest_does_not_depend_on_the_angle():
    _, state, fields = random_case(2)
    at = lambda deg, terms: make_law(LawSpec(reading="L1a", terms=terms)).forces(state, fields, Frame.si(deg)).composer  # noqa: E731
    assert np.allclose(at(90.0, "-2 cos q E"), 0.0, atol=1e-14 * np.abs(at(10.0, "-2 cos q E")).max())
    for terms in ("s q g", "-m E/s"):
        assert np.array_equal(at(10.0, terms), at(170.0, terms))
    assert not np.allclose(at(10.0, "all"), at(170.0, "all"))  # L1a as a whole depends on the angle
