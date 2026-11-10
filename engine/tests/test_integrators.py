"""E3.1, E3.3, E6a: velocity Verlet, the adaptive step and reflecting walls, under the
reference law L0 (compendium, "Cosmic Field Law and law of motion"). The energy used
here is Newton plus Coulomb: kinetic energy minus G_N Re(mu_a mu_b) / r for each pair."""

import math

import numpy as np
import pytest

from cosmicmass.constants import SCALED, SI
from cosmicmass.environment import Environment
from cosmicmass.fields import DirectSum
from cosmicmass.frame import Frame
from cosmicmass.integrators import REST_FRACTION, StepContext, VelocityVerlet, nearest_distances
from cosmicmass.laws import Law, make_law
from cosmicmass.matter import preset
from cosmicmass.settings import SettingsError, load_settings
from cosmicmass.state import State

VV = VelocityVerlet()


def energy(s: State, c=SCALED) -> float:
    kinetic = 0.5 * np.sum(s.m * np.sum(s.v**2, axis=1))
    potential = 0.0
    for i in range(s.n):
        for j in range(i + 1, s.n):
            strength = c.G_N * s.m[i] * s.m[j] - c.K * s.q[i] * s.q[j]  # G_N Re(mu_i mu_j)
            potential -= strength / np.linalg.norm(s.x[i] - s.x[j])
    return float(kinetic + potential)


def context(constants=SCALED, **kw) -> StepContext:
    frame = Frame.scaled() if constants.units == "scaled" else Frame.si()
    return StepContext(make_law("L0"), DirectSum(constants=constants), frame, **kw)


def kepler(m1=1.0, m2=1e-3, a=1.0, e=0.0, q=(0.0, 0.0), c=SCALED):
    """Two bodies at apocentre in their centre-of-mass frame; returns the state and the period."""
    strength = c.G_N * m1 * m2 - c.K * q[0] * q[1]
    M, reduced = m1 + m2, m1 * m2 / (m1 + m2)
    k = strength / reduced  # the relative motion obeys r'' = -k r / r^3
    r, speed = a * (1 + e), math.sqrt(k * (1 - e) / (a * (1 + e)))
    x = np.array([[-m2 / M * r, 0.0], [m1 / M * r, 0.0]])
    v = np.array([[0.0, -m2 / M * speed], [0.0, m1 / M * speed]])
    state = State(m=[m1, m2], q=list(q), x=x, v=v, fixed=[False, False], alive=[True, True])
    return state, 2 * math.pi * math.sqrt(a**3 / k)


def orbit(state, period, periods, ctx, adaptive=False, steps_per_period=6283):
    """Advance and return (final state, relative energy change at the end, largest relative
    energy excursion on the way, measured period, number of steps)."""
    c = ctx.solver.constants
    e0, worst, angle, steps = energy(state, c), 0.0, 0.0, 0
    previous = math.atan2(*(state.x[1] - state.x[0])[::-1])
    while state.t < periods * period:
        dt = VV.next_dt(state, ctx) if adaptive else period / steps_per_period
        state = VV.step(state, dt, ctx)
        steps += 1
        now = math.atan2(*(state.x[1] - state.x[0])[::-1])
        angle += (now - previous + math.pi) % (2 * math.pi) - math.pi
        previous = now
        if steps % 97 == 0:
            worst = max(worst, abs(energy(state, c) - e0) / abs(e0))
    drift = abs(energy(state, c) - e0) / abs(e0)
    return state, drift, max(worst, drift), state.t * 2 * math.pi / abs(angle), steps


def test_kepler_orbit_keeps_its_energy_and_its_period():
    state, period = kepler()
    _, drift, excursion, measured, _ = orbit(state, period, 5, context())
    assert excursion < 1e-6 and drift < 1e-6  # E3.1, shortened to 5 periods
    assert measured == pytest.approx(period, rel=1e-3)


@pytest.mark.slow
def test_kepler_orbit_for_1000_periods():
    # E3.1 as written: 1000 periods, energy drift below 1e-6, period within 0.1%
    state, period = kepler()
    final, drift, excursion, measured, steps = orbit(state, period, 1000, context())
    print(f"\n1000 periods: steps={steps}, energy change at the end={drift:.3e}, "
          f"largest excursion={excursion:.3e}, period error={measured / period - 1:.3e}")
    assert excursion < 1e-6 and measured == pytest.approx(period, rel=1e-3)


def test_the_adaptive_step_follows_D10():
    state, period = kepler()
    ctx = context()
    dt = VV.next_dt(state, ctx)
    a = np.linalg.norm(ctx.acceleration(state), axis=1)
    speed = np.linalg.norm(state.v, axis=1)
    assert dt == pytest.approx(1e-3 * np.min(speed / a))  # velocity change per step: 1e-3 of the speed
    assert dt == pytest.approx(0.9995e-3, rel=1e-4)
    assert VV.next_dt(state, context(dt_max=dt / 10)) == dt / 10 and VV.next_dt(state, context(dt_min=dt * 10)) == dt * 10
    assert VV.next_dt(state, context(eta=1e-4)) == pytest.approx(dt / 10)


def test_eccentric_orbit_with_the_adaptive_step():
    # E3.3: an eccentric orbit (e = 0.9) keeps its energy drift below 1e-6.
    # "Drift" is read as the change of the energy at the end of the run. On the way the energy
    # makes an excursion of a few 1e-6 at each pericentre, where the step is 1.9e-3 of the local
    # turning time; it comes back. With eta = 3e-4 the excursion too stays below 1e-6.
    state, period = kepler(e=0.9)
    final, drift, excursion, measured, steps = orbit(state, period, 1, context(), adaptive=True)
    assert drift < 1e-6 and excursion < 1e-5 and measured == pytest.approx(period, rel=1e-3) and steps > 6283
    state, period = kepler(e=0.9)
    assert orbit(state, period, 1, context(eta=3e-4), adaptive=True)[2] < 1e-6


def test_hydrogen_like_orbit_in_SI_has_the_bohr_period():
    # T0.1(b) as a unit test: electron and proton at the Bohr radius, period 1.5198e-16 s within 0.1%
    e, p = preset("electron"), preset("proton")
    state, period = kepler(p.m, e.m, a=5.29177e-11, q=(p.q, e.q), c=SI)
    _, drift, excursion, measured, _ = orbit(state, period, 3, context(SI))
    assert measured == pytest.approx(1.5198e-16, rel=1e-3) and excursion < 1e-6


def test_bodies_starting_at_rest_get_a_free_fall_step():
    # two protons at rest (T0.1c): no speed, so the step comes from the floor of the speed
    p = preset("proton")
    state = State.from_bodies([p, p], [[0.0, 0.0], [1e-10, 0.0]])
    ctx = context(SI)
    a = np.linalg.norm(ctx.acceleration(state), axis=1)[0]
    dt = VV.next_dt(state, ctx)
    assert dt == pytest.approx(1e-3 * REST_FRACTION * math.sqrt(1e-10 / a))
    e0, x0 = energy(state, SI), state.x.copy()
    for _ in range(2000):
        state = VV.step(state, VV.next_dt(state, ctx), ctx)
    assert abs(energy(state, SI) - e0) / abs(e0) < 1e-6
    assert state.x[1, 0] - state.x[0, 0] > 1e-10  # they repel
    assert abs(np.sum(state.m[:, None] * state.v)) < 1e-12 * np.sum(state.m * np.linalg.norm(state.v, axis=1))
    assert np.allclose(state.x.mean(axis=0), x0.mean(axis=0), rtol=0, atol=1e-22)  # the centre does not move


def test_fixed_step_verlet_is_time_reversible_and_keeps_momentum():
    rng = np.random.default_rng(3)
    state = State(m=rng.uniform(1, 2, 5), q=rng.normal(0, 0.3, 5), x=rng.uniform(-1, 1, (5, 3)) * 3,
                  v=rng.normal(0, 0.2, (5, 3)), fixed=np.zeros(5, bool), alive=np.ones(5, bool))
    start, ctx = state.copy(), context()
    p0 = np.sum(state.m[:, None] * state.v, axis=0)
    for _ in range(300):
        state = VV.step(state, 1e-3, ctx)
    assert np.allclose(np.sum(state.m[:, None] * state.v, axis=0), p0, rtol=0, atol=1e-13)
    assert abs(energy(state) - energy(start)) / abs(energy(start)) < 1e-6
    state.v *= -1
    ctx.invalidate()
    for _ in range(300):
        state = VV.step(state, 1e-3, ctx)
    assert np.allclose(state.x, start.x, rtol=0, atol=1e-11) and state.t == pytest.approx(0.6)


def test_fixed_and_absorbed_bodies_do_not_move():
    state = State(m=[1.0, 1.0, 1.0], q=[0.0, 0.0, 0.0], x=[[0.0, 0.0], [1.0, 0.0], [0.0, 2.0]],
                  v=[[0.0, 0.0], [0.0, 1.0], [5.0, 5.0]], fixed=[True, False, False], alive=[True, True, False])
    ctx = context()
    assert np.array_equal(ctx.acceleration(state)[[0, 2]], np.zeros((2, 2)))
    assert ctx.acceleration(state)[1] == pytest.approx([-1.0, 0.0])  # pulled by the fixed body only
    assert nearest_distances(state)[2] == np.inf and nearest_distances(state)[0] == 1.0
    new = VV.step(state, 0.01, ctx)
    assert np.array_equal(new.x[[0, 2]], state.x[[0, 2]]) and np.array_equal(new.v[2], [5.0, 5.0])
    assert not np.array_equal(new.x[1], state.x[1]) and new is not state and state.t == 0.0


def test_no_derivable_step_and_zero_inertia_are_refused():
    alone = State.from_bodies([preset("neutron")], [[0.0, 0.0]])
    with pytest.raises(ValueError, match="dt_max"):
        VV.next_dt(alone, context())
    assert VV.next_dt(alone, context(dt_max=2.0)) == 2.0
    # at rest in a uniform background: no neighbour and no speed, so the step is dt_max
    falling = StepContext(make_law("L0"), DirectSum(constants=SCALED, background={"g": [0.0, -1.0]}), Frame.scaled(), dt_max=0.5)
    assert VV.next_dt(alone, falling) == 0.5
    # thrown upwards, it must pass its turning point: D10's rule alone would shrink the step to
    # nothing there, so a body with no neighbour does not limit the step
    thrown = State.from_bodies([preset("neutron")], [[0.0, 0.0]], v=[[0.0, 1.0]])
    falling.dt_max = 0.01
    for _ in range(300):
        assert VV.next_dt(thrown, falling) == 0.01
        thrown = VV.step(thrown, 0.01, falling)
    assert thrown.t == pytest.approx(3.0) and thrown.v[0, 1] == pytest.approx(-2.0)  # v = 1 - g t
    assert thrown.x[0, 1] == pytest.approx(3.0 - 0.5 * 9.0)  # exact in a uniform field

    class Weightless(Law):
        reading = "L0"

        def forces(self, state, fields, frame):
            return make_law("L0").forces(state, fields, frame)

        def inertia(self, state, frame):
            return np.zeros(state.n)

    state, _ = kepler()
    with pytest.raises(ZeroDivisionError, match="zero inertia"):
        StepContext(Weightless(), DirectSum(constants=SCALED), Frame.scaled()).acceleration(state)


def test_reflecting_walls_turn_bodies_back_without_changing_their_energy():
    state = State(m=[1.0, 2.0, 3.0, 4.0], q=[0.0] * 4, x=[[1.5, 0.0], [-1.2, 0.5], [0.0, 0.0], [1.5, 0.0]],
                  v=[[2.0, 1.0], [-1.0, 0.0], [3.0, 3.0], [-2.0, 0.0]], fixed=[False] * 4, alive=[True] * 4)
    before, speeds = state.copy(), np.linalg.norm(state.v, axis=1)
    events = Environment("reflecting", [-1.0, -1.0], [1.0, 1.0]).apply(state, 0.1)
    assert np.array_equal(state.v, [[-2.0, 1.0], [1.0, 0.0], [3.0, 3.0], [-2.0, 0.0]])  # body 3 was already returning
    assert np.array_equal(state.x, before.x) and np.array_equal(np.linalg.norm(state.v, axis=1), speeds)
    assert events == [{"wall": (0, "lower"), "bodies": 1, "impulse": -4.0}, {"wall": (0, "upper"), "bodies": 1, "impulse": 4.0}]
    assert Environment().apply(state, 0.1) == []
    state.fixed[0], state.v[0, 0] = True, 2.0
    assert Environment("reflecting", [-1.0, -1.0], [1.0, 1.0]).apply(state, 0.1) == []  # fixed bodies are left alone


def test_a_cloud_in_a_reflecting_box_keeps_its_energy():
    rng = np.random.default_rng(8)
    state = State(m=np.full(6, 1e-2), q=rng.normal(0, 2e-2, 6), x=rng.uniform(-0.8, 0.8, (6, 2)),
                  v=rng.normal(0, 0.7, (6, 2)), fixed=np.zeros(6, bool), alive=np.ones(6, bool))
    ctx, box = context(), Environment("reflecting", [-1.0, -1.0], [1.0, 1.0])
    e0, hits = energy(state), 0
    for _ in range(6000):
        dt = VV.next_dt(state, ctx)
        state = VV.step(state, dt, ctx)
        if box.apply(state, dt):
            hits += 1
            ctx.invalidate()
    assert hits >= 3 and abs(energy(state) - e0) / abs(e0) < 1e-6
    assert np.all(np.abs(state.x) < 1.0 + 0.05)  # bodies overshoot a wall by one step at most


def test_bad_boxes_are_refused():
    for bad in (dict(boundary="absorbing"), dict(boundary="reflecting"), dict(boundary="reflecting", lower=[0.0], upper=[0.0]),
                dict(boundary="reflecting", lower=[0.0, 0.0], upper=[1.0])):
        with pytest.raises(ValueError):
            Environment(**bad)
    with pytest.raises(ValueError):
        Environment("reflecting", [0.0] * 3, [1.0] * 3).apply(kepler()[0], 0.1)
    base = {"scenario": "x", "bodies": [{"preset": "proton", "x": [0.0, 0.0]}]}
    s = load_settings({**base, "universe": {"boundary": "reflecting", "box": {"lower": [-1, -1], "upper": [1, 1]}},
                       "numerics": {"eta": 1e-4, "dt_max": 1.0}})
    assert s.universe["box"] == {"lower": [-1, -1], "upper": [1, 1]} and s.numerics["eta"] == 1e-4
    for box in (None, {"lower": [0, 0]}, {"lower": [0], "upper": [1]}, {"lower": [1, 1], "upper": [0, 0]}, {"size": 2}):
        with pytest.raises(SettingsError, match="box"):
            load_settings({**base, "universe": {"boundary": "reflecting", "box": box}})
