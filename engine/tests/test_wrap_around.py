"""Decision D36: the wrap-around boundary (setting U4) with the nearest image. Software checks of the
engine parts; the comparison with the alternatives is scenario D36 (`scenarios/D36/run.py`)."""

import json

import numpy as np
import pytest

from cosmicmass.constants import scaled
from cosmicmass.energy import energy_l0, interaction_energy
from cosmicmass.environment import Environment
from cosmicmass.fields import DirectSum, nearest_image, smallest_separation
from cosmicmass.frame import Frame
from cosmicmass.integrators import nearest_distances, nearest_neighbours, pair_approach_time
from cosmicmass.laws import make_law
from cosmicmass.runner import WRAP_LIMIT, run
from cosmicmass.settings import SettingsError, load_settings
from cosmicmass.stamp import build_stamp
from cosmicmass.state import State
from cosmicmass.universe import Universe

BOX = {"lower": [0.0, 0.0], "upper": [1.0, 2.0]}
PERIOD = np.array([1.0, 2.0])
FRAME, UNITS = Frame.scaled(), scaled()


def gas(n=12, seed=3):
    rng = np.random.default_rng(seed)
    return State(
        m=rng.uniform(0.5, 2.0, n),
        q=rng.normal(0.0, 0.3, n),
        x=rng.uniform(0.0, 1.0, (n, 2)) * PERIOD,
        v=rng.normal(0.0, 0.2, (n, 2)),
        fixed=np.zeros(n, dtype=bool),
        alive=np.ones(n, dtype=bool),
    )


def forces(state, period, softening=0.01):
    solver = DirectSum(constants=UNITS, softening=softening, period=period)
    return make_law("L0").forces(state, solver.at_bodies(state), FRAME).total()


def test_nearest_image():
    r = np.array([[0.9, 0.1], [-0.6, 1.5], [0.2, -0.3], [3.4, -4.1]])
    out = nearest_image(r, PERIOD)
    assert np.allclose(out, [[-0.1, 0.1], [0.4, -0.5], [0.2, -0.3], [0.4, -0.1]])
    assert np.all(np.abs(out) <= PERIOD / 2) and nearest_image(r, None) is r
    assert np.allclose((r - out) / PERIOD, np.round((r - out) / PERIOD))  # moved by whole periods only


def test_the_walls_wrap_the_bodies_and_touch_nothing_else():
    env = Environment("wrap", BOX["lower"], BOX["upper"])
    assert list(env.period) == [1.0, 2.0] and Environment("reflecting", [0, 0], [1, 1]).period is None
    assert Environment().period is None
    state = State(
        m=[1.0] * 5,
        q=[0.0] * 5,
        x=[[1.25, 0.5], [-0.25, 2.5], [0.5, 1.0], [3.5, -4.5], [1.5, 0.5]],
        v=[[1.0, 0.0], [-1.0, 1.0], [0.0, 0.0], [2.0, -2.0], [1.0, 1.0]],
        fixed=[False, False, False, False, True],
        alive=[True, True, True, False, True],
    )
    v = state.v.copy()
    assert env.apply(state, 0.1) == []  # no wall is hit: no event, nothing given to a wall
    assert np.allclose(state.x[:3], [[0.25, 0.5], [0.75, 0.5], [0.5, 1.0]])
    assert list(state.x[3]) == [3.5, -4.5] and list(state.x[4]) == [1.5, 0.5]  # absorbed and fixed bodies stay
    assert np.array_equal(state.v, v) and env.collected["bodies"] == 0
    edge = State(m=[1.0], q=[0.0], x=[[1.0, 2.0]], v=[[0.0, 0.0]], fixed=[False], alive=[True])
    env.apply(edge, 0.1)
    assert list(edge.x[0]) == [0.0, 0.0]  # the upper face is the lower face
    with pytest.raises(ValueError, match="needs its"):
        Environment("wrap")


def test_a_pair_across_the_edge_feels_its_true_separation():
    pair = State(
        m=[1.0, 2.0],
        q=[0.3, -0.2],
        x=[[0.02, 1.0], [0.98, 1.0]],
        v=np.zeros((2, 2)),
        fixed=[False] * 2,
        alive=[True] * 2,
    )
    near = State(
        m=[1.0, 2.0],
        q=[0.3, -0.2],
        x=[[0.02, 1.0], [-0.02, 1.0]],
        v=np.zeros((2, 2)),
        fixed=[False] * 2,
        alive=[True] * 2,
    )
    wrapped, open_space = forces(pair, PERIOD, 0.0), forces(near, None, 0.0)
    assert np.allclose(wrapped, open_space, rtol=1e-13) and wrapped[0, 0] < 0  # pulled through the edge
    assert np.allclose(wrapped[0], -wrapped[1], rtol=1e-13)
    far = forces(pair, None, 0.0)  # without the period: the long way round
    assert far[0, 0] > 0 and abs(far[0, 0] / wrapped[0, 0]) == pytest.approx((0.04 / 0.96) ** 2, rel=1e-12)
    inside = State(
        m=[1.0, 1.0], q=[0.0, 0.0], x=[[0.4, 1.0], [0.6, 1.2]], v=np.zeros((2, 2)), fixed=[False] * 2, alive=[True] * 2
    )
    assert np.array_equal(forces(inside, PERIOD), forces(inside, None))  # nothing changes inside the box
    assert smallest_separation(pair, PERIOD) == pytest.approx(0.04) and smallest_separation(pair) == pytest.approx(0.96)
    assert interaction_energy(pair, FRAME, UNITS, period=PERIOD) == pytest.approx(
        interaction_energy(near, FRAME, UNITS)
    )
    solver = DirectSum.with_default_softening(pair, constants=UNITS, period=PERIOD)
    assert solver.softening == pytest.approx(0.04e-3) and list(solver.period) == [1.0, 2.0]


def test_forces_are_pairwise_and_do_not_depend_on_where_the_box_is_cut():
    state = gas()
    f = forces(state, PERIOD)
    assert np.linalg.norm(f.sum(axis=0)) < 1e-12 * np.sum(np.linalg.norm(f, axis=1))  # Newton's third law
    energy = energy_l0(state, FRAME, UNITS, 0.01, period=PERIOD)["total"]
    env = Environment("wrap", BOX["lower"], BOX["upper"])
    for shift in ([0.37, 0.0], [0.0, -1.21], [5.3, 7.9]):
        moved = state.copy()
        moved.x = moved.x + shift
        env.apply(moved, 0.0)  # the same bodies, the box cut elsewhere
        assert np.all((moved.x >= 0) & (moved.x < PERIOD))
        assert np.allclose(forces(moved, PERIOD), f, rtol=1e-9, atol=1e-12)
        assert energy_l0(moved, FRAME, UNITS, 0.01, period=PERIOD)["total"] == pytest.approx(energy, rel=1e-12)
    assert not np.allclose(forces(state, None), f)  # open space is another problem


def test_the_step_rule_measures_distances_through_the_edge():
    state = State(
        m=[1.0, 1.0, 1.0],
        q=[0.0, 0.0, 0.0],
        x=[[0.02, 1.0], [0.98, 1.0], [0.5, 1.0]],
        v=[[1.0, 0.0], [-1.0, 0.0], [0.0, 0.0]],
        fixed=[False] * 3,
        alive=[True] * 3,
    )
    reach, partner = nearest_neighbours(state, PERIOD)
    assert np.allclose(reach, [0.04, 0.04, 0.48]) and list(partner[:2]) == [1, 0]
    assert np.allclose(nearest_distances(state, PERIOD), reach) and np.allclose(
        nearest_distances(state), [0.48, 0.48, 0.48]
    )
    assert pair_approach_time(state, PERIOD) == pytest.approx(0.04 / 2.0) and pair_approach_time(
        state
    ) == pytest.approx(0.48)


def test_bad_periods_and_universes_are_refused():
    with pytest.raises(ValueError, match="period"):
        DirectSum(period=[1.0, -1.0])
    with pytest.raises(ValueError, match="period"):
        DirectSum(period=1.0)
    with pytest.raises(ValueError, match="3 sides"):
        DirectSum(period=[1.0, 1.0, 1.0]).at_bodies(gas())
    u = Universe.from_dict({"boundary": "wrap", "box": BOX})
    assert u.period == [1.0, 2.0] and Universe.from_dict({"boundary": "reflecting", "box": BOX}).period is None
    assert u.environment().boundary == "wrap"
    inflow = {"preset": "neutron", "rate": 1.0, "velocity": [1.0, 0.0], "face": {"axis": 0, "side": "lower"}}
    with pytest.raises(SettingsError, match="wrap-around box has none"):
        load_settings({"scenario": "x", "universe": {"boundary": "wrap", "box": BOX, "inflow": inflow}, "bodies": []})
    with pytest.raises(SettingsError, match="universe.box"):
        load_settings(
            {"scenario": "x", "universe": {"boundary": "wrap"}, "bodies": [{"m": 1.0, "q": 0.0, "x": [0.0, 0.0]}]}
        )


def pair_settings(boundary):
    a, speed, drift = 0.02, 10.0, 3.0  # not 5: with it the kinetic and the binding energy cancel exactly
    block = {"units": "scaled", "boundary": boundary}
    if boundary == "wrap":
        block["box"] = {"lower": [0.0, 0.0], "upper": [1.0, 1.0]}
    return {
        "scenario": "pair",
        "universe": block,
        "bodies": [
            {"m": 1.0, "q": 0.0, "x": [0.95 - a / 2, 0.5], "v": [drift, -speed / 2]},
            {"m": 1.0, "q": 0.0, "x": [0.95 + a / 2, 0.5], "v": [drift, speed / 2]},
        ],
        "numerics": {"softening": 0.0, "t_end": 0.02, "eta": 0.01, "record_every": 20},
        "kind": "mathematics",
    }


def test_a_run_in_a_wrap_around_box(tmp_path):
    """A bound pair carried through the edge by the runner: it does not notice the edge."""
    wrapped, free = run(pair_settings("wrap"), out_root=tmp_path), run(pair_settings("open"), out_root=None)
    assert wrapped.summary["steps"] == free.summary["steps"] and wrapped.state.t == free.state.t == 0.02
    assert np.all((wrapped.state.x >= 0.0) & (wrapped.state.x < 1.0)) and np.any(free.state.x > 1.0)  # it crossed
    relative = nearest_image(wrapped.state.x[1] - wrapped.state.x[0], np.array([1.0, 1.0]))
    assert np.allclose(relative, free.state.x[1] - free.state.x[0], rtol=0, atol=1e-9)
    assert np.allclose(wrapped.state.v, free.state.v, rtol=0, atol=1e-7)
    assert wrapped.summary["final"]["energy_re"] == pytest.approx(free.summary["final"]["energy_re"], rel=1e-9)
    assert wrapped.summary["conservation"]["energy_re"] < 1e-5  # no jump at the crossing
    assert wrapped.summary["conservation"]["momentum"] < 1e-13 and wrapped.summary["conservation"]["mass"] == 0.0
    # the stated limit of D36 is written into the result; a wrap-around box has no wall pressure
    note = wrapped.summary["wrap_around"]
    assert (
        note == {"field": "nearest image (D36)", "period": [1.0, 1.0], "limit": WRAP_LIMIT}
        and "half the box" in WRAP_LIMIT
    )
    assert wrapped.summary["walls"] is None and wrapped.summary["wall_hits"] == 0
    assert free.summary["wrap_around"] is None
    stamp = json.loads((wrapped.folder / "stamp.json").read_text(encoding="utf-8"))
    assert stamp["universe"]["boundary"] == "wrap" and stamp["universe"]["box"]["upper"] == [1.0, 1.0]
    assert stamp == json.loads(
        json.dumps(build_stamp(load_settings(pair_settings("wrap")), 0.0, stamp["run_id"], stamp["engine_commit"]))
    )


def test_the_default_softening_uses_the_nearest_image():
    settings = pair_settings("wrap")
    settings["bodies"][0]["x"], settings["bodies"][1]["x"] = [0.01, 0.5], [0.99, 0.5]
    settings["numerics"] = {"steps": 1}
    assert run(settings, out_root=None).stamp["numerics"]["softening"] == pytest.approx(1e-3 * 0.02)
