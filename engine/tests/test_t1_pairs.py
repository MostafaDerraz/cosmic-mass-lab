"""T1.1 to T1.6 (`scenarios/T1/`): the cases, the event that ends a run, the two passes through the
unchanged runner, the measurements, and the files of a test. Pre-registration: `team/notes/T1.md`.
The runs here use a coarse step (they check the protocol, not its numbers)."""

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from cosmicmass.constants import SI
from cosmicmass.runner import run
from cosmicmass.settings import load_settings

ROOT = Path(__file__).resolve().parents[2]


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


pairs = _load("pairs", ROOT / "scenarios" / "T1" / "pairs.py")
t1 = _load("t1_run", ROOT / "scenarios" / "T1" / "run.py")
validate = _load("validate_results", ROOT / "tools" / "validate_results.py")

COARSE = 1e-2


def test_ten_cases_in_six_tests_and_their_scales():
    assert list(pairs.TESTS) == ["T1.1", "T1.2", "T1.3", "T1.4", "T1.5", "T1.6"]
    assert [case for cases in pairs.TESTS.values() for case in cases] == list(pairs.CASES)
    moon = pairs.scales("T1.1")
    assert moon["T_ref"] == pytest.approx(2.357430e6, rel=1e-6) and moon["close_event"] == "plunge"
    assert moon["close"] == pytest.approx(3.844e7)
    earth = pairs.scales("T1.2a")  # D46: the proton reaches the Earth's surface before r0 / 10
    assert earth["close_event"] == "contact" and earth["close"] == pytest.approx(6.371e6)
    assert earth["a_ref"] == pytest.approx(8.134473, rel=1e-6)
    hydrogen = pairs.scales("T1.4")
    assert hydrogen["v_ref"] == pytest.approx(2.188287e6, rel=1e-6)
    a, b = pairs.bodies("T1.5")  # on the diagonal: Re(mu_a mu_b) = m_a m_b - s^2 q_a q_b = 0 to rounding
    assert abs(SI.G_N * a.m * b.m - SI.K * a.q * b.q) < 1e-14 * SI.K * a.q * b.q


@pytest.mark.parametrize("case", list(pairs.CASES))
def test_the_settings_of_every_case_load(case):
    loaded = load_settings(pairs.settings(case, "exploratory", producer="test"))
    state, sc = loaded.state, pairs.scales(case)
    assert loaded.law_spec.reading == "L1a" and loaded.numerics["softening"] == 0.0
    assert np.linalg.norm(state.x[1] - state.x[0]) == pytest.approx(sc["r0"])
    assert np.linalg.norm(state.v[1] - state.v[0]) == pytest.approx(sc["beta"] * sc["v_ref"])
    momentum = np.sum(state.m[:, None] * state.v, axis=0)
    assert np.allclose(momentum, 0.0, atol=1e-12 * np.sum(state.m * np.linalg.norm(state.v, axis=1)) + 1e-300)
    assert list(state.radius) == [pairs.RADII[body.name] for body in pairs.bodies(case)]  # D46
    assert pairs.settings(case, steps=7)["numerics"]["steps"] == 7
    assert pairs.settings(case, eta=5e-4)["numerics"]["steps"] == 2 * pairs.STEP_CAP


def test_a_law_block_a_name_or_a_definition_of_t2():
    block = {"id": "L1b", "sign": "-", "direction": "perpendicular"}
    assert pairs.settings("T1.3a", block)["law"] == block
    assert pairs.settings("T1.3a", {"id": "L1b-perpendicular.Re", "law": block})["law"] == block
    assert pairs.settings("T1.3a", "reference")["kind"] == "mathematics"
    assert pairs.settings("T1.3a", block)["kind"] == "finding"
    assert pairs.law_label("reference") == "L0" and pairs.law_label("exploratory") == "L1a"
    assert pairs.law_label(block) == "L1b (sign -, direction perpendicular)"


def test_the_first_acceleration_pull_push_and_swirl():
    assert pairs.first_acceleration("T1.3a")["what"] == "push"
    assert pairs.first_acceleration("T1.3a")["radial"] == pytest.approx(1.0)
    share = pairs.first_acceleration("T1.3b", "exploratory")  # the charge share: 1 - 2 cos(phi_bar)
    assert share["what"] == "pull" and share["radial"] == pytest.approx(-(1 - 2 * math.cos(math.radians(79.86))))
    proton = pairs.first_acceleration("T1.2a", "exploratory")  # s q g: s e / m_p times gravity
    assert proton["what"] == "pull" and -proton["radial"] == pytest.approx(SI.s * SI.e / 1.67262192e-27, rel=1e-9)
    assert pairs.first_acceleration("T1.2c", "exploratory")["what"] == "push"
    assert pairs.first_acceleration("T1.5")["what"] == "none"
    assert pairs.first_acceleration("T1.5", "exploratory")["radial"] == pytest.approx(-1.176053996, rel=1e-8)
    turned = pairs.first_acceleration("T1.5", {"id": "L1a", "direction": "perpendicular"})
    assert turned["what"] == "swirl" and abs(turned["tangential"]) == pytest.approx(1.176053996, rel=1e-8)
    # T1.1: two pure masses; the Composer part adds nothing at any angle
    assert pairs.first_acceleration("T1.1", "exploratory", 45.0) == pairs.first_acceleration("T1.1", "reference", 45.0)


def _row(t, x, y, vx=0.0, vy=0.0):
    return {"t": t, "x0_0": 0.0, "x0_1": 0.0, "x1_0": x, "x1_1": y, "v0_0": 0.0, "v0_1": 0.0, "v1_0": vx, "v1_1": vy}


class _Start:
    t = 0.0
    x = np.array([[0.0, 0.0], [1.0, 0.0]])
    v = np.zeros((2, 2))


def _watch(rows):
    watch = pairs.Watch(close=0.1, close_event="plunge", far=10.0)
    watch.start(_Start)
    for row in rows:
        assert not watch.stop()
        watch(row)
    return watch


def test_the_watch_names_the_event():
    assert _watch([_row(1.0, 0.5, 0.0), _row(2.0, 0.09, 0.0)]).event == "plunge"
    assert _watch([_row(1.0, 5.0, 0.0), _row(2.0, 10.0, 0.0)]).event == "escape"
    assert _watch([_row(1.0, math.nan, 0.0)]).event == "not finite"
    assert _watch([_row(1.0, 0.5, 0.0, vx=3.0)]).fastest == 3.0
    # a circle, in either sense: the event comes when the turn passes three revolutions
    for sense in (1, -1):
        angles = np.arange(1, 400) * 0.05
        rows = [_row(float(a), math.cos(a), sense * math.sin(a)) for a in angles]
        watch = pairs.Watch(close=0.1, close_event="plunge", far=10.0)
        watch.start(_Start)
        for row in rows:
            watch(row)
            if watch.stop():
                break
        assert watch.event == "three revolutions"
        assert sense * watch.track[-1][5] == pytest.approx(6 * math.pi, abs=0.05)
        assert sense * watch.track[-1][5] >= 6 * math.pi


def test_a_fall_to_contact_in_two_passes_through_the_plain_runner(tmp_path):
    out = pairs.run_case("T1.2a", "reference", eta=COARSE, out_root=tmp_path, producer="test")
    m = out.measured
    assert m["event"] == "contact" and m["t_event"] == pytest.approx(387.2847, rel=1e-3)
    assert m["first_acceleration"]["what"] == "pull" and m["revolutions"] == []
    assert m["separation_end_over_r0"] <= pairs.scales("T1.2a")["close"] / 7e6 < m["separation_max_over_r0"]
    assert m["third_law_residual_start"] < 1e-12 and m["spin_end"] == [0.0, 0.0]
    assert m["conservation"]["mass"] == 0 and m["conservation"]["momentum"] < 1e-9
    # the settings of the result are a plain run: the runner alone gives the same final state again
    again = run(out.settings, out_root=None)
    assert again.summary["steps"] == m["steps"] == out.settings["numerics"]["steps"]
    assert np.array_equal(again.state.x, out.result.state.x) and np.array_equal(again.state.v, out.result.state.v)
    # and the folder is an ordinary run folder (what the classifier of T2 reads)
    assert out.result.folder.parent == tmp_path / "T1.2a"
    assert validate.validate_folder(out.result.folder) == []
    assert {"stamp.json", "summary.json", "timeseries.csv"} <= {p.name for p in out.result.folder.iterdir()}


def test_an_orbit_turns_three_times_and_closes():
    out = pairs.run_case("T1.6", "reference", eta=COARSE)
    m = out.measured
    assert m["event"] == "three revolutions" and len(m["revolutions"]) == 3
    assert m["revolutions"][0]["t"] == pytest.approx(2.127206e-16, rel=1e-3)
    assert m["t_event"] == pytest.approx(m["revolutions"][2]["t"])
    assert max(turn["closure_error"] for turn in m["revolutions"]) < 2e-2
    apsides = m["apsides"]
    assert apsides["defined"] and apsides["inner_over_r0"][0] == pytest.approx(0.324503, rel=1e-3)
    assert apsides["largest_precession_rad"] < 2e-2 and m["angular_momentum_drift"] < 1e-9
    assert m["binding_energy_start_eV"] == pytest.approx(10.87173, rel=1e-5)


def test_no_force_no_event_and_the_composer_position_under_l1e():
    still = pairs.run_case("T1.5", "reference", eta=COARSE)
    assert still.measured["event"] == "no event" and still.measured["first_acceleration"]["what"] == "none"
    assert abs(still.measured["separation_max_over_r0"] - 1) < 1e-9
    assert still.measured["energy_change_over_reference"] < 1e-12
    two = pairs.run_case("T1.3a", {"id": "L1e"}, eta=COARSE)  # D49: xc reported beside x
    place = two.measured["composer_position_end"]
    assert np.any(np.array(place["xc"]) != 0) and np.array(place["x_minus_s2_q_over_m_xc"]).shape == (2, 2)
    assert "composer_position_end" not in still.measured


def test_the_torque_turns_the_spin_of_a_body_with_a_radius():
    out = pairs.run_case("T1.3a", {"id": "L1a", "direction": "torque"}, eta=COARSE)  # D46: protons have a radius
    assert all(spin != 0 for spin in out.measured["spin_end"])
    base = pairs.run_case("T1.3a", "reference", eta=COARSE)
    assert np.array_equal(out.result.state.x, base.result.state.x)  # the torque pushes nothing: L0's motion


def test_a_test_is_judged_and_written(tmp_path):
    lines = []
    outcome = t1.run_test("T1.3", producer="test", report=lines.append, eta=COARSE, only=("T1.3a",))
    assert set(outcome["cases"]) == {"T1.3a; L0", "T1.3a; L1a"} and len(lines) == 2
    for case in outcome["cases"].values():
        flags, block = case["flags"], case["convergence"]
        assert flags["event as predicted"] and flags["mass and charge exact"]
        assert flags["third law: residual below 1e-12"] and flags["momentum drift below 1e-6"]
        assert block["converged"] and block["eta"] == [COARSE, COARSE / 2, COARSE / 4]
        assert block["steps"][0] < block["steps"][1] < block["steps"][2]
        # at these coarse steps the escape is not yet second order (the ratio is 4 from eta = 2.5e-3 down)
        assert block["energy_drift_ratio_per_halving"][1] > 2.0
    folder = t1.write_test(outcome, tmp_path, "test")
    assert validate.validate_folder(folder) == []
    names = {p.name for p in folder.iterdir()}
    assert names == {"stamp.json", "summary.json", "timeseries-T1.3a__L0.csv", "timeseries-T1.3a__L1a.csv"}
    summary = json.loads((folder / "summary.json").read_text(encoding="utf-8"))
    assert summary["verdict"] in ("finding", "fail") and summary["question"] == t1.QUESTIONS["T1.3"]
    stamp = json.loads((folder / "stamp.json").read_text(encoding="utf-8"))
    assert stamp["cases"]["T1.3a; L1a"]["law"]["id"] == "L1a" and stamp["cases"]["T1.3a; L0"]["producer"] == "test"
    assert stamp["cases"]["T1.3a; L0"]["numerics"]["steps"] == summary["cases"]["T1.3a; L0"]["steps"]


def test_the_judge_marks_what_fails():
    measured = {
        "event": "plunge",
        "t_event": 1.0,
        "steps": 10,
        "eta": 1e-3,
        "first_acceleration": {"radial": -1.0, "tangential": 0.0, "what": "pull"},
        "revolutions": [],
        "apsides": {"defined": False},
        "conservation": {"mass": 0.0, "charge": 0.0, "momentum": 1e-3, "energy_total": 1e-9},
        "energy_used": "energy_total",
        "energy_change_over_reference": 0.0,
        "third_law_residual_start": 0.0,
        "third_law_residual_end": None,
        "separation_end_over_r0": 0.1,
        "binding_energy_start_eV": 1.0,
    }
    flags = t1.judge("T1.3b", "L0", [measured] * 3)
    assert flags["event as predicted"] and flags["first acceleration as predicted"]
    assert not flags["momentum drift below 1e-6"] and not flags["time of the event within 1e-4"]
    assert flags["third law: residual below 1e-12"] and flags["energy drift below 1e-6 at eta/4"]
    assert set(t1.judge("T1.3b", "L1b", [measured] * 3)) == {  # a law with no prediction: conservation only
        "mass and charge exact",
        "momentum drift below 1e-6",
        "energy drift below 1e-6 at eta/4",
        "third law: residual below 1e-12",
    }
    assert set(t1.PREDICTED) == {(case, law) for case in pairs.CASES for law in ("L0", "L1a")}


def test_the_script_refuses_an_unknown_test():
    with pytest.raises(SystemExit):
        t1.main(["T9.9"])
