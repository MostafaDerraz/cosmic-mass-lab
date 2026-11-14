"""T2.3: the classifier of the survival table (decisions D43, D44; pre-registration in team/notes/T2.md).

These tests check the method on cases whose answer is known in advance; they do not test the classes of
the grid, which are the result of T2."""

import json

import pytest

from cosmicmass.analysis import survival as sv
from cosmicmass.laws import LawSpec

L1A = {"id": "L1a", "sign": "+", "lambda": 1.0, "direction": "radial", "terms": "all", "B_cal": "off",
       "inertia": "Re", "range": "1/r^2"}  # fmt: skip
TURNED = dict(L1A, direction="perpendicular")


def test_d23_scope():
    assert sv.in_d23_scope(L1A)
    assert not sv.in_d23_scope(TURNED) and not sv.in_d23_scope(dict(L1A, inertia="Sh"))
    assert not sv.in_d23_scope(dict(L1A, id="L1b")) and not sv.in_d23_scope({"id": "L0"})


def test_the_loop_test_passes_a_gradient_and_fails_a_curl():
    # L0 is conservative (compendium); a central field turned by 90 degrees is not: its work round a loop
    # is of the order of the work along one side
    assert sv.closed_loop_ratio({"id": "L0"}, 79.86) <= 1e-12
    assert sv.closed_loop_ratio(L1A, 45.0) <= 1e-12
    assert sv.closed_loop_ratio(TURNED, 45.0) > 1e-3
    assert sv.closed_loop_ratio(dict(L1A, direction="velocity"), 45.0) > 1e-3


def record(setup, angle, mass=0.0, momentum=0.0, energy=0.0, x=1.0, definition="d"):
    summary = {"conservation": {"mass": mass, "charge": 0.0, "momentum": momentum, "energy_total": energy},
               "momentum_with_walls_drift": momentum,
               "final": {"energy_re": -1.0, "kinetic": 0.5, "momentum": [0.0, 0.0]}}  # fmt: skip
    return {"definition": definition, "setup": setup, "angle": angle, "summary": summary,
            "final": {"x0_0": x, "v0_0": 2.0}}  # fmt: skip


def runs(**change):
    return [record(s, a, **change) for s in sv.SETUPS for a in sv.ANGLES]


LOOPS_PASS = {a: 1e-15 for a in sv.ANGLES}


def test_classes_from_hand_made_runs():
    assert sv.classify("d", L1A, runs())["class"] == "A"  # D23 scope: the energy of the runs decides
    assert sv.classify("d", TURNED, runs(), LOOPS_PASS)["class"] == "A"  # outside: the loop test decides
    tilted = runs()
    for r in tilted:
        if r["setup"] == "T1.4" and r["angle"] == 90.0:
            r["final"]["x0_0"] = 1.0 + 1e-6  # one place moves with the angle
    assert sv.classify("d", L1A, tilted)["class"] == "B"
    out = sv.classify("d", L1A, runs(momentum=1e-5))
    assert out["class"] == "C" and all(f.startswith("(b)") for f in out["failed"])
    assert sv.classify("d", L1A, runs(energy=1e-4))["failed"][0].startswith("(c) energy")
    assert sv.classify("d", L1A, runs(mass=1e-12))["failed"][0].startswith("(a)")
    curl = {**LOOPS_PASS, 45.0: 0.5}
    out = sv.classify("d", TURNED, runs(), curl)
    assert out["class"] == "C" and out["failed"] == ["(c) closed-loop work 0.5, at 45 degrees"]
    assert sv.classify("d", TURNED, runs(), None)["class"] == "C"  # no loop ratios: not shown to conserve
    assert sv.classify("d", L1A, runs()[:-1])["class"] == "incomplete"


def test_load_records_checks_each_run_against_its_definition(tmp_path):
    law = LawSpec.from_dict(L1A).to_dict()
    definitions = {"definitions": [{"id": "L1a+radial.Re", "law": L1A}]}
    for run_id, angle_value in (("20261114-000001-x", 45.0), ("20261114-000002-x", 45.0)):
        folder = tmp_path / "L1a+radial.Re" / "T1.1" / "phi45" / run_id
        folder.mkdir(parents=True)
        (folder / "stamp.json").write_text(json.dumps({"law": law, "universe": {"phi_bar_deg": angle_value}}))
        summary = record("T1.1", 45.0)["summary"]
        summary["marker"] = run_id
        (folder / "summary.json").write_text(json.dumps(summary))
        (folder / "timeseries.csv").write_text("step,x0_0,v0_0\n0,0.0,1.0\n9,1.5,2.5\n")
    (records,) = sv.load_records(tmp_path, definitions)
    assert records["summary"]["marker"] == "20261114-000002-x"  # the newest run counts
    assert records["final"] == {"step": 9.0, "x0_0": 1.5, "v0_0": 2.5} and records["angle"] == 45.0
    wrong = tmp_path / "L1a+radial.Re" / "T1.1" / "phi45" / "20261114-000003-x"
    wrong.mkdir()
    (wrong / "stamp.json").write_text(json.dumps({"law": dict(law, sign="-"), "universe": {"phi_bar_deg": 45.0}}))
    (wrong / "summary.json").write_text("{}")
    with pytest.raises(ValueError, match="not definition"):
        sv.load_records(tmp_path, definitions)
