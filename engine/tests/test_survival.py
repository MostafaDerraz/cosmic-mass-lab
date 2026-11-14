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


def write_run(folder, law, scenario, angle, marker):
    folder.mkdir(parents=True)
    stamp = {"law": law, "scenario": scenario, "universe": {"phi_bar_deg": angle}}
    (folder / "stamp.json").write_text(json.dumps(stamp))
    summary = record("T1.1", angle)["summary"]
    summary["marker"] = marker
    (folder / "summary.json").write_text(json.dumps(summary))
    (folder / "timeseries.csv").write_text("step,x0_0,v0_0\n0,0.0,1.0\n9,1.5,2.5\n")


def test_load_records_reads_runs_by_their_stamps_whatever_the_layout(tmp_path):
    law = LawSpec.from_dict(L1A).to_dict()
    other = LawSpec.from_dict(TURNED).to_dict()
    definitions = {"definitions": [{"id": "L1a+radial.Re", "law": L1A}, {"id": "L1a+perpendicular.Re", "law": TURNED}]}
    # the layout of T2.2 (#95): results/<scene>/<definition>/<run-id>/, scenario "<scene>/<definition>"
    write_run(tmp_path / "T1.1@45" / "L1a+radial.Re" / "20261114-000001-x", law, "T1.1@45/L1a+radial.Re", 45.0, "old")
    write_run(tmp_path / "T1.1@45" / "L1a+radial.Re" / "20261114-000002-x", law, "T1.1@45/L1a+radial.Re", 45.0, "new")
    # the layout first pre-registered: results/T2/<definition>/<setup>/phi<angle>/<run-id>/
    write_run(
        tmp_path / "T2" / "L1a+perpendicular.Re" / "T1.1" / "phi90" / "20261114-000003-x", other, "T1.1", 90.0, "b"
    )
    # runs that are not T2's are left out: another law, another setup, an angle that is not one of the six
    write_run(tmp_path / "T0.1" / "20261110-000000-x", {"id": "L0"}, "T0.1", 79.86, "gate")
    write_run(tmp_path / "x" / "20261114-000004-x", law, "Q5", 45.0, "other setup")
    write_run(tmp_path / "y" / "20261114-000005-x", law, "T1.1", 30.0, "other angle")
    records = sv.load_records(tmp_path, definitions)
    assert [(r["definition"], r["setup"], r["angle"], r["summary"]["marker"]) for r in records] == [
        ("L1a+perpendicular.Re", "T1.1", 90.0, "b"),
        ("L1a+radial.Re", "T1.1", 45.0, "new"),  # the newest run counts
    ]
    assert records[1]["final"] == {"step": 9.0, "x0_0": 1.5, "v0_0": 2.5}
    # a folder that names a definition its stamp's law is not: refused
    write_run(tmp_path / "T1.2" / "L1a+radial.Re" / "20261114-000006-x", other, "T1.2/L1a+radial.Re", 45.0, "bad")
    with pytest.raises(ValueError, match="the folder names"):
        sv.load_records(tmp_path, definitions)
