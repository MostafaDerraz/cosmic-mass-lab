"""T2.3: the closed-loop work test of every coupling definition, at the six angles of D44 (milestone M2).

The test and its predictions are pre-registered in `team/notes/T2.md` ("T2.3: the classifier"), committed
before this script was written or run. It writes `results/T2.3/<run-id>/` (one case per definition: its
ratio at each angle, whether it passes the 1e-9 of D43, and the prediction), through the gate writer, so
that the stamp of every case is valid.

The verdict is about the **predictions**: "pass" if every definition behaved as predicted, "fail" if one
did not. Which definitions pass the loop test is the result, whatever the verdict.

    python scenarios/T2/loop_test.py
"""

import json
import sys
from pathlib import Path

from cosmicmass.analysis.survival import ANGLES, LOOP_TOLERANCE, SOURCES, TEST_BODY, closed_loop_ratio
from cosmicmass.gates.t0 import write_gate_result

ROOT = Path(__file__).resolve().parents[2]
PRODUCER = "Claude Opus 5.5 in Claude Code, local (sponsor @MostafaDerraz)"


def predicted_pass(law: dict) -> bool:
    """The pre-registered prediction: radial L1a, L1b, L1c; every torque; radial L1e pass; the rest fail."""
    reading, direction = law["id"], law["direction"]
    if direction == "torque":
        return True
    if direction == "radial":
        return reading in ("L1a", "L1b", "L1c", "L1e")
    return False


def settings_of(definition: dict) -> dict:
    """The settings of the test, for the stamp: the two fixed sources and the test body (scaled units)."""
    bodies = [{"m": m, "q": q, "x": list(place), "fixed": True} for m, q, place in SOURCES]
    bodies.append({"m": TEST_BODY[0], "q": TEST_BODY[1], "x": [1.0, 1.5]})
    return {
        "scenario": "T2.3",
        "universe": {"units": "scaled", "dimension": 2, "phi_bar_deg": 79.86},
        "law": definition["law"],
        "bodies": bodies,
        "numerics": {"steps": 1, "softening": 0.0},
        "kind": "mathematics",
        "preregistered": True,
    }


def main() -> dict:
    definitions = json.loads((ROOT / "scenarios" / "T2" / "definitions.json").read_text(encoding="utf-8"))
    cases, runs = {}, {}
    for definition in definitions["definitions"]:
        ratios = {f"{a:g}": closed_loop_ratio(definition["law"], a) for a in ANGLES}
        passed = all(r <= LOOP_TOLERANCE for r in ratios.values())
        expected = predicted_pass(definition["law"])
        cases[definition["id"]] = {
            "pass": passed == expected,  # the case "passes" when the prediction held
            "loop_test": "pass" if passed else "fail",
            "predicted": "pass" if expected else "fail",
            "largest_ratio": max(ratios.values()),
            "ratio_by_angle": ratios,
        }
        runs[definition["id"]] = settings_of(definition)
    held = sum(c["pass"] for c in cases.values())
    return {"id": "T2.3", "pass": held == len(cases), "cases": cases, "runs": runs, "held": held}


if __name__ == "__main__":
    outcome = main()
    folder = write_gate_result(outcome, producer=PRODUCER, kind="mathematics",
                               law="the 102 coupling definitions of T2.1 (hypotheses, D1)")  # fmt: skip
    cases = outcome["cases"]
    passing = sorted(k for k, c in cases.items() if c["loop_test"] == "pass")
    print(f"loop test: {len(passing)} of {len(cases)} definitions pass (predicted: 48)")
    print(f"predictions held: {outcome['held']} of {len(cases)}")
    for key, case in cases.items():
        if not case["pass"]:
            print(f"  not as predicted: {key}: {case['loop_test']} (predicted {case['predicted']}), "
                  f"largest ratio {case['largest_ratio']:.3e}")  # fmt: skip
    print(f"results in {folder}")
    sys.exit(0)
