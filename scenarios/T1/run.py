"""T1.1 to T1.6: pairs of cosmic masses under L0 and the default exploratory law (phase 1 of M2).

Pre-registration, predictions and analysis plan: team/notes/T1.md. L0 is the reference law; its runs
are judged against the two-body problem (mathematics). L1a is a hypothesis (decision D1), excluded at
lambda = 1 by measurement (T9.1): its runs are findings of the virtual universe, compared with the
predictions written before any run.

    python scenarios/T1/run.py                 # the six tests; writes results/T1.x/<run-id>/
    python scenarios/T1/run.py T1.3 T1.5       # some of them
    python scenarios/T1/run.py T1.4 --law '{"id": "L1b"}'    # one more law: measured, no prediction
    python scenarios/T1/run.py T1.3 --dry-run  # runs, prints, writes nothing
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pairs  # noqa: E402
from cosmicmass.constants import SI  # noqa: E402
from cosmicmass.fields import DirectSum  # noqa: E402
from cosmicmass.laws import third_law_residual  # noqa: E402
from cosmicmass.settings import load_settings  # noqa: E402
from cosmicmass.stamp import ResultsError, engine_commit, make_run_id, validate_outputs  # noqa: E402

PRODUCER = "the OpenCode agent / Claude Fable 5.2, local (sponsor @MostafaDerraz)"
LAWS = ("reference", "exploratory")
L0, L1A = "L0", "L1a"

# tolerances of the analysis plan (team/notes/T1.md; team/RESULTS-STANDARD.md, section 5)
TIME, CLOSURE, DRIFT, THIRD_LAW, STILL, CONVERGED = 1e-4, 1e-3, 1e-6, 1e-12, 1e-9, 1e-2

QUESTIONS = {
    "T1.1": "Does the Composer part stay silent for pure masses?",
    "T1.2": "`s q g`: pull, push or swirl? Does the orbit close?",
    "T1.3": "How does the charge share `-2 cos(phi_bar) q E` change the electric force?",
    "T1.4": "Bound? Precession? Spiral? Binding energy?",
    "T1.5": "No interaction under L0. What under L1 to L8?",
    "T1.6": "Newton's third law: force on 1 = -(force on 2)?",
}

# The predictions of the pre-registration, with the digits written there. `f = k / k_ref`: positive is
# a pull. `t`: the time of the event (s); `period`: the time of one revolution (s).
_ORBIT = {"event": "three revolutions", "closes": True}
PREDICTED = {
    ("T1.1", L0): {**_ORBIT, "f": 1.0, "period": 2.357430e6},
    ("T1.2a", L0): {"event": "contact", "f": 1.0, "t": 387.2847},
    ("T1.2b", L0): {**_ORBIT, "f": 1.0, "period": 5.828599e3},
    ("T1.2c", L0): {"event": "contact", "f": 1.0, "t": 387.2847},
    ("T1.2d", L0): {**_ORBIT, "f": 1.0, "period": 5.828599e3},
    ("T1.3a", L0): {"event": "escape", "f": -1.0, "t": 1.522016e-14},
    ("T1.3b", L0): {"event": "plunge", "f": 1.0, "t": 6.880887e-17},
    ("T1.4", L0): {**_ORBIT, "f": 1.0, "period": 1.519415e-16, "binding_eV": 13.60570},
    ("T1.5", L0): {"event": "no event", "f": 0.0, "still": True},
    ("T1.6", L0): {**_ORBIT, "f": 1.0, "period": 2.127206e-16, "inner_apsis": 0.324503},
    ("T1.1", L1A): {**_ORBIT, "f": 1.0, "period": 2.357430e6, "same_as_L0": True},
    ("T1.2a", L1A): {"event": "contact", "f": 1.111554e18, "t": 3.673374e-7},
    ("T1.2b", L1A): {"event": "contact", "f": 1.111554e18, "t": 3.6734e-7, "closes": False},
    ("T1.2c", L1A): {"event": "escape", "f": -2.040982e21, "t": 1.641463e-7},
    ("T1.2d", L1A): {"event": "escape", "f": -2.040982e21, "t": 1.6415e-7, "closes": False},
    ("T1.3a", L1A): {"event": "escape", "f": -0.647892008, "t": 1.890897e-14},
    ("T1.3b", L1A): {"event": "plunge", "f": 0.647892008, "t": 8.548563e-17},
    ("T1.4", L1A): {
        **_ORBIT,
        "f": 0.647892008,
        "period": 6.119516e-16,
        "outer_apsis": 3.380845,
        "binding_eV": 4.024348,
    },
    ("T1.5", L1A): {"event": "plunge", "f": 1.176053996, "t": 1.433635e9},
    ("T1.6", L1A): {**_ORBIT, "f": 0.647892008, "period": 3.535491e-16, "inner_apsis": 0.608103},
}


def _near(value, reference, tolerance) -> bool:
    return bool(abs(value - reference) <= tolerance * abs(reference))


def _drift(measured: dict) -> float:
    """The relative drift of the energy the law conserves; where the energy starts at zero and the
    runner reports no relative drift, its change in units of `k_ref / r0`."""
    drift = measured["conservation"][measured["energy_used"]]
    return measured["energy_change_over_reference"] if drift is None else drift


def converge(series: list) -> dict:
    """The convergence block: the same quantities at eta, eta/2 and eta/4, and their changes."""
    full, half, quarter = series
    out = {
        "eta": [m["eta"] for m in series],
        "event": [m["event"] for m in series],
        "steps": [m["steps"] for m in series],
        "t_event": [m["t_event"] for m in series],
        "separation_end_over_r0": [m["separation_end_over_r0"] for m in series],
        "energy_drift": [_drift(m) for m in series],
        "energy_drift_below_1e-6_at_eta": bool(_drift(full) < DRIFT),  # reported; judged at eta/4
        "momentum_drift": [m["conservation"]["momentum"] for m in series],
        "largest_closure_error": [max((r["closure_error"] for r in m["revolutions"]), default=None) for m in series],
    }
    changes = {}
    for key in ("t_event", "separation_end_over_r0"):
        a, b, c = out[key]
        changes[key] = {
            "eta_to_half": abs(a - b) / abs(b) if b else abs(a - b),
            "half_to_quarter": abs(b - c) / abs(c) if c else abs(b - c),
        }
    out["changes"] = changes
    drifts = out["energy_drift"]
    out["energy_drift_ratio_per_halving"] = [
        drifts[i] / drifts[i + 1] if drifts[i + 1] else None for i in range(2)
    ]  # 4 for a second-order method
    out["converged"] = bool(
        full["event"] == half["event"] == quarter["event"]
        and all(change["half_to_quarter"] < CONVERGED for change in changes.values())
    )
    return out


def judge(case: str, label: str, series: list, reference: dict | None = None) -> dict:
    """The flags of the analysis plan for one case under one law: `{flag: True or False}`. `series`:
    the measurements at eta, eta/2 and eta/4; `reference`: the L0 measurement of the same case (T1.1)."""
    measured, predicted = series[0], PREDICTED.get((case, label))
    conservation = measured["conservation"]
    flags = {
        "mass and charge exact": conservation["mass"] == 0 and conservation["charge"] == 0,
        "momentum drift below 1e-6": conservation["momentum"] < DRIFT,
        "energy drift below 1e-6 at eta/4": _drift(series[2]) < DRIFT,
        "third law: residual below 1e-12": all(
            residual is None or residual < THIRD_LAW
            for residual in (measured["third_law_residual_start"], measured["third_law_residual_end"])
        ),
    }
    if predicted is None:
        return {key: bool(value) for key, value in flags.items()}
    first, turns = measured["first_acceleration"], measured["revolutions"]
    flags["event as predicted"] = measured["event"] == predicted["event"]
    if predicted["f"] == 0:
        flags["first acceleration: none"] = first["what"] == "none"
    else:
        flags["first acceleration as predicted"] = _near(-first["radial"], predicted["f"], TIME) and (
            "swirl" not in first["what"]
        )
    if "t" in predicted:
        flags["time of the event within 1e-4"] = _near(measured["t_event"], predicted["t"], TIME)
    if "period" in predicted:
        flags["time of each revolution within 1e-4"] = len(turns) == pairs.REVOLUTIONS and all(
            _near(turn["t"], turn["revolution"] * predicted["period"], TIME) for turn in turns
        )
    if predicted.get("closes") is True:
        flags["the orbit closes (closure error below 1e-3)"] = len(turns) == pairs.REVOLUTIONS and all(
            turn["closure_error"] < CLOSURE for turn in turns
        )
    if predicted.get("closes") is False:
        flags["the orbit does not close (no revolution)"] = len(turns) == 0
    for kind in ("inner", "outer"):
        if f"{kind}_apsis" in predicted:
            found = measured["apsides"].get(f"{kind}_over_r0") or []
            others = [r for r in found if not _near(r, 1.0, TIME)]  # the start is the other apsis
            flags[f"{kind} apsis within 1e-4"] = bool(others) and all(
                _near(r, predicted[f"{kind}_apsis"], TIME) for r in others
            )
    if "binding_eV" in predicted:
        flags["binding energy within 1e-4"] = _near(measured["binding_energy_start_eV"], predicted["binding_eV"], TIME)
    if predicted.get("still"):
        flags["the separation changes by less than 1e-9"] = (
            abs(measured["separation_max_over_r0"] - 1) < STILL and abs(measured["separation_min_over_r0"] - 1) < STILL
        )
    if predicted.get("same_as_L0") and reference is not None:
        flags["identical to L0 (same steps, same final state)"] = (
            measured["steps"] == reference["steps"] and measured["final_state"] == reference["final_state"]
        )
    return {key: bool(value) for key, value in flags.items()}


def third_law_at_start(case: str, law) -> float | None:
    """T1.6: the third-law residual of a case at its first state (no run)."""
    loaded = load_settings(pairs.settings(case, law))
    solver = DirectSum(constants=SI, softening=0.0)
    return third_law_residual(loaded.make_law(), loaded.state, solver, loaded.frame)["relative"]


def answers(test: str, cases: dict, laws: dict) -> dict:
    """The answer to the test's question under each law, from the measured numbers only:
    yes, no, impossible or unclear, with the numbers behind it. `laws`: `{label: law}`."""
    out = {}
    for label, law in laws.items():
        got = {case: cases[f"{case}; {label}"] for case in pairs.TESTS[test]}
        if test == "T1.1":
            one, base = got["T1.1"], cases[f"T1.1; {L0}"]
            same = one["steps"] == base["steps"] and one["final_state"] == base["final_state"]
            out[label] = {"answer": "yes" if same else "no", "identical_to_L0_bit_for_bit": same}
        elif test == "T1.2":
            closes = {
                case: bool(got[case]["revolutions"])
                and all(t["closure_error"] < CLOSURE for t in got[case]["revolutions"])
                for case in ("T1.2b", "T1.2d")
            }
            out[label] = {
                "proton": got["T1.2a"]["first_acceleration"],
                "electron": got["T1.2c"]["first_acceleration"],
                "the orbit closes": {case: "yes" if ok else "no" for case, ok in closes.items()},
                "events": {case: got[case]["event"] for case in got},
            }
        elif test == "T1.3":
            ratio = {
                case: got[case]["first_acceleration"]["radial"] / cases[f"{case}; {L0}"]["first_acceleration"]["radial"]
                for case in got
            }
            out[label] = {"answer": "the electric force is multiplied by", "factor": ratio}
        elif test == "T1.4":
            one = got["T1.4"]
            turns, apsides = one["revolutions"], one["apsides"]
            bound = one["event"] == "three revolutions"
            closes = bool(turns) and all(t["closure_error"] < CLOSURE for t in turns)
            out[label] = {
                "bound": "yes" if bound else ("no" if one["event"] == "escape" else "unclear"),
                "precession": "no" if closes else ("yes" if apsides.get("largest_precession_rad") else "unclear"),
                "spiral": "no" if closes else "unclear",
                "binding_energy_eV": one["binding_energy_start_eV"],
                "largest_closure_error": max((t["closure_error"] for t in turns), default=None),
                "largest_precession_rad": apsides.get("largest_precession_rad"),
            }
        elif test == "T1.5":
            one = got["T1.5"]
            out[label] = {
                "answer": "no interaction" if one["first_acceleration"]["what"] == "none" else "interaction",
                "first_acceleration": one["first_acceleration"],
                "event": one["event"],
            }
        else:
            start = {case: third_law_at_start(case, law) for case in pairs.CASES}
            run_values = [got["T1.6"]["third_law_residual_start"], got["T1.6"]["third_law_residual_end"]]
            worst = max((v for v in list(start.values()) + run_values if v is not None), default=None)
            out[label] = {
                "answer": "yes" if worst is None or worst < THIRD_LAW else "no",
                "largest_residual": worst,
                "residual_at_the_start_of_every_case": start,
                "momentum_drift_of_the_run": got["T1.6"]["conservation"]["momentum"],
            }
    return out


def run_test(test: str, laws=LAWS, producer: str = PRODUCER, report=print, eta: float = pairs.ETA, only=None) -> dict:
    """One test under each law, with its convergence series. Returns the outcome for `write_test`.
    `eta` and `only` (some of the test's cases) are for the tests of this script: a production run
    uses the pre-registered step and every case."""
    cases, runs, labels, failed_l0 = {}, {}, {}, False
    for law in laws:
        label = pairs.law_label(law)
        labels[label] = law
        for case in only or pairs.TESTS[test]:
            series, first = [], None
            for step in (eta, eta / 2, eta / 4):
                one = pairs.run_case(case, law, eta=step, producer=producer)
                series.append(one.measured)
                first = one if first is None else first
            name = f"{case}; {label}"
            measured = first.measured
            measured["final_state"] = {"x": first.result.state.x.tolist(), "v": first.result.state.v.tolist()}
            measured["convergence"] = converge(series)
            measured["flags"] = judge(case, label, series, cases.get(f"{case}; {L0}"))
            measured["holds"] = all(measured["flags"].values())
            measured["predicted"] = PREDICTED.get((case, label))
            measured["settings"] = first.settings  # with the plain runner, they give this run again
            measured["pass"] = measured["holds"]
            cases[name], runs[name] = measured, first.result
            failed_l0 = failed_l0 or (label == L0 and not measured["holds"])
            report(f"{name}: {measured['event']} at t = {measured['t_event']:.6e} s; " + _failures(measured))
    predicted_fail = [name for name, case in cases.items() if not case["holds"]]
    return {
        "id": test,
        "question": QUESTIONS[test],
        "verdict": "fail" if failed_l0 else "finding",
        "cases": cases,
        "runs": runs,
        "answers": answers(test, cases, labels) if only is None else {},
        "not_holding": predicted_fail,
    }


def _failures(measured: dict) -> str:
    bad = [flag for flag, ok in measured["flags"].items() if not ok]
    return "every flag holds" if not bad else "FAILS: " + "; ".join(bad)


def write_test(outcome: dict, out_root="results", producer: str = PRODUCER, now=None) -> Path:
    """Write one test to `results/<ID>/<run-id>/`: `stamp.json` (the runner's stamp of every case),
    `summary.json` (the measurements, flags, convergence and answers) and `timeseries-<case>.csv`."""
    commit = engine_commit()
    run_id = make_run_id(commit, now)
    stamps = {name: dict(result.stamp, run_id=run_id, producer=producer) for name, result in outcome["runs"].items()}
    stamp = {
        "scenario": outcome["id"],
        "run_id": run_id,
        "engine_commit": commit,
        "kind": "finding",
        "preregistered": True,
        "producer": producer,
        "cases": stamps,
    }
    bad = outcome["not_holding"]
    summary = {
        "scenario": outcome["id"],
        "verdict": outcome["verdict"],
        "question": outcome["question"],
        "answers": outcome["answers"],
        "cases": outcome["cases"],
        "prediction": "team/notes/T1.md, 'Predictions' (committed before any run); each case carries its own",
        "explanation": f"{outcome['id']}, {outcome['question']} "
        + (
            "Every case meets its pre-registered flags."
            if not bad
            else "Cases with a flag that does not hold: " + ", ".join(bad) + "."
        )
        + " L0 is mathematics; L1a is a hypothesis (D1), excluded at lambda = 1 by measurement (T9.1).",
    }
    problems = validate_outputs(stamp, summary)
    if problems:
        raise ResultsError("the outputs break the results standard: " + "; ".join(problems))
    folder = Path(out_root) / outcome["id"] / run_id
    folder.mkdir(parents=True, exist_ok=False)
    for name, data in (("stamp.json", stamp), ("summary.json", summary)):
        (folder / name).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for name, result in outcome["runs"].items():
        safe = "".join(c if c.isalnum() or c in ".-=," else "_" for c in name)
        with open(folder / f"timeseries-{safe}.csv", "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(result.timeseries[0]))
            writer.writeheader()
            writer.writerows(result.timeseries)
    return folder


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="T1.1 to T1.6 under L0 and the default exploratory law.")
    parser.add_argument("tests", nargs="*", default=list(pairs.TESTS), help="T1.1 ... T1.6 (default: all)")
    parser.add_argument("--law", action="append", default=[], help="one more law block, as JSON (no prediction)")
    parser.add_argument("--producer", default=PRODUCER)
    parser.add_argument("--out", default="results")
    parser.add_argument("--dry-run", action="store_true", help="run and print, write nothing")
    args = parser.parse_args(argv)
    unknown = [test for test in args.tests if test not in pairs.TESTS]
    if unknown:
        parser.error(f"unknown test(s) {unknown}; known: {list(pairs.TESTS)}")
    laws = LAWS + tuple(json.loads(block) for block in args.law)
    failed = False
    for test in args.tests:
        outcome = run_test(test, laws, args.producer)
        failed = failed or outcome["verdict"] == "fail"
        print(f"{test}: {outcome['verdict']}; answers: {json.dumps(outcome['answers'], default=_plain)}")
        if not args.dry_run:
            print(f"results in {write_test(outcome, args.out, args.producer)}")
    return 1 if failed else 0


def _plain(value):
    return value.tolist() if isinstance(value, np.ndarray) else str(value)


if __name__ == "__main__":
    sys.exit(main())
