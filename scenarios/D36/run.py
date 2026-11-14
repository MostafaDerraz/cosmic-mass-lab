"""D36.1 to D36.4: the wrap-around boundary with the nearest image (decision D36), against its alternatives.

Pre-registration: team/notes/D36.md (committed before this file was written). Every case runs under
the reference law L0, so every statement is mathematics: a check of a numerical method.

A: the nearest image (the engine: `Environment("wrap", ...)` and `DirectSum(period=...)`).
B: positions wrap, forces ignore the copies (written here, not in the engine).
C: the full periodic sum, for a pair with no net charge, as a plain sum over the copies (written here).

    python scenarios/D36/run.py            # all four, writes results/D36.1/ ... results/D36.4/
    python scenarios/D36/run.py 3          # only D36.3
"""

import math
import sys

import numpy as np

from cosmicmass.constants import scaled
from cosmicmass.energy import energy_l0
from cosmicmass.environment import Environment
from cosmicmass.fields import DirectSum, nearest_image
from cosmicmass.gates.t0 import write_gate_result
from cosmicmass.integrators import StepContext, VelocityVerlet
from cosmicmass.runner import run
from cosmicmass.settings import load_settings

PRODUCER = "Claude chat session on claude.ai / Claude Fable 5.2 (sponsor @MostafaDerraz)"
LOWER, UPPER, SIDE = [0.0, 0.0], [1.0, 1.0], 1.0
PERIOD = np.array([SIDE, SIDE])
LATTICE_SUM = 9.0336  # sum of 1/R^3 over the square lattice without the origin: 4 zeta(3/2) beta(3/2)
ETA = 1e-3


def universe(boundary: str) -> dict:
    block = {"units": "scaled", "dimension": 2, "phi_bar_deg": 90.0, "boundary": boundary}
    if boundary != "open":
        block["box"] = {"lower": LOWER, "upper": UPPER}
    return block


def settings(scenario: str, bodies: list, boundary: str, numerics: dict) -> dict:
    return {
        "scenario": scenario,
        "universe": universe(boundary),
        "law": {"id": "L0"},
        "bodies": bodies,
        "numerics": numerics,
        "kind": "mathematics",
        "preregistered": True,
    }


def forces(config: dict, variant: str) -> np.ndarray:
    """The force on every body, evaluated once: variant "A" (nearest image) or "B" (open space)."""
    s = load_settings(config)
    solver = DirectSum(constants=scaled(), softening=s.numerics["softening"], period=PERIOD if variant == "A" else None)
    return s.make_law().forces(s.state, solver.at_bodies(s.state), s.frame).total()


def simulate(config: dict, variant: str) -> dict:
    """Integrate `config` in the wrap-around box with the engine's parts: variant "A" (nearest image
    in the forces, the step rule and the energy) or "B" (positions wrap, everything else as in open
    space); "open" runs the same bodies with no box at all."""
    s = load_settings(config)
    period = PERIOD if variant == "A" else None
    softening = s.numerics["softening"]
    solver = DirectSum(constants=scaled(), softening=softening, period=period)
    environment = Environment("open" if variant == "open" else "wrap", LOWER, UPPER)
    ctx = StepContext(s.make_law(), solver, s.frame, eta=ETA)
    integrator, state = VelocityVerlet(), s.state.copy()
    t_end = math.inf if s.numerics["t_end"] is None else s.numerics["t_end"]
    max_steps = math.inf if s.numerics["steps"] is None else s.numerics["steps"]

    def energy(st):
        return energy_l0(st, s.frame, scaled(), softening, period=period)["total"]

    def momentum(st):
        return np.sum(st.m[:, None] * st.v, axis=0)

    e0, p0 = energy(state), momentum(state)
    # the scale of the energy drift: the interaction energy at the start (change to the plan, 2026-11-14:
    # in D36.3 the kinetic and the binding energy cancel, so the total energy is zero and cannot be the scale)
    u0 = abs(e0 - 0.5 * float(np.sum(state.m * np.sum(state.v**2, axis=1))))
    scale = float(np.sum(state.m * np.linalg.norm(state.v, axis=1)))
    steps, crossings, excursion, separations = 0, 0, 0.0, []
    while steps < max_steps and state.t < t_end:
        dt = integrator.next_dt(state, ctx)
        last = state.t + dt >= t_end
        state = integrator.step(state, min(dt, t_end - state.t), ctx)
        if last:
            state.t = float(t_end)
        steps += 1
        before = state.x.copy()
        environment.apply(state, dt)
        moved = np.any(state.x != before, axis=1)
        if np.any(moved):
            crossings += int(moved.sum())
            ctx.invalidate()
        excursion = max(excursion, abs(energy(state) - e0) / u0)
        if state.n == 2:
            separations.append(float(np.linalg.norm(nearest_image(state.x[1] - state.x[0], PERIOD))))
    return {
        "state": state,
        "steps": steps,
        "crossings": crossings,
        "energy_end_drift": abs(energy(state) - e0) / u0,
        "energy_at_the_start": e0,
        "energy_scale_interaction_at_the_start": u0,
        "energy_largest_excursion": excursion,
        "momentum_drift": float(np.linalg.norm(momentum(state) - p0) / scale) if scale else 0.0,
        "separations": separations,
    }


# ------------------------------------------------------------------ D36.1
def pair_bodies(x_plus: float, x_minus: float) -> list:
    return [
        {"name": "plus", "m": 1e-12, "q": 1.0, "x": [x_plus, 0.5]},
        {"name": "minus", "m": 1e-12, "q": -1.0, "x": [x_minus, 0.5]},
    ]


def full_sum(d: float, shells: int) -> float:
    """C: the attraction on the positive charge (towards the negative one, along x) from the negative
    charge, every copy of it and every copy of the positive charge, in the shells |i|, |j| <= shells.
    Coulomb only (K = 1, charges +1 and -1): Newton's part is 1e-24 of it."""
    i, j = np.meshgrid(np.arange(-shells, shells + 1), np.arange(-shells, shells + 1), indexing="ij")
    shift_x, shift_y = (SIDE * i).ravel(), (SIDE * j).ravel()
    own = (shift_x == 0) & (shift_y == 0)
    # the positive charge at the origin; the negative one at (-d, 0) and its copies pull, the copies of the positive push
    rx, ry = -d + shift_x, shift_y  # from the positive charge to each copy of the negative one
    pull = np.sum(rx / (rx * rx + ry * ry) ** 1.5)  # x component of the attraction towards the negatives
    sx, sy = shift_x[~own], shift_y[~own]  # the copies of the positive charge itself
    push = np.sum(sx / (sx * sx + sy * sy) ** 1.5)  # zero by symmetry; kept for the record
    return float(-(pull - push))  # positive: towards the negative charge (at -x)


def d36_1() -> dict:
    cases, runs = {}, {}
    straddle = settings("D36.1", pair_bodies(0.02, 0.98), "wrap", {"softening": 0.0, "steps": 0})
    fa, fb = forces(straddle, "A"), forces(straddle, "B")
    cases["a: straddling the edge"] = {
        "true_separation": 0.04,
        "expected_attraction": 1 / 0.04**2,
        "A_attraction_on_plus": float(-fa[0, 0]),  # the partner is at -0.04 through the edge
        "B_attraction_on_plus": float(fb[0, 0]),  # B sees it at +0.96
        "A_relative_error": abs(-fa[0, 0] * 0.04**2 - 1),
        "B_fraction_of_the_right_force": float(fb[0, 0] * 0.04**2),
        "third_law_A": float(abs(fa[0] + fa[1]).max() / abs(fa[0]).max()),
    }
    cases["a: straddling the edge"]["pass"] = bool(
        cases["a: straddling the edge"]["A_relative_error"] < 1e-12
        and abs(cases["a: straddling the edge"]["B_fraction_of_the_right_force"] / (0.04 / 0.96) ** 2 - 1) < 1e-9
    )
    runs["a"] = straddle
    for d in (0.05, 0.1, 0.2, 0.4):
        inside = settings("D36.1", pair_bodies(0.5 + d / 2, 0.5 - d / 2), "wrap", {"softening": 0.0, "steps": 0})
        fa, fb = forces(inside, "A"), forces(inside, "B")
        attraction = float(-fa[0, 0])
        sums = {n: full_sum(d, n) for n in (50, 100, 200)}
        excess = (attraction - sums[200]) / attraction
        leading = LATTICE_SUM / 2 * d**3
        case = {
            "d_over_L": d,
            "A_attraction": attraction,
            "A_equals_B": bool(np.array_equal(fa, fb)),
            "C_attraction_shells_50_100_200": [sums[50], sums[100], sums[200]],
            "A_excess_over_C": excess,
            "predicted_leading_excess": leading,
            "excess_over_prediction": excess / leading,
            "C_convergence_100_to_200_of_the_correction": abs(sums[200] - sums[100]) / abs(attraction - sums[200]),
        }
        tolerance = {0.05: 0.2, 0.1: 0.2, 0.2: 1.0}.get(d)  # within 20%, 20%, a factor 2; 0.4: reported only
        judged = (
            True if tolerance is None else (excess > 0 and (1 / (1 + tolerance)) <= excess / leading <= 1 + tolerance)
        )
        case["judged"] = tolerance is not None
        case["pass"] = bool(case["A_equals_B"] and judged and case["C_convergence_100_to_200_of_the_correction"] < 1e-3)
        cases[f"b: inside the box, d = {d}"] = case
        runs[f"b, d = {d}"] = inside
    return {"id": "D36.1", "pass": all(c["pass"] for c in cases.values()), "cases": cases, "runs": runs}


# ------------------------------------------------------------------ the gas of D36.2 and D36.4
def gas_bodies() -> list:
    rng = np.random.default_rng(11)
    x, v = rng.uniform(0.0, SIDE, (40, 2)), rng.normal(0.0, 0.3, (40, 2))
    return [
        {"m": 1.0, "q": 0.05 if k < 20 else -0.05, "x": [float(c) for c in x[k]], "v": [float(c) for c in v[k]]}
        for k in range(40)
    ]


def gas_settings(scenario: str, steps: int = 4000) -> dict:
    return settings(scenario, gas_bodies(), "wrap", {"softening": 0.02, "steps": steps, "eta": ETA, "seed": 11})


def d36_2() -> dict:
    config, cases = gas_settings("D36.2", steps=0), {}
    for variant in ("A", "B"):
        f = forces(config, variant)
        residual = float(np.linalg.norm(f.sum(axis=0)) / np.sum(np.linalg.norm(f, axis=1)))
        cases[variant] = {"third_law_residual": residual, "pass": bool(residual < 1e-12)}
    cases["note"] = {"text": "both are pairwise: the third law does not tell A from B", "pass": True}
    return {"id": "D36.2", "pass": all(c["pass"] for c in cases.values()), "cases": cases, "runs": {"gas": config}}


# ------------------------------------------------------------------ D36.3
def bound_pair(boundary: str) -> dict:
    a, speed, drift = 0.02, 10.0, 5.0
    bodies = [
        {"name": "one", "m": 1.0, "q": 0.0, "x": [0.5 - a / 2, 0.5], "v": [drift, -speed / 2]},
        {"name": "two", "m": 1.0, "q": 0.0, "x": [0.5 + a / 2, 0.5], "v": [drift, speed / 2]},
    ]
    return settings("D36.3", bodies, boundary, {"softening": 0.0, "t_end": 0.35, "eta": ETA})


def d36_3() -> dict:
    a = 0.02
    out = {name: simulate(bound_pair("wrap" if name != "open" else "open"), name) for name in ("A", "B", "open")}
    relative = {k: nearest_image(v["state"].x[1] - v["state"].x[0], PERIOD) for k, v in out.items()}
    relative["open"] = out["open"]["state"].x[1] - out["open"]["state"].x[0]
    cases = {}
    for name, result in out.items():
        seps = np.array(result.pop("separations"))
        result.pop("state")
        result["separation_largest_change"] = float(np.max(np.abs(seps / a - 1)))
        result["separation_at_the_end_over_a"] = float(np.linalg.norm(relative[name]) / a)
        cases[name] = result
    cases["A"]["relative_position_against_open_over_a"] = float(np.linalg.norm(relative["A"] - relative["open"]) / a)
    cases["A"]["steps_against_open"] = cases["A"]["steps"] / cases["open"]["steps"]
    flags = {
        "A: the separation stays a within 1e-5": cases["A"]["separation_largest_change"] < 1e-5,
        "A: energy end drift below 1e-6": cases["A"]["energy_end_drift"] < 1e-6,
        "A: energy largest excursion below 1e-6": cases["A"]["energy_largest_excursion"] < 1e-6,
        "A: final relative position within 1e-6 a of the open-space run": cases["A"][
            "relative_position_against_open_over_a"
        ]
        < 1e-6,
        "A: steps within 1% of the open-space run": abs(cases["A"]["steps_against_open"] - 1) < 0.01,
        "A: the pair crossed the edge (4 body crossings)": cases["A"]["crossings"] == 4,
        "B: the pair is torn apart (separation above 10 a at the end)": cases["B"]["separation_at_the_end_over_a"] > 10,
    }
    cases["predictions"] = {"flags": {k: bool(v) for k, v in flags.items()}, "pass": all(flags.values())}
    for case in cases.values():
        case.setdefault("pass", True)
    return {"id": "D36.3", "pass": cases["predictions"]["pass"], "cases": cases, "runs": {"pair": bound_pair("wrap")}}


# ------------------------------------------------------------------ D36.4
def d36_4() -> dict:
    config = gas_settings("D36.4")
    stamped = run(config, out_root=None)  # the engine's runner, under A
    loops = {variant: simulate(config, variant) for variant in ("A", "B")}
    cases = {}
    for variant, result in loops.items():
        state = result.pop("state")
        result.pop("separations")
        result["inside_the_box"] = bool(np.all((state.x >= 0.0) & (state.x < SIDE)))
        result["final_x"] = state.x
        cases[f"{variant} (the script's loop)"] = result
    a, b = cases["A (the script's loop)"], cases["B (the script's loop)"]
    agreement = float(np.max(np.abs(nearest_image(a.pop("final_x") - stamped.state.x, PERIOD))))
    b.pop("final_x")
    conservation = stamped.summary["conservation"]
    cases["A (the engine's runner)"] = {
        "steps": stamped.summary["steps"],
        "mass_drift": conservation["mass"],
        "charge_drift": conservation["charge"],
        "energy_end_drift": conservation["energy_re"],
        "momentum_drift": conservation["momentum"],
        "inside_the_box": bool(np.all((stamped.state.x >= 0.0) & (stamped.state.x < SIDE))),
        "wrap_around": stamped.summary["wrap_around"],
        "agreement_with_the_loop": agreement,
    }
    runner = cases["A (the engine's runner)"]
    flags = {
        "A: mass and charge exactly constant": runner["mass_drift"] == 0.0 and runner["charge_drift"] == 0.0,
        "A: every body inside the box at the end": runner["inside_the_box"] and a["inside_the_box"],
        "A: momentum drift below 1e-12": runner["momentum_drift"] < 1e-12 and a["momentum_drift"] < 1e-12,
        "A: at least 20 edge crossings": a["crossings"] >= 20,
        "A: energy end drift below 1e-3": runner["energy_end_drift"] < 1e-3 and a["energy_end_drift"] < 1e-3,
        "B: energy end drift at least 10 times that of A": b["energy_end_drift"] >= 10 * a["energy_end_drift"],
        "the loop and the runner agree for A to 1e-9": agreement < 1e-9,
        "the summary carries the stated limit of D36": bool(
            runner["wrap_around"] and "half the box" in runner["wrap_around"]["limit"]
        ),
    }
    cases["predictions"] = {"flags": {k: bool(v) for k, v in flags.items()}, "pass": all(flags.values())}
    for case in cases.values():
        case.setdefault("pass", True)
    return {"id": "D36.4", "pass": cases["predictions"]["pass"], "cases": cases, "runs": {"gas": stamped}}


TESTS = {"1": d36_1, "2": d36_2, "3": d36_3, "4": d36_4}

if __name__ == "__main__":
    chosen = [a for a in sys.argv[1:] if a in TESTS] or list(TESTS)
    failed = False
    for key in chosen:
        outcome = TESTS[key]()
        folder = write_gate_result(
            outcome,
            producer=PRODUCER,
            kind="mathematics",
            law="the reference law L0 (a test of the wrap-around boundary)",
        )
        for name, case in outcome["cases"].items():
            if "flags" in case:
                for flag, ok in case["flags"].items():
                    print(f"  {'holds' if ok else 'FAILS'}  {flag}")
            else:
                shown = {
                    k: (f"{v:.4g}" if isinstance(v, float) else v)
                    for k, v in case.items()
                    if not isinstance(v, (dict, list))
                }
                print(f"  {name}: {shown}")
        print(f"{outcome['id']}: {'pass' if outcome['pass'] else 'FAIL'}; results in {folder}")
        failed = failed or not outcome["pass"]
    sys.exit(1 if failed else 0)
