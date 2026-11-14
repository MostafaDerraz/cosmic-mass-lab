"""T1: pairs of cosmic masses (T1.1 to T1.6): the cases, the protocol and the measurements.

Setups and questions: `docs/05-scenarios-and-tests.md`, section 2, T1. Pre-registration:
`team/notes/T1.md`. Decisions: D46 (the bodies carry their physical radii), D49 (under L1e the
measurements use the Real position `x`), D21 (two-body cases run with softening 0).

One protocol for every law, so that the batch runner of T2 can run any `LawSpec` on the same cases:

- the reference scales of a case (`scales`) are built from `k_ref = G_N m_a m_b + K |q_a q_b|`, not
  from the law under test;
- a run ends at an **event** (`Watch`): contact, plunge, escape, three revolutions, or a time cap;
- the engine's runner is used unchanged: a first pass with every step recorded finds the step `n`
  of the event, and the result is the plain run of the same settings with `numerics.steps = n`
  (`run_case`). The two passes must end in the same state, bit for bit.

This file is scenario code, outside `engine/cosmicmass/`: the engine and its `t0_status` are untouched.

    from pairs import CASES, settings, run_case      # with scenarios/T1 on sys.path
"""

import math
from dataclasses import dataclass, field

import numpy as np

from cosmicmass.constants import SI
from cosmicmass.fields import DirectSum
from cosmicmass.integrators import StepContext
from cosmicmass.laws import third_law_residual
from cosmicmass.matter import CosmicMass, preset
from cosmicmass.runner import Interrupted, RunResult, run
from cosmicmass.settings import NAMED_LAWS, load_settings

__all__ = [
    "CASES",
    "TESTS",
    "RADII",
    "PHI",
    "ETA",
    "STEP_CAP",
    "bodies",
    "scales",
    "settings",
    "Watch",
    "find_event",
    "run_case",
    "first_acceleration",
    "law_label",
]

PHI = 79.86  # degrees: the default frame angle
ETA = 1e-3  # D10
STEP_CAP = 60_000  # steps at eta = 1e-3; in proportion at a smaller eta
TIME_CAP = 20.0  # reference periods
DT_MAX = 1.0 / 500.0  # of a reference period: only so that a pair with no force has a step
PLUNGE, ESCAPE, REVOLUTIONS = 0.1, 10.0, 3
SILENT = 1e-12  # a first acceleration below this, in units of a_ref, is "none"

# D46: Earth, moon and proton have radii (m); the electron stays a point. The Stoney body is a point
# (a labelled choice of the pre-registration: D46 names no radius for it).
RADII = {"Earth": 6.371e6, "moon": 1.7374e6, "proton": 0.84e-15, "electron": 0.0, "Stoney body": 0.0}
STONEY = CosmicMass(SI.stoney_mass, SI.e, "Stoney body")  # m_S + e G: Re(mu mu) = 0 to rounding

# case: (first body, second body, separation r0 in m, relative speed in units of v_ref)
CASES = {
    "T1.1": ("Earth", "moon", 3.844e8, 1.0),
    "T1.2a": ("Earth", "proton", 7e6, 0.0),
    "T1.2b": ("Earth", "proton", 7e6, 1.0),
    "T1.2c": ("Earth", "electron", 7e6, 0.0),
    "T1.2d": ("Earth", "electron", 7e6, 1.0),
    "T1.3a": ("proton", "proton", 1e-10, 0.0),
    "T1.3b": ("proton", "electron", 1e-10, 0.0),
    "T1.4": ("proton", "electron", 5.29177e-11, 1.0),
    "T1.5": ("Stoney body", "Stoney body", 1.0, 0.0),
    "T1.6": ("proton", "electron", 1e-10, 0.7),
}
TESTS = {
    "T1.1": ("T1.1",),
    "T1.2": ("T1.2a", "T1.2b", "T1.2c", "T1.2d"),
    "T1.3": ("T1.3a", "T1.3b"),
    "T1.4": ("T1.4",),
    "T1.5": ("T1.5",),
    "T1.6": ("T1.6",),
}


def _body(name: str) -> CosmicMass:
    return STONEY if name == STONEY.name else preset(name)


def bodies(case: str) -> tuple:
    """The two cosmic masses of a case."""
    first, second, _, _ = CASES[case]
    return _body(first), _body(second)


def scales(case: str) -> dict:
    """The reference scales of a case, the same for every law: `k_ref = G_N m_a m_b + K |q_a q_b|`,
    the reduced mass, `a_ref`, `v_ref`, `T_ref`, and the separation at which a run ends as "contact"
    (`R_a + R_b`, D46) or "plunge" (`r0 / 10`), whichever is larger."""
    a, b = bodies(case)
    r0, beta = CASES[case][2], CASES[case][3]
    k_ref = SI.G_N * a.m * b.m + SI.K * abs(a.q * b.q)
    reduced = a.m * b.m / (a.m + b.m)
    touch = RADII[a.name] + RADII[b.name]
    return {
        "r0": r0,
        "beta": beta,
        "k_ref": k_ref,
        "reduced_mass": reduced,
        "a_ref": k_ref / (reduced * r0**2),
        "v_ref": math.sqrt(k_ref / (reduced * r0)),
        "T_ref": 2 * math.pi * math.sqrt(reduced * r0**3 / k_ref),
        "close": max(touch, PLUNGE * r0),
        "close_event": "contact" if touch > PLUNGE * r0 else "plunge",
    }


def _law_block(law) -> dict:
    """A `law` block from a block, a name ("reference", "exploratory") or an entry of
    `scenarios/T2/definitions.json` (`{"id": ..., "law": {...}}`)."""
    if isinstance(law, str):
        return dict(NAMED_LAWS[law])
    return dict(law["law"]) if isinstance(law.get("law"), dict) else dict(law)


def law_label(law) -> str:
    """A short name of a law: its reading and every option that is not the default exploratory one."""
    spec = load_settings({"scenario": "-", "law": _law_block(law), "bodies": [{"m": 1.0, "q": 0.0, "x": [0.0, 0.0]}]})
    block = spec.law_spec.to_dict()
    if block["id"] == "L0":
        return "L0"
    extra = [
        f"{key} {value}" for key, value in block.items() if key != "id" and value != NAMED_LAWS["exploratory"][key]
    ]
    return block["id"] + ("" if not extra else " (" + ", ".join(extra) + ")")


def settings(
    case: str,
    law="reference",
    phi_bar_deg: float = PHI,
    eta: float = ETA,
    steps: int | None = None,
    record_every: int = 1,
    producer: str = "",
    scenario: str | None = None,
) -> dict:
    """The settings of a case for `runner.run`. `steps=None`: the step cap of the first pass; `steps=n`:
    the plain run that ends at step `n`. The bodies start on the x axis in their centre-of-mass frame,
    the first at negative x, with the relative velocity `beta v_ref` along +y."""
    a, b = bodies(case)
    sc = scales(case)
    total, r0, speed = a.m + b.m, sc["r0"], sc["beta"] * sc["v_ref"]
    block = _law_block(law)
    listed = [
        {"name": a.name, "m": a.m, "q": a.q, "x": [-b.m / total * r0, 0.0], "v": [0.0, -b.m / total * speed]},
        {"name": b.name, "m": b.m, "q": b.q, "x": [a.m / total * r0, 0.0], "v": [0.0, a.m / total * speed]},
    ]
    for entry in listed:
        if RADII[entry["name"]] > 0:
            entry["radius"] = RADII[entry["name"]]
    return {
        "scenario": case if scenario is None else scenario,
        "universe": {"phi_bar_deg": phi_bar_deg, "dimension": 2, "units": "SI"},
        "law": block,
        "bodies": listed,
        "numerics": {
            "softening": 0.0,
            "eta": eta,
            "dt_max": DT_MAX * sc["T_ref"],
            "t_end": TIME_CAP * sc["T_ref"],
            # the cap is never below the steps of the time cap at dt_max (it matters only above eta = 6e-3)
            "steps": max(int(round(STEP_CAP * ETA / eta)), int(TIME_CAP / DT_MAX) + 1) if steps is None else int(steps),
            "record_every": record_every,
        },
        "kind": "mathematics" if block["id"] == "L0" else "finding",
        "preregistered": True,
        "producer": producer,
    }


@dataclass
class Watch:
    """Follows the relative position of the pair, row by row of the runner's time series (given to
    `runner.run` as `progress`), and names the event that ends the run."""

    close: float
    close_event: str
    far: float
    turn: float = 2 * math.pi * REVOLUTIONS
    event: str | None = None
    track: list = field(default_factory=list)  # (t, x, y, vx, vy of the second body relative to the first, angle)
    fastest: float = 0.0  # the largest speed of a body

    def start(self, state) -> None:
        rel, vel = state.x[1] - state.x[0], state.v[1] - state.v[0]
        self.track = [(float(state.t), float(rel[0]), float(rel[1]), float(vel[0]), float(vel[1]), 0.0)]
        self.fastest = float(np.max(np.linalg.norm(state.v, axis=1)))

    def __call__(self, row: dict) -> None:
        x, y = row["x1_0"] - row["x0_0"], row["x1_1"] - row["x0_1"]
        vx, vy = row["v1_0"] - row["v0_0"], row["v1_1"] - row["v0_1"]
        _, px, py, _, _, angle = self.track[-1]
        angle += math.atan2(px * y - py * x, px * x + py * y)  # the turn since the last row
        self.track.append((row["t"], x, y, vx, vy, angle))
        self.fastest = max(self.fastest, math.hypot(row["v0_0"], row["v0_1"]), math.hypot(row["v1_0"], row["v1_1"]))
        r = math.hypot(x, y)
        if not math.isfinite(r):
            self.event = "not finite"
        elif r <= self.close:
            self.event = self.close_event
        elif r >= self.far:
            self.event = "escape"
        elif abs(angle) >= self.turn:
            self.event = "three revolutions"

    def stop(self) -> bool:
        return self.event is not None


def find_event(case: str, law="reference", phi_bar_deg: float = PHI, eta: float = ETA) -> tuple:
    """The first pass: the runner with every step recorded, stopped at the event. Returns
    `(steps, watch)`; `watch.event` is the event, "no event" at the time cap or "step cap"."""
    sc = scales(case)
    watch = Watch(sc["close"], sc["close_event"], ESCAPE * sc["r0"])
    loaded = load_settings(settings(case, law, phi_bar_deg, eta))
    watch.start(loaded.state)
    try:
        result = run(loaded, out_root=None, progress=watch, should_stop=watch.stop)
    except Interrupted as stopped:
        return stopped.step, watch
    watch.event = "no event" if result.summary["stopped_by"] == "t_end" else "step cap"
    return result.summary["steps"], watch


def first_acceleration(case: str, law="reference", phi_bar_deg: float = PHI) -> dict:
    """The acceleration of the second body relative to the first at the start, in units of `a_ref`:
    `radial` (negative: a pull; positive: a push), `tangential` (a swirl; positive towards +90 degrees)
    and the name of what it is. A part counts when it is above 1e-12 of `a_ref`, and above 1e-12 of
    the whole acceleration (so that the rounding of a very large pull is not read as a swirl)."""
    loaded = load_settings(settings(case, law, phi_bar_deg))
    state, sc = loaded.state, scales(case)
    ctx = StepContext(loaded.make_law(), DirectSum(constants=SI, softening=0.0), loaded.frame)
    a = ctx.acceleration(state)
    rel = a[1] - a[0]
    unit = (state.x[1] - state.x[0]) / np.linalg.norm(state.x[1] - state.x[0])
    radial = float(rel @ unit) / sc["a_ref"]
    tangential = float(-unit[1] * rel[0] + unit[0] * rel[1]) / sc["a_ref"]
    floor = SILENT * max(1.0, math.hypot(radial, tangential))
    names = []
    if abs(radial) > floor:
        names.append("pull" if radial < 0 else "push")
    if abs(tangential) > floor:
        names.append("swirl")
    return {"radial": radial, "tangential": tangential, "what": " and ".join(names) or "none"}


def _crossing(track, quantity, level: float, last_only: bool = False):
    """The first place where `quantity(row)` passes `level` between two rows of the track, by linear
    interpolation: `(t, row values interpolated)`, or None. With `last_only`, only the last pair of rows."""
    start = max(1, len(track) - 1) if last_only else 1
    for i in range(start, len(track)):
        before, after = quantity(track[i - 1]) - level, quantity(track[i]) - level
        if before == 0 and i == 1:
            continue  # the start itself is not a crossing
        if before * after <= 0 and before != after:
            w = before / (before - after)
            return tuple(p + w * (q - p) for p, q in zip(track[i - 1], track[i]))
    return None


def _radius(row) -> float:
    return math.hypot(row[1], row[2])


def _radial_speed(row) -> float:
    return (row[1] * row[3] + row[2] * row[4]) / math.hypot(row[1], row[2])


def _event_time(watch: Watch, sc: dict) -> float:
    """The time of the event: where the separation (or the angle) passes its threshold inside the last
    step, by linear interpolation; the time of the last row for "no event" and "step cap"."""
    if watch.event in ("contact", "plunge"):
        hit = _crossing(watch.track, _radius, sc["close"], last_only=True)
    elif watch.event == "escape":
        hit = _crossing(watch.track, _radius, ESCAPE * sc["r0"], last_only=True)
    elif watch.event == "three revolutions":
        hit = _crossing(watch.track, lambda row: abs(row[5]), watch.turn, last_only=True)
    else:
        hit = None
    return watch.track[-1][0] if hit is None else hit[0]


def _revolutions(watch: Watch, sc: dict) -> list:
    """For each completed revolution `k`: its time and the closure error
    `sqrt((r/r0 - 1)^2 + (r'/v_ref)^2)` where the relative position has turned by `2 pi k`."""
    out = []
    for k in range(1, REVOLUTIONS + 1):
        hit = _crossing(watch.track, lambda row: abs(row[5]), 2 * math.pi * k)
        if hit is None:
            break
        r, speed = _radius(hit), _radial_speed(hit)
        out.append(
            {
                "revolution": k,
                "t": hit[0],
                "separation_over_r0": r / sc["r0"],
                "radial_speed_over_v_ref": speed / sc["v_ref"],
                "closure_error": math.hypot(r / sc["r0"] - 1, speed / sc["v_ref"]),
            }
        )
    return out


def _apsides(watch: Watch, sc: dict) -> dict:
    """The extremes of the separation (where the radial speed changes sign) and the apsidal precession:
    the angle turned between two apsides of the same kind, minus a full turn. Not defined for an orbit
    whose separation changes by less than 1e-3 of `r0` (a circle has no apsides)."""
    radii = [_radius(row) for row in watch.track]
    if (max(radii) - min(radii)) / sc["r0"] < 1e-3:
        return {"defined": False, "why": "the separation changes by less than 1e-3 of r0"}
    track, found = watch.track, []
    first = next((_radial_speed(row) for row in track[1:] if _radial_speed(row) != 0), 0.0)
    if _radial_speed(track[0]) == 0 and first != 0:  # the start is itself an apsis
        found.append({"kind": "inner" if first > 0 else "outer", "t": track[0][0], "angle": 0.0, "r": radii[0]})
    for i in range(2, len(track)):
        before, after = _radial_speed(track[i - 1]), _radial_speed(track[i])
        if before * after < 0:
            w = before / (before - after)
            row = tuple(p + w * (q - p) for p, q in zip(track[i - 1], track[i]))
            found.append({"kind": "inner" if after > 0 else "outer", "t": row[0], "angle": row[5], "r": _radius(row)})
    precession = []
    for kind in ("inner", "outer"):
        same = [a for a in found if a["kind"] == kind]
        precession += [abs(q["angle"] - p["angle"]) - 2 * math.pi for p, q in zip(same, same[1:])]
    inner = [a["r"] for a in found if a["kind"] == "inner"]
    outer = [a["r"] for a in found if a["kind"] == "outer"]
    return {
        "defined": True,
        "count": len(found),
        "inner_over_r0": [r / sc["r0"] for r in inner],
        "outer_over_r0": [r / sc["r0"] for r in outer],
        "radial_period": [q["t"] - p["t"] for p, q in zip(found, found[2:])][:1] or None,
        "precession_rad_per_radial_period": precession,
        "largest_precession_rad": max((abs(p) for p in precession), default=None),
    }


def _angular_momentum(state) -> float:
    return float(np.sum(state.m * (state.x[:, 0] * state.v[:, 1] - state.x[:, 1] * state.v[:, 0])))


@dataclass
class PairRun:
    """One T1 run: the runner's result (`stamp`, `summary`, `timeseries`, final `state`), the
    measurements of the pre-registration, and the settings that reproduce it with the plain runner."""

    result: RunResult
    measured: dict
    settings: dict
    watch: Watch


def run_case(
    case: str,
    law="reference",
    phi_bar_deg: float = PHI,
    eta: float = ETA,
    out_root=None,
    producer: str = "",
    scenario: str | None = None,
    rows: int = 500,
    now=None,
) -> PairRun:
    """A T1 case under a law: the first pass finds the step of the event (`find_event`), the second is
    the plain run to that step (written to `out_root/<scenario>/<run-id>/` if `out_root` is given,
    with about `rows` rows of time series). Returns the run and its measurements."""
    steps, watch = find_event(case, law, phi_bar_deg, eta)
    plain = settings(case, law, phi_bar_deg, eta, steps, max(1, steps // rows), producer, scenario)
    loaded = load_settings(plain)
    start = loaded.state.copy()
    result = run(loaded, out_root=out_root, now=now)
    end = result.state
    last = watch.track[-1]
    same = (float(end.x[1, 0] - end.x[0, 0]), float(end.x[1, 1] - end.x[0, 1])) == (last[1], last[2])
    if result.summary["steps"] != steps or not same:
        raise RuntimeError(f"{case}: the plain run to step {steps} does not end where the first pass did")

    sc = scales(case)
    law_object = loaded.make_law()
    solver = DirectSum(constants=SI, softening=0.0)
    radii = [_radius(row) for row in watch.track]
    summary = result.summary
    energy_key = "energy_re" if summary["initial"].get("energy_total") is None else "energy_total"
    unit_l = sc["reduced_mass"] * sc["r0"] * sc["v_ref"]
    measured = {
        "case": case,
        "law": law_label(law),
        "phi_bar_deg": phi_bar_deg,
        "eta": eta,
        "first_acceleration": first_acceleration(case, law, phi_bar_deg),
        "event": watch.event,
        "t_event": _event_time(watch, sc),
        "t_event_over_T_ref": _event_time(watch, sc) / sc["T_ref"],
        "steps": steps,
        "t_end": summary["t_end"],
        "separation_end_over_r0": radii[-1] / sc["r0"],
        "separation_min_over_r0": min(radii) / sc["r0"],
        "separation_max_over_r0": max(radii) / sc["r0"],
        "turned_revolutions": watch.track[-1][5] / (2 * math.pi),
        "fastest_speed_over_c": watch.fastest / SI.c,
        "revolutions": _revolutions(watch, sc),
        "apsides": _apsides(watch, sc),
        "conservation": dict(summary["conservation"], momentum_with_walls=summary["momentum_with_walls_drift"]),
        "energy_used": energy_key,
        # in units of k_ref / r0: defined also where the energy starts at zero and has no relative drift
        "energy_change_over_reference": abs(summary["final"][energy_key] - summary["initial"][energy_key])
        / (sc["k_ref"] / sc["r0"]),
        "angular_momentum_drift": abs(_angular_momentum(end) - _angular_momentum(start)) / unit_l,
        "binding_energy_start_eV": -summary["initial"][energy_key] / SI.e,
        "binding_energy_end_eV": -summary["final"][energy_key] / SI.e,
        "third_law_residual_start": third_law_residual(law_object, start, solver, loaded.frame)["relative"],
        "third_law_residual_end": third_law_residual(law_object, end, solver, loaded.frame)["relative"],
        "spin_end": [float(value) for value in end.spin],
        "reference": {key: sc[key] for key in ("r0", "T_ref", "v_ref", "a_ref", "close", "close_event")},
    }
    if loaded.law_spec.reading == "L1e" or end.xc.any():  # D49: the Composer position beside the Real one
        ratio = loaded.frame.s**2 * end.q / end.m
        measured["composer_position_end"] = {
            "x": end.x.tolist(),
            "xc": end.xc.tolist(),
            "x_minus_s2_q_over_m_xc": (end.x - ratio[:, None] * end.xc).tolist(),
        }
    return PairRun(result, measured, plain, watch)
