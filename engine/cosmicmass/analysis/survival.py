"""The classifier of the survival table (task T2.3, milestone M2; decisions D43 and D44).

Every rule here is the pre-registration of `team/notes/T2.md` ("T2.3: the classifier"), committed before
this code was written. A coupling definition of `scenarios/T2/definitions.json` is classed:

- **A**: conserving and angle-free;
- **B**: conserving, angle-dependent;
- **C**: not conserving (the table says which condition failed, and where);
- **incomplete**: a run it needs is missing.

"Conserving" (D43), for every run of the definition (T0.2 and T1.1-T1.6, at every angle of D44):
(a) the mass and charge drifts are exactly 0; (b) the momentum with the walls' share drifts by at most
1e-6; (c) inside D23's scope the energy that the law conserves (E9's `energy_total`) drifts by at most
1e-6, and outside it the closed-loop work test (`closed_loop_ratio`) gives at most 1e-9 at every angle.

"Angle-free" (D44): the Real-coordinate results of every setup are the same, to 1e-9 relative, at the six
angles of T0.3.

    from cosmicmass.analysis.survival import closed_loop_ratio, classify, load_records
"""

import csv
import json
from pathlib import Path

import numpy as np

from ..constants import scaled
from ..energy import D23_SCOPE
from ..fields.direct import DirectSum
from ..frame import Frame
from ..laws import LawSpec, make_law
from ..state import State

__all__ = [
    "ANGLES",
    "SETUPS",
    "MOMENTUM_TOLERANCE",
    "ENERGY_TOLERANCE",
    "LOOP_TOLERANCE",
    "ANGLE_TOLERANCE",
    "in_d23_scope",
    "closed_loop_ratio",
    "load_records",
    "classify",
]

ANGLES = (10.0, 45.0, 79.86, 90.0, 135.0, 170.0)  # D44: the angles of T0.3
SETUPS = ("T0.2", "T1.1", "T1.2", "T1.3", "T1.4", "T1.5", "T1.6")
MOMENTUM_TOLERANCE = 1e-6  # D43 (b)
ENERGY_TOLERANCE = 1e-6  # D43 (c), inside D23's scope
LOOP_TOLERANCE = 1e-9  # D43 (c), outside it
ANGLE_TOLERANCE = 1e-9  # D44

# The closed-loop work test, as pre-registered: scaled units, 2D, two fixed sources, one test body.
SOURCES = ((1.0, 0.3, (0.0, 0.0)), (0.5, -0.4, (3.0, 1.0)))  # (m, q, place)
TEST_BODY = (0.2, 0.15)  # (m, q): a charge-to-mass ratio unlike both sources
LOOPS = (
    ((1.0, 1.5), (2.0, 1.5), (2.0, 2.5), (1.0, 2.5)),  # encloses no source
    ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)),  # encloses S1
)
NODES = 64  # Gauss-Legendre points per side


def in_d23_scope(law: dict) -> bool:
    """True if a `law` block is inside the scope of decision D23 (where E9 defines the energy)."""
    spec = LawSpec.from_dict(law)
    return all(getattr(spec, key) == value for key, value in D23_SCOPE.items())


def _loop_work(law, frame: Frame, solver: DirectSum, sources: State, corners) -> tuple[float, float]:
    """`(W_loop, max |W_side|)` of the law's force on the test body carried round one loop."""
    t, w = np.polynomial.legendre.leggauss(NODES)
    sides = []
    for k in range(4):
        start, end = np.array(corners[k], float), np.array(corners[(k + 1) % 4], float)
        length = float(np.linalg.norm(end - start))
        tangent = (end - start) / length
        points = start + np.outer((t + 1) / 2, end - start)
        m, q = TEST_BODY
        bodies = State(
            m=np.full(NODES, m),
            q=np.full(NODES, q),
            x=points,
            v=np.tile(tangent, (NODES, 1)),  # along the path, unit speed (the direction "velocity" reads it)
            fixed=np.zeros(NODES, bool),
            alive=np.ones(NODES, bool),
        )
        # the field of the sources only: the test bodies do not act on each other
        force = law.forces(bodies, solver.sample(sources, points), frame).total()
        sides.append(float(np.sum(w * (force @ tangent)) * length / 2))
    return float(np.sum(sides)), float(np.max(np.abs(sides)))


def closed_loop_ratio(law_block: dict, angle_deg: float) -> float:
    """The pre-registered closed-loop work test of one definition at one frame angle: the larger, over the
    two loops, of `|W_loop| / max(|W_side|)`; 0 where the force vanishes on a loop."""
    law, frame = make_law(LawSpec.from_dict(law_block)), Frame.scaled(angle_deg)
    solver = DirectSum(constants=scaled(), softening=0.0)
    sources = State(
        m=[s[0] for s in SOURCES],
        q=[s[1] for s in SOURCES],
        x=[s[2] for s in SOURCES],
        v=[[0.0, 0.0]] * len(SOURCES),
        fixed=[True] * len(SOURCES),
        alive=[True] * len(SOURCES),
    )
    worst = 0.0
    for corners in LOOPS:
        loop, side = _loop_work(law, frame, solver, sources, corners)
        if side > 0:
            worst = max(worst, abs(loop) / side)
    return worst


def _last_row(folder: Path) -> dict:
    path = folder / "timeseries.csv"
    if not path.exists():
        return {}
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {k: float(v) for k, v in rows[-1].items() if v not in ("", None)} if rows else {}


def _setup_of(scenario: str) -> str | None:
    """The setup (T0.2, T1.1 ... T1.6) a run's scenario name begins with (T2.2 names a job `<scene>/<definition>`,
    and a scene may carry its angle, for example `T1.1@45`); None for a run of something else."""
    head = scenario.split("/")[0].split("@")[0]
    return head if head in SETUPS else None


def _angle_of(stamp: dict) -> float | None:
    """The angle of D44 a run was made at (from its stamp), or None if it is not one of the six."""
    value = float(stamp["universe"]["phi_bar_deg"])
    return next((a for a in ANGLES if abs(a - value) <= 1e-9), None)


def load_records(root, definitions: dict) -> list[dict]:
    """The runs of T2.2 anywhere under `root`, **read by their stamps**, whatever the folder layout: the
    definition by the stamp's `law` block (each of the definitions has its own), the setup by the scenario name,
    the angle by `universe.phi_bar_deg`. The newest run (by run id) of each definition, setup and angle counts.
    A run whose folder path names a definition that its stamp's law is not raises ValueError; runs of other
    laws, setups or angles are left out. Returns records `{"definition", "setup", "angle", "summary", "final"}`."""
    by_law = {
        json.dumps(LawSpec.from_dict(d["law"]).to_dict(), sort_keys=True): d["id"] for d in definitions["definitions"]
    }
    ids = set(by_law.values())
    newest = {}
    root = Path(root)
    for stamp_path in sorted(root.rglob("stamp.json")) if root.exists() else []:
        run = stamp_path.parent
        if not (run / "summary.json").exists():
            continue
        stamp = json.loads(stamp_path.read_text(encoding="utf-8"))
        definition = by_law.get(json.dumps(stamp.get("law", {}), sort_keys=True))
        named = ids.intersection(run.relative_to(root).parts)  # a definition the folder path names, if any
        if named and named != {definition}:
            raise ValueError(
                f"{run}: the folder names {sorted(named)}, but the stamp's law is {definition or 'none of them'}"
            )
        setup, angle = _setup_of(str(stamp.get("scenario", ""))), _angle_of(stamp)
        if definition is None or setup is None or angle is None:
            continue
        key = (definition, setup, angle)
        if key not in newest or run.name > newest[key].name:  # run ids sort by time: the newest counts
            newest[key] = run
    records = []
    for (definition, setup, angle), run in sorted(newest.items()):
        summary = json.loads((run / "summary.json").read_text(encoding="utf-8"))
        records.append({"definition": definition, "setup": setup, "angle": angle, "summary": summary,
                        "final": _last_row(run)})  # fmt: skip
    return records


def _angle_values(record: dict) -> dict:
    """The Real-coordinate results compared across angles (D44): final places and velocities for T1, the
    final energy, kinetic energy and momentum for T0.2."""
    if record["setup"] == "T0.2":
        final = record["summary"]["final"]
        values = {"energy_re": final["energy_re"], "kinetic": final["kinetic"]}
        values.update({f"p{k}": p for k, p in enumerate(final["momentum"])})
        return values
    return {k: v for k, v in record["final"].items() if k[:1] in ("x", "v") and "_" in k}


def classify(definition_id: str, law_block: dict, records: list[dict], loop_ratios: dict | None = None) -> dict:
    """The class of one definition from its runs (`load_records`) and, outside D23's scope, its closed-loop
    ratios by angle (`closed_loop_ratio`). Returns `{"id", "class", "failed", "numbers"}`: "failed" lists
    each condition of D43 that failed, with the setup and angle where it failed."""
    mine = [r for r in records if r["definition"] == definition_id]
    have = {(r["setup"], r["angle"]) for r in mine}
    missing = [(s, a) for s in SETUPS for a in ANGLES if (s, a) not in have]
    if missing:
        return {"id": definition_id, "class": "incomplete", "failed": [], "numbers": {"missing": missing}}

    failed, numbers = [], {"momentum": 0.0, "energy": 0.0, "loop": 0.0, "angle": 0.0}
    d23 = in_d23_scope(law_block)
    for r in mine:
        drift, where = r["summary"]["conservation"], f"{r['setup']} at {r['angle']:g} degrees"
        if drift["mass"] != 0 or drift["charge"] != 0:
            failed.append(f"(a) mass or charge, {where}")
        momentum = r["summary"].get("momentum_with_walls_drift")
        momentum = drift["momentum"] if momentum is None else momentum
        numbers["momentum"] = max(numbers["momentum"], momentum or 0.0)
        if momentum is None or momentum > MOMENTUM_TOLERANCE:
            failed.append(f"(b) momentum {momentum}, {where}")
        if d23:
            energy = drift.get("energy_total")
            numbers["energy"] = max(numbers["energy"], energy or 0.0)
            if energy is None or energy > ENERGY_TOLERANCE:
                failed.append(f"(c) energy {energy}, {where}")
    if not d23:
        ratios = loop_ratios or {}
        for angle in ANGLES:
            ratio = ratios.get(angle)
            numbers["loop"] = max(numbers["loop"], ratio if ratio is not None else np.inf)
            if ratio is None or ratio > LOOP_TOLERANCE:
                failed.append(f"(c) closed-loop work {ratio}, at {angle:g} degrees")

    for setup in SETUPS:
        by_angle = [_angle_values(r) for r in mine if r["setup"] == setup]
        for key in set().union(*by_angle):
            column = [v.get(key) for v in by_angle]
            if any(c is None for c in column):  # a quantity missing at one angle: not the same
                spread = np.inf
            else:
                scale = max(abs(c) for c in column)
                spread = (max(column) - min(column)) / scale if scale > 0 else 0.0
            numbers["angle"] = max(numbers["angle"], spread)
    if failed:
        return {"id": definition_id, "class": "C", "failed": failed, "numbers": numbers}
    return {"id": definition_id, "class": "A" if numbers["angle"] <= ANGLE_TOLERANCE else "B", "failed": [],
            "numbers": numbers}  # fmt: skip
