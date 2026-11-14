#!/usr/bin/env python3
"""T2.1: the coupling definitions of the survival table (milestone M2, `docs/M2-PLAN.md`).

The grid of `docs/05-scenarios-and-tests.md`, T2: the readings L1a to L1e x the sign (L2: +, -) x the
direction (L4: radial, perpendicular, velocity, torque) x the inertia (L7: Re, Sh, abs) = 120 cells.

- Decision **D45**: L1e has no Composer force to direct, so its 18 non-radial cells are **not applicable**
  (the engine refuses them); the grid has **102 definitions**.
- Decision **D47** and the default exploratory law of `docs/05`: every definition has the strength
  lambda = 1, all terms (L5), the moving-body term off (L6) and the range 1/r^2 (L8). The lambda bounds
  of T2.4 vary the strength later; the definition stays the same.

Each definition has a stable ID, `<reading><sign><direction>.<inertia>` (for example `L1a+radial.Re`,
`L1e-radial.abs`), and a `law` block that goes unchanged into the `law` of a `settings.json`.

Writes `scenarios/T2/definitions.json`, in the fixed order of the grid, and checks every entry against the
engine: each definition builds its law; each cell that is not applicable is refused.

    python tools/t2_definitions.py            # write the file
    python tools/t2_definitions.py --check    # exit 1 if the file is missing or out of date
"""

import argparse
import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "engine"))

from cosmicmass.laws import LawSpec, make_law  # noqa: E402

OUT = ROOT / "scenarios" / "T2" / "definitions.json"
READINGS = ("L1a", "L1b", "L1c", "L1d", "L1e")
SIGNS = ("+", "-")
DIRECTIONS = ("radial", "perpendicular", "velocity", "torque")
INERTIAS = ("Re", "Sh", "abs")
FIXED = {"lambda": 1.0, "terms": "all", "B_cal": "off", "range": "1/r^2"}  # D47, default exploratory law
NOT_APPLICABLE = "D45: L1e has no Composer force to direct; only the radial direction exists for it"


def definition_id(reading: str, sign: str, direction: str, inertia: str) -> str:
    return f"{reading}{sign}{direction}.{inertia}"


def law_block(reading: str, sign: str, direction: str, inertia: str) -> dict:
    return {"id": reading, "sign": sign, "direction": direction, "inertia": inertia, **FIXED}


def grid() -> dict:
    """The definitions and the cells that are not applicable, in the fixed order of the grid."""
    definitions, not_applicable = [], []
    for reading, sign, direction, inertia in itertools.product(READINGS, SIGNS, DIRECTIONS, INERTIAS):
        key = definition_id(reading, sign, direction, inertia)
        if reading == "L1e" and direction != "radial":
            not_applicable.append({"id": key, "reason": NOT_APPLICABLE})
        else:
            definitions.append({"id": key, "law": law_block(reading, sign, direction, inertia)})
    return {
        "about": "T2.1 (docs/M2-PLAN.md): the coupling definitions of the survival table; written by "
        "tools/t2_definitions.py, do not edit by hand",
        "grid": {"readings": READINGS, "signs": SIGNS, "directions": DIRECTIONS, "inertias": INERTIAS},
        "fixed": FIXED,
        "count": {"definitions": len(definitions), "not_applicable": len(not_applicable)},
        "definitions": definitions,
        "not_applicable": not_applicable,
    }


def problems(data: dict) -> list[str]:
    """Every definition must build its law in the engine; every cell that is not applicable must be refused."""
    out = []
    for entry in data["definitions"]:
        try:
            make_law(LawSpec.from_dict(entry["law"]))
        except (ValueError, TypeError, NotImplementedError) as exc:
            out.append(f"{entry['id']}: the engine refuses it ({exc})")
    for entry in data["not_applicable"]:
        reading, rest = entry["id"][:3], entry["id"][3:]
        sign, direction_inertia = rest[0], rest[1:]
        direction, inertia = direction_inertia.split(".")
        try:
            make_law(LawSpec.from_dict(law_block(reading, sign, direction, inertia)))
            out.append(f"{entry['id']}: marked not applicable, but the engine accepts it")
        except (ValueError, TypeError, NotImplementedError):
            pass
    return out


def text(data: dict) -> str:
    return json.dumps(data, indent=2, ensure_ascii=False) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="fail if the file is missing or out of date")
    args = parser.parse_args(argv)
    data = grid()
    errors = problems(data)
    if errors:
        print("Definitions the engine does not agree with:")
        for error in errors:
            print("  -", error)
        return 1
    if args.check:
        if not OUT.exists() or OUT.read_text(encoding="utf-8").replace("\r\n", "\n") != text(data):
            print(f"{OUT.relative_to(ROOT).as_posix()} is out of date; run: python tools/t2_definitions.py")
            return 1
        print(
            f"T2 definitions: {data['count']['definitions']} definitions, {data['count']['not_applicable']} not applicable, up to date"
        )
        return 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text(data), encoding="utf-8", newline="\n")
    print(f"wrote {OUT.relative_to(ROOT).as_posix()}: {data['count']['definitions']} definitions, "
          f"{data['count']['not_applicable']} not applicable")  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
