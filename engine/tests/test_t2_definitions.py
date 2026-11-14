"""T2.1: the coupling definitions of the survival table (docs/M2-PLAN.md; decisions D45, D47)."""

import importlib.util
import json
import re
from pathlib import Path

from cosmicmass.laws import LawSpec, make_law
from cosmicmass.settings import load_settings

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("t2_definitions", ROOT / "tools" / "t2_definitions.py")
t2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t2)


def test_the_committed_file_is_up_to_date_and_the_engine_agrees(capsys):
    assert t2.main(["--check"]) == 0, capsys.readouterr().out


def test_102_definitions_and_18_cells_not_applicable():
    data = t2.grid()
    ids = [d["id"] for d in data["definitions"]] + [d["id"] for d in data["not_applicable"]]
    assert len(data["definitions"]) == 102 and len(data["not_applicable"]) == 18  # D45
    assert len(ids) == len(set(ids)) == 120  # every cell of the grid, once
    assert all(d["id"].startswith("L1e") and ".radial" not in d["id"] for d in data["not_applicable"])
    assert all(re.fullmatch(r"L1[a-e][+-](radial|perpendicular|velocity|torque)\.(Re|Sh|abs)", i) for i in ids)


def test_every_definition_is_a_law_block_a_run_accepts():
    data = t2.grid()
    for entry in data["definitions"]:
        law = entry["law"]
        assert law["lambda"] == 1.0 and law["terms"] == "all" and law["B_cal"] == "off"  # D47
        spec_ = LawSpec.from_dict(law)
        assert t2.definition_id(spec_.reading, spec_.sign, spec_.direction, spec_.inertia) == entry["id"]
        make_law(spec_)
    one = next(d for d in data["definitions"] if d["id"] == "L1b-torque.Sh")
    settings = {"scenario": "x", "law": one["law"], "bodies": [{"preset": "proton", "x": [0.0, 0.0]}]}
    assert load_settings(settings).law_spec.direction == "torque"


def test_the_ids_are_stable():
    # the order and the names are part of the survival table: changing them changes every reference to them
    data = json.loads((ROOT / "scenarios" / "T2" / "definitions.json").read_text(encoding="utf-8"))
    ids = [d["id"] for d in data["definitions"]]
    assert ids[:3] == ["L1a+radial.Re", "L1a+radial.Sh", "L1a+radial.abs"]
    assert ids[-1] == "L1e-radial.abs" and "L1e+radial.Re" in ids
