"""T2.2: the batch runner. Every definition on every scene, stamped, resumable, in parallel.
Software checks on small scenes; nothing here is a result."""

import copy
import json
import shutil
import signal

import pytest

from cosmicmass.analysis import batch
from cosmicmass.analysis.batch import (
    BatchError,
    Job,
    index,
    load_definitions,
    load_scenes,
    main,
    make_jobs,
    run_batch,
    timing,
)
from cosmicmass.runner import run
from cosmicmass.stamp import validate_outputs

PAIR = {
    "scenario": "pair",
    "universe": {"units": "scaled", "dimension": 2},
    "law": {"id": "L0"},
    "bodies": [
        {"name": "heavy", "m": 1.0, "q": 1e-3, "x": [0.0, 0.0], "v": [0.0, -5e-4]},
        {"name": "light", "m": 1e-3, "q": -1e-3, "x": [1.0, 0.0], "v": [0.0, 0.9]},
    ],
    "numerics": {"steps": 60, "softening": 0.0, "eta": 0.02, "record_every": 20},
    "kind": "mathematics",
    "producer": "unit test",
}
TRIO = copy.deepcopy(PAIR)
TRIO["bodies"].append({"name": "third", "m": 1e-3, "q": 0.0, "x": [0.0, 2.0], "v": [-0.6, 0.0]})
TRIO["numerics"] = {"t_end": 0.5, "softening": 0.01, "eta": 0.02, "record_every": 20}
SCENES = {"scenes": {"S1": PAIR, "S2": TRIO}}
DEFINITIONS = {
    "definitions": [
        {"id": "L0", "law": {"id": "L0"}},
        {"id": "L1a+rad-Re", "law": {"id": "L1a", "lambda": 1e-3}},
        {
            "id": "L1b-perp-Sh",
            "law": {"id": "L1b", "lambda": 1e-3, "sign": "-", "direction": "perpendicular", "inertia": "Sh"},
        },
    ]
}
REFUSED = {"id": "not-a-law", "law": {"id": "L9"}}


def jobs(extra=()):
    return make_jobs({"definitions": DEFINITIONS["definitions"] + list(extra)}, SCENES)


def test_definitions_and_scenes_are_read_and_checked(tmp_path):
    assert [name for name, _ in load_definitions(DEFINITIONS)] == ["L0", "L1a+rad-Re", "L1b-perp-Sh"]
    assert load_definitions(DEFINITIONS["definitions"]) == load_definitions(DEFINITIONS)
    (tmp_path / "pair.json").write_text(json.dumps(PAIR), encoding="utf-8")
    scenes_file = tmp_path / "scenes.json"
    scenes_file.write_text(
        json.dumps({"scenes": {"S1": "pair.json", "S2": TRIO}}), encoding="utf-8"
    )  # a path, read beside it
    definitions_file = tmp_path / "definitions.json"
    definitions_file.write_text(json.dumps(DEFINITIONS), encoding="utf-8")
    assert load_scenes(scenes_file) == load_scenes(SCENES) and load_scenes(SCENES["scenes"]) == load_scenes(SCENES)
    assert make_jobs(definitions_file, scenes_file) == jobs()
    for bad, named in (
        ([], "non-empty list"),
        ({"definitions": "L0"}, "non-empty list"),
        ([{"id": "a"}], r"definitions\[0\]"),
        ([{"id": "a", "law": {}, "colour": 1}], r"definitions\[0\]"),
        ([{"id": "a/b", "law": {}}], "names a folder"),
        ([{"id": "", "law": {}}], "names a folder"),
        ([{"id": ".hidden", "law": {}}], "names a folder"),
        ([{"id": 3, "law": {}}], "names a folder"),
        ([{"id": "a", "law": {}}, {"id": "a", "law": {}}], "used twice"),
    ):
        with pytest.raises(BatchError, match=named):
            load_definitions(bad)
    for bad, named in (({}, "non-empty object"), ({"a b": PAIR}, "names a folder"), ({"S": 3}, "settings object")):
        with pytest.raises(BatchError, match=named):
            load_scenes(bad)


def test_jobs_are_every_definition_on_every_scene_in_a_fixed_order():
    all_jobs = jobs()
    assert [j.name for j in all_jobs] == [
        "S1/L0",
        "S1/L1a+rad-Re",
        "S1/L1b-perp-Sh",
        "S2/L0",
        "S2/L1a+rad-Re",
        "S2/L1b-perp-Sh",
    ]
    job = all_jobs[4]
    assert job.settings["law"] == {"id": "L1a", "lambda": 1e-3} and job.settings["scenario"] == "S2/L1a+rad-Re"
    assert job.settings["bodies"] == TRIO["bodies"] and TRIO["law"] == {"id": "L0"}  # the scene is not changed
    assert len({j.key for j in all_jobs}) == 6 and job.key == jobs()[4].key  # a key per job, stable
    shards = [make_jobs(DEFINITIONS, SCENES, shard=(k, 4)) for k in range(4)]
    assert sorted(j.name for s in shards for j in s) == sorted(j.name for j in all_jobs)  # the shards cover it once
    assert [j.name for j in shards[1]] == ["S1/L1a+rad-Re", "S2/L1b-perp-Sh"]
    for bad in ((4, 4), (-1, 2), (0, 0), ("0", 2)):
        with pytest.raises(BatchError, match="shard"):
            make_jobs(DEFINITIONS, SCENES, shard=bad)


def test_a_batch_runs_stamps_and_is_not_run_twice(tmp_path):
    out, state = tmp_path / "results", tmp_path / "state"
    seen, handler = [], signal.getsignal(signal.SIGINT)
    first = run_batch(jobs([REFUSED]), out_root=out, state_dir=state, progress=seen.append)
    assert first["counts"] == {"done": 6, "refused": 2} and first["skipped"] == 0 and len(seen) == 8
    for job in jobs():
        record = first["records"][job.name]
        folder = out / job.scene / job.definition / record["run_id"]  # results/<scene>/<definition>/<run-id>/
        assert record["status"] == "done" and record["folder"] == folder.as_posix() and folder.is_dir()
        stamp = json.loads((folder / "stamp.json").read_text(encoding="utf-8"))
        summary = json.loads((folder / "summary.json").read_text(encoding="utf-8"))
        assert validate_outputs(stamp, summary) == [] and stamp["scenario"] == job.name
        assert stamp["law"]["id"] == job.settings["law"]["id"] and stamp["numerics"]["checkpoint_every"] == 1000
        assert record["conservation"] == summary["conservation"] and record["steps"] == summary["steps"]
        assert not (folder.parent / "checkpoint.npz").exists() and record["seconds"] > 0
        alone = run(job.settings, out_root=None)  # the batch runs what a single run would
        assert alone.summary == summary
    refused = first["records"]["S1/not-a-law"]
    assert refused["status"] == "refused" and "SettingsError" in refused["error"] and "folder" not in refused

    assert signal.getsignal(signal.SIGINT) is handler  # Ctrl-C is given back to the caller
    again = run_batch(jobs([REFUSED]), out_root=out, state_dir=state, progress=seen.append)
    assert again["counts"] == first["counts"] and again["skipped"] == 8 and len(seen) == 8  # nothing ran
    assert len(list(out.rglob("stamp.json"))) == 6

    changed = copy.deepcopy(SCENES)  # other settings for one scene: its jobs run again, the others do not
    changed["scenes"]["S1"]["numerics"]["steps"] = 40
    third = run_batch(make_jobs(DEFINITIONS, changed), out_root=out, state_dir=state)
    assert third["skipped"] == 3 and third["counts"] == {"done": 6} and third["records"]["S1/L0"]["steps"] == 40

    table = index(jobs(), state)  # the state now belongs to the changed batch for S1
    assert table["complete"] is False and table["jobs"]["S1/L0"] == {"status": "not started"}
    assert table["jobs"]["S2/L0"]["status"] == "done" and "key" not in table["jobs"]["S2/L0"]
    assert index(make_jobs(DEFINITIONS, changed), state)["complete"] is True

    lost = third["records"]["S2/L0"]["folder"]  # a finished job whose results are gone runs again
    shutil.rmtree(lost)
    assert run_batch(make_jobs(DEFINITIONS, changed), out_root=out, state_dir=state)["skipped"] == 5


def test_a_stopped_batch_is_taken_up_again_with_the_same_results(tmp_path, monkeypatch):
    whole = run_batch(jobs(), out_root=tmp_path / "whole", state_dir=tmp_path / "state-whole")
    calls = {"n": 0}

    def stop_after_100_steps():  # as a time limit or Ctrl-C would: in the middle of the second job
        calls["n"] += 1
        return calls["n"] > 100

    monkeypatch.setattr(batch, "_should_stop", stop_after_100_steps)
    out, state = tmp_path / "cut", tmp_path / "state-cut"
    stopped = run_batch(jobs(), out_root=out, state_dir=state, checkpoint_every=25)
    assert stopped["counts"] == {"done": 1, "interrupted": 1, "not started": 4}
    cut = stopped["records"]["S1/L1a+rad-Re"]
    assert cut["step"] == 39 and (out / "S1" / "L1a+rad-Re" / "checkpoint.npz").exists()
    monkeypatch.undo()

    resumed = run_batch(jobs(), out_root=out, state_dir=state, checkpoint_every=25)
    assert resumed["counts"] == {"done": 6} and resumed["skipped"] == 1
    assert resumed["records"]["S1/L1a+rad-Re"]["resumed_at_steps"] == [39]
    assert resumed["records"]["S1/L1a+rad-Re"]["seconds"] > cut["seconds"]  # the time of both parts
    for name, record in whole["records"].items():
        other = resumed["records"][name]
        assert other["conservation"] == record["conservation"] and other["steps"] == record["steps"]
        a = (tmp_path / "whole" / name / record["run_id"] / "timeseries.csv").read_bytes()
        assert a == (out / name / other["run_id"] / "timeseries.csv").read_bytes()  # byte for byte
    assert not list(out.rglob("checkpoint.npz"))


def test_a_checkpoint_of_another_run_does_not_block_a_job(tmp_path, monkeypatch):
    out, state = tmp_path / "results", tmp_path / "state"
    monkeypatch.setattr(batch, "_should_stop", lambda: True)
    assert run_batch(jobs()[:1], out_root=out, state_dir=state, workers=1)["counts"] == {"not started": 1}
    calls = {"n": 0}

    def once():
        calls["n"] += 1
        return calls["n"] == 10

    monkeypatch.setattr(batch, "_should_stop", once)
    assert run_batch(jobs()[:1], out_root=out, state_dir=state)["counts"] == {"interrupted": 1}
    monkeypatch.undo()
    other = copy.deepcopy(SCENES)
    other["scenes"]["S1"]["numerics"]["eta"] = 0.01  # the same job name, other numerics: the checkpoint is not its own
    done = run_batch(make_jobs(DEFINITIONS, other)[:1], out_root=out, state_dir=state)
    assert done["counts"] == {"done": 1} and done["records"]["S1/L0"]["resumed_at_steps"] == []


def test_several_processes_give_the_same_results(tmp_path):
    one = run_batch(jobs([REFUSED]), out_root=tmp_path / "one", state_dir=tmp_path / "state-one")
    two = run_batch(jobs([REFUSED]), out_root=tmp_path / "two", state_dir=tmp_path / "state-two", workers=2)
    assert two["counts"] == one["counts"] == {"done": 6, "refused": 2}
    for name, record in one["records"].items():
        if record["status"] == "done":
            assert two["records"][name]["conservation"] == record["conservation"]
            assert two["records"][name]["steps"] == record["steps"]
    limited = run_batch(
        jobs(), out_root=tmp_path / "late", state_dir=tmp_path / "state-late", workers=2, time_limit=0.0
    )
    assert set(limited["counts"]) == {"interrupted"} and limited["counts"]["interrupted"] == 6  # out of time at once


def test_a_failure_that_is_not_a_refusal_is_recorded_as_failed(tmp_path, monkeypatch):
    def broken(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(batch, "run", broken)
    outcome = run_batch(jobs()[:1], out_root=tmp_path / "r", state_dir=tmp_path / "s")
    assert outcome["counts"] == {"failed": 1} and "disk full" in outcome["records"]["S1/L0"]["error"]
    monkeypatch.undo()
    assert run_batch(jobs()[:1], out_root=tmp_path / "r", state_dir=tmp_path / "s")["counts"] == {
        "done": 1
    }  # run again


def test_the_timing_run_sizes_the_budget_and_writes_nothing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    report = timing(jobs([REFUSED]), steps=5)
    assert report["refused"] == ["S1/not-a-law", "S2/not-a-law"] and report["without_estimate"] == []
    pair, trio = report["jobs"]["S1/L1a+rad-Re"], report["jobs"]["S2/L1a+rad-Re"]
    assert pair["steps_planned"] == 60 and pair["bodies"] == 2 and pair["seconds_per_step"] > 0
    assert pair["seconds_estimated"] == pytest.approx(60 * pair["seconds_per_step"])
    assert trio["steps_planned"] is None and trio["steps_guessed"] > 5  # ends at a time: guessed from the sample
    assert report["seconds_estimated"] == pytest.approx(
        sum(j.get("seconds_estimated", 0.0) for j in report["jobs"].values())
    )
    assert "SettingsError" in report["jobs"]["S1/not-a-law"]["refused"] and not list(tmp_path.iterdir())
    open_ended = copy.deepcopy(SCENES)
    open_ended["scenes"]["S2"]["numerics"]["t_end"] = 0.0  # no time to scale by
    assert timing(make_jobs(DEFINITIONS, open_ended)[3:4], steps=3)["without_estimate"] == ["S2/L0"]


def test_command_line(tmp_path, capsys):
    definitions, scenes = tmp_path / "d.json", tmp_path / "s.json"
    definitions.write_text(json.dumps(DEFINITIONS), encoding="utf-8")
    scenes.write_text(json.dumps(SCENES), encoding="utf-8")
    common = [str(definitions), str(scenes), "--out", str(tmp_path / "results"), "--state", str(tmp_path / "state")]
    assert main(common + ["--index"]) == 1 and '"complete": false' in capsys.readouterr().out
    assert main(common + ["--timing", "3"]) == 0 and "seconds_per_step" in capsys.readouterr().out
    assert main(common + ["--shard", "0/2", "--checkpoint-every", "0"]) == 0
    printed = capsys.readouterr().out
    assert "3 jobs: 3 done; 0 were already finished" in printed and printed.count("done        S") == 3
    assert (
        main(common + ["--workers", "2"]) == 0 and "6 jobs: 6 done; 3 were already finished" in capsys.readouterr().out
    )
    assert main(common + ["--index"]) == 0 and '"complete": true' in capsys.readouterr().out
    assert main(common + ["--time-limit", "0", "--state", str(tmp_path / "other")]) == 1  # nothing could finish
    assert "not started" in capsys.readouterr().out
    assert main([str(tmp_path / "missing.json"), str(scenes)]) == 2 and "error:" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        main(common + ["--shard", "half"])
    capsys.readouterr()
    assert Job("S", "d", {"a": 1}).name == "S/d"
