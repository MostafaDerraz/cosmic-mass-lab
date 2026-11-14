"""The batch runner (task T2.2): every coupling definition on every scene, stamped, resumable, in parallel.

A **definition** is a law block with a stable id (`{"id": "L1b+perp-Sh", "law": {...}}`, from the
generator of T2.1). A **scene** is a complete settings file with an id (`"T1.1"`, from the scenarios
of T1, or T0.2). A **job** is one definition on one scene: the scene's settings with the law block
replaced. Each job is an ordinary run of the engine (`runner.run`), stamped like any other, in

    <out>/<scene>/<definition id>/<run-id>/        (the layout `results/T1.4/L0/` already uses)

What the batch adds:

- **Resumable.** A finished job is recorded in the state folder and skipped the next time, as long as
  its settings and the engine commit are the same. A job that was stopped resumes from its checkpoint
  (E7.1), and gives the same files as if it had not been stopped.
- **Parallel.** `workers` processes on one machine; `shard=(k, n)` takes every n-th job, for n machines
  (a CI matrix). The order of the jobs is fixed: scenes in the order given, definitions in the order given.
- **A job the engine refuses does not stop the batch.** Its error is recorded with the job
  (`status: "refused"`), for example a moving body with zero inertia under an inertia reading.
- **A time limit.** After `time_limit` seconds the running jobs stop at their next step with a
  checkpoint (`status: "interrupted"`), and the next call continues them.
- **A timing run** (`timing`): a few steps of every job, to size the compute budget before the batch.

The state folder holds one small JSON file per job. It is not a result and is not committed; the
results are the stamped run folders. `index` gathers them into one table.

    python -m cosmicmass.analysis.batch definitions.json scenes.json --out results --workers 4
    python -m cosmicmass.analysis.batch definitions.json scenes.json --timing 20
"""

import argparse
import copy
import hashlib
import json
import os
import signal
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from ..checkpoint import CheckpointError
from ..runner import CHECKPOINT_NAME, Interrupted, run
from ..settings import SettingsError, load_settings
from ..stamp import ResultsError, engine_commit

__all__ = ["Job", "BatchError", "load_definitions", "load_scenes", "make_jobs", "run_batch", "timing", "index", "main"]

_SAFE = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.+-=,_")
DONE, REFUSED, INTERRUPTED, FAILED = "done", "refused", "interrupted", "failed"
# what the engine raises when it will not run a definition on a scene: recorded, not fatal
_REFUSALS = (SettingsError, ZeroDivisionError, NotImplementedError, ValueError)


class BatchError(ValueError):
    """The definitions or the scenes of a batch are invalid; the message names what is wrong."""


@dataclass(frozen=True)
class Job:
    scene: str
    definition: str
    settings: dict  # the scene's settings with the definition's law; scenario = "<scene>/<definition>"

    @property
    def name(self) -> str:
        return f"{self.scene}/{self.definition}"

    @property
    def key(self) -> str:
        """What a finished job is recognised by: its settings (the engine commit is checked beside it)."""
        text = json.dumps(self.settings, sort_keys=True)
        return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _check_id(value, what: str) -> str:
    if not isinstance(value, str) or not value or set(value) - _SAFE or value.startswith("."):
        raise BatchError(f"{what} {value!r} must be a name made of letters, digits and . + - = , _ (it names a folder)")
    return value


def _read(source, key: str):
    data = source
    if not isinstance(source, (dict, list)):
        data = json.loads(Path(source).read_text(encoding="utf-8"))
    return data[key] if isinstance(data, dict) and key in data else data


def load_definitions(source) -> list:
    """The definitions of a batch: a list of `{"id", "law"}`, or `{"definitions": [...]}`, or a path
    to such a file. Returns `[(id, law), ...]` in the order given; ids must be distinct."""
    entries = _read(source, "definitions")
    if not isinstance(entries, list) or not entries:
        raise BatchError("the definitions must be a non-empty list of {'id', 'law'}")
    out, seen = [], set()
    for k, entry in enumerate(entries):
        if not isinstance(entry, dict) or set(entry) - {"id", "law", "note"} or "law" not in entry:
            raise BatchError(f"definitions[{k}] must be an object with 'id' and 'law'")
        name = _check_id(entry.get("id"), f"definitions[{k}].id")
        if name in seen:
            raise BatchError(f"definitions[{k}].id {name!r} is used twice")
        seen.add(name)
        out.append((name, entry["law"]))
    return out


def load_scenes(source) -> list:
    """The scenes of a batch: `{"scenes": {id: settings or path}}`, or that inner object, or a path
    to such a file. A path given for a scene is read relative to the file that names it. Returns
    `[(id, settings), ...]` in the order given."""
    base = Path(source).parent if not isinstance(source, (dict, list)) else Path()
    entries = _read(source, "scenes")
    if not isinstance(entries, dict) or not entries:
        raise BatchError("the scenes must be a non-empty object {id: settings}")
    out = []
    for name, settings in entries.items():
        _check_id(name, "the scene id")
        if isinstance(settings, str):
            settings = json.loads((base / settings).read_text(encoding="utf-8"))
        if not isinstance(settings, dict):
            raise BatchError(f"scene {name!r} must be a settings object or a path to one")
        out.append((name, settings))
    return out


def make_jobs(definitions, scenes, shard: tuple | None = None) -> list:
    """Every definition on every scene, scenes first: `[Job, ...]`. With `shard=(k, n)`: the jobs
    number k, k + n, k + 2n, ... of that list (k from 0), so that n calls cover it once."""
    definitions = definitions if _is_pairs(definitions) else load_definitions(definitions)
    scenes = scenes if _is_pairs(scenes) else load_scenes(scenes)
    jobs = []
    for scene, settings in scenes:
        for name, law in definitions:
            config = copy.deepcopy(settings)
            config["law"] = copy.deepcopy(law)
            config["scenario"] = f"{scene}/{name}"
            jobs.append(Job(scene, name, config))
    if shard is not None:
        k, n = shard
        if not (isinstance(k, int) and isinstance(n, int) and 0 <= k < n):
            raise BatchError("a shard is (k, n) with 0 <= k < n")
        jobs = jobs[k::n]
    return jobs


def _is_pairs(value) -> bool:
    return isinstance(value, list) and all(isinstance(item, tuple) for item in value) and bool(value)


# ------------------------------------------------------------------ one job
_STOP = {"asked": False, "deadline": None}


def _should_stop() -> bool:
    return _STOP["asked"] or (_STOP["deadline"] is not None and time.monotonic() >= _STOP["deadline"])


def _worker_start(deadline) -> None:
    """In each worker: Ctrl-C asks the run to stop after its step (it then writes its checkpoint)."""
    _STOP["asked"], _STOP["deadline"] = False, deadline
    try:
        signal.signal(signal.SIGINT, lambda *_: _STOP.__setitem__("asked", True))
    except ValueError:  # not in the main thread of its process
        pass


def _state_file(state_dir: Path, job: Job) -> Path:
    return Path(state_dir) / job.scene / f"{job.definition}.json"


def _load_state(state_dir: Path, job: Job) -> dict | None:
    path = _state_file(state_dir, job)
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _save_state(state_dir: Path, job: Job, record: dict) -> None:
    path = _state_file(state_dir, job)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".writing")
    temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _run_job(job: Job, out_root: str, state_dir: str, checkpoint_every: int | None) -> dict:
    """Run one job, resuming its checkpoint if it has one; returns and stores its record."""
    out_root, started = Path(out_root), time.perf_counter()
    record = {"scene": job.scene, "definition": job.definition, "key": job.key, "engine_commit": engine_commit()}
    settings = copy.deepcopy(job.settings)
    if checkpoint_every:
        settings.setdefault("numerics", {})["checkpoint_every"] = checkpoint_every
    checkpoint = out_root / job.scene / job.definition / CHECKPOINT_NAME
    previous = _load_state(state_dir, job) or {}
    spent = previous.get("seconds", 0.0) if previous.get("status") == INTERRUPTED else 0.0
    try:
        loaded = load_settings(settings)
        try:
            result = run(loaded, out_root=out_root, resume=checkpoint.exists(), should_stop=_should_stop)
        except CheckpointError:  # a checkpoint of other settings or of another engine: start again
            checkpoint.unlink(missing_ok=True)
            result = run(loaded, out_root=out_root, should_stop=_should_stop)
        except FileExistsError:  # a run of this job ended in this very second: its run id is taken
            time.sleep(1.1)
            result = run(loaded, out_root=out_root, should_stop=_should_stop)
        conservation = result.summary["conservation"]
        record.update(
            status=DONE,
            folder=result.folder.as_posix(),
            run_id=result.stamp["run_id"],
            steps=result.summary["steps"],
            t_end=result.summary["t_end"],
            conservation=conservation,
            resumed_at_steps=result.stamp.get("resumed_at_steps", []),
        )
    except Interrupted as stop:
        record.update(
            status=INTERRUPTED, step=stop.step, checkpoint=None if stop.checkpoint is None else str(stop.checkpoint)
        )
    except _REFUSALS as exc:
        record.update(status=REFUSED, error=f"{type(exc).__name__}: {exc}")
    except (ResultsError, OSError, ArithmeticError) as exc:  # not a refusal: something went wrong
        record.update(status=FAILED, error=f"{type(exc).__name__}: {exc}")
    record["seconds"] = spent + time.perf_counter() - started
    _save_state(Path(state_dir), job, record)
    return record


def _finished(record: dict | None, job: Job, commit: str) -> bool:
    """A job is not run again when it was done or refused with the same settings on the same engine."""
    return bool(
        record
        and record.get("status") in (DONE, REFUSED)
        and record.get("key") == job.key
        and record.get("engine_commit") == commit
        and (record["status"] == REFUSED or Path(record.get("folder", "")).is_dir())
    )


# ------------------------------------------------------------------ the batch
def run_batch(
    jobs: list,
    out_root="results",
    state_dir=".batch",
    workers: int = 1,
    checkpoint_every: int | None = 1000,
    time_limit: float | None = None,
    progress=None,
) -> dict:
    """Run the jobs that are not finished yet. Returns `{"records": {job name: record}, "counts":
    {status: number}, "skipped": number already finished, "seconds": wall time}`. `progress(record)`
    is called as each job ends. With `workers=1` the jobs run in this process, one after another."""
    started, commit = time.perf_counter(), engine_commit()
    out_root, state_dir = Path(out_root), Path(state_dir)
    records, todo = {}, []
    for job in jobs:
        record = _load_state(state_dir, job)
        if _finished(record, job, commit):
            records[job.name] = record
        else:
            todo.append(job)
    skipped = len(records)
    deadline = None if time_limit is None else time.monotonic() + time_limit
    arguments = (str(out_root), str(state_dir), checkpoint_every)

    def note(record):
        records[f"{record['scene']}/{record['definition']}"] = record
        if progress is not None:
            progress(record)

    if workers <= 1:
        try:
            before = signal.getsignal(signal.SIGINT)
        except ValueError:
            before = None
        _worker_start(deadline)
        try:
            for job in todo:
                if _should_stop():
                    break  # asked to stop, or out of time: the rest wait for the next call
                note(_run_job(job, *arguments))
        finally:  # this process is the caller's: give Ctrl-C back
            if before is not None:
                try:
                    signal.signal(signal.SIGINT, before)
                except ValueError:
                    pass
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=_worker_start, initargs=(deadline,)) as pool:
            futures = [pool.submit(_run_job, job, *arguments) for job in todo]
            for future in as_completed(futures):
                note(future.result())
    counts = {}
    for record in records.values():
        counts[record["status"]] = counts.get(record["status"], 0) + 1
    waiting = len(jobs) - len(records)
    if waiting:
        counts["not started"] = waiting
    return {"records": records, "counts": counts, "skipped": skipped, "seconds": time.perf_counter() - started}


def timing(jobs: list, steps: int = 20) -> dict:
    """A timing run: `steps` steps of every job, nothing written. For each job the seconds per step
    and, where the settings give a number of steps, the estimate for the whole run; a refused job is
    listed with its error. Returns `{"jobs": {name: {...}}, "seconds_estimated", "refused": [...],
    "without_estimate": [...]}`. The estimate is rough: the first steps of a run are not its average."""
    out, total, refused, unknown = {}, 0.0, [], []
    for job in jobs:
        settings = copy.deepcopy(job.settings)
        numerics = settings.setdefault("numerics", {})
        planned = numerics.get("steps")
        numerics.update(steps=steps, t_end=None, checkpoint_every=None)
        started = time.perf_counter()
        try:
            result = run(settings, out_root=None)
        except _REFUSALS as exc:
            out[job.name] = {"refused": f"{type(exc).__name__}: {exc}"}
            refused.append(job.name)
            continue
        per_step = (time.perf_counter() - started) / max(1, result.summary["steps"])
        entry = {"seconds_per_step": per_step, "steps_planned": planned, "bodies": int(result.state.n)}
        if planned:
            entry["seconds_estimated"] = per_step * planned
            total += entry["seconds_estimated"]
        else:  # the run ends at a time, not at a step count: the elapsed time of the sample gives a guess
            t_end, elapsed = job.settings.get("numerics", {}).get("t_end"), result.summary["t_end"]
            if t_end and elapsed > 0:
                entry["steps_guessed"] = int(steps * t_end / elapsed)
                entry["seconds_estimated"] = per_step * entry["steps_guessed"]
                total += entry["seconds_estimated"]
            else:
                unknown.append(job.name)
        out[job.name] = entry
    return {"jobs": out, "seconds_estimated": total, "refused": refused, "without_estimate": unknown}


def index(jobs: list, state_dir=".batch") -> dict:
    """The table of a batch from its state folder: for every job its status and, when done, its run
    folder, steps and conservation drifts; `"complete"` is True when every job is done or refused."""
    rows, complete = {}, True
    for job in jobs:
        record = _load_state(Path(state_dir), job)
        if record is None or record.get("key") != job.key:
            rows[job.name] = {"status": "not started"}
            complete = False
            continue
        rows[job.name] = {k: v for k, v in record.items() if k not in ("key", "scene", "definition")}
        complete = complete and record["status"] in (DONE, REFUSED)
    return {"jobs": rows, "complete": complete}


# ------------------------------------------------------------------ command line
def _shard(text: str) -> tuple:
    try:
        k, n = text.split("/")
        return int(k), int(n)
    except ValueError:
        raise argparse.ArgumentTypeError("a shard is written k/n, for example 0/8") from None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m cosmicmass.analysis.batch", description=__doc__.split("\n")[0])
    parser.add_argument("definitions", help="JSON file: {'definitions': [{'id', 'law'}, ...]}")
    parser.add_argument("scenes", help="JSON file: {'scenes': {id: settings or path}}")
    parser.add_argument("--out", default="results", help="root folder of the results (default: results)")
    parser.add_argument("--state", default=".batch", help="folder of the job records (default: .batch)")
    parser.add_argument("--workers", type=int, default=1, help="processes on this machine (default: 1)")
    parser.add_argument("--shard", type=_shard, default=None, metavar="K/N", help="take every N-th job, starting at K")
    parser.add_argument("--checkpoint-every", type=int, default=1000, help="steps between checkpoints (0: none)")
    parser.add_argument(
        "--time-limit", type=float, default=None, help="seconds after which the jobs stop with a checkpoint"
    )
    parser.add_argument("--timing", type=int, default=None, metavar="STEPS", help="only time STEPS steps of every job")
    parser.add_argument("--index", action="store_true", help="only print the table of the batch from its state folder")
    args = parser.parse_args(argv)
    try:
        jobs = make_jobs(args.definitions, args.scenes, shard=args.shard)
    except (BatchError, OSError, json.JSONDecodeError) as exc:
        print(f"error: {exc}")
        return 2
    if args.timing is not None:
        report = timing(jobs, args.timing)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.index:
        table = index(jobs, args.state)
        print(json.dumps(table, indent=2, sort_keys=True))
        return 0 if table["complete"] else 1

    def line(record):
        extra = record.get("folder") or record.get("error") or f"at step {record.get('step')}"
        print(
            f"{record['status']:11s} {record['scene']}/{record['definition']}  ({record['seconds']:.1f} s)  {extra}",
            flush=True,
        )

    outcome = run_batch(
        jobs,
        out_root=args.out,
        state_dir=args.state,
        workers=args.workers,
        checkpoint_every=args.checkpoint_every or None,
        time_limit=args.time_limit,
        progress=line,
    )
    counts = ", ".join(f"{n} {status}" for status, n in sorted(outcome["counts"].items()))
    print(f"{len(jobs)} jobs: {counts}; {outcome['skipped']} were already finished; {outcome['seconds']:.1f} s")
    unfinished = sum(n for status, n in outcome["counts"].items() if status not in (DONE, REFUSED))
    return 0 if unfinished == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
