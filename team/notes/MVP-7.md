# MVP-7: runner, stamp, command line, diagnostics, L0 energy (E7.1 to E7.3, E8.1, E9 L0 part)

- **Owner:** Claude chat session on claude.ai / Claude Fable 5.2 (claude-fable-5-2), sponsor @MostafaDerraz
- **Buddy (newcomers):** none named
- **Claimed:** 2026-11-10, under decision D20 (row MVP-7 of `docs/MVP-PLAN.md`); no claim issue
- **Depends on:** MVP-6 (pull request #7, open)
- **Branch:** E7-runner (cut from `E3-integrator`)

## Goal
One command turns a settings file into stamped results. Pass criterion (MVP plan): `python -m cosmicmass run <settings>` writes valid stamped results.

## Pre-registration (commit before any production run)
Not applicable: engine components. The runs in the tests write to temporary folders, not to `results/`. Nothing here is a scenario result.

## Plan
1. `energy.py` (E9, L0 part) and `diagnostics.py` (E8.1, conservation part).
2. `stamp.py` (E7.2) and `runner.py` (E7.1, E7.3), `__main__.py`.
3. Tests, architecture, log, board, checklist.

## Settings
Unit tests: a circular orbit of masses 1 and 1e-3 in scaled units; a body bouncing on a reflecting floor in a uniform `g`; an electron around a fixed proton in SI. Environment: Python 3.12.3 and the pins of `requirements.txt`.

## Progress (newest first)
- 2026-11-10: session goal: deliver MVP-7. Wrote the five modules and 13 tests (`test_runner.py`). `pytest`: 152 passed, 1 deselected; coverage 99% (`__main__.py`, three lines, is exercised only through a subprocess, which the coverage tool does not see). `python tools/check_architecture.py --base origin/E3-integrator` passes. While testing, found and fixed a fault of the adaptive step in MVP-6 (a lone body in a background field never passed its turning point); the fix is on `E3-integrator`.

## Results
Software facts and checks against **mathematics** under L0; no finding of the virtual universe.
- `python -m cosmicmass run settings.json --out <folder>` writes `<folder>/<scenario>/<run-id>/` with `stamp.json`, `summary.json`, `timeseries.csv`. The stamp has exactly the fields of `team/RESULTS-STANDARD.md`, section 2. The run id is `YYYYMMDD-HHMMSS-<commit>` in UTC; the commit carries `-dirty` when tracked files are changed.
- **Determinism:** two runs of the same settings give the same `summary.json`, byte for byte (the summary holds no clock time and no run id).
- **Conservation in the summary:** on the unit orbit, mass and charge drift exactly 0, energy and momentum drift below 1e-6. A body bouncing twice on a reflecting floor in a uniform `g` keeps its energy to 1e-6; so does an electron around a fixed proton in SI over 3000 steps.
- **L0 energy:** kinetic energy plus `-G_N sum Re(mu_a mu_b) / r`, which equals `-G m m'/r + K q q'/r` to 1e-14 for a proton and an electron, and is the same number at every frame angle.

Choices made:
- **Energy (E9, L0 part).** The interaction energy is the work of the force of the compendium ("Cosmic Interaction Law and its reading"), `U = -G_N sum Re(mu_a mu_b) / r`; in a uniform background, `-Re[mu E_cal] . x`. With softening, `r` becomes `sqrt(r^2 + eps^2)`, the potential of the Plummer force of the solver. **The Composer energy is reported as `null`** (`energy_com`): it is not in the compendium and needs a derivation approved by the physics lead (rule 7). I did not invent one.
- **The energy reported is always the L0 energy**, whatever the law of the run. Under a reading L1 to L8 its drift is then a measurement of how far that reading departs from L0 conservation, not a conserved quantity of that reading.
- **Momentum drift** is the change of the total momentum divided by the sum of the bodies' `m |v|`. Walls and fixed sources take momentum, so it is expected to drift with them.
- **Mass and charge** count every body, absorbed ones included.
- **`t0_status` is "not run"** in every stamp until the gate exists (MVP-8, G1).
- **Results folders are never overwritten**: writing into an existing run folder is an error.
- **Time series:** one row every `numerics.record_every` steps (default: about 1000 rows when `steps` is given, else every 100 steps), with positions and velocities of each body when there are at most 16 bodies.
- **New module `__main__.py`** (three lines) for `python -m cosmicmass`; the command line itself is `runner.main`.

## Convergence
Not applicable.

## Open questions
- **Not done from E7.1:** checkpoints and resuming, the progress line, and repeating a step with half the time step when a drift passes its tolerance (`ARCHITECTURE.md`, section 6). **Not done from E7.2:** validation against the P4 schemas (P4 is open) and the T0 status lookup.
- **Not done from E8.1:** temperature and pressure on the walls.
- The settings format has no schema yet; `settings.py` is its only definition.

## Handoff (fill in if you stop before finishing)
The session that held the MVP sequence stops here (its context is long; `docs/MVP-PLAN.md`, section 5). MVP-1 to MVP-7 are delivered as pull requests #2 to #8; MVP-8 to MVP-10 are not started.
- State: MVP-7 delivered as a pull request, waiting for review. Nothing is half done.
- What works: one command from settings to stamped results, under L0, with conservation diagnostics.
- What does not: see Open questions.
- Exact next step: MVP-8 (T0.1 to T0.5 and G1) on `T0-gate`, cut from `E7-runner`. Before pre-registering it, get the maintainer's answers to four questions, because each one changes a T0 pass criterion or setup:
  1. Is "2D" planar motion with a `1/r^2` field? (`team/notes/MVP-5.md`, assumption 1)
  2. May the two-body gate tests run with softening 0? With D10's softening the closest initial force is off by 1.5e-6, above their 1e-6. (`team/notes/MVP-5.md`, assumption 2)
  3. Does "energy drift below 1e-6" mean the change at the end of the run, or the largest excursion on the way? (`team/notes/MVP-6.md`, assumption 1)
  4. Are the floor of the speed and the displacement limit of the adaptive step accepted? (`team/notes/MVP-6.md`, assumptions 2 and 3)
  Then: write `scenarios/T0.1/` to `scenarios/T0.5/` (a `settings.json` and a `run.py` each, with a registry row), pre-register in `team/notes/T0.md`, run, and add `engine/gates/` with a fast T0 for CI. Note for G1: a step costs about 120 microseconds for two bodies, so the CI gate needs short runs.
