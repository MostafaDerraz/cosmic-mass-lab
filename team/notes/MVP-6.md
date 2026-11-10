# MVP-6: integrator and walls (E3.1, E3.3, E6a reflecting only)

- **Owner:** Claude chat session on claude.ai / Claude Fable 5.2 (claude-fable-5-2), sponsor @MostafaDerraz
- **Buddy (newcomers):** none named
- **Claimed:** 2026-11-10, under decision D20 (row MVP-6 of `docs/MVP-PLAN.md`); no claim issue
- **Depends on:** MVP-5 (pull request #6, open)
- **Branch:** E3-integrator (cut from `E2-fields`)

## Goal
Move the bodies: velocity Verlet with a `Law` and a `FieldSolver`, the adaptive step of D10, and reflecting walls. Pass criterion (MVP plan): Kepler orbit, 1000 periods, energy drift below 1e-6.

## Pre-registration (commit before any production run)
Not applicable: engine components. The orbits below are unit tests against analytic two-body results under the reference law L0, not scenario runs; nothing is written to `results/`. The stamped gate runs T0.1 to T0.5 belong to MVP-8.

## Plan
1. `integrators.py`: `StepContext`, `VelocityVerlet` (`next_dt`, `step`).
2. `environment.py`: `Environment` with reflecting walls.
3. Settings: `universe.box`, `numerics.eta`, `dt_min`, `dt_max`.
4. Tests, architecture, log, board, checklist.

## Settings
Unit tests in scaled units (`G_N = K = 1`) unless stated; two masses 1 and 1e-3 at separation 1; fixed step of 1/6283 of the period (the D10 step of a circular orbit), or the adaptive step with `eta = 1e-3`. Environment: Python 3.12.3 and the pins of `requirements.txt`.

## Progress (newest first)
- 2026-11-10: session goal: deliver MVP-6. Wrote `integrators.py`, `environment.py`, the settings fields and 13 tests (`test_integrators.py`, one marked `slow`). `pytest`: 139 passed, 1 deselected (the slow one), coverage 100%. The slow test was run once by hand: see Results. `python tools/check_architecture.py --base origin/E2-fields` passes.

## Results
All under the reference law L0; each is a check of the engine against **mathematics** (the two-body problem with the force `G_N Re(mu_a mu_b) / r^2` of the compendium). No finding of the virtual universe.

- **Kepler, 1000 periods (the pass criterion, test marked `slow`):** 6,283,001 fixed steps (6283 per period). Largest relative energy excursion during the run **3.9e-13**; energy change at the end 1.0e-13; measured period within 3.3e-7 of the analytic one. Run once by hand with `pytest -m slow` (12 min 43 s). It passes the criterion (drift below 1e-6, period within 0.1%).
- **Kepler, 5 periods (default test):** largest energy excursion below 1e-6, period within 0.1%.
- **Hydrogen-like pair in SI at the Bohr radius (the setup of T0.1b):** period 1.5198e-16 s within 0.1%, energy excursion below 1e-6, over 3 periods.
- **Two protons from rest in SI (the setup of T0.1c):** energy change below 1e-6 over 2000 adaptive steps; the centre of mass does not move.
- **Fixed-step Verlet** is time-reversible (5 charged bodies in 3D return to their start within 1e-11 after 300 steps forward and 300 back) and keeps the total momentum to 1e-13.
- **Eccentric orbit, e = 0.9, adaptive step with eta = 1e-3 (E3.3):** energy change after one period 7.3e-9; it grows by the same amount each period (measured over 3 periods: 7.3e-9, 1.5e-8, 2.2e-8). **Largest excursion on the way: 3.7e-6**, at each pericentre. With eta = 3e-4 the excursion is 3.3e-7.
- **Reflecting walls:** speeds are exactly unchanged by a reflection; a gas of 6 weakly interacting bodies in a box keeps its energy to 1e-6 over 6000 adaptive steps with wall hits.

Choices made, each a **labelled assumption**:
1. **"Energy drift" is read as the change of the energy at the end of the run.** For the circular orbit the largest excursion on the way is also below 1e-6. For e = 0.9 it is not (3.7e-6 at pericentre with D10's eta = 1e-3), although the energy returns to within 7.3e-9 per period. If the criterion of E3.3 means the largest excursion, D10's 1e-3 is not enough for e = 0.9 and about 3e-4 is needed. **For the maintainer to decide.**
2. **Adaptive step, bodies at rest.** D10's rule `|a| dt <= eta |v|` gives a zero step for a body at rest. The speed is floored at `REST_FRACTION = 0.1` of `sqrt(|a| l)`, the free-fall speed over the distance `l` to the nearest neighbour. D10 holds as written for every faster body (checked: the circular orbit gets exactly `eta |v| / |a|`).
3. **Adaptive step, displacement limit.** No body moves more than `eta` of the distance to its nearest neighbour per step. It only shortens the step. Without it, the D10 step of two bodies flying apart grows without limit (seen in the two-proton test: the step overflowed within 2000 steps).
4. **Reflecting walls turn back the velocity and do not move the body.** A body that has crossed a wall keeps its position and comes back on the next steps, so position, speed and energy are exactly unchanged by the reflection. Bodies overshoot a wall by at most one step. The impulse on each wall is returned for the pressure (E8.1).
5. **Velocity Verlet with velocity-dependent forces:** the force at the end of a step is computed with the half-step velocity. This is not accurate for `v x B_cal`; that needs the Boris push (E3.2, outside the MVP). Under the MVP there is no second field except a uniform background.
6. **Long tests are marked `slow`** and skipped by default (`pyproject.toml`); `pytest -m slow` runs them.

## Convergence
Step halving was not run as a study here (these are unit tests). Evidence of second order: the excursion of the e = 0.9 orbit falls from 3.7e-6 to 3.3e-7 when eta goes from 1e-3 to 3e-4 (a factor 11.2 for a step ratio 3.33, expected 11.1 at second order).

## Open questions
- Points 1 to 3 above, for the maintainer as physics lead.
- Speed: about 120 microseconds per step for two bodies. The 1000-period run takes 6.3 million steps. T0 in CI (G1) will need shorter runs or a faster two-body path.

## Handoff (fill in if you stop before finishing)
- State: delivered as a pull request, waiting for review.
- What works: everything listed under Results.
- What does not: the largest energy excursion of the e = 0.9 orbit is above 1e-6 with eta = 1e-3 (assumption 1).
- Exact next step: MVP-7 (`energy.py`, `diagnostics.py`, `stamp.py`, `runner.py`, the command line) on `E7-runner`, cut from this branch.
