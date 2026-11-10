# MVP plan: the first working virtual universe

**Goal of the MVP (version 0.1.0):** anyone can type one command and watch two cosmic masses move under the reference law L0 and under one reading of the Composer coordinate (L1a). The engine must pass the consistency gate T0 and reproduce the measured limits T9.1 and T9.3, with stamped results and figures.

This is a cut of `docs/WORK-PLAN.md`: the smallest set of tasks that gives a working, trustworthy universe. Everything else in the work plan comes after version 0.1.0.

> **Status board of the MVP:** the checklist in section 4. The agent updates it after every pull request, so that any new session can resume from here.

---

## 1. Inside the MVP, and outside it

| Inside (0.1.0) | Outside (after 0.1.0) |
|---|---|
| algebra core: constants, frame, General Numbers (scalar and arrays), readings, matter presets | `B_cal` and moving bodies (E2.2, Boris, L6) |
| state, minimal settings loader, law interface | laws L1b-L1e and the modifiers L2-L8 |
| L0, and L1a (radial, sign +, lambda as a setting) | populations, inflow, heat bath (E5.4, E6b) |
| direct-sum static field `g`, `e` with softening; uniform backgrounds | fast solver, super-particles (E2.5, E10) |
| velocity Verlet with an adaptive step; reflecting walls | T1 (except T1.4), T2-T8 |
| runner, `stamp.json`, `summary.json`, `timeseries.csv`, command line | the interactive web universe (UI1-UI4) |
| conservation diagnostics; energy under L0 | Composer energy for L1-L8 (needs an approved derivation, E9) |
| gates T0.1-T0.5 (T0 fast run in CI) | |
| T9.1 (neutron 9.0e-19 e) and T9.3 (1.14e11 V/m) | T9.2, T9.4 |
| demo: T1.4 hydrogen-like pair under L0 and L1a (labelled exploratory), with figures | |
| CI runs the tests (P2); release notes; tag proposal v0.1.0 | |

---

## 2. Decision D20: the MVP sprint (stacked pull requests)

For the MVP only:

- **One agent may hold the whole sequence below** as one claim (an exception to D15).
- **Stacked pull requests:** each pull request's branch is cut from the previous one in the sequence, and its **base** on GitHub is that previous branch, not `main`. So the agent never waits for a merge, and the maintainer later reviews and merges them **in order** (MVP-1, then MVP-2, ...). After a merge, GitHub moves the next pull request's base to `main`.
- **Never merge, never push to `main`.** Every pull request is reviewed separately (D14, D17).
- Every pull request still follows every other rule: tests, architecture (D18), log, notebook, labels.

---

## 3. The sequence

Each row is one pull request, of about 400 changed lines at most (tests included). Branch names are given in full. Task details and acceptance criteria are in `docs/WORK-PLAN.md`.

| # | Pull request | Tasks (WORK-PLAN) | Branch, cut from | Done when |
|---|---|---|---|---|
| MVP-1 | constants, frame, General Number | E1.1, E1.2 | `E1-algebra-core` from `main` | Rule of the Square exact at 50 digits; inverse; tests pass |
| MVP-2 | arrays, readings, matter, compendium tests | E1.3, E1.4, E1.5, E1.6 | `E1b-arrays-readings` from MVP-1 | arrays agree with scalars to 1e-15; `mu E_cal` product identity; `Re(mu_a mu_b)` angle-free |
| MVP-3 | CI runs pytest | P2 | `P2-ci-tests` from **`main`** (independent) | a failing test turns the pull request red |
| MVP-4 | state, settings, law interface, L0 | E5.1, E5.3 (minimal), E4.1, E4.2 | `E4-l0` from MVP-2 | two bodies: L0 force = `G_N m m'/r^2 - K q q'/r^2` to 1e-15 |
| MVP-5 | field solver (static) and backgrounds | E2.1, E2.3, E2.4 | `E2-fields` from MVP-4 | T0.5 analytic check within 1e-6 as a unit test |
| MVP-6 | integrator and walls | E3.1, E3.3, E6a (reflecting only) | `E3-integrator` from MVP-5 | Kepler 1000 periods: energy drift < 1e-6 |
| MVP-7 | runner, stamp, CLI, diagnostics, L0 energy | E7.1-E7.3, E8.1, E9 (L0 part only) | `E7-runner` from MVP-6 | `python -m cosmicmass run <settings>` writes valid stamped results |
| MVP-8 | the consistency gate | T0.1-T0.5, G1 | `T0-gate` from MVP-7 | all T0 pass criteria met; fast T0 in CI |
| MVP-9 | L1a and the measured limits | E4.3 (L1a only), T9.1, T9.3 | `T9-limits` from MVP-8 | 9.0e-19 e and 1.14e11 V/m reproduced |
| MVP-10 | demo and release notes | T1.4 (L0 and L1a, exploratory), README "Try it", `CHANGELOG.md` for 0.1.0 | `MVP-demo` from MVP-9 | one command reproduces the demo figures; the tag is proposed, not created |

Rules for every row:

1. Pre-register in `team/notes/<ID>.md` before any production run (MVP-8 to MVP-10).
2. Update `ARCHITECTURE.md` (D18): registry status, remove `:::planned`, contracts if they changed, change log.
3. Log entry in `team/LOG.md`, with the `Architecture:` line.
4. Update the board (`docs/OPEN-SCENARIOS.md`) and the checklist below.
5. Run `python tools/check_architecture.py --base origin/<base branch>` and `pytest` before pushing.

---

## 4. MVP checklist (updated by the agent after every pull request)

| # | State | Pull request | Notes |
|---|---|---|---|
| MVP-1 | PR open | [#2](https://github.com/MostafaDerraz/cosmic-mass-lab/pull/2) | E1.1 and E1.2 done 2026-11-10: 28 tests pass |
| MVP-2 | PR open | [#3](https://github.com/MostafaDerraz/cosmic-mass-lab/pull/3) | E1.3 to E1.6 done 2026-11-10: 64 tests pass |
| MVP-3 | PR open | [#4](https://github.com/MostafaDerraz/cosmic-mass-lab/pull/4) | P2 done 2026-11-10 on its own branch `P2-ci-tests` (base `main`); its notebook and log entry are on that branch |
| MVP-4 | PR open | [#5](https://github.com/MostafaDerraz/cosmic-mass-lab/pull/5) | E5.1, minimal E5.3, E4.1, E4.2 done 2026-11-10: 115 tests pass |
| MVP-5 | PR open | [#6](https://github.com/MostafaDerraz/cosmic-mass-lab/pull/6) | E2.1, E2.3, E2.4 done 2026-11-10: 128 tests pass |
| MVP-6 | PR open | [#7](https://github.com/MostafaDerraz/cosmic-mass-lab/pull/7) | E3.1, E3.3, E6a (reflecting) done 2026-11-10: 139 tests pass, 1 slow test run by hand; three assumptions for the maintainer in `team/notes/MVP-6.md` |
| MVP-7 | PR open | [#8](https://github.com/MostafaDerraz/cosmic-mass-lab/pull/8) | E7.1 to E7.3, E8.1, E9 (L0 part) done 2026-11-10: 152 tests pass; no checkpoints yet |
| MVP-8 | PR open | [#9](https://github.com/MostafaDerraz/cosmic-mass-lab/pull/9) | T0.1 to T0.5 pass on engine commit `91ba4ae` (status `done`, to be reproduced); fast gate in `engine/gates/`; the four questions were answered by the maintainer (D21) |
| MVP-9 | PR open | [#10](https://github.com/MostafaDerraz/cosmic-mass-lab/pull/10) | L1a (hypothesis); T9.1 gives -9.009e-19 e, T9.3 gives 1.1380e11 V/m, on engine commit `2baf36f` (status `done`, to be reproduced) |
| MVP-10 | PR open | [#11](https://github.com/MostafaDerraz/cosmic-mass-lab/pull/11) | demo `python scenarios/T1.4/run.py` (exploratory), README "Try it", `CHANGELOG.md`, version 0.1.0; the tag `v0.1.0` is proposed, not created |

States: `open`, `in progress`, `PR open` (with link), `merged`, `blocked` (with reason).

---

## 5. When to stop and wait for the maintainer

Stop, write the handoff in the notebook and the log, and update the checklist, if:

- a verification script or a T0 criterion fails and the cause is not a bug in your own new code;
- a formula is needed that is not in the compendium (it is a derivation for the physics lead, `AGENTS.md`, rule 7);
- an interface contract of `ARCHITECTURE.md`, section 10, must change in a way that breaks earlier pull requests;
- your session is ending or its context is getting long. Finish the current pull request cleanly first.

**Resuming:** a new session reads this file (section 4), the newest entries of `team/LOG.md`, and the notebook of the open item, then continues with the first item that is not `PR open` or `merged`.

---

## 6. Review and merge (maintainer)

Review and merge in order: MVP-1, MVP-2, MVP-3 (independent), then MVP-4 to MVP-10. An independent review (another model or a human) runs the tests on a clean checkout before each merge (D14). After MVP-10 is merged, create the tag `v0.1.0` and a GitHub release from `CHANGELOG.md`.
