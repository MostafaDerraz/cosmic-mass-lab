# MVP-5: field solver (static) and backgrounds (E2.1, E2.3, E2.4)

- **Owner:** Claude chat session on claude.ai / Claude Fable 5.2 (claude-fable-5-2), sponsor @MostafaDerraz
- **Buddy (newcomers):** none named
- **Claimed:** 2026-11-10, under decision D20 (row MVP-5 of `docs/MVP-PLAN.md`); no claim issue
- **Depends on:** MVP-4 (pull request #5, open)
- **Branch:** E2-fields (cut from `E4-l0`)

## Goal
The static cosmic field of point cosmic masses, with softening, no self-field and uniform backgrounds. Pass criterion (MVP plan): the T0.5 analytic check within 1e-6, as a unit test.

## Pre-registration (commit before any production run)
Not applicable: engine component, no production run, no finding. The T0.5 check here is a unit test against an analytic formula; the stamped T0.5 run belongs to MVP-8.

## Plan
1. `fields/direct.py`: `DirectSum` with `sample(state, points)` and `at_bodies(state)`.
2. Backgrounds in the solver and in the settings loader.
3. Tests, architecture, log, board, checklist.

## Settings
No run. Environment: Python 3.12.3 and the pins of `requirements.txt`.

## Progress (newest first)
- 2026-11-10: session goal: deliver MVP-5. Wrote `fields/direct.py`, extended `settings.py` (backgrounds, boundary names) and 13 tests (`test_fields.py`). `pytest`: 128 passed, coverage 100%. `python tools/check_architecture.py --base origin/E4-l0` passes.

## Results
Mathematics (compendium, "Cosmic Field Law and law of motion": the field of a source at rest is `E_cal = -G mu_s e_r / r^2`), checked by the tests:
- **T0.5 as a unit test:** one proton gives `E = -s^2 e = K e / r^2` to 1e-12 relative (asked: 1e-6), in 2D and 3D, at radii from 1e-11 m to 1 m; its mass gives `g = -G_N M / r^2`.
- The field at each of 30 bodies equals the sum over the others taken one by one (1e-10); there is no self-field and no infinity.
- Under L0 the internal forces of a cloud of 60 bodies sum to zero (Newton's third law) to 1e-12 of their size.
- Uniform backgrounds `g`, `e`, `B_g`, `b` appear unchanged at every point.
- Software fact: 1,000 bodies in 2D, one field evaluation in well under 0.5 s on the machine used (the test fails above 0.5 s).

Choices made, each a **labelled assumption** for the physics lead:
1. **2D means planar motion of point bodies of three-dimensional space**: the field falls as `1/r^2` in 2D too. Reason: T0.1 asks for the 27.3-day month and the Bohr period in 2D. A field law solved in a true two-dimensional space would fall as `1/r`; that would be a different universe and needs a decision.
2. **Softening has the Plummer form** `r / (r^2 + eps^2)^(3/2)`. D10 fixes the length (1e-3 of the smallest initial separation) but not the form. At the separation `r` it changes the force by about `1.5 (eps/r)^2`: **1.5e-6 at the closest initial pair**. This is above the 1e-6 of T0.1(c) and T0.5 at that one distance, so two-body gate tests should set `softening` to 0 explicitly (the settings allow it). To raise with the maintainer before MVP-8.
3. **Static part only:** `B_g` and `b` come from the background alone. E2.2 is outside the MVP.
4. **Sources are the alive bodies**, fixed ones included (U6: a pinned source is a body with `fixed = true`). Absorbed bodies do not act.
5. **Backgrounds are given in the theory's coordinates** (`g`, `e`, `B_g`, `b`). A Maxwell field `E` is entered as `e = -E / s^2`.

## Convergence
Not applicable.

## Open questions
- Points 1 and 2 above.
- `at_bodies(state)` is an addition to the solver interface of section 10 (`sample(state, points)` cannot know which point is which body). The fast solver (E2.5) must provide it too.

## Handoff (fill in if you stop before finishing)
- State: delivered as a pull request, waiting for review.
- What works: everything listed under Results.
- What does not: nothing known.
- Exact next step: MVP-6, `integrators.py` and reflecting walls in `environment.py`, on `E3-integrator` cut from this branch.
