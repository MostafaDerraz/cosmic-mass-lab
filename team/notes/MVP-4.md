# MVP-4: state, settings, law interface, L0 (E5.1, E5.3 minimal, E4.1, E4.2)

- **Owner:** Claude chat session on claude.ai / Claude Fable 5.2 (claude-fable-5-2), sponsor @MostafaDerraz
- **Buddy (newcomers):** none named
- **Claimed:** 2026-11-10, under decision D20 (row MVP-4 of `docs/MVP-PLAN.md`); no claim issue
- **Depends on:** MVP-2 (pull request #3, open)
- **Branch:** E4-l0 (cut from `E1b-arrays-readings`)

## Goal
The data and the reference law: a `State` of bodies, a minimal settings loader, the `Law` interface, and L0. Pass criterion (MVP plan): two bodies feel the L0 force `G_N m m'/r^2 - K q q'/r^2` to 1e-15.

## Pre-registration (commit before any production run)
Not applicable: engine components, no production run, no finding.

## Plan
1. `fields/sample.py`: `FieldSample` and the `FieldSolver` interface (needed by the `Law` signature).
2. `state.py` (E5.1), `laws/base.py` (E4.1), `laws/l0.py` (E4.2), `settings.py` (E5.3, minimal).
3. Tests, architecture, log, board, checklist.

## Settings
No run. Environment: Python 3.12.3 and the pins of `requirements.txt`.

## Progress (newest first)
- 2026-11-10: session goal: deliver MVP-4. Wrote the five modules and 51 tests (`test_state_laws.py`, `test_settings.py`). `pytest`: 115 passed, coverage 100% of `engine/cosmicmass`. `python tools/check_architecture.py --base origin/E1b-arrays-readings` passes.

## Results
Mathematics (compendium, section "Cosmic Mass"), checked by the tests:
- **Two bodies under L0:** with the field of a source at rest, `E_cal = -G mu_s e_r / r^2`, the force is `(G_N m m' - K q q') / r^2` to 1e-15 relative, for the 100 pairs of presets, in 2D and 3D, at 4 frame angles, with equal and opposite forces on the two bodies.
- **Angle-free:** the L0 force on moving bodies is the same array, bit for bit, at the 33 frame angles from 10 to 170 degrees.
- **Newton and Lorentz:** `Re[mu (E_cal + v x B_cal)] = m (g + v x B_g) + q (E + v x B)` with `E = -s^2 e`, `B = -s^2 b`, on random fields in 2D and 3D, to 1e-12.

Choices made:
- **New module `fields/sample.py`.** `Law.forces(state, fields, frame)` needs the type `FieldSample` before the solver exists (MVP-5). It is registered in `ARCHITECTURE.md` with its node and edges.
- **`Law` is a base class, not a Protocol**, with the members of the contract (`id`, `spec`, `forces`, `inertia`). Laws register themselves by reading (`@register`), and `make_law(spec)` builds them.
- **`ForcePair(real, composer)`**: `composer` is the force that the law reads from the Composer coordinate; it is zero under L0. The force applied is `real + composer`.
- **`LawSpec` uses the keys of the stamp** (`id`, `sign`, `lambda`, `direction`, `terms`, `B_cal`, `inertia`, `range`). `moving` (L6) is a string `off`, `on`, `only`, because L6 has three options (the class diagram had a Boolean). `lambda` may be a number or `"cos_phi_bar"` (L3).
- **L0 takes no option:** an L0 spec with any L2 to L8 option is refused.
- **L0 includes `v x B_cal`**, as the formula does, although the MVP solver gives no second field yet. In 2D the second field points out of the plane and has shape `(N,)`.
- **Inertia:** `Re(mu) = m`. The readings `Sh` and `abs` (L7) raise `NotImplementedError` until E4.4.
- **Settings (minimal):** JSON with `scenario`, `universe`, `law`, `bodies`, `numerics`, `kind`, `preregistered`, `producer`. Bodies are a preset name (SI only, optional charge) or `m` and `q`. `universe.boundary` and `universe.background` are accepted and stored but not acted on yet (MVP-5, MVP-6). Default frame angle: 79.86 degrees in SI, 90 degrees in scaled units.

## Convergence
Not applicable.

## Open questions
- `Universe` (E5.2) does not exist yet; the loader keeps the universe as a validated dictionary. E5.2 and the P4 schemas should replace it after the MVP.
- Neither the law nor the state treats absorbed or fixed bodies specially: the integrator (MVP-6) must skip them, and the solver (MVP-5) must not count absorbed bodies as sources.

## Handoff (fill in if you stop before finishing)
- State: delivered as a pull request, waiting for review.
- What works: everything listed under Results.
- What does not: nothing known.
- Exact next step: MVP-5, `fields/direct.py`, on `E2-fields` cut from this branch.
