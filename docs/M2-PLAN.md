# M2 plan: what can the coupling be? (T1, T2)

Stage S4 of `docs/WORK-PLAN.md`; tests T1 and T2 of `docs/05-scenarios-and-tests.md`. Proposed by Claude Opus 5.5 in Claude Code (sponsor @MostafaDerraz) on 2026-11-14, after M1 closed (#78). **The maintainer accepted the plan and all seven decisions of section 3 on 2026-11-14: decisions D43 to D49.**

## 1. The goal and the exit gate

M2 asks of every reading of the Composer coordinate: does a universe built on it conserve what it should, does it depend on the frame angle, and how strong may it be before the measured limits (T9.1) exclude it?

- **T1.1-T1.6:** six pairs of cosmic masses, each run under L0, the default exploratory law, and then every definition of the T2 grid.
- **T2:** the survival table. Each definition is classed as **A** (conserving, angle-free), **B** (conserving, angle-dependent) or **C** (non-conserving), with the largest strength λ that T9.1 allows. Only A and B go on to M3.
- **Exit gate (S4):** the survival table confirmed: every definition in class A, B or C, the table reproduced by a contributor of a different model or a person (T2.5), and the M2 review closed (M2-R).

## 2. What M1 already gives M2

Every engine piece T1 and T2 need is merged:
- **The laws:** all the readings L1a-L1e, and the options L2-L8.
- **Integrators:** Boris, and the adaptive step with the D29 guard.
- **Energy:** the Composer energy where D23 defines it.
- **Runs:** checkpoints and resume (for a long batch), the convergence helper, the standard plots, and the third-law utility (E4.5).
- **T0 passes on the current engine** (#90), so every M2 result is stamped `t0_status: pass`.

## 3. Decisions before T2 (all accepted: D43 to D49)

Each is a proposal with a recommendation. The pre-registration of T2 records the answers before any production run (D12).

| # | Question | Proposal |
|---|---|---|
| M2-1 | **What does "conserving" mean for the classifier?** E9 defines an energy only within D23's scope (L1a, radial, inertia `Re`); most of the grid has no defined energy. | A definition is **conserving** if, on T0.2 and on each of T1.1-T1.6: (a) mass and charge are exact (always true); (b) the momentum, with the walls' share, drifts by at most 1e-6; and (c) an energy is conserved. For (c), where an approved energy exists (D23), its drift is at most 1e-6. Elsewhere, a **closed-loop work test** decides: the work done by the force on a test body carried round a closed loop in the field of fixed sources is at most 1e-9 of the work along one side. A force that does work round a loop has no energy, so the definition is class C. Thresholds as in `team/RESULTS-STANDARD.md` where it has them. |
| M2-2 | **What does "angle-free" mean?** | Every result of the Real coordinate is identical (to 1e-9 relative) at the six frame angles of T0.3 (10, 45, 79.86, 90, 135, 170 degrees). |
| M2-3 | **L1e with directions other than radial** (18 of the 120). L1e has no Composer force to direct (#69 refuses them). | **Not applicable.** The grid has 102 definitions, and the table says why 18 cells are empty. |
| M2-4 | **The direction "torque" on point bodies.** It turns a spin only through a radius (D37), so on point bodies its 30 definitions move exactly as L0. | Give the T1 bodies their **physical radii** (Earth, moon, proton 0.84 fm; the electron stays a point), so that the torque acts. Report the spin as the observable. |
| M2-5 | **The strength for classification**, and the λ bound (T2.4). | Classify at **λ = 1**, as the test plan says. Then bisect λ on the L3 sweep, 1e-24 to 1, for the largest value T9.1 allows, with its uncertainty. |
| M2-6 | **Approve the L1a plug-in and what "radial" means** (questions 2 and 3; proposed answers on the board). L1a is in the grid. | Approve both as proposed. |
| M2-7 | **Where is a body under L1e?** (question 25). T1 measures orbits. | Measure at the Real position `x`, and report `xc` and the combination `x − s²(q/m) xc` beside it, so the choice can change later without a re-run. |

Not needed for M2, and can wait for M3: question 6 (`c` and `e` in scaled units; M2's pairs use no second field), and question 7 (how to judge T0.1).

## 4. The tasks, in order

| Phase | Task | What it delivers | Proposed owner | Needs | Effort (plan) |
|---|---|---|---|---|---|
| 0 | Decisions M2-1 to M2-7 | recorded in `team/DECISIONS.md` | maintainer | this plan | one session |
| 1 | **T1.1-T1.6 scenarios** | `scenarios/T1.x/`: settings for each pair, scripts that run any `LawSpec`, pre-registration of the six questions; first runs under L0 and the default exploratory law | Claude Opus 5.5 in OpenCode | M2-4, M2-7 | ~15 h |
| 1 | **T2.1 generator** | the 102 (or 120) definitions as `LawSpec` files, each with a stable ID (for example `L1b+perp-Sh`) | Claude Opus 5.5 in Claude Code | M2-3 | 2 h |
| 1 | **T2.2 batch runner** | every definition on T0.2 and T1.1-T1.6 (714 or 840 runs), resumable with E7's checkpoints, parallel in CI or Codespaces; a first timing run to size the budget | Claude Fable 5.2 | T2.1, T1 | 6 h |
| 1 | **T2.3 classifier** | the pre-registered thresholds of M2-1 and M2-2; the closed-loop work test; class A, B or C for each definition, with its numbers | Claude Opus 5.5 in Claude Code | M2-1, M2-2 | 4 h |
| 2 | **The production runs** | all runs stamped, plus the angle sweep for definitions that conserve | batch runner (compute) | phase 1 | depends on the timing run |
| 2 | **T2.4 λ bounds** | the largest λ that T9.1 allows, for each definition, by bisection | Claude Fable 5.2 | T2.2 | 4 h |
| 3 | **T2.5 survival table** | the report, a figure, rows in `team/FINDINGS.md`, and a journal entry (D35) | Claude Opus 5.5 in Claude Code | T2.3, T2.4 | 4 h |
| 3 | **Reproduction of the table** | the table re-derived independently | **GPT-6** (a different model and provider, D42) or a person | T2.5 | — |
| 3 | **M2-R** | the milestone review: which definitions continue to M3 | maintainer closes it | all of the above | 2 h |

The owners are proposals, following the split that worked for M1 (D40). Each agent holds one claim at a time (D15), and every result needs a reproduction by a different model (D14, D27).

## 5. Risks

- **The compute budget is unknown** until the first timing run. The pairs are small (2 bodies), but T0.2 has 200 bodies and the angle sweep multiplies the count by 6. If it is too large: run the sweep only on definitions that pass the conservation test, and use a shorter T0.2.
- **The classifier is the delicate part.** A loose threshold puts non-conserving readings into A or B; a tight one fails good readings on rounding. Thresholds are pre-registered (M2-1) and the table reports the numbers, not only the classes.
- **At λ = 1 most definitions are excluded by T9.1**, as L1a is. That is the expected result, not a failure. The useful content of the table is the λ bounds (T2.4) and which readings conserve at all.
- **Results by the same model** do not count as confirmed (D27). The table needs GPT-6, Claude Fable 5.2 for results produced by Claude Opus 5.5 (and the reverse), or a person.
