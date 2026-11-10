# Q1 to Q3 and T9.4: open questions turned into scenarios

- **Owner:** Claude chat session on claude.ai / Claude Fable 5.2 (claude-fable-5-2), sponsor @MostafaDerraz
- **Buddy (newcomers):** none named
- **Claimed:** 2026-11-10, on the sponsor's instruction in the sponsoring chat ("all questions turn it as scenarios"); no claim issue
- **Depends on:** T0 (`confirmed`), T9.1 and T9.3 (on `main`)
- **Branch:** E4.5-T9.4 (cut from `main` after the merge of pull request #14)

## Goal
The MVP left questions for the maintainer. Following the example of the D21 scenarios, each question that a run can answer is run, with its alternatives, instead of being settled by opinion. Three questions qualify, plus the open test T9.4. The utility E4.5 (Newton's third law check) is written for Q1.

Questions that a run cannot answer stay questions, listed under "Open questions".

## Pre-registration (commit before any production run)
Committed before any code of these scenarios was written or run. No exploratory run was made. One thing I knew beforehand: the demo T1.4 (exploratory) showed the L0 energy swinging by 44% and the momentum constant to 3e-15 under L1a.

**A correction that motivates Q1.** In `team/notes/T1.4.md` (and in the description of pull request #11) I wrote that under L1a "the Coupling term does not obey Newton's third law". Reasoning it through again, I now expect the opposite for L1a with all its terms, and Q1 tests it. If Q1 confirms, that sentence was wrong and is corrected in the same pull request.

### T9.4: particle against antiparticle (gate, hypothesis L1a)
- **Question (docs/05, T9.4):** L1a, lambda = 1; a neutral Earth; a proton, then an antiproton (`m_p - e G`), at 7e6 m: what is the sign of the extra push from `s q g`?
- **Settings:** SI, 2D, `phi_bar` = 79.86 deg, softening 0, Earth fixed at the origin, each body alone at (7e6 m, 0), at rest; forces evaluated once.
- **Prediction (analytic):** the Earth is neutral, so `e = 0` and the Composer force is `s q g` alone. For the proton it points along `g` (towards the Earth); for the antiproton it is exactly opposite. Size `s e g` = 1.16043e10 x 1.602177e-19 x 8.13447 = 1.5124e-8 N, which is `s e / m_p` = 1.11e18 times the weight. The Real force (the weight) is the same for both.
- **Pass:** the two Composer forces are opposite (their sum is zero to 1e-12 of their size), the proton's points towards the Earth, and each equals `s q g` to 1e-12. Kind: measurement check in the sense of T9 (it reproduces the book's statement "particle against antiparticle"; no measurement is compared here).

### Q1: does L1a obey Newton's third law?
- **Question:** for a pair of bodies, is the force on a from b equal and opposite to the force on b from a, under L1a with each choice of terms (option L5)?
- **Law(s):** L0 (reference); L1a (hypothesis), lambda = 1, sign +, with terms `all`, `s q g`, `-m E/s`, `-2 cos q E`.
- **Settings:** SI, 2D, 79.86 deg, softening 0, forces evaluated once. Pairs at 1e-10 m: proton and electron; proton and neutron; two protons; and Earth with a proton at 7e6 m.
- **Prediction and why:** the field of a body at rest is `E_cal = -G mu_b e_r / r^2`, so `mu_a E_cal = -G (mu_a mu_b) e_r / r^2`, and the product of two General Numbers commutes. Therefore **with all terms the pair forces are equal and opposite** (residual at rounding, below 1e-12 of the force), in both coordinates. Term by term: the charge share is symmetric (`q_a q_b`), so it also obeys the law alone. The two halves of the Coupling term do not, taken alone: `s q g` gives `-s G q_a m_b` on a and `-s G q_b m_a` on b, so the residual is `s G |q_a m_b - q_b m_a| / r^2`; `-m E/s` alone gives the same residual with the other sign. They cancel only together. For two protons every option is symmetric (`q_a m_b = q_b m_a`).
- **Analysis:** `third_law_residual` (E4.5) gives the largest `|F_ab + F_ba|` over the pairs and the largest force; their ratio is compared with the prediction for each pair and each option.

### Q2: which energy does L1a conserve?
- **Question:** the L0 energy is not conserved under L1a (demo: 44%). Is there a quantity that is? This is the question that blocks the survival table T2 and the Composer part of E9.
- **Candidates:**
  - A: the L0 energy, `kinetic - G sum Re(mu_a mu_b) / r` (what the engine reports today);
  - B: `kinetic - G sum [Re(mu_a mu_b) + lambda s Com(mu_a mu_b)] / r`, the same sum read with the reading of L1a. It is the work of the pair force predicted in Q1 (central, symmetric, `1/r^2`).
- **Law(s):** L1a (hypothesis), all terms, inertia `Re(mu)`.
- **Settings:** (i) the demo pair: proton and electron on the Bohr circle of L0, SI, 79.86 deg, lambda = 1, softening 0, 5 L0 periods; (ii) the same with lambda = 0.3 and sign -; (iii) the same at `phi_bar` = 120 deg (where the charge share adds to the pull); (iv) 12 bodies in scaled units (masses 1, charges drawn with seed 7 from a normal law of width 0.5, positions uniform in a square of side 4, at rest, 79.86 deg, lambda = 0.2), 3000 adaptive steps, softening of D10 used in both the force and the candidate energies.
- **Prediction:** A swings by tens of percent in (i) to (iii) and visibly in (iv). B is conserved like any energy under velocity Verlet: end drift below 1e-6 in all four, largest excursion below 1e-5 (the orbits are eccentric, compare D21.3).
- **Analysis:** end drift and largest excursion of A and of B over each run (D21 as refined: both numbers). B is computed by the scenario script from the stored positions and velocities; it is **not** added to the engine.
- **What a pass would and would not mean:** it would show numerically that B is the conserved energy of L1a with all terms for bodies at rest or moving slowly without `B_cal`. It would not make B an approved formula: B is not in the compendium, so adopting it in `energy.py` remains a derivation for the physics lead (`AGENTS.md`, rule 7). It says nothing about the single-term options, which Q1 predicts to be non-reciprocal.

### Q3: does the choice of `K` (CODATA 2018 or 2022) change any verdict?
- **Question:** `docs/05` quotes `K = 8.9875517923e9` (CODATA 2018); the engine uses `1/(4 pi eps0)` with the CODATA 2022 `eps0`, 8.9875517862e9 (`team/notes/E1.md`, open question 1). Does it matter?
- **Settings:** the quantities of T0.1b (analytic period at the Bohr radius), T0.1c (closed-form separation), T9.1 (`q_eff`) and T9.3 (`s g`), each computed with both values of `K`, through the engine's own formulas (`Constants`, `Frame`, `DirectSum`, the laws).
- **Prediction:** `K` changes by 6.8e-10, so `s` by 3.4e-10, the Bohr period by 3.4e-10, `q_eff` and `s g` by 3.4e-10. Every one is thousands of times below the tightest criterion (1e-6). **No verdict changes.**
- **Pass:** every relative change is below 1e-8 and no pass or fail flips.

- **Changes to the plan (with date and reason):**
  - 2026-11-10, Q1, after the first production run (engine commit `e1a3165`): for the pair proton and electron, the predicted residual of the single-term options (2.1e-26 N) is 18 orders of magnitude below the Coulomb force (2.3e-8 N), so it is below the rounding of the total force, and the planned analysis (the residual of the total force) could neither confirm nor refute it: it returned exactly 0. My script nevertheless marked those two cases "pass", because its tolerance was tied to the size of the force. That was a flaw of the script, not a result. I added the residual of the **Composer forces alone** to every case, made "pass" require that it equals the prediction, and re-ran Q1 on the new commit. The first Q1 output was discarded uncommitted. The settings and the predictions are unchanged.
  - 2026-11-10, T9.4: I ran `t9_4` once while writing it, before committing its code (a check that the function worked). It gave the numbers of the prediction. The stamped run is the one made on the committed code.

## Plan
1. This pre-registration (first commit).
2. `third_law_residual` in `laws/base.py` (E4.5) with tests; `t9_4` in `gates/t9.py`; `scenarios/T9.4/run.py`; `scenarios/Q/run.py`.
3. Run on the committed code; commit the stamped results; correct the sentence in `team/notes/T1.4.md` if Q1 confirms.

## Settings
See the pre-registration; stamped in the results.

## Progress (newest first)

## Results
(to be filled after the production runs)

## Convergence
(to be filled)

## Open questions

## Handoff (fill in if you stop before finishing)
- State:
- What works:
- What does not:
- Exact next step:
