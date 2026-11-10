# General note: the virtual universe, before we build it

## 1. What the simulator is for

It is a laboratory for one question the theory leaves open. **The product of two cosmic masses has two coordinates.** The real part is Newton plus Coulomb, which is measured and verified. The Composer coordinate (the Coupling term `m_a q_b + q_a m_b` and the charge share `2 s cos(phi_bar) q_a q_b`) has never been read.

The simulator lets us give the Composer coordinate a meaning (laws L1 to L8), switch it on in a universe built only from cosmic masses, and watch what happens. Everything stays in the language of the theory: `mu`, `G`, the Rule of the Square, `E_cal`, `B_cal`, and the readings `Re`, `Com`, `Sh`, `Ht`, `abs`.

## 2. What is solid, and what is a hypothesis

- **Solid, and used as engine inputs:**
  - every formula in `04-formulas-compendium`, checked again for this package (106 + 204 + 253 + 144 laws + 9,304 + 79 checks, all passing);
  - the paper's law L0 with the Cosmic Field Law.
- **Hypotheses, the subject of the experiments:**
  - every reading of the Composer coordinate (L1 to L8);
  - any physical claim of the book's Chapters 7, 8, 10 and 11.
- **Known failures of the physics, even under L0 (stated in the paper):**
  - LAGEOS frame dragging comes out at 1/4 of the measured value;
  - no light bending;
  - negative gravitational-wave energy;
  - no mass-to-energy conversion;
  - slow bodies only.

  The simulator works in the slow-body world, where these do not arise, and it must not be used to argue against them.

## 3. What the book already settled

The book's Chapter 10 has done careful work, and the simulator must start from it:

- **A force along the line is ruled out at full strength.** A reading of the Coupling term as a force along the line at full strength would make a neutron respond as if it carried 9.0e-19 e. Measured: (-0.2 +/- 0.8)e-21 e. So such a reading is limited to below about **0.002** of its natural size.
- **The charge share as a force is excluded.** It would bind hydrogen at 5.71 eV; measured: 13.598 eV.
- **Readings of other kinds stay open.** Measurement says nothing yet about the Composer coordinate as something other than a force along the line: a mass, a property of a pair, a phase, a spin, a motion in a second coordinate (L1e).

So the most interesting experiments are not "turn the coupling on at full strength and watch planets fly apart". That world is already excluded. They are:

1. **the readings of other kinds** (L1e motion in two coordinates, L4 torque or perpendicular, L7 inertia readings);
2. **weak couplings below the limit** (lambda < 0.002), which might accumulate in large, long-lived systems: planetary charge, rings, stellar plasma (T4).

An earlier rough estimate in our conversation, from electrons falling near the Earth, gave a much smaller limit. The book is right that this test is hidden by metal screening, so the **neutron limit is the reliable one**. Test T9 reproduces the book's numbers before anything else runs.

## 4. Observations from the mathematics (worth testing first)

- **The charge share vanishes at `phi_bar` = 90 deg.** That is exactly the angle of natural units (Planck or Stoney). The frame in which the theory is cleanest is also the frame in which the only angle-dependent term switches off (T3.2).
- **If the law of motion holds in both coordinates (L1e), every body accelerates by the field itself.** The charge-to-mass ratio drops out. Bodies would move in a two-coordinate space: the real coordinate with gravity, the Composer coordinate with the electric field. This is a different kind of universe, not excluded by the neutron limit, and worth seeing (T1, T2).
- **Self-charging (T4.3).** A massive body that collects charge from a surrounding population, under a radial coupling, settles at `Q = lambda M / s`. The real parts of two such bodies then shrink as `G M1 M2 (1 - lambda^2)`, plus coupling terms of order `lambda^2`. Measured planetary and solar charges limit lambda independently of the neutron. The simulator can show how fast and how fully such charging happens.
- **The diagonal `G m^2 = K q^2`** (`Re(mu1 mu2) = 0` for two equal bodies) is where the real force vanishes. The book links it to extremal black holes. It is a natural test case for every reading (T1.5, T8.3).

## 5. Design principles

1. **T0 is a gate.** The engine must reproduce the paper's law exactly before any coupling experiment counts.
2. **One switch per idea.** Every law option and every environment module can be turned on or off on its own, so any effect can be traced to its cause.
3. **Cosmic language only, in the engine and in the outputs.** Temperature and pressure are expressed through populations of `mu` and their momentum flux.
4. **Both coordinates are always reported**, `Re` and `Com`, for forces, energy, momentum and each `mu(t)`.
5. **SI and scaled units side by side.** Scaled units make effects visible, and SI keeps them honest.
6. **Every run is stamped** with its settings and the T0 and T9 status, so results can be compared and reproduced.
7. **Findings are recorded as findings.** A result that contradicts measurement is reported as excluded, not tuned away.

## 6. Risks and how the plan handles them

| Risk | Handling |
|---|---|
| Effects 10^21 times apart in SI cannot be integrated together | super-particles, scaled units, and separate runs per scale (T7.1) |
| Integrator error looks like new physics | conservation checks every step, and time-step halving (T7.4) |
| A reading breaks energy conservation | classified in T2 as class C, and reported, not hidden |
| 2D distorts the 1/r^2 laws | 3D runs for the key scenes (T7.2) |
| Over-interpreting a simulation | the report states what is mathematics, what is measured, and what is a finding of the virtual universe only |

## 7. Next steps and effort

| Step | Work | Estimated working time |
|---|---|---|
| 1 | engine and T0 | ~2 h |
| 2 | T9 (book limits) | ~1 h |
| 3 | T1, T2 (survival table) | 2-3 h |
| 4 | T3, T4 | 3-4 h |
| 5 | T5, T6, T7, T8 | 3-4 h |
| 6 | interactive universe | 3-4 h |

**Total:** about 14-18 hours of working time, in stages. Each stage ends with a findings note.
