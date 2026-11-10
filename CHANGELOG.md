# Changelog

## 0.1.0 (proposed; the tag is created by the maintainer after the MVP pull requests are merged)

The first working virtual universe of cosmic masses (`docs/MVP-PLAN.md`).

### What it does

- **Algebra of the compendium:** constants (CODATA 2022), the leaning frame, General Numbers as scalars (with an exact mode at any precision) and as arrays, the readings `Re`, `Com`, `Sh`, `Ht`, `abs`, the matter presets.
- **The reference law L0** (`Re(mu) dv/dt = Re[mu (E_cal + v x B_cal)]`): Newton plus Coulomb from one product, the same at every frame angle.
- **One reading of the Composer coordinate, L1a** (a hypothesis): force `Re(mu E_cal) + lambda s Com(mu E_cal)`, with the options sign, strength and terms.
- **The static cosmic field** of point bodies by direct sum, in 2D and 3D, with softening and uniform backgrounds.
- **Motion:** velocity Verlet with the adaptive step of D10 (as completed by D21), reflecting walls.
- **Runs:** `python -m cosmicmass run settings.json` writes `stamp.json`, `summary.json` and `timeseries.csv`, with conservation diagnostics and the L0 energy.
- **Gates:** T0.1 to T0.5 pass; T9.1 (a neutron under L1a responds as -9.009e-19 e) and T9.3 (1.138e11 V/m) reproduce the book's numbers. A fast gate runs with the unit tests.
- **Demo:** `python scenarios/T1.4/run.py`, a hydrogen-like pair under L0 and under L1a, with figures (exploratory).

### What it does not do yet

- The second field `B_cal` of moving bodies, the Boris push, and the readings L1b to L1e and options L4, L6, L7, L8.
- Populations, inflow, heat bath, absorbing walls; the fast solver for many bodies.
- The energy of the Composer coordinate (it needs a derivation approved by the physics lead; reported as `null`).
- Checkpoints and resuming; JSON schemas for the settings and the outputs; the stamp's lookup of the T0 status.
- T9.2, T9.4, and all of T1 to T8 except this demo.

### Status of the results

T0.1-T0.5, T9.1 and T9.3 are `confirmed`: reproduced independently by a second model (decision D14). D21.1-D21.4 are `done` and await reproduction. The demo is exploratory and is not a result. L1a is a hypothesis, and at lambda = 1 it is excluded by the measured charge of the neutron (T9.1).

### Conventions to know

Decision D21 (approved as refined after scenarios D21.1-D21.4): 2D is planar motion with a `1/r^2` field; softening has the Plummer form and is 0 in checks against closed forms; every result reports both the energy change at the end of a run and its largest excursion; the adaptive step has a floor for bodies at rest (0.1) and a displacement limit. The scope of the field law is slow bodies; it does not replace general relativity (`CONTRIBUTING.md`, rule 14).
