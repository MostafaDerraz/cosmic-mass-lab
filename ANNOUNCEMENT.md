# Call for contributors: build a virtual universe of cosmic masses

**An open, collective project for developers, physicists, mathematicians and AI agents.**

## The idea

Newton's law of gravitation and Coulomb's law have the same shape, yet physics describes every body by two numbers, its mass and its charge. The *cosmic mass* describes a body by one number,

    mu = m + q G,

where `G` is the second unit of a *General Number* plane: a number plane whose second axis leans at an angle and has a length of its own, `s = sqrt(K/G_N) = 1.16043e10 kg/C`. Products of cosmic masses follow one rule, the Rule of the Square, `G^2 = 2 s cos(phi) G - s^2`.

The **real part** of the product of two cosmic masses is exactly Newton plus Coulomb. A field law written with these numbers reproduces Maxwell's equations, the Lorentz force and Newtonian gravity with its magnetic part. All of this is verified: 106 symbolic and numerical checks for the physics, and over 10,000 for the mathematics.

The product also has a **second coordinate**, which mixes mass and charge: `m_a q_b + q_a m_b + 2 s cos(phi) q_a q_b`. No measurement has ever read it. Read at full strength as a force along the line, it is already ruled out: a neutron would respond as if it carried 9.0e-19 e, and the measured limit is about 1e-21 e. Read as anything else, it is an open question:

- a motion in a second coordinate;
- a torque;
- a change of inertia;
- a property of pairs;
- a weak force acting over astronomical time.

**We are building a simulator, a virtual universe made only of cosmic masses, to find out what each reading would do.**

## What we offer

- A verified foundation:
  - a formulas compendium;
  - a report of what is correct, corrected or excluded;
  - verification scripts that run automatically on every change.
- A complete, open test plan: about 50 scenarios (`docs/OPEN-SCENARIOS.md`), each with a setup, measurements and pass criteria, from the consistency gate T0 to large plasma scenes.
- An architecture for the engine (`docs/06-engine-architecture.mmd`).
- Honest framing. This is exploration, not a claim that physics is wrong. Every result, positive or negative, is recorded as a finding.

## What we are looking for

- **Developers:** the N-body engine, the field solver, integrators (Verlet, Boris), and the interactive web universe.
- **Physicists and mathematicians:** derivations of each reading, analytic predictions to check simulations against, and links to measurements.
- **AI agents:** running scenarios, writing reproducible scripts, checking each other's results, under the rules for agents in `AGENTS.md`.
- **Reviewers:** people who try to break results. A finding that survives an attempt to break it is worth more.

## How to join

1. Read `CONTRIBUTING.md` (the rules) and, for agents, `AGENTS.md`.
2. Pick an open scenario from `docs/OPEN-SCENARIOS.md`, and claim it with the "Claim a scenario" issue template.
3. Submit your work as a pull request with a reproducible script, the settings stamp and the output.
4. Your contribution is credited in `CONTRIBUTORS.md` and in any publication of the findings.

**First milestone:** an engine that passes the consistency gate T0. The engine must reproduce Newton plus Coulomb exactly, conserve mass, charge, energy and momentum, and give results that do not depend on the frame angle.

Contact: the project maintainer, through the repository's issues.
