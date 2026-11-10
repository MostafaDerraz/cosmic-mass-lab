# Virtual universe: scenarios and tests with detailed setups

Everything is written in the language of the cosmic mass and General Numbers (vocabulary in section 0). The verified formulas are in `04-formulas-compendium/`, and their status is in `03-report-verified-formulas.md`.

---

## 0. Vocabulary and constants

| Symbol | Value or meaning |
|---|---|
| `G_N` | 6.67430e-11 m^3 kg^-1 s^-2 |
| `K` | 8.9875517923e9 N m^2 C^-2 |
| `c` | 299792458 m/s |
| `s = sqrt(K/G_N)` | 1.16043e10 kg/C, the length of the Composer Unit `G` |
| `phi_bar` | frame angle of `G` (default 79.86 deg; also 90 deg, 10.14 deg, sweeps) |
| `mu = m + q G` | cosmic mass |
| `G^2 = 2 s cos(phi_bar) G - s^2` | Rule of the Square |
| `E_cal = g + e G` | cosmic field (Maxwell `E = -s^2 e`) |
| `B_cal = B_g + b G` | second field (Maxwell `B = -s^2 b`) |
| `rho_mu`, `J_mu` | density and current of cosmic mass |
| readings | `Re`, `Com` (Composer coordinate), `Sh`, `Ht`, `abs` |
| `m_S = e s` | 1.85921e-9 kg (Stoney mass) |

**Full product of a cosmic mass with the cosmic field** (both coordinates):

`mu E_cal = [m g + q E] + [q g - m E/s^2 - (2 cos(phi_bar)/s) q E] G`

### Matter presets (all as `mu = m + q G`)

| Name | m (kg) | q (C) |
|---|---|---|
| electron | 9.1093837e-31 | -1.602176634e-19 |
| proton | 1.67262192e-27 | +1.602176634e-19 |
| neutron | 1.67492750e-27 | 0 |
| alpha | 6.6446573e-27 | +3.204353268e-19 |
| dust grain | 1e-15 | +/- 1e-16 (about 600 e) |
| moon | 7.342e22 | 0 |
| Earth | 5.972e24 | 0 (or measured ~ -5e5) |
| Sun | 1.989e30 | 0 (or estimated ~ +100) |
| Saturn | 5.683e26 | 0 |

### Unit systems

| System | Definition | What changes |
|---|---|---|
| **SI** | values above | nothing |
| **Scaled** | `G_N = K = 1`, so `s = 1` | the drawing length becomes `A = 1`, and the natural frame angle is 90 deg. `phi_bar` stays a free setting for the readings L1b, L1c and L5c |

---

## 1. Settings

### U: universe

| Code | Setting | Choices |
|---|---|---|
| U1 | frame angle `phi_bar` | 79.86 deg / 90 deg / 10.14 deg / sweep 10..170 deg in 5 deg steps |
| U2 | dimension | 2D / 3D / 2D slice of 3D |
| U3 | units | SI / scaled / both side by side |
| U4 | boundary | open / reflecting walls / absorbing walls / wrap-around |
| U5 | background fields | none / uniform `g` / uniform `e` (i.e. `E`) / uniform `B_g` / uniform `b` (i.e. `B`) |
| U6 | fixed sources | list of pinned `mu` |
| U7 | inflow | rate, `mu` type and velocity of a stream entering from one side |
| U8 | population state | number density per `mu` type, velocity spread (temperature) per type, bulk flow; pressure is computed, not entered |

### L: laws

| Code | Law |
|---|---|
| L0 | `Re(mu) dv/dt = Re[mu(E_cal + v x B_cal)]` (paper; reference) |
| L1a | force = `Re(mu E_cal) + lambda s Com(mu E_cal)` (two forces along the field) |
| L1b | force = `Sh(mu E_cal) = Re + s cos(phi_bar) Com` |
| L1c | force = `Ht(mu E_cal) = s sin(phi_bar) Com` |
| L1d | force = `abs(mu E_cal)` along the field |
| L1e | `mu a = mu (E_cal + v x B_cal)` as a General Number equation, so `a = E_cal + v x B_cal` (motion in two coordinates) |
| L2 | sign of the Composer part: + or - |
| L3 | strength lambda: 1 / sweep 1e-24..1 (log, 25 steps) / cos(phi_bar) |
| L4 | direction of the Composer force: radial / perpendicular (+90 deg in the plane) / torque (spin) / along the velocity |
| L5 | terms: `s q g` only / `-m E/s` only / `-2 cos(phi_bar) q E` only / all |
| L6 | moving bodies: off / `Com` of `mu (v x B_cal)` on / that term only |
| L7 | inertia: `Re(mu) = m` / `Sh(mu) = m + q s cos(phi_bar)` / `abs(mu)` |
| L8 | range: 1/r^2 / screened by the surrounding population (length from the population state) |

**Default exploratory law:** L1a, sign +, lambda = 1, radial, all terms, L6 off, inertia `Re(mu)`, 1/r^2.

**Numerics (all runs):**

- integrator: velocity Verlet, or Boris for `v x B`;
- time step: adaptive, with `dt` so that the largest velocity change per step is below 1e-3 of the speed;
- softening length: 1e-3 of the smallest initial separation.

---

## 2. Tests

Each test lists its setup, what it measures and how it passes. "Pass" means agreement with the stated reference within the tolerance. "Finding" means there is no reference: we record what happens.

### T0: consistency gate (must pass before anything else)

| Test | Setup | Measured | Pass |
|---|---|---|---|
| T0.1 | L0, SI, 2D. (a) Earth and moon at 3.844e8 m, circular orbit speed. (b) electron and proton at the Bohr radius 5.29177e-11 m, circular speed. (c) two protons at 1e-10 m, at rest | period, radius, energy | (a) period 27.3 days within 0.1%; (b) period 1.5198e-16 s within 0.1%; (c) the motion matches Coulomb's analytic escape within 1e-6 |
| T0.2 | L0, 200 random `mu` (electrons, protons, neutrons), 2D box with reflecting walls, 1e4 steps | total `rho_mu` (both coordinates), energy `u`, momentum | mass and charge exactly constant; energy and momentum drift < 1e-6 |
| T0.3 | T0.1 repeated at `phi_bar` = 10, 45, 79.86, 90, 135, 170 deg | every `Re` result | identical to machine precision |
| T0.4 | L0, neutral planet (Earth); test bodies: electron, proton, neutron, dust, each at 7e6 m from rest | fall acceleration | all equal `G_N M/r^2` within 1e-12 |
| T0.5 | field solver: static `rho_mu` of one proton | `g` and `e` on a grid | Composer part gives `E = -s^2 e = K e/r^2` within 1e-6 |

### T9: reproduce the measured limits of the book (Chapter 10) before exploring

| Test | Setup | Measured | Pass |
|---|---|---|---|
| T9.1 | L1a radial, lambda = 1; a neutron at 1e-10 m from a fixed proton | the neutron's response, expressed as an effective charge `q_eff = -m_n/s` | `q_eff = 9.0e-19 e` (book); the bound "share below 0.002" against the measured (-0.2 +/- 0.8)e-21 e |
| T9.2 | charge-share reading only (L5c radial); hydrogen pair | binding energy | 5.71 eV (book), against 13.598 eV measured |
| T9.3 | L5a (`s q g` only); electron in Earth's field | equivalent electric field | 1.14e11 V/m (book) |
| T9.4 | L1a, lambda = 1; neutral planet (Earth); a proton, then an antiproton (`m_p - e G`), at 7e6 m | sign of the extra push from `s q g` | opposite pushes, as the book states ("particle against antiparticle") |

### T1: pairs of cosmic masses (run under L0 and every L1 to L8 variant)

| Test | Pair (`mu`) | Initial state | Question |
|---|---|---|---|
| T1.1 | Earth with moon (both `q = 0`) | circular orbit | Does the Composer part stay silent for pure masses? |
| T1.2 | Earth (`q = 0`) with proton; Earth with electron | at 7e6 m, at rest and in orbit | `s q g`: pull, push or swirl? Does the orbit close? |
| T1.3 | two protons; proton with electron | at 1e-10 m | How does the charge share `-2 cos(phi_bar) q E` change the electric force? |
| T1.4 | electron with proton (hydrogen-like) | Bohr orbit | Bound? Precession? Spiral? Binding energy? |
| T1.5 | two `mu` with `Re(mu1 mu2) = 0` (for example `m1 m2 = s^2 q1 q2`, the diagonal) | at 1 m, at rest | No interaction under L0. What under L1 to L8? |
| T1.6 | any pair above | | Newton's third law: force on 1 = -(force on 2)? |

### T2: survival of each coupling definition

- **Grid:** L1 (a-e) x L2 (+/-) x L4 (4 directions) x L7 (3 inertias) = 120 definitions, on T0.2 and T1.1-T1.6.
- **Classify each definition:**
  - (A) conserving and angle-free;
  - (B) conserving, angle-dependent;
  - (C) non-conserving.
- **Also record:** whether it passes T9 at lambda = 1, and the largest lambda that respects T9.1.
- **Output:** a survival table. Only classes A and B go on to T3-T7.

### T3: the angle

| Test | Setup | Question |
|---|---|---|
| T3.1 | T1.2-T1.4 with the `phi_bar` sweep | Which quantities change with the angle? Plot each against `phi_bar`. |
| T3.2 | `phi_bar` = 90 deg | Does the charge share vanish, giving the cleanest universe? |
| T3.3 | mirror pair 10.14 and 79.86 deg (same `s`) | Same, opposite or different worlds? |
| T3.4 | `phi_bar` -> 1 deg and 179 deg | How does the universe behave as the frame closes? |

### T4: many cosmic masses (SI and scaled)

| Test | Setup | Measured | Question |
|---|---|---|---|
| T4.1 | Sun (`q = 0`) fixed; 400 particles: 200 protons, 200 electrons, at 1 AU in a shell, velocity spread at 1e5 K | radial distribution by sign of `Com(mu)`, over 10 orbital times (scaled) | Does the cloud separate by charge sign? |
| T4.2 | Saturn (`q = 0`); 500 dust grains, `m = 1e-15` kg, `q = +/- 600 e` (equal numbers), ring at 1.2e8 m, Keplerian | ring radius, thickness and tilt by sign | Does the ring split, tilt or spread? |
| T4.3 | two planets (Earth `mu`, at 1e9 m); a population of 2000 protons and electrons (super-particles) around each; absorbing planetary surfaces (charge is collected into the planet's `q`) | each planet's `mu(t)`; `Re(mu1 mu2)`; their mutual force | Do the planets gain a Composer part, and does their attraction fall toward zero (`Q -> lambda M / s`)? |
| T4.4 | Sun; inflow of protons and electrons at 400 km/s, 5 per cm^3 (super-particles) | charge layer, separation | Does a separated `Com` layer form? |
| T4.5 | box, 1000 neutral pairs (bound electron + proton) at 300 K | stability, pressure | Do the pairs stay bound under each L? |

### T5: moving bodies (`B_cal`)

| Test | Setup | Question |
|---|---|---|
| T5.1 | proton crossing uniform `b` (`B` = 5 nT) and uniform `B_g` | Circling in the Real reading against the Composer reading (L6 on). |
| T5.2 | two parallel streams of protons, 1e5 m/s, 1 m apart | Real attraction or repulsion; what the `Com` part of `mu (v x B_cal)` adds. |
| T5.3 | rotating planet `mu` (Earth spin) with charged test `mu` nearby | Does rotation make a `B_cal` that acts through the Composer part? |

### T6: populations (temperature, pressure, screening)

| Test | Setup | Measured | Pass or question |
|---|---|---|---|
| T6.1 | L0; box, 1000 protons + 1000 electrons, 1e4 K | pressure from `Re` momentum flux on the walls | Pass: `p = n k T` within 2% |
| T6.2 | T6.1 under each surviving L | pressure, temperature, drift | Does the Composer part change them? |
| T6.3 | test `mu` (proton) inside the T6.1 population | its field and its coupling against distance | Is `E_cal` screened (Debye length)? Is the Composer coupling screened too? |
| T6.4 | T4.1 heated: 1e4 K, then 1e6 K | separation against temperature | Does heat destroy the separation? |

### T7: robustness

| Test | Check | Pass |
|---|---|---|
| T7.1 | T4.2 and T4.3 in SI and scaled units | same qualitative finding |
| T7.2 | T4.1 and T4.2 in 2D and 3D | same qualitative finding |
| T7.3 | particle counts 50 / 500 / 5000 | results converge |
| T7.4 | time step halved twice | results converge (change < 1%) |

### T8: limits of the language

| Test | Matter | Question |
|---|---|---|
| T8.1 | pure charge `mu = q G` (`m = 0`) | With inertia `Re(mu) = 0`, the motion is undefined. With L7 `Sh` or `abs`, what is its motion? |
| T8.2 | electron (large `Com/Re`) against planet (tiny `Com/Re`) in one run | Do the extremes stay numerically stable? |
| T8.3 | `mu` on the diagonal `G_N m^2 = K q^2` (two equal bodies) | Balance: no Real force between them. What under L1 to L8? |

---

## 3. Recorded in every run

- `Re` and `Com` of the total force on each body, over time.
- `Re` and `Com` of the total energy and momentum (`u`, `S`, stress), and their drift.
- `mu(t)` of each body: its readings `Re`, `Com`, `Sh`, `Ht`, `abs`, and its angle in the plane.
- Separation of the populations by the sign of `Com(mu)`.
- Dependence on `phi_bar` (for sweeps).
- Pass or fail of T0 and T9, with the settings stamped on the output.

## 4. Order of work

1. Engine and T0.
2. T9.
3. T1.
4. T2.
5. T3.
6. T4.
7. T5.
8. T6.
9. T7 and T8.
10. Interactive universe.

Each stage ends with a findings note. A stage starts only if T0 still passes.
