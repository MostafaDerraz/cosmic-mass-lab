# Virtual Universe of Cosmic Masses: final test and scenario list

Everything is stated in the language of the cosmic mass and General Numbers only.

## Vocabulary

| Symbol | Meaning |
|---|---|
| `G` | Composer Unit, length `s = sqrt(K/G_N) = 1.16043e10 kg/C`, frame angle `phi_bar` |
| `mu = m + q G` | cosmic mass of a body |
| `G^2 = 2 s cos(phi_bar) G - s^2` | Rule of the Square |
| `E_cal = g + e G` | cosmic field (Maxwell field `E = -s^2 e`) |
| `B_cal = B_g + b G` | second field (`B = -s^2 b`) |
| `rho_mu = rho_m + rho_q G`, `J_mu = J_m + J_q G` | density and current of cosmic mass |
| `Re`, `Com`, `Sh`, `Ht`, `abs` | readings of a General Number: real part, Composer coordinate, shadow, height, magnitude |
| F1-F4, M | Cosmic Field Law and law of motion (Cosmic Mass paper, Postulate 6.2) |

Full product: `mu E_cal = [m g + q E] + [q g - m E/s^2 - (2 cos(phi_bar)/s) q E] G`

---

## Part U: universe settings

- **U1, frame angle `phi_bar`:** 79.86 deg, 90 deg, the mirror root 10.14 deg, and a sweep from 10 to 170 deg.
- **U2, dimension:** 2D, 3D, or a 2D slice of a 3D universe.
- **U3, scale:** SI, scaled (all effects visible), or both side by side.
- **U4, boundary:** open, closed walls, wrap-around.
- **U5, background fields:** uniform `E_cal` (g and/or e), uniform `B_cal` (B_g and/or b), or none.
- **U6, fixed sources:** cosmic masses `mu` pinned in place, such as a star or a charged plate.
- **U7, inflow:** a stream of cosmic masses entering from one side (a wind of `mu`).
- **U8, state of a population, all derived from the cosmic masses:**
  - density `rho_mu` and current `J_mu`;
  - temperature = spread of the velocities of a population of `mu`;
  - pressure = momentum of the cosmic masses delivered to the walls, read by `Re`.

## Part M: matter (cosmic masses only)

- **M1, types:**
  - pure mass `mu = m`;
  - pure charge `mu = qG`, idealised;
  - electron-like `mu = m_e - e G`;
  - proton-like `mu = m_p + e G`;
  - neutral pair;
  - dust `mu`;
  - planet `mu`;
  - star `mu`.
- **M2, populations:** N copies of a `mu`, with positions, velocities and a velocity spread.
- **M3, readout per body:**
  - `mu` itself;
  - `Re(mu)`, `Com(mu)`, `Sh(mu)`, `Ht(mu)`, `abs(mu)`;
  - the angle of `mu` in the plane;
  - how `mu` changes over time (when a body gathers charge).

## Part L: laws (what the coupling is)

**L0, paper law (reference):**
`Re(mu) dv/dt = Re[mu (E_cal + v x B_cal)]`

**L1, reading of the product used as the force:**
- L1a: `Re + s Com` (two forces; law A)
- L1b: `Sh(mu E_cal)`
- L1c: `Ht(mu E_cal)`
- L1d: `abs(mu E_cal)`, along the field
- L1e: full General Number motion, `mu a = mu (E_cal + v x B_cal)` (law B; acceleration with two coordinates)

**L2, sign of the Composer part:** `+`, as the algebra gives, or `-`.

**L3, strength:** lambda = 1, a lambda sweep from 1e-24 to 1, or lambda = cos(phi_bar).

**L4, direction of the Composer force:**
- radial;
- perpendicular to the line between the bodies;
- torque (spin);
- along the velocity.

**L5, terms of the Composer part:**
- `s q g` alone;
- `-m E / s` alone;
- `-2 cos(phi_bar) q E` alone;
- all three.

**L6, moving bodies:**
- no `B_cal`;
- `Com` of `mu (v x B_cal)` included;
- that term alone.

**L7, inertia:** `Re(mu) = m`, `Sh(mu) = m + q s cos(phi_bar)`, or `abs(mu)`.

**L8, range:** `1/r^2` (from F1-F4), or screened by the surrounding population of `mu`.

---

## Part T: tests

### T0: consistency of the virtual universe (must pass first)

| Test | Check | Pass when |
|---|---|---|
| T0.1 | L0 with two bodies at rest | force = `-G_N Re(mu1 mu2)/r^2`, i.e. Newton and Coulomb exactly |
| T0.2 | conservation of `rho_mu` (both coordinates) | total mass and total charge constant |
| T0.3 | balance of energy and momentum (paper's `u`, `S`, stress `T`) | drift below the integrator tolerance |
| T0.4 | angle invariance under L0 | every `Re` result unchanged for every `phi_bar` |
| T0.5 | free fall under L0 in the field of a pure-mass source | every `mu` falls the same way |
| T0.6 | field split (F1-F4 in two coordinates) | the `Com` set reproduces Maxwell's field `E = -s^2 e` |

### T1: pairs of cosmic masses (every L-definition)

| Test | Pair | Question |
|---|---|---|
| T1.1 | pure mass with pure mass | does the Composer part stay silent? |
| T1.2 | pure mass with charged `mu` | the coupling `s q g`: pull, push, or swirl? |
| T1.3 | charged with charged | the angle term `-2 cos(phi_bar) q E`: how does it change the electric force? |
| T1.4 | electron-like with proton-like | does the bound pair stay bound? Does it precess or spiral? |
| T1.5 | two equal `mu` with `Re(mu1 mu2) = 0` | no interaction under L0; what under L1-L8? |
| T1.6 | Newton's third law | is the force on 1 equal and opposite to the force on 2? |

### T2: survival of each coupling definition

Run T0 and T1 for every combination L1 x L2 x L4 x L7. Classify each definition:

- **conserving and angle-free:** compatible with the paper's structure;
- **conserving but angle-dependent:** the frame becomes measurable;
- **non-conserving:** energy or momentum created or destroyed (a new kind of universe, or an inconsistency).

The definitions that survive go on to T3-T7.

### T3: the angle

| Test | Setting | Question |
|---|---|---|
| T3.1 | sweep `phi_bar` over T1.2-T1.4 | which effects depend on the angle? |
| T3.2 | `phi_bar` = 90 deg | does the angle term vanish, giving the cleanest universe? |
| T3.3 | the mirror pair 10.14 and 79.86 deg (same `s`) | identical, opposite, or different worlds? |
| T3.4 | `phi_bar` near 0 or 180 deg | how does the universe behave as the frame closes? |

### T4: many cosmic masses

| Test | Scene | Question |
|---|---|---|
| T4.1 | cloud of mixed `mu` (both signs of `Com`) around a pure-mass star | does the cloud separate by the sign of `Com(mu)`? |
| T4.2 | ring of charged dust `mu` around a planet `mu` | does the ring split, tilt, or spread? |
| T4.3 | two planets in a population of free charged `mu` | does each planet's `mu` gain a Composer part, and does `Re(mu1 mu2)`, and with it their attraction, fall toward zero? |
| T4.4 | star with inflowing wind of `mu` | does a layer of separated `Com` form? |
| T4.5 | gas of neutral pairs in a box | is it stable? |

### T5: moving bodies (`B_cal`)

| Test | Scene | Question |
|---|---|---|
| T5.1 | `mu` crossing a uniform `B_cal` (both coordinates) | how does it circle in the Real and in the Composer reading? |
| T5.2 | two parallel streams of `mu` | Real part attracts or repels; what does the `Com` part of `mu (v x B_cal)` add? |
| T5.3 | rotating planet `mu` | does its motion make a `B_cal` that acts on nearby charged `mu`? |

### T6: populations (temperature, pressure, screening)

| Test | Scene | Question |
|---|---|---|
| T6.1 | hot population of `mu` in a box (L0) | pressure from `Re` momentum flux behaves like `n k T` |
| T6.2 | the same under each surviving L | does the Composer part change pressure or temperature? |
| T6.3 | test `mu` inside a population | is its `E_cal` screened, and is its Composer coupling screened too? |
| T6.4 | heating and cooling of a population | does the separation of T4.1 survive heat? |

### T7: robustness of the virtual universe

| Test | Check | Question |
|---|---|---|
| T7.1 | SI against scaled units on T4.2 and T4.3 | same qualitative answer? |
| T7.2 | 2D against 3D on T4.1 and T4.2 | same answer? |
| T7.3 | particle count (50 / 500 / thousands) | do results converge? |
| T7.4 | time step halved | do results converge? |

### T8: limits of the language

| Test | Matter | Question |
|---|---|---|
| T8.1 | pure charge `mu = qG` (no `Re`) | what is its motion when inertia is `Re(mu)` = 0? (Under L7 with `Sh` or `abs`, does it move?) |
| T8.2 | `mu` with very large `Com/Re` (electron-like) against very small (planet) | do extremes stay stable? |
| T8.3 | `mu` whose angle equals the frame angle (`m = 0`) | special behaviour? |

---

## Measured in every run (cosmic-mass language)

- `Re` and `Com` of the total force on each body.
- `Re` and `Com` of the total energy and momentum.
- `mu` of each body over time, including its angle in the plane.
- separation of populations by the sign of `Com(mu)`.
- dependence of every result on `phi_bar`.
- pass or fail of the T0 checks.

## Order of work

1. Engine and T0.
2. T1 and T2 (which couplings survive).
3. T3.
4. T4.
5. T5.
6. T6.
7. T7 and T8.
8. Interactive virtual universe with all U, M, L settings and a compare mode.

Each step ends with a short findings note.
