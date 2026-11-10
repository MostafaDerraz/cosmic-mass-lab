# Report: which formulas and laws are correct

Work preparation package, 10 November 2026. This report covers the four sources in `01-sources/`:

- the book *The Cosmic Mass, Volume 1* (clean sources of 9 November 2026);
- the Cosmic Mass paper (the version posted on Zenodo, doi:10.5281/zenodo.24094855);
- the General Numbers paper (Zenodo version, doi:10.5281/zenodo.24091305, and the Gazette version);
- the General Sinus paper (revision 3, the version submitted to *Elemente der Mathematik*).

The status of every item rests on one of three things:

- **a verification script** in `02-verification-scripts/`, re-run for this package on the copied sources (results below);
- **an audit done in this conversation** (numerical, at 40-50 significant digits);
- **a decision recorded in this conversation**, such as excluding a chapter's claims from the papers.

Nothing in this report is a judgement about whether the physics describes nature. It says whether the mathematics is right, and where the sources themselves say a law fails a measurement.

---

## 1. Verification runs for this package

| Script | Checks | Result | Covers |
|---|---|---|---|
| `cosmic-mass/verify_cmt.py` | 106 (41 symbolic, 9 numerical, 55 printed values, 1 coverage) | **106 / 106 pass** | every numbered statement and computed number of the Cosmic Mass paper |
| `general-numbers/verify_paper.py` | 204 | **204 / 204 pass** | every numbered statement of the General Numbers paper (51 statements, 62,387 numerical comparisons in 75 planes) |
| `general-numbers/book-algebra/verify_general_number_algebra.py` | 253 | **253 / 253 pass** | the algebra of the book's Chapter 4 (180 symbolic checks for every angle and length) |
| `general-numbers/book-laws-test/test_general_number_laws.py` | 144 laws, 684,167 evaluations | **all laws hold** (largest error 9.5e-13) | the laws of the book's Chapter 4 and Appendix F, with 14 wrong laws and 5 planted errors caught as negative controls |
| `general-numbers/gazette-version/check_gazette.py` | 9,304 | **all pass** | every identity and printed number of the Gazette article |
| `general-sinus/verify_gsp.py` | 79 | **79 / 79 pass** | every statement of General Sinus revision 3, including the new closed form and rate |

Each script's output is saved next to it as `verify_output.txt`.

---

## 2. General Sinus (oblique triangle)

### 2.1 Correct (verified)

Correct, and checked by `verify_gsp.py` at 50 digits:

- **Definitions:** `sin(a,p) = sin a / sin(a+p)`, `cos(a,p) = sin p / sin(a+p)`, `tan(a,p) = sin a / sin p`, and the reciprocals.
- **Coordinate reading:** the general cosine and sine are the coordinates of the unit circle in the oblique frame `(e1, u)`, with `u` at the angle `a+p`.
- **Oblique Pythagorean identity:** `cos^2 + sin^2 + 2 cos sin cos(a+p) = 1`.
- **Exchange, reciprocals, signs and pi-shifts, special values** (`p = pi/2`, `a = p`, `pi/3`, the golden ratio at `pi/5`).
- **Addition in the studied angle** (all four forms in the book's Appendix A, and the closed forms for `cos(a+b,p)` and `sec(a+b,p)`), **addition in the projection angle**, and the double- and half-angle formulas listed in revision 3 and in the book's Appendix A. All reduce to Ptolemy's identity.
- **Addition in a fixed frame** (new in revision 3): `C(a+b) = C(a)C(b) - S(a)S(b)`, `S(a+b) = C(a)S(b) + S(a)C(b) + 2 cos(pb) S(a)S(b)`, and De Moivre's formula.
- **First partial derivatives**, and the six relations among them in the book's Appendix A.
- **Second partial derivatives, corrected forms:** see 2.2.
- **Quarter-turn shift:** `D^n f(x,y) = f(x + n pi/2, y - n pi/2)` for `f = sin, cos`, where `D = d/dx - d/dy` (frame held fixed). In particular `D^2 f = -f`. Primitives follow by the inverse shift.
- **Product identities for every N >= 1**, not only `N = 2^n`.
- **General Sinus Transform:** exact closed form `S_N = (rho(x/N)/rho((x+y)/N))^N sin(N delta(x/N)) / sin(N delta((x+y)/N))`, and the rate `S_N = sin(x,y) + c(x,y)/N^2 + O(N^-3)`.

### 2.2 Found false and corrected

These are in the Research Square preprint (rs-9525842 v1) and the rev1 file, **not** in the book, whose Appendix A holds them back with a placeholder.

| Original statement | Status | Correct form |
|---|---|---|
| `tan(2a,p) - sec(2a,p) = sin(a,p)` | **false** | deleted (no correct form of this shape) |
| `tan(a,p) = sin(a/2,p) + sec(a,p)` | **false** | deleted |
| `d2/dy2 sin(x,y) = sin(x,y)(2cot^2(x+y) - 1)` | **false** | `sin(x,y)(1 + 2cot^2(x+y))` |
| `d2/dx2 cos(x,y) = cos(x,y)(2cot^2(x+y) - 1)` | **false** | `cos(x,y)(1 + 2cot^2(x+y))` |
| `d2/dy2 tan(x,y) = tan(x,y)(2cot^2(x+y) - 1)` | **false** | `tan(x,y)(1 + 2cot^2 y)` |
| total derivative `sin(x,y)'' = +sin(x,y)`, `cos'' = +cos` ("no sign change") | **false** | `D^2 sin = -sin`, `D^2 cos = -cos` |
| n-th derivative and primitive by the operator `omega(z,n)` | **false for n = 2, 3 (mod 4)** | `D^n f(x,y) = f(x + n pi/2, y - n pi/2)` |
| `tan^(n) = n! tan' tan(pi/2+y, y)^(n-1)` | **true only for n = 1, 2** | `D^n tan(x,y) = sin(x+y) (-d/dy)^(n-1) csc^2 y` |
| recurrence `u_{n,k}`, `v_{n,k}` "exactly" `sin`, `cos` | **false** (approximation only; error accumulates) | use exact `sin(pi k/N)`, `cos(pi k/N)` |
| "five-digit accuracy for n >= 20 at every point" | **false with the recurrence** (error 3e-4 at N = 2^18 near `x+y -> pi`) | rate `c(x,y)/N^2` with explicit `c`; proved |
| heat-equation analogy (fit `R^2 = 0.83`) | **withdrawn** (not evidence) | removed from the paper |

**For the book:** the placeholder at the end of Appendix A ("Still to come") can now be filled with the corrected forms above.

---

## 3. General Numbers (relative planes)

Correct and verified: 204 checks on the paper, plus 253 checks and 144 laws on the book's Chapter 4.

- **Axiom:** `|G| = |cos pb + sin pb|`, `G = |G| u`, with `u = (cos pb, sin pb)`.
- **Two angles for one length:** `pb1 = (1/2) arcsin(A^2 - 1)` and `pb2 = pi/2 - pb1`. For `A = 1.16043` they are 10.14 and 79.86 degrees.
- **Rule of the Square:** `G^2 = 2|G| cos pb G - |G|^2`.
- **Product:** `zw = (ac - |G|^2 bd) + (ad + bc + 2|G| cos pb bd) G`. Its real part depends on the angle only through `|G|`.
- **Conjugate, skew flip, division, the norm `a^2 + 2ab|G|cos pb + b^2|G|^2`** (multiplicative).
- **Turning Unit:** `J = -cot pb + G/(|G| sin pb)`, with `J^2 = -1`. The isomorphism `z -> Sh z + i Ht z` onto C.
- **Laws of the parts:** `Com(zw) = Sh(z)Com(w) + Sh(w)Com(z)`, `Re(conj z) = 2 Sh z - Re z`, and the full habits table.
- **Euler's formula in a relative plane:** `e^{theta J} = C_theta + S_theta u`. Also De Moivre, the powers of G (Chebyshev coefficients), the logarithmic spiral, and the n-th roots.
- **Geometry of one magnitude** (ellipse with axes on the diagonals, ratio `tan(pb/2)`), the polar-like frame, and the Cauchy-Riemann equations in a leaning frame.

**Caveats** (stated in the paper itself, not errors):

- Every true relative plane is the complex plane C written in an oblique frame. The algebra is that of the elliptic generalized complex numbers. What is new is the frame, its trigonometry and the laws written in it.
- The Axiom's choice of length has a geometric meaning, but no principle singles it out (an open question).

---

## 4. Cosmic Mass paper (Zenodo version)

Correct and verified: 106 / 106.

- **Cosmic mass:** `mu = m + q G`, with `|G| = s = sqrt(K/G) = 1.16043e10 kg/C`, and the drawing length `A = s/u = 1.16043` for `u = 1e10 kg/C`.
- **Cosmic Interaction Law:** read by the real part, `F = G Re(mu_a mu_b)/r^2 = (G m_a m_b - K q_a q_b)/r^2`. This is Newton and Coulomb, with the minus sign supplied by the Rule of the Square.
- **Angle invariance:** every measured quantity depends on s only. The angle is a convention tied to the unit. In Planck or Stoney units, `A = 1` and `pb = 90` degrees.
- **Stoney mass:** `m_S = e s = 1.85921e-9 kg = sqrt(alpha) m_P`.
- **Cosmic Field Law (F1-F4) with the law of motion (M).** It splits into Newton plus Heaviside gravitomagnetism (Real set) and Maxwell (Composer set, with `E = -s^2 e` and `B = -s^2 b`). The rest of the field law:
  - the Lorentz force;
  - waves at c in both coordinates;
  - conservation of `rho_mu`;
  - the energy and momentum balances;
  - energy density `u = -(g^2 + c^2 B_g^2)/(8 pi G) + (E^2 + c^2 B^2)/(8 pi K)`;
  - the action.

**Physical limits** (stated in the paper's Table 5, and true):

- The gravitomagnetic frame dragging on LAGEOS comes out at **1/4 of the measured value**.
- Gravity does not act on light: **no light bending** and no Shapiro delay.
- **Gravity waves carry negative energy**, which contradicts the observed inspirals.
- Rest mass is conserved, so the law **cannot describe fusion** or other mass-to-energy processes.
- Valid for **slow bodies** only.

**Decision recorded in this conversation:** the paper reads only the real part, and says the Coupling term "exerts no force". The second coordinate is left unread, with the conditions a reading must meet.

---

## 5. The book, chapter by chapter

| Chapter | Content | Status |
|---|---|---|
| 1. One Number for Two Laws | introduction | prose; its formulas are those of the papers |
| 2. Old Laws, New Results (`Alignment.tex`) | old laws in the new language | consistent with the papers; not separately audited |
| 3. The Oblique Triangle, and Appendix A | General Sinus | **verified** (section 2 above); the placeholder can now be filled |
| 4. General Numbers, and Appendix F | the algebra | **verified** (253 checks, 144 laws) |
| 5. Prologue | prose | none needed |
| 6. The Interaction Law | the Cosmic Interaction Law and the real reading | **verified** through the paper's script for the parts the paper reproduces; the remaining remarks are not separately audited |
| 7. The Cosmic Mass Standard Model | identity cards and the frame | definitions are consistent mathematics. The physical claims are **not verified**, and were kept out of the paper by decision in this conversation |
| 8. The Closure Law | an error law refitted on measured masses | **empirical fit, not derived**; kept out of the paper by decision |
| 9. The Unified Cosmic Law | the field law | **verified** through the paper's script; the same physical limits as in section 4 |
| 10. What Is the Coupling Force? | readings of the Composer coordinate | **open by design.** Its measurement table (section 6 below) already **excludes every reading as a force along the line at full size** |
| 11. What the Frame Allows | antimatter, negative mass, black holes | "allowed, not measured"; black-hole facts are quoted from general relativity, not derived |
| Appendices B-E | tables and supporting material | not separately audited |

---

## 6. What measurement already says about the Coupling term (book, Chapter 10, Table)

| Reading | What it would do | Measurement | Status |
|---|---|---|---|
| Coupling term along the line | a neutron acts as if charged, 9.0e-19 e | neutron charge (-0.2 +/- 0.8)e-21 e | **excluded; any share below ~0.002** |
| charge share along the line | hydrogen holds its electron with 5.71 eV | 13.598 eV | **excluded** |
| a charge in gravity | every charge as if in 1.14e11 V/m | electrons in a copper tube (Witteborn-Fairbank) | hidden by metal; strength limited by the neutron |
| Coupling term inside the atom | 9.0e-19 of Coulomb's pull in hydrogen | atomic spectra | too small to see; no test |
| particle against antiparticle | opposite pushes near a neutral mass | antihydrogen falls at (0.75 +/- 0.13 +/- 0.16) g | neutral, so no test |
| products as bodies | charge not kept | every fusion | excluded as composition |
| the diagonal `G m^2 = K q^2` | line of balance, extremal black holes | none needed | holds as mathematics; no new prediction |

**Consequence for the simulator:** any coupling reading used "as a force along the line" must have a strength **below about 0.002 of its natural size**, or be something other than a force along the line. The simulator must reproduce these numbers before it explores anything else (test family T9 in `05-scenarios-and-tests.md`).

---

## 7. What the simulator may take as established

**Established inputs:**

- every formula in sections 2.1, 3 and 4;
- the paper's law L0: `Re(mu) dv/dt = Re[mu (E + v x B)]`, with the field law F1-F4.

**Exploratory hypotheses** (the subject of the experiments, never established inputs):

- every reading of the Composer coordinate as a force or a motion, laws L1-L8 of the test plan;
- any physical claim of the book's Chapters 7, 8, 10 and 11.
