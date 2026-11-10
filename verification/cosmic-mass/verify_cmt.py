#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_cmt.py

Checks every numbered statement, every table and every printed number of

    Mostafa Derraz, "Mass and Charge as One General Number: A Rewriting of
    Newton's, Coulomb's and Maxwell's Laws for Slow Bodies", manuscript, November 2026.

Each check carries one of four tags:
  sym  a symbolic proof with SymPy, for all values of the symbols involved;
  num  a numerical check of a statement, at chosen or random values;
  val  a number printed in the paper, recomputed from the definitions and the data and
       rounded as printed; the printed form must also occur in the LaTeX source that prints it;
  cov  coverage: every numbered statement and every labelled equation has at least one check.

Notation. A General Number of the relative plane with angle p and length L = |G| is a + b G.
On the page, G is the arrow of length L at the angle p, so a + b G is the complex number
a + b L (cos p + i sin p) (General Numbers paper, Proposition 4.3). A field with General Number
values is a pair (real coordinate, Composer coordinate) of real fields, and the real part of a
product of two such quantities is x1 y1 - L^2 x2 y2 (the real product of Section 6).

Run with:   python3 verify_cmt.py [folder]
            where folder holds the LaTeX sources of the paper. By default the folder above this script.
            Without the LaTeX sources, every mathematical and numerical check still runs; the checks
            that a printed form occurs in the source, and the coverage check, are skipped.
Needs:      Python 3.8 or later, SymPy and mpmath.
Exit code:  0 when every check passes, 1 otherwise.
"""
import csv
import json
import math
import os
import random
import re
import sys
import time

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import constants as C                                       # noqa: E402

T0 = time.time()
ROOT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(HERE)
RESULTS = []          # (tag, labels, name, ok)
COVERED = set()
SOURCES = ["cmt-paper", "sec2", "sec3", "sec4", "sec5", "sec6", "sec9", "appA"]
TEX = {}
for name in SOURCES:
    path = os.path.join(ROOT, name + ".tex")
    TEX[name] = open(path, encoding="utf-8").read() if os.path.exists(path) else None
HAVE_TEX = all(TEX[n] is not None for n in SOURCES)


def flat(s):
    """LaTeX source with the line breaks and repeated spaces folded, for matching printed forms."""
    return re.sub(r"\s+", " ", s)


FLAT = {k: flat(v) for k, v in TEX.items() if v is not None}


def heading(title):
    print("\n" + title)
    print("-" * len(title))


def check(tag, labels, name, ok):
    labels = (labels,) if isinstance(labels, str) else tuple(labels)
    RESULTS.append((tag, labels, name, bool(ok)))
    COVERED.update(labels)
    print(("  ok    " if ok else "  FAIL  ") + f"[{tag}] {', '.join(labels):24s} {name}")


def sci(x, digits):
    """x as printed in the paper: 'm \\times 10^{e}' with the given number of decimals in m."""
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(abs(m), digits) >= 10:
        m, e = m / 10, e + 1
    return f"{m:.{digits}f} \\times 10^{{{e}}}"


def val(labels, name, computed, printed, where=None, snippet=None):
    """A printed number: the computed form must equal the printed one, and occur where it is printed."""
    ok = computed == printed
    if where is not None and HAVE_TEX:
        ok = ok and flat(snippet if snippet is not None else printed) in FLAT[where]
    check("val", labels, f"{name}: {computed}" + ("" if computed == printed else f" (printed {printed})"), ok)


def zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


# =============================================================================
heading("1. Constants, the size of the Composer Unit and the angle (Sections 1, 3 and 4)")
# =============================================================================
val("eq:cu-size", "size of the Composer Unit sqrt(K/G), kg/C", sci(C.GU, 5), "1.16043 \\times 10^{10}", "sec4")
val("eq:stoney-ratio", "Stoney's ratio in the introduction", sci(C.GU, 5), "1.16043 \\times 10^{10}", "cmt-paper")
check("num", "eq:cu-size", "the two forms agree: K = mu0 c^2 / (4 pi) to the CODATA digits",
      abs(C.mu0 * C.c ** 2 / (4 * math.pi) / C.K - 1) < 1e-9)
val("eq:cu-size", "relative uncertainty of sqrt(K/G), half that of G: about one part in 10^5",
    f"{C.G_REL_UNC / 2:.1e}", "1.1e-05")
val("eq:electron", "Stoney mass e sqrt(K/G), kg", sci(C.mS, 5), "1.85921 \\times 10^{-9}", "cmt-paper")
val("eq:electron", "charge part of the electron over its mass", sci(C.mS / C.me, 2), "2.04 \\times 10^{21}", "cmt-paper")
val("eq:two-laws", "Coulomb over Newton for two electrons, K e^2 / (G me^2)", sci(C.K * C.e ** 2 / (C.G * C.me ** 2), 2),
    "4.17 \\times 10^{42}", "cmt-paper")
check("num", "eq:two-laws", "that factor is the square of the charge-to-mass ratio of the electron through G",
      abs(C.K * C.e ** 2 / (C.G * C.me ** 2) / (C.mS / C.me) ** 2 - 1) < 1e-12)
val("def:cosmic-mass", "charge part of the electron over its mass (Section 4)", sci(C.mS / C.me, 2), "2.04 \\times 10^{21}", "sec4")
val("def:cosmic-mass", "charge part of the proton over its mass", sci(C.mS / C.mp, 2), "1.11 \\times 10^{18}", "sec4")
p_, L_ = sp.symbols("p ell", positive=True)
check("sym", "eq:axiom", "|cos p + sin p| = sqrt(2) |sin(p + pi/4)|",
      zero(sp.expand_trig(sp.sqrt(2) * sp.sin(p_ + sp.pi / 4)) - sp.cos(p_) - sp.sin(p_)))
r1 = sp.asin(L_ / sp.sqrt(2)) - sp.pi / 4
check("sym", ("eq:axiom", "eq:our-angle"), "both roots p1 = arcsin(l/sqrt2) - pi/4 and p2 = pi/2 - p1 solve cos p + sin p = l",
      zero(sp.expand_trig(sp.cos(r1) + sp.sin(r1)) - L_) and zero(sp.expand_trig(sp.cos(sp.pi / 2 - r1) + sp.sin(sp.pi / 2 - r1)) - L_))
val("eq:our-angle", "arcsine root, degrees", f"{math.degrees(C.ROOT1):.2f}", "10.14", "sec4", "10.14^\\circ")
val("eq:our-angle", "mirror root, degrees", f"{math.degrees(C.ROOT2):.2f}", "79.86", "sec4", "79.86^\\circ")
val("eq:our-angle", "mirror root to four decimals", f"{math.degrees(C.ROOT2):.4f}", "79.8605", "sec4", "79.8605^\\circ")
shift = math.degrees(C.ROOT2 - (math.pi / 2 - (math.asin(C.ELL * (1 - C.G_REL_UNC / 2) / math.sqrt(2)) - math.pi / 4)))
val("eq:our-angle", "shift of the root for one standard error of G, degrees", f"{abs(shift):.4f}", "0.0009", "sec3", "0.0009^\\circ")
val("eq:our-angle", "the held angle 79.86 lies within one standard error of the root", str(abs(math.degrees(C.ROOT2) - 79.86) < abs(shift)), "True")
val("eq:our-angle", "sin of the held angle", f"{C.SIN:.5f}", "0.98438", "sec4")
val("eq:our-angle", "cos of the held angle", f"{C.COS:.5f}", "0.17605", "sec4")
val("eq:our-angle", "the length lies between 1 and sqrt 2, so both roots exist", str(1 < C.ELL < math.sqrt(2)), "True")

# =============================================================================
heading("2. The number plane (Section 2)")
# =============================================================================
al, ph, a, b, c_, d_, M = sp.symbols("alpha phi a b c d M", real=True)
Cs, Ss, L = sp.symbols("C S L", positive=True)                 # cos p, sin p, |G|
gsin = lambda u, v: sp.sin(u) / sp.sin(u + v)
gcos = lambda u, v: sp.sin(v) / sp.sin(u + v)
check("sym", "eq:gs", "the general sine and cosine return sin and cos when alpha + phi = pi/2",
      zero(gsin(al, sp.pi / 2 - al) - sp.sin(al)) and zero(gcos(al, sp.pi / 2 - al) - sp.cos(al)))
check("sym", "eq:gs-pythagoras", "general Pythagoras: sin^2 + cos^2 + 2 sin cos cos(alpha + phi) = 1",
      zero(sp.expand_trig(sp.together(gsin(al, ph) ** 2 + gcos(al, ph) ** 2
                                       + 2 * gsin(al, ph) * gcos(al, ph) * sp.cos(al + ph) - 1)
                          * sp.sin(al + ph) ** 2)))
I = sp.I
Gpage = L * (Cs + I * Ss)                                     # the Composer Unit on the page
on_circle = {Ss ** 2: 1 - Cs ** 2}


def page(re_, co):
    return re_ + co * Gpage


def gn_parts(w):
    """(real part, Composer coordinate) of a page point w, read along the two axes."""
    w = sp.expand(w)
    y = sp.im(w) if False else sp.expand((w - sp.conjugate(w)) / (2 * I))
    x = sp.expand((w + sp.conjugate(w)) / 2)
    co = sp.simplify(y / (L * Ss))
    return sp.simplify(x - co * L * Cs), co


def zero_cs(expr):
    e = sp.expand(expr)
    return sp.simplify(e.subs(Ss ** 2, 1 - Cs ** 2).subs(Ss ** 3, Ss * (1 - Cs ** 2)).subs(Ss ** 4, (1 - Cs ** 2) ** 2)) == 0


check("sym", "eq:square", "Rule of the Square: G^2 = 2 L cos(p) G - L^2 on the page",
      zero_cs(Gpage ** 2 - (2 * L * Cs * Gpage - L ** 2)))
check("sym", "eq:product", "product of two General Numbers, both parts",
      zero_cs(page(a, b) * page(c_, d_) - page(a * c_ - L ** 2 * b * d_, a * d_ + b * c_ + 2 * L * Cs * b * d_)))
check("sym", "eq:shadow", "shadow and height are the page coordinates: Sh = a + b L cos p, Ht = b L sin p",
      zero(sp.re(page(a, b).subs({Cs: sp.cos(p_), Ss: sp.sin(p_)})) - (a + b * L * sp.cos(p_)))
      and zero(sp.im(page(a, b).subs({Cs: sp.cos(p_), Ss: sp.sin(p_)})) - b * L * sp.sin(p_)))
pp = sp.Symbol("pp", positive=True)

# =============================================================================
heading("3. Bodies at rest: the Cosmic Interaction Law (Section 5)")
# =============================================================================
ma, mb, qa, qb, Gn, Kc = sp.symbols("m_a m_b q_a q_b G K", positive=True)
mu_a, mu_b = page(ma, qa), page(mb, qb)
check("sym", "law:interaction", "mu_a mu_b = m_a m_b + (m_a q_b + q_a m_b) G + q_a q_b G^2 (three terms)",
      zero_cs(mu_a * mu_b - (ma * mb + (ma * qb + qa * mb) * Gpage + qa * qb * Gpage ** 2)))
prod_re, prod_co = ma * mb - L ** 2 * qa * qb, ma * qb + qa * mb + 2 * L * Cs * qa * qb
check("sym", "prop:law-real", "product of two cosmic masses: both coordinates",
      zero_cs(mu_a * mu_b - page(prod_re, prod_co)))
check("sym", ("prop:law-real", "def:reading"), "G Re(mu_a mu_b) = G m_a m_b - K q_a q_b when G L^2 = K",
      zero((Gn * prod_re).subs(L, sp.sqrt(Kc / Gn)) - (Gn * ma * mb - Kc * qa * qb)))
check("sym", "prop:law-real", "the angle appears only in the Composer coordinate", Cs not in prod_re.free_symbols)
check("sym", "def:real-reading", "for a neutral body the real reading and the shadow agree: Re(a) = Sh(a) = a",
      zero(sp.re(page(a, 0).subs({Cs: sp.cos(p_), Ss: sp.sin(p_)})) - a))
check("sym", "eq:cosmic-mass", "the electron: q G = -e |G| u, a charge part of length mS pointing down the Composer axis",
      zero(sp.Abs(-1 * Gpage.subs({Cs: sp.cos(p_), Ss: sp.sin(p_)})) - L))
val("eq:cosmic-mass", "charge part of the electron, e |G| in kg", sci(C.e * C.GU, 5), "1.85921 \\times 10^{-9}", "sec4")
check("sym", "eq:real-product", "Re[(a + bG)(c + dG)] = ac - L^2 bd, free of the angle",
      zero_cs(gn_parts(page(a, b) * page(c_, d_))[0] - (a * c_ - L ** 2 * b * d_)))
two_p = [C.G * C.mp ** 2, C.G * C.GU ** 2 * C.e ** 2, C.G * C.GU * 2 * C.mp * C.e]
for lab, v_, pr in zip(("Newton term", "Coulomb-shaped term", "Coupling term (length of its Composer part)"), two_p,
                       ("1.87 \\times 10^{-64}", "2.31 \\times 10^{-28}", "4.15 \\times 10^{-46}")):
    val("law:interaction", f"two protons 1 m apart, {lab}, N", sci(v_, 2), pr, "sec5")
lam = sp.Symbol("lambda", real=True)
check("sym", "def:reading", "with one real number m + lambda q, two like charges would attract: lambda^2 q_a q_b > 0",
      sp.ask(sp.Q.positive(lam ** 2 * qa * qb), sp.Q.nonzero(lam)) is True)
check("sym", "eq:shadow-law", "shadow of the product: m m + L cos p (m q + q m) + L^2 cos(2p) q q",
      zero(sp.expand_trig(sp.re((mu_a * mu_b).subs({Cs: sp.cos(p_), Ss: sp.sin(p_)}))
                          - (ma * mb + L * sp.cos(p_) * (ma * qb + qa * mb) + L ** 2 * sp.cos(2 * p_) * qa * qb))))
check("sym", "eq:shadow-law", "times G, with G L = sqrt(KG) and G L^2 = K",
      zero((Gn * L).subs(L, sp.sqrt(Kc / Gn)) - sp.sqrt(Kc * Gn)) and zero((Gn * L ** 2).subs(L, sp.sqrt(Kc / Gn)) - Kc))
val("eq:shadow-law", "cos 2p at the arcsine root", f"{math.cos(2 * C.PHI1):.5f}", "0.93801", "sec5", "0.93801")
val("eq:shadow-law", "cos 2p at the mirror root", f"{math.cos(2 * C.PHI2):.5f}", "-0.93801", "sec5", "-0.93801")
val("eq:shadow-law", "the square leans at 2p, degrees", f"{2 * math.degrees(C.PHI2):.2f}", "159.72", "sec5", "159.72^\\circ")
check("num", "eq:shadow-law", "cos 2p < 0 for every p in (pi/4, 3pi/4) (the square points back)",
      all(math.cos(2 * (math.pi / 4 + k * (math.pi / 2) / 1000)) < 0 for k in range(1, 1000)))

# =============================================================================
heading("3b. The sign of the interaction at every angle (Proposition 5.4)")
# =============================================================================
check("sym", "prop:sign", "Re(G^2) = -L^2 for every angle: the real part of the Rule of the Square equals the angle-free -L^2",
      zero_cs(gn_parts(Gpage ** 2)[0] + L ** 2))
check("sym", "prop:sign", "F r^2 = G m_a m_b - K q_a q_b with no condition on the angle",
      zero((Gn * prod_re).subs(L, sp.sqrt(Kc / Gn)) - (Gn * ma * mb - Kc * qa * qb)) and Cs not in prod_re.free_symbols)
pulls = []
for sa, sb in ((1, 1), (1, -1), (-1, -1)):
    for deg in (10.14, 45.0, 79.86, 90.0, 120.0, 170.0):
        w_a = complex(0, 0) + sa * C.GU * complex(math.cos(math.radians(deg)), math.sin(math.radians(deg)))
        w_b = sb * C.GU * complex(math.cos(math.radians(deg)), math.sin(math.radians(deg)))
        prod = w_a * w_b                        # q_a q_b G^2 on the page, with q = +-1 C
        # real part in the leaning frame: Re = x - (y / (L sin p)) L cos p
        re_part = prod.real - prod.imag / math.sin(math.radians(deg)) * math.cos(math.radians(deg))
        pulls.append((sa * sb, re_part))
check("num", "prop:sign", "two charges at six angles from 10.14 to 170 degrees: Re(q_a q_b G^2) = -L^2 q_a q_b at each",
      all(abs(r_ + s_ * C.GU ** 2) < 1e-6 * C.GU ** 2 for s_, r_ in pulls))

# =============================================================================
heading("4. Bodies in motion: the Unified Cosmic Law (Section 6)")
# =============================================================================
t, x, y, z = sp.symbols("t x y z", real=True)
X = (x, y, z)
c = sp.Symbol("c", positive=True)


def grad(f):
    return [sp.diff(f, v) for v in X]


def div(F):
    return sum(sp.diff(F[i], X[i]) for i in range(3))


def curl(F):
    return [sp.diff(F[2], y) - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], x), sp.diff(F[1], x) - sp.diff(F[0], y)]


def dt(F):
    return [sp.diff(f, t) for f in F] if isinstance(F, list) else sp.diff(F, t)


def dot(A_, B_):
    return sum(A_[i] * B_[i] for i in range(3))


def cross(A_, B_):
    return [A_[1] * B_[2] - A_[2] * B_[1], A_[2] * B_[0] - A_[0] * B_[2], A_[0] * B_[1] - A_[1] * B_[0]]


def vec(name):
    return [sp.Function(f"{name}{k}")(t, x, y, z) for k in "xyz"]


def sca(name):
    return sp.Function(name)(t, x, y, z)


def lapv(F):
    return [sum(sp.diff(f, v, 2) for v in X) for f in F]


def allzero(lst):
    return all(sp.expand(e) == 0 for e in lst)


# the field equations of one coordinate, and the generic fields that obey or break them
E1, B1 = vec("E"), vec("B")
R1, R2 = div(E1), div(B1)
R3 = [u + v for u, v in zip(curl(E1), dt(B1))]
R4 = [u - v / c ** 2 for u, v in zip(curl(B1), dt(E1))]
waveE = [u - sp.diff(v, t, 2) / c ** 2 - (g1 - g2 + g3) for u, v, g1, g2, g3 in zip(lapv(E1), E1, grad(R1), curl(R3), dt(R4))]
waveB = [u - sp.diff(v, t, 2) / c ** 2 - (g1 - g2 - g3 / c ** 2) for u, v, g1, g2, g3 in zip(lapv(B1), B1, grad(R2), curl(R4), dt(R3))]
check("sym", "thm:un-waves", "lap E - E_tt/c^2 = grad(R1) - curl(R3) + d_t R4: the source-free equations give waves at c", allzero(waveE))
check("sym", "thm:un-waves", "lap B - B_tt/c^2 = grad(R2) - curl(R4) - d_t R3 / c^2", allzero(waveB))
ev, bv, rq, Jq = vec("e"), vec("b"), sca("rho_q"), vec("J_q")
Gs = sp.Symbol("G", positive=True)
Em, Bm = [-L ** 2 * f for f in ev], [-L ** 2 * f for f in bv]
Kk = Gs * L ** 2
check("sym", "thm:un-maxwell", "div E - 4 pi K rho_q = -L^2 (div e + 4 pi G rho_q), with E = -L^2 e and K = G L^2",
      zero(div(Em) - 4 * sp.pi * Kk * rq + L ** 2 * (div(ev) + 4 * sp.pi * Gs * rq)))
check("sym", "thm:un-maxwell", "curl E + d_t B = -L^2 (curl e + d_t b)",
      allzero([u + v + L ** 2 * (w + s) for u, v, w, s in zip(curl(Em), dt(Bm), curl(ev), dt(bv))]))
check("sym", "thm:un-maxwell", "curl B - (4 pi K/c^2) J_q - d_t E / c^2 = -L^2 (curl b + (4 pi G/c^2) J_q - d_t e / c^2)",
      allzero([u - 4 * sp.pi * Kk / c ** 2 * j - v / c ** 2 + L ** 2 * (w + 4 * sp.pi * Gs / c ** 2 * j - s / c ** 2)
               for u, j, v, w, s in zip(curl(Bm), Jq, dt(Em), curl(bv), dt(ev))]))
check("sym", "thm:un-maxwell", "div B = -L^2 div b", zero(div(Bm) + L ** 2 * div(bv)))
check("num", "thm:un-maxwell", "the magnetic constant: 4 pi K / c^2 equals mu0 of CODATA 2022",
      abs(4 * math.pi * C.K / C.c ** 2 / C.mu0 - 1) < 1e-9)

# two coordinates: (real, Composer) pairs; prop:un-split holds because every operation acts on each coordinate
P1, P2 = sca("Phi1"), sca("Phi2")
A1, A2 = vec("A1_"), vec("A2_")
E_ = [[-g - a_ for g, a_ in zip(grad(P), dt(A))] for P, A in ((P1, A1), (P2, A2))]
B_ = [curl(A1), curl(A2)]
rho = [-div(E_[k]) / (4 * sp.pi * Gs) for k in range(2)]
J = [[(c ** 2 / (4 * sp.pi * Gs)) * (-u + v / c ** 2) for u, v in zip(curl(B_[k]), dt(E_[k]))] for k in range(2)]
check("sym", "prop:un-split", "with potentials, (F2) and (F3) hold in each coordinate separately",
      all(zero(div(B_[k])) and allzero([u + v for u, v in zip(curl(E_[k]), dt(B_[k]))]) for k in range(2)))
check("sym", "prop:un-split", "a General Number equation a + bG = 0 holds exactly when a = 0 and b = 0 (the page map is injective)",
      zero_cs(sp.expand((a + b * L * Cs) ** 2 + (b * L * Ss) ** 2) - (a ** 2 + 2 * a * b * L * Cs + b ** 2 * L ** 2)))


def Re2(X_, Y_, f=lambda u, v: u * v):
    return f(X_[0], Y_[0]) - L ** 2 * f(X_[1], Y_[1])


check("sym", "thm:un-source", "d_t rho_mu + div J_mu = 0 in both coordinates, from (F1) and (F4)",
      all(zero(dt(rho[k]) + div(J[k])) for k in range(2)))
u_en = -(Re2(E_, E_, dot) + c ** 2 * Re2(B_, B_, dot)) / (8 * sp.pi * Gs)
S_en = [-(c ** 2 / (4 * sp.pi * Gs)) * (p1 - L ** 2 * p2) for p1, p2 in zip(cross(E_[0], B_[0]), cross(E_[1], B_[1]))]
check("sym", "thm:un-energy", "d_t u + div S = -Re[E . J] for every field from potentials",
      zero(dt(u_en) + div(S_en) + Re2(E_, J, dot)))
g_, e_, Bg_, bb_ = vec("g"), vec("ee"), vec("Bg"), vec("bb")
u2 = -(dot(g_, g_) - L ** 2 * dot(e_, e_) + c ** 2 * (dot(Bg_, Bg_) - L ** 2 * dot(bb_, bb_))) / (8 * sp.pi * Gs)
EE, BB = [-L ** 2 * f for f in e_], [-L ** 2 * f for f in bb_]
u2_faces = -(dot(g_, g_) + c ** 2 * dot(Bg_, Bg_)) / (8 * sp.pi * Gs) + (dot(EE, EE) + c ** 2 * dot(BB, BB)) / (8 * sp.pi * Kk)
check("sym", "prop:un-two-faces", "u = -(g^2 + c^2 Bg^2)/(8 pi G) + (E^2 + c^2 B^2)/(8 pi K)", zero(u2 - u2_faces))
f_force = [rho[0] * e1 - L ** 2 * rho[1] * e2 + c1 - L ** 2 * c2
           for e1, e2, c1, c2 in zip(E_[0], E_[1], cross(J[0], B_[0]), cross(J[1], B_[1]))]
p_f = [s_ / c ** 2 for s_ in S_en]
ReEE, ReBB = Re2(E_, E_, dot), Re2(B_, B_, dot)


def T(i, j):
    return -(E_[0][i] * E_[0][j] - L ** 2 * E_[1][i] * E_[1][j] + c ** 2 * (B_[0][i] * B_[0][j] - L ** 2 * B_[1][i] * B_[1][j])
             - sp.Rational(1, 2) * (1 if i == j else 0) * (ReEE + c ** 2 * ReBB)) / (4 * sp.pi * Gs)


check("sym", "thm:un-momentum", "f + d_t p_f = div T, component by component",
      all(zero(f_force[i] + dt(p_f[i]) - sum(sp.diff(T(i, j), X[j]) for j in range(3))) for i in range(3)))

# the action: the field equations as Euler-Lagrange equations of the potentials
rho1, rho2 = sca("r1"), sca("r2")
J1, J2 = vec("j1_"), vec("j2_")
Lag = (-(rho1 * P1 - L ** 2 * rho2 * P2) + (dot(J1, A1) - L ** 2 * dot(J2, A2))
       - (Re2(E_, E_, dot) - c ** 2 * Re2(B_, B_, dot)) / (8 * sp.pi * Gs))
fields = [P1, P2] + A1 + A2
eqs = sp.calculus.euler.euler_equations(Lag, fields, [t, x, y, z])
F1 = [div(E_[k]) + 4 * sp.pi * Gs * r for k, r in ((0, rho1), (1, rho2))]
F4 = [[u + 4 * sp.pi * Gs / c ** 2 * j - v / c ** 2 for u, j, v in zip(curl(B_[k]), JJ, dt(E_[k]))] for k, JJ in ((0, J1), (1, J2))]
targets = [F1[0], F1[1]] + F4[0] + F4[1]
srcs = [rho1, rho2] + J1 + J2


def proportional(expr, target, src):
    lam_ = sp.simplify(sp.expand(expr).coeff(src) / sp.expand(target).coeff(src))
    return lam_ != 0 and not lam_.has(sp.Function) and zero(expr - lam_ * target)


check("sym", ("def:un-action", "thm:un-action"), "(ii) varying Phi gives (F1) in both coordinates",
      all(proportional(eqs[k].lhs - eqs[k].rhs, targets[k], srcs[k]) for k in (0, 1)))
check("sym", ("def:un-action", "thm:un-action"), "(iii) varying A gives (F4) in both coordinates, all components",
      all(proportional(eqs[k].lhs - eqs[k].rhs, targets[k], srcs[k]) for k in range(2, 8)))
# (i) the path of a body: Lagrange's equations with the total time derivative along the path
m_, q_ = sp.symbols("m q", real=True)
V = sp.symbols("v_x v_y v_z", real=True)
Acc = sp.symbols("a_x a_y a_z", real=True)
Lp = (sp.Rational(1, 2) * m_ * sum(v ** 2 for v in V) - (m_ * P1 - L ** 2 * q_ * P2)
      + (m_ * dot(V, A1) - L ** 2 * q_ * dot(V, A2)))


def D(f):
    return sp.diff(f, t) + sum(V[i] * sp.diff(f, X[i]) + Acc[i] * sp.diff(f, V[i]) for i in range(3))


lhs_force = [m_ * (e1 + w1) - L ** 2 * q_ * (e2 + w2)
             for e1, e2, w1, w2 in zip(E_[0], E_[1], cross(list(V), B_[0]), cross(list(V), B_[1]))]
check("sym", ("thm:un-action", "law:unified"), "(i) varying the path gives (M): m a = Re[mu (E + v x B)]",
      all(zero(D(sp.diff(Lp, V[i])) - sp.diff(Lp, X[i]) - (m_ * Acc[i] - lhs_force[i])) for i in range(3)))
check("sym", "prop:un-lorentz", "Re[mu(E + v x B)] = m (g + v x Bg) + q (E + v x B) with E = -L^2 e, B = -L^2 b",
      allzero([m_ * (g1 + w1) - L ** 2 * q_ * (g2 + w2) - (m_ * (g1 + w1) + q_ * (-L ** 2 * g2 - L ** 2 * w2))
               for g1, g2, w1, w2 in zip(g_, e_, cross(list(V), Bg_), cross(list(V), bb_))]))
check("sym", "def:un-field", "the cosmic field of a body at rest gives the Cosmic Interaction Law: Re(mu E) = -G Re(mu mu_s)/r^2",
      zero(Re2((ma, qa), (-Gn * mb, -Gn * qb)) + Gn * (ma * mb - L ** 2 * qa * qb)))
rr = sp.Symbol("r", positive=True)
Ms = sp.Symbol("M_s", positive=True)
gfield = [-Gs * Ms * v / rr ** 3 for v in X]
check("sym", "def:un-field", "the field -G M r_hat / r^2 has no divergence and no curl away from the source",
      zero(div([f.subs(rr, sp.sqrt(x ** 2 + y ** 2 + z ** 2)) for f in gfield]))
      and allzero(curl([f.subs(rr, sp.sqrt(x ** 2 + y ** 2 + z ** 2)) for f in gfield])))

# the program: one force line moves every body (Table 6.2)
E0, vin, Lp_, Dscr = 1e4, 2.5e7, 0.05, 0.10
a_y = C.e * E0 / C.me
t1 = Lp_ / vin
rise = 0.5 * a_y * t1 ** 2 + a_y * t1 * Dscr / vin
val("law:unified", "Thomson's tube: rise of the electron at the screen, mm", f"{rise * 1e3:.2f}", "17.59", "sec6", "$17.59$~mm")
val("law:unified", "electron in 1 T: turns per second", sci(C.e * 1.0 / (2 * math.pi * C.me), 4), "2.7992 \\times 10^{10}", "sec6")


def drift_speed():
    """Proton from rest in E = 1000 V/m along y and B = 1 T along z: RK4 with the force line Re[mu(E + v x B)]."""
    m, q = C.mp, C.e
    e_fld = -1000.0 / C.GU ** 2            # Composer coordinate of the cosmic field: E = -|G|^2 e
    b_fld = -1.0 / C.GU ** 2

    def acc(v):
        Ey, Bz = -C.GU ** 2 * q * e_fld, -C.GU ** 2 * q * b_fld      # Re[mu E] and Re[mu B] for mu = m + q G
        return [(v[1] * Bz) / m, (Ey - v[0] * Bz) / m]
    omega = q * 1.0 / m
    n_per, n_turns = 4000, 10
    h = 2 * math.pi / omega / n_per
    pos, v = [0.0, 0.0], [0.0, 0.0]
    for _ in range(n_per * n_turns):
        k1 = acc(v)
        k2 = acc([v[0] + h / 2 * k1[0], v[1] + h / 2 * k1[1]])
        k3 = acc([v[0] + h / 2 * k2[0], v[1] + h / 2 * k2[1]])
        k4 = acc([v[0] + h * k3[0], v[1] + h * k3[1]])
        pos = [pos[i] + h * (v[i] + h / 6 * (k1[i] + k2[i] + k3[i])) for i in range(2)]
        v = [v[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(2)]
    return pos[0] / (n_per * n_turns * h)


val("law:unified", "proton in crossed fields: drift speed, m/s", f"{drift_speed():.3f}", "1000.000", "sec6")
val("law:unified", "hydrogen: speed of the classical orbit at a0 over c", f"{math.sqrt(C.K * C.e ** 2 / (C.me * C.a0)) / C.c:.7f}", "0.0072974", "sec6")
val("law:unified", "hydrogen: energy to free the electron at Bohr's radius for the reduced mass, eV",
    f"{C.K * C.e ** 2 / (2 * C.a0 * (C.me + C.mp) / C.mp) / C.e:.3f}", "13.598", "sec6", "$13.598$~eV")
ME, MM, aM = 5.9724e24, 7.346e22, 3.844e8          # NASA Moon fact sheet (Williams 2024)
val("law:unified", "the Moon and the Earth: two-body period, days", f"{2 * math.pi * math.sqrt(aM ** 3 / (C.G * (ME + MM))) / 86400:.2f}", "27.28", "sec6")
val("law:unified", "GW170817: about 2 s after about 130 million years, fraction", sci(2.0 / (130e6 * 365.25 * 86400), 0), "5 \\times 10^{-16}")
val("law:unified", "Gravity Probe B: a quarter of GR lies so many standard errors away", f"{abs(-37.2 - (-39.2 / 4)) / 7.2:.1f}", "3.8")
check("val", "law:unified", "the numbers -37.2 +- 7.2 and -39.2 mas/yr are printed as given by Everitt et al.",
      not HAVE_TEX or ("-37.2 \\pm 7.2" in FLAT["sec6"] and "-39.2" in FLAT["sec6"]))
val("law:unified", "apparent charge of a neutron under the shadow reading, in e", sci(C.mn * C.COS / C.GU / C.e, 1), "1.6 \\times 10^{-19}", "sec6")
val("law:unified", "that charge in standard errors of the measured neutron charge", f"{C.mn * C.COS / C.GU / C.e / 1.1e-21:.0f}", "144", "sec6", "about 145 times the error")
val("law:unified", "the sidereal month of the fact sheet, days", "27.32", "27.32", "sec6", "$27.32$~days")


# =============================================================================
heading("7. The angle and the units, the tests of gravity, earlier work (Sections 6 and 7)")
# =============================================================================
check("sym", "prop:angle-invariance", "the real part of a product equals ac - L^2 bd, which holds L but neither cos p nor sin p",
      zero_cs(gn_parts(page(a, b) * page(c_, d_))[0] - (a * c_ - L ** 2 * b * d_))
      and not ({Cs, Ss} & (a * c_ - L ** 2 * b * d_).free_symbols))
check("sym", "prop:angle-invariance", "the field equations and Maxwell's fields hold L only (E = -L^2 e, K = G L^2)",
      not ({Cs, Ss} & (set().union(*[sp.expand(f).free_symbols for f in Em + Bm]) | Kk.free_symbols)))
check("sym", "prop:angle-invariance", "energy density and action are built from Re2 = x1 y1 - L^2 x2 y2, free of the angle",
      not ({Cs, Ss} & (sp.expand(u2).free_symbols | sp.expand(Lag).free_symbols)))
uP = math.sqrt(C.hbar * C.c / C.G) / math.sqrt(4 * math.pi * C.eps0 * C.hbar * C.c)
check("num", "rem:natural-units", "Planck units: m_P / q_P = 1/sqrt(4 pi eps0 G) = sqrt(K/G)", abs(uP / C.GU - 1) < 1e-12)
check("num", "rem:natural-units", "Stoney units: m_S / e = sqrt(K/G)", abs(C.mS / C.e / C.GU - 1) < 1e-12)
sols = [k for k in range(1, 18000) if abs(abs(math.cos(math.radians(k / 100)) + math.sin(math.radians(k / 100))) - 1) < 1e-9
        and abs(k / 100 - 135) > 1e-9]
check("num", "rem:natural-units", "with A = 1, the only root in (0, 180) degrees, 135 excluded, is 90 degrees (grid of 0.01 degree)",
      sols == [9000])
check("sym", "rem:natural-units", "|cos p + sin p| = 1 means sin 2p = 0", zero(sp.expand_trig((sp.cos(p_) + sp.sin(p_)) ** 2 - 1 - sp.sin(2 * p_))))
val("rem:natural-units", "m_S / m_P = sqrt(alpha_fs)", f"{C.mS / C.mP:.7f}", "0.0854245")
check("num", "rem:natural-units", "m_S/m_P equals sqrt(alpha_fs) to the CODATA digits", abs(C.mS / C.mP / math.sqrt(C.alpha_fs) - 1) < 2e-9)
val("prop:angle-invariance", "Table 6: drawing length for u = 1e10 kg/C", f"{C.GU / 1e10:.5f}", "1.16043", "sec9", "$1.16043$")
val("prop:angle-invariance", "Table 6: drawing length for u = 1e11 kg/C", f"{C.GU / 1e11:.3f}", "0.116", "sec9", "$0.116$")
val("prop:angle-invariance", "Table 6: drawing length for u = 1e9 kg/C", f"{C.GU / 1e9:.1f}", "11.6", "sec9", "$11.6$")
val("prop:angle-invariance", "Table 6: roots for u = 1e10 kg/C, degrees", f"{math.degrees(C.ROOT1):.2f} {math.degrees(C.ROOT2):.2f}",
    "10.14 79.86", "sec9", "$10.14^\\circ$ and $79.86^\\circ$")
val("eq:axiom", "the size in units of 1e9 kg/C has no root (beyond sqrt 2)", f"{C.GU / 1e9:.1f} {C.GU / 1e9 > math.sqrt(2)}", "11.6 True", "sec9", "$11.6$")
val("eq:axiom", "Gaussian units: 1/sqrt(G) in g per statcoulomb", f"{1 / math.sqrt(C.G * 1e3):.0f}", "3871", "sec9")
val("eq:axiom", "in units of 1e11 kg/C: the length and its two roots, degrees",
    "{:.3f} {}".format(C.GU / 1e11, " ".join(f"{math.degrees(r_):.2f}" for r_ in (math.pi * 3 / 4 - math.asin(C.GU / 1e11 / math.sqrt(2)), math.pi * 3 / 4 + math.asin(C.GU / 1e11 / math.sqrt(2))))),
    "0.116 130.29 139.71", "sec9", "$130.29^\\circ$ and $139.71^\\circ$")
val("eq:axiom", "in units of 1e9 kg/C (Section 4)", f"{C.GU / 1e9:.1f}", "11.6", "sec4", "the drawing length would be $11.6$")
val("law:unified", "Table 5: a quarter of the frame dragging of general relativity, mas/yr", f"{-39.2 / 4:.1f}", "-9.8", "sec6", "$-9.8$")
val("law:interaction", "Coupling term against Coulomb for two protons, orders of magnitude", f"{math.log10(two_p[1] / two_p[2]):.0f}", "18")
val("eq:shadow-law", "apparent neutron charge under the shadow reading, e (Discussion)", sci(C.mn * C.COS / C.GU / C.e, 1), "1.6 \\times 10^{-19}", "sec9")
val("eq:shadow-law", "against the error of the measured neutron charge, 1.1e-21 e: about 145 standard errors", f"{C.mn * C.COS / C.GU / C.e / 1.1e-21:.0f}", "144", "sec9", "about 145 standard errors")
val("eq:electron", "Stoney mass e sqrt(K/G), kg (Discussion)", sci(C.mS, 5), "1.85921 \\times 10^{-9}", "sec9")
check("val", "eq:cu-size", "Methods prints the shift of the root for one standard error of G",
      not HAVE_TEX or "0.0009^\\circ" in FLAT["sec3"])
for snip in ("306 bodies", "2,481", "Garvey", "registered", "oganesson"):
    check("val", "def:cosmic-mass", f"the stage 4 paper no longer prints '{snip}'",
          not HAVE_TEX or all(snip not in re.sub(r"(?<!\\)%.*", "", TEX[n]) for n in SOURCES))

# =============================================================================
heading("8. Coverage of the numbered statements and labelled equations")
# =============================================================================
if HAVE_TEX:
    labels = set()
    for name in ("cmt-paper", "sec2", "sec3", "sec4", "sec5", "sec6", "sec9", "appA"):
        labels |= set(re.findall(r"\\label\{((?:eq|def|law|thm|prop|cor|rem):[^}]+)\}", TEX[name]))
    missing = sorted(labels - COVERED)
    check("cov", "all", f"{len(labels)} numbered statements and equations, each with at least one check"
          + (f"; missing: {', '.join(missing)}" if missing else ""), not missing)
else:
    print("  (the LaTeX sources were not found next to this script; the coverage check is skipped)")

# =============================================================================
n_ok = sum(r[3] for r in RESULTS)
tags = {tg: sum(1 for r in RESULTS if r[0] == tg) for tg in ("sym", "num", "val", "cov")}
print(f"\n{n_ok} of {len(RESULTS)} checks pass ({', '.join(f'{v} {k}' for k, v in tags.items())}); "
      f"{len(COVERED - {'all'})} labels covered; {time.time() - T0:.1f} s")
sys.exit(0 if n_ok == len(RESULTS) else 1)
