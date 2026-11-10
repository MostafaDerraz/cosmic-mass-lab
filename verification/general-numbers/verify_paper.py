#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_paper.py

Checks every numbered statement, every table row and every printed number of

    Mostafa Derraz, "General Numbers: The Number Plane Beyond the Right Angle",
    preprint, November 2026.

Each check carries one of four tags:
  sym  a symbolic proof with SymPy, valid for every angle of the Composer axis and,
       for the arithmetic, for every positive length of the Composer Unit;
  num  arrows multiplied directly on the page, as complex numbers, in relative
       planes at many angles and lengths, with random coefficients;
  ctl  a negative control: a habit of the complex plane (Table 1), copied into a
       leaning frame, must FAIL there and hold at pi/2;
  val  a number printed in the paper, recomputed from the definitions and rounded
       as printed;
  cov  coverage: every numbered statement of the paper has at least one check.

Notation. A General Number of the relative plane with angle p (phi-bar) and length
L = |G| is z = a + b G. On the page, G is the arrow of length L at the angle p, so z is
the complex number a + b L (cos p + i sin p). In the symbolic checks C = cos p and
S = sin p, with C^2 + S^2 = 1, and L is any positive length.

Run with:   python3 verify_paper.py [path/to/gn-paper.tex]
            (the path is optional; with it, the script also checks that every
            numbered statement of the paper is covered by at least one check)
Needs:      Python 3.8 or later, and SymPy (which brings mpmath).
Exit code:  0 when every check passes, 1 otherwise.
"""
import cmath
import math
import os
import random
import re
import sys
import time

import mpmath as mp
import sympy as sp

T0 = time.time()
RESULTS = []           # (tag, label, name, ok)
EVALS = [0]            # number of numerical comparisons made


def heading(title):
    print("\n" + title)
    print("-" * len(title))


def check(tag, label, name, ok):
    RESULTS.append((tag, label, name, bool(ok)))
    print(("  ok    " if ok else "  FAIL  ") + f"[{tag}] {label:22s} {name}")


# =============================================================================
# Symbolic tools
# =============================================================================
a, b, c, d, x, y, t = sp.symbols("a b c d x y t", real=True)
C, S = sp.symbols("C S", real=True)              # cos(phibar), sin(phibar)
L, M = sp.symbols("L M", positive=True)          # |G| (any length), a magnitude
I = sp.I
al, be, ph, th, pp, hh = sp.symbols("alpha beta phi theta p h", real=True)


def split_ri(expr):
    """Real and imaginary parts of an expanded expression in real symbols."""
    e = sp.expand(expr)
    im = sp.expand(e.coeff(I))
    re_ = sp.expand(e - I * im)
    return re_, im


def zero_CS(expr):
    """True if expr vanishes identically under C^2 + S^2 = 1 (denominators non-zero)."""
    num, _ = sp.fraction(sp.together(sp.expand(expr)))
    for part in split_ri(num):
        if part == 0:
            continue
        gens = [S, C] + sorted(part.free_symbols - {S, C}, key=lambda s_: s_.name)
        rem = sp.reduced(part, [S**2 + C**2 - 1], *gens, order="lex")[1]
        if sp.expand(rem) != 0:
            return False
    return True


def trig_zero(expr, angles):
    """True if a rational expression in sines and cosines of `angles` vanishes identically."""
    e = sp.expand_trig(sp.expand(expr))
    sub, rels, gens = {}, [], []
    for t_ in angles:
        s_t = sp.Symbol("s_" + t_.name, real=True)
        c_t = sp.Symbol("c_" + t_.name, real=True)
        sub[sp.sin(t_)] = s_t
        sub[sp.cos(t_)] = c_t
        rels.append(s_t**2 + c_t**2 - 1)
        gens += [s_t, c_t]
    e = e.subs(sub)
    num, _ = sp.fraction(sp.together(e))
    for part in split_ri(num):
        if part == 0:
            continue
        if part.has(sp.sin, sp.cos, sp.tan, sp.cot):
            raise ValueError("trigonometric term left unexpanded: " + str(part)[:120])
        others = sorted(part.free_symbols - set(gens), key=lambda s_: s_.name)
        rem = sp.reduced(part, rels, *gens, *others, order="lex")[1]
        if sp.expand(rem) != 0:
            return False
    return True


# General Number arithmetic, written in the leaning frame (the laws of the paper)
U = C + I * S                 # the direction u on the page
G = L * U                     # the Composer Unit on the page


def page(z):
    return z[0] + z[1] * G


def mul(z, w):                # Proposition 4.2
    (a1, b1), (a2, b2) = z, w
    return (a1 * a2 - L**2 * b1 * b2, a1 * b2 + b1 * a2 + 2 * L * C * b1 * b2)


def add(z, w):
    return (z[0] + w[0], z[1] + w[1])


def scal(k, z):
    return (k * z[0], k * z[1])


def conj(z):                  # Definition 4.4
    return (z[0] + 2 * z[1] * L * C, -z[1])


def Sh(z):                    # Definition 3.10
    return z[0] + z[1] * L * C


def Ht(z):
    return z[1] * L * S


def mag2(z):                  # Theorem 3.8(1)
    return z[0]**2 + 2 * z[0] * z[1] * L * C + z[1]**2 * L**2


def inv(z):                   # Theorem 4.7(2)
    return scal(1 / mag2(z), conj(z))


def quo(z, w):                # Theorem 4.7(3)
    den = mag2(w)
    a1, b1 = z
    a2, b2 = w
    return ((a1 * a2 + 2 * a1 * b2 * L * C + b1 * b2 * L**2) / den, (b1 * a2 - a1 * b2) / den)


Jc = (-C / S, 1 / (L * S))    # Definition 4.8
z_, w_ = (a, b), (c, d)
ONE, GU = (sp.Integer(1), sp.Integer(0)), (sp.Integer(0), sp.Integer(1))


def same(z, w):
    return zero_CS(z[0] - w[0]) and zero_CS(z[1] - w[1])


def gs_sin(al_, ph_):
    return sp.sin(al_) / sp.sin(al_ + ph_)


def gs_cos(al_, ph_):
    return sp.sin(ph_) / sp.sin(al_ + ph_)


def gs_tan(al_, ph_):
    return sp.sin(al_) / sp.sin(ph_)


def Cth(theta, p=pp):         # C_theta = cos(theta, p - theta)
    return gs_cos(theta, p - theta)


def Sth(theta, p=pp):         # S_theta = sin(theta, p - theta)
    return gs_sin(theta, p - theta)


# =============================================================================
# Numerical tools
# =============================================================================
rng = random.Random(20261109)


def Lax(p):
    return abs(math.cos(p) + math.sin(p))


class Plane:
    """A relative plane with angle p and length L (by default the Axiom's length)."""

    def __init__(self, p, Lv=None):
        self.p = p
        self.L = Lax(p) if Lv is None else Lv
        self.c, self.s = math.cos(p), math.sin(p)
        self.u = complex(self.c, self.s)
        self.G = self.L * self.u

    def page(self, z):
        return z[0] + z[1] * self.G

    def coords(self, X):
        yy = X.imag / (self.L * self.s)
        return (X.real - yy * self.L * self.c, yy)

    def mul(self, z, w):
        (a1, b1), (a2, b2) = z, w
        return (a1 * a2 - self.L**2 * b1 * b2, a1 * b2 + b1 * a2 + 2 * self.L * self.c * b1 * b2)

    def conj(self, z):
        return (z[0] + 2 * z[1] * self.L * self.c, -z[1])

    def Sh(self, z):
        return z[0] + z[1] * self.L * self.c

    def Ht(self, z):
        return z[1] * self.L * self.s

    def mag(self, z):
        return math.sqrt(z[0]**2 + 2 * z[0] * z[1] * self.L * self.c + z[1]**2 * self.L**2)

    def inv(self, z):
        m2 = self.mag(z)**2
        zc = self.conj(z)
        return (zc[0] / m2, zc[1] / m2)

    def J(self):
        return (-self.c / self.s, 1 / (self.L * self.s))

    def C_(self, th_):
        return math.sin(self.p - th_) / self.s

    def S_(self, th_):
        return math.sin(th_) / self.s

    def unit(self, th_):      # C + S u = C + (S/L) G
        return (self.C_(th_), self.S_(th_) / self.L)

    def matrix(self, z):      # Appendix A.7
        a1, b1 = z
        return mp.matrix([[a1, -self.L**2 * b1], [b1, a1 + 2 * b1 * self.L * self.c]])

    def exp(self, z):         # the exponential, computed natively from the matrix
        E = mp.expm(self.matrix(z))
        return (float(E[0, 0]), float(E[1, 0]))

    def power(self, z, n):
        r = (1.0, 0.0)
        base = z if n >= 0 else self.inv(z)
        for _ in range(abs(n)):
            r = self.mul(r, base)
        return r


def close(u, v, tol=1e-9):
    EVALS[0] += 1
    return abs(u - v) <= tol * (1 + abs(u) + abs(v))


def closez(z, w, tol=1e-9):
    return close(z[0], w[0], tol) and close(z[1], w[1], tol)


def rnum(lo=-3.0, hi=3.0):
    return (rng.uniform(lo, hi), rng.uniform(lo, hi))


DEG = math.pi / 180
ANGLES = [k * math.pi / 40 for k in range(1, 40) if k != 30]          # 3pi/4 excluded
ANGLES += [79.86 * DEG, 10.14 * DEG, math.pi / 3, 0.02, math.pi - 0.02,              # pi/2 is k = 20 above
           3 * math.pi / 4 - 0.05, 3 * math.pi / 4 + 0.05]
PLANES = [Plane(p) for p in ANGLES]                                   # lengths of the Axiom
PLANES_ANY = PLANES + [Plane(p, 0.37) for p in ANGLES[::3]] + [Plane(p, 2.6) for p in ANGLES[1::3]]
NUM_PLANES = len(PLANES_ANY)


def num_check(label, name, fn, planes=None, trials=6):
    planes = PLANES_ANY if planes is None else planes
    ok = True
    for P in planes:
        for _ in range(trials):
            if not fn(P):
                ok = False
                break
        if not ok:
            break
    check("num", label, name, ok)


# =============================================================================
heading("Section 2. The General Sinus")
# =============================================================================
lab = "def:gs, thm:gs-analytic"


def triangle_ok():
    al_ = rng.uniform(0.05, 3.0)
    ph_ = rng.uniform(0.05, math.pi - al_ - 0.05)
    # O at 0, Z on the second side at distance 1, projection line at the angle phi with ZO
    Z = cmath.exp(1j * al_)
    direc = cmath.exp(1j * (al_ + math.pi + ph_))
    tt = -Z.imag / direc.imag
    Zb = Z + tt * direc
    oppo, adja, base = abs(Zb - Z), abs(Zb), abs(Z)
    sin_g, cos_g, tan_g = oppo / base, adja / base, oppo / adja
    sin_a = math.sin(al_) / math.sin(al_ + ph_)
    cos_a = math.sin(ph_) / math.sin(al_ + ph_)
    tan_a = math.sin(al_) / math.sin(ph_)
    return abs(Zb.imag) < 1e-12 and close(sin_g, sin_a) and close(cos_g, cos_a) and close(tan_g, tan_a)


check("num", "def:gs", "the geometric ratios of 2000 oblique triangles", all(triangle_ok() for _ in range(2000)))
check("num", "thm:gs-analytic", "the analytic form equals the geometric form (law of sines)",
      all(triangle_ok() for _ in range(2000)))
check("sym", "cor:gs-classical", "classical recovery at alpha + phi = pi/2",
      trig_zero(gs_sin(al, sp.pi / 2 - al) - sp.sin(al), [al]) and
      trig_zero(gs_cos(al, sp.pi / 2 - al) - sp.cos(al), [al]) and
      trig_zero(gs_tan(al, sp.pi / 2 - al) - sp.sin(al) / sp.cos(al), [al]))
check("sym", "prop:gs-coordinates", "e^{i alpha} = cos(alpha, phi) + sin(alpha, phi) u, phi = p - alpha",
      trig_zero(sp.cos(al) + I * sp.sin(al) - (Cth(al) + Sth(al) * (sp.cos(pp) + I * sp.sin(pp))), [al, pp]))
check("sym", "prop:gs-basic", "reciprocity sin(a, f) = cos(f, a), tan(a, f) = cot(f, a)",
      trig_zero(gs_sin(al, ph) - gs_cos(ph, al), [al, ph]) and trig_zero(gs_tan(al, ph) - 1 / gs_tan(ph, al), [al, ph]))
check("sym", "prop:gs-basic", "symmetry identities",
      trig_zero(gs_sin(-al, -ph) - gs_sin(al, ph), [al, ph]) and trig_zero(gs_cos(-al, -ph) - gs_cos(al, ph), [al, ph])
      and trig_zero(gs_sin(-al, ph) - gs_sin(al, -ph), [al, ph]))
ok = True
for p_ in range(-2, 3):
    for q_ in range(-2, 3):
        ok &= trig_zero(gs_sin(al + p_ * sp.pi, ph + q_ * sp.pi) - (-1)**q_ * gs_sin(al, ph), [al, ph])
        ok &= trig_zero(gs_cos(al + p_ * sp.pi, ph + q_ * sp.pi) - (-1)**p_ * gs_cos(al, ph), [al, ph])
check("sym", "prop:gs-basic", "pi-addition identities for p, q = -2..2", ok)
check("num", "prop:gs-basic", "domains: sin, cos finite off alpha+phi in pi Z, tan off phi in pi Z",
      all(math.isfinite(math.sin(u_) / math.sin(u_ + v_)) and math.isfinite(math.sin(u_) / math.sin(v_))
          for u_, v_ in [(rng.uniform(-9, 9), rng.uniform(-9, 9)) for _ in range(500)]))
check("sym", "prop:gs-pythagoras", "sin^2 + cos^2 + 2 sin cos cos(alpha + phi) = 1",
      trig_zero(gs_sin(al, ph)**2 + gs_cos(al, ph)**2 + 2 * gs_sin(al, ph) * gs_cos(al, ph) * sp.cos(al + ph) - 1, [al, ph]))
xx, yy_, zz = sp.symbols("xx yy zz", real=True)
check("sym", "thm:gs-addition", "Ptolemy's identity",
      trig_zero(sp.sin(xx + yy_) * sp.sin(yy_ + zz) - sp.sin(xx) * sp.sin(zz) - sp.sin(yy_) * sp.sin(xx + yy_ + zz), [xx, yy_, zz]))
check("sym", "thm:gs-addition", "sin(a+b, f) = sin(a, f) + cos(a, f) sin(b, f + a)  [GSP (5.13)]",
      trig_zero(gs_sin(al + be, ph) - gs_sin(al, ph) - gs_cos(al, ph) * gs_sin(be, ph + al), [al, be, ph]))
check("sym", "thm:gs-addition", "cos(a+b, f) = cos cos + sin sin cos(a+b, f)  [GSP (5.15)]",
      trig_zero(gs_cos(al + be, ph) - gs_cos(al, ph) * gs_cos(be, ph) - gs_sin(al, ph) * gs_sin(be, ph) * gs_cos(al + be, ph), [al, be, ph]))
check("sym", "thm:gs-addition", "sin(a, f + t) = sin(a, f) sin(a + f, t)  [GSP (5.17)]",
      trig_zero(gs_sin(al, ph + th) - gs_sin(al, ph) * gs_sin(al + ph, th), [al, ph, th]))
check("sym", "thm:gs-addition", "cos(a, f + t) = cos(a, f) + sin(a, f) cos(a + f, t)  [GSP (5.18)]",
      trig_zero(gs_cos(al, ph + th) - gs_cos(al, ph) - gs_sin(al, ph) * gs_cos(al + ph, th), [al, ph, th]))
check("sym", "prop:gs-double", "sin(2a, f) = sin(a, f)[1 + cos(2a, f)]",
      trig_zero(gs_sin(2 * al, ph) - gs_sin(al, ph) * (1 + gs_cos(2 * al, ph)), [al, ph]))
check("sym", "prop:gs-double", "cos(2a, f) = cos^2(a, f) + sin^2(a, f) cos(2a, f)",
      trig_zero(gs_cos(2 * al, ph) - gs_cos(al, ph)**2 - gs_sin(al, ph)**2 * gs_cos(2 * al, ph), [al, ph]))
check("sym", "prop:gs-double", "the four first-order partial derivatives",
      trig_zero(sp.diff(gs_sin(al, ph), al) - gs_cos(al, ph) / sp.sin(al + ph), [al, ph]) and
      trig_zero(sp.diff(gs_sin(al, ph), ph) + gs_sin(al, ph) * sp.cos(al + ph) / sp.sin(al + ph), [al, ph]) and
      trig_zero(sp.diff(gs_cos(al, ph), al) + gs_cos(al, ph) * sp.cos(al + ph) / sp.sin(al + ph), [al, ph]) and
      trig_zero(sp.diff(gs_cos(al, ph), ph) - gs_sin(al, ph) / sp.sin(al + ph), [al, ph]))

mp.mp.dps = 40


def msin(u_, v_):
    return mp.sin(u_) / mp.sin(u_ + v_)


def mcos(u_, v_):
    return mp.sin(v_) / mp.sin(u_ + v_)


def mtan(u_, v_):
    return mp.sin(u_) / mp.sin(v_)


def products_ok(N, xv, yv):
    pi = mp.pi
    ps = mp.fprod(msin((xv + pi * k) / N, yv / N) for k in range(N))
    pc = mp.fprod(mcos(xv / N, (yv + pi * k) / N) for k in range(N))
    pt = mp.fprod(mtan((xv + pi * k) / N, (yv + pi * k) / N) for k in range(N))
    tol = mp.mpf(10)**-25
    return (abs(ps - msin(xv, yv)) < tol * (1 + abs(ps)) and abs(pc - mcos(xv, yv)) < tol * (1 + abs(pc))
            and abs(pt - mtan(xv, yv)) < tol * (1 + abs(pt)))


def random_xy():
    while True:
        xv, yv = mp.mpf(rng.uniform(-7, 7)), mp.mpf(rng.uniform(-7, 7))
        if abs(mp.sin(xv + yv)) > 0.05 and abs(mp.sin(yv)) > 0.05:
            return xv, yv


check("num", "thm:gs-transform", "the 2^n products for n = 1..5 (sin, cos, tan), 40 digits",
      all(products_ok(2**n, *random_xy()) for n in range(1, 6) for _ in range(20)))
check("num", "thm:every-n", "the products for every N = 1..24 (sin, cos, tan), 40 digits",
      all(products_ok(N, *random_xy()) for N in range(1, 25) for _ in range(20)))
check("num", "thm:every-n", "the classical product of sines, N = 1..24",
      all(abs(mp.fprod(mp.sin(tv + mp.pi * k / N) for k in range(N)) - mp.sin(N * tv) / 2**(N - 1)) < mp.mpf(10)**-30
          for N in range(1, 25) for tv in [mp.mpf(rng.uniform(-5, 5)) for _ in range(10)]))
Xs, Ys = sp.symbols("X Y", real=True)
ok = True
for N in (1, 2):          # x = N X, y = N Y; the factors have rational trigonometric constants
    ok &= trig_zero(gs_sin(N * Xs, N * Ys) - sp.Mul(*[gs_sin(Xs + sp.pi * k / N, Ys) for k in range(N)]), [Xs, Ys])
    ok &= trig_zero(gs_cos(N * Xs, N * Ys) - sp.Mul(*[gs_cos(Xs, Ys + sp.pi * k / N) for k in range(N)]), [Xs, Ys])
    ok &= trig_zero(gs_tan(N * Xs, N * Ys) - sp.Mul(*[gs_tan(Xs + sp.pi * k / N, Ys + sp.pi * k / N) for k in range(N)]), [Xs, Ys])
check("sym", "thm:every-n", "the three products for N = 1 and 2, symbolically", ok)
ok = True
for N in range(1, 13):
    for _ in range(20):
        xv, yv = random_xy()
        dens = [mp.sin((xv + mp.pi * k) / N + yv / N) for k in range(N)]
        ok &= abs(mp.fprod(dens) - mp.sin(xv + yv) / 2**(N - 1)) < mp.mpf(10)**-30 and all(abs(dv) > 0 for dv in dens)
check("num", "thm:every-n", "no factor leaves its domain (product of denominators)", ok)

# =============================================================================
heading("Section 3. The General Number Plane")
# =============================================================================
check("sym", "ax:gn", "|cos p + sin p| = sqrt2 |sin(p + pi/4)| (second form)",
      trig_zero(sp.cos(pp) + sp.sin(pp) - sp.sqrt(2) * sp.sin(pp + sp.pi / 4), [pp]))
grid = [k * math.pi / 20000 for k in range(20000)]
check("num", "ax:gn", "length 1 at 0 and pi/2, largest sqrt2 at pi/4, 0 at 3pi/4, back towards 1 near pi",
      close(Lax(0), 1) and close(Lax(math.pi / 2), 1) and close(Lax(math.pi / 4), math.sqrt(2))
      and abs(Lax(3 * math.pi / 4)) < 1e-15 and max(Lax(g_) for g_ in grid) <= math.sqrt(2) + 1e-15
      and abs(Lax(math.pi - 1e-7) - 1) < 1e-6)
P0 = Plane(0.0)
check("num", "ax:gn", "at the angle 0, 1 + 2 G and 2 + G are the same number, 3",
      close(P0.page((1, 2)), 3) and close(P0.page((2, 1)), 3))
num_check("def:relative-plane", "the sum is taken part by part",
          lambda P: (lambda z, w: close(P.page((z[0] + w[0], z[1] + w[1])), P.page(z) + P.page(w)))(rnum(), rnum()))
num_check("rem:page", "every arrow of the page has exactly one pair (a, b) in a true plane",
          lambda P: (lambda X: closez(P.coords(P.page(P.coords(X))), P.coords(X)) and close(P.page(P.coords(X)), X))(
              complex(*rnum())))
check("num", "thm:classical-recovery", "at pi/2 the Composer Unit is i", close(Plane(math.pi / 2).G, 1j))
num_check("prop:ladder", "the real numbers lie on the real axis of every plane",
          lambda P: (lambda r_: close(P.page((r_, 0.0)), r_))(rng.uniform(-9, 9)))

ok = True
for k in range(1, 400):
    A = 1 + (math.sqrt(2) - 1) * k / 400
    p1 = math.asin(A / math.sqrt(2)) - math.pi / 4
    p2 = math.pi / 2 - p1
    ok &= 0 < p1 < math.pi / 4 < p2 < math.pi / 2
    ok &= close(Lax(p1), A) and close(Lax(p2), A)
    ok &= close(math.cos(p2), math.sin(p1)) and close(math.sin(p2), math.cos(p1))
    ok &= close(math.sin(2 * p1), A * A - 1) and close(math.sin(2 * p2), A * A - 1) and close(p1, 0.5 * math.asin(A * A - 1))
    ok &= math.cos(2 * p1) > 0 and close(math.cos(2 * p2), -math.cos(2 * p1))
    vals = [Lax(g_) - A for g_ in grid]
    changes = sum(1 for i_ in range(len(vals) - 1) if vals[i_] == 0 or vals[i_] * vals[i_ + 1] < 0)
    ok &= changes == 2
check("num", "thm:two-roots", "399 lengths in (1, sqrt2): exactly two roots, statements (1) to (3)", ok)
ok = True
for k in range(1, 200):
    A = k / 200                                   # 0 < A < 1
    r1_, r2_ = 3 * math.pi / 4 - math.asin(A / math.sqrt(2)), 3 * math.pi / 4 + math.asin(A / math.sqrt(2))
    ok &= close(Lax(r1_), A) and close(Lax(r2_), A) and math.pi / 2 < r1_ < 3 * math.pi / 4 < r2_ < math.pi
    vals = [Lax(g_) - A for g_ in grid]
    ok &= sum(1 for i_ in range(len(vals) - 1) if vals[i_] * vals[i_ + 1] < 0) == 2
fine = [k * math.pi / 200000 for k in range(200000)]
ones = [g_ for g_ in fine if abs(Lax(g_) - 1) < 1e-12]
ok &= close(Lax(0), 1) and close(Lax(math.pi / 2), 1)
ok &= all(Lax(g_) < math.sqrt(2) - 1e-9 for g_ in fine if abs(g_ - math.pi / 4) > 1e-3)
ok &= all(Lax(g_) > 1e-9 for g_ in fine if abs(g_ - 3 * math.pi / 4) > 1e-3)
ok &= all(abs(Lax(g_) - 1) > 1e-9 for g_ in fine if min(abs(g_), abs(g_ - math.pi / 2)) > 1e-3)
check("num", "thm:two-roots", "the other lengths: A < 1 at 3pi/4 -+ arcsin(A/sqrt2); 1 at 0 and pi/2; sqrt2 at pi/4; 0 at 3pi/4", ok)
ok = all((Lax(g_)**2 * math.cos(2 * g_) < 0) == (math.pi / 4 < g_ < 3 * math.pi / 4) for g_ in grid
         if abs(g_ - 3 * math.pi / 4) > 1e-9 and abs(g_ - math.pi / 4) > 1e-9)
ok &= all(close(Plane(p_).G**2, Lax(p_)**2 * cmath.exp(2j * p_)) for p_ in ANGLES)
check("num", "thm:square-band", "the square is the arrow L^2 at 2p; the backward band is (pi/4, 3pi/4)", ok)
ok = True
for k in range(1, 400):
    A = 1 + (math.sqrt(2) - 1) * k / 400
    p1 = math.asin(A / math.sqrt(2)) - math.pi / 4
    p2 = math.pi / 2 - p1
    ok &= math.pi / 4 < p2 < 3 * math.pi / 4 and not (math.pi / 4 < p1 < 3 * math.pi / 4)
    ok &= close(A * A * math.cos(2 * p2), -A * A * math.cos(2 * p1))
check("num", "thm:square-band", "the mirror root is in the band, the arcsine root outside, shadows +-A^2 cos 2p1", ok)

check("sym", "thm:triangle", "(1) |z|^2 = a^2 + 2ab|G| cos p + b^2 |G|^2",
      zero_CS(sp.expand(page(z_) * page(z_).subs(I, -I)) - mag2(z_)))


def triangle_num(P):
    z = rnum()
    Z = P.page(z)
    m = abs(Z)
    alpha = cmath.phase(Z) + 2 * math.pi * rng.randint(-2, 2)          # any branch
    phi = P.p - alpha
    ok_ = close(z[0], m * math.sin(phi) / math.sin(alpha + phi)) and close(z[1] * P.L, m * math.sin(alpha) / math.sin(alpha + phi))
    eta = P.L * z[1] / z[0]
    if abs(1 + eta * P.c) > 1e-6:
        beta = math.atan(eta * P.s / (1 + eta * P.c))
        kk = (alpha - beta) / math.pi
        ok_ &= abs(kk - round(kk)) < 1e-9
    return ok_


num_check("thm:triangle", "(2) a = |z| cos(alpha, phi), b|G| = |z| sin(alpha, phi); (3) the phase angle", triangle_num, PLANES)


def corner_num(P):
    a1, b1 = rnum()
    v1, v2 = -a1, b1 * P.G                       # the two sides that meet at the corner
    ang = math.acos(max(-1, min(1, (v1 * v2.conjugate()).real / (abs(v1) * abs(v2)))))
    return close(ang, math.pi - P.p if a1 * b1 > 0 else P.p)


num_check("thm:triangle", "the corner angle is pi - p (same signs) or p (opposite signs)", corner_num, PLANES)
num_check("def:parts", "z is the sum of its two parts; equal numbers have equal parts",
          lambda P: (lambda z: close(P.page(z), z[0] + P.page((0, z[1]))) and closez(P.coords(P.page(z)), z))(rnum()))
num_check("def:shadow-height", "Sh and Ht are the perpendicular readings on the page",
          lambda P: (lambda z: close(P.Sh(z), P.page(z).real) and close(P.Ht(z), P.page(z).imag))(rnum()))

# =============================================================================
heading("Section 4. Arithmetic in a leaning frame")
# =============================================================================
check("sym", "prop:square-rule", "G^2 = 2|G| cos p G - |G|^2 and u^2 = 2 cos p u - 1 (any length)",
      zero_CS(G**2 - (2 * L * C * G - L**2)) and zero_CS(U**2 - (2 * C * U - 1)))
check("sym", "prop:product", "(1) the product formula equals Wessel's product on the page",
      zero_CS(page(mul(z_, w_)) - page(z_) * page(w_)))
check("sym", "prop:product", "(2) the real part a c - |G|^2 b d contains neither cos p nor sin p",
      not (mul(z_, w_)[0].has(C) or mul(z_, w_)[0].has(S)))
X_ = (x, y)
check("sym", "prop:product", "(3) Re(1 X) = x and Re(G X) = -|G|^2 y, so Re(wX) = 0 for all w forces X = 0",
      sp.expand(mul(ONE, X_)[0] - x) == 0 and sp.expand(mul(GU, X_)[0] + L**2 * y) == 0)


def tworoot_real(P):
    if not (1 < P.L < math.sqrt(2)) or P.p >= math.pi / 2:
        return True
    Q2 = Plane(math.pi / 2 - P.p)                 # the other root of the same length
    z, w = rnum(), rnum()
    return close(P.mul(z, w)[0], Q2.mul(z, w)[0])


num_check("prop:product", "(2) the real part of zw is the same in the two planes of one length", tworoot_real, PLANES)


def field_num(P):
    z, w, v = rnum(), rnum(), rnum()
    ok_ = closez(P.mul(z, w), P.mul(w, z))
    ok_ &= closez(P.mul(P.mul(z, w), v), P.mul(z, P.mul(w, v)))
    ok_ &= closez(P.mul(z, (w[0] + v[0], w[1] + v[1])), tuple(p_ + q_ for p_, q_ in zip(P.mul(z, w), P.mul(z, v))))
    ok_ &= closez(P.mul((1.0, 0.0), z), z) and closez((z[0] - z[0], z[1] - z[1]), (0, 0))
    ok_ &= closez(P.mul(z, P.inv(z)), (1.0, 0.0))
    return ok_


num_check("thm:field", "order, grouping, distribution, unit, opposite and inverse", field_num)
check("sym", "def:conjugate, thm:conjugate", "(1) the conjugate is the mirror image of z in the real axis",
      zero_CS(page(conj(z_)) - page(z_).subs(I, -I)))
check("sym", "thm:conjugate", "(2) conj conj z = z, conj(z + w), conj(z w) = conj z conj w",
      same(conj(conj(z_)), z_) and same(conj(add(z_, w_)), add(conj(z_), conj(w_)))
      and same(conj(mul(z_, w_)), mul(conj(z_), conj(w_))))
check("sym", "thm:conjugate", "(3) z + conj z = 2 Sh(z)", same(add(z_, conj(z_)), (2 * Sh(z_), 0)))
check("sym", "thm:conjugate", "(4) z conj z = |z|^2, a positive definite form",
      same(mul(z_, conj(z_)), (mag2(z_), 0)) and zero_CS((2 * L * C)**2 - 4 * L**2 + 4 * L**2 * S**2))
check("sym", "thm:conjugate", "(5) G + conj G = 2|G| cos p, G conj G = |G|^2, roots of x^2 - 2|G|cos p x + |G|^2",
      same(add(GU, conj(GU)), (2 * L * C, 0)) and same(mul(GU, conj(GU)), (L**2, 0))
      and zero_CS(G**2 - 2 * L * C * G + L**2) and zero_CS(G.subs(I, -I)**2 - 2 * L * C * G.subs(I, -I) + L**2))
check("sym", "rem:skew-flip", "|a - bG|^2 = |z|^2 - 4ab|G| cos p and (a + bG)(a - bG)",
      zero_CS(mag2((a, -b)) - (mag2(z_) - 4 * a * b * L * C))
      and same(mul(z_, (a, -b)), (a**2 + b**2 * L**2, -2 * b**2 * L * C)))
check("num", "rem:skew-flip", "at pi/2 the skew flip is the conjugate",
      closez(Plane(math.pi / 2).conj((1.3, 0.7)), (1.3, -0.7)))
p_, q_ = mul(z_, w_)
check("sym", "thm:division", "(1) the two-squares identity and |z w| = |z| |w|",
      zero_CS(mag2(z_) * mag2(w_) - (p_**2 + 2 * p_ * q_ * L * C + q_**2 * L**2)))
check("sym", "thm:division", "(2) z (1/z) = 1 and 1/G = (2|G| cos p - G)/|G|^2",
      same(mul(z_, inv(z_)), ONE) and same(inv(GU), (2 * C / L, -1 / L**2)))
check("sym", "thm:division", "(3) w (z/w) = z", same(mul(w_, quo(z_, w_)), z_))
check("sym", "def:turning-unit", "(u - cos p)/sin p = -cot p + G/(|G| sin p)",
      zero_CS((U - C) / S - (-C / S + G / (L * S))))
check("sym", "thm:turning-unit", "(1) J is the arrow i of the page: length 1, a right angle",
      zero_CS(page(Jc) - I))
check("sym", "thm:turning-unit", "(2) J^2 = -1, and x^2 + 1 = (x - J)(x + J), so only +-J",
      same(mul(Jc, Jc), (-1, 0)) and same(mul(add(X_, scal(-1, Jc)), add(X_, Jc)), add(mul(X_, X_), ONE)))
check("sym", "thm:turning-unit", "(3) J z is z turned by a right angle, same length",
      zero_CS(page(mul(Jc, z_)) - I * page(z_)))
check("num", "thm:turning-unit", "(4) at pi/2, J = G = i",
      closez(Plane(math.pi / 2).J(), (0, 1)) and close(Plane(math.pi / 2).G, 1j))
check("sym", "prop:shadow-height", "z = Sh + Ht J, a = Sh - Ht cot p, b = Ht/(|G| sin p)",
      same(add((Sh(z_), 0), scal(Ht(z_), Jc)), z_) and zero_CS(Sh(z_) - Ht(z_) * C / S - a)
      and zero_CS(Ht(z_) / (L * S) - b))
check("sym", "prop:shadow-height", "z - conj z = 2 Ht J and |z|^2 = Sh^2 + Ht^2",
      same(add(z_, scal(-1, conj(z_))), scal(2 * Ht(z_), Jc)) and zero_CS(mag2(z_) - Sh(z_)**2 - Ht(z_)**2))


def Phi(z):
    return Sh(z) + I * Ht(z)


xc, yc = sp.symbols("xc yc", real=True)
Phi_inv = (xc - yc * C / S, yc / (L * S))
check("sym", "cor:bridge", "Phi is additive and multiplicative; Phi(G) = |G|e^{ip}, Phi(J) = i",
      zero_CS(Phi(add(z_, w_)) - Phi(z_) - Phi(w_)) and zero_CS(Phi(mul(z_, w_)) - Phi(z_) * Phi(w_))
      and zero_CS(Phi(GU) - G) and zero_CS(Phi(Jc) - I))
check("sym", "cor:bridge", "Phi keeps magnitudes and conjugates; the inverse formula",
      zero_CS(Phi(z_) * Phi(z_).subs(I, -I) - mag2(z_)) and zero_CS(Phi(conj(z_)) - Phi(z_).subs(I, -I))
      and zero_CS(Phi(Phi_inv) - (xc + I * yc)))
check("sym", "prop:pairs", "(1) the parts from the readings, and back",
      zero_CS(Sh(z_) - Ht(z_) * C / S - a) and zero_CS(Ht(z_) / (L * S) - b) and zero_CS(Sh(z_) - (a + b * L * C)))
kk_ = sp.symbols("kk", real=True)
check("sym", "prop:pairs", "(2) the four maps are additive and commute with real factors",
      zero_CS(Sh(add(z_, w_)) - Sh(z_) - Sh(w_)) and zero_CS(Ht(add(z_, w_)) - Ht(z_) - Ht(w_))
      and zero_CS(Sh(scal(kk_, z_)) - kk_ * Sh(z_)) and zero_CS(Ht(scal(kk_, z_)) - kk_ * Ht(z_)))
Gb = G.subs(I, -I)
zp, zbp = page(z_), page(z_).subs(I, -I)
check("sym", "prop:pairs", "(3) conj z = a + b conj G; Re, Com, Sh, Ht J read from z and conj z",
      zero_CS(zbp - (a + b * Gb)) and zero_CS((G * zbp - Gb * zp) - a * (G - Gb))
      and zero_CS(G * (zp - zbp) - b * G * (G - Gb)) and zero_CS((zp + zbp) / 2 - Sh(z_))
      and zero_CS((zp - zbp) / 2 - Ht(z_) * I) and zero_CS(G - Gb - 2 * L * S * I))


def pairs4(P):
    r_, s_ = rng.uniform(-3, 3), rng.uniform(-3, 3)
    real, pure = (r_, 0.0), (0.0, s_)
    gen = rnum()
    ok_ = closez(P.conj(real), real) and not closez(P.conj(gen), gen)
    Gp_, Gbp = P.G, P.G.conjugate()
    ok_ &= close(Gp_ * P.page(P.conj(pure)), Gbp * P.page(pure))
    ok_ &= not close(Gp_ * P.page(P.conj(gen)), Gbp * P.page(gen))
    return ok_


num_check("prop:pairs", "(4) real exactly when conj z = z; pure Composer exactly when G conj z = conj G z", pairs4, PLANES)
check("sym", "thm:parts-product", "(1) Re(zw) = Re Re - |G|^2 bd; Com(zw) = Sh(z)Com(w) + Sh(w)Com(z)",
      zero_CS(mul(z_, w_)[0] - (a * c - L**2 * b * d)) and zero_CS(mul(z_, w_)[1] - (Sh(z_) * d + Sh(w_) * b)))
check("sym", "thm:parts-product", "(2) Sh(zw) = Sh Sh - Ht Ht; Ht(zw) = Sh Ht + Sh Ht",
      zero_CS(Sh(mul(z_, w_)) - (Sh(z_) * Sh(w_) - Ht(z_) * Ht(w_)))
      and zero_CS(Ht(mul(z_, w_)) - (Sh(z_) * Ht(w_) + Sh(w_) * Ht(z_))))
z2 = mul(z_, z_)
check("sym", "thm:parts-product", "(3) the four laws for a square",
      zero_CS(z2[0] - (a**2 - L**2 * b**2)) and zero_CS(z2[1] - 2 * Sh(z_) * b)
      and zero_CS(Sh(z2) - (Sh(z_)**2 - Ht(z_)**2)) and zero_CS(Ht(z2) - 2 * Sh(z_) * Ht(z_)))
zc_, zi_ = conj(z_), inv(z_)
check("sym", "prop:parts-conjugate", "the parts, shadow and height of conj z and of 1/z",
      zero_CS(zc_[0] - (2 * Sh(z_) - a)) and zero_CS(zc_[1] + b) and zero_CS(Sh(zc_) - Sh(z_)) and zero_CS(Ht(zc_) + Ht(z_))
      and zero_CS(zi_[0] - (2 * Sh(z_) - a) / mag2(z_)) and zero_CS(zi_[1] + b / mag2(z_))
      and zero_CS(Sh(zi_) - Sh(z_) / mag2(z_)) and zero_CS(Ht(zi_) + Ht(z_) / mag2(z_)))
zw_ = mul(conj(z_), w_)
check("sym", "prop:signed-area", "the four laws of conj(z) w",
      zero_CS(zw_[1] - (a * d - b * c)) and zero_CS(Ht(zw_) - (Sh(z_) * Ht(w_) - Ht(z_) * Sh(w_)))
      and zero_CS(zw_[0] - (conj(z_)[0] * c + L**2 * b * d)) and zero_CS(Sh(zw_) - (Sh(z_) * Sh(w_) + Ht(z_) * Ht(w_))))
zpg, wpg = page(z_), page(w_)
check("sym", "prop:signed-area", "Ht(conj z w) = (ad - bc)|G| sin p = the cross product of the arrows",
      zero_CS(Ht(zw_) - (a * d - b * c) * L * S)
      and zero_CS(Ht(zw_) - (split_ri(zpg)[0] * split_ri(wpg)[1] - split_ri(zpg)[1] * split_ri(wpg)[0])))


def area_num(P):
    z, w = rnum(), rnum()
    Z, W = P.page(z), P.page(w)
    zw = P.mul(P.conj(z), w)
    ok_ = close(P.Sh(zw), abs(Z) * abs(W) * math.cos(cmath.phase(W) - cmath.phase(Z)))
    ok_ &= (P.Ht(zw) > 0) == ((cmath.phase(W / Z)) > 0)
    tt = rng.uniform(-2, 2)
    col = (tt * z[0], tt * z[1])
    ok_ &= close(z[0] * col[1] - z[1] * col[0], 0)
    perp = P.mul(P.J(), (tt * z[0], tt * z[1]))
    ok_ &= close(P.Sh(P.mul(P.conj(z), perp)), 0)
    return ok_


num_check("prop:signed-area", "Sh(conj z w) = |z||w| cos; sign of the area; on one line; at a right angle", area_num, PLANES)


def longer_num(P):
    z = rnum()
    m = P.mag(z)
    ok_ = abs(z[0]) <= m / P.s * (1 + 1e-12) and abs(z[1]) * P.L <= m / P.s * (1 + 1e-12)
    ok_ &= abs(P.Sh(z)) <= m * (1 + 1e-12) and abs(P.Ht(z)) <= m * (1 + 1e-12)
    zr = P.coords(cmath.exp(1j * (P.p - math.pi / 2)))         # at a right angle to the Composer axis
    ok_ &= close(abs(zr[0]), 1 / P.s)
    Jn = P.J()
    ok_ &= close(abs(Jn[1]) * P.L, 1 / P.s) and close(P.mag(Jn), 1)
    return ok_


num_check("rem:longer", "|Re z|, |Com z| <= |z|/sin p, both bounds reached (the second by J); readings <= |z|", longer_num, PLANES)

# ---------- Table 1: negative controls ----------
TZ = [(1.3, 0.7), (-0.4, 1.1), (2.0, -1.5), (0.0, 1.0), (1.0, 1.0)]


def habits(P):
    """For each habit: (does the habit, copied with i -> G and Im -> b, hold for every test number?)"""
    out = {}
    out["i^2 = -1"] = close(P.G**2, -1)
    out["Re conj z = Re z"] = all(close(P.conj(z)[0], z[0]) for z in TZ)
    out["z + conj z = 2 Re z"] = all(close(P.page(z) + P.page(P.conj(z)), 2 * z[0]) for z in TZ)
    out["z conj z = a^2 + b^2"] = all(close(P.mag(z)**2, z[0]**2 + z[1]**2) for z in TZ)
    out["(a + bG)(a - bG) is real"] = all(close(P.mul(z, (z[0], -z[1]))[1], 0) for z in TZ)
    out["Com coordinate of zw = ad + bc"] = all(close(P.mul(z, w)[1], z[0] * w[1] + z[1] * w[0]) for z in TZ for w in TZ)
    zr = P.coords(cmath.exp(1j * (P.p - math.pi / 2)))
    out["|Re z| <= |z|"] = all(abs(z[0]) <= P.mag(z) * (1 + 1e-12) for z in TZ + [zr])
    out["the Composer part is as long as the height"] = all(close(abs(z[1]) * P.L, abs(P.Ht(z))) for z in TZ)
    out["e^{tG} = cos t + G sin t"] = all(close(cmath.exp(tv * P.G), math.cos(tv) + P.G * math.sin(tv)) for tv in (0.3, 1.1, 2.5))
    return out


def laws(P):
    z, w = TZ[0], TZ[1]
    zr = P.coords(cmath.exp(1j * (P.p - math.pi / 2)))
    tv = 0.9
    J = P.J()
    eJ = P.exp((tv * J[0], tv * J[1]))
    return all([
        close(P.G**2, 2 * P.L * P.c * P.G - P.L**2) and closez(P.mul(J, J), (-1, 0)),
        close(P.conj(z)[0], 2 * P.Sh(z) - z[0]),
        close(P.page(z) + P.page(P.conj(z)), 2 * P.Sh(z)),
        close(P.mag(z)**2, P.Sh(z)**2 + P.Ht(z)**2),
        close(P.mul(z, (z[0], -z[1]))[1], -2 * z[1]**2 * P.L * P.c),
        close(P.mul(z, w)[1], P.Sh(z) * w[1] + P.Sh(w) * z[1]),
        abs(zr[0]) <= P.mag(zr) / P.s * (1 + 1e-12),
        close(z[1] * P.G, P.Ht(z) / P.s * P.u),
        closez(eJ, P.unit(tv)),
    ])


lean = [P for P in PLANES if abs(P.c) > 1e-9]
P90 = Plane(math.pi / 2)
H90 = habits(P90)
for name in H90:
    check("ctl", "tab:habits", f"habit fails in all {len(lean)} leaning planes, holds at pi/2: {name}",
          H90[name] and all(not habits(P)[name] for P in lean))
check("num", "tab:habits", "every law of the right-hand column holds in every plane", all(laws(P) for P in PLANES))

# =============================================================================
heading("Section 5. Turning, powers and roots")
# =============================================================================
check("sym", "thm:triangle, prop:gs-pythagoras", "C^2 + 2 cos p C S + S^2 = 1 for C_t, S_t",
      trig_zero(Cth(th)**2 + 2 * sp.cos(pp) * Cth(th) * Sth(th) + Sth(th)**2 - 1, [th, pp]))
check("sym", "thm:turning", "C_{a+b} = C_a C_b - S_a S_b",
      trig_zero(Cth(al + be) - (Cth(al) * Cth(be) - Sth(al) * Sth(be)), [al, be, pp]))
check("sym", "thm:turning", "S_{a+b} = C_a S_b + S_a C_b + 2 cos p S_a S_b",
      trig_zero(Sth(al + be) - (Cth(al) * Sth(be) + Sth(al) * Cth(be) + 2 * sp.cos(pp) * Sth(al) * Sth(be)), [al, be, pp]))


def turning_num(P):
    t1, t2 = rng.uniform(-7, 7), rng.uniform(-7, 7)
    r1, r2 = rng.uniform(0.2, 3), rng.uniform(0.2, 3)
    z = (r1 * P.unit(t1)[0], r1 * P.unit(t1)[1])
    w = (r2 * P.unit(t2)[0], r2 * P.unit(t2)[1])
    pr = P.mul(z, w)
    qu = P.mul(z, P.inv(w))
    return (closez(pr, (r1 * r2 * P.unit(t1 + t2)[0], r1 * r2 * P.unit(t1 + t2)[1]))
            and closez(qu, (r1 / r2 * P.unit(t1 - t2)[0], r1 / r2 * P.unit(t1 - t2)[1])))


num_check("thm:turning", "products and quotients in polar form", turning_num)


def euler_num(P):
    tv = rng.uniform(-8, 8)
    J = P.J()
    ok_ = closez(P.exp((tv * J[0], tv * J[1])), P.unit(tv), 1e-8)
    z = rnum()
    m = P.mag(z)
    alpha = cmath.phase(P.page(z))
    ok_ &= closez(z, tuple(m * v for v in P.exp((alpha * J[0], alpha * J[1]))), 1e-8)
    ok_ &= closez(P.conj(z), tuple(m * v for v in P.exp((-alpha * J[0], -alpha * J[1]))), 1e-8)
    ok_ &= closez((0.0, 1.0), tuple(P.L * v for v in P.exp((P.p * J[0], P.p * J[1]))), 1e-8)
    ok_ &= closez(J, P.exp((math.pi / 2 * J[0], math.pi / 2 * J[1])), 1e-8)
    ok_ &= closez(P.exp((math.pi * J[0], math.pi * J[1])), (-1.0, 0.0), 1e-8)
    return ok_


mp.mp.dps = 30
num_check("thm:euler", "e^{tJ} = C_t + S_t u (matrix exponential); z = |z|e^{aJ}; G, J, e^{pi J} = -1", euler_num, PLANES, trials=3)
check("sym", "rem:euler-native", "dS/dt = C/sin p + S cot p and dC/dt = -C cot p - S/sin p",
      trig_zero(sp.diff(Sth(th), th) - (Cth(th) / sp.sin(pp) + Sth(th) * sp.cos(pp) / sp.sin(pp)), [th, pp])
      and trig_zero(sp.diff(Cth(th), th) - (-Cth(th) * sp.cos(pp) / sp.sin(pp) - Sth(th) / sp.sin(pp)), [th, pp]))
Ls = sp.symbols("Ls", positive=True)
Jt = (-sp.cos(pp) / sp.sin(pp), 1 / (Ls * sp.sin(pp)))
wt = (Cth(th), Sth(th) / Ls)
Jw = (Jt[0] * wt[0] - Ls**2 * Jt[1] * wt[1], Jt[0] * wt[1] + Jt[1] * wt[0] + 2 * Ls * sp.cos(pp) * Jt[1] * wt[1])
check("sym", "rem:euler-native", "w(t) = C_t + S_t u obeys dw/dt = J w (in the frame {1, G})",
      trig_zero(sp.diff(wt[0], th) - Jw[0], [th, pp]) and trig_zero(sp.diff(wt[1], th) - Jw[1], [th, pp]))


def demoivre_num(P):
    tv = rng.uniform(-4, 4)
    n = rng.randint(-7, 7)
    return closez(P.power(P.unit(tv), n), P.unit(n * tv), 1e-8)


num_check("thm:de-moivre", "(C_t + S_t u)^n = C_{nt} + S_{nt} u for n = -7..7", demoivre_num)
ok = True
for n in range(-3, 8):
    Gn = (L * (sp.cos(pp) + I * sp.sin(pp)))**n if n >= 0 else (sp.cos(pp) - I * sp.sin(pp))**(-n) / L**(-n)
    formula = -L**n * sp.sin((n - 1) * pp) / sp.sin(pp) + L**(n - 1) * sp.sin(n * pp) / sp.sin(pp) * L * (sp.cos(pp) + I * sp.sin(pp))
    ok &= trig_zero(Gn - formula, [pp])
check("sym", "cor:powers", "G^n = -|G|^n sin(n-1)p/sin p + |G|^{n-1} sin np/sin p G, n = -3..7", ok)
ok = all(zero_CS(G**(n + 1) - (2 * L * C * G**n - L**2 * G**(n - 1))) for n in range(1, 8))
check("sym", "cor:powers", "the step rule, n = 1..7", ok)
ok = all(trig_zero(sp.sin(n * pp) / sp.sin(pp) - sp.chebyshevu(n - 1, sp.cos(pp)), [pp]) for n in range(1, 9))
check("sym", "cor:powers", "sin np/sin p = U_{n-1}(cos p), n = 1..8", ok)


def spiral_num(P):
    """Proposition 5.6: the powers of G lie on rho = |G|^(theta/pbar), at the angles n*pbar."""
    ok_ = True
    for n in range(-5, 11):
        gn = P.power((0.0, 1.0), n)
        Z = P.page(gn)
        ok_ &= close(Z, P.L**n * cmath.exp(1j * n * P.p), 1e-8)
        ok_ &= close(abs(Z), P.L**((n * P.p) / P.p), 1e-8)                 # on the spiral at theta = n pbar
    return ok_


num_check("prop:spiral", "G^n = |G|^n e^{i n pbar} for n = -5..10, on the spiral rho = |G|^(theta/pbar)", spiral_num, PLANES, trials=1)
fine_ = [k * math.pi / 100000 for k in range(1, 100000)]
ok = all((Lax(g_) > 1) == (0 < g_ < math.pi / 2) for g_ in fine_ if abs(g_ - math.pi / 2) > 1e-9)
ok &= close(Lax(math.pi / 2), 1) and all(Lax(g_) < 1 for g_ in fine_ if g_ > math.pi / 2 + 1e-9)
check("num", "prop:spiral", "(1)-(3) |G| > 1 exactly on (0, pi/2), = 1 at pi/2, < 1 on (pi/2, pi)", ok)
ok = True
for pdeg, N in ((60, 6), (45, 8), (30, 12), (90, 4), (120, 3), (105, 24), (72, 5), (150, 12)):
    Pq = Plane(math.radians(pdeg))
    gN = Pq.power((0.0, 1.0), N)
    ok &= abs(gN[1]) < 1e-9 * (1 + abs(gN[0])) and gN[0] > 0 and close(gN[0], Pq.L**N, 1e-9)
check("num", "prop:spiral", "(4) N pbar = 2 pi k gives a positive real G^N (eight planes)", ok)
check("sym", "prop:spiral", "exact returns: G^6 = (26 + 15 sqrt3)/8 at pi/3; G^8 = 16 at pi/4; |G| = sqrt2/2 at 7pi/12",
      sp.simplify(sp.expand(((sp.cos(sp.pi / 3) + sp.sin(sp.pi / 3)) * (sp.cos(sp.pi / 3) + I * sp.sin(sp.pi / 3)))**6) - (26 + 15 * sp.sqrt(3)) / 8) == 0
      and sp.simplify(sp.expand(((sp.cos(sp.pi / 4) + sp.sin(sp.pi / 4)) * (sp.cos(sp.pi / 4) + I * sp.sin(sp.pi / 4)))**8) - 16) == 0
      and sp.simplify(sp.Abs(sp.cos(7 * sp.pi / 12) + sp.sin(7 * sp.pi / 12)) - sp.sqrt(2) / 2) == 0)


# ---------- Remark 5.7: the golden ratio at pi/5 ----------
mp.mp.dps = 40
tau = (1 + mp.sqrt(5)) / 2
p5 = mp.pi / 5
L5 = mp.cos(p5) + mp.sin(p5)
u5 = mp.expj(p5)
G5 = L5 * u5
ok = abs(2 * mp.cos(p5) - tau) < mp.mpf(10)**-35 and abs(u5**2 - (tau * u5 - 1)) < mp.mpf(10)**-35
ok &= abs(G5**2 - (tau * L5 * G5 - L5**2)) < mp.mpf(10)**-35
cn = [mp.sin(n * p5) / mp.sin(p5) for n in range(1, 11)]
ok &= all(abs(x_ - y_) < mp.mpf(10)**-35 for x_, y_ in zip(cn, [1, tau, tau, 1, 0, -1, -tau, -tau, -1, 0]))
ok &= abs(G5**2 - (-L5**2 + tau * L5 * G5)) < mp.mpf(10)**-35
ok &= abs(G5**3 - (-tau * L5**3 + tau * L5**2 * G5)) < mp.mpf(10)**-35
ok &= abs(G5**4 - (-tau * L5**4 + L5**3 * G5)) < mp.mpf(10)**-35
ok &= abs(G5**5 + L5**5) < mp.mpf(10)**-35 and abs(G5**10 - L5**10) < mp.mpf(10)**-33
check("num", "rem:golden", "2cos(pi/5) = tau; u^2 = tau u - 1; coefficients 1, tau, tau, 1, 0; G^2..G^5; G^10 real (40 digits)", ok)
q5 = 2 * mp.pi / 5
C5, S5 = mp.sin(q5 - p5) / mp.sin(q5), mp.sin(p5) / mp.sin(q5)
P25 = Plane(float(q5))
eJ = P25.exp((float(p5) * P25.J()[0], float(p5) * P25.J()[1]))
ok = abs(C5 - 1 / tau) < mp.mpf(10)**-35 and abs(S5 - 1 / tau) < mp.mpf(10)**-35
ok &= closez(eJ, (float(1 / tau), float(1 / tau) / P25.L), 1e-8)
ok &= abs(abs(1 + mp.expj(q5)) - tau) < mp.mpf(10)**-35
ok &= abs(1 / C5 - tau) < mp.mpf(10)**-35 and abs(1 / S5 - tau) < mp.mpf(10)**-35 and abs(1 / tau - (tau - 1)) < mp.mpf(10)**-35
ok &= abs(sum(mp.expj(q5)**k for k in range(5))) < mp.mpf(10)**-35      # the sides 1, u, ..., u^4 close a regular pentagon
check("num", "rem:golden", "sin(pi/5, pi/5) = cos(pi/5, pi/5) = 1/tau = tau - 1; csc = sec = tau; e^{(pi/5)J} = (1 + u)/tau at 2pi/5; "
      "|1 + u| = tau, a diagonal of the regular pentagon with the sides 1, u, ..., u^4", ok)
mp.mp.dps = 30


def twotrig_num(P):
    z = rnum()
    m = P.mag(z)
    alpha = cmath.phase(P.page(z))
    n = rng.randint(-6, 6)
    zn = P.power(z, n)
    ok_ = close(zn[0], m**n * P.C_(n * alpha), 1e-8) and close(zn[1] * P.L, m**n * P.S_(n * alpha), 1e-8)
    ok_ &= close(P.Sh(zn), m**n * math.cos(n * alpha), 1e-8) and close(P.Ht(zn), m**n * math.sin(n * alpha), 1e-8)
    zz_ = (rng.uniform(-1.5, 1.5), rng.uniform(-1.5, 1.5))
    ez = P.exp(zz_)
    sh, ht = P.Sh(zz_), P.Ht(zz_)
    ok_ &= closez(ez, (math.exp(sh) * P.C_(ht), math.exp(sh) * P.S_(ht) / P.L), 1e-8)
    return ok_


num_check("cor:two-trig", "the parts and readings of z^n; e^z = e^{Sh}(C_{Ht} + S_{Ht} u)", twotrig_num, PLANES, trials=3)


def roots_num(P):
    z = rnum()
    m = P.mag(z)
    alpha = cmath.phase(P.page(z))
    n = rng.randint(1, 8)
    ws = [(m**(1 / n) * P.unit((alpha + 2 * math.pi * k) / n)[0], m**(1 / n) * P.unit((alpha + 2 * math.pi * k) / n)[1])
          for k in range(n)]
    ok_ = all(closez(P.power(w, n), z, 1e-8) for w in ws)
    ok_ &= all(abs(P.page(ws[i_]) - P.page(ws[j_])) > 1e-6 for i_ in range(n) for j_ in range(i_))
    return ok_


num_check("thm:roots", "the n roots w_k satisfy w_k^n = z and are distinct, n = 1..8", roots_num)
check("sym", "cor:sqrt", "(1 + u)^2 = 4 cos^2(p/2) u, and (sqrt|G|/(2cos(p/2)))^2 (1 + u)^2 = G  (p = 2h)",
      trig_zero((1 + sp.cos(2 * hh) + I * sp.sin(2 * hh))**2 - 4 * sp.cos(hh)**2 * (sp.cos(2 * hh) + I * sp.sin(2 * hh)), [hh])
      and trig_zero(L / (4 * sp.cos(hh)**2) * (1 + sp.cos(2 * hh) + I * sp.sin(2 * hh))**2 - L * (sp.cos(2 * hh) + I * sp.sin(2 * hh)), [hh]))
check("num", "cor:sqrt", "the first root lies on the ray alpha = p/2",
      all(close(cmath.phase(1 + P.u), P.p / 2) for P in PLANES))

# =============================================================================
heading("Section 6. The geometry of one magnitude")
# =============================================================================


def circle_num(P):
    alpha = rng.uniform(-math.pi, math.pi)
    z = P.coords(rng.uniform(0.3, 3) * cmath.exp(1j * alpha))
    inside = ((P.p - alpha) % (2 * math.pi)) < math.pi
    return abs(z[0]) < 1e-12 or (z[0] > 0) == inside


num_check("prop:circle", "Re(z) > 0 exactly on the half of the circle that holds the positive real axis", circle_num, PLANES, trials=40)
aa = M * sp.sin(pp - al) / sp.sin(pp)
BB = M * sp.sin(al) / sp.sin(pp)
check("sym", "thm:ellipse", "the parametric point satisfies a^2 + B^2 + 2 cos p aB = M^2",
      trig_zero(aa**2 + BB**2 + 2 * sp.cos(pp) * aa * BB - M**2, [pp, al]))
check("sym", "thm:ellipse", "a + B cos p = M cos(alpha) and B sin p = M sin(alpha)",
      trig_zero(aa + BB * sp.cos(pp) - M * sp.cos(al), [pp, al]) and trig_zero(BB * sp.sin(pp) - M * sp.sin(al), [pp, al]))
Q = sp.Matrix([[1, C], [C, 1]])
check("sym", "thm:ellipse", "the axes: eigenvalues 1 +- cos p on the diagonals (1, +-1)",
      sp.simplify(Q * sp.Matrix([1, 1]) - (1 + C) * sp.Matrix([1, 1])) == sp.zeros(2, 1)
      and sp.simplify(Q * sp.Matrix([1, -1]) - (1 - C) * sp.Matrix([1, -1])) == sp.zeros(2, 1))
check("sym", "thm:ellipse", "the ratio of the semi-axes is tan(p/2)  (p = 2h)",
      trig_zero((1 - sp.cos(2 * hh)) / (1 + sp.cos(2 * hh)) - sp.sin(hh)**2 / sp.cos(hh)**2, [hh]))


def ellipse_num(P):
    c_ = P.c
    s1, s2 = 1 / math.sqrt(1 + c_), 1 / math.sqrt(1 - c_)
    ok_ = close(math.pi * s1 * s2, math.pi / P.s)
    short, long_ = min(s1, s2), max(s1, s2)
    ok_ &= close(1 - (short / long_)**2, 2 * abs(c_) / (1 + abs(c_)))
    # one turn: the winding number of (a(alpha), B(alpha)) about 0 is 1
    pts = [(math.sin(P.p - al_) / P.s, math.sin(al_) / P.s) for al_ in [2 * math.pi * k / 720 for k in range(721)]]
    wind = sum(math.atan2(x1 * y2 - y1 * x2, x1 * x2 + y1 * y2) for (x1, y1), (x2, y2) in zip(pts, pts[1:]))
    ok_ &= close(wind, 2 * math.pi)
    if P.p < math.pi / 2 - 1e-9:
        for (aa_, BB_) in pts[1:-1]:
            if abs(aa_ * BB_) > 1e-9:
                ok_ &= (aa_ * aa_ + BB_ * BB_ < 1) == (aa_ * BB_ > 0)
    return ok_


num_check("thm:ellipse", "area pi M^2/sin p, eccentricity 2|c|/(1+|c|) at every angle, one turn, inside/outside", ellipse_num, PLANES, trials=1)
check("sym", "thm:ellipse", "at p = 0 the curve is the pair of lines a + B = +-M",
      sp.expand((a + b)**2 - (a**2 + b**2 + 2 * 1 * a * b)) == 0)
num_check("def:rings-rays", "a number of magnitude M and angle alpha is the point (M C_alpha, M S_alpha) of the chart",
          lambda P: (lambda m, al_: closez((P.coords(m * cmath.exp(1j * al_))[0], P.coords(m * cmath.exp(1j * al_))[1] * P.L),
                                           (m * P.C_(al_), m * P.S_(al_))))(rng.uniform(0.2, 3), rng.uniform(-3, 3)), PLANES)


def frame_readings(P):
    """Remark 6.4: the axes are rays, the uneven spread of the rays, the axes of the rings."""
    def chart(M, al_):
        return (M * math.sin(P.p - al_) / P.s, M * math.sin(al_) / P.s)
    def ang(al_):
        x_, y_ = chart(1.0, al_)
        return math.atan2(y_, x_)
    ok_ = closez(chart(1.0, 0.0), (1.0, 0.0)) and closez(chart(1.0, P.p), (0.0, 1.0))
    grid_ = [P.p * k / 50 for k in range(1, 50)]                       # 0 < alpha < pbar: first quadrant
    ok_ &= all(0 < ang(g_) < math.pi / 2 for g_ in grid_)
    ok_ &= all(ang(g1) < ang(g2) for g1, g2 in zip(grid_, grid_[1:]))
    grid2 = [P.p + (math.pi - P.p) * k / 50 for k in range(1, 50)]     # pbar < alpha < pi: second quadrant
    ok_ &= all(math.pi / 2 < ang(g_) < math.pi for g_ in grid2)
    M = rng.uniform(0.3, 3)
    xr, _ = chart(M, 0.0)
    _, yc = chart(M, P.p)
    ok_ &= close(abs(xr), M) and close(abs(yc), M)                     # M = |a| and M = |B| on the axes
    for al_ in (P.p / 2, P.p / 2 + math.pi):                           # B = a: the axis 1/sqrt(1 + c)
        x_, y_ = chart(M, al_)
        ok_ &= close(x_, y_) and close(math.hypot(x_, y_), M / math.sqrt(1 + P.c))
    for al_ in (P.p / 2 + math.pi / 2, P.p / 2 - math.pi / 2):        # B = -a: the axis 1/sqrt(1 - c)
        x_, y_ = chart(M, al_)
        ok_ &= close(x_, -y_) and close(math.hypot(x_, y_), M / math.sqrt(1 - P.c))
    return ok_


num_check("rem:frame-readings", "axes are the rays 0 and pbar; quadrants hold pbar and pi - pbar of alpha; ring axes on pbar/2 (+pi/2)",
          frame_readings, PLANES, trials=3)


def near_num(P):
    a1 = rng.uniform(0.1, 2)
    bsign = rng.choice([-1, 1])
    b1 = bsign * a1 * rng.uniform(1.01, 50) / P.L
    eta = abs(b1) * P.L / a1
    Z = P.page((a1, b1))
    axis = P.u if b1 > 0 else -P.u
    delta = abs(cmath.phase(Z / axis))
    formula = math.atan(P.s / (eta + P.c)) if b1 > 0 else math.atan(P.s / (eta - P.c))
    ok_ = close(delta, formula)
    if b1 > 0:
        ok_ &= close(delta, P.p - cmath.phase(Z))
    big = P.page((1.0, 1e4 / P.L))
    ok_ &= abs(abs(cmath.phase(big / P.u)) - P.s / 1e4) < 1e-8
    return ok_


num_check("prop:near-axis", "delta = arctan(sin p/(|eta| +- cos p)); delta = phi for b > 0; delta ~ sin p/|eta|", near_num, PLANES, trials=20)

# =============================================================================
heading("Appendix A. The formulary (rows not stated in the body)")
# =============================================================================


def change_num(P):
    Qp = PLANES[rng.randrange(len(PLANES))]
    z = rnum()
    a2 = P.Sh(z) - P.Ht(z) * Qp.c / Qp.s
    b2 = P.Ht(z) / (Qp.L * Qp.s)
    return close(Qp.page((a2, b2)), P.page(z))


num_check("app:units", "to another true relative plane", change_num, PLANES)
check("sym", "app:plane", "the opposite angle: C_{-t} = C_t + 2 cos p S_t, S_{-t} = -S_t",
      trig_zero(Cth(-th) - Cth(th) - 2 * sp.cos(pp) * Sth(th), [th, pp]) and trig_zero(Sth(-th) + Sth(th), [th, pp]))
check("sym", "app:plane", "reciprocity: C_{p-t} = S_t and S_{p-t} = C_t",
      trig_zero(Cth(pp - th) - Sth(th), [th, pp]) and trig_zero(Sth(pp - th) - Cth(th), [th, pp]))
check("sym", "app:plane", "the mirror z -> u conj z sends x + y u to y + x u",
      zero_CS(U * (x + y * U).subs(I, -I) - (y + x * U)))


def plane_rows(P):
    z = rnum()
    Z = P.page(z)
    tv = cmath.phase(Z)
    m = abs(Z)
    ok_ = close(math.tan(tv), P.Ht(z) / P.Sh(z), 1e-7)
    ok_ &= closez(P.conj(z), (m * P.C_(-tv), m * P.S_(-tv) / P.L))
    return ok_


num_check("app:plane", "the angle tan t = Ht/Sh; the conjugate in polar form", plane_rows, PLANES)
check("sym", "app:turning", "subtraction: C_{a-b} and S_{a-b}",
      trig_zero(Cth(al - be) - (Cth(al) * Cth(be) + Sth(al) * Sth(be) + 2 * sp.cos(pp) * Cth(al) * Sth(be)), [al, be, pp])
      and trig_zero(Sth(al - be) - (Sth(al) * Cth(be) - Cth(al) * Sth(be)), [al, be, pp]))
check("sym", "app:turning", "double angle: C_{2a} = C^2 - S^2, S_{2a} = 2S(C + cos p S)",
      trig_zero(Cth(2 * al) - (Cth(al)**2 - Sth(al)**2), [al, pp])
      and trig_zero(Sth(2 * al) - 2 * Sth(al) * (Cth(al) + sp.cos(pp) * Sth(al)), [al, pp]))
check("sym", "app:turning", "a turn by p, a quarter turn and a half turn",
      trig_zero(Cth(th + pp) + Sth(th), [th, pp]) and trig_zero(Sth(th + pp) - Cth(th) - 2 * sp.cos(pp) * Sth(th), [th, pp])
      and trig_zero(Cth(th + sp.pi / 2) - (-Cth(th) * sp.cos(pp) / sp.sin(pp) - Sth(th) / sp.sin(pp)), [th, pp])
      and trig_zero(Sth(th + sp.pi / 2) - (Cth(th) / sp.sin(pp) + Sth(th) * sp.cos(pp) / sp.sin(pp)), [th, pp])
      and trig_zero(Cth(th + sp.pi) + Cth(th), [th, pp]) and trig_zero(Sth(th + sp.pi) + Sth(th), [th, pp]))
check("sym", "app:euler", "sin np/sin p for n = 1..4: 1, 2c, 4c^2 - 1, 8c^3 - 4c",
      all(trig_zero(sp.sin(n * pp) / sp.sin(pp) - e_, [pp]) for n, e_ in
          [(1, 1), (2, 2 * sp.cos(pp)), (3, 4 * sp.cos(pp)**2 - 1), (4, 8 * sp.cos(pp)**3 - 4 * sp.cos(pp))]))


def log_num(P):
    z = rnum()
    Z = P.page(z)
    tv = cmath.phase(Z)
    lg = (math.log(abs(Z)) - tv * P.c / P.s, tv / (P.L * P.s))
    return closez(P.exp(lg), z, 1e-8)


num_check("app:euler", "the logarithm: e^{ln z} = z", log_num, PLANES, trials=3)
num_check("app:roots", "the roots of unity",
          lambda P: (lambda n: all(closez(P.power(P.unit(2 * math.pi * k / n), n), (1.0, 0.0), 1e-8) for k in range(n)))(rng.randint(1, 9)))


def quad_num(P):
    Pq, Qq = rnum(), rnum()
    disc = (P.mul(Pq, Pq)[0] - 4 * Qq[0], P.mul(Pq, Pq)[1] - 4 * Qq[1])
    m = P.mag(disc)
    alpha = cmath.phase(P.page(disc))
    r = (math.sqrt(m) * P.unit(alpha / 2)[0], math.sqrt(m) * P.unit(alpha / 2)[1])
    ok_ = True
    for sg in (1, -1):
        v = ((-Pq[0] + sg * r[0]) / 2, (-Pq[1] + sg * r[1]) / 2)
        val_ = P.mul(v, v)
        pv = P.mul(Pq, v)
        ok_ &= closez((val_[0] + pv[0] + Qq[0], val_[1] + pv[1] + Qq[1]), (0, 0), 1e-8)
    return ok_


num_check("app:roots", "the quadratic equation v^2 + P v + Q = 0", quad_num)


def Mz(z):
    return sp.Matrix([[z[0], -L**2 * z[1]], [z[1], z[0] + 2 * z[1] * L * C]])


check("sym", "app:matrix", "the matrix of zw is the product of the matrices",
      all(zero_CS(e_) for e_ in (Mz(mul(z_, w_)) - Mz(z_) * Mz(w_))))
check("sym", "app:matrix", "det = |z|^2 and trace = 2 Sh(z); z and conj z are the roots of x^2 - 2Sh x + |z|^2",
      zero_CS(Mz(z_).det() - mag2(z_)) and zero_CS(Mz(z_).trace() - 2 * Sh(z_))
      and same(add(add(mul(z_, z_), scal(-2 * Sh(z_), z_)), (mag2(z_), 0)), (0, 0))
      and same(add(add(mul(conj(z_), conj(z_)), scal(-2 * Sh(z_), conj(z_))), (mag2(z_), 0)), (0, 0)))
MG = Mz(GU)
check("sym", "app:matrix", "the matrix of G obeys the Rule of the Square",
      all(zero_CS(e_) for e_ in (MG * MG - 2 * L * C * MG + L**2 * sp.eye(2))))
qz = quo(z_, w_)
zwb = mul(z_, conj(w_))
check("sym", "app:parts", "the quotient by its top z conj w",
      zero_CS(qz[0] - zwb[0] / mag2(w_)) and zero_CS(qz[1] - zwb[1] / mag2(w_))
      and zero_CS(zwb[0] - (a * conj(w_)[0] + L**2 * b * d)) and zero_CS(zwb[1] - (Sh(w_) * b - Sh(z_) * d))
      and zero_CS(Sh(zwb) - (Sh(z_) * Sh(w_) + Ht(z_) * Ht(w_))) and zero_CS(Ht(zwb) - (Sh(w_) * Ht(z_) - Sh(z_) * Ht(w_))))
Jz = mul(Jc, z_)
check("sym", "app:parts", "a quarter turn: the parts and readings of J z",
      zero_CS(Jz[0] - (-Ht(z_) - Sh(z_) * C / S)) and zero_CS(Jz[1] - Sh(z_) / (L * S))
      and zero_CS(Sh(Jz) + Ht(z_)) and zero_CS(Ht(Jz) - Sh(z_)))
Gz = mul(GU, z_)
check("sym", "app:parts", "the product by G: Re = -|G|^2 b, Com = Re(conj z) G, Sh, Ht",
      zero_CS(Gz[0] + L**2 * b) and zero_CS(Gz[1] - conj(z_)[0])
      and zero_CS(Sh(Gz) - L * (Sh(z_) * C - Ht(z_) * S)) and zero_CS(Ht(Gz) - L * (Sh(z_) * S + Ht(z_) * C)))


def gpowers_num(P):
    n = rng.randint(-5, 8)
    gn = P.power((0.0, 1.0), n)
    return (close(gn[0], -P.L**n * math.sin((n - 1) * P.p) / P.s, 1e-8) and close(gn[1], P.L**(n - 1) * math.sin(n * P.p) / P.s, 1e-8)
            and close(P.Sh(gn), P.L**n * math.cos(n * P.p), 1e-8) and close(P.Ht(gn), P.L**n * math.sin(n * P.p), 1e-8))


num_check("app:parts", "the parts and readings of G^n", gpowers_num, PLANES)


def criteria_num(P):
    r_ = rng.uniform(-3, 3)
    ok_ = closez(P.conj((r_, 0.0)), (r_, 0.0)) and close(P.Ht((r_, 0.0)), 0)
    Jm = (r_ * P.J()[0], r_ * P.J()[1])                         # a multiple of J
    ok_ &= close(P.Sh(Jm), 0) and closez(P.conj(Jm), (-Jm[0], -Jm[1]))
    return ok_


num_check("app:parts", "Ht = 0 exactly when conj z = z; Sh = 0 exactly when conj z = -z", criteria_num, PLANES)
SIX = {"1": (ONE, (1, 0, 1, 0)),
       "G": (GU, (0, 1, L * C, L * S)),
       "conj G": (conj(GU), (2 * L * C, -1, L * C, -L * S)),
       "G^2": (mul(GU, GU), (-L**2, 2 * L * C, L**2 * (2 * C**2 - 1), L**2 * 2 * S * C)),
       "1/G": (inv(GU), (2 * C / L, -1 / L**2, C / L, -S / L)),
       "J": (Jc, (-C / S, 1 / (L * S), 0, 1))}
ok = all(zero_CS(n_[0] - e_[0]) and zero_CS(n_[1] - e_[1]) and zero_CS(Sh(n_) - e_[2]) and zero_CS(Ht(n_) - e_[3])
         for n_, e_ in SIX.values())
check("sym", "app:parts", "the table of the six numbers that recur", ok)

# =============================================================================
heading("Introduction and Summary: prior work and open questions")
# =============================================================================
pq, qq = -L**2, 2 * L * C                         # i^2 = p + q i with these p and q
check("sym", "sec:intro", "Rule of the Square = generalized complex unit with p = -|G|^2, q = 2|G| cos p (elliptic)",
      zero_CS(G**2 - (pq + qq * G)) and zero_CS(qq**2 + 4 * pq + 4 * L**2 * S**2))
check("sym", "sec:intro", "conjugate = (G -> q - G); |z|^2 = a^2 + qab - pb^2; J = (G - q/2)/sqrt(-p - q^2/4)",
      same(conj(z_), (a + b * qq, -b)) and zero_CS(mag2(z_) - (a**2 + qq * a * b - pq * b**2))
      and zero_CS((G - qq / 2) / (L * S) - page(Jc)) and zero_CS((L * S)**2 - (-pq - qq**2 / 4)))
ok = True
for _ in range(2000):
    qv = rng.uniform(-5, 5)
    pv = -qv * qv / 4 - rng.uniform(1e-3, 5)          # an elliptic pair: q^2 + 4p < 0
    Lv = math.sqrt(-pv)
    ang = math.acos(qv / (2 * Lv))
    ok &= 0 < ang < math.pi and close(-Lv**2, pv) and close(2 * Lv * math.cos(ang), qv)
check("num", "sec:intro", "every elliptic pair (p, q) is a leaning frame: length sqrt(-p), angle arccos(q/(2 sqrt(-p)))", ok)
Zp = a + b * G
ok = True
for n in (2, 3, 4):
    F = sp.expand(Zp**n)
    Fr, Fi = split_ri(F)
    V_ = Fi / (L * S)
    U_ = Fr - V_ * L * C
    ok &= zero_CS(sp.diff(U_, b) + L**2 * sp.diff(V_, a))
    ok &= zero_CS(sp.diff(V_, b) - sp.diff(U_, a) - 2 * L * C * sp.diff(V_, a))
check("sym", "sec:summary", "Cauchy-Riemann equations of a leaning frame, for f = Z^2, Z^3, Z^4", ok)

# =============================================================================
heading("Printed numbers (Sections 1, 3, 4, 6 and 7, Figures 1, 4 and 5, Table 2)")
# =============================================================================
mp.mp.dps = 50
pm = mp.mpf("79.86") * mp.pi / 180
Lm = mp.cos(pm) + mp.sin(pm)
cm_, sm_ = mp.cos(pm), mp.sin(pm)


def mmul(z, w):
    (a1, b1), (a2, b2) = z, w
    return (a1 * a2 - Lm**2 * b1 * b2, a1 * b2 + b1 * a2 + 2 * Lm * cm_ * b1 * b2)


def mconj(z):
    return (z[0] + 2 * z[1] * Lm * cm_, -z[1])


def mSh(z):
    return z[0] + z[1] * Lm * cm_


def mHt(z):
    return z[1] * Lm * sm_


def mmag(z):
    return mp.sqrt(mSh(z)**2 + mHt(z)**2)


def minv(z):
    m2 = mmag(z)**2
    zc = mconj(z)
    return (zc[0] / m2, zc[1] / m2)


def fmt(v, printed):
    """v rounded to as many decimals as the printed value has."""
    n = len(printed.split(".")[1]) if "." in printed else 0
    s_ = f"{float(mp.mpf(v)):.{n}f}"
    return s_[1:] if s_.startswith("-") and float(s_) == 0 else s_


def val(name, v, printed):
    got = fmt(v, printed)
    ok_ = got == printed
    check("val", "printed", f"{name}: {printed}" + ("" if ok_ else f"  (computed {got})"), ok_)


z1 = (mp.mpf(1), mp.mpf(1))
Jm = (-cm_ / sm_, 1 / (Lm * sm_))
r1 = mp.asin(Lm / mp.sqrt(2)) - mp.pi / 4
val("length of the Composer Unit at 79.86 deg", Lm, "1.16043")
val("Table 2: the same length to seven decimals", Lm, "1.1604345")
val("cos 79.86 deg", cm_, "0.17605")
val("sin 79.86 deg", sm_, "0.98438")
val("arcsine root of that length, degrees", r1 * 180 / mp.pi, "10.14")
val("mirror root of that length, degrees", 90 - r1 * 180 / mp.pi, "79.86")
val("arcsine root of the printed length 1.16043, degrees",
    (mp.asin(mp.mpf("1.16043") / mp.sqrt(2)) - mp.pi / 4) * 180 / mp.pi, "10.14")
G2 = mmul((0, 1), (0, 1))
G3 = mmul(G2, (0, 1))
val("G^2, real part", G2[0], "-1.34661")
val("G^2, Composer coordinate", G2[1], "0.40860")
val("G^3, real part", G3[0], "-0.55022")
val("G^3, Composer coordinate", G3[1], "-1.17966")
val("J, real part (= -cot p)", Jm[0], "-0.17885")
val("J, Composer coordinate", Jm[1], "0.87542")
k_sq = mp.sqrt(Lm) / (2 * mp.cos(pm / 2))
val("square root of G, real part", k_sq, "0.70240")
val("square root of G, Composer coordinate", k_sq / Lm, "0.60529")
sqG = (k_sq, k_sq / Lm)
check("val", "printed", "the square of that root is G", abs(mmul(sqG, sqG)[0]) < 1e-40 and abs(mmul(sqG, sqG)[1] - 1) < 1e-40)
val("Sh(1 + G)", mSh(z1), "1.20430")
val("Ht(1 + G)", mHt(z1), "1.14231")
val("|1 + G|", mmag(z1), "1.65988")
val("|1 + G|^2", mmag(z1)**2, "2.75521")
alpha1 = mp.atan2(mHt(z1), mSh(z1))
val("studied angle of 1 + G, degrees", alpha1 * 180 / mp.pi, "43.49")
val("projection angle of 1 + G, degrees", (pm - alpha1) * 180 / mp.pi, "36.37")
zc1 = mconj(z1)
val("conjugate, real part", zc1[0], "1.40860")
check("val", "printed", "conjugate, Composer coordinate: -1", zc1[1] == -1)
val("Section 7: 2(1.20430) - 1", 2 * mp.mpf("1.20430") - 1, "1.40860")
z2m = mmul(z1, z1)
val("z^2, real part", z2m[0], "-0.34661")
val("z^2, Composer coordinate (= 2 Sh(z))", z2m[1], "2.40860")
val("2 Sh(z)", 2 * mSh(z1), "2.40860")
z3m = mmul(z2m, z1)
val("z^3, real part", z3m[0], "-3.59005")
val("z^3, Composer coordinate", z3m[1], "3.04614")
zi = minv(z1)
val("1/z, real part", zi[0], "0.51125")
val("1/z, Composer coordinate", zi[1], "-0.36295")
Jz1 = mmul(Jm, z1)
val("J z, real part", Jz1[0], "-1.35770")
val("J z, Composer coordinate", Jz1[1], "1.05427")
half = alpha1 / 2
rt = (mp.sqrt(mmag(z1)) * mp.sin(pm - half) / sm_, mp.sqrt(mmag(z1)) * mp.sin(half) / sm_ / Lm)
val("square root of z, real part", rt[0], "1.11134")
val("square root of z, Composer coordinate", rt[1], "0.41782")
val("square root of z, angle in degrees", half * 180 / mp.pi, "21.74")
check("val", "printed", "that root squared is z", max(abs(mmul(rt, rt)[0] - 1), abs(mmul(rt, rt)[1] - 1)) < 1e-40)
sk = mmul(z1, (1, -1))
val("(1 + G)(1 - G), real part", sk[0], "2.34661")
val("(1 + G)(1 - G), Composer coordinate", sk[1], "-0.40860")
val("signed area of 1 + G and 2 - G", mHt(mmul(mconj(z1), (2, -1))), "-3.42693")
val("semi-axis on B = a, per unit of M", 1 / mp.sqrt(1 + cm_), "0.92212")
val("semi-axis on B = -a, per unit of M", 1 / mp.sqrt(1 - cm_), "1.10167")
val("their ratio", mp.sqrt((1 - cm_) / (1 + cm_)), "0.83702")
val("their ratio = tan(p/2)", mp.tan(pm / 2), "0.83702")
p3 = mp.pi / 3
L3 = mp.cos(p3) + mp.sin(p3)
m3 = mp.sqrt(1 + 2 * L3 * mp.cos(p3) + L3**2)
a3 = mp.atan2(L3 * mp.sin(p3), 1 + L3 * mp.cos(p3))
val("Remark 5.7: |G| at pi/5", mp.cos(mp.pi / 5) + mp.sin(mp.pi / 5), "1.39680")
val("Remark 5.7: G^10 at pi/5", (mp.cos(mp.pi / 5) + mp.sin(mp.pi / 5))**10, "28.27153")
val("Figure 6(a): |G| at pi/3", L3, "1.36603")
val("Figure 6(a): G^6 at pi/3", L3**6, "6.49760")
val("Figure 6(f): G^12 at pi/6", (mp.cos(mp.pi / 6) + mp.sin(mp.pi / 6))**12, "42.21874")
val("Figure 6(d): the stretch per step at 79.86 deg, per cent", (Lm - 1) * 100, "16")
val("Figure 8(b): the ring of 1 + G at pi/3", m3, "2.05719")
val("Figure 8(b): the ray of 1 + G at pi/3, degrees", a3 * 180 / mp.pi, "35.10")
val("a part longer than its number, at most (per cent)", (1 / sm_ - 1) * 100, "1.6")
for nm, ang, pr in [("pi/6", mp.pi / 6, "1.366"), ("pi/4", mp.pi / 4, "1.414"), ("pi/3", mp.pi / 3, "1.366"),
                    ("2pi/3", 2 * mp.pi / 3, "0.366"), ("5pi/6", 5 * mp.pi / 6, "0.366"), ("pi/2", mp.pi / 2, "1")]:
    val(f"Figures 1 and 4: length of the Composer Unit at {nm}", abs(mp.cos(ang) + mp.sin(ang)), pr)
K_coul = mp.mpf("8.9875517862e9")       # CODATA 2022, 1/(4 pi eps0) in N m^2 C^-2
G_grav = mp.mpf("6.67430e-11")          # CODATA 2022, N m^2 kg^-2
val("Section 1.6: sqrt(K/G) in units of 10^10 kg/C", mp.sqrt(K_coul / G_grav) / mp.mpf(10)**10, "1.16043")
check("val", "printed", "Section 1.1: (2 +- i)^3 = 2 +- 11i; 4^3 = 15*4 + 4; 5 +- sqrt(-15) has sum 10 and product 40",
      (2 + 1j)**3 == (2 + 11j) and (2 - 1j)**3 == (2 - 11j) and 4**3 == 15 * 4 + 4
      and (5 + cmath.sqrt(-15)) + (5 - cmath.sqrt(-15)) == 10 and abs((5 + cmath.sqrt(-15)) * (5 - cmath.sqrt(-15)) - 40) < 1e-12)

# =============================================================================
heading("Coverage of the numbered statements")
# =============================================================================
covered = set()
for _, lab_, _, _ in RESULTS:
    for part in re.split(r"[,\s]+", lab_):
        if part:
            covered.add(part)
tex = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gn-paper.tex")
statements = []
if os.path.exists(tex):
    src = open(tex, encoding="utf-8").read()
    statements = re.findall(r"\\begin\{(?:theorem|proposition|corollary|lemma|axiom|definition|remark)\}(?:\[(?:[^\[\]]|\[[^\]]*\])*\])?\s*\\label\{([^}]+)\}", src)
    missing = [s_ for s_ in statements if s_ not in covered]
    check("cov", "coverage", f"all {len(statements)} numbered statements of the paper are checked"
          + ("" if not missing else f"; missing: {missing}"), not missing)
else:
    print("  (gn-paper.tex not found: coverage not checked)")

# =============================================================================
heading("Summary")
# =============================================================================
tags = {}
for tag, _, _, ok_ in RESULTS:
    tags.setdefault(tag, [0, 0])
    tags[tag][0 if ok_ else 1] += 1
for tag in ("sym", "num", "ctl", "val", "cov"):
    if tag in tags:
        print(f"  {tag}: {tags[tag][0]} passed, {tags[tag][1]} failed")
npass = sum(1 for r_ in RESULTS if r_[3])
nfail = len(RESULTS) - npass
print(f"\n  {npass} of {len(RESULTS)} checks passed, {nfail} failed.")
print(f"  {len(statements)} numbered statements covered; {EVALS[0]} numerical comparisons in "
      f"{NUM_PLANES} relative planes ({len(ANGLES)} angles; the Axiom's length and two other lengths).")
print(f"  Run time {time.time() - T0:.1f} s.")
sys.exit(0 if nfail == 0 else 1)
