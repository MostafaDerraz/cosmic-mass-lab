#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_general_number_laws.py  (version 2, 6 November 2026)

A digital test of the algebra of General Numbers: every law of Section 4.9, with its last
part, Subsection 4.9.7 (the real part and the Composer part beside the shadow and the
height), and every law of Appendix F, with its Section F.8, in "The Cosmic Mass", Volume 1.
It also answers one question on its own: are the complex numbers still inside, as the
special case phibar = pi/2?

HOW IT TESTS
Each law is computed in two independent ways, and the two results are compared.

  1. In the leaning frame, as the book computes. A General Number z = a + bG is the
     pair (a, b), and pairs are multiplied with the Rule of the Square,
     G^2 = 2|G| cos(phibar) G - |G|^2 (Propositions 4.22 and 4.23).
  2. On the page, as arrows. z is the point a + b|G| (cos phibar + i sin phibar),
     and arrows are multiplied by Wessel's rule (lengths multiply, angles add),
     which is the multiplication of complex numbers that Python does by itself.
     The real part and the Composer coordinate of an arrow are read back from the page by
     solving two linear equations; its shadow and height are its two coordinates.

The values of the General Sinus always come from the definitions of Appendix A,
sin(alpha, phi) = sin(alpha)/sin(alpha + phi) and cos(alpha, phi) = sin(phi)/sin(alpha + phi),
never from the law under test. Further independent checks: the exponential from its
power series, derivatives from finite differences, and all roots of a polynomial from an
independent root finder (Durand-Kerner).

Every law is tried on many random numbers in many relative planes:
  - the book's plane, phibar = 79.86 degrees, with |G| = |cos phibar + sin phibar| = 1.16043;
  - the plane phibar = pi/3, where cos, sin and |G| are exact numbers;
  - the complex plane, phibar = pi/2, where the Axiom gives |G| = 1;
  - two planes next to it, 89.9 and 90.1 degrees;
  - 60 random angles, with the length that the Axiom gives;
  - 20 random angles with a free length (the laws hold for any length of the Composer Unit).
Negative controls check that the test catches wrong laws: each wrong law must fail.

THE COMPLEX CASE
  - Every law, old and new, is evaluated at phibar = pi/2 and compared with its law in
    Appendix C, computed with Python's own complex numbers.
  - Exact checks, with no rounding at all: at pi/2 and at pi/3 the numbers are written
    p + q*sqrt(3) with fractions p and q, and the laws are checked exactly. At pi/2 the
    arithmetic of the leaning frame must be exactly that of the complex numbers.
  - The wrong laws of the controls are the habits of the complex plane. Each must fail
    where the Composer axis leans and hold where it stands upright.
  - A table follows the habits as the angle approaches 90 degrees.

Only the Python standard library is needed (Python 3.8 or later). Run

    python3 test_general_number_laws.py

The script prints a summary and writes two reports next to itself:
test-report.html (open it in a browser; it prints well to PDF) and test-report.txt.
"""
import cmath
import datetime
import html
import math
import os
import platform
import random
import sys
import time
from fractions import Fraction

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

SEED = 20261106          # the random numbers are the same at every run
SAMPLES = 100            # random cases per law and per plane
TOL = 1e-9               # largest relative error accepted for an exact law
TOL_DIFF = 1e-7          # for laws checked with finite differences
NEGATIVE_MIN = 1e-3      # a wrong law must miss by at least this much

# symbols for the reports
GU = "\U0001D4A2"        # the Composer Unit (script G)
JU = "\U0001D4A5"        # the Turning Unit (script J)
PB = "φ̄"      # phi with a bar


# =====================================================================================
# 1. Relative planes and General Numbers
# =====================================================================================
class Plane:
    """A relative plane: the angle phibar of the Composer axis and the length |G|."""

    def __init__(self, phibar, length=None, label="", kind="axiom"):
        self.p = phibar
        self.c, self.s = math.cos(phibar), math.sin(phibar)
        self.L = abs(self.c + self.s) if length is None else length   # Axiom 4.1
        self.U = complex(self.c, self.s)          # the direction of G on the page
        self.Gpage = self.L * self.U              # the Composer Unit as an arrow
        self.label = label
        self.kind = kind

    def deg(self):
        return math.degrees(self.p)


class GN:
    """A General Number a + bG of one relative plane, computed in the leaning frame."""

    __slots__ = ("a", "b", "P")

    def __init__(self, a, b, P):
        self.a, self.b, self.P = a, b, P

    def lift(self, o):
        return o if isinstance(o, GN) else GN(float(o), 0.0, self.P)

    def __add__(self, o):
        o = self.lift(o)
        return GN(self.a + o.a, self.b + o.b, self.P)

    __radd__ = __add__

    def __neg__(self):
        return GN(-self.a, -self.b, self.P)

    def __sub__(self, o):
        o = self.lift(o)
        return GN(self.a - o.a, self.b - o.b, self.P)

    def __rsub__(self, o):
        return self.lift(o) - self

    def __mul__(self, o):
        if not isinstance(o, GN):                       # a real number only scales
            return GN(self.a * o, self.b * o, self.P)
        P = self.P                                      # Proposition 4.23
        return GN(self.a * o.a - P.L ** 2 * self.b * o.b,
                  self.a * o.b + self.b * o.a + 2 * P.L * P.c * self.b * o.b, P)

    __rmul__ = __mul__

    def conj(self):                                     # Definition 4.37
        P = self.P
        return GN(self.a + 2 * self.b * P.L * P.c, -self.b, P)

    def norm2(self):                                    # Theorem 4.9: |z|^2
        P = self.P
        return self.a ** 2 + 2 * self.a * self.b * P.L * P.c + self.b ** 2 * P.L ** 2

    def inv(self):                                      # Theorem 4.41
        n2 = self.norm2()
        zc = self.conj()
        return GN(zc.a / n2, zc.b / n2, self.P)

    def shadow(self):                                   # Definition 4.19
        return self.a + self.b * self.P.L * self.P.c

    def height(self):                                   # Proposition 4.44
        return self.b * self.P.L * self.P.s

    def page(self):
        """The same number as an arrow of the page: a steps of the arrow 1 and b steps of the arrow G."""
        return self.a * 1 + self.b * self.P.Gpage

    def power(self, n):
        out = GN(1.0, 0.0, self.P)
        base = self if n >= 0 else self.inv()
        for _ in range(abs(n)):
            out = out * base
        return out


def from_page(w, P):
    """The pair (a, b) with a + b G = w: two linear equations, solved by Cramer's rule."""
    gx, gy = P.Gpage.real, P.Gpage.imag
    det = 1.0 * gy - 0.0 * gx
    return GN((w.real * gy - w.imag * gx) / det, (1.0 * w.imag - 0.0 * w.real) / det, P)


def parts_page(w, P):
    """Read an arrow of the page in both frames: its real part, Composer coordinate, shadow and height."""
    z = from_page(w, P)
    return z.a, z.b, w.real, w.imag


# The General Sinus, from its definitions (Appendix A)
def gsin(al, ph):
    return math.sin(al) / math.sin(al + ph)


def gcos(al, ph):
    return math.sin(ph) / math.sin(al + ph)


def C(t, P):
    """cos(t, phibar - t): the general cosine at the studied angle t of the plane."""
    return gcos(t, P.p - t)


def S(t, P):
    """sin(t, phibar - t): the general sine at the studied angle t of the plane."""
    return gsin(t, P.p - t)


def unit(t, P):
    """The General Number C_t + S_t u, with u = G/|G| (the claimed unit number at the angle t)."""
    return GN(C(t, P), S(t, P) / P.L, P)


def G_of(P):
    return GN(0.0, 1.0, P)


def u_of(P):
    return GN(0.0, 1.0 / P.L, P)


def J_of(P):
    """The Turning Unit of Definition 4.42: -cot(phibar) + G/(|G| sin(phibar))."""
    return GN(-P.c / P.s, 1.0 / (P.L * P.s), P)


def gn_exp(w):
    """e^w from the power series, summed in the leaning frame (scaled down, then squared back)."""
    m, size = 0, abs(w.page())
    while size > 0.25:
        size /= 2
        m += 1
    x = w * (1.0 / 2 ** m)
    term = GN(1.0, 0.0, w.P)
    total = GN(1.0, 0.0, w.P)
    for k in range(1, 30):
        term = (term * x) * (1.0 / k)
        total = total + term
    for _ in range(m):
        total = total * total
    return total


def poly_roots(coeffs):
    """All roots of w^n + c[n-1] w^(n-1) + ... + c[0], by the Durand-Kerner iteration."""
    n = len(coeffs)

    def p(w):
        acc = 1 + 0j
        for ci in reversed(coeffs):
            acc = acc * w + ci
        return acc

    roots = [(0.4 + 0.9j) ** k for k in range(n)]
    for _ in range(5000):
        new = []
        for i, r in enumerate(roots):
            den = 1 + 0j
            for j, s_ in enumerate(roots):
                if i != j:
                    den *= (r - s_)
            new.append(r - p(r) / den)
        done = max(abs(x - y) for x, y in zip(new, roots)) < 1e-15
        roots = new
        if done:
            break
    return roots


# Errors
def err_num(x, y):
    return abs(x - y) / max(1.0, abs(y))


def err_gn(z, w):
    w = z.lift(w) if not isinstance(w, GN) else w
    return max(abs(z.a - w.a), abs(z.b - w.b)) / max(1.0, abs(w.a), abs(w.b))


def rgn(P, rng, lo=-4.0, hi=4.0):
    """A random General Number, kept away from 0."""
    while True:
        z = GN(rng.uniform(lo, hi), rng.uniform(lo, hi), P)
        if abs(z.page()) > 0.05:
            return z


def rangle(rng, lo=-math.pi, hi=math.pi):
    return rng.uniform(lo, hi)


# =====================================================================================
# 2. The laws
# =====================================================================================
GROUPS = []          # (key, title, note)
LAWS = []            # registered tests


def group(key, title, note=""):
    GROUPS.append((key, title, note))


def law(grp, ref, name, formula, method, planes="all", samples=SAMPLES, tol=TOL, negative=False):
    def deco(fn):
        LAWS.append(dict(group=grp, ref=ref, name=name, formula=formula, method=method,
                         planes=planes, samples=samples, tol=tol, negative=negative, fn=fn))
        return fn
    return deco


FRAME_PAGE = "frame against page"

# ---------------------------------------------------------------- 2.1 foundations
group("found", "Foundations: the leaning frame against the page",
      "The frame arithmetic of the book is checked against arrows multiplied on the page.")


@law("found", "Prop. 4.22", "The Rule of the Square",
     f"{GU}² = 2|{GU}| cos {PB} {GU} − |{GU}|²", FRAME_PAGE, samples=1)
def t_square(P, rng):
    G = G_of(P)
    return max(err_num((G * G).page(), P.Gpage ** 2),
               err_gn(G * G, GN(-P.L ** 2, 2 * P.L * P.c, P)))


@law("found", "Def. 4.17", "Sum, part by part", f"(a + b{GU}) + (c + d{GU}) = (a + c) + (b + d){GU}", FRAME_PAGE)
def t_sum(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    return err_num((z1 + z2).page(), z1.page() + z2.page())


@law("found", "Prop. 4.23", "Product of two General Numbers",
     f"(ac − |{GU}|²bd) + (ad + bc + 2|{GU}| cos {PB} bd){GU}", FRAME_PAGE)
def t_product(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    return err_num((z1 * z2).page(), z1.page() * z2.page())


@law("found", "Def. 4.19", "Shadow", f"Sh(z) = a + b|{GU}| cos {PB}", "frame against page")
def t_shadow(P, rng):
    z = rgn(P, rng)
    return err_num(z.shadow(), z.page().real)


@law("found", "Thm. 4.9", "Magnitude", f"|z|² = a² + 2ab|{GU}| cos {PB} + b²|{GU}|²", FRAME_PAGE)
def t_magnitude(P, rng):
    z = rgn(P, rng)
    return err_num(z.norm2(), abs(z.page()) ** 2)


# ---------------------------------------------------------------- 2.2 field
group("field", "The rules of arithmetic (Theorem 4.36)",
      "Every true relative plane is a field: each rule is checked in the leaning frame.")


@law("field", "Thm. 4.36 (2)", "Order of a product", "z₁z₂ = z₂z₁", "leaning frame")
def t_comm(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    return err_gn(z1 * z2, z2 * z1)


@law("field", "Thm. 4.36 (2)", "Grouping of a product", "(z₁z₂)z₃ = z₁(z₂z₃)", "leaning frame")
def t_assoc(P, rng):
    z1, z2, z3 = rgn(P, rng), rgn(P, rng), rgn(P, rng)
    return err_gn((z1 * z2) * z3, z1 * (z2 * z3))


@law("field", "Thm. 4.36 (3)", "Distribution", "z₁(z₂ + z₃) = z₁z₂ + z₁z₃", "leaning frame")
def t_dist(P, rng):
    z1, z2, z3 = rgn(P, rng), rgn(P, rng), rgn(P, rng)
    return err_gn(z1 * (z2 + z3), z1 * z2 + z1 * z3)


@law("field", "Thm. 4.36 (4), (5)", "Zero, one and the opposite", "z + 0 = z, 1·z = z, z + (−z) = 0", "leaning frame")
def t_neutral(P, rng):
    z = rgn(P, rng)
    return max(err_gn(z + 0, z), err_gn(GN(1.0, 0.0, P) * z, z), err_gn(z + (-z), GN(0.0, 0.0, P)))


@law("field", "Thm. 4.36 (6)", "Every number but 0 has an inverse", "z · (1/z) = 1", "leaning frame")
def t_inverse(P, rng):
    z = rgn(P, rng)
    return err_gn(z * z.inv(), GN(1.0, 0.0, P))


# ---------------------------------------------------------------- 2.3 conjugate
group("conj", "The conjugate (Definition 4.37, Theorem 4.38, Remark 4.39)")


@law("conj", "Thm. 4.38 (1)", "The conjugate is the mirror image in the real axis",
     f"z̄ = (a + 2b|{GU}| cos {PB}) − b{GU}", FRAME_PAGE)
def t_mirror(P, rng):
    z = rgn(P, rng)
    return err_num(z.conj().page(), z.page().conjugate())


@law("conj", "Thm. 4.38 (2)", "Conjugate of the conjugate", "z̄̄ = z", "leaning frame")
def t_conj2(P, rng):
    z = rgn(P, rng)
    return err_gn(z.conj().conj(), z)


@law("conj", "Thm. 4.38 (2)", "Conjugate of a sum and of a product",
     "conj(z₁ + z₂) = z̄₁ + z̄₂,  conj(z₁z₂) = z̄₁ z̄₂", "leaning frame")
def t_conj_ops(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    return max(err_gn((z1 + z2).conj(), z1.conj() + z2.conj()), err_gn((z1 * z2).conj(), z1.conj() * z2.conj()))


@law("conj", "Thm. 4.38 (3)", "A number plus its conjugate", "z + z̄ = 2 Sh(z)", "leaning frame")
def t_conj_sum(P, rng):
    z = rgn(P, rng)
    return err_gn(z + z.conj(), GN(2 * z.shadow(), 0.0, P))


@law("conj", "Thm. 4.38 (4)", "A number times its conjugate",
     f"z z̄ = a² + 2ab|{GU}| cos {PB} + b²|{GU}|² = |z|², a real number", "leaning frame")
def t_conj_prod(P, rng):
    z = rgn(P, rng)
    return err_gn(z * z.conj(), GN(z.norm2(), 0.0, P))


@law("conj", "Thm. 4.38 (4)", "The square of the magnitude is positive",
     "|z|² > 0 for z ≠ 0, and |0|² = 0", "sign check")
def t_positive(P, rng):
    z = rgn(P, rng)
    return 0.0 if (z.norm2() > 0 and GN(0.0, 0.0, P).norm2() == 0) else 1.0


@law("conj", "Thm. 4.38 (5)", "The Composer Unit and its conjugate",
     f"{GU} + {GU}̄ = 2|{GU}| cos {PB},  {GU}{GU}̄ = |{GU}|²", "leaning frame", samples=1)
def t_G_conj(P, rng):
    G = G_of(P)
    return max(err_gn(G + G.conj(), GN(2 * P.L * P.c, 0.0, P)), err_gn(G * G.conj(), GN(P.L ** 2, 0.0, P)))


@law("conj", "Thm. 4.38 (5)", "They are the two roots of one equation",
     f"x² − 2|{GU}| cos {PB} x + |{GU}|² = 0 at x = {GU} and x = {GU}̄", "leaning frame", samples=1)
def t_G_roots(P, rng):
    G = G_of(P)
    return max(err_gn(x * x - x * (2 * P.L * P.c) + P.L ** 2, GN(0.0, 0.0, P)) / max(1, P.L ** 2) for x in (G, G.conj()))


@law("conj", "Rem. 4.39", "The skew flip changes the length",
     f"|a − b{GU}|² = |z|² − 4ab|{GU}| cos {PB}", "leaning frame")
def t_skew_len(P, rng):
    z = rgn(P, rng)
    return err_num(GN(z.a, -z.b, P).norm2(), z.norm2() - 4 * z.a * z.b * P.L * P.c)


@law("conj", "Rem. 4.39", "The skew flip gives no real product",
     f"(a + b{GU})(a − b{GU}) = a² + b²|{GU}|² − 2b²|{GU}| cos {PB} {GU}", "leaning frame")
def t_skew_prod(P, rng):
    z = rgn(P, rng)
    return err_gn(z * GN(z.a, -z.b, P), GN(z.a ** 2 + z.b ** 2 * P.L ** 2, -2 * z.b ** 2 * P.L * P.c, P))


# ---------------------------------------------------------------- 2.4 division
group("div", "Magnitude and division (Theorems 4.40 and 4.41)")


@law("div", "Thm. 4.40", "The magnitude of a product", "|z₁z₂| = |z₁| |z₂|", "leaning frame and page")
def t_mag_prod(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    return max(err_num(math.sqrt((z1 * z2).norm2()), math.sqrt(z1.norm2() * z2.norm2())),
               err_num(math.sqrt((z1 * z2).norm2()), abs(z1.page()) * abs(z2.page())))


@law("div", "Thm. 4.40", "The general two-squares identity",
     f"(a₁² + 2a₁b₁|{GU}|cos {PB} + b₁²|{GU}|²)(a₂² + …) = c² + 2cd|{GU}|cos {PB} + d²|{GU}|²", "numbers")
def t_two_squares(P, rng):
    a1, b1, a2, b2 = (rng.uniform(-4, 4) for _ in range(4))
    L, c_ = P.L, P.c
    c = a1 * a2 - L ** 2 * b1 * b2
    d = a1 * b2 + b1 * a2 + 2 * L * c_ * b1 * b2
    lhs = (a1 ** 2 + 2 * a1 * b1 * L * c_ + b1 ** 2 * L ** 2) * (a2 ** 2 + 2 * a2 * b2 * L * c_ + b2 ** 2 * L ** 2)
    return err_num(lhs, c ** 2 + 2 * c * d * L * c_ + d ** 2 * L ** 2)


@law("div", "Thm. 4.41", "The inverse", f"1/z = [(a + 2b|{GU}| cos {PB}) − b{GU}] / |z|²", FRAME_PAGE)
def t_inv_page(P, rng):
    z = rgn(P, rng)
    return err_num(z.inv().page(), 1 / z.page())


@law("div", "Thm. 4.41", "The quotient",
     f"z₁/z₂ = [(a₁a₂ + 2a₁b₂|{GU}|cos {PB} + b₁b₂|{GU}|²) + (b₁a₂ − a₁b₂){GU}] / |z₂|²", FRAME_PAGE)
def t_quotient(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    L, c_ = P.L, P.c
    den = z2.norm2()
    q = GN((z1.a * z2.a + 2 * z1.a * z2.b * L * c_ + z1.b * z2.b * L ** 2) / den, (z1.b * z2.a - z1.a * z2.b) / den, P)
    return err_num(q.page(), z1.page() / z2.page())


@law("div", "Thm. 4.41", f"The inverse of {GU}", f"1/{GU} = (2|{GU}| cos {PB} − {GU}) / |{GU}|²", FRAME_PAGE, samples=1)
def t_inv_G(P, rng):
    return err_gn(GN(2 * P.c / P.L, -1 / P.L ** 2, P), from_page(1 / P.Gpage, P))


# ---------------------------------------------------------------- 2.5 Turning Unit
group("turnunit", "The Turning Unit and the two frames (Definition 4.42, Theorem 4.43, Proposition 4.44)")


@law("turnunit", "Thm. 4.43 (1)", "The Turning Unit stands upright, with length 1",
     f"{JU} = −cot {PB} + {GU}/(|{GU}| sin {PB}) is the arrow i", FRAME_PAGE, samples=1)
def t_J_page(P, rng):
    return err_num(J_of(P).page(), 1j)


@law("turnunit", "Def. 4.42", "The two forms of the Turning Unit agree",
     f"({GU}/|{GU}| − cos {PB}) / sin {PB} = −cot {PB} + {GU}/(|{GU}| sin {PB})", "leaning frame", samples=1)
def t_J_forms(P, rng):
    return err_gn((u_of(P) - P.c) * (1 / P.s), J_of(P))


@law("turnunit", "Thm. 4.43 (2)", "Its square is −1", f"{JU}² = −1", "leaning frame", samples=1)
def t_J_square(P, rng):
    J = J_of(P)
    return err_gn(J * J, GN(-1.0, 0.0, P))


@law("turnunit", "Thm. 4.43 (2)", "No other number squares to −1",
     f"the roots of w² + 1 = 0 are {JU} and −{JU} only", "independent root finder", samples=1)
def t_J_only(P, rng):
    J = J_of(P)
    found = [from_page(r, P) for r in poly_roots([1 + 0j, 0j])]
    best = min(max(err_gn(found[0], J), err_gn(found[1], -J)), max(err_gn(found[0], -J), err_gn(found[1], J)))
    return best


@law("turnunit", "Thm. 4.43 (3)", "Multiplying by the Turning Unit is a quarter turn",
     f"page({JU}z) = i · page(z), the same length", FRAME_PAGE)
def t_J_turn(P, rng):
    z = rgn(P, rng)
    w = J_of(P) * z
    return max(err_num(w.page(), 1j * z.page()), err_num(abs(w.page()), abs(z.page())))


@law("turnunit", "Thm. 4.43 (4)", "At π/2 the two units are one", f"{JU} = {GU} = i", "complex plane", planes="complex", samples=1)
def t_J_complex(P, rng):
    return max(err_gn(J_of(P), G_of(P)), err_num(G_of(P).page(), 1j))


@law("turnunit", "Prop. 4.44", "Shadow and height", f"z = Sh(z) + Ht(z) {JU}", "leaning frame")
def t_sh_ht(P, rng):
    z = rgn(P, rng)
    return err_gn(GN(z.shadow(), 0.0, P) + J_of(P) * z.height(), z)


@law("turnunit", "Prop. 4.44", "Back to the leaning frame",
     f"a = Sh − Ht cot {PB},  b = Ht/(|{GU}| sin {PB})", "page to frame")
def t_back(P, rng):
    w = complex(rng.uniform(-4, 4), rng.uniform(-4, 4))
    z = from_page(w, P)
    return err_gn(GN(w.real - w.imag * P.c / P.s, w.imag / (P.L * P.s), P), z)


@law("turnunit", "Prop. 4.44", "A number minus its conjugate", f"z − z̄ = 2 Ht(z) {JU}", "leaning frame")
def t_minus_conj(P, rng):
    z = rgn(P, rng)
    return err_gn(z - z.conj(), J_of(P) * (2 * z.height()))


@law("turnunit", "Prop. 4.44", "Pythagoras for the shadow and the height", "|z|² = Sh(z)² + Ht(z)²", "numbers")
def t_pyth_sh(P, rng):
    z = rgn(P, rng)
    return err_num(z.norm2(), z.shadow() ** 2 + z.height() ** 2)


# ---------------------------------------------------------------- 2.6 polar form
group("polar", "The polar form with one angle (Theorem 4.45, Theorem 4.10)")


@law("polar", "Thm. 4.45", "Polar form",
     f"z = |z| [cos(α, {PB} − α) + sin(α, {PB} − α) u], α the angle of the arrow", FRAME_PAGE)
def t_polar(P, rng):
    z = rgn(P, rng)
    al = cmath.phase(z.page())
    return err_gn(unit(al, P) * abs(z.page()), z)


@law("polar", "Thm. 4.45", "The parts in polar form",
     f"a = |z| cos(α, {PB} − α),  b|{GU}| = |z| sin(α, {PB} − α)", FRAME_PAGE)
def t_parts(P, rng):
    z = rgn(P, rng)
    al, M = cmath.phase(z.page()), abs(z.page())
    return max(err_num(M * C(al, P), z.a), err_num(M * S(al, P), z.b * P.L))


@law("polar", "Thm. 4.45", "The general Pythagoras identity",
     f"cos(α, {PB}−α)² + 2 cos {PB} cos(α, {PB}−α) sin(α, {PB}−α) + sin(α, {PB}−α)² = 1", "numbers")
def t_gen_pyth(P, rng):
    al = rangle(rng, -10, 10)
    return err_num(C(al, P) ** 2 + 2 * P.c * C(al, P) * S(al, P) + S(al, P) ** 2, 1.0)


@law("polar", "Thm. 4.10", "The phase angle",
     f"tan α = η sin {PB} / (1 + η cos {PB}),  η = |{GU}| b/a", FRAME_PAGE)
def t_phase(P, rng):
    while True:
        z = rgn(P, rng)
        eta = P.L * z.b / z.a if abs(z.a) > 0.05 else None
        if eta is not None and abs(1 + eta * P.c) > 0.05 and abs(math.cos(cmath.phase(z.page()))) > 0.05:
            break
    return err_num(math.tan(cmath.phase(z.page())), eta * P.s / (1 + eta * P.c))


# ---------------------------------------------------------------- 2.7 turning
group("turn", "Multiplying means turning (Theorem 4.46, Remark 4.47)",
      f"C_α = cos(α, {PB} − α) and S_α = sin(α, {PB} − α), computed from the definitions of Appendix A.")


@law("turn", "Thm. 4.46", "Product in polar form", "|z₁|(C_α + S_α u) · |z₂|(C_β + S_β u) = |z₁||z₂|(C_{α+β} + S_{α+β} u)", "leaning frame")
def t_turn_prod(P, rng):
    al, be = rangle(rng), rangle(rng)
    r1, r2 = rng.uniform(0.2, 3), rng.uniform(0.2, 3)
    return err_gn((unit(al, P) * r1) * (unit(be, P) * r2), unit(al + be, P) * (r1 * r2))


@law("turn", "Thm. 4.46", "Quotient in polar form", "z₁/z₂ = (|z₁|/|z₂|)(C_{α−β} + S_{α−β} u)", "leaning frame")
def t_turn_quot(P, rng):
    al, be = rangle(rng), rangle(rng)
    r1, r2 = rng.uniform(0.2, 3), rng.uniform(0.2, 3)
    return err_gn((unit(al, P) * r1) * (unit(be, P) * r2).inv(), unit(al - be, P) * (r1 / r2))


@law("turn", "Thm. 4.46", "Addition formula, cosine", "C_{α+β} = C_α C_β − S_α S_β", "numbers")
def t_add_cos(P, rng):
    al, be = rangle(rng, -6, 6), rangle(rng, -6, 6)
    return err_num(C(al + be, P), C(al, P) * C(be, P) - S(al, P) * S(be, P))


@law("turn", "Thm. 4.46", "Addition formula, sine", f"S_{{α+β}} = C_α S_β + S_α C_β + 2 cos {PB} S_α S_β", "numbers")
def t_add_sin(P, rng):
    al, be = rangle(rng, -6, 6), rangle(rng, -6, 6)
    return err_num(S(al + be, P), C(al, P) * S(be, P) + S(al, P) * C(be, P) + 2 * P.c * S(al, P) * S(be, P))


def _remark_angles(P, rng):
    while True:
        al, be = rangle(rng, -3, 3), rangle(rng, -3, 3)
        ph = P.p - al - be
        if min(abs(C(al, P)), abs(C(be, P)), abs(math.sin(al + ph)), abs(math.sin(be + ph)),
               abs(1 - gsin(al, ph) * gsin(be, ph)), abs(C(al, P) * C(be, P) - S(al, P) * S(be, P))) > 0.05:
            return al, be, ph


@law("turn", "Rem. 4.47 (1)", "Values at the projection angle of α + β",
     "sin(α, φ) = S_α/C_β,  cos(α, φ) = C_α − S_α S_β/C_β, and the same for β", "numbers")
def t_rem_1(P, rng):
    al, be, ph = _remark_angles(P, rng)
    Ca, Sa, Cb, Sb = C(al, P), S(al, P), C(be, P), S(be, P)
    return max(err_num(gsin(al, ph), Sa / Cb), err_num(gcos(al, ph), Ca - Sa * Sb / Cb),
               err_num(gsin(be, ph), Sb / Ca), err_num(gcos(be, ph), Cb - Sa * Sb / Ca))


@law("turn", "Rem. 4.47 (2)", "The sine formula from Appendix A",
     f"S_{{α+β}} = C_α S_β + S_α (1 − S_β²)/C_β,  (1 − S_β²)/C_β = C_β + 2 cos {PB} S_β", "numbers")
def t_rem_2(P, rng):
    al, be, ph = _remark_angles(P, rng)
    Ca, Sa, Cb, Sb = C(al, P), S(al, P), C(be, P), S(be, P)
    return max(err_num(gsin(al + be, ph), gsin(al, ph) + gcos(al, ph) * gsin(be, ph + al)),
               err_num(S(al + be, P), Ca * Sb + Sa * (1 - Sb ** 2) / Cb),
               err_num((1 - Sb ** 2) / Cb, Cb + 2 * P.c * Sb))


@law("turn", "Rem. 4.47 (3)", "The cosine formula from Appendix A",
     "cos(α+β, φ) = cos(α,φ) cos(β,φ) / (1 − sin(α,φ) sin(β,φ)) = C_α C_β − S_α S_β", "numbers")
def t_rem_3(P, rng):
    al, be, ph = _remark_angles(P, rng)
    Ca, Sa, Cb, Sb = C(al, P), S(al, P), C(be, P), S(be, P)
    v = (Ca - Sa * Sb / Cb) * (Cb - Sa * Sb / Ca) / (1 - (Sa / Cb) * (Sb / Ca))
    return max(err_num(gcos(al + be, ph), gcos(al, ph) * gcos(be, ph) / (1 - gsin(al, ph) * gsin(be, ph))),
               err_num(v, Ca * Cb - Sa * Sb), err_num(C(al + be, P), v))


# ---------------------------------------------------------------- 2.8 Euler and de Moivre
group("euler", "Euler and de Moivre (Theorems 4.48 and 4.50, Remark 4.49, Corollary 4.51)",
      "The exponential is summed from its power series in the leaning frame.")


@law("euler", "Thm. 4.48", "Euler's formula, from the power series",
     f"e^(θ{JU}) = cos(θ, {PB} − θ) + sin(θ, {PB} − θ) u", "power series in the frame", samples=25)
def t_euler_series(P, rng):
    th = rangle(rng, -4, 4)
    return err_gn(gn_exp(J_of(P) * th), unit(th, P))


@law("euler", "Thm. 4.48", "Euler's formula on the page", f"page(e^(θ{JU})) = e^(iθ)", FRAME_PAGE)
def t_euler_page(P, rng):
    th = rangle(rng, -10, 10)
    return err_num(unit(th, P).page(), cmath.exp(1j * th))


@law("euler", "Thm. 4.48 (3)", "Euler's identity", f"e^(π{JU}) + 1 = 0", "power series and formula", samples=1)
def t_euler_id(P, rng):
    return max(err_gn(gn_exp(J_of(P) * math.pi) + 1, GN(0.0, 0.0, P)), err_gn(unit(math.pi, P) + 1, GN(0.0, 0.0, P)))


@law("euler", "Thm. 4.48 (2)", "The two units as exponentials",
     f"{GU} = |{GU}| e^({PB}{JU}),  {JU} = e^((π/2){JU})", "power series and formula", samples=1)
def t_units_exp(P, rng):
    J = J_of(P)
    return max(err_gn(gn_exp(J * P.p) * P.L, G_of(P)), err_gn(unit(P.p, P) * P.L, G_of(P)),
               err_gn(gn_exp(J * (math.pi / 2)), J), err_gn(unit(math.pi / 2, P), J))


@law("euler", "Thm. 4.48 (1)", "Every number and its conjugate as exponentials",
     f"z = |z| e^(α{JU}),  z̄ = |z| e^(−α{JU})", "power series in the frame", samples=25)
def t_exp_form(P, rng):
    z = rgn(P, rng)
    al, M = cmath.phase(z.page()), abs(z.page())
    J = J_of(P)
    return max(err_gn(gn_exp(J * al) * M, z), err_gn(gn_exp(J * (-al)) * M, z.conj()))


@law("euler", "Rem. 4.49", "The derivatives along the Composer axis",
     f"dS/dθ = C/sin {PB} + S cot {PB},  dC/dθ = −C cot {PB} − S/sin {PB}", "finite differences", tol=TOL_DIFF)
def t_derivs(P, rng):
    th, h = rangle(rng), 1e-5
    dS = (S(th + h, P) - S(th - h, P)) / (2 * h)
    dC = (C(th + h, P) - C(th - h, P)) / (2 * h)
    return max(err_num(dS, C(th, P) / P.s + S(th, P) * P.c / P.s), err_num(dC, -C(th, P) * P.c / P.s - S(th, P) / P.s))


@law("euler", "Rem. 4.49", "A unit number turns at the rate of its quarter turn",
     f"dw/dθ = {JU} w,  w(θ) = C_θ + S_θ u", "finite differences", tol=TOL_DIFF)
def t_dw(P, rng):
    th, h = rangle(rng), 1e-5
    d = (unit(th + h, P) - unit(th - h, P)) * (1 / (2 * h))
    return err_gn(d, J_of(P) * unit(th, P))


@law("euler", "Thm. 4.50", "De Moivre's formula, n = −8 … 12", "(C_θ + S_θ u)ⁿ = C_{nθ} + S_{nθ} u", "repeated products in the frame", samples=10)
def t_demoivre(P, rng):
    th = rangle(rng)
    w = unit(th, P)
    return max(err_gn(w.power(n), unit(n * th, P)) for n in range(-8, 13))


@law("euler", "Cor. 4.51", "The powers of the Composer Unit, n = −6 … 12",
     f"{GU}ⁿ = |{GU}|ⁿ cos(n{PB}, (1−n){PB}) + |{GU}|ⁿ⁻¹ sin(n{PB}, (1−n){PB}) {GU}", "repeated products in the frame", samples=1)
def t_powers(P, rng):
    G, L = G_of(P), P.L
    worst = 0.0
    for n in range(-6, 13):
        g = G.power(n)
        f1 = GN(L ** n * gcos(n * P.p, (1 - n) * P.p), L ** (n - 1) * gsin(n * P.p, (1 - n) * P.p), P)
        f2 = GN(-L ** n * math.sin((n - 1) * P.p) / P.s, L ** (n - 1) * math.sin(n * P.p) / P.s, P)
        worst = max(worst, err_gn(g, f1), err_gn(g, f2))
    return worst


@law("euler", "Cor. 4.51", "The step rule", f"{GU}ⁿ⁺¹ = 2|{GU}| cos {PB} {GU}ⁿ − |{GU}|² {GU}ⁿ⁻¹", "repeated products in the frame", samples=1)
def t_step(P, rng):
    G, L = G_of(P), P.L
    return max(err_gn(G.power(n + 1), G.power(n) * (2 * L * P.c) - G.power(n - 1) * L ** 2) for n in range(-5, 12))


@law("euler", "Cor. 4.51", "n = 2 is the Rule of the Square", f"|{GU}|² cos(2{PB}, −{PB}) = −|{GU}|²,  |{GU}| sin(2{PB}, −{PB}) = 2|{GU}| cos {PB}", "numbers", samples=1)
def t_n2(P, rng):
    L = P.L
    return err_gn(GN(L ** 2 * gcos(2 * P.p, -P.p), L * gsin(2 * P.p, -P.p), P), GN(-L ** 2, 2 * L * P.c, P))


@law("euler", "Thm. 4.1", f"On the page, {GU}ⁿ is the arrow |{GU}|ⁿ at the angle n{PB}", f"page({GU}ⁿ) = |{GU}|ⁿ e^(in{PB})", FRAME_PAGE, samples=1)
def t_thm41(P, rng):
    G = G_of(P)
    return max(err_num(G.power(n).page(), P.L ** n * cmath.exp(1j * n * P.p)) for n in range(1, 13))


# ---------------------------------------------------------------- 2.9 roots
group("roots", "Roots (Theorem 4.52, Corollary 4.53)")


@law("roots", "Thm. 4.52", "Each root, raised to the power n, gives z back",
     "w_k = |z|^(1/n) (C_{θk} + S_{θk} u),  θ_k = (α + 2πk)/n,  w_kⁿ = z", "repeated products in the frame", samples=10)
def t_roots(P, rng):
    z = rgn(P, rng)
    al, M = cmath.phase(z.page()), abs(z.page())
    worst = 0.0
    for n in range(1, 9):
        for k in range(n):
            w = unit((al + 2 * math.pi * k) / n, P) * M ** (1 / n)
            worst = max(worst, err_gn(w.power(n), z))
    return worst


@law("roots", "Thm. 4.52", "Exactly n roots, all different",
     "the n formula roots are the n roots that a root finder finds for wⁿ − z", "independent root finder", samples=5)
def t_roots_all(P, rng):
    z = rgn(P, rng)
    al, M = cmath.phase(z.page()), abs(z.page())
    worst = 0.0
    for n in range(2, 8):
        ours = [(unit((al + 2 * math.pi * k) / n, P) * M ** (1 / n)).page() for k in range(n)]
        found = poly_roots([-z.page()] + [0j] * (n - 1))
        match = [min(range(n), key=lambda j: abs(f - ours[j])) for f in found]
        if len(set(match)) != n:
            return 1.0
        worst = max(worst, max(abs(f - ours[j]) / max(1, abs(ours[j])) for f, j in zip(found, match)))
        if min(abs(ours[i] - ours[j]) for i in range(n) for j in range(i)) < 1e-6:
            return 1.0
    return worst


@law("roots", "Thm. 4.52", "The roots of unity lie on the ring of magnitude 1",
     "C² + 2 cos φ̄ C S + S² = 1 and (C + S u)ⁿ = 1 at θ = 2πk/n", "frame and numbers", samples=1)
def t_unity(P, rng):
    worst = 0.0
    for n in range(1, 13):
        for k in range(n):
            t = 2 * math.pi * k / n
            worst = max(worst, err_num(C(t, P) ** 2 + 2 * P.c * C(t, P) * S(t, P) + S(t, P) ** 2, 1.0),
                        err_gn(unit(t, P).power(n), GN(1.0, 0.0, P)))
    return worst


@law("roots", "Cor. 4.53", "The square root of the Composer Unit",
     f"[√|{GU}| / (2 cos({PB}/2))] (1 + u), squared, is {GU}", "leaning frame", samples=1)
def t_sqrtG(P, rng):
    k = math.sqrt(P.L) / (2 * math.cos(P.p / 2))
    w = GN(k, k / P.L, P)
    return max(err_gn(w * w, G_of(P)), err_gn(unit(P.p / 2, P) * math.sqrt(P.L), w))


@law("roots", "Cor. 4.53", "It lies on the ray of equal parts",
     f"real part = length of the Composer part, angle = {PB}/2", FRAME_PAGE, samples=1)
def t_sqrtG_ray(P, rng):
    k = math.sqrt(P.L) / (2 * math.cos(P.p / 2))
    w = GN(k, k / P.L, P)
    return max(err_num(w.a, w.b * P.L), err_num(cmath.phase(w.page()), P.p / 2))


# ---------------------------------------------------------------- 2.10 the parts beside the shadow and the height
group("parts", "The parts beside the shadow and the height (Subsection 4.9.7 and Appendix F.8)",
      f"Each law gives the real part Re(z) = a and the Composer part Com(z) = b{GU} beside the shadow Sh(z) and the height Ht(z). "
      f"The formula of the book is computed in the leaning frame; the other side is read on the page, from the arrow itself.")


@law("parts", "Prop. 4.54 (1)", "The parts from the shadow and the height",
     f"Re(z) = Sh(z) − Ht(z) cot {PB},  Com(z) = Ht(z) {GU}/(|{GU}| sin {PB})", FRAME_PAGE)
def t_p_dict(P, rng):
    z = rgn(P, rng)
    w = z.page()
    sh, ht = w.real, w.imag                       # read on the page
    return max(err_num(sh - ht * P.c / P.s, z.a), err_num(ht / (P.L * P.s), z.b))


@law("parts", "Prop. 4.54 (1)", "The shadow and the height from the parts",
     f"Sh(z) = Re(z) + Sh(Com z) = Re(z) + b|{GU}| cos {PB},  Ht(z) = Ht(Com z) = b|{GU}| sin {PB}", FRAME_PAGE)
def t_p_dict2(P, rng):
    z = rgn(P, rng)
    w, com = z.page(), GN(0.0, z.b, P).page()
    return max(err_num(z.a + com.real, w.real), err_num(com.imag, w.imag),
               err_num(z.shadow(), w.real), err_num(z.height(), w.imag))


@law("parts", "Prop. 4.54 (2)", "All four are additive and keep real factors",
     "Re, Com, Sh and Ht of z₁ + z₂ and of kz, k real", FRAME_PAGE)
def t_p_lin(P, rng):
    z1, z2, k = rgn(P, rng), rgn(P, rng), rng.uniform(-3, 3)
    w1, w2 = z1.page(), z2.page()
    s, ks = parts_page(w1 + w2, P), parts_page(k * w1, P)
    r1, r2 = parts_page(w1, P), parts_page(w2, P)
    return max(max(err_num(s[i], r1[i] + r2[i]) for i in range(4)), max(err_num(ks[i], k * r1[i]) for i in range(4)),
               err_gn(z1 + z2, from_page(w1 + w2, P)), err_gn(z1 * k, from_page(k * w1, P)))


@law("parts", "Prop. 4.54 (3)", "The conjugate keeps both coordinates and mirrors the unit",
     f"z̄ = a + b {GU}̄", FRAME_PAGE)
def t_p_conjunit(P, rng):
    z = rgn(P, rng)
    return max(err_num(z.conj().page(), z.a + z.b * P.Gpage.conjugate()),
               err_gn(z.conj(), GN(z.a, 0.0, P) + G_of(P).conj() * z.b))


@law("parts", "Prop. 4.54 (3)", "The parts read from z and z̄",
     f"Re(z) = ({GU}z̄ − {GU}̄z)/({GU} − {GU}̄),  Com(z) = {GU}(z − z̄)/({GU} − {GU}̄),  {GU} − {GU}̄ = 2|{GU}| sin {PB} {JU}",
     "page and leaning frame")
def t_p_fromconj(P, rng):
    z = rgn(P, rng)
    w, G = z.page(), P.Gpage
    Gb = G.conjugate()
    re_page = (G * w.conjugate() - Gb * w) / (G - Gb)              # on the page, with Python's complex numbers
    com_page = G * (w - w.conjugate()) / (G - Gb)
    g = G_of(P)
    gb = g.conj()
    d = (g - gb).inv()
    re_frame = (g * z.conj() - gb * z) * d                          # in the leaning frame
    com_frame = g * (z - z.conj()) * d
    return max(err_num(re_page, z.a), err_num(com_page, z.b * G), err_num(G - Gb, 2 * P.L * P.s * 1j),
               err_gn(re_frame, GN(z.a, 0.0, P)), err_gn(com_frame, GN(0.0, z.b, P)),
               err_gn(g - gb, J_of(P) * (2 * P.L * P.s)))


@law("parts", "Prop. 4.54 (4)", "Real numbers and pure Composer numbers",
     f"z − z̄ = b({GU} − {GU}̄): z̄ = z exactly when b = 0;  {GU}z̄ − {GU}̄z = a({GU} − {GU}̄): {GU}z̄ = {GU}̄z exactly when Re(z) = 0",
     "leaning frame and page")
def t_p_pure(P, rng):
    z = rgn(P, rng)
    g = G_of(P)
    gb = g.conj()
    e1 = err_gn(z - z.conj(), (g - gb) * z.b)
    e2 = err_gn(g * z.conj() - gb * z, (g - gb) * z.a)
    w = rng.uniform(-3, 3) * P.Gpage                               # a pure Composer number, on the page
    e3 = err_num(P.Gpage * w.conjugate(), P.Gpage.conjugate() * w)
    t = rng.uniform(-3, 3)
    e4 = err_gn((J_of(P) * t).conj(), -(J_of(P) * t))               # on the upright axis z̄ = −z (Appendix F.8)
    return max(e1, e2, e3, e4)


@law("parts", "Thm. 4.55 (1)", "The parts of a product",
     f"Re(z₁z₂) = Re(z₁)Re(z₂) − |{GU}|² b₁b₂,  Com(z₁z₂) = Sh(z₁)Com(z₂) + Sh(z₂)Com(z₁)", FRAME_PAGE)
def t_p_prod(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    pa, pb_, _, _ = parts_page(z1.page() * z2.page(), P)
    return max(err_num(z1.a * z2.a - P.L ** 2 * z1.b * z2.b, pa),
               err_num(z1.shadow() * z2.b + z2.shadow() * z1.b, pb_))


@law("parts", "Thm. 4.55 (2)", "The shadow and the height of a product",
     "Sh(z₁z₂) = Sh(z₁)Sh(z₂) − Ht(z₁)Ht(z₂),  Ht(z₁z₂) = Sh(z₁)Ht(z₂) + Sh(z₂)Ht(z₁)", FRAME_PAGE)
def t_p_prod_sh(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    w = z1.page() * z2.page()
    s1, s2, h1, h2 = z1.shadow(), z2.shadow(), z1.height(), z2.height()
    return max(err_num(s1 * s2 - h1 * h2, w.real), err_num(s1 * h2 + s2 * h1, w.imag))


@law("parts", "Thm. 4.55 (3)", "The parts, the shadow and the height of a square",
     f"Re(z²) = Re(z)² − |{GU}|²b²,  Com(z²) = 2 Sh(z) Com(z);  Sh(z²) = Sh(z)² − Ht(z)²,  Ht(z²) = 2 Sh(z) Ht(z)", FRAME_PAGE)
def t_p_square(P, rng):
    z = rgn(P, rng)
    w = z.page() ** 2
    pa, pb_, _, _ = parts_page(w, P)
    s, h = z.shadow(), z.height()
    return max(err_num(z.a ** 2 - P.L ** 2 * z.b ** 2, pa), err_num(2 * s * z.b, pb_),
               err_num(s * s - h * h, w.real), err_num(2 * s * h, w.imag))


@law("parts", "Prop. 4.56 (1)", "The parts, the shadow and the height of the conjugate",
     "Re(z̄) = 2 Sh(z) − Re(z),  Com(z̄) = −Com(z);  Sh(z̄) = Sh(z),  Ht(z̄) = −Ht(z)", FRAME_PAGE)
def t_p_conj(P, rng):
    z = rgn(P, rng)
    w = z.page().conjugate()                       # the mirror image, on the page
    pa, pb_, _, _ = parts_page(w, P)
    s, h = z.shadow(), z.height()
    return max(err_num(2 * s - z.a, pa), err_num(-z.b, pb_), err_num(s, w.real), err_num(-h, w.imag))


@law("parts", "Prop. 4.56 (1)", "The shadow lies halfway between the two real parts",
     "Sh(z) = [Re(z) + Re(z̄)]/2, and z is as far from Re(z) as from Re(z̄): the isosceles triangle", "page")
def t_p_mid(P, rng):
    z = rgn(P, rng)
    w = z.page()
    r1 = from_page(w, P).a
    r2 = from_page(w.conjugate(), P).a
    return max(err_num((r1 + r2) / 2, w.real), err_num(abs(w - r1), abs(w - r2)), err_num(abs(w - r1), abs(z.b) * P.L))


@law("parts", "Prop. 4.56 (2)", "The inverse",
     "Re(1/z) = [2 Sh(z) − Re(z)]/|z|²,  Com(1/z) = −Com(z)/|z|²;  Sh(1/z) = Sh(z)/|z|²,  Ht(1/z) = −Ht(z)/|z|²", FRAME_PAGE)
def t_p_inv(P, rng):
    z = rgn(P, rng)
    w = 1 / z.page()
    pa, pb_, _, _ = parts_page(w, P)
    n2, s, h = z.norm2(), z.shadow(), z.height()
    return max(err_num((2 * s - z.a) / n2, pa), err_num(-z.b / n2, pb_), err_num(s / n2, w.real), err_num(-h / n2, w.imag))


@law("parts", "Prop. 4.56 (3)", "The quotient",
     f"Re(z₁/z₂) = [Re(z₁)Re(z̄₂) + |{GU}|²b₁b₂]/|z₂|²,  Com(z₁/z₂) = [Sh(z₂)Com(z₁) − Sh(z₁)Com(z₂)]/|z₂|²;  "
     "Sh(z₁/z₂) = [Sh(z₁)Sh(z₂) + Ht(z₁)Ht(z₂)]/|z₂|²,  Ht(z₁/z₂) = [Sh(z₂)Ht(z₁) − Sh(z₁)Ht(z₂)]/|z₂|²", FRAME_PAGE)
def t_p_quot(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    w = z1.page() / z2.page()
    pa, pb_, _, _ = parts_page(w, P)
    n2 = z2.norm2()
    s1, s2, h1, h2 = z1.shadow(), z2.shadow(), z1.height(), z2.height()
    return max(err_num((z1.a * z2.conj().a + P.L ** 2 * z1.b * z2.b) / n2, pa),
               err_num((s2 * z1.b - s1 * z2.b) / n2, pb_),
               err_num((s1 * s2 + h1 * h2) / n2, w.real), err_num((s2 * h1 - s1 * h2) / n2, w.imag))


@law("parts", "Cor. 4.57 (1)", f"A quarter turn, the product by {JU}",
     f"Re({JU}z) = −Ht(z) − Sh(z) cot {PB},  Com({JU}z) = Sh(z) {GU}/(|{GU}| sin {PB});  Sh({JU}z) = −Ht(z),  Ht({JU}z) = Sh(z)",
     FRAME_PAGE)
def t_p_J(P, rng):
    z = rgn(P, rng)
    w = 1j * z.page()                              # a quarter turn, on the page
    pa, pb_, _, _ = parts_page(w, P)
    s, h = z.shadow(), z.height()
    want = GN(-h - s * P.c / P.s, s / (P.L * P.s), P)
    return max(err_num(want.a, pa), err_num(want.b, pb_), err_num(-h, w.real), err_num(s, w.imag),
               err_gn(J_of(P) * z, want))


@law("parts", "Cor. 4.57 (2)", "The product by the Composer Unit",
     f"Re({GU}z) = −|{GU}|²b,  Com({GU}z) = Re(z̄){GU};  Sh({GU}z) = |{GU}|(Sh cos {PB} − Ht sin {PB}),  Ht({GU}z) = |{GU}|(Sh sin {PB} + Ht cos {PB})",
     FRAME_PAGE)
def t_p_G(P, rng):
    z = rgn(P, rng)
    w = P.Gpage * z.page()
    pa, pb_, _, _ = parts_page(w, P)
    s, h = z.shadow(), z.height()
    return max(err_num(-P.L ** 2 * z.b, pa), err_num(z.conj().a, pb_),
               err_num(P.L * (s * P.c - h * P.s), w.real), err_num(P.L * (s * P.s + h * P.c), w.imag),
               err_gn(G_of(P) * z, GN(-P.L ** 2 * z.b, z.conj().a, P)))


@law("parts", "Prop. 4.58", "The signed area",
     f"Com(z̄₁z₂) = (a₁b₂ − b₁a₂){GU},  Ht(z̄₁z₂) = Sh(z₁)Ht(z₂) − Ht(z₁)Sh(z₂) = (a₁b₂ − b₁a₂)|{GU}| sin {PB}, the area on the page",
     FRAME_PAGE)
def t_p_area(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    w1, w2 = z1.page(), z2.page()
    w = w1.conjugate() * w2
    _, pb_, _, _ = parts_page(w, P)
    det = z1.a * z2.b - z1.b * z2.a
    cross = w1.real * w2.imag - w1.imag * w2.real                  # the signed area of the parallelogram, on the page
    return max(err_num(det, pb_), err_num(z1.shadow() * z2.height() - z1.height() * z2.shadow(), w.imag),
               err_num(det * P.L * P.s, cross), err_num(cross, w.imag))


@law("parts", "Prop. 4.58", "The dot product",
     f"Re(z̄₁z₂) = Re(z̄₁)Re(z₂) + |{GU}|²b₁b₂,  Sh(z̄₁z₂) = Sh(z₁)Sh(z₂) + Ht(z₁)Ht(z₂)", FRAME_PAGE)
def t_p_dot(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    w1, w2 = z1.page(), z2.page()
    w = w1.conjugate() * w2
    pa, _, _, _ = parts_page(w, P)
    dot = w1.real * w2.real + w1.imag * w2.imag
    return max(err_num(z1.conj().a * z2.a + P.L ** 2 * z1.b * z2.b, pa),
               err_num(z1.shadow() * z2.shadow() + z1.height() * z2.height(), dot), err_num(dot, w.real))


@law("parts", "Prop. 4.58", "One line through 0, and a right angle",
     f"z₂ = t z₁ gives a₁b₂ = b₁a₂;  z₂ = t {JU} z₁ gives Sh(z̄₁z₂) = 0", "page")
def t_p_lines(P, rng):
    z1, t = rgn(P, rng), rng.uniform(-3, 3)
    z2 = from_page(t * z1.page(), P)
    z3 = from_page(t * 1j * z1.page(), P)
    w = z1.page().conjugate() * z3.page()
    return max(abs(z1.a * z2.b - z1.b * z2.a) / max(1.0, abs(z1.a * z2.b), abs(z1.b * z2.a)), abs(w.real) / max(1.0, abs(w)))


@law("parts", "Prop. 4.23, Prop. 4.58", "Two laws carry no angle",
     f"with |{GU}| held fixed, Re(z₁z₂) and the Composer coordinate of z̄₁z₂ do not change when {PB} changes", "page, in two planes")
def t_p_noangle(P, rng):
    Q = Plane(rng.uniform(0.1, math.pi - 0.1), P.L)                # another plane with the same length of the Composer Unit
    a, b, c, d = (rng.uniform(-4, 4) for _ in range(4))
    vals = []
    for R in (P, Q):
        w1, w2 = GN(a, b, R).page(), GN(c, d, R).page()
        vals.append((from_page(w1 * w2, R).a, from_page(w1.conjugate() * w2, R).b))
    return max(err_num(vals[0][0], vals[1][0]), err_num(vals[0][1], vals[1][1]),
               err_num(vals[0][0], a * c - P.L ** 2 * b * d), err_num(vals[0][1], a * d - b * c))


@law("parts", "Cor. 4.59", "Two trigonometries: the powers, n = −4 … 7",
     "Re(zⁿ) = |z|ⁿ C_{nα},  Com(zⁿ) = |z|ⁿ S_{nα} u;  Sh(zⁿ) = |z|ⁿ cos nα,  Ht(zⁿ) = |z|ⁿ sin nα", "page", samples=20)
def t_p_pow(P, rng):
    z = rgn(P, rng, -1.5, 1.5)
    w = z.page()
    r, al = abs(w), cmath.phase(w)
    worst = 0.0
    for n in range(-4, 8):
        wn = w ** n
        pa, pb_, _, _ = parts_page(wn, P)
        worst = max(worst, err_num(r ** n * C(n * al, P), pa), err_num(r ** n * S(n * al, P) / P.L, pb_),
                    err_num(r ** n * math.cos(n * al), wn.real), err_num(r ** n * math.sin(n * al), wn.imag))
    return worst


@law("parts", "Cor. 4.59", "Two trigonometries: the exponential",
     "Re(e^(z)) = exp(Sh z) C_(Ht z),  Com(e^(z)) = exp(Sh z) S_(Ht z) u;  Sh(e^(z)) = exp(Sh z) cos(Ht z),  Ht(e^(z)) = exp(Sh z) sin(Ht z)",
     "page and power series")
def t_p_exp(P, rng):
    z = rgn(P, rng, -1.5, 1.5)
    w = cmath.exp(z.page())
    pa, pb_, _, _ = parts_page(w, P)
    s, h = z.shadow(), z.height()
    es = math.exp(s)
    f = gn_exp(z)                                  # the power series, summed in the leaning frame
    return max(err_num(es * C(h, P), pa), err_num(es * S(h, P) / P.L, pb_),
               err_num(es * math.cos(h), w.real), err_num(es * math.sin(h), w.imag),
               err_num(f.a, pa), err_num(f.b, pb_))


@law("parts", "Rem. 4.60", "Parts longer than the number",
     f"|Sh(z)|, |Ht(z)| ≤ |z|, but |Re(z)|, |Com(z)| ≤ |z|/sin {PB}, reached at α = {PB} − π/2 and at α = π/2", "page")
def t_p_long(P, rng):
    z = rgn(P, rng)
    w = z.page()
    m = abs(w)
    pa, pb_, _, _ = parts_page(w, P)
    over = max(0.0, abs(w.real) - m, abs(w.imag) - m, abs(pa) - m / P.s, abs(pb_) * P.L - m / P.s) / max(1.0, m)
    e1 = cmath.exp(1j * (P.p - math.pi / 2))       # at a right angle to the Composer axis: the real part is longest
    e2 = 1j                                        # at a right angle to the real axis, as the Turning Unit: the Composer part is longest
    reach = max(err_num(abs(from_page(e1, P).a), 1 / P.s), err_num(abs(from_page(e2, P).b) * P.L, 1 / P.s))
    return max(over, reach)


@law("parts", "App. F.8", "The magnitude in the two pairs",
     f"|z|² = Re(z)² + 2 cos {PB} Re(z) b|{GU}| + b²|{GU}|² = Sh(z)² + Ht(z)²", FRAME_PAGE)
def t_p_mag(P, rng):
    z = rgn(P, rng)
    m2 = abs(z.page()) ** 2
    return max(err_num(z.a ** 2 + 2 * P.c * z.a * z.b * P.L + z.b ** 2 * P.L ** 2, m2),
               err_num(z.shadow() ** 2 + z.height() ** 2, m2))


@law("parts", "App. F.8", "The powers of the Composer Unit, n = −4 … 9",
     f"Re({GU}ⁿ) = −|{GU}|ⁿ sin((n − 1){PB})/sin {PB},  Com({GU}ⁿ) = |{GU}|ⁿ⁻¹ (sin n{PB}/sin {PB}) {GU};  "
     f"Sh({GU}ⁿ) = |{GU}|ⁿ cos n{PB},  Ht({GU}ⁿ) = |{GU}|ⁿ sin n{PB}", FRAME_PAGE, samples=1)
def t_p_Gn(P, rng):
    worst = 0.0
    for n in range(-4, 10):
        w = P.Gpage ** n
        pa, pb_, _, _ = parts_page(w, P)
        f = G_of(P).power(n)
        worst = max(worst, err_num(-P.L ** n * math.sin((n - 1) * P.p) / P.s, pa),
                    err_num(P.L ** (n - 1) * math.sin(n * P.p) / P.s, pb_),
                    err_num(P.L ** n * math.cos(n * P.p), w.real), err_num(P.L ** n * math.sin(n * P.p), w.imag),
                    err_gn(f, GN(pa, pb_, P)))
    return worst


@law("parts", "App. F.8", "The table of six numbers",
     f"1, {GU}, {GU}̄, {GU}², 1/{GU} and {JU}: their real part, Composer coordinate, shadow and height", FRAME_PAGE, samples=1)
def t_p_table(P, rng):
    L, c, s = P.L, P.c, P.s
    G = P.Gpage
    rows = [(1 + 0j, (1, 0, 1, 0)), (G, (0, 1, L * c, L * s)), (G.conjugate(), (2 * L * c, -1, L * c, -L * s)),
            (G * G, (-L * L, 2 * L * c, L * L * math.cos(2 * P.p), L * L * math.sin(2 * P.p))),
            (1 / G, (2 * c / L, -1 / L ** 2, c / L, -s / L)), (1j, (-c / s, 1 / (L * s), 0, 1))]
    worst = 0.0
    for w, want in rows:
        got = parts_page(w, P)
        worst = max(worst, max(err_num(g_, w_) for g_, w_ in zip(got, want)))
    return worst


# ---------------------------------------------------------------- 2.11 Appendix F
group("appf", "More laws of Appendix F")


@law("appf", "App. F.3", "The opposite angle", f"C_{{−θ}} = C_θ + 2 cos {PB} S_θ,  S_{{−θ}} = −S_θ", "numbers")
def t_opp(P, rng):
    t = rangle(rng, -6, 6)
    return max(err_num(C(-t, P), C(t, P) + 2 * P.c * S(t, P)), err_num(S(-t, P), -S(t, P)))


@law("appf", "App. F.3", "The conjugate in polar form", "z̄ = |z| (C_{−α} + S_{−α} u)", "leaning frame")
def t_conj_polar(P, rng):
    z = rgn(P, rng)
    al, M = cmath.phase(z.page()), abs(z.page())
    return err_gn(z.conj(), unit(-al, P) * M)


@law("appf", "App. F.3", "The mirror in the bisector of the two axes",
     "u z̄ sends x + y u to y + x u", "leaning frame")
def t_mirror2(P, rng):
    x, y = rng.uniform(-4, 4), rng.uniform(-4, 4)
    z = GN(x, y / P.L, P)
    return err_gn(u_of(P) * z.conj(), GN(y, x / P.L, P))


@law("appf", "App. F.3", "Reciprocity in the plane", f"C_{{{PB}−θ}} = S_θ,  S_{{{PB}−θ}} = C_θ", "numbers")
def t_recip_plane(P, rng):
    t = rangle(rng, -6, 6)
    return max(err_num(C(P.p - t, P), S(t, P)), err_num(S(P.p - t, P), C(t, P)))


@law("appf", "App. F.3", "The angle", f"tan θ = Ht(z)/Sh(z) = b|{GU}| sin {PB} / (a + b|{GU}| cos {PB})", FRAME_PAGE)
def t_angle(P, rng):
    while True:
        z = rgn(P, rng)
        if abs(z.shadow()) > 0.05:
            break
    return err_num(math.tan(cmath.phase(z.page())), z.b * P.L * P.s / (z.a + z.b * P.L * P.c))


@law("appf", "App. F.4", "Subtraction, cosine", f"C_{{α−β}} = C_α C_β + S_α S_β + 2 cos {PB} C_α S_β", "numbers")
def t_sub_cos(P, rng):
    al, be = rangle(rng, -6, 6), rangle(rng, -6, 6)
    return err_num(C(al - be, P), C(al, P) * C(be, P) + S(al, P) * S(be, P) + 2 * P.c * C(al, P) * S(be, P))


@law("appf", "App. F.4", "Subtraction, sine", "S_{α−β} = S_α C_β − C_α S_β", "numbers")
def t_sub_sin(P, rng):
    al, be = rangle(rng, -6, 6), rangle(rng, -6, 6)
    return err_num(S(al - be, P), S(al, P) * C(be, P) - C(al, P) * S(be, P))


@law("appf", "App. F.4", "Double angle", f"C_{{2α}} = C_α² − S_α²,  S_{{2α}} = 2 S_α (C_α + cos {PB} S_α)", "numbers")
def t_double(P, rng):
    al = rangle(rng, -6, 6)
    Ca, Sa = C(al, P), S(al, P)
    return max(err_num(C(2 * al, P), Ca ** 2 - Sa ** 2), err_num(S(2 * al, P), 2 * Sa * (Ca + P.c * Sa)))


@law("appf", "App. F.4", f"A turn by {PB}, the product by u", f"C_{{θ+{PB}}} = −S_θ,  S_{{θ+{PB}}} = C_θ + 2 cos {PB} S_θ", "numbers and frame")
def t_turn_phi(P, rng):
    t = rangle(rng)
    return max(err_num(C(t + P.p, P), -S(t, P)), err_num(S(t + P.p, P), C(t, P) + 2 * P.c * S(t, P)),
               err_gn(u_of(P) * unit(t, P), unit(t + P.p, P)))


@law("appf", "App. F.4", f"A quarter turn, the product by {JU}",
     f"C_{{θ+π/2}} = −C_θ cot {PB} − S_θ/sin {PB},  S_{{θ+π/2}} = C_θ/sin {PB} + S_θ cot {PB}", "numbers and frame")
def t_turn_quarter(P, rng):
    t = rangle(rng)
    Ct, St = C(t, P), S(t, P)
    return max(err_num(C(t + math.pi / 2, P), -Ct * P.c / P.s - St / P.s),
               err_num(S(t + math.pi / 2, P), Ct / P.s + St * P.c / P.s),
               err_gn(J_of(P) * unit(t, P), unit(t + math.pi / 2, P)))


@law("appf", "App. F.4", "A half turn, the product by −1", "C_{θ+π} = −C_θ,  S_{θ+π} = −S_θ", "numbers")
def t_turn_half(P, rng):
    t = rangle(rng)
    return max(err_num(C(t + math.pi, P), -C(t, P)), err_num(S(t + math.pi, P), -S(t, P)))


@law("appf", "App. F.5", "The exponential of a General Number",
     "e^z = e^(Sh(z)) (C_{Ht(z)} + S_{Ht(z)} u)", "power series and page", samples=25)
def t_exp(P, rng):
    z = GN(rng.uniform(-2, 2), rng.uniform(-2, 2), P)
    f = unit(z.height(), P) * math.exp(z.shadow())
    return max(err_gn(gn_exp(z), f), err_num(f.page(), cmath.exp(z.page())))


@law("appf", "App. F.5", "The logarithm", f"ln z = (ln|z| − θ cot {PB}) + θ/(|{GU}| sin {PB}) {GU}, and e^(ln z) = z", "power series and page", samples=25)
def t_log(P, rng):
    z = rgn(P, rng)
    th, M = cmath.phase(z.page()), abs(z.page())
    lg = GN(math.log(M) - th * P.c / P.s, th / (P.L * P.s), P)
    return max(err_gn(gn_exp(lg), z), err_num(lg.page(), cmath.log(z.page())))


@law("appf", "App. F.5", "The numbers sin n" + PB + "/sin " + PB + " (Chebyshev)",
     f"sin n{PB}/sin {PB} = U_(n−1)(cos {PB}): 1, 2cos {PB}, 4cos²{PB} − 1, 8cos³{PB} − 4cos {PB}, …", "numbers", samples=1)
def t_cheb(P, rng):
    x = P.c
    U = [1.0, 2 * x]
    for _ in range(12):
        U.append(2 * x * U[-1] - U[-2])
    worst = max(err_num(math.sin(n * P.p) / P.s, U[n - 1]) for n in range(1, 13))
    explicit = [1.0, 2 * x, 4 * x ** 2 - 1, 8 * x ** 3 - 4 * x]
    return max(worst, max(err_num(math.sin(n * P.p) / P.s, explicit[n - 1]) for n in range(1, 5)))


@law("appf", "App. F.6", "One step along each axis, squared", f"(1 + u)² = 4 cos²({PB}/2) u", "leaning frame", samples=1)
def t_one_plus_u(P, rng):
    w = GN(1.0, 0.0, P) + u_of(P)
    return err_gn(w * w, u_of(P) * (4 * math.cos(P.p / 2) ** 2))


@law("appf", "App. F.6", "The quadratic equation",
     "v² + Pv + Q = 0 at v = (−P ± √(P² − 4Q))/2, square root from Thm. 4.52", "frame and root finder", samples=25)
def t_quadratic(P, rng):
    Pq, Qq = rgn(P, rng), rgn(P, rng)
    D = Pq * Pq - Qq * 4
    al, M = cmath.phase(D.page()), abs(D.page())
    r = unit(al / 2, P) * math.sqrt(M)
    vs = [(-Pq + r) * 0.5, (-Pq - r) * 0.5]
    scale = max(1.0, abs(Pq.page()) ** 2, abs(Qq.page()))
    worst = max(err_gn(v * v + Pq * v + Qq, GN(0.0, 0.0, P)) / scale for v in vs)
    found = poly_roots([Qq.page(), Pq.page()])
    pg = [v.page() for v in vs]
    worst = max(worst, min(max(abs(found[0] - pg[0]), abs(found[1] - pg[1])),
                           max(abs(found[0] - pg[1]), abs(found[1] - pg[0]))) / max(1.0, abs(pg[0]), abs(pg[1])))
    return worst


def _mat(z):
    P = z.P
    return [[z.a, -P.L ** 2 * z.b], [z.b, z.a + 2 * P.L * z.b * P.c]]


def _mmul(A, B):
    return [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
            [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]]


def _merr(A, B):
    return max(abs(A[i][j] - B[i][j]) for i in range(2) for j in range(2)) / max(1.0, max(abs(B[i][j]) for i in range(2) for j in range(2)))


@law("appf", "App. F.7", "The matrix of a General Number",
     f"[[a, −|{GU}|²b], [b, a + 2b|{GU}| cos {PB}]]: products, det = |z|², trace = 2 Sh(z)", "matrices")
def t_matrix(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    A = _mat(z1)
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    return max(_merr(_mmul(_mat(z1), _mat(z2)), _mat(z1 * z2)), err_num(det, z1.norm2()),
               err_num(A[0][0] + A[1][1], 2 * z1.shadow()))


@law("appf", "App. F.7", "The matrix of the Composer Unit obeys the Rule of the Square",
     f"M² = 2|{GU}| cos {PB} M − |{GU}|²", "matrices", samples=1)
def t_matrix_G(P, rng):
    M = _mat(G_of(P))
    M2 = _mmul(M, M)
    R = [[2 * P.L * P.c * M[i][j] - (P.L ** 2 if i == j else 0.0) for j in range(2)] for i in range(2)]
    return _merr(M2, R)


@law("appf", "App. F.1", "To another relative plane",
     f"b′ = Ht(z)/(|{GU}_ψ̄| sin ψ̄),  a′ = Sh(z) − Ht(z) cot ψ̄: the same point", FRAME_PAGE)
def t_change(P, rng):
    z = rgn(P, rng)
    while True:
        q = rng.uniform(0.1, math.pi - 0.1)
        if abs(q - 3 * math.pi / 4) > 0.1:
            break
    Q = Plane(q, rng.uniform(0.3, 3.0))
    w = GN(z.shadow() - z.height() * Q.c / Q.s, z.height() / (Q.L * Q.s), Q)
    return err_num(w.page(), z.page())


@law("appf", "App. F.1", "From the complex plane",
     f"x + iy = (x − y cot {PB}) + y/(|{GU}| sin {PB}) {GU}", FRAME_PAGE)
def t_from_c(P, rng):
    x, y = rng.uniform(-4, 4), rng.uniform(-4, 4)
    return err_num(GN(x - y * P.c / P.s, y / (P.L * P.s), P).page(), complex(x, y))


# ---------------------------------------------------------------- 2.12 Appendix A
group("appa", "The formulas of Appendix A that the section uses",
      "Random angles; cases near a pole of a formula are skipped.")


def _angles(rng, k=4):
    return [rng.uniform(-3, 3) for _ in range(k)]


def _safe(*vals):
    return all(abs(v) > 0.05 for v in vals)


def _appa(fn):
    def run(rng):
        while True:
            al, be, ph, th = _angles(rng)
            out = fn(al, be, ph, th)
            if out is not None:
                return out
    return run


def appa_law(name, formula, tol=TOL, method="numbers"):
    def deco(fn):
        LAWS.append(dict(group="appa", ref="App. A", name=name, formula=formula, method=method,
                         planes="none", samples=2000, tol=tol, negative=False, fn=_appa(fn)))
        return fn
    return deco


@appa_law("Tangent", "tan(α, φ) = sin(α, φ)/cos(α, φ) = sin α/sin φ")
def a_tan(al, be, ph, th):
    if not _safe(math.sin(al + ph), math.sin(ph)):
        return None
    return err_num(gsin(al, ph) / gcos(al, ph), math.sin(al) / math.sin(ph))


@appa_law("Reciprocity", "sin(α, φ) = cos(φ, α)")
def a_recip(al, be, ph, th):
    if not _safe(math.sin(al + ph)):
        return None
    return err_num(gsin(al, ph), gcos(ph, al))


@appa_law("The right-angle case", "α + φ = π/2:  sin(α, φ) = sin α,  cos(α, φ) = cos α")
def a_special(al, be, ph, th):
    return max(err_num(gsin(al, math.pi / 2 - al), math.sin(al)), err_num(gcos(al, math.pi / 2 - al), math.cos(al)))


@appa_law("The isosceles case", "sin(θ, θ) = cos(θ, θ) = 1/(2 cos θ)")
def a_iso(al, be, ph, th):
    if not _safe(math.cos(th), math.sin(2 * th)):
        return None
    return max(err_num(gsin(th, th), 0.5 / math.cos(th)), err_num(gcos(th, th), 0.5 / math.cos(th)))


@appa_law("π-addition", "sin(α + pπ, φ + qπ) = (−1)^q sin(α, φ),  cos(α + pπ, φ + qπ) = (−1)^p cos(α, φ)")
def a_pi(al, be, ph, th):
    if not _safe(math.sin(al + ph)):
        return None
    worst = 0.0
    for p_ in range(-2, 3):
        for q_ in range(-2, 3):
            worst = max(worst, err_num(gsin(al + p_ * math.pi, ph + q_ * math.pi), (-1) ** q_ * gsin(al, ph)),
                        err_num(gcos(al + p_ * math.pi, ph + q_ * math.pi), (-1) ** p_ * gcos(al, ph)))
    return worst


@appa_law("Studied-angle addition, four lines", "sin(α+β, φ) = sin(α,φ) + cos(α,φ) sin(β,φ+α), and the three other lines")
def a_add(al, be, ph, th):
    if not _safe(math.sin(al + ph), math.sin(be + ph), math.sin(al + be + ph), math.sin(al + ph + be), math.sin(be + ph + al)):
        return None
    s, c = gsin, gcos
    return max(err_num(s(al + be, ph), s(al, ph) + c(al, ph) * s(be, ph) * c(al, ph + be)),
               err_num(s(al + be, ph), s(al, ph) + c(al, ph) * s(be, ph + al)),
               err_num(c(al + be, ph), c(al, ph) * c(be, ph) + c(al, ph) * s(be, ph) * s(al, ph + be)),
               err_num(c(al + be, ph), c(al, ph) * c(be, ph) + s(al, ph) * s(be, ph) * c(al + be, ph)))


@appa_law("Projection-angle addition", "sin(α, φ+θ) = sin(α,φ) sin(α+φ, θ),  cos(α, φ+θ) = cos(α,φ) + sin(α,φ) cos(α+φ, θ)")
def a_proj(al, be, ph, th):
    if not _safe(math.sin(al + ph), math.sin(al + ph + th)):
        return None
    return max(err_num(gsin(al, ph + th), gsin(al, ph) * gsin(al + ph, th)),
               err_num(gcos(al, ph + th), gcos(al, ph) + gsin(al, ph) * gcos(al + ph, th)))


@appa_law("A consequence of the addition formulas", "cos(α+β, φ) = cos(α,φ) cos(β,φ) / (1 − sin(α,φ) sin(β,φ))")
def a_conseq(al, be, ph, th):
    if not _safe(math.sin(al + ph), math.sin(be + ph), math.sin(al + be + ph), 1 - gsin(al, ph) * gsin(be, ph) if _safe(math.sin(al + ph), math.sin(be + ph)) else 0):
        return None
    return err_num(gcos(al + be, ph), gcos(al, ph) * gcos(be, ph) / (1 - gsin(al, ph) * gsin(be, ph)))


@appa_law("Double angle", "sin(2α, φ) = sin(α,φ)(cos(2α,φ) + 1),  cos(2α, φ) = cos(α,φ)²/(1 − sin(α,φ)²)")
def a_double(al, be, ph, th):
    if not _safe(math.sin(al + ph), math.sin(2 * al + ph)):
        return None
    if not _safe(1 - gsin(al, ph) ** 2):
        return None
    return max(err_num(gsin(2 * al, ph), gsin(al, ph) * (gcos(2 * al, ph) + 1)),
               err_num(gcos(2 * al, ph), gcos(al, ph) ** 2 + gsin(al, ph) ** 2 * gcos(2 * al, ph)),
               err_num(gcos(2 * al, ph), gcos(al, ph) ** 2 / (1 - gsin(al, ph) ** 2)))


@appa_law("Partial derivatives", "sin(x,y)′ₓ = cos(x,y)/sin(x+y), sin(x,y)′ᵧ = −sin(x,y) cot(x+y), and the two for cos",
          tol=TOL_DIFF, method="finite differences")
def a_partial(al, be, ph, th):
    x, y, h = al, ph, 1e-5
    if not _safe(math.sin(x + y)) or abs(math.sin(x + y)) < 0.2:
        return None
    sx = (gsin(x + h, y) - gsin(x - h, y)) / (2 * h)
    sy = (gsin(x, y + h) - gsin(x, y - h)) / (2 * h)
    cx = (gcos(x + h, y) - gcos(x - h, y)) / (2 * h)
    cy = (gcos(x, y + h) - gcos(x, y - h)) / (2 * h)
    cot = math.cos(x + y) / math.sin(x + y)
    e = max(err_num(sx, gcos(x, y) / math.sin(x + y)), err_num(sy, -gsin(x, y) * cot),
            err_num(cx, -gcos(x, y) * cot), err_num(cy, gsin(x, y) / math.sin(x + y)))
    return e


# ---------------------------------------------------------------- 2.13 the complex case
group("complex", f"The complex case: each law at {PB} = π/2 against Appendix C",
      "In the complex plane the Axiom gives |" + GU + "| = 1. Each general formula is evaluated there and compared with the formula of Appendix C. "
      "The lines for 4.9.7 and F.8 compare with Python's own complex numbers.")


def claw(name, formula, samples=SAMPLES):
    def deco(fn):
        LAWS.append(dict(group="complex", ref="App. C", name=name, formula=formula, method="general formula against Appendix C",
                         planes="complex", samples=samples, tol=TOL, negative=False, fn=fn))
        return fn
    return deco


@claw("The Axiom gives the length 1", "|cos(π/2) + sin(π/2)| = 1", samples=1)
def c_length(P, rng):
    return err_num(P.L, 1.0)


@claw("Product", "(ac − bd) + (ad + bc) i")
def c_product(P, rng):
    a, b, c, d = (rng.uniform(-4, 4) for _ in range(4))
    return err_gn(GN(a, b, P) * GN(c, d, P), GN(a * c - b * d, a * d + b * c, P))


@claw("Conjugate, and the skew flip is the conjugate here", "z̄ = a − bi")
def c_conj(P, rng):
    a, b = rng.uniform(-4, 4), rng.uniform(-4, 4)
    return max(err_gn(GN(a, b, P).conj(), GN(a, -b, P)), err_num(GN(a, -b, P).page(), GN(a, b, P).page().conjugate()))


@claw("A number times its conjugate", "z z̄ = a² + b²")
def c_norm(P, rng):
    a, b = rng.uniform(-4, 4), rng.uniform(-4, 4)
    z = GN(a, b, P)
    return err_gn(z * z.conj(), GN(a * a + b * b, 0.0, P))


@claw("Division", "(a + bi)/(c + di) = [(ac + bd) + (bc − ad) i]/(c² + d²)")
def c_div(P, rng):
    a, b, c, d = (rng.uniform(-4, 4) for _ in range(4))
    if c * c + d * d < 0.01:
        return 0.0
    L, c_ = P.L, P.c
    den = c * c + 2 * c * d * L * c_ + d * d * L * L
    q = GN((a * c + 2 * a * d * L * c_ + b * d * L * L) / den, (b * c - a * d) / den, P)
    return err_gn(q, GN((a * c + b * d) / (c * c + d * d), (b * c - a * d) / (c * c + d * d), P))


@claw("Modulus, and the two-squares identity", "|z| = √(a² + b²);  (a² + b²)(c² + d²) = (ac − bd)² + (ad + bc)²")
def c_mod(P, rng):
    a, b, c, d = (rng.uniform(-4, 4) for _ in range(4))
    z, w = GN(a, b, P), GN(c, d, P)
    return max(err_num(math.sqrt(z.norm2()), math.hypot(a, b)),
               err_num(z.norm2() * w.norm2(), (a * c - b * d) ** 2 + (a * d + b * c) ** 2))


@claw("The Composer Unit, its direction and the Turning Unit are all i", f"{GU} = u = {JU} = i", samples=1)
def c_units(P, rng):
    return max(err_num(G_of(P).page(), 1j), err_num(u_of(P).page(), 1j), err_num(J_of(P).page(), 1j),
               err_gn(J_of(P), GN(0.0, 1.0, P)))


@claw("Shadow and height are the real and imaginary parts", "Sh(z) = a,  Ht(z) = b")
def c_shht(P, rng):
    a, b = rng.uniform(-4, 4), rng.uniform(-4, 4)
    z = GN(a, b, P)
    return max(err_num(z.shadow(), a), err_num(z.height(), b))


@claw("Polar form and Euler's formula", "C_θ + S_θ u = cos θ + i sin θ = e^(iθ)")
def c_euler(P, rng):
    t = rangle(rng, -6, 6)
    return max(err_gn(unit(t, P), GN(math.cos(t), math.sin(t), P)),
               err_gn(gn_exp(J_of(P) * t), GN(math.cos(t), math.sin(t), P)))


@claw("Euler's identity", "e^(iπ) + 1 = 0", samples=1)
def c_eid(P, rng):
    return err_gn(gn_exp(J_of(P) * math.pi) + 1, GN(0.0, 0.0, P))


@claw("Addition formulas and Pythagoras", "cos(α+β) = cos α cos β − sin α sin β, sin(α+β) = sin α cos β + cos α sin β, cos² + sin² = 1")
def c_add(P, rng):
    al, be = rangle(rng), rangle(rng)
    return max(err_num(C(al, P) * C(be, P) - S(al, P) * S(be, P), math.cos(al + be)),
               err_num(C(al, P) * S(be, P) + S(al, P) * C(be, P) + 2 * P.c * S(al, P) * S(be, P), math.sin(al + be)),
               err_num(C(al, P) ** 2 + 2 * P.c * C(al, P) * S(al, P) + S(al, P) ** 2, math.cos(al) ** 2 + math.sin(al) ** 2))


@claw("De Moivre's formula", "(cos θ + i sin θ)ⁿ = cos nθ + i sin nθ", samples=20)
def c_dm(P, rng):
    t = rangle(rng)
    return max(err_gn(unit(t, P).power(n), GN(math.cos(n * t), math.sin(n * t), P)) for n in range(-6, 11))


@claw("The powers of i and the step rule", "iⁿ = cos(nπ/2) + i sin(nπ/2),  iⁿ⁺¹ = −iⁿ⁻¹", samples=1)
def c_ipow(P, rng):
    G = G_of(P)
    return max(max(err_gn(G.power(n), GN(math.cos(n * math.pi / 2), math.sin(n * math.pi / 2), P)) for n in range(-4, 13)),
               max(err_gn(G.power(n + 1), -G.power(n - 1)) for n in range(-3, 12)))


@claw("Roots of unity and the square root of i", "e^(2πik/n);  √i = (1 + i)/√2", samples=1)
def c_roots(P, rng):
    worst = 0.0
    for n in range(1, 10):
        for k in range(n):
            t = 2 * math.pi * k / n
            worst = max(worst, err_gn(unit(t, P), GN(math.cos(t), math.sin(t), P)))
    k_ = math.sqrt(P.L) / (2 * math.cos(P.p / 2))
    return max(worst, err_gn(GN(k_, k_ / P.L, P), GN(1 / math.sqrt(2), 1 / math.sqrt(2), P)))


@claw("Exponential and logarithm", "e^(a + bi) = e^(a) (cos b + i sin b);  ln z = ln r + iθ")
def c_explog(P, rng):
    a, b = rng.uniform(-2, 2), rng.uniform(-2, 2)
    z = GN(a, b, P)
    f = unit(z.height(), P) * math.exp(z.shadow())
    th, M = cmath.phase(z.page()), abs(z.page())
    lg = GN(math.log(M) - th * P.c / P.s, th / (P.L * P.s), P)
    return max(err_gn(f, GN(math.exp(a) * math.cos(b), math.exp(a) * math.sin(b), P)), err_gn(lg, GN(math.log(M), th, P)))


@claw("The matrix", "[[a, −b], [b, a]]")
def c_mat(P, rng):
    a, b = rng.uniform(-4, 4), rng.uniform(-4, 4)
    return _merr(_mat(GN(a, b, P)), [[a, -b], [b, a]])


# the laws of Subsection 4.9.7 at pi/2, against the complex numbers that Python computes by itself
def _cz(rng):
    a, b = rng.uniform(-4, 4), rng.uniform(-4, 4)
    while a * a + b * b < 0.01:
        a, b = rng.uniform(-4, 4), rng.uniform(-4, 4)
    return a, b


@claw("4.9.7: the two pairs are one", "Re(z) = Sh(z) = Re z,  Com(z) = i Im z,  Ht(z) = Im z")
def c_pairs(P, rng):
    a, b = _cz(rng)
    z, w = GN(a, b, P), complex(a, b)              # w is the complex number a + bi, made by Python
    pa, pb_, sh, ht = parts_page(z.page(), P)
    return max(err_num(z.a, w.real), err_num(z.shadow(), w.real), err_num(z.b * P.Gpage, 1j * w.imag),
               err_num(z.height(), w.imag), err_num(pa, w.real), err_num(pb_, w.imag), err_num(sh, w.real), err_num(ht, w.imag))


@claw("4.9.7: the parts read from z and z̄ (Prop. 4.54)", "Re z = (z + z̄)/2,  i Im z = (z − z̄)/2;  z̄ = −z for imaginary numbers")
def c_fromconj(P, rng):
    a, b = _cz(rng)
    w, G = complex(a, b), P.Gpage
    Gb = G.conjugate()
    re_ = (G * w.conjugate() - Gb * w) / (G - Gb)
    com = G * (w - w.conjugate()) / (G - Gb)
    t = rng.uniform(-3, 3)
    return max(err_num(re_, (w + w.conjugate()) / 2), err_num(com, (w - w.conjugate()) / 2),
               err_num(re_, w.real), err_num(com, 1j * w.imag),
               err_num(Gb / G * complex(0, t), -complex(0, t)))      # the test of Prop. 4.54 (4): z̄ = (𝒢̄/𝒢) z


@claw("4.9.7: the parts of a product and of a square (Thm. 4.55)",
      "Re(zw) = Re z Re w − Im z Im w,  Im(zw) = Re z Im w + Re w Im z;  Re(z²) = a² − b²,  Im(z²) = 2ab")
def c_prod_parts(P, rng):
    (a, b), (c, d) = _cz(rng), _cz(rng)
    z, w = GN(a, b, P), GN(c, d, P)
    zw, sq = complex(a, b) * complex(c, d), complex(a, b) ** 2
    return max(err_num(z.a * w.a - P.L ** 2 * z.b * w.b, zw.real), err_num(z.shadow() * w.b + w.shadow() * z.b, zw.imag),
               err_num(z.shadow() * w.shadow() - z.height() * w.height(), zw.real),
               err_num(z.shadow() * w.height() + w.shadow() * z.height(), zw.imag),
               err_num(z.a ** 2 - P.L ** 2 * z.b ** 2, sq.real), err_num(2 * z.shadow() * z.b, sq.imag),
               err_num(z.shadow() ** 2 - z.height() ** 2, sq.real), err_num(2 * z.shadow() * z.height(), sq.imag))


@claw("4.9.7: the conjugate, the inverse and the quotient (Prop. 4.56)",
      "Re z̄ = Re z,  Im z̄ = −Im z;  Re(1/z) = a/(a² + b²),  Im(1/z) = −b/(a² + b²);  the parts of z/w")
def c_conj_parts(P, rng):
    (a, b), (c, d) = _cz(rng), _cz(rng)
    z, w = GN(a, b, P), GN(c, d, P)
    zc, zi, q = complex(a, b).conjugate(), 1 / complex(a, b), complex(a, b) / complex(c, d)
    n2, m2 = z.norm2(), w.norm2()
    s1, s2, h1, h2 = z.shadow(), w.shadow(), z.height(), w.height()
    return max(err_num(2 * s1 - z.a, zc.real), err_num(-z.b, zc.imag),
               err_num((2 * s1 - z.a) / n2, zi.real), err_num(-z.b / n2, zi.imag),
               err_num(s1 / n2, zi.real), err_num(-h1 / n2, zi.imag),
               err_num((z.a * w.conj().a + P.L ** 2 * z.b * w.b) / m2, q.real), err_num((s2 * z.b - s1 * w.b) / m2, q.imag),
               err_num((s1 * s2 + h1 * h2) / m2, q.real), err_num((s2 * h1 - s1 * h2) / m2, q.imag))


@claw("4.9.7: the two turns (Cor. 4.57)", f"Re(iz) = −Im z,  Im(iz) = Re z, for {JU}z and for {GU}z")
def c_turns(P, rng):
    a, b = _cz(rng)
    z, iz = GN(a, b, P), 1j * complex(a, b)
    s, h = z.shadow(), z.height()
    return max(err_num(-h - s * P.c / P.s, iz.real), err_num(s / (P.L * P.s), iz.imag), err_num(-h, iz.real), err_num(s, iz.imag),
               err_num(-P.L ** 2 * z.b, iz.real), err_num(z.conj().a, iz.imag),
               err_num(P.L * (s * P.c - h * P.s), iz.real), err_num(P.L * (s * P.s + h * P.c), iz.imag))


@claw("4.9.7: the signed area and the dot product (Prop. 4.58)", "z̄w = (ac + bd) + (ad − bc) i")
def c_area(P, rng):
    (a, b), (c, d) = _cz(rng), _cz(rng)
    z, w = GN(a, b, P), GN(c, d, P)
    v = complex(a, b).conjugate() * complex(c, d)
    return max(err_num(a * d - b * c, v.imag), err_num(z.shadow() * w.height() - z.height() * w.shadow(), v.imag),
               err_num((a * d - b * c) * P.L * P.s, v.imag), err_num(z.conj().a * w.a + P.L ** 2 * z.b * w.b, v.real),
               err_num(z.shadow() * w.shadow() + z.height() * w.height(), v.real))


@claw("4.9.7: the two trigonometries are one (Cor. 4.59)",
      "Re(zⁿ) = rⁿ cos nα,  Im(zⁿ) = rⁿ sin nα;  Re(e^(z)) = e^(a) cos b,  Im(e^(z)) = e^(a) sin b", samples=20)
def c_trig(P, rng):
    a, b = rng.uniform(-1.5, 1.5), rng.uniform(-1.5, 1.5)
    w = complex(a, b)
    r, al = abs(w), cmath.phase(w)
    worst = 0.0
    for n in range(-4, 8):
        wn = w ** n
        worst = max(worst, err_num(r ** n * C(n * al, P), wn.real), err_num(r ** n * S(n * al, P) / P.L, wn.imag))
    e = cmath.exp(w)
    return max(worst, err_num(math.exp(a) * C(b, P), e.real), err_num(math.exp(a) * S(b, P) / P.L, e.imag))


@claw("4.9.7: the size of the parts (Rem. 4.60)", f"|z|/sin {PB} = |z|:  |Re z| ≤ |z|,  |Im z| ≤ |z|")
def c_size(P, rng):
    a, b = _cz(rng)
    m = abs(complex(a, b))
    return max(err_num(m / P.s, m), max(0.0, abs(a) - m, abs(b) - m) / max(1.0, m))


@claw("App. F.8: the table of six numbers", f"1, {GU}, {GU}̄, {GU}², 1/{GU} and {JU} are 1, i, −i, −1, −i and i", samples=1)
def c_table(P, rng):
    G = G_of(P)
    nums = [GN(1.0, 0.0, P), G, G.conj(), G * G, G.inv(), J_of(P)]
    want = [1, 1j, -1j, -1, -1j, 1j]
    return max(err_num(z.page(), w) for z, w in zip(nums, want))


@claw("The habits of the complex plane hold here",
      "the wrong laws of the negative controls below are all right at π/2: there the parts and the readings are one pair")
def c_habits(P, rng):
    return max(evaluate(fn, P, rng) for fn in (n_p_com, n_p_reconj, n_p_mix, n_p_comht, n_p_bound, n_p_J, n_p_pow, n_p_inv))


# ---------------------------------------------------------------- 2.14 negative controls
group("neg", "Negative controls: wrong laws must fail",
      f"Each line is a law that is wrong where the Composer axis leans. It is tried in every plane that leans by more than "
      f"about 3° from the upright (|cos {PB}| > 0.05), and in each of them it must miss by more than {NEGATIVE_MIN:g}. "
      f"The table gives the smallest of these misses. If a wrong law passed, the test above could not see an error. "
      f"The last eight are habits of the complex plane used for the parts: right at π/2, wrong where the axis leans.")


def nlaw(name, formula):
    def deco(fn):
        LAWS.append(dict(group="neg", ref="control", name=name, formula=formula, method="must fail",
                         planes="leaning", samples=SAMPLES, tol=TOL, negative=True, fn=fn))
        return fn
    return deco


@nlaw("The skew flip as the conjugate", f"a − b{GU} would be the mirror image of a + b{GU}")
def n_skew(P, rng):
    z = rgn(P, rng)
    return err_num(GN(z.a, -z.b, P).page(), z.page().conjugate())


@nlaw("The addition formula without its extra term", "S_{α+β} = C_α S_β + S_α C_β")
def n_add(P, rng):
    al, be = rangle(rng), rangle(rng)
    return err_num(S(al + be, P), C(al, P) * S(be, P) + S(al, P) * C(be, P))


@nlaw("Euler's formula with the ordinary sine and cosine", "e^(θ" + JU + ") = cos θ + sin θ u")
def n_euler(P, rng):
    t = rangle(rng)
    return err_num(GN(math.cos(t), math.sin(t) / P.L, P).page(), cmath.exp(1j * t))


@nlaw("De Moivre with the ordinary sine and cosine", "(cos θ + sin θ u)ⁿ = cos nθ + sin nθ u, n = 2 … 5")
def n_dm(P, rng):
    t = rangle(rng)
    w = GN(math.cos(t), math.sin(t) / P.L, P)
    return max(err_gn(w.power(n), GN(math.cos(n * t), math.sin(n * t) / P.L, P)) for n in range(2, 6))


@nlaw("The square of the Composer Unit as a real number", f"{GU}² = −|{GU}|²")
def n_square(P, rng):
    G = G_of(P)
    return err_gn(G * G, GN(-P.L ** 2, 0.0, P))


@nlaw("Pythagoras in the leaning frame", "|z|² = a² + b²|" + GU + "|²")
def n_pyth(P, rng):
    z = rgn(P, rng)
    return err_num(z.a ** 2 + z.b ** 2 * P.L ** 2, abs(z.page()) ** 2)


# the habits of the complex plane, used for the parts: right at pi/2, wrong where the axis leans
@nlaw("The Composer part of a product with the real parts in place of the shadows", "Com(z₁z₂) = Re(z₁)Com(z₂) + Re(z₂)Com(z₁)")
def n_p_com(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    return err_num(z1.a * z2.b + z2.a * z1.b, from_page(z1.page() * z2.page(), P).b)


@nlaw("The real part of the conjugate kept, as in the complex plane", "Re(z̄) = Re(z)")
def n_p_reconj(P, rng):
    z = rgn(P, rng)
    return err_num(z.a, from_page(z.page().conjugate(), P).a)


@nlaw("The real part of a product from the shadows and the heights", "Re(z₁z₂) = Sh(z₁)Sh(z₂) − Ht(z₁)Ht(z₂)")
def n_p_mix(P, rng):
    z1, z2 = rgn(P, rng), rgn(P, rng)
    return err_num(z1.shadow() * z2.shadow() - z1.height() * z2.height(), from_page(z1.page() * z2.page(), P).a)


@nlaw("The Composer part as the height along the Composer axis", f"Com(z) = Ht(z) u, without the factor 1/sin {PB}")
def n_p_comht(P, rng):
    z = GN(rng.uniform(-4, 4), rng.choice((-1, 1)) * rng.uniform(1, 4), P)
    return err_num(z.height() / P.L, from_page(z.page(), P).b)


@nlaw("The parts never longer than the number, as in the complex plane", "|Re(z)| ≤ |z|")
def n_p_bound(P, rng):
    w = cmath.exp(1j * (P.p - math.pi / 2)) * rng.uniform(1, 3)    # a number at a right angle to the Composer axis
    return max(0.0, abs(from_page(w, P).a) - abs(w)) / max(1.0, abs(w))


@nlaw(f"The quarter turn of the parts without its cot term", f"Re({JU}z) = −Ht(z)")
def n_p_J(P, rng):
    z = rgn(P, rng)
    return err_num(-z.height(), from_page(1j * z.page(), P).a)


@nlaw("The powers of the parts with the ordinary cosine", "Re(zⁿ) = |z|ⁿ cos nα, n = 2 … 5")
def n_p_pow(P, rng):
    z = rgn(P, rng, -1.5, 1.5)
    w = z.page()
    r, al = abs(w), cmath.phase(w)
    return max(err_num(r ** n * math.cos(n * al), from_page(w ** n, P).a) for n in range(2, 6))


@nlaw("The real part of the inverse, as in the complex plane", "Re(1/z) = Re(z)/|z|²")
def n_p_inv(P, rng):
    z = rgn(P, rng)
    return err_num(z.a / z.norm2(), from_page(1 / z.page(), P).a)


# =====================================================================================
# 3. Worked examples of Section 4.9: the printed numbers
# =====================================================================================
def worked_examples():
    P = Plane(math.radians(79.86))
    c, s, L = P.c, P.s, P.L
    G = G_of(P)
    z = GN(1.0, 1.0, P)
    zb = z.conj()
    inv = z.inv()
    num_a = 2 + 2 * 2 * L * c - L * L
    q = GN(2.0, -1.0, P) * z.inv()
    skew = GN(1.0, -1.0, P)
    J = J_of(P)
    Jz = J * z
    al = cmath.phase(z.page())
    ap = al + P.p
    prod = z * G
    G2, G3 = G * G, G * G * G
    k = math.sqrt(L) / (2 * math.cos(P.p / 2))
    sq = z * z                                    # Subsection 4.9.7: the square of 1 + G
    Gz = G * z
    top = zb * GN(2.0, -1.0, P)                   # the top of (2 - G)/(1 + G)
    w1, w2 = z.page(), GN(2.0, -1.0, P).page()
    cross = w1.real * w2.imag - w1.imag * w2.real # the signed area on the page
    deg = math.degrees
    rows = [
        # (where, what, printed, computed, decimals)
        ("4.9.2 example", "cos 79.86°", 0.17605, c, 5),
        ("4.9.2 example", "sin 79.86°", 0.98438, s, 5),
        ("4.9.2 example", f"|{GU}| = |cos {PB} + sin {PB}|", 1.16043, L, 5),
        ("4.9.2 example", f"|{GU}| cos {PB}, one step across", 0.20430, L * c, 5),
        ("4.9.2 example", f"|{GU}| sin {PB}, one step up", 1.14231, L * s, 5),
        ("4.9.2 example", f"|{GU}|²", 1.34661, L * L, 5),
        ("4.9.2 example", f"2|{GU}| cos {PB}", 0.40860, 2 * L * c, 5),
        ("4.9.2 example", f"tip of 1 + {GU}, across (its shadow)", 1.20430, z.page().real, 5),
        ("4.9.2 example", f"real part of the conjugate of 1 + {GU}", 1.40860, zb.a, 5),
        ("4.9.2 example", "tip of the conjugate, across", 1.20430, zb.page().real, 5),
        ("4.9.2 example", "tip of the conjugate, up", -1.14231, zb.page().imag, 5),
        ("4.9.2 example", f"|1 + {GU}|²", 2.75521, z.norm2(), 5),
        ("4.9.2 example", f"|1 + {GU}|", 1.65988, math.sqrt(z.norm2()), 5),
        ("4.9.2 example", f"1/(1 + {GU}), real part", 0.51125, inv.a, 5),
        ("4.9.2 example", f"1/(1 + {GU}), Composer coordinate", -0.36295, inv.b, 5),
        ("4.9.2 example", f"(2 − {GU})/(1 + {GU}), real part of the top", 1.47059, num_a, 5),
        ("4.9.2 example", f"(2 − {GU})/(1 + {GU}), real part", 0.53375, q.a, 5),
        ("4.9.2 example", f"(2 − {GU})/(1 + {GU}), Composer coordinate", -1.08885, q.b, 5),
        ("4.9.2 example", f"|1 − {GU}|", 1.39212, math.sqrt(skew.norm2()), 5),
        ("4.9.2 example", f"(1 + {GU})(1 − {GU}), real part", 2.34661, (z * skew).a, 5),
        ("4.9.2 example", f"(1 + {GU})(1 − {GU}), Composer coordinate", -0.40860, (z * skew).b, 5),
        ("4.9.3 example", f"cot {PB}", 0.17885, c / s, 5),
        ("4.9.3 example", f"1/(|{GU}| sin {PB})", 0.87542, 1 / (L * s), 5),
        ("4.9.3 example", f"{JU}, real part", -0.17885, J.a, 5),
        ("4.9.3 example", f"{JU}, Composer coordinate", 0.87542, J.b, 5),
        ("4.9.3 example", "0.17885² (exactly cot²)", 0.03199, (c / s) ** 2, 5),
        ("4.9.3 example", "1.34661 × 0.87542² (exactly 1/sin²)", 1.03199, 1 / s ** 2, 5),
        ("4.9.3 example", "2 × 0.17885 × 0.87542 (exact value)", 0.31313, 2 * (c / s) / (L * s), 5),
        ("4.9.3 example", "0.40860 × 0.87542² (exact value)", 0.31313, 2 * L * c / (L * s) ** 2, 5),
        ("4.9.3 example", f"{JU}², real part", -1.0, (J * J).a, 5),
        ("4.9.3 example", f"{JU}², Composer coordinate", 0.0, (J * J).b, 5),
        ("4.9.3 example", f"1 + {GU} = Sh + Ht {JU}: Sh", 1.20430, z.shadow(), 5),
        ("4.9.3 example", f"1 + {GU} = Sh + Ht {JU}: Ht", 1.14231, z.height(), 5),
        ("4.9.3 example", f"{JU}(1 + {GU}), real part", -1.35770, Jz.a, 5),
        ("4.9.3 example", f"{JU}(1 + {GU}), Composer coordinate", 1.05427, Jz.b, 5),
        ("4.9.3 example", f"tip of {JU}(1 + {GU}), across", -1.14231, Jz.page().real, 5),
        ("4.9.3 example", f"tip of {JU}(1 + {GU}), up", 1.20430, Jz.page().imag, 5),
        ("4.9.4 example", f"(1 + {GU}){GU}, real part", -1.34661, prod.a, 5),
        ("4.9.4 example", f"(1 + {GU}){GU}, Composer coordinate", 1.40860, prod.b, 5),
        ("4.9.4 example", f"angle α of 1 + {GU}, degrees", 43.487, deg(al), 3),
        ("4.9.4 example", "C_α", 0.60245, C(al, P), 5),
        ("4.9.4 example", "S_α", 0.69911, S(al, P), 5),
        ("4.9.4 example", "|z| C_α (the real part 1)", 1.0, math.sqrt(z.norm2()) * C(al, P), 5),
        ("4.9.4 example", f"|z| S_α (the length {GU}, 1.16043)", 1.16043, math.sqrt(z.norm2()) * S(al, P), 5),
        ("4.9.4 example", f"C_β = cos({PB}, 0)", 0.0, gcos(P.p, 0.0), 5),
        ("4.9.4 example", f"S_β = sin({PB}, 0)", 1.0, gsin(P.p, 0.0), 5),
        ("4.9.4 example", "α + β, degrees", 123.347, deg(ap), 3),
        ("4.9.4 example", "C_{α+β}", -0.69911, C(ap, P), 5),
        ("4.9.4 example", f"2 cos {PB}", 0.35211, 2 * c, 5),
        ("4.9.4 example", "S_{α+β}", 0.84861, S(ap, P), 5),
        ("4.9.4 example", "magnitude of the product", 1.92618, math.sqrt(z.norm2()) * L, 5),
        ("4.9.4 example", "real part, by turning", -1.34661, math.sqrt(z.norm2()) * L * C(ap, P), 5),
        ("4.9.4 example", "Composer coordinate, by turning", 1.40860, math.sqrt(z.norm2()) * L * S(ap, P) / L, 5),
        ("4.9.4 example", f"extra term 2 cos {PB} S_α S_β", 0.24616, 2 * c * S(al, P), 5),
        ("4.9.5 example", f"{GU}², real part", -1.34661, G2.a, 5),
        ("4.9.5 example", f"{GU}², Composer coordinate", 0.40860, G2.b, 5),
        ("4.9.5 example", f"{GU}³, real part", -0.55022, G3.a, 5),
        ("4.9.5 example", f"{GU}³, Composer coordinate", -1.17966, G3.b, 5),
        ("4.9.5 example", f"|{GU}|³", 1.56265, L ** 3, 5),
        ("4.9.5 example", f"sin 3{PB}/sin {PB}", -0.87602, math.sin(3 * P.p) / s, 5),
        ("4.9.5 example", f"real part by the formula, −|{GU}|³ × 2 cos {PB}", -0.55022, -L ** 3 * 2 * c, 5),
        ("4.9.5 example", f"Composer coordinate by the formula, |{GU}|² sin 3{PB}/sin {PB}", -1.17966, L ** 2 * math.sin(3 * P.p) / s, 5),
        ("4.9.5 example", f"3 × 79.86°", 239.58, 3 * 79.86, 2),
        ("4.9.5 example", f"tip of {GU}³, across", -0.79122, G3.page().real, 5),
        ("4.9.5 example", f"tip of {GU}³, down", -1.34753, G3.page().imag, 5),
        ("4.9.6 example", f"√|{GU}|", 1.07723, math.sqrt(L), 5),
        ("4.9.6 example", f"{PB}/2, degrees", 39.93, 79.86 / 2, 2),
        ("4.9.6 example", f"cos({PB}/2)", 0.76683, math.cos(P.p / 2), 5),
        ("4.9.6 example", f"√{GU} = ±0.70240 (1 + u): the factor", 0.70240, k, 5),
        ("4.9.6 example", f"√{GU}, Composer coordinate", 0.60529, k / L, 5),
        ("4.9.6 example", "2 × 0.70240 × 0.60529 (exact value)", 0.85030, 2 * k * k / L, 5),
        ("4.9.6 example", "0.40860 × 0.60529² (exact value)", 0.14970, 2 * L * c * (k / L) ** 2, 5),
        ("4.9.6 example", f"(√{GU})², real part", 0.0, (GN(k, k / L, P) * GN(k, k / L, P)).a, 5),
        ("4.9.6 example", f"(√{GU})², Composer coordinate", 1.0, (GN(k, k / L, P) * GN(k, k / L, P)).b, 5),
        ("4.9.6 example", "sixth root w₁, real part", 0.34511, C(math.pi / 3, P), 5),
        ("4.9.6 example", "sixth roots w₁, w₂: Composer coordinate", 0.75814, S(math.pi / 3, P) / L, 5),
        ("4.9.6 example", "sixth root w₂, real part", -0.65489, C(2 * math.pi / 3, P), 5),
        ("4.9.6 example", "w₁ − w₂, real parts: exactly 1", 1.0, C(math.pi / 3, P) - C(2 * math.pi / 3, P), 5),
        ("4.9.6 example", f"sin 60°/(|{GU}| sin {PB})", 0.75814, math.sin(math.pi / 3) / (L * s), 5),
        ("4.9.7 example", f"Re(z̄) = 2 Sh(z) − Re(z), z = 1 + {GU}", 1.40860, 2 * z.shadow() - z.a, 5),
        ("4.9.7 example", "the same, read on the page from the mirror image", 1.40860, from_page(z.page().conjugate(), P).a, 5),
        ("4.9.7 example", "the shadow, halfway between 1 and 1.40860", 1.20430, (z.a + zb.a) / 2, 5),
        ("4.9.7 example", "Re(z²) = 1 − 1.34661", -0.34661, sq.a, 5),
        ("4.9.7 example", "Com(z²): the coordinate 2 Sh(z) = 2 × 1.20430", 2.40860, sq.b, 5),
        ("4.9.7 example", "the same, 2 Sh(z) computed", 2.40860, 2 * z.shadow(), 5),
        ("4.9.7 example", "Sh(z²) = 1.20430² − 1.14231² (exact value)", 0.14547, z.shadow() ** 2 - z.height() ** 2, 5),
        ("4.9.7 example", "Ht(z²) = 2 × 1.20430 × 1.14231 (exact value)", 2.7514, 2 * z.shadow() * z.height(), 4),
        ("4.9.7 example", "tip of z², across: −0.34661 + 2.40860 × 0.20430 (exact value)", 0.14547, sq.page().real, 5),
        ("4.9.7 example", "tip of z², up: 2.40860 × 1.14231 (exact value)", 2.7514, sq.page().imag, 4),
        ("4.9.7 example", "Re(1/z) = 1.40860/2.75521", 0.51125, (2 * z.shadow() - z.a) / z.norm2(), 5),
        ("4.9.7 example", "Com(1/z): −1/2.75521", -0.36295, -z.b / z.norm2(), 5),
        ("4.9.7 example", "Sh(1/z) = 1.20430/2.75521", 0.43710, (1 / z.page()).real, 5),
        ("4.9.7 example", "Ht(1/z) = −1.14231/2.75521", -0.41460, (1 / z.page()).imag, 5),
        ("4.9.7 example", f"{GU}z, real part −|{GU}|² b", -1.34661, Gz.a, 5),
        ("4.9.7 example", f"{GU}z, Composer coordinate Re(z̄)", 1.40860, Gz.b, 5),
        ("4.9.7 example", f"{JU}z, real part −1.14231 − 1.20430 × 0.17885 (exact value)", -1.35770, -z.height() - z.shadow() * c / s, 5),
        ("4.9.7 example", f"{JU}z, Composer coordinate 1.20430 × 0.87542 (exact value)", 1.05427, z.shadow() / (L * s), 5),
        ("4.9.7 example", f"z̄(2 − {GU}), real part", 1.47059, top.a, 5),
        ("4.9.7 example", f"z̄(2 − {GU}), Composer coordinate 1 × (−1) − 1 × 2", -3.0, top.b, 5),
        ("4.9.7 example", f"z̄(2 − {GU}), height −3 × 1.14231 = the signed area", -3.42693, top.height(), 5),
        ("4.9.7 example", "the signed area, from the two arrows on the page", -3.42693, cross, 5),
        ("4.9.7 example", f"{JU} = −0.17885 + 0.87542 {GU}: length of its Composer part, 1/sin {PB}", 1.0159, abs(J.b) * L, 4),
        ("Rem. 4.60", f"1/sin {PB} − 1, in percent", 1.6, (1 / s - 1) * 100, 1),
        ("Fig. 4.20", f"the base of the triangle, 2|{GU}| cos {PB}", 0.40860, zb.a - z.a, 5),
        ("Fig. 4.20", "the middle of the dashed leg, measured from Re(z) = 1", 0.30645, 1.5 * L * c, 5),
        ("Fig. 4.20", f"each leg of the triangle is |{GU}| long", 1.16043, abs(z.page() - zb.a), 5),
        ("Fig. 4.20", "tip of z², up (as drawn)", 2.75136, sq.page().imag, 5),
    ]
    out = []
    for where, what, printed, computed, d in rows:
        ok = abs(computed - printed) <= 0.5 * 10 ** (-d) + 1e-12
        out.append(dict(where=where, what=what, printed=printed, computed=computed, decimals=d, ok=ok))
    return out


# =====================================================================================
# 4. Exact checks, with no rounding at all
# =====================================================================================
class Q3:
    """An exact number p + q*sqrt(3), with p and q fractions."""

    __slots__ = ("p", "q")

    def __init__(self, p, q=0):
        self.p, self.q = Fraction(p), Fraction(q)

    @staticmethod
    def of(x):
        return x if isinstance(x, Q3) else Q3(x)

    def __add__(self, o):
        o = Q3.of(o)
        return Q3(self.p + o.p, self.q + o.q)

    __radd__ = __add__

    def __neg__(self):
        return Q3(-self.p, -self.q)

    def __sub__(self, o):
        return self + (-Q3.of(o))

    def __rsub__(self, o):
        return Q3.of(o) - self

    def __mul__(self, o):
        o = Q3.of(o)
        return Q3(self.p * o.p + 3 * self.q * o.q, self.p * o.q + self.q * o.p)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = Q3.of(o)
        n = o.p * o.p - 3 * o.q * o.q
        return self * Q3(o.p / n, -o.q / n)

    def __rtruediv__(self, o):
        return Q3.of(o) / self

    def __eq__(self, o):
        o = Q3.of(o)
        return self.p == o.p and self.q == o.q

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash((self.p, self.q))

    def is_zero(self):
        return self.p == 0 and self.q == 0


class EPlane:
    """A relative plane whose cos, sin and length are exact numbers p + q*sqrt(3)."""

    def __init__(self, key, label, c, s, L, c2, s2):
        self.key, self.label = key, label
        self.c, self.s, self.L, self.c2, self.s2 = c, s, L, c2, s2      # cos, sin, |G|, cos 2phibar, sin 2phibar
        self.G = EC(L * c, L * s)                                       # the Composer Unit as an arrow of the page


class EC:
    """An exact arrow of the page, x + iy, multiplied by Wessel's rule."""

    __slots__ = ("x", "y")

    def __init__(self, x, y=0):
        self.x, self.y = Q3.of(x), Q3.of(y)

    @staticmethod
    def of(v):
        return v if isinstance(v, EC) else EC(v, 0)

    def __add__(self, o):
        o = EC.of(o)
        return EC(self.x + o.x, self.y + o.y)

    def __sub__(self, o):
        o = EC.of(o)
        return EC(self.x - o.x, self.y - o.y)

    def __mul__(self, o):
        o = EC.of(o)
        return EC(self.x * o.x - self.y * o.y, self.x * o.y + self.y * o.x)

    def conj(self):
        return EC(self.x, -self.y)

    def __truediv__(self, o):
        o = EC.of(o)
        n = o.x * o.x + o.y * o.y
        return self * EC(o.x / n, -o.y / n)

    def __eq__(self, o):
        o = EC.of(o)
        return self.x == o.x and self.y == o.y


class EGN:
    """An exact General Number a + bG, computed in the leaning frame with the Rule of the Square."""

    __slots__ = ("a", "b", "P")

    def __init__(self, a, b, P):
        self.a, self.b, self.P = Q3.of(a), Q3.of(b), P

    def lift(self, o):
        return o if isinstance(o, EGN) else EGN(o, 0, self.P)

    def __add__(self, o):
        o = self.lift(o)
        return EGN(self.a + o.a, self.b + o.b, self.P)

    def __sub__(self, o):
        o = self.lift(o)
        return EGN(self.a - o.a, self.b - o.b, self.P)

    def __mul__(self, o):
        if not isinstance(o, EGN):
            return EGN(self.a * o, self.b * o, self.P)
        P = self.P
        return EGN(self.a * o.a - P.L * P.L * self.b * o.b, self.a * o.b + self.b * o.a + 2 * P.L * P.c * self.b * o.b, P)

    def conj(self):
        P = self.P
        return EGN(self.a + 2 * self.b * P.L * P.c, -self.b, P)

    def norm2(self):
        P = self.P
        return self.a * self.a + 2 * self.a * self.b * P.L * P.c + self.b * self.b * P.L * P.L

    def inv(self):
        n, zc = self.norm2(), self.conj()
        return EGN(zc.a / n, zc.b / n, self.P)

    def shadow(self):
        return self.a + self.b * self.P.L * self.P.c

    def height(self):
        return self.b * self.P.L * self.P.s

    def page(self):
        return EC(self.a) + self.P.G * EC(self.b)

    def __eq__(self, o):
        o = self.lift(o)
        return self.a == o.a and self.b == o.b


def e_from_page(w, P):
    """The exact pair (a, b) with a + bG = w: the leg along the Composer axis rises b|G| sin(phibar)."""
    b = w.y / (P.L * P.s)
    return EGN(w.x - b * P.L * P.c, b, P)


E90 = EPlane("π/2", "π/2, the complex plane", Q3(0), Q3(1), Q3(1), Q3(-1), Q3(0))
E60 = EPlane("π/3", "π/3", Q3(Fraction(1, 2)), Q3(0, Fraction(1, 2)), Q3(Fraction(1, 2), Fraction(1, 2)),
             Q3(Fraction(-1, 2)), Q3(0, Fraction(1, 2)))
EXACT_PLANES = {"π/2": E90, "π/3": E60}
EXACT_SAMPLES = 300


def rq(rnd):
    """A random exact number p + q*sqrt(3), with small fractions; q = 0 half of the time."""
    p = Fraction(rnd.randint(-9, 9), rnd.randint(1, 6))
    q = Fraction(rnd.randint(-9, 9), rnd.randint(1, 6)) if rnd.random() < 0.5 else Fraction(0)
    return Q3(p, q)


def regn(P, rnd):
    while True:
        z = EGN(rq(rnd), rq(rnd), P)
        if not (z.a.is_zero() and z.b.is_zero()):
            return z


XCHECKS = []


def xcheck(ref, name, planes=("π/2", "π/3")):
    def deco(fn):
        XCHECKS.append(dict(ref=ref, name=name, planes=planes, fn=fn))
        return fn
    return deco


@xcheck("Prop. 4.22, 4.23", "The product in the leaning frame is the product of the arrows on the page")
def x_prod(P, rnd):
    z, w = regn(P, rnd), regn(P, rnd)
    return (z * w).page() == z.page() * w.page()


@xcheck("Thm. 4.38 (1)", "The conjugate is the mirror image in the real axis")
def x_conj(P, rnd):
    z = regn(P, rnd)
    return z.conj().page() == z.page().conj()


@xcheck("Prop. 4.54 (1)", "The parts follow from the shadow and the height, and back")
def x_dict(P, rnd):
    z = regn(P, rnd)
    w = z.page()
    return (e_from_page(w, P) == z and w.x == z.shadow() and w.y == z.height()
            and z.a == w.x - w.y * P.c / P.s and z.b == w.y / (P.L * P.s))


@xcheck("Prop. 4.54 (3)", "The parts read from z and z̄, and z̄ = a + b" + GU + "̄")
def x_fromconj(P, rnd):
    z = regn(P, rnd)
    w, G = z.page(), P.G
    Gb = G.conj()
    return ((G * w.conj() - Gb * w) / (G - Gb) == EC(z.a) and G * (w - w.conj()) / (G - Gb) == G * EC(z.b)
            and z.conj() == EGN(z.a, 0, P) + EGN(0, 1, P).conj() * z.b
            and G - Gb == EC(0, 1) * EC(2 * P.L * P.s))


@xcheck("Thm. 4.55", "The parts, the shadow and the height of a product and of a square")
def x_prod_parts(P, rnd):
    z, w = regn(P, rnd), regn(P, rnd)
    pw = z.page() * w.page()
    p = e_from_page(pw, P)
    sq = e_from_page(z.page() * z.page(), P)
    s1, s2, h1, h2 = z.shadow(), w.shadow(), z.height(), w.height()
    return (p.a == z.a * w.a - P.L * P.L * z.b * w.b and p.b == s1 * w.b + s2 * z.b
            and pw.x == s1 * s2 - h1 * h2 and pw.y == s1 * h2 + s2 * h1
            and sq.a == z.a * z.a - P.L * P.L * z.b * z.b and sq.b == 2 * s1 * z.b)


@xcheck("Prop. 4.56", "The conjugate, the inverse and the quotient")
def x_cinvq(P, rnd):
    z, w = regn(P, rnd), regn(P, rnd)
    s1, s2, h1, h2 = z.shadow(), w.shadow(), z.height(), w.height()
    n1, n2 = z.norm2(), w.norm2()
    c_ = e_from_page(z.page().conj(), P)
    zi = EC(1) / z.page()
    i_ = e_from_page(zi, P)
    q = z.page() / w.page()
    q_ = e_from_page(q, P)
    return (c_.a == 2 * s1 - z.a and c_.b == -z.b and 2 * s1 == z.a + c_.a
            and i_.a == (2 * s1 - z.a) / n1 and i_.b == -z.b / n1 and zi.x == s1 / n1 and zi.y == -h1 / n1
            and q_.a == (z.a * w.conj().a + P.L * P.L * z.b * w.b) / n2 and q_.b == (s2 * z.b - s1 * w.b) / n2
            and q.x == (s1 * s2 + h1 * h2) / n2 and q.y == (s2 * h1 - s1 * h2) / n2)


@xcheck("Cor. 4.57", "The products by " + JU + " and by " + GU)
def x_turns(P, rnd):
    z = regn(P, rnd)
    s, h = z.shadow(), z.height()
    jz = e_from_page(EC(0, 1) * z.page(), P)
    gz = e_from_page(P.G * z.page(), P)
    J = EGN(-P.c / P.s, 1 / (P.L * P.s), P)
    return (jz.a == -h - s * P.c / P.s and jz.b == s / (P.L * P.s) and J * z == jz
            and gz.a == -P.L * P.L * z.b and gz.b == z.conj().a and EGN(0, 1, P) * z == gz)


@xcheck("Prop. 4.58", "The signed area and the dot product")
def x_area(P, rnd):
    z, w = regn(P, rnd), regn(P, rnd)
    w1, w2 = z.page(), w.page()
    v = w1.conj() * w2
    v_ = e_from_page(v, P)
    det = z.a * w.b - z.b * w.a
    return (v_.b == det and v.y == det * P.L * P.s and v.y == w1.x * w2.y - w1.y * w2.x
            and v.x == w1.x * w2.x + w1.y * w2.y and v_.a == z.conj().a * w.a + P.L * P.L * z.b * w.b)


@xcheck("App. F.8", "The magnitude in the two pairs, and the table of six numbers")
def x_table(P, rnd):
    z = regn(P, rnd)
    w = z.page()
    m2 = w.x * w.x + w.y * w.y
    if not (m2 == z.a * z.a + 2 * P.c * z.a * z.b * P.L + z.b * z.b * P.L * P.L and m2 == z.shadow() * z.shadow() + z.height() * z.height()):
        return False
    G, L, c, s = P.G, P.L, P.c, P.s
    six = [(EC(1), (1, 0, 1, 0)), (G, (0, 1, L * c, L * s)), (G.conj(), (2 * L * c, -1, L * c, -L * s)),
           (G * G, (-L * L, 2 * L * c, L * L * P.c2, L * L * P.s2)), (EC(1) / G, (2 * c / L, -1 / (L * L), c / L, -s / L)),
           (EC(0, 1), (-c / s, 1 / (L * s), 0, 1))]
    for v, want in six:
        e = e_from_page(v, P)
        if not all(Q3.of(g_) == Q3.of(w_) for g_, w_ in zip((e.a, e.b, v.x, v.y), want)):
            return False
    return True


@xcheck("App. C", "At π/2 the leaning frame is the complex plane: a + b" + GU + " is a + bi, and every operation agrees", planes=("π/2",))
def x_complex(P, rnd):
    z, w = regn(P, rnd), regn(P, rnd)
    zc, wc = EC(z.a, z.b), EC(w.a, w.b)            # the complex numbers a + bi and c + di
    return (z.page() == zc and (z * w).page() == zc * wc and z.conj().page() == zc.conj()
            and z.inv().page() == EC(1) / zc and (z * w.inv()).page() == zc / wc
            and z.shadow() == z.a and z.height() == z.b and EGN(-P.c / P.s, 1 / (P.L * P.s), P) == EGN(0, 1, P))


@xcheck("App. C", "At π/2 the laws of 4.9.7 are the laws of the real and imaginary parts", planes=("π/2",))
def x_complex_parts(P, rnd):
    z, w = regn(P, rnd), regn(P, rnd)
    zc, wc = EC(z.a, z.b), EC(w.a, w.b)
    p = zc * wc
    u = zc.conj() * wc
    return (p.x == zc.x * wc.x - zc.y * wc.y and p.y == zc.x * wc.y + wc.x * zc.y
            and (z * w).a == p.x and (z * w).b == p.y
            and z.conj().a == z.a and z.conj().b == -z.b
            and e_from_page(EC(0, 1) * zc, P) == EGN(-z.b, z.a, P) and EGN(0, 1, P) * z == EGN(-z.b, z.a, P)
            and (z.conj() * w).b == z.a * w.b - z.b * w.a and (z.conj() * w).a == z.a * w.a + z.b * w.b
            and u.y == z.a * w.b - z.b * w.a and u.x == z.a * w.a + z.b * w.b)


@xcheck("Controls", "At π/2 the habits of the complex plane are exact", planes=("π/2",))
def x_habits(P, rnd):
    z, w = regn(P, rnd), regn(P, rnd)
    return ((z * w).b == z.a * w.b + w.a * z.b and z.conj().a == z.a
            and (z * w).a == z.shadow() * w.shadow() - z.height() * w.height()
            and z.height() / P.L == z.b and e_from_page(EC(0, 1) * z.page(), P).a == -z.height()
            and z.inv().a == z.a / z.norm2())


@xcheck("Controls", "At π/3 the same habits are wrong, for every number with b ≠ 0", planes=("π/3",))
def x_habits_wrong(P, rnd):
    z = EGN(rq(rnd), Q3(rnd.choice((-1, 1)) * rnd.randint(1, 9)), P)
    w = EGN(rq(rnd), Q3(rnd.choice((-1, 1)) * rnd.randint(1, 9)), P)
    return ((z * w).b != z.a * w.b + w.a * z.b and z.conj().a != z.a and z.height() / P.L != z.b
            and z.inv().a != z.a / z.norm2())


def exact_checks():
    rows = []
    for X in XCHECKS:
        for key in X["planes"]:
            P = EXACT_PLANES[key]
            rnd = random.Random(f"{SEED}-exact-{X['name']}-{key}")
            ok, n = True, 0
            for _ in range(EXACT_SAMPLES):
                try:
                    good = bool(X["fn"](P, rnd))
                except ZeroDivisionError:
                    continue
                ok &= good
                n += 1
            rows.append(dict(ref=X["ref"], name=X["name"], plane=key, cases=n, ok=ok and n > 0))
    return rows


# =====================================================================================
# 5. Approaching the complex plane
# =====================================================================================
APPROACH_DEG = (79.86, 85.0, 89.0, 89.9, 89.99, 89.999, 90.0)


def approach():
    """The habits of the complex plane, and the true laws of 4.9.7, as the Composer axis rises to the upright."""
    habits = (n_p_com, n_p_reconj, n_p_mix, n_p_comht, n_p_bound, n_p_J, n_p_pow, n_p_inv)
    true_laws = (t_p_prod, t_p_prod_sh, t_p_conj, t_p_inv, t_p_quot, t_p_J, t_p_G, t_p_area, t_p_dot, t_p_pow, t_p_exp)
    rows = []
    for deg in APPROACH_DEG:
        P = Plane(math.pi / 2 if deg == 90.0 else math.radians(deg))
        rng = random.Random(f"{SEED}-approach")
        miss = max(evaluate(fn, P, rng) for fn in habits for _ in range(100))
        rng = random.Random(f"{SEED}-approach-laws")
        err = max(evaluate(fn, P, rng) for fn in true_laws for _ in range(30))
        rows.append(dict(deg=deg, cos=abs(P.c), miss=miss, err=err))
    ok = (all(r["err"] <= TOL for r in rows) and rows[-1]["miss"] <= 1e-12
          and all(rows[i + 1]["miss"] <= rows[i]["miss"] for i in range(len(rows) - 1)))
    return dict(rows=rows, ok=ok)


# =====================================================================================
# 6. Running everything
# =====================================================================================
def build_planes(rng):
    planes = [Plane(math.radians(79.86), label="79.86°, the book's angle", kind="book"),
              Plane(math.pi / 3, label="π/3, where cos, sin and |" + GU + "| are exact numbers", kind="axiom"),
              Plane(math.pi / 2, label="π/2, the complex plane", kind="complex")]
    while len(planes) < 63:
        p = rng.uniform(0.1, math.pi - 0.1)
        if abs(p - 3 * math.pi / 4) < 0.1 or abs(p - math.pi / 2) < 1e-3:
            continue
        planes.append(Plane(p, label="random angle, the Axiom's length", kind="axiom"))
    while len(planes) < 83:
        p = rng.uniform(0.1, math.pi - 0.1)
        if abs(p - math.pi / 2) < 1e-3:
            continue
        planes.append(Plane(p, rng.uniform(0.3, 3.0), label="random angle, free length", kind="free"))
    planes.append(Plane(math.radians(89.9), label="89.9°, next to the complex plane", kind="near"))
    planes.append(Plane(math.radians(90.1), label="90.1°, next to the complex plane", kind="near"))
    return planes


def evaluate(fn, *args):
    try:
        e = fn(*args)
    except (ZeroDivisionError, OverflowError, ValueError):
        return float("inf")
    return float("inf") if e != e else e          # a result that is not a number counts as a failure


def run_laws(planes):
    results = []
    for L_ in LAWS:
        if L_["planes"] == "all":
            use = planes
        elif L_["planes"] == "complex":
            use = [P for P in planes if P.kind == "complex"]
        elif L_["planes"] == "leaning":
            use = [P for P in planes if abs(P.c) > 0.05]
        else:
            use = [None]
        lrng = random.Random(f"{SEED}-{L_['name']}-{L_['ref']}")
        if L_["negative"]:
            # a wrong law: its largest miss in each plane, and the smallest of these must still be large
            per_plane, cases = [], 0
            for P in use:
                m = 0.0
                for _ in range(L_["samples"]):
                    m = max(m, evaluate(L_["fn"], P, lrng))
                    cases += 1
                per_plane.append((m, P))
            smallest, Pmin = min(per_plane, key=lambda t: t[0])
            results.append(dict(L_, cases=cases, worst=smallest, worst_at=f"{Pmin.deg():.2f}°",
                                passed=smallest > NEGATIVE_MIN, nplanes=len(use)))
            continue
        worst, worst_at, cases = 0.0, "", 0
        for P in use:
            for _ in range(L_["samples"]):
                e = evaluate(L_["fn"], P, lrng) if P is not None else evaluate(L_["fn"], lrng)
                cases += 1
                if e > worst:
                    worst = e
                    worst_at = f"{P.deg():.2f}°" if P is not None else ""
        results.append(dict(L_, cases=cases, worst=worst, worst_at=worst_at, passed=worst <= L_["tol"],
                            nplanes=len(use) if use != [None] else 0))
    return results


def sensitivity(planes):
    """Plant small errors in the laws, one at a time, and count how many laws then fail."""
    g = globals()
    plants = []

    def conj_wrong(self):
        return GN(self.a + self.b * self.P.L * self.P.c, -self.b, self.P)

    def mul_wrong(self, o):
        if not isinstance(o, GN):
            return GN(self.a * o, self.b * o, self.P)
        P = self.P
        return GN(self.a * o.a - P.L ** 2 * self.b * o.b,
                  self.a * o.b + self.b * o.a + 2 * P.L * P.c * (1 + 1e-6) * self.b * o.b, P)

    def shadow_wrong(self):
        return self.a + self.b * self.P.L * self.P.c * (1 + 1e-6)

    cases = [
        ("The conjugate without its factor 2: (a + b|" + GU + "| cos " + PB + ") − b" + GU, ("GN", "conj", conj_wrong)),
        ("The Rule of the Square with 2|" + GU + "| cos " + PB + " off by one part in a million", ("GN", "__mul__", mul_wrong)),
        ("The Turning Unit with its real part off by 0.1%", ("g", "J_of", lambda P: GN(-P.c / P.s * 1.001, 1.0 / (P.L * P.s), P))),
        ("The general cosine off by one part in a million", ("g", "gcos", lambda al, ph: math.sin(ph) / math.sin(al + ph) * (1 + 1e-6))),
        ("The shadow with its term b|" + GU + "| cos " + PB + " off by one part in a million", ("GN", "shadow", shadow_wrong)),
    ]
    for label, (where, name, fake) in cases:
        if where == "GN":
            old = getattr(GN, name)
            setattr(GN, name, fake)
            if name == "__mul__":
                GN.__rmul__ = fake
        else:
            old = g[name]
            g[name] = fake
        try:
            res = run_laws(planes)
        finally:
            if where == "GN":
                setattr(GN, name, old)
                if name == "__mul__":
                    GN.__rmul__ = old
            else:
                g[name] = old
        failing = [r for r in res if not r["negative"] and not r["passed"]]
        new_failing = [r for r in failing if r["group"] == "parts" or (r["group"] == "complex" and r["name"].startswith("4.9.7"))]
        plants.append(dict(label=label, failing=len(failing), total=len([r for r in res if not r["negative"]]),
                           new=len(new_failing), examples=[r["name"] for r in failing[:3]], found=len(failing) > 0))
    return plants


def complex_summary(results, xs, appr):
    """The answer to the question: are the complex numbers still inside, as the special case pi/2?"""
    cl = [r for r in results if r["group"] == "complex"]
    cl_new = [r for r in cl if r["name"].startswith("4.9.7") or r["name"].startswith("App. F.8") or r["name"].startswith("The habits")]
    habits = [r for r in results if r["negative"] and r["fn"].__name__.startswith("n_p_")]
    xc = [x for x in xs if x["plane"] == "π/2"]
    near = [r for r in results if not r["negative"] and r["group"] == "parts"]
    ok = (all(r["passed"] for r in cl) and all(x["ok"] for x in xc) and all(r["passed"] for r in habits) and appr["ok"]
          and all(r["passed"] for r in near))
    return dict(ok=ok, n_cl=len(cl), ok_cl=sum(r["passed"] for r in cl), n_new=len(cl_new), ok_new=sum(r["passed"] for r in cl_new),
                worst_cl=max(r["worst"] for r in cl), n_xc=len(xc), ok_xc=sum(x["ok"] for x in xc),
                cases_xc=sum(x["cases"] for x in xc), n_hab=len(habits), ok_hab=sum(r["passed"] for r in habits))


def run():
    t0 = time.time()
    rng = random.Random(SEED)
    planes = build_planes(rng)
    results = run_laws(planes)
    wex = worked_examples()
    xs = exact_checks()
    appr = approach()
    plants = sensitivity(planes)
    cs = complex_summary(results, xs, appr)
    return planes, results, wex, xs, appr, plants, cs, time.time() - t0


# =====================================================================================
# 7. Reports
# =====================================================================================
def fmt_err(e):
    if e == 0:
        return "0"
    return f"{e:.1e}"


def text_report(planes, results, wex, xs, appr, plants, cs, secs):
    lines = []
    W = lines.append
    exact = [r for r in results if not r["negative"]]
    neg = [r for r in results if r["negative"]]
    ok_laws = sum(r["passed"] for r in exact)
    ok_neg = sum(r["passed"] for r in neg)
    W("DIGITAL TEST OF THE ALGEBRA OF GENERAL NUMBERS (version 2)")
    W("The Cosmic Mass, Volume 1: Section 4.9 with Subsection 4.9.7, and Appendix F with F.8 (draft of 6 November 2026)")
    W(f"Run on {datetime.datetime.now():%Y-%m-%d %H:%M}, Python {platform.python_version()}, seed {SEED}, {secs:.1f} s")
    W("")
    W(f"Laws:              {ok_laws} of {len(exact)} passed ({sum(r['group'] == 'parts' for r in exact)} of them new in 4.9.7 and F.8)")
    W(f"Evaluations:       {sum(r['cases'] for r in exact):,}")
    W(f"Largest error:     {max(r['worst'] for r in exact if r['tol'] == TOL):.1e} for the exact laws, "
      f"{max(r['worst'] for r in exact if r['tol'] != TOL):.1e} with finite differences")
    W(f"Negative controls: {ok_neg} of {len(neg)} wrong laws caught")
    W(f"Complex case:      {cs['ok_cl']} of {cs['n_cl']} laws at pi/2 match Appendix C; exact: {cs['ok_xc']} of {cs['n_xc']} checks hold")
    W(f"Exact checks:      {sum(x['ok'] for x in xs)} of {len(xs)} hold, with no rounding ({sum(x['cases'] for x in xs):,} exact cases)")
    W(f"Worked examples:   {sum(r['ok'] for r in wex)} of {len(wex)} printed numbers match")
    W(f"Planted errors:    {sum(p['found'] for p in plants)} of {len(plants)} found")
    W(f"Planes:            {len(planes)} (79.86°, π/3, π/2, 89.9°, 90.1°, 60 random angles with the Axiom's length, 20 with a free length)")
    W("")
    W("IS THE COMPLEX PLANE STILL INSIDE, AS THE SPECIAL CASE pi/2?  " + ("YES" if cs["ok"] else "NOT SHOWN: see below"))
    W(f"  - All {cs['n_cl']} laws of the complex case, {cs['n_new']} of them for 4.9.7 and F.8, hold at pi/2 against Appendix C, "
      f"computed with Python's complex numbers; largest error {cs['worst_cl']:.1e}.")
    W(f"  - Exact checks at pi/2, no rounding: {cs['ok_xc']} of {cs['n_xc']} hold on {cs['cases_xc']:,} exact numbers. There the arithmetic "
      f"of the leaning frame is exactly the arithmetic of the complex numbers.")
    W(f"  - The {cs['n_hab']} habits of the complex plane, used for the parts, are caught as wrong in every leaning plane "
      f"({cs['ok_hab']} of {cs['n_hab']}), and they hold at pi/2.")
    r0, r1 = appr["rows"][0], appr["rows"][-1]
    W(f"  - Approaching 90°, their largest miss falls from {r0['miss']:.1e} at {r0['deg']}° to {r1['miss']:.1e} at 90°, "
      f"while the true laws hold at every step.")
    W("")
    for key, title, note in GROUPS:
        rs = [r for r in results if r["group"] == key]
        if not rs:
            continue
        W(title)
        W("-" * len(title))
        for r in rs:
            tag = ("CAUGHT" if r["passed"] else "MISSED") if r["negative"] else ("PASS" if r["passed"] else "FAIL")
            W(f"  {tag:6}  {r['ref']:22} {r['name']}")
            W(f"          {r['formula']}")
            W(f"          {r['method']}; {r['cases']:,} cases; {'smallest miss' if r['negative'] else 'largest error'} {fmt_err(r['worst'])}")
        W("")
    W("Exact checks, with no rounding (numbers p + q√3 with fractions p and q)")
    W("-----------------------------------------------------------------------")
    for x in xs:
        W(f"  {'PASS' if x['ok'] else 'FAIL':6}  {x['plane']:5} {x['ref']:18} {x['name']} ({x['cases']} cases)")
    W("")
    W("Approaching the complex plane: the habits of the complex plane, and the true laws of 4.9.7")
    W("-------------------------------------------------------------------------------------------")
    W(f"  {'angle':>9}  {'|cos φ̄|':>9}  {'miss of the habits':>19}  {'error of the true laws':>23}")
    for r in appr["rows"]:
        W(f"  {r['deg']:>8}°  {r['cos']:9.2e}  {r['miss']:19.2e}  {r['err']:23.2e}")
    W("  " + ("PASS" if appr["ok"] else "FAIL") + "  the misses fall toward 90° and vanish there; the true laws hold at every angle")
    W("")
    W("Worked examples: printed numbers against exact computation")
    W("-----------------------------------------------------------")
    for r in wex:
        W(f"  {'PASS' if r['ok'] else 'FAIL':6}  {r['where']:14} {r['what']}: printed {r['printed']:.{r['decimals']}f}, "
          f"computed {r['computed']:.{r['decimals'] + 3}f}")
    W("")
    W("Planted errors: the test must find each of them")
    W("-----------------------------------------------")
    for pl in plants:
        W(f"  {'FOUND' if pl['found'] else 'MISSED':6}  {pl['label']}: {pl['failing']} of {pl['total']} laws fail "
          f"({pl['new']} of them laws of 4.9.7 and F.8), for example {', '.join(pl['examples'])}")
    return "\n".join(lines) + "\n"


CSS = """
:root { --ink:#1d2433; --muted:#5b6475; --line:#d9dee7; --soft:#f4f6fa; --ok:#1f7a4d; --oksoft:#e7f4ec;
        --bad:#b42318; --badsoft:#fdecea; --accent:#2f4f8f; }
* { box-sizing: border-box; }
body { margin: 0; background: #ffffff; color: var(--ink);
       font-family: "Palatino Linotype", Palatino, "Book Antiqua", "TeX Gyre Pagella", "Latin Modern Roman", Georgia,
       "Cambria Math", "Latin Modern Math", "DejaVu Math TeX Gyre", serif; font-size: 10.5pt; line-height: 1.45; }
main { max-width: 980px; margin: 0 auto; padding: 28px 32px 48px; }
h1 { font-size: 22pt; margin: 0 0 4px; font-weight: normal; letter-spacing: .01em; }
h2 { font-size: 13.5pt; margin: 28px 0 6px; font-weight: normal; color: var(--accent);
     border-bottom: 1px solid var(--line); padding-bottom: 4px; break-after: avoid; page-break-after: avoid; }
h2 + p.note, h2 + table, p.note + table { break-before: avoid; page-break-before: avoid; }
thead { display: table-header-group; }
table.laws { table-layout: fixed; }
sub, sup { font-size: 75%; line-height: 0; }
p { margin: 6px 0; }
.sub { color: var(--muted); margin: 0 0 18px; }
.note { color: var(--muted); font-size: 9.5pt; margin: 2px 0 8px; }
.summary { display: grid; grid-template-columns: repeat(4, 1fr); gap: 10px; margin: 14px 0 8px; }
.tile { background: var(--soft); border: 1px solid var(--line); border-radius: 6px; padding: 10px 12px; }
.tile b { display: block; font-size: 17pt; font-weight: normal; }
.tile span { color: var(--muted); font-size: 9.5pt; }
.verdict { margin: 12px 0 4px; padding: 10px 14px; border-radius: 6px; font-size: 12pt; }
.verdict.ok { background: var(--oksoft); color: var(--ok); border: 1px solid #b9dfc8; }
.verdict.bad { background: var(--badsoft); color: var(--bad); border: 1px solid #f3c1bb; }
table { width: 100%; border-collapse: collapse; margin: 6px 0 4px; font-size: 9.5pt; }
th { text-align: left; font-weight: normal; color: var(--muted); font-size: 8.5pt; text-transform: uppercase;
     letter-spacing: .04em; border-bottom: 1px solid var(--line); padding: 4px 6px; }
td { border-bottom: 1px solid #edf0f5; padding: 5px 6px; vertical-align: top; }
tr { break-inside: avoid; }
td.num { text-align: right; white-space: nowrap; font-variant-numeric: tabular-nums; }
td.ref { color: var(--muted); }
td.f { color: #2b3446; }
.pass { color: var(--ok); white-space: nowrap; }
.fail { color: var(--bad); white-space: nowrap; font-weight: bold; }
ul { margin: 4px 0 8px 18px; padding: 0; }
li { margin: 2px 0; }
.ov { text-decoration: overline; text-decoration-thickness: 0.07em; }
footer { margin-top: 28px; color: var(--muted); font-size: 9pt; border-top: 1px solid var(--line); padding-top: 8px; }
@page { size: A4; margin: 14mm 12mm 16mm; }
@media print { main { padding: 0; max-width: none; } h2 { margin-top: 20px; } }
"""


def math_html(text):
    """Escape a formula for HTML and set x_a, x_{ab}, x_(ab) as subscripts and e^(ab) as superscripts."""
    import re
    t = html.escape(text)
    t = t.replace(GU + "\u0304", "<span class='ov'>" + GU + "</span>").replace(JU + "\u0304", "<span class='ov'>" + JU + "</span>")
    t = re.sub(r"_\{([^}]*)\}", r"<sub>\1</sub>", t)
    t = re.sub(r"_\(([^)]*)\)", r"<sub>\1</sub>", t)
    t = re.sub(r"\^\(([^)]*)\)", r"<sup>\1</sup>", t)
    t = re.sub(r"_([A-Za-z0-9αβθ]+)", r"<sub>\1</sub>", t)
    return t


def html_report(planes, results, wex, xs, appr, plants, cs, secs):
    E = html.escape
    M = math_html
    exact = [r for r in results if not r["negative"]]
    neg = [r for r in results if r["negative"]]
    ok_laws = sum(r["passed"] for r in exact)
    ok_neg = sum(r["passed"] for r in neg)
    n_new = sum(r["group"] == "parts" for r in exact)
    all_ok = (ok_laws == len(exact) and ok_neg == len(neg) and all(r["ok"] for r in wex) and all(x["ok"] for x in xs)
              and appr["ok"] and all(pl["found"] for pl in plants))
    n_eval = sum(r["cases"] for r in exact)
    worst = max(r["worst"] for r in exact if r["tol"] == TOL)
    worst_fd = max(r["worst"] for r in exact if r["tol"] != TOL)
    out = []
    W = out.append
    W("<!doctype html><html lang='en'><head><meta charset='utf-8'>")
    W("<meta name='viewport' content='width=device-width, initial-scale=1'>")
    W("<title>General Number Laws: Test Report</title><style>" + CSS + "</style></head><body><main>")
    W("<h1>The Algebra of General Numbers: a digital test</h1>")
    W(f"<p class='sub'>The Cosmic Mass, Volume 1 &middot; Section 4.9 with Subsection 4.9.7, and Appendix F with F.8 "
      f"(draft of 6 November 2026) &middot; version 2 of the test &middot; run on {datetime.datetime.now():%d %B %Y, %H:%M}</p>")
    if all_ok:
        W(f"<div class='verdict ok'>Every law held: {ok_laws} of {len(exact)} laws passed on {n_eval:,} random cases "
          f"({n_new} of them the new laws of 4.9.7 and F.8), all {len(neg)} wrong laws were caught, every planted error was found, "
          f"all {len(xs)} exact checks hold, and every printed number checks out.</div>")
    else:
        W(f"<div class='verdict bad'>Not everything held: {len(exact) - ok_laws} laws failed, {len(neg) - ok_neg} controls were missed, "
          f"{sum(not x['ok'] for x in xs)} exact checks failed. See the lines marked FAIL or MISSED.</div>")
    W("<div class='summary'>")
    tiles = [(f"{ok_laws} / {len(exact)}", f"laws passed, {n_new} of them new in 4.9.7 and F.8"),
             (f"{n_eval:,}", "random cases computed"),
             (f"{worst:.1e}", f"largest error of an exact law ({worst_fd:.1e} with finite differences)"),
             (f"{ok_neg} / {len(neg)}", "wrong laws caught"),
             (f"{cs['ok_cl']} / {cs['n_cl']}", "laws at π/2 that match Appendix C"),
             (f"{sum(x['ok'] for x in xs)} / {len(xs)}", "exact checks, with no rounding"),
             (f"{sum(r['ok'] for r in wex)} / {len(wex)}", "printed numbers of the worked examples"),
             (f"{sum(p['found'] for p in plants)} / {len(plants)}", "planted errors found")]
    for big, small in tiles:
        W(f"<div class='tile'><b>{E(big)}</b><span>{E(small)}</span></div>")
    W("</div>")

    r0, r1 = appr["rows"][0], appr["rows"][-1]
    W("<h2>Is the complex plane still inside, as the special case " + PB + " = π/2?</h2>")
    W(f"<div class='verdict {'ok' if cs['ok'] else 'bad'}'>{'Yes.' if cs['ok'] else 'Not shown: see the lines marked FAIL.'} "
      f"At {PB} = π/2 every law of the section, the new laws of the parts included, becomes its law of Appendix C, "
      f"exactly and not only in the limit.</div><ul>")
    W(f"<li><b>Against Python's complex numbers.</b> All {cs['n_cl']} laws of the complex case hold at π/2, "
      f"{cs['n_new']} of them for 4.9.7 and F.8; the largest error is {cs['worst_cl']:.1e}, the size of rounding. "
      f"There the real part and the shadow are both Re z, and the Composer part is i Im z.</li>")
    W(f"<li><b>Exactly, with no rounding.</b> With numbers p + q√3 whose p and q are fractions, {cs['ok_xc']} of {cs['n_xc']} "
      f"checks hold at π/2 on {cs['cases_xc']:,} exact cases: there the arithmetic of the leaning frame is the arithmetic of the "
      f"complex numbers, digit for digit, and the laws of 4.9.7 are the laws of the real and imaginary parts.</li>")
    W(f"<li><b>The habits of the complex plane.</b> {cs['n_hab']} rules that are true for complex numbers, such as "
      f"Re(z̄) = Re(z) or Com(zw) = Re(z)Com(w) + Re(w)Com(z), are used for the parts as negative controls. "
      f"All {cs['ok_hab']} fail in every leaning plane and hold at π/2. So the new laws are needed exactly where the axis leans, "
      f"and they reduce to the old ones exactly where it stands upright.</li>")
    W(f"<li><b>Approaching 90°.</b> The largest miss of these habits falls from {r0['miss']:.1e} at {r0['deg']}° to "
      f"{r1['miss']:.1e} at 90°, in step with |cos {PB}|, while the true laws hold at every angle on the way "
      f"(table below). The complex plane is the upright point of one family, not a separate rule.</li></ul>")

    W("<h2>How the test works</h2>")
    W("<p>Each law is computed in two independent ways, and the two results must agree.</p><ul>")
    W(f"<li><b>In the leaning frame</b>, as the book computes: z = a + b{GU} is the pair (a, b), and pairs are multiplied "
      f"with the Rule of the Square, {GU}² = 2|{GU}| cos {PB} {GU} − |{GU}|² (Propositions 4.22 and 4.23).</li>")
    W(f"<li><b>On the page</b>, as arrows: z is the point a + b|{GU}|(cos {PB} + i sin {PB}), and arrows are multiplied by "
      f"Wessel's rule (lengths multiply, angles add), which is the multiplication of complex numbers that Python does by itself. "
      f"The real part and the Composer coordinate of an arrow are read back from the page by solving two linear equations; "
      f"its shadow and its height are its two coordinates.</li>")
    W(f"<li>The General Sinus always comes from its definitions in Appendix A, sin(α, φ) = sin α / sin(α + φ) and "
      f"cos(α, φ) = sin φ / sin(α + φ), never from the law under test. Here C<sub>α</sub> = cos(α, {PB} − α), "
      f"S<sub>α</sub> = sin(α, {PB} − α) and u = {GU}/|{GU}|, the direction of the Composer Unit.</li>")
    W("<li>Further independent checks: the exponential is summed from its power series in the leaning frame; derivatives "
      "come from finite differences; all roots of a polynomial come from an independent root finder (Durand-Kerner); "
      "and the exact checks use numbers p + q√3, where p and q are fractions, in the planes π/2 and π/3, whose cosine, sine "
      f"and |{GU}| are such numbers.</li>")
    W(f"<li>A law passes when its largest relative error, over all cases, is below {TOL:g} ({TOL_DIFF:g} for the laws checked with "
      f"finite differences). Ordinary computer arithmetic carries about 16 digits, so a true law shows errors near 1e-15, "
      f"and a wrong one misses by far more.</li></ul>")

    W("<h2>The planes tested</h2><table><thead><tr><th>planes</th><th>angle of the Composer axis</th><th>length of the Composer Unit</th></tr></thead>")
    W(f"<tr><td>1</td><td>79.86°, the book's angle</td><td>{planes[0].L:.5f}, from the Axiom</td></tr>")
    W(f"<tr><td>1</td><td>π/3, where cos {PB}, sin {PB} and |{GU}| are exact numbers</td><td>{planes[1].L:.5f}, from the Axiom</td></tr>")
    W("<tr><td>1</td><td>π/2, the complex plane</td><td>1, from the Axiom</td></tr>")
    W("<tr><td>2</td><td>89.9° and 90.1°, next to the complex plane</td><td>from the Axiom</td></tr>")
    W("<tr><td>60</td><td>random, between 0.1 and π − 0.1, away from 3π/4</td><td>from the Axiom, |cos " + PB + " + sin " + PB + "|</td></tr>")
    W("<tr><td>20</td><td>random, between 0.1 and π − 0.1</td><td>free, between 0.3 and 3</td></tr></table>")
    rand = ", ".join(f"{P.deg():.2f}°" + (f" ({P.L:.3f})" if P.kind == "free" else "") for P in planes[3:83])
    W(f"<p class='note'>The random planes, by the angle of the Composer axis (with the free length in brackets): {E(rand)}.</p>")
    W(f"<p class='note'>Angles too close to 0, π or 3π/4 are left out: there the relative plane closes into a line (Remark 4.13). "
      f"Each law takes {SAMPLES} random numbers in each plane, unless the table says otherwise; the random numbers are the same at every run (seed {SEED}).</p>")

    for key, title, note in GROUPS:
        rs = [r for r in results if r["group"] == key]
        if not rs:
            continue
        W(f"<h2>{E(title)}</h2>")
        if note:
            W(f"<p class='note'>{M(note)}</p>")
        neg_group = key == "neg"
        W("<table class='laws'><colgroup><col style='width:12%'><col style='width:19%'><col style='width:34%'>"
          "<col style='width:12%'><col style='width:7%'><col style='width:8%'><col style='width:8%'></colgroup>"
          "<thead><tr><th>in the book</th><th>law</th><th>formula</th><th>method</th><th class='num'>cases</th>"
          f"<th class='num'>{'smallest miss' if neg_group else 'largest error'}</th><th>result</th></tr></thead>")
        for r in rs:
            if r["negative"]:
                res = "<span class='pass'>&#10003; caught</span>" if r["passed"] else "<span class='fail'>&#10007; missed</span>"
            else:
                res = "<span class='pass'>&#10003; pass</span>" if r["passed"] else "<span class='fail'>&#10007; FAIL</span>"
            W(f"<tr><td class='ref'>{E(r['ref'])}</td><td>{M(r['name'])}</td><td class='f'>{M(r['formula'])}</td>"
              f"<td>{E(r['method'])}</td><td class='num'>{r['cases']:,}</td><td class='num'>{fmt_err(r['worst'])}</td><td>{res}</td></tr>")
        W("</table>")

    W("<h2>Exact checks, with no rounding</h2>")
    W(f"<p class='note'>In the planes π/2 and π/3 the cosine, the sine and the length |{GU}| are exact numbers p + q√3, with "
      f"fractions p and q (at π/3: 1/2, √3/2 and (1 + √3)/2). Each check takes {EXACT_SAMPLES} random exact numbers, computes the "
      f"law in the leaning frame and on the page with exact fractions, and asks for equality to the last digit. "
      f"At π/2 the arithmetic of the leaning frame must be exactly that of the complex numbers.</p>")
    W("<table><thead><tr><th>in the book</th><th>check</th><th>plane</th><th class='num'>cases</th><th>result</th></tr></thead>")
    for x in xs:
        res = "<span class='pass'>&#10003; exact</span>" if x["ok"] else "<span class='fail'>&#10007; FAIL</span>"
        W(f"<tr><td class='ref'>{E(x['ref'])}</td><td>{M(x['name'])}</td><td>{E(x['plane'])}</td>"
          f"<td class='num'>{x['cases']}</td><td>{res}</td></tr>")
    W("</table>")

    W("<h2>Approaching the complex plane</h2>")
    W(f"<p class='note'>The {cs['n_hab']} habits of the complex plane of the negative controls, and eleven true laws of 4.9.7, "
      f"in planes whose Composer axis rises toward the upright. The miss of the habits shrinks with |cos {PB}| and vanishes at 90°; "
      f"the true laws hold at every angle.</p>")
    W("<table><thead><tr><th>angle of the Composer axis</th><th class='num'><span style='text-transform:none'>|cos " + PB + "|</span></th>"
      "<th class='num'>largest miss of the habits</th><th class='num'>largest error of the true laws</th></tr></thead>")
    for r in appr["rows"]:
        W(f"<tr><td>{r['deg']}°</td><td class='num'>{r['cos']:.2e}</td><td class='num'>{r['miss']:.2e}</td>"
          f"<td class='num'>{r['err']:.2e}</td></tr>")
    W("</table>")
    W("<p class='note'>" + ("&#10003; The misses fall toward 90° and vanish there, and the true laws hold at every angle."
                            if appr["ok"] else "&#10007; The approach does not behave as expected.") + "</p>")

    W("<h2>Planted errors: the test must find each of them</h2>")
    W("<p class='note'>To be sure that the test can see a mistake, the script plants one small error at a time in its own "
      "copy of the laws, runs every law again, and counts the laws that then fail. Each planted error must make at least one law fail.</p>")
    W("<table><thead><tr><th>planted error</th><th class='num'>laws that fail</th><th class='num'>of them in 4.9.7, F.8</th>"
      "<th>for example</th><th>result</th></tr></thead>")
    for pl in plants:
        res = "<span class='pass'>&#10003; found</span>" if pl["found"] else "<span class='fail'>&#10007; missed</span>"
        W(f"<tr><td>{M(pl['label'])}</td><td class='num'>{pl['failing']} of {pl['total']}</td><td class='num'>{pl['new']}</td>"
          f"<td>{E('; '.join(pl['examples']))}</td><td>{res}</td></tr>")
    W("</table>")

    W("<h2>The worked examples of Section 4.9</h2>")
    W("<p class='note'>Every number printed in the worked examples, and in Figure 4.20, recomputed from scratch at 79.86° and "
      "compared with the printed value to the number of decimals printed.</p>")
    W("<table><thead><tr><th>where</th><th>quantity</th><th class='num'>printed</th><th class='num'>computed</th><th>result</th></tr></thead>")
    for r in wex:
        res = "<span class='pass'>&#10003; match</span>" if r["ok"] else "<span class='fail'>&#10007; differs</span>"
        W(f"<tr><td class='ref'>{E(r['where'])}</td><td>{M(r['what'])}</td><td class='num'>{r['printed']:.{r['decimals']}f}</td>"
          f"<td class='num'>{r['computed']:.{r['decimals'] + 3}f}</td><td>{res}</td></tr>")
    W("</table>")

    W(f"<footer>Generated by test_general_number_laws.py (version 2) with Python {E(platform.python_version())} on "
      f"{E(platform.system())}; {secs:.1f} seconds. Only the Python standard library is used. "
      f"Run <code>python3 test_general_number_laws.py</code> to reproduce this report.</footer>")
    W("</main></body></html>")
    return "\n".join(out)


def main():
    planes, results, wex, xs, appr, plants, cs, secs = run()
    txt = text_report(planes, results, wex, xs, appr, plants, cs, secs)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "test-report.txt"), "w", encoding="utf-8") as f:
        f.write(txt)
    with open(os.path.join(here, "test-report.html"), "w", encoding="utf-8") as f:
        f.write(html_report(planes, results, wex, xs, appr, plants, cs, secs))
    exact = [r for r in results if not r["negative"]]
    neg = [r for r in results if r["negative"]]
    print("\n\n".join(txt.split("\n\n")[:3]))
    bad = [r for r in exact if not r["passed"]] + [r for r in neg if not r["passed"]]
    bad_w = [r for r in wex if not r["ok"]]
    bad_x = [x for x in xs if not x["ok"]]
    bad_p = [p for p in plants if not p["found"]]
    if bad or bad_w or bad_x or bad_p or not appr["ok"]:
        print("\nPROBLEMS:")
        for r in bad:
            print(f"  {r['ref']} {r['name']}: {fmt_err(r['worst'])} at {r['worst_at']}")
        for r in bad_w:
            print(f"  worked example {r['where']} {r['what']}: printed {r['printed']}, computed {r['computed']}")
        for x in bad_x:
            print(f"  exact check {x['ref']} {x['name']} at {x['plane']}")
        for p in bad_p:
            print(f"  planted error not found: {p['label']}")
        if not appr["ok"]:
            print("  the approach to the complex plane does not behave as expected")
    else:
        print("\nAll laws held. Reports written: test-report.html and test-report.txt")
    return 0 if not (bad or bad_w or bad_x or bad_p or not appr["ok"]) else 1


if __name__ == "__main__":
    sys.exit(main())
