#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_general_number_algebra.py

Checks for "The Algebra of General Numbers": Section 4.9 and Appendix F of
The Cosmic Mass, Volume 1 (draft of 6 November 2026), including Subsection 4.9.7
and Appendix F.8, the laws of the parts beside the shadow and the height.

The script checks six things:
  1. every law, proved with symbolic algebra for any angle phibar of the Composer
     axis and any length L = |G| of the Composer Unit;
  2. every law on arrows multiplied directly on the page (as complex numbers),
     at 79.86 degrees, pi/3 and other angles;
  3. that phibar = pi/2, where the Axiom gives L = 1, turns each law into its
     formula in Appendix C;
  4. the formulas of Appendix A that the section uses;
  5. every number printed in the worked examples of Section 4.9;
  6. every law of Subsection 4.9.7 for the real part, the Composer part, the shadow
     and the height, symbolically and on arrows of the page (sections 15 and 16).

Run with:  python3 verify_general_number_algebra.py
Needs:     Python 3 and SymPy (pip install sympy)

Notation: a General Number of the relative plane with angle p and length L is
z = a + b G. On the page, G is the arrow of length L at the angle p, so z is the
complex number a + b L (cos p + i sin p). The product uses the Rule of the
Square, G^2 = 2 L cos(p) G - L^2 (Proposition 4.22).
"""
import cmath
import math
import random

import sympy as sp

PASSED, FAILED = [], []


def report(name, ok, how=""):
    (PASSED if ok else FAILED).append(name)
    tag = "  ok    " if ok else "  FAIL  "
    print(tag + name + (f"   [{how}]" if how and ok else ""))


def heading(title):
    print("\n" + title)
    print("-" * len(title))


# --------------------------------------------------------------------------
# Symbols and the General Number arithmetic
# --------------------------------------------------------------------------
a, b, c, d, x, y = sp.symbols("a b c d x y", real=True)
al, be, th = sp.symbols("alpha beta theta", real=True)
p, q = sp.symbols("phibar psibar", real=True)  # angles of two Composer axes
L, K = sp.symbols("L K", positive=True)        # lengths of the two Composer Units
r, s_ = sp.symbols("r s", positive=True)
I = sp.I


def zero(expr):
    """True if expr is identically zero. Symbolic first; if SymPy cannot decide,
    a numeric test at 40 random points is reported instead."""
    e = sp.expand(expr)
    for f in (sp.simplify, lambda t: sp.simplify(sp.expand_trig(t)),
              lambda t: sp.simplify(t.rewrite(sp.exp)),
              lambda t: sp.trigsimp(sp.expand_trig(t))):
        try:
            if f(e) == 0:
                return True, "symbolic"
        except Exception:
            pass
    syms = sorted(e.free_symbols, key=lambda t: t.name)
    fn = sp.lambdify(syms, e, "mpmath")
    rnd = random.Random(11)
    for _ in range(40):
        vals = []
        for t in syms:
            if t in (L, K, r, s_):
                vals.append(rnd.uniform(0.3, 2.0))
            elif t in (p, q):
                vals.append(rnd.uniform(0.2, 2.9))
            else:
                vals.append(rnd.uniform(-2.5, 2.5))
        try:
            v = complex(fn(*vals))
        except ZeroDivisionError:
            continue
        if abs(v) > 1e-9:
            return False, "numeric"
    return True, "numeric only"


class GN:
    """A General Number a + b G of the relative plane with angle p and length L."""

    def __init__(self, a_, b_, p_=p, L_=L):
        self.a, self.b, self.p, self.L = sp.sympify(a_), sp.sympify(b_), p_, L_

    def lift(self, o):
        return o if isinstance(o, GN) else GN(o, 0, self.p, self.L)

    def __add__(self, o):
        o = self.lift(o)
        return GN(self.a + o.a, self.b + o.b, self.p, self.L)

    __radd__ = __add__

    def __neg__(self):
        return GN(-self.a, -self.b, self.p, self.L)

    def __sub__(self, o):
        return self + (-self.lift(o))

    def __rsub__(self, o):
        return self.lift(o) - self

    def __mul__(self, o):
        # Proposition 4.23: (a1 + b1 G)(a2 + b2 G)
        #   = (a1 a2 - L^2 b1 b2) + (a1 b2 + b1 a2 + 2 L cos p b1 b2) G
        o = self.lift(o)
        k = 2 * self.L * sp.cos(self.p)
        return GN(self.a * o.a - self.L**2 * self.b * o.b,
                  self.a * o.b + self.b * o.a + k * self.b * o.b, self.p, self.L)

    __rmul__ = __mul__

    def conj(self):
        # Definition: the conjugate (a + 2 L b cos p) - b G
        return GN(self.a + 2 * self.L * self.b * sp.cos(self.p), -self.b, self.p, self.L)

    def norm2(self):
        # |z|^2 = a^2 + 2 a b L cos p + b^2 L^2 (Theorem 4.9)
        return self.a**2 + 2 * self.L * self.a * self.b * sp.cos(self.p) + self.L**2 * self.b**2

    def inv(self):
        n2 = self.norm2()
        zc = self.conj()
        return GN(zc.a / n2, zc.b / n2, self.p, self.L)

    def __truediv__(self, o):
        return self * self.lift(o).inv()

    def __pow__(self, n):
        out = GN(1, 0, self.p, self.L)
        base = self if n >= 0 else self.inv()
        for _ in range(abs(n)):
            out = out * base
        return out

    def shadow(self):
        return self.a + self.b * self.L * sp.cos(self.p)

    def height(self):
        return self.b * self.L * sp.sin(self.p)

    def page(self):
        return self.a + self.b * self.L * (sp.cos(self.p) + I * sp.sin(self.p))

    def subs(self, *args, **kw):
        return GN(self.a.subs(*args, **kw), self.b.subs(*args, **kw), self.p, self.L)


def same(z, w):
    """Two General Numbers are equal when both parts are (Definition 4.12)."""
    ok1, h1 = zero(z.a - w.a)
    ok2, h2 = zero(z.b - w.b)
    how = "symbolic" if (h1 == h2 == "symbolic") else "numeric only"
    return ok1 and ok2, how


def page_eq(e1, e2):
    ok1, h1 = zero(sp.re(sp.expand(e1 - e2)))
    ok2, h2 = zero(sp.im(sp.expand(e1 - e2)))
    how = "symbolic" if (h1 == h2 == "symbolic") else "numeric only"
    return ok1 and ok2, how


# The General Sinus (Appendix A, Definitions)
def gsin(alpha, phi):
    return sp.sin(alpha) / sp.sin(alpha + phi)


def gcos(alpha, phi):
    return sp.sin(phi) / sp.sin(alpha + phi)


# Along the Composer axis the two angles add up to phibar: phi = phibar - alpha
def Cg(t, p_=p):
    return gcos(t, p_ - t)       # = sin(phibar - t) / sin(phibar)


def Sg(t, p_=p):
    return gsin(t, p_ - t)       # = sin(t) / sin(phibar)


G = GN(0, 1)                     # the Composer Unit
u = GN(0, 1 / L)                 # its direction, G/|G|, a General Number of length 1
J = GN(-sp.cos(p) / sp.sin(p), 1 / (L * sp.sin(p)))   # the Turning Unit


def E(t, p_=p, L_=L):
    """The unit General Number at the angle t: cos(t, p - t) + sin(t, p - t) G/|G|."""
    return GN(Cg(t, p_), Sg(t, p_) / L_, p_, L_)


z1, z2, z3 = GN(a, b), GN(c, d), GN(x, y)

# --------------------------------------------------------------------------
heading("1. The arithmetic agrees with arrows on the page")
report("Rule of the Square: page(G*G) = page(G)^2", *page_eq((G * G).page(), G.page() ** 2))
report("Product of Proposition 4.23 is Wessel's product on the page",
       *page_eq((z1 * z2).page(), z1.page() * z2.page()))

heading("2. The rules of arithmetic (Theorem: a true relative plane is a field)")
report("order of a product: z1 z2 = z2 z1", *same(z1 * z2, z2 * z1))
report("grouping: (z1 z2) z3 = z1 (z2 z3)", *same((z1 * z2) * z3, z1 * (z2 * z3)))
report("distribution: z1 (z2 + z3) = z1 z2 + z1 z3", *same(z1 * (z2 + z3), z1 * z2 + z1 * z3))
report("neutral numbers: z + 0 = z and 1 z = z", same(z1 + 0, z1)[0] and same(GN(1, 0) * z1, z1)[0], "symbolic")
report("inverse: z (1/z) = 1", *same(z1 * z1.inv(), GN(1, 0)))

heading("3. The conjugate")
report("the conjugate is the mirror image in the real axis",
       *page_eq(z1.conj().page(), sp.conjugate(z1.page())))
report("conjugate of the conjugate is z", *same(z1.conj().conj(), z1))
report("parts of the conjugate (Definition 4.11): Re(conj z) = a + 2b|G|cos(phibar), Com(conj z) = -b G = -Com(z)",
       zero(z1.conj().a - (a + 2 * b * L * sp.cos(p)))[0] and zero(z1.conj().b + b)[0], "symbolic")
report("  Re(conj z) is Re(z) moved by twice the shadow of the Composer part b G",
       *zero(z1.conj().a - (a + 2 * GN(0, b).shadow())))
report("  the real part stays only where cos(phibar) = 0: at pi/2, Re(conj z) = a",
       *zero(z1.conj().a.subs(p, sp.pi / 2) - a))
report("conjugate of a sum", *same((z1 + z2).conj(), z1.conj() + z2.conj()))
report("conjugate of a product", *same((z1 * z2).conj(), z1.conj() * z2.conj()))
report("z + conj(z) = 2 Sh(z)", *same(z1 + z1.conj(), GN(2 * z1.shadow(), 0)))
report("z conj(z) = |z|^2, a real number", *same(z1 * z1.conj(), GN(z1.norm2(), 0)))
report("|z|^2 is the squared length on the page", *zero(z1.norm2() - sp.Abs(z1.page()) ** 2))
report("|z|^2 = Sh(z)^2 + Ht(z)^2", *zero(z1.norm2() - z1.shadow() ** 2 - z1.height() ** 2))
report("G + conj(G) = 2 L cos(phibar)", *same(G + G.conj(), GN(2 * L * sp.cos(p), 0)))
report("G conj(G) = L^2", *same(G * G.conj(), GN(L**2, 0)))
report("G and conj(G) are the roots of x^2 - 2L cos(phibar) x + L^2",
       same(G * G - 2 * L * sp.cos(p) * G + L**2, GN(0, 0))[0]
       and same(G.conj() * G.conj() - 2 * L * sp.cos(p) * G.conj() + L**2, GN(0, 0))[0], "symbolic")
skew = GN(a, -b)
prod = z1 * skew
report("the skew flip a - bG is not the mirror: (a+bG)(a-bG) has Composer coordinate -2L cos(phibar) b^2",
       *zero(prod.b + 2 * L * sp.cos(p) * b**2))
report("the skew flip changes the length: |a - bG|^2 - |a + bG|^2 = -4 a b L cos(phibar)",
       *zero(skew.norm2() - z1.norm2() + 4 * a * b * L * sp.cos(p)))

heading("4. Magnitude and division")
report("|z1 z2|^2 = |z1|^2 |z2|^2 (the general two-squares identity)",
       *zero((z1 * z2).norm2() - z1.norm2() * z2.norm2()))
N2 = c**2 + 2 * c * d * L * sp.cos(p) + d**2 * L**2
quot = GN((a * c + 2 * a * d * L * sp.cos(p) + L**2 * b * d) / N2, (b * c - a * d) / N2)
report("division formula: [(ac + 2ad L cos + L^2 bd) + (bc - ad) G] / (c^2 + 2cd L cos + d^2 L^2)",
       *same(quot * z2, z1))
report("the division formula is z1 conj(z2) / |z2|^2", *same(quot, GN((z1 * z2.conj()).a / N2, (z1 * z2.conj()).b / N2)))
report("1/G = (2 L cos(phibar) - G) / L^2", *same(G.inv(), GN(2 * sp.cos(p) / L, -1 / L**2)))

heading("5. The Turning Unit, shadow and height")
report("J is the arrow of length 1 straight up: page(J) = i", *page_eq(J.page(), I))
report("J^2 = -1", *same(J * J, GN(-1, 0)))
report("J = (G/|G| - cos(phibar)) / sin(phibar)", *same(J, GN(-sp.cos(p) / sp.sin(p), 0) + GN(0, 1 / L) * (1 / sp.sin(p))))
report("multiplying by J turns by a right angle: page(J z) = i page(z)", *page_eq((J * z1).page(), I * z1.page()))
report("z = Sh(z) + Ht(z) J", *same(z1, GN(z1.shadow(), 0) + J * z1.height()))
report("z - conj(z) = 2 Ht(z) J", *same(z1 - z1.conj(), J * (2 * z1.height())))
Sh_, Ht_ = sp.symbols("Sh Ht", real=True)
back = GN(Sh_ - Ht_ * sp.cos(p) / sp.sin(p), Ht_ / (L * sp.sin(p)))
report("back from the upright frame: a = Sh - Ht cot(phibar), b = Ht/(L sin(phibar))",
       *same(back, GN(Sh_, 0) + J * Ht_))
# the only square roots of -1
xs, ys = sp.symbols("xs ys", real=True)
w = GN(xs, ys)
sol = sp.solve([sp.expand((w * w).a + 1), sp.expand((w * w).b)], [xs, ys], dict=True)
roots_ok = len(sol) == 2 and all(
    zero(sv[xs] - sg * J.a)[0] and zero(sv[ys] - sg * J.b)[0] for sv, sg in
    zip(sorted(sol, key=lambda t: float(t[ys].subs({p: 1.2, L: 1.1}))), (-1, 1)))
report("+J and -J are the only General Numbers whose square is -1", roots_ok, "symbolic")

heading("6. The polar form with one angle")
report("page(E(alpha)) = cos(alpha) + i sin(alpha): the unit number at the angle alpha",
       *page_eq(E(al).page(), sp.cos(al) + I * sp.sin(al)))
report("general Pythagoras: C^2 + 2 cos(phibar) C S + S^2 = 1",
       *zero(Cg(al) ** 2 + 2 * sp.cos(p) * Cg(al) * Sg(al) + Sg(al) ** 2 - 1))
zr = GN(r * Cg(al), r * Sg(al) / L)
report("an arrow of length r at the angle alpha has a = r cos(alpha, phibar-alpha), bL = r sin(alpha, phibar-alpha)",
       *page_eq(zr.page(), r * (sp.cos(al) + I * sp.sin(al))))
report("its magnitude is r", *zero(zr.norm2() - r**2))

heading("7. Multiplying means turning")
report("E(alpha) E(beta) = E(alpha + beta)", *same(E(al) * E(be), E(al + be)))
report("cos(a+b, pb-a-b) = cos(a, pb-a) cos(b, pb-b) - sin(a, pb-a) sin(b, pb-b)",
       *zero(Cg(al + be) - (Cg(al) * Cg(be) - Sg(al) * Sg(be))))
report("sin(a+b, pb-a-b) = cos(a,pb-a) sin(b,pb-b) + sin(a,pb-a) cos(b,pb-b) + 2cos(pb) sin(a,pb-a) sin(b,pb-b)",
       *zero(Sg(al + be) - (Cg(al) * Sg(be) + Sg(al) * Cg(be) + 2 * sp.cos(p) * Sg(al) * Sg(be))))
report("quotient: E(alpha) / E(beta) = E(alpha - beta)", *same(E(al) / E(be), E(al - be)))
report("product in polar form: r E(alpha) times s E(beta) = r s E(alpha + beta)",
       *same((E(al) * r) * (E(be) * s_), E(al + be) * (r * s_)))
# The same formulas from the lines of Appendix A (Remark in Section 4.9.4)
phi0 = p - al - be  # the projection angle of the sum
report("Appendix A, projection angle + reciprocity: sin(a, pb-a-b) cos(b, pb-b) = sin(a, pb-a)",
       *zero(gsin(al, phi0) * gsin(p - be, be) - Sg(al)))
report("  and reciprocity: sin(pb-b, b) = cos(b, pb-b)", *zero(gsin(p - be, be) - Cg(be)))
report("Appendix A, projection angle + reciprocity: cos(a, pb-a-b) = cos(a, pb-a) - sin(a, pb-a) sin(b, pb-b)/cos(b, pb-b)",
       *zero(gcos(al, phi0) - (Cg(al) - Sg(al) * Sg(be) / Cg(be))))
report("Appendix A, studied angle (second line): sin(a+b, pb-a-b) = sin(a, pb-a-b) + cos(a, pb-a-b) sin(b, pb-b)",
       *zero(gsin(al + be, phi0) - (gsin(al, phi0) + gcos(al, phi0) * gsin(be, phi0 + al))))
report("Appendix A, direct consequence: cos(a+b, phi) = cos(a,phi) cos(b,phi) / (1 - sin(a,phi) sin(b,phi))",
       *zero(gcos(al + be, phi0) * (1 - gsin(al, phi0) * gsin(be, phi0)) - gcos(al, phi0) * gcos(be, phi0)))
report("the chain gives the sine law (with general Pythagoras)",
       *zero(Sg(al) / Cg(be) + (Cg(al) - Sg(al) * Sg(be) / Cg(be)) * Sg(be)
             - (Cg(al) * Sg(be) + Sg(al) * Cg(be) + 2 * sp.cos(p) * Sg(al) * Sg(be))))
P_, Q_ = Sg(al) * Sg(be), Cg(al) * Cg(be)
report("the chain gives the cosine law",
       *zero((Cg(al) - P_ / Cg(be)) * (Cg(be) - P_ / Cg(al)) / (1 - P_ / Q_) - (Q_ - P_)))
report("double angle: cos(2a, pb-2a) = cos(a, pb-a)^2 - sin(a, pb-a)^2",
       *zero(Cg(2 * al) - (Cg(al) ** 2 - Sg(al) ** 2)))
report("double angle: sin(2a, pb-2a) = 2 sin(a, pb-a) [cos(a, pb-a) + cos(pb) sin(a, pb-a)]",
       *zero(Sg(2 * al) - 2 * Sg(al) * (Cg(al) + sp.cos(p) * Sg(al))))

heading("8. Euler's formula and de Moivre's formula")
report("cos(theta) + J sin(theta) = cos(theta, pb-theta) + sin(theta, pb-theta) G/|G|",
       *same(GN(sp.cos(th), 0) + J * sp.sin(th), E(th)))
# d/dtheta along x + y = phibar, with the partial derivatives of Appendix A
X, Y = sp.symbols("X Y", real=True)
sin_x = sp.sin(Y) / sp.sin(X + Y) ** 2            # sin(x, y)'_x
sin_y = -sp.sin(X) * sp.cos(X + Y) / sp.sin(X + Y) ** 2
cos_x = -sp.sin(Y) * sp.cos(X + Y) / sp.sin(X + Y) ** 2
cos_y = sp.sin(X) / sp.sin(X + Y) ** 2
sub = {X: th, Y: p - th}
report("Appendix A partial derivatives are right (sin(x,y)'_x and '_y)",
       zero(sp.diff(gsin(X, Y), X) - sin_x)[0] and zero(sp.diff(gsin(X, Y), Y) - sin_y)[0], "symbolic")
report("Appendix A partial derivatives are right (cos(x,y)'_x and '_y)",
       zero(sp.diff(gcos(X, Y), X) - cos_x)[0] and zero(sp.diff(gcos(X, Y), Y) - cos_y)[0], "symbolic")
dE = GN((cos_x - cos_y).subs(sub), ((sin_x - sin_y).subs(sub)) / L)
report("d/dtheta E(theta) = J E(theta), from the Appendix A derivatives", *same(dE, J * E(th)))
report("E(0) = 1", *same(E(0), GN(1, 0)))
report("e^{pi J} + 1 = 0: E(pi) = -1", *same(E(sp.pi), GN(-1, 0)))
report("G = |G| e^{phibar J}: |G| E(phibar) = G", *same(E(p) * L, G))
report("J = e^{(pi/2) J}", *same(E(sp.pi / 2), J))
report("conj(E(alpha)) = E(-alpha)", *same(E(al).conj(), E(-al)))
for n in (2, 3, -1, -2):
    report(f"de Moivre, n = {n}: E(theta)^{n} = E({n} theta)", *same(E(th) ** n, E(n * th)))
# numeric de Moivre and the exponential series, at several planes
ok_dm, ok_series = True, True
for pv in (math.radians(79.86), math.pi / 3, 2.4):
    Lv = abs(math.cos(pv) + math.sin(pv))
    Uv = cmath.exp(1j * pv)
    Jv = (Uv - math.cos(pv)) / math.sin(pv)
    for n in range(-6, 10):
        for tv in (0.3, 1.7, -2.2):
            lhs = (math.sin(pv - tv) / math.sin(pv) + math.sin(tv) / math.sin(pv) * Uv) ** n
            rhs = math.sin(pv - n * tv) / math.sin(pv) + math.sin(n * tv) / math.sin(pv) * Uv
            ok_dm &= abs(lhs - rhs) < 1e-9
    for tv in (0.5, 2.0, -1.3):
        ser = sum((tv * Jv) ** k / math.factorial(k) for k in range(60))
        ok_series &= abs(ser - (math.sin(pv - tv) / math.sin(pv) + math.sin(tv) / math.sin(pv) * Uv)) < 1e-12
report("de Moivre for n = -6 ... 9 on the page, at 79.86 deg, pi/3 and 2.4 rad", ok_dm, "numeric")
report("the exponential series of theta J sums to E(theta)", ok_series, "numeric")

heading("9. The powers of the Composer Unit in its own frame")
for n in range(-3, 9):
    if n == 0:
        continue
    rhs = GN(L**n * Cg(n * p), L ** (n - 1) * Sg(n * p))
    report(f"G^{n} = L^{n} cos({n}pb, {1-n}pb) + L^{n-1} sin({n}pb, {1-n}pb) G", *same(G ** n, rhs))
n_ = 5
report("step rule: G^(n+1) = 2 L cos(phibar) G^n - L^2 G^(n-1)  (n = 5)",
       *same(G ** (n_ + 1), G ** n_ * (2 * L * sp.cos(p)) - G ** (n_ - 1) * L**2))
report("n = 2 is the Rule of the Square", *same(G ** 2, GN(-L**2, 2 * L * sp.cos(p))))
report("on the page, G^n is the arrow L^n at the angle n phibar (Theorem 4.1), n = 4",
       *page_eq((G ** 4).page(), L**4 * (sp.cos(4 * p) + I * sp.sin(4 * p))))
cheb = all(zero(sp.sin(k * p) / sp.sin(p) - sp.chebyshevu(k - 1, sp.cos(p)))[0] for k in range(1, 9))
report("sin(n phibar)/sin(phibar) is the Chebyshev polynomial U_{n-1}(cos phibar), n = 1 ... 8", cheb, "symbolic")

heading("10. Roots")
rnd = random.Random(5)
ok_roots = True
for pv in (math.radians(79.86), math.pi / 3, 1.0):
    for n in (2, 3, 5, 7):
        zv = complex(rnd.uniform(-2, 2), rnd.uniform(-2, 2))
        Mz, az = abs(zv), cmath.phase(zv)
        ws = []
        for k in range(n):
            tk = (az + 2 * math.pi * k) / n
            wv = Mz ** (1 / n) * (math.sin(pv - tk) / math.sin(pv) + math.sin(tk) / math.sin(pv) * cmath.exp(1j * pv))
            ok_roots &= abs(wv**n - zv) < 1e-9
            ws.append(wv)
        ok_roots &= min(abs(ws[i] - ws[j]) for i in range(n) for j in range(i)) > 1e-6
report("the n roots w_k = |z|^(1/n) E((alpha + 2 pi k)/n) all satisfy w^n = z and are distinct", ok_roots, "numeric")
ok_unity = True
for pv in (math.radians(79.86), math.pi / 3):
    for n in (3, 4, 6, 8):
        for k in range(n):
            tk = 2 * math.pi * k / n
            Cv, Sv = math.sin(pv - tk) / math.sin(pv), math.sin(tk) / math.sin(pv)
            ok_unity &= abs((Cv + Sv * cmath.exp(1j * pv)) ** n - 1) < 1e-9
            ok_unity &= abs(Cv**2 + 2 * math.cos(pv) * Cv * Sv + Sv**2 - 1) < 1e-12
report("roots of unity E(2 pi k / n) satisfy z^n = 1 and lie on the unit ring", ok_unity, "numeric")
sqrtG = GN(sp.sqrt(L) / (2 * sp.cos(p / 2)), sp.sqrt(L) / (2 * sp.cos(p / 2)) / L)
report("square root of G: [sqrt(L) / (2 cos(phibar/2))] (1 + G/|G|), squared, is G", *same(sqrtG * sqrtG, G))
report("  and it is E(phibar/2) sqrt(L), on the ray of equal parts", *same(sqrtG, E(p / 2) * sp.sqrt(L)))
report("  Appendix A isosceles case: sin(t, t) = cos(t, t) = 1/(2 cos t)",
       zero(gsin(th, th) - 1 / (2 * sp.cos(th)))[0] and zero(gcos(th, th) - 1 / (2 * sp.cos(th)))[0], "symbolic")
report("(1 + G/|G|)^2 = 4 cos^2(phibar/2) G/|G|", *same((GN(1, 0) + u) ** 2, u * (4 * sp.cos(p / 2) ** 2)))

heading("11. Exponential, logarithm, matrix, change of plane, second mirror (Appendix F)")
ex_rhs = E(z1.height()) * sp.exp(z1.shadow())
report("e^z = e^{Sh(z)} [cos(Ht, pb-Ht) + sin(Ht, pb-Ht) G/|G|]  (on the page)",
       *page_eq(ex_rhs.page(), sp.exp(z1.page())))
lg = GN(sp.log(r) - al * sp.cos(p) / sp.sin(p), al / (L * sp.sin(p)))
report("log of r e^{alpha J}: (ln r - alpha cot(phibar)) + alpha/(L sin(phibar)) G  =  ln r + alpha J",
       *same(lg, GN(sp.log(r), 0) + J * al))
report("  and its page value is ln r + i alpha", *page_eq(lg.page(), sp.log(r) + I * al))


def Mat(z):
    return sp.Matrix([[z.a, -z.L**2 * z.b], [z.b, z.a + 2 * z.L * z.b * sp.cos(z.p)]])


report("matrix of a product is the product of the matrices",
       all(zero(t)[0] for t in (Mat(z1) * Mat(z2) - Mat(z1 * z2))), "symbolic")
report("det = |z|^2 and trace = 2 Sh(z)",
       zero(Mat(z1).det() - z1.norm2())[0] and zero(Mat(z1).trace() - 2 * z1.shadow())[0], "symbolic")
MG = Mat(G)
report("the matrix of G obeys the Rule of the Square: M^2 = 2L cos(phibar) M - L^2",
       all(zero(t)[0] for t in (MG * MG - 2 * L * sp.cos(p) * MG + L**2 * sp.eye(2))), "symbolic")
bq = b * L * sp.sin(p) / (K * sp.sin(q))
aq = a + b * L * sp.cos(p) - bq * K * sp.cos(q)
report("moving a number to another relative plane: same shadow, same height",
       *page_eq(aq + bq * K * (sp.cos(q) + I * sp.sin(q)), z1.page()))
report("from the complex plane: x + iy = (x - y cot(phibar)) + y/(L sin(phibar)) G",
       *page_eq(GN(x - y * sp.cos(p) / sp.sin(p), y / (L * sp.sin(p))).page(), x + I * y))
uu = sp.cos(p) + I * sp.sin(p)
report("second mirror: z -> (G/|G|) conj(z) swaps the two unit coordinates",
       *page_eq(uu * sp.conjugate(x + y * uu), y + x * uu))
report("reciprocity in the frame: E(phibar - theta) = sin(theta, pb-theta) + cos(theta, pb-theta) G/|G|",
       *same(E(p - th), GN(Sg(th), Cg(th) / L)))

heading("11b. More laws listed in Appendix F")
report("negative angle: cos(-t, pb+t) = cos(t, pb-t) + 2cos(pb) sin(t, pb-t)",
       *zero(Cg(-th) - (Cg(th) + 2 * sp.cos(p) * Sg(th))))
report("negative angle: sin(-t, pb+t) = -sin(t, pb-t)", *zero(Sg(-th) + Sg(th)))
report("subtraction: cos(a-b, pb-a+b) = C(a)C(b) + S(a)S(b) + 2cos(pb) C(a)S(b)",
       *zero(Cg(al - be) - (Cg(al) * Cg(be) + Sg(al) * Sg(be) + 2 * sp.cos(p) * Cg(al) * Sg(be))))
report("subtraction: sin(a-b, pb-a+b) = S(a)C(b) - C(a)S(b)",
       *zero(Sg(al - be) - (Sg(al) * Cg(be) - Cg(al) * Sg(be))))
report("multiplying by G/|G| turns by phibar: C(t + pb) = -S(t), S(t + pb) = C(t) + 2cos(pb) S(t)",
       zero(Cg(th + p) + Sg(th))[0] and zero(Sg(th + p) - (Cg(th) + 2 * sp.cos(p) * Sg(th)))[0], "symbolic")
report("multiplying by J turns by pi/2: C(t + pi/2) = -C cot(pb) - S/sin(pb), S(t + pi/2) = C/sin(pb) + S cot(pb)",
       zero(Cg(th + sp.pi / 2) - (-Cg(th) * sp.cos(p) / sp.sin(p) - Sg(th) / sp.sin(p)))[0]
       and zero(Sg(th + sp.pi / 2) - (Cg(th) / sp.sin(p) + Sg(th) * sp.cos(p) / sp.sin(p)))[0], "symbolic")
report("conjugate in polar form: conj(r E(t)) = r [cos(-t, pb+t) + sin(-t, pb+t) G/|G|]",
       *same((E(th) * r).conj(), E(-th) * r))
report("angle: tan(theta) = Ht(z)/Sh(z) (Theorem 4.10 with eta = |G| b/a)",
       *zero(sp.tan(sp.atan2(z1.height(), z1.shadow())) - z1.height() / z1.shadow()))
report("Theorem 4.10: eta sin(pb)/(1 + eta cos(pb)) = Ht/Sh with eta = |G| b / a",
       *zero((L * b / a) * sp.sin(p) / (1 + (L * b / a) * sp.cos(p)) - z1.height() / z1.shadow()))
report("Chebyshev: sin(2pb)/sin(pb) = 2cos(pb), sin(3pb)/sin(pb) = 4cos^2(pb) - 1, sin(4pb)/sin(pb) = 8cos^3 - 4cos",
       zero(sp.sin(2 * p) / sp.sin(p) - 2 * sp.cos(p))[0]
       and zero(sp.sin(3 * p) / sp.sin(p) - (4 * sp.cos(p) ** 2 - 1))[0]
       and zero(sp.sin(4 * p) / sp.sin(p) - (8 * sp.cos(p) ** 3 - 4 * sp.cos(p)))[0], "symbolic")
report("multiplying by J: page(J z) = i page(z), and J z = -(a cos(pb) + b|G|)/sin(pb) + Sh(z)/(|G| sin(pb)) G",
       *same(J * z1, GN(-(a * sp.cos(p) + b * L) / sp.sin(p), z1.shadow() / (L * sp.sin(p)))))
ok_quad = True
rq = random.Random(9)
for pv_ in (math.radians(79.86), math.pi / 3, 2.0):
    Uq = cmath.exp(1j * pv_)
    for _ in range(50):
        Pq = complex(rq.uniform(-3, 3), rq.uniform(-3, 3))
        Qq = complex(rq.uniform(-3, 3), rq.uniform(-3, 3))
        Dq = Pq * Pq - 4 * Qq
        tq = cmath.phase(Dq) / 2       # the square root of Theorem (Roots), k = 0
        sq_ = math.sqrt(abs(Dq)) * (math.sin(pv_ - tq) / math.sin(pv_) + math.sin(tq) / math.sin(pv_) * Uq)
        for sg in (1, -1):
            wq = (-Pq + sg * sq_) / 2
            ok_quad &= abs(wq * wq + Pq * wq + Qq) < 1e-9
report("quadratic formula: w^2 + P w + Q = 0 has the roots (-P +/- sqrt(P^2 - 4Q))/2, square root from the roots theorem",
       ok_quad, "numeric")

heading("12. At phibar = pi/2 every law is the law of Appendix C")
hp = {p: sp.pi / 2, L: 1}


def at90(z):
    return GN(sp.simplify(z.a.subs(hp)), sp.simplify(z.b.subs(hp)), sp.pi / 2, 1)


def c_eq(z, re_, im_):
    return zero(at90(z).a - re_)[0] and zero(at90(z).b - im_)[0]


report("the Axiom gives L = |cos + sin| = 1 at pi/2", sp.Abs(sp.cos(sp.pi / 2) + sp.sin(sp.pi / 2)) == 1, "symbolic")
report("product -> (ac - bd) + (ad + bc) i", c_eq(z1 * z2, a * c - b * d, a * d + b * c), "symbolic")
report("conjugate -> a - bi", c_eq(z1.conj(), a, -b), "symbolic")
report("z conj(z) -> a^2 + b^2", c_eq(z1 * z1.conj(), a**2 + b**2, 0), "symbolic")
report("division -> [(ac + bd) + (bc - ad) i] / (c^2 + d^2)",
       c_eq(quot, (a * c + b * d) / (c**2 + d**2), (b * c - a * d) / (c**2 + d**2)), "symbolic")
report("two-squares identity -> (a^2+b^2)(c^2+d^2) = (ac-bd)^2 + (ad+bc)^2",
       zero(((z1 * z2).norm2() - z1.norm2() * z2.norm2()).subs(hp))[0]
       and zero((z1.norm2() * z2.norm2()).subs(hp) - (a**2 + b**2) * (c**2 + d**2))[0], "symbolic")
report("Turning Unit -> i, and G/|G| -> i", c_eq(J, 0, 1) and c_eq(u, 0, 1), "symbolic")
report("shadow and height -> real and imaginary parts",
       zero(z1.shadow().subs(hp) - a)[0] and zero(z1.height().subs(hp) - b)[0], "symbolic")
report("unit number -> cos(theta) + i sin(theta)", c_eq(E(th), sp.cos(th), sp.sin(th)), "symbolic")
report("addition formulas -> cos(a+b) = cos a cos b - sin a sin b, sin(a+b) = sin a cos b + cos a sin b",
       zero((Cg(al) * Cg(be) - Sg(al) * Sg(be)).subs(hp) - sp.cos(al + be))[0]
       and zero((Cg(al) * Sg(be) + Sg(al) * Cg(be) + 2 * sp.cos(p) * Sg(al) * Sg(be)).subs(hp) - sp.sin(al + be))[0],
       "symbolic")
report("general Pythagoras -> cos^2 + sin^2 = 1",
       zero((Cg(al) ** 2 + 2 * sp.cos(p) * Cg(al) * Sg(al) + Sg(al) ** 2).subs(hp) - 1)[0], "symbolic")
report("powers of G -> i^n = cos(n pi/2) + i sin(n pi/2), n = 1 ... 8",
       all(c_eq(GN(L**n * Cg(n * p), L ** (n - 1) * Sg(n * p)), sp.cos(n * sp.pi / 2), sp.sin(n * sp.pi / 2))
           for n in range(1, 9)), "symbolic")
report("step rule -> i^(n+1) = -i^(n-1)", zero((2 * L * sp.cos(p)).subs(hp))[0] and (L**2).subs(hp) == 1, "symbolic")
report("square root of G -> (1 + i)/sqrt(2)", c_eq(sqrtG, 1 / sp.sqrt(2), 1 / sp.sqrt(2)), "symbolic")
report("exponential -> e^a (cos b + i sin b)",
       zero(z1.shadow().subs(hp) - a)[0] and zero(z1.height().subs(hp) - b)[0], "symbolic")
report("logarithm -> ln r + i alpha", c_eq(lg, sp.log(r), al), "symbolic")
report("matrix -> [[a, -b], [b, a]]", all(zero(t)[0] for t in (Mat(z1).subs(hp) - sp.Matrix([[a, -b], [b, a]]))), "symbolic")
report("from the complex plane -> a = x, b = y", c_eq(GN(x - y * sp.cos(p) / sp.sin(p), y / (L * sp.sin(p))), x, y), "symbolic")

heading("13. The formulas of Appendix A that the section uses (random tests)")
rnd = random.Random(3)


def f_s(al_, ph_):
    return math.sin(al_) / math.sin(al_ + ph_)


def f_c(al_, ph_):
    return math.sin(ph_) / math.sin(al_ + ph_)


def close(u1, u2):
    return abs(u1 - u2) <= 1e-8 * max(1.0, abs(u1), abs(u2))


tests = {
    "reciprocity sin(a,f) = cos(f,a)": lambda A, B, F, T: close(f_s(A, F), f_c(F, A)),
    "special case a + f = pi/2": lambda A, B, F, T: close(f_s(A, math.pi / 2 - A), math.sin(A))
    and close(f_c(A, math.pi / 2 - A), math.cos(A)),
    "isosceles sin(t,t) = cos(t,t) = sec(t)/2": lambda A, B, F, T: close(f_s(A, A), 0.5 / math.cos(A)),
    "pi-addition sin(a+p pi, f+q pi) = (-1)^q sin(a,f)": lambda A, B, F, T: close(f_s(A + math.pi, F - math.pi), -f_s(A, F))
    and close(f_s(A + 2 * math.pi, F - math.pi), -f_s(A, F)),
    "pi-addition cos(a+p pi, f+q pi) = (-1)^p cos(a,f)": lambda A, B, F, T: close(f_c(A + math.pi, F - math.pi), -f_c(A, F))
    and close(f_c(A + 2 * math.pi, F - 2 * math.pi), f_c(A, F)),
    "studied-angle addition, line 1": lambda A, B, F, T: close(f_s(A + B, F), f_s(A, F) + f_c(A, F) * f_s(B, F) * f_c(A, F + B)),
    "studied-angle addition, line 2": lambda A, B, F, T: close(f_s(A + B, F), f_s(A, F) + f_c(A, F) * f_s(B, F + A)),
    "studied-angle addition, line 3": lambda A, B, F, T: close(f_c(A + B, F), f_c(A, F) * f_c(B, F) + f_c(A, F) * f_s(B, F) * f_s(A, F + B)),
    "studied-angle addition, line 4": lambda A, B, F, T: close(f_c(A + B, F), f_c(A, F) * f_c(B, F) + f_s(A, F) * f_s(B, F) * f_c(A + B, F)),
    "projection-angle addition, sine": lambda A, B, F, T: close(f_s(A, F + T), f_s(A, F) * f_s(A + F, T)),
    "projection-angle addition, cosine": lambda A, B, F, T: close(f_c(A, F + T), f_c(A, F) + f_s(A, F) * f_c(A + F, T)),
    "consequence cos(a+b,f) = cos cos / (1 - sin sin)": lambda A, B, F, T: close(f_c(A + B, F), f_c(A, F) * f_c(B, F) / (1 - f_s(A, F) * f_s(B, F))),
    "double angle, sine": lambda A, B, F, T: close(f_s(2 * A, F), f_s(A, F) * (f_c(2 * A, F) + 1)),
    "double angle, cosine (solved form)": lambda A, B, F, T: close(f_c(2 * A, F), f_c(A, F) ** 2 / (1 - f_s(A, F) ** 2)),
}
for name, fn in tests.items():
    ok, n_ok = True, 0
    while n_ok < 3000:
        A, B, F, T = (rnd.uniform(-3, 3) for _ in range(4))
        args = [A + F, A + B + F, 2 * A + F, A + F + T, B + F, A + B + 2 * F]
        if min(abs(math.sin(t)) for t in args) < 0.05 or abs(math.cos(A)) < 0.05 or abs(math.sin(F)) < 0.05:
            continue
        if abs(1 - f_s(A, F) * f_s(B, F)) < 0.05 or abs(1 - f_s(A, F) ** 2) < 0.05:
            continue
        ok &= fn(A, B, F, T)
        n_ok += 1
    report("Appendix A: " + name, ok, "numeric, 3000 points")

# --------------------------------------------------------------------------
heading("14. Numbers of the worked examples, at phibar = 79.86 degrees")
pv = math.radians(79.86)
cv, sv = math.cos(pv), math.sin(pv)
Lv = abs(cv + sv)
f5 = lambda t: f"{t:.5f}"
print(f"  cos = {f5(cv)}, sin = {f5(sv)}, |G| = {f5(Lv)}, 2|G|cos = {f5(2*Lv*cv)}, |G|^2 = {f5(Lv**2)}")
print(f"  |G| cos = {f5(Lv*cv)}, |G| sin = {f5(Lv*sv)}, cot = {f5(cv/sv)}, 1/(|G| sin) = {f5(1/(Lv*sv))}")


class N:
    """Numeric General Number at 79.86 degrees."""

    def __init__(self, a_, b_):
        self.a, self.b = a_, b_

    def __mul__(self, o):
        return N(self.a * o.a - Lv**2 * self.b * o.b, self.a * o.b + self.b * o.a + 2 * Lv * cv * self.b * o.b)

    def conj(self):
        return N(self.a + 2 * Lv * self.b * cv, -self.b)

    def n2(self):
        return self.a**2 + 2 * Lv * self.a * self.b * cv + Lv**2 * self.b**2

    def page(self):
        return complex(self.a + self.b * Lv * cv, self.b * Lv * sv)

    def __repr__(self):
        return f"{f5(self.a)} + ({f5(self.b)}) G"


zA = N(1, 1)
print("\n  Worked example: the conjugate, the magnitude and the inverse of z = 1 + G")
print("    page of z      :", f5(zA.page().real), f5(zA.page().imag))
print("    conjugate      :", zA.conj(), "  page:", f5(zA.conj().page().real), f5(zA.conj().page().imag))
print("    |z|^2 = z conj :", f5(zA.n2()), " |z| =", f5(math.sqrt(zA.n2())), " check on page:", f5(abs(zA.page()) ** 2))
inv = N(zA.conj().a / zA.n2(), zA.conj().b / zA.n2())
print("    1/z            :", inv, " check z (1/z):", zA * inv)
sk = N(1, -1)
print("    skew flip 1 - G: (1+G)(1-G) =", zA * sk, " |1-G| =", f5(math.sqrt(sk.n2())))
zB, zC = N(2, -1), N(1, 1)
num = zB * zC.conj()
qq = N(num.a / zC.n2(), num.b / zC.n2())
print("    (2 - G)/(1 + G): numerator", num, " denominator", f5(zC.n2()), " quotient", qq, " check q(1+G):", qq * zC)
print("    on the page    :", f5((zB.page() / zC.page()).real), f5((zB.page() / zC.page()).imag), " vs", f5(qq.page().real), f5(qq.page().imag))

print("\n  Worked example: the Turning Unit, and 1 + G in the upright frame")
JN = N(-cv / sv, 1 / (Lv * sv))
print("    J              :", JN, " J^2 =", JN * JN, " page of J:", f5(JN.page().real), f5(JN.page().imag))
Sh_A, Ht_A = zA.a + zA.b * Lv * cv, zA.b * Lv * sv
print("    1 + G = Sh + Ht J with Sh =", f5(Sh_A), " Ht =", f5(Ht_A))
chk = N(Sh_A, 0)
tmp = N(Ht_A * JN.a, Ht_A * JN.b)
print("    check Sh + Ht J:", N(chk.a + tmp.a, chk.b + tmp.b))
print("    J (1 + G)      :", JN * zA, " page:", f5((JN * zA).page().real), f5((JN * zA).page().imag), " = i page(z):", f5((1j * zA.page()).real), f5((1j * zA.page()).imag))

print("\n  Worked example: (1 + G) G, two ways")
GN_ = N(0, 1)
pr = zA * GN_
print("    by Proposition 4.23:", pr)
m1, a1_ = math.sqrt(zA.n2()), math.atan2(zA.page().imag, zA.page().real)
m2, a2_ = Lv, pv
mp, ap = m1 * m2, a1_ + a2_
print("    |1+G| =", f5(m1), " angle =", f"{math.degrees(a1_):.3f}", " ; |G| =", f5(m2), " angle =", f"{math.degrees(a2_):.3f}")
print("    product: length", f5(mp), " angle", f"{math.degrees(ap):.3f}")
Cp, Sp = math.sin(pv - ap) / sv, math.sin(ap) / sv
print("    cos(ap, pb-ap) =", f5(Cp), " sin(ap, pb-ap) =", f5(Sp))
print("    real part =", f5(mp * Cp), "  Composer coordinate =", f5(mp * Sp / Lv))
# also the General Sinus values of the two factors
C1, S1 = math.sin(pv - a1_) / sv, math.sin(a1_) / sv
print("    1+G: cos(a1, pb-a1) =", f5(C1), " sin(a1, pb-a1) =", f5(S1), " -> |z| cos =", f5(m1 * C1), " |z| sin / |G| =", f5(m1 * S1 / Lv))

print("\n  Worked example: G^3, three ways")
G3 = GN_ * GN_ * GN_
print("    by multiplying :", G3)
print("    by the formula : L^3 cos(3pb, -2pb) =", f5(Lv**3 * math.sin(pv - 3 * pv) / sv), "; L^2 sin(3pb, -2pb) =", f5(Lv**2 * math.sin(3 * pv) / sv))
G2 = GN_ * GN_
print("    G^2            :", G2)
print("    step rule      : 2|G|cos G^2 - |G|^2 G =", N(2 * Lv * cv * G2.a, 2 * Lv * cv * G2.b - Lv**2))
print("    on the page    : L^3 =", f5(Lv**3), " at", f"{3*79.86:.2f}", "deg:", f5(Lv**3 * math.cos(3 * pv)), f5(Lv**3 * math.sin(3 * pv)))
print("                     read back: b =", f5(Lv**3 * math.sin(3 * pv) / (Lv * sv)), " a =", f5(Lv**3 * math.cos(3 * pv) - Lv**3 * math.sin(3 * pv) / (Lv * sv) * Lv * cv))
print("    sin(3pb)/sin(pb) =", f5(math.sin(3 * pv) / sv), " = 3 - 4 sin^2 =", f5(3 - 4 * sv**2), " ; sin(2pb)/sin(pb) = 2cos =", f5(2 * cv))

print("\n  Worked example: the square root of G and the sixth roots of unity")
hv = pv / 2
k_ = math.sqrt(Lv) / (2 * math.cos(hv))
sq = N(k_, k_ / Lv)
print("    sqrt(|G|) =", f5(math.sqrt(Lv)), " cos(pb/2) =", f5(math.cos(hv)), " factor =", f5(k_))
print("    sqrt(G) = +/-", sq, "  squared:", sq * sq)
print("    angle of sqrt(G):", f"{math.degrees(math.atan2(sq.page().imag, sq.page().real)):.3f}", "deg; length", f5(abs(sq.page())))
print("    sixth roots of unity (k, angle, real part, Composer coordinate, page x, page y):")
for k in range(6):
    tk = 2 * math.pi * k / 6
    Cv, Sv = math.sin(pv - tk) / sv, math.sin(tk) / sv
    zk = N(Cv, Sv / Lv)
    pg = zk.page()
    zk6 = N(1, 0)
    for _ in range(6):
        zk6 = zk6 * zk
    print(f"      {k}  {60*k:3d}  {f5(Cv):>9}  {f5(Sv/Lv):>9}   page {f5(pg.real):>9} {f5(pg.imag):>9}   z^6 = {zk6}")

# --------------------------------------------------------------------------
heading("15. The parts beside the shadow and the height (Subsection 4.9.7 and Appendix F.8)")
k_ = sp.symbols("k", real=True)


def Sh(z):
    return z.shadow()


def Ht(z):
    return z.height()


def com(z):
    """The Composer part b G of z, itself a General Number (Definition 4.11)."""
    return GN(0, z.b, z.p, z.L)


def allz(*exprs):
    res = [zero(e) for e in exprs]
    how = "symbolic" if all(h == "symbolic" for _, h in res) else "numeric only"
    return all(ok for ok, _ in res), how


ct = sp.cos(p) / sp.sin(p)
S1, S2, H1, H2 = Sh(z1), Sh(z2), Ht(z1), Ht(z2)
n1, n2 = z1.norm2(), z2.norm2()
Gb = G.conj()                        # the conjugate of the Composer Unit
zc = z1.conj()

# Proposition 4.54, From One Pair to the Other
report("Prop 4.54 (1): Re(z) = Sh(z) - Ht(z) cot(pb), Com(z) = Ht(z)/(|G| sin pb) G",
       *allz(z1.a - (S1 - H1 * ct), z1.b - H1 / (L * sp.sin(p))))
report("Prop 4.54 (1): Sh(z) = Re(z) + b|G| cos(pb) = Re(z) + Sh(Com z), Ht(z) = b|G| sin(pb) = Ht(Com z)",
       *allz(S1 - (z1.a + b * L * sp.cos(p)), S1 - (z1.a + Sh(com(z1))),
             H1 - b * L * sp.sin(p), H1 - Ht(com(z1))))
report("Prop 4.54 (2): Re, Com, Sh and Ht are additive and keep real factors",
       *allz(*[e for f in (lambda t: t.a, lambda t: t.b, Sh, Ht)
               for e in (f(z1 + z2) - f(z1) - f(z2), f(z1 * k_) - k_ * f(z1))]))
report("Prop 4.54 (3): the conjugate keeps the coordinates and mirrors the unit, conj(z) = a + b conj(G)",
       *same(zc, GN(a, 0) + Gb * b))
report("Prop 4.54 (3): G - conj(G) = 2 |G| sin(pb) J", *same(G - Gb, J * (2 * L * sp.sin(p))))
report("Prop 4.54 (3): Re(z) = (G conj(z) - conj(G) z) / (G - conj(G))", *same((G * zc - Gb * z1) / (G - Gb), GN(a, 0)))
report("Prop 4.54 (3): Com(z) = G (z - conj(z)) / (G - conj(G))", *same(G * (z1 - zc) / (G - Gb), com(z1)))
report("Prop 4.54 (3): Sh(z) = (z + conj z)/2 and Ht(z) J = (z - conj z)/2",
       same((z1 + zc) * sp.Rational(1, 2), GN(S1, 0))[0] and same((z1 - zc) * sp.Rational(1, 2), J * H1)[0], "symbolic")
report("Prop 4.54 (4): z - conj(z) = b (G - conj G), zero exactly when b = 0", *same(z1 - zc, (G - Gb) * b))
report("Prop 4.54 (4): G conj(z) - conj(G) z = a (G - conj G), zero exactly when Re(z) = 0", *same(G * zc - Gb * z1, (G - Gb) * a))
report("F.8: a number t J on the upright axis has conj = -(t J), and z + conj z = 2 Sh(z)",
       same((J * x).conj(), -(J * x))[0] and same(z1 + zc, GN(2 * S1, 0))[0], "symbolic")

# Theorem 4.55, The Parts of a Product
pr12 = z1 * z2
report("Thm 4.55 (1): Re(z1 z2) = Re(z1) Re(z2) - |G|^2 b1 b2", *zero(pr12.a - (a * c - L**2 * b * d)))
report("Thm 4.55 (1): Com(z1 z2) = Sh(z1) Com(z2) + Sh(z2) Com(z1)", *same(com(pr12), com(z2) * S1 + com(z1) * S2))
report("Thm 4.55 (2): Sh(z1 z2) = Sh(z1) Sh(z2) - Ht(z1) Ht(z2), Ht(z1 z2) = Sh(z1) Ht(z2) + Sh(z2) Ht(z1)",
       *allz(Sh(pr12) - (S1 * S2 - H1 * H2), Ht(pr12) - (S1 * H2 + S2 * H1)))
report("Thm 4.55 (2) on the page: Sh + i Ht of the product is the product of the arrows",
       *page_eq(Sh(pr12) + I * Ht(pr12), z1.page() * z2.page()))
sq1 = z1 * z1
report("Thm 4.55 (3): Re(z^2) = Re(z)^2 - |G|^2 b^2, Com(z^2) = 2 Sh(z) Com(z)",
       zero(sq1.a - (a**2 - L**2 * b**2))[0] and same(com(sq1), com(z1) * (2 * S1))[0], "symbolic")
report("Thm 4.55 (3): Sh(z^2) = Sh(z)^2 - Ht(z)^2, Ht(z^2) = 2 Sh(z) Ht(z)",
       *allz(Sh(sq1) - (S1**2 - H1**2), Ht(sq1) - 2 * S1 * H1))

# Proposition 4.56, The Parts of the Conjugate, the Inverse and the Quotient
report("Prop 4.56 (1): Re(conj z) = 2 Sh(z) - Re(z), Com(conj z) = -Com(z), Sh(conj z) = Sh(z), Ht(conj z) = -Ht(z)",
       *allz(zc.a - (2 * S1 - a), zc.b + b, Sh(zc) - S1, Ht(zc) + H1))
report("Prop 4.56 (1): the shadow is the midpoint of Re(z) and Re(conj z)", *zero(S1 - (z1.a + zc.a) / 2))
report("Prop 4.56 (1) on the page: z is as far from Re(z) as from Re(conj z), |b| |G| (isosceles triangle)",
       *allz((z1 - GN(z1.a, 0)).norm2() - b**2 * L**2, (z1 - GN(zc.a, 0)).norm2() - b**2 * L**2))
zi = z1.inv()
report("Prop 4.56 (2): Re(1/z) = (2 Sh - Re)/|z|^2, Com(1/z) = -Com(z)/|z|^2, Sh(1/z) = Sh/|z|^2, Ht(1/z) = -Ht/|z|^2",
       *allz(zi.a - (2 * S1 - a) / n1, zi.b + b / n1, Sh(zi) - S1 / n1, Ht(zi) + H1 / n1))
qz, z2c = z1 / z2, z2.conj()
report("Prop 4.56 (3): Re(z1/z2) = [Re(z1) Re(conj z2) + |G|^2 b1 b2]/|z2|^2",
       *zero(qz.a - (a * z2c.a + L**2 * b * d) / n2))
report("Prop 4.56 (3): Com(z1/z2) = [Sh(z2) Com(z1) - Sh(z1) Com(z2)]/|z2|^2", *zero(qz.b - (S2 * b - S1 * d) / n2))
report("Prop 4.56 (3): Sh(z1/z2) = [Sh1 Sh2 + Ht1 Ht2]/|z2|^2, Ht(z1/z2) = [Sh2 Ht1 - Sh1 Ht2]/|z2|^2",
       *allz(Sh(qz) - (S1 * S2 + H1 * H2) / n2, Ht(qz) - (S2 * H1 - S1 * H2) / n2))
report("F.8: the top z1 conj(z2) of the quotient has these parts and readings",
       *allz((z1 * z2c).a - (a * z2c.a + L**2 * b * d), (z1 * z2c).b - (S2 * b - S1 * d),
             Sh(z1 * z2c) - (S1 * S2 + H1 * H2), Ht(z1 * z2c) - (S2 * H1 - S1 * H2)))

# Corollary 4.57, Turning by the Two Units
Jz, Gz = J * z1, G * z1
report("Cor 4.57 (1): Re(J z) = -Ht(z) - Sh(z) cot(pb), Com(J z) = Sh(z)/(|G| sin pb) G, Sh(J z) = -Ht(z), Ht(J z) = Sh(z)",
       *allz(Jz.a - (-H1 - S1 * ct), Jz.b - S1 / (L * sp.sin(p)), Sh(Jz) + H1, Ht(Jz) - S1))
report("Cor 4.57 (2): Re(G z) = -|G|^2 b, Com(G z) = Re(conj z) G",
       *allz(Gz.a + L**2 * b, Gz.b - zc.a))
report("Cor 4.57 (2): Sh(G z) = |G|(Sh cos pb - Ht sin pb), Ht(G z) = |G|(Sh sin pb + Ht cos pb)",
       *allz(Sh(Gz) - L * (S1 * sp.cos(p) - H1 * sp.sin(p)), Ht(Gz) - L * (S1 * sp.sin(p) + H1 * sp.cos(p))))
report("Cor 4.57: G z is the second column of the matrix of z (Appendix F.7)",
       *allz(Mat(z1)[0, 1] - Gz.a, Mat(z1)[1, 1] - Gz.b))

# Proposition 4.58, The Signed Area
A_ = zc * z2
report("Prop 4.58: Com(conj(z1) z2) = (a1 b2 - b1 a2) G", *zero(A_.b - (a * d - b * c)))
report("Prop 4.58: Ht(conj(z1) z2) = Sh1 Ht2 - Ht1 Sh2 = (a1 b2 - b1 a2) |G| sin(pb)",
       *allz(Ht(A_) - (S1 * H2 - H1 * S2), Ht(A_) - (a * d - b * c) * L * sp.sin(p)))
report("Prop 4.58: Re(conj(z1) z2) = Re(conj z1) Re(z2) + |G|^2 b1 b2, Sh(conj(z1) z2) = Sh1 Sh2 + Ht1 Ht2",
       *allz(A_.a - (zc.a * c + L**2 * b * d), Sh(A_) - (S1 * S2 + H1 * H2)))
report("Prop 4.58 on the page: Sh + i Ht of conj(z1) z2 is conj(arrow 1) times arrow 2 (dot product + i cross product)",
       *page_eq(Sh(A_) + I * Ht(A_), sp.conjugate(z1.page()) * z2.page()))
report("Prop 4.58: two numbers on one line through 0, z2 = t z1, have a1 b2 - b1 a2 = 0", *zero((zc * (z1 * x)).b))
report("two laws carry neither cos(pb) nor sin(pb): Re(z1 z2) and Com(conj(z1) z2)",
       p not in sp.expand(pr12.a).free_symbols and p not in sp.expand(A_.b).free_symbols, "symbolic")

# Corollary 4.59, Two Trigonometries
zp = E(al) * r
ok_pw, how_pw = True, "symbolic"
for n in (-2, -1, 1, 2, 3, 4):
    zn = zp ** n
    ok_n, h_n = allz(zn.a - r**n * Cg(n * al), zn.b * L - r**n * Sg(n * al),
                     Sh(zn) - r**n * sp.cos(n * al), Ht(zn) - r**n * sp.sin(n * al))
    ok_pw &= ok_n
    how_pw = how_pw if h_n == "symbolic" else "numeric only"
report("Cor 4.59, n = -2 ... 4: Re(z^n) = r^n C(n alpha), Com(z^n) = r^n S(n alpha) u, Sh = r^n cos, Ht = r^n sin", ok_pw, how_pw)
ex2 = E(H1) * sp.exp(S1)
report("Cor 4.59: e^z = e^Sh(z) [C(Ht z) + S(Ht z) u] is exp of the arrow on the page", *page_eq(ex2.page(), sp.exp(z1.page())))
report("Cor 4.59: Sh(e^z) = e^Sh(z) cos Ht(z), Ht(e^z) = e^Sh(z) sin Ht(z)",
       *allz(Sh(ex2) - sp.exp(S1) * sp.cos(H1), Ht(ex2) - sp.exp(S1) * sp.sin(H1)))

# Remark 4.60, Parts Longer than the Number
ok_b, top_ratio = True, 0.0
rb = random.Random(17)
for _ in range(20000):
    pq = rb.uniform(0.05, math.pi - 0.05)
    Lq = rb.uniform(0.2, 2.5)
    A0, B0 = rb.uniform(-3, 3), rb.uniform(-3, 3)
    wv = A0 + B0 * Lq * cmath.exp(1j * pq)
    m = abs(wv)
    if m < 1e-9:
        continue
    top_ratio = max(top_ratio, abs(A0) * math.sin(pq) / m, abs(B0) * Lq * math.sin(pq) / m)
    ok_b &= abs(wv.real) <= m * (1 + 1e-12) and abs(wv.imag) <= m * (1 + 1e-12)
report("Rem 4.60: |Sh|, |Ht| <= |z| and |Re|, |Com| <= |z|/sin(pb), on 20000 random numbers and planes",
       ok_b and top_ratio <= 1 + 1e-9, "numeric")
report("Rem 4.60: the bounds are reached, Re at alpha = pb - pi/2 and Com at alpha = pi/2",
       *allz(Cg(p - sp.pi / 2) - 1 / sp.sin(p), Sg(sp.pi / 2) - 1 / sp.sin(p)))
report("Rem 4.60: the Composer part of J is 1/sin(pb) long", *zero(J.b * L - 1 / sp.sin(p)))

# Appendix F.8: the rows not proved above
report("F.8 magnitude: |z|^2 = Re(z)^2 + 2 cos(pb) Re(z) b|G| + b^2|G|^2 = Sh(z)^2 + Ht(z)^2",
       *allz(n1 - (a**2 + 2 * sp.cos(p) * a * b * L + b**2 * L**2), n1 - (S1**2 + H1**2)))
ok_g, how_g = True, "symbolic"
for n in (-2, -1, 1, 2, 3, 4, 5, 6):
    Gn = G ** n
    o, h = allz(Gn.a + L**n * sp.sin((n - 1) * p) / sp.sin(p), Gn.b - L ** (n - 1) * sp.sin(n * p) / sp.sin(p),
                Sh(Gn) - L**n * sp.cos(n * p), Ht(Gn) - L**n * sp.sin(n * p))
    ok_g &= o
    how_g = how_g if h == "symbolic" else "numeric only"
report("F.8 powers of G, n = -2 ... 6: Re = -L^n sin((n-1)pb)/sin pb, Com = L^(n-1) (sin n pb/sin pb) G, Sh = L^n cos n pb, Ht = L^n sin n pb",
       ok_g, how_g)
cs_, sn_ = sp.cos(p), sp.sin(p)
units = [("1", GN(1, 0), (1, 0, 1, 0)),
         ("G", G, (0, 1, L * cs_, L * sn_)),
         ("conj(G)", Gb, (2 * L * cs_, -1, L * cs_, -L * sn_)),
         ("G^2", G * G, (-L**2, 2 * L * cs_, L**2 * sp.cos(2 * p), L**2 * sp.sin(2 * p))),
         ("1/G", G.inv(), (2 * cs_ / L, -1 / L**2, cs_ / L, -sn_ / L)),
         ("J", J, (-ct, 1 / (L * sn_), 0, 1))]
for name, zu, (ra, cb, shv, htv) in units:
    report(f"F.8 table of six numbers, {name}: real part, Composer coordinate, shadow, height",
           *allz(zu.a - ra, zu.b - cb, Sh(zu) - shv, Ht(zu) - htv))
report("F.8 table at pi/2: the six numbers are 1, i, -i, -1, -i and i",
       all(zero(sp.simplify(zu.page().subs(hp)) - v)[0] for (_, zu, _), v in zip(units, (1, I, -I, -1, -I, I))), "symbolic")

# At phibar = pi/2 the parts and the readings are one pair
report("pi/2: Re(z) = Sh(z) = a and the Composer coordinate is Ht(z) = b", *allz((z1.a - S1).subs(hp), (z1.b - H1).subs(hp)))
report("pi/2: Prop 4.54 (3) gives Re z = (z + conj z)/2 and Com(z) = (z - conj z)/2 = i Im z",
       c_eq((G * zc - Gb * z1) / (G - Gb), a, 0) and c_eq(G * (z1 - zc) / (G - Gb), 0, b)
       and c_eq((z1 - zc) * sp.Rational(1, 2), 0, b), "symbolic")
report("pi/2: Re(conj z) = Re(z) = a (Prop 4.56)", *zero((2 * S1 - a).subs(hp) - a))
report("pi/2: Thm 4.55 gives Re(zw) = ac - bd and Com(zw) = (ad + bc) i", c_eq(pr12, a * c - b * d, a * d + b * c), "symbolic")
report("pi/2: Cor 4.57 gives Re(iz) = -b and Im(iz) = a for both units", c_eq(Jz, -b, a) and c_eq(Gz, -b, a), "symbolic")
report("pi/2: Prop 4.58 gives conj(z) w = (ac + bd) + (ad - bc) i", c_eq(A_, a * c + b * d, a * d - b * c), "symbolic")
report("pi/2: the bound |z|/sin(pb) of Remark 4.60 is |z|", (1 / sp.sin(p)).subs(hp) == 1, "symbolic")

# laws on arrows multiplied directly on the page (complex numbers), independent of the class GN
ok_pg, cases, planes = True, 0, 0
rq = random.Random(23)


def parts_of(wv, pq, Lq):
    """Real part, Composer coordinate, shadow and height of the arrow wv (Proposition 4.44)."""
    return wv.real - wv.imag / math.tan(pq), wv.imag / (Lq * math.sin(pq)), wv.real, wv.imag


def cl(u1, u2):
    return abs(u1 - u2) <= 1e-9 * max(1.0, abs(u1), abs(u2))


for pq in [math.radians(79.86), math.pi / 3, math.pi / 2] + [rq.uniform(0.1, 3.0) for _ in range(40)]:
    if abs(pq - 3 * math.pi / 4) < 0.1:
        continue
    planes += 1
    Lq = abs(math.cos(pq) + math.sin(pq)) if rq.random() < 0.5 else rq.uniform(0.3, 2.0)
    Gq = Lq * cmath.exp(1j * pq)
    for _ in range(60):
        A1, B1, A2, B2 = (rq.uniform(-3, 3) for _ in range(4))
        w1, w2 = A1 + B1 * Gq, A2 + B2 * Gq
        if abs(w1) < 1e-3 or abs(w2) < 1e-3:
            continue
        R1, C1, S1v, H1v = parts_of(w1, pq, Lq)
        R2, C2, S2v, H2v = parts_of(w2, pq, Lq)
        Rp, Cp_, Sp_, Hp = parts_of(w1 * w2, pq, Lq)
        Rq, Cq, Sq, Hq = parts_of(w1 / w2, pq, Lq)
        Rc, Cc, Sc, Hc = parts_of(w1.conjugate(), pq, Lq)
        Ra, Ca, Sa, Ha = parts_of(w1.conjugate() * w2, pq, Lq)
        RJ, CJ, SJ, HJ = parts_of(1j * w1, pq, Lq)
        RG, CG, SG, HG = parts_of(Gq * w1, pq, Lq)
        m2 = abs(w2) ** 2
        ok_pg &= cl(Rp, R1 * R2 - Lq**2 * B1 * B2) and cl(Cp_, S1v * C2 + S2v * C1)
        ok_pg &= cl(Sp_, S1v * S2v - H1v * H2v) and cl(Hp, S1v * H2v + S2v * H1v)
        ok_pg &= cl(Rc, 2 * S1v - R1) and cl(Cc, -C1) and cl(Sc, S1v) and cl(Hc, -H1v)
        ok_pg &= cl(Rq, (R1 * (2 * S2v - R2) + Lq**2 * B1 * B2) / m2) and cl(Cq, (S2v * C1 - S1v * C2) / m2)
        ok_pg &= cl(Sq, (S1v * S2v + H1v * H2v) / m2) and cl(Hq, (S2v * H1v - S1v * H2v) / m2)
        ok_pg &= cl(Ca, A1 * B2 - B1 * A2) and cl(Ha, (A1 * B2 - B1 * A2) * Lq * math.sin(pq))
        ok_pg &= cl(Ra, (2 * S1v - R1) * R2 + Lq**2 * B1 * B2) and cl(Sa, S1v * S2v + H1v * H2v)
        ok_pg &= cl(RJ, -H1v - S1v / math.tan(pq)) and cl(CJ, S1v / (Lq * math.sin(pq)))
        ok_pg &= cl(RG, -Lq**2 * B1) and cl(CG, 2 * S1v - R1)
        ok_pg &= cl(R1, A1) and cl(C1, B1)
        cases += 1
report(f"the laws of Subsection 4.9.7 on arrows multiplied on the page: {cases} pairs in {planes} planes", ok_pg, "numeric")

# --------------------------------------------------------------------------
heading("16. Numbers of the worked example of Subsection 4.9.7 and of Figure 4.20, at 79.86 degrees")
f4 = lambda t: f"{t:.4f}"


def num(label, value, printed, fmt=f5):
    shown = fmt(value)
    report(f"{label}: {printed}" + ("" if shown == printed else f"  (computed {shown})"), shown == printed)


zA = N(1, 1)
ShA, HtA, n2A = zA.a + zA.b * Lv * cv, zA.b * Lv * sv, zA.n2()
num("|G|", Lv, "1.16043")
num("|G|^2", Lv**2, "1.34661")
num("cot(pb)", cv / sv, "0.17885")
num("1/(|G| sin pb)", 1 / (Lv * sv), "0.87542")
num("|G| cos(pb), one step across", Lv * cv, "0.20430")
num("Sh(1 + G)", ShA, "1.20430")
num("Ht(1 + G)", HtA, "1.14231")
num("Re(conj z)", zA.conj().a, "1.40860")
num("  2 x 1.20430 - 1", 2 * 1.20430 - 1, "1.40860")
num("  the shadow is halfway: (1 + 1.40860)/2", (1 + 1.40860) / 2, "1.20430")
sqA = zA * zA
num("Re(z^2)", sqA.a, "-0.34661")
num("  1 - 1.34661", 1 - 1.34661, "-0.34661")
num("Composer coordinate of z^2", sqA.b, "2.40860")
num("  2 x 1.20430", 2 * 1.20430, "2.40860")
num("Sh(z^2)", sqA.a + sqA.b * Lv * cv, "0.14547")
num("  1.20430^2 - 1.14231^2", 1.20430**2 - 1.14231**2, "0.14547")
num("Ht(z^2), four places", sqA.b * Lv * sv, "2.7514", f4)
num("  2 x 1.20430 x 1.14231, four places", 2 * 1.20430 * 1.14231, "2.7514", f4)
num("  the walk across: -0.34661 + 2.40860 x 0.20430", -0.34661 + 2.40860 * 0.20430, "0.14547")
num("  the walk up: 2.40860 x 1.14231, four places", 2.40860 * 1.14231, "2.7514", f4)
num("|z|^2", n2A, "2.75521")
ziA = N(zA.conj().a / n2A, zA.conj().b / n2A)
num("Re(1/z)", ziA.a, "0.51125")
num("  1.40860/2.75521", 1.40860 / 2.75521, "0.51125")
num("Composer coordinate of 1/z", ziA.b, "-0.36295")
num("  -1/2.75521", -1 / 2.75521, "-0.36295")
num("Sh(1/z)", ziA.a + ziA.b * Lv * cv, "0.43710")
num("  1.20430/2.75521", 1.20430 / 2.75521, "0.43710")
num("Ht(1/z)", ziA.b * Lv * sv, "-0.41460")
num("  -1.14231/2.75521", -1.14231 / 2.75521, "-0.41460")
GzA = N(0, 1) * zA
num("Re(G z)", GzA.a, "-1.34661")
num("Composer coordinate of G z, which is Re(conj z)", GzA.b, "1.40860")
JzA = N(-cv / sv, 1 / (Lv * sv)) * zA
num("Re(J z)", JzA.a, "-1.35770")
num("  -1.14231 - 1.20430 x 0.17885", -1.14231 - 1.20430 * 0.17885, "-1.35770")
num("Composer coordinate of J z", JzA.b, "1.05427")
num("  1.20430 x 0.87542", 1.20430 * 0.87542, "1.05427")
topA = zA.conj() * N(2, -1)
num("Re of conj(z) (2 - G), the top of (2 - G)/(1 + G)", topA.a, "1.47059")
num("Composer coordinate of conj(z) (2 - G)", topA.b, "-3.00000")
num("  1 x (-1) - 1 x 2", 1 * (-1) - 1 * 2, "-3.00000")
pA, pB = zA.page(), N(2, -1).page()
crossA = pA.real * pB.imag - pA.imag * pB.real
num("height of conj(z) (2 - G) = the signed area on the page", crossA, "-3.42693")
num("  -3 x 1.14231", -3 * 1.14231, "-3.42693")
num("  the height from the walk: -3 |G| sin(pb)", topA.b * Lv * sv, "-3.42693")
report("  2 - G lies clockwise from 1 + G (the area is negative)", crossA < 0)
num("J = -0.17885 + 0.87542 G, real part", -cv / sv, "-0.17885")
num("the Composer part of J, 1/sin(pb), four places", 1 / sv, "1.0159", f4)
num("  0.87542 x 1.16043, four places", 0.87542 * 1.16043, "1.0159", f4)
report("Remark 4.60: 1/sin(pb) - 1 is about 1.6 percent: " + f"{(1 / sv - 1) * 100:.1f}", f"{(1 / sv - 1) * 100:.1f}" == "1.6")
num("Figure 4.20: |G| sin(pb), one step up", Lv * sv, "1.14231")
num("Figure 4.20: 2 |G| cos(pb), the base of the triangle", 2 * Lv * cv, "0.40860")
num("Figure 4.20: the middle of the dashed leg, (0.20430 + 0.40860)/2", (0.20430 + 0.40860) / 2, "0.30645")
num("Figure 4.20: the leg from Re(conj z) to z is |G| long", abs(pA - zA.conj().a), "1.16043")
num("Figure 4.20: the tip of z^2, across", sqA.a + sqA.b * Lv * cv, "0.14547")
num("Figure 4.20: the tip of z^2, up", sqA.b * Lv * sv, "2.75136")

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
print("\n" + "=" * 60)
print(f"  {len(PASSED)} checks passed, {len(FAILED)} failed")
if FAILED:
    print("  Failed:")
    for f in FAILED:
        print("   -", f)
print("=" * 60)

import sys
sys.exit(1 if FAILED else 0)
