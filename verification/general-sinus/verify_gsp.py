"""
verify_gsp.py -- checks every numbered statement of
"The General Sinus: Trigonometric Functions of an Oblique Projection".

Identities are tested at 50 significant digits at many random points of their domains.
Derivatives are computed symbolically with SymPy and compared at 50 digits.
Requires: Python 3, SymPy (which installs mpmath). Run: python verify_gsp.py
"""
import math
import random
import sympy as sp
from mpmath import mp, mpf, sin, cos, tan, cot, pi, nstr

mp.dps = 50
random.seed(2026)
TOL = mpf(10) ** -40

S = lambda a, p: sin(a) / sin(a + p)
C = lambda a, p: sin(p) / sin(a + p)
T = lambda a, p: sin(a) / sin(p)
CSC = lambda a, p: 1 / S(a, p)
SEC = lambda a, p: 1 / C(a, p)
COT = lambda a, p: 1 / T(a, p)

passed, failed = 0, []


def rnd():
    return mpf(random.uniform(-3.1, 3.1))


def ok_point(*vals):
    """Reject points too close to a pole (a sine of any argument below 1e-3)."""
    return all(abs(sin(v)) > 1e-3 for v in vals)


def check(label, fn, args=4, guard=lambda *v: v, n=300, tol=TOL):
    global passed
    worst = mpf(0)
    tested = 0
    while tested < n:
        v = [rnd() for _ in range(args)]
        if not ok_point(*guard(*v)):
            continue
        lhs, rhs = fn(*v)
        worst = max(worst, abs(lhs - rhs) / (1 + abs(lhs) + abs(rhs)))
        tested += 1
    good = worst < tol
    passed += good
    if not good:
        failed.append(label)
    print(f"  {'ok  ' if good else 'FAIL'} {label:<40} worst rel. error {nstr(worst, 3)}")


print("Section 2: definitions")
# Theorem 2.2 by the law of sines on an explicit triangle
def triangle(a, p):
    a, p = abs(a) % 1.5 + mpf("0.05"), abs(p) % 1.5 + mpf("0.05")
    Z = (cos(a), sin(a))                         # [OZ] = 1
    # foot: line through Z making angle p with ZO meets the x-axis
    direction = a + pi + p                        # rotate the ray Z->O by -p ... choose the side that hits the axis
    d1 = (cos(pi + a - p), sin(pi + a - p))
    d2 = (cos(pi + a + p), sin(pi + a + p))
    for d in (d1, d2):
        if d[1] < 0:
            t = -Z[1] / d[1]
            foot = (Z[0] + t * d[0], mpf(0))
            if foot[0] > 0:
                break
    zz = math.hypot(float(Z[0] - foot[0]), float(Z[1]))
    oz_bar = float(foot[0])
    return (mpf(zz), S(a, p)), (mpf(oz_bar), C(a, p))
def thm_analytic(a, p, *_):
    (x1, y1), (x2, y2) = triangle(a, p)
    return x1 + x2, y1 + y2
check("thm:analytic (law of sines, geometry)", thm_analytic, tol=mpf(10) ** -12, guard=lambda a, p, b, t: (a + p, a, p))
check("cor:classical", lambda a, *_: (S(a, pi / 2 - a) + C(a, pi / 2 - a) + T(a, pi / 2 - a),
                                       sin(a) + cos(a) + tan(a)), guard=lambda a, p, b, t: (a, a + pi / 2))

print("Section 3: coordinates")
def coords(a, p, *_):
    pb = a + p
    vx = C(a, p) * 1 + S(a, p) * cos(pb)
    vy = S(a, p) * sin(pb)
    return vx + vy, cos(a) + sin(a)
check("prop:coordinates", coords, guard=lambda a, p, b, t: (a + p,))
check("prop:pythagoras", lambda a, p, *_: (C(a, p) ** 2 + S(a, p) ** 2 + 2 * C(a, p) * S(a, p) * cos(a + p), 1),
      guard=lambda a, p, b, t: (a + p,))

print("Section 4: symmetries and special values")
g = lambda a, p, b, t: (a, p, a + p, a - p)
check("prop:exchange (sin=cos swapped)", lambda a, p, *_: (S(a, p), C(p, a)), guard=g)
check("prop:exchange (tan=cot swapped)", lambda a, p, *_: (T(a, p), COT(p, a)), guard=g)
check("prop:exchange (csc=sec swapped)", lambda a, p, *_: (CSC(a, p), SEC(p, a)), guard=g)
check("prop:exchange (csc=sin(a+p,-p))", lambda a, p, *_: (CSC(a, p), S(a + p, -p)), guard=g)
check("prop:exchange (sec=tan(a+p,p))", lambda a, p, *_: (SEC(a, p), T(a + p, p)), guard=g)
check("prop:exchange (cot=cos(a-p,p))", lambda a, p, *_: (COT(a, p), C(a - p, p)), guard=g)
check("prop:symmetry (sin(-a,-p))", lambda a, p, *_: (S(-a, -p), S(a, p)), guard=g)
check("prop:symmetry (cos(-a,-p))", lambda a, p, *_: (C(-a, -p), C(a, p)), guard=g)
check("prop:symmetry (sin(-a,p)=sin(a,-p))", lambda a, p, *_: (S(-a, p), S(a, -p)), guard=g)
for P_ in range(-2, 3):
    for Q_ in range(-2, 3):
        pass
def shifts(a, p, *_):
    tot_l = tot_r = mpf(0)
    for P_ in range(-2, 3):
        for Q_ in range(-2, 3):
            tot_l += S(a + P_ * pi, p + Q_ * pi) + C(a + P_ * pi, p + Q_ * pi) + T(a + P_ * pi, p + Q_ * pi)
            tot_r += (-1) ** Q_ * S(a, p) + (-1) ** P_ * C(a, p) + (-1) ** (P_ + Q_) * T(a, p)
    return tot_l, tot_r
check("prop:symmetry (pi shifts)", shifts, guard=g)
check("prop:special (phi=pi/2)", lambda a, *_: (S(a, pi / 2) + C(a, pi / 2) + T(a, pi / 2),
                                                tan(a) + 1 / cos(a) + sin(a)), guard=lambda a, *_: (a, a + pi / 2))
check("prop:special (theta,theta)", lambda a, *_: (S(a, a) + C(a, a) + T(a, a), 1 / cos(a) + 1),
      guard=lambda a, *_: (a, 2 * a))
check("special value pi/3 all six = 1", lambda *_: (S(pi/3, pi/3) + C(pi/3, pi/3) + T(pi/3, pi/3)
                                                   + CSC(pi/3, pi/3) + SEC(pi/3, pi/3) + COT(pi/3, pi/3), 6), n=1)
tau = (1 + mp.sqrt(5)) / 2
check("special value pi/5 = 1/tau", lambda *_: (S(pi/5, pi/5), 1 / tau), n=1)

print("Section 5: addition formulas")
check("lem:ptolemy", lambda x, y, z, _: (sin(x + y) * sin(y + z), sin(x) * sin(z) + sin(y) * sin(x + y + z)))
ga = lambda a, p, b, t: (a, b, p, a + p, b + p, a + b + p, a + b)
check("eq:add-sin", lambda a, p, b, t: (S(a + b, p), S(a, p) + C(a, p) * S(b, p + a)), guard=ga)
check("eq:add-cos", lambda a, p, b, t: (C(a + b, p), C(a, p) * C(b, p) + S(a, p) * S(b, p) * C(a + b, p)), guard=ga)
check("eq:add-tan", lambda a, p, b, t: (T(a + b, p), T(a, p) * cos(b) + T(b, p) * cos(a)), guard=ga)
check("eq:add-closed (cos)", lambda a, p, b, t: (C(a + b, p), C(a, p) * C(b, p) / (1 - S(a, p) * S(b, p))), guard=ga)
check("eq:add-closed (sec)", lambda a, p, b, t: (SEC(a + b, p), SEC(a, p) * SEC(b, p) - T(a, p) * T(b, p)), guard=ga)
gp = lambda a, p, b, t: (a, p, t, a + p, p + t, a + p + t)
check("eq:addp-sin", lambda a, p, b, t: (S(a, p + t), S(a, p) * S(a + p, t)), guard=gp)
check("eq:addp-cos", lambda a, p, b, t: (C(a, p + t), C(a, p) + S(a, p) * C(a + p, t)), guard=gp)
gd = lambda a, p, b, t: (a, p, a + p, 2 * a + p, a / 2 + p, a / 2)
check("cor:double (sin)", lambda a, p, *_: (S(2 * a, p), S(a, p) * (1 + C(2 * a, p))), guard=gd)
check("cor:double (cos)", lambda a, p, *_: (C(2 * a, p), C(a, p) ** 2 / (1 - S(a, p) ** 2)), guard=gd)
check("cor:double (sec)", lambda a, p, *_: (SEC(2 * a, p), SEC(a, p) ** 2 - T(a, p) ** 2), guard=gd)
check("cor:double (tan)", lambda a, p, *_: (T(2 * a, p), 2 * T(a, p) * cos(a)), guard=gd)
check("cor:double (half angle)", lambda a, p, *_: (C(a, p), C(a / 2, p) ** 2 / (1 - S(a / 2, p) ** 2)), guard=gd)

print("Section 6: derivatives (SymPy, compared at 50 digits)")
x, y, t = sp.symbols("x y t", real=True)
s_, c_, t_ = sp.sin(x) / sp.sin(x + y), sp.sin(y) / sp.sin(x + y), sp.sin(x) / sp.sin(y)
def num(expr, a, p):
    return sp.N(expr.subs({x: sp.Float(str(a), 50), y: sp.Float(str(p), 50)}), 50)
claims = {
    "prop:partials d_x sin": (sp.diff(s_, x), c_ / sp.sin(x + y)),
    "prop:partials d_y sin": (sp.diff(s_, y), -s_ * sp.cot(x + y)),
    "prop:partials d_x cos": (sp.diff(c_, x), -c_ * sp.cot(x + y)),
    "prop:partials d_y cos": (sp.diff(c_, y), s_ / sp.sin(x + y)),
    "prop:partials d_x tan": (sp.diff(t_, x), sp.cos(x) / sp.sin(y)),
    "prop:partials d_y tan": (sp.diff(t_, y), -t_ * sp.cot(y)),
    "prop:partials d_xx sin": (sp.diff(s_, x, 2), -2 * sp.cot(x + y) * sp.diff(s_, x)),
    "prop:partials d_yy sin": (sp.diff(s_, y, 2), s_ * (1 + 2 * sp.cot(x + y) ** 2)),
    "prop:partials d_xx cos": (sp.diff(c_, x, 2), c_ * (1 + 2 * sp.cot(x + y) ** 2)),
    "prop:partials d_yy cos": (sp.diff(c_, y, 2), -2 * sp.cot(x + y) * sp.diff(c_, y)),
    "prop:partials d_xx tan": (sp.diff(t_, x, 2), -t_),
    "prop:partials d_yy tan": (sp.diff(t_, y, 2), t_ * (1 + 2 * sp.cot(y) ** 2)),
}
path = lambda e: e.subs({x: x + t, y: y - t}, simultaneous=True)
gen = lambda f, a, b: f.subs({x: a, y: b}, simultaneous=True)
for n in range(1, 7):
    Dn_s = sp.diff(path(s_), t, n).subs(t, 0)
    Dn_c = sp.diff(path(c_), t, n).subs(t, 0)
    Dn_t = sp.diff(path(t_), t, n).subs(t, 0)
    claims[f"thm:shift D^{n} sin"] = (Dn_s, gen(s_, x + n * sp.pi / 2, y - n * sp.pi / 2))
    claims[f"thm:shift D^{n} cos"] = (Dn_c, gen(c_, x + n * sp.pi / 2, y - n * sp.pi / 2))
    claims[f"D^{n} tan formula"] = (Dn_t, sp.sin(x + y) * sp.diff(sp.csc(y) ** 2, y, n - 1) * (-1) ** (n - 1))
    # primitives: D^n of the shifted function returns the function
    F_s = gen(s_, x - n * sp.pi / 2, y + n * sp.pi / 2)
    claims[f"cor:primitives n={n} (sin)"] = (sp.diff(path(F_s), t, n).subs(t, 0), s_)
for label, (lhs, rhs) in claims.items():
    worst = mpf(0)
    k = 0
    while k < 25:
        a, p = random.uniform(-1.4, 1.4), random.uniform(0.2, 2.9)
        if min(abs(math.sin(a + p)), abs(math.sin(p)), abs(math.sin(a))) < 0.05:
            continue
        l, r = num(lhs, a, p), num(rhs, a, p)
        worst = max(worst, mpf(str(abs(l - r) / (1 + abs(l) + abs(r)))))
        k += 1
    good = worst < TOL
    passed += good
    if not good:
        failed.append(label)
    print(f"  {'ok  ' if good else 'FAIL'} {label:<40} worst rel. error {nstr(worst, 3)}")

print("Section 7: product identities and the General Sinus Transform")
def products(a, p, *_):
    lhs = rhs = mpf(0)
    for N in range(1, 41):
        ps = pc = pt = mpf(1)
        for k in range(N):
            ps *= S((a + pi * k) / N, p / N)
            pc *= C(a / N, (p + pi * k) / N)
            pt *= T((a + pi * k) / N, (p + pi * k) / N)
        lhs += ps + pc + pt
        rhs += S(a, p) + C(a, p) + T(a, p)
    return lhs, rhs
check("thm:products, N = 1..40", products, n=40, tol=mpf(10) ** -35, guard=lambda a, p, b, t: (a + p, p))

def P(N, k, z):
    u, v = math.sin(math.pi * k / N), math.cos(math.pi * k / N)
    return u * z * z - 2 * N * v * z - 2 * N * N * u
def S_N(N, a, p):
    return math.prod(P(N, k, a) / P(N, k, a + p) for k in range(N))
pts = [(0.5, 0.7), (1.0, 1.0), (2.0, 0.9), (0.3, 2.7)]
printed = {64: "1.6e-02", 1024: "6.3e-05", 16384: "2.5e-07"}
for N, txt in printed.items():
    err = max(abs(S_N(N, a, p) - math.sin(a) / math.sin(a + p)) for a, p in pts)
    good = f"{err:.1e}" == txt
    passed += good
    if not good:
        failed.append(f"rem:accuracy N={N}")
    print(f"  {'ok  ' if good else 'FAIL'} rem:accuracy N={N:<5} printed {txt}, computed {err:.1e}, N^2*err = {N*N*err:.1f}")
for N in (100, 3000):
    err = max(abs(S_N(N, a, p) - math.sin(a) / math.sin(a + p)) for a, p in pts)
    good = f"{N*N*err:.1f}" == "66.1"
    passed += good
    if not good:
        failed.append(f"rem:accuracy N^2 err N={N}")
    print(f"  {'ok  ' if good else 'FAIL'} rem:accuracy N^2*err = 66.1 at N={N:<5} computed {N*N*err:.1f}")

print("Section 6: fixed frame (Proposition 6.2)")
from mpmath import mpc, exp as mexp, atan, sqrt as msqrt, chebyu
def fixed_frame(pb, a, b, *_):
    pb = abs(pb) % 3 + mpf("0.05")
    Cf = lambda t: sin(pb - t) / sin(pb)
    Sf = lambda t: sin(t) / sin(pb)
    u = mexp(mpc(0, 1) * pb)
    lhs = Cf(a + b) + Sf(a + b)
    rhs = (Cf(a) * Cf(b) - Sf(a) * Sf(b)) + (Cf(a) * Sf(b) + Sf(a) * Cf(b) + 2 * cos(pb) * Sf(a) * Sf(b))
    tot = abs(mexp(mpc(0, 1) * a) - (Cf(a) + Sf(a) * u))          # oblique Euler formula
    tot += abs(u * u - (2 * cos(pb) * u - 1))                       # u^2 = 2 cos(pb) u - 1
    for n in range(-3, 6):
        tot += abs((Cf(a) + Sf(a) * u) ** n - (Cf(n * a) + Sf(n * a) * u))   # De Moivre
    for n in range(2, 8):
        tot += abs(u ** n - (-chebyu(n - 2, cos(pb)) + chebyu(n - 1, cos(pb)) * u))  # Chebyshev powers
    tot += abs(Cf(a)**2 + 2*cos(pb)*Cf(a)*Sf(a) + Sf(a)**2 - 1)        # Pythagoras in C,S notation
    return lhs + tot, rhs
check("prop:fixed-frame, eq:oblique-euler, De Moivre", fixed_frame, guard=lambda a, b, c, d: (1,))

print("Section 7: closed form and rate of the General Sinus Transform")
rho = lambda e: msqrt(1 + e**4 / 4)
dlt = lambda e: atan(e / (1 - e**2 / 2))
def Pm(N, k, z):
    return sin(pi * k / N) * z * z - 2 * N * cos(pi * k / N) * z - 2 * N * N * sin(pi * k / N)
def SNm(N, a, b): return mp.fprod(Pm(N, k, a) / Pm(N, k, a + b) for k in range(N))
def CNm(N, a, b): return mp.fprod(Pm(N, k, b) / Pm(N, k, a + b) for k in range(N))
def TNm(N, a, b): return mp.fprod(Pm(N, k, a) / Pm(N, k, b) for k in range(N))
def closed(a, b, *_):
    a, b = a / 2, b / 2 + mpf("0.2")
    lhs = rhs = mpf(0)
    for N in (4, 7, 16, 33, 100):
        e = lambda z: z / N
        lhs += SNm(N, a, b) + CNm(N, a, b) + TNm(N, a, b)
        rhs += ((rho(e(a)) / rho(e(a + b))) ** N * sin(N * dlt(e(a))) / sin(N * dlt(e(a + b)))
                + (rho(e(b)) / rho(e(a + b))) ** N * sin(N * dlt(e(b))) / sin(N * dlt(e(a + b)))
                + (rho(e(a)) / rho(e(b))) ** N * sin(N * dlt(e(a))) / sin(N * dlt(e(b))))
    return lhs, rhs
check("thm:GST-closed (N = 4, 7, 16, 33, 100)", closed, n=40, tol=mpf(10) ** -30,
      guard=lambda a, b, c, d: (a / 2 + b / 2 + mpf("0.2"), b / 2 + mpf("0.2")))
cc = lambda a, b: (a**3 * cos(a) * sin(a + b) - (a + b)**3 * sin(a) * cos(a + b)) / (6 * sin(a + b)**2)
def rate(a, b, *_):
    a, b = a / 2, b / 2 + mpf("0.2")
    N = 10**6                       # via the closed form (exact), N^2 err -> c with O(1/N)
    e = lambda z: z / N
    SN = (rho(e(a)) / rho(e(a + b))) ** N * sin(N * dlt(e(a))) / sin(N * dlt(e(a + b)))
    return N * N * (SN - S(a, b)), cc(a, b)
check("cor:GST-rate coefficient c(x,y)", rate, n=60, tol=mpf(10) ** -4,
      guard=lambda a, b, c, d: (a / 2 + b / 2 + mpf("0.2"), a / 2))
x0, y0 = mpf("0.3"), mpf("2.7")
quoted = [("c(0.3,2.7) = 66.14", cc(x0, y0), "66.14"),
          ("N^2 err at N=64 = 66.28", 64**2 * (SNm(64, x0, y0) - S(x0, y0)), "66.28"),
          ("N^2 err at N=1000 = 66.12", 1000**2 * (SNm(1000, x0, y0) - S(x0, y0)), "66.12")]
for label, val, txt in quoted:
    good = mp.nstr(val, 4) == txt
    passed += good
    if not good:
        failed.append(label)
    print(f"  {'ok  ' if good else 'FAIL'} rem:accuracy {label:<28} computed {mp.nstr(val, 6)}")

print()
print(f"{passed} checks passed, {len(failed)} failed" + (": " + ", ".join(failed) if failed else "."))

import sys
sys.exit(1 if failed else 0)
