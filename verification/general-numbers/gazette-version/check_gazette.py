"""Check the identities and printed numbers of gn-gazette.tex against complex arithmetic."""
import cmath, math, random

random.seed(1)
ok = 0
def check(name, u, v, tol=1e-9):
    global ok
    if abs(u - v) > tol * (1 + abs(u) + abs(v)):
        raise SystemExit(f"FAIL {name}: {u} vs {v}")
    ok += 1

def frame(phi, L):
    G = L * cmath.exp(1j * phi)
    def coords(z):  # z = a + b G
        b = z.imag / G.imag
        return z.real - b * G.real, b
    return G, coords

for _ in range(300):
    phi = random.uniform(0.05, math.pi - 0.05)
    L = random.uniform(0.2, 3)
    G, co = frame(phi, L)
    u = cmath.exp(1j * phi)
    a, b, c, d = (random.uniform(-3, 3) for _ in range(4))
    z, w = a + b * G, c + d * G
    sh = lambda x: x.real
    ht = lambda x: x.imag
    check("shadow", sh(z), a + b * L * math.cos(phi)); check("height", ht(z), b * L * math.sin(phi))
    check("square", G * G, 2 * L * math.cos(phi) * G - L**2)
    s, t = co(z * w)
    check("prod re", s, a * c - L**2 * b * d); check("prod com", t, a * d + b * c + 2 * L * math.cos(phi) * b * d)
    check("com product", t * G, sh(z) * d * G + sh(w) * b * G)
    check("conj", z.conjugate(), (a + 2 * b * L * math.cos(phi)) - b * G)
    check("skew", (a + b * G) * (a - b * G), a*a + b*b*L*L - 2*b*b*L*math.cos(phi)*G)
    check("z+zbar", z + z.conjugate(), 2 * sh(z))
    check("norm", abs(z)**2, a*a + 2*a*b*L*math.cos(phi) + b*b*L*L)
    J = (u - math.cos(phi)) / math.sin(phi)
    check("J", J, 1j); check("J alt", J, -1/math.tan(phi) + G / (L * math.sin(phi)))
    check("bound", 0, max(0, abs(a) - abs(z) / math.sin(phi)) )
    th, al, be = (random.uniform(-4, 4) for _ in range(3))
    C = lambda x: math.sin(phi - x) / math.sin(phi); S = lambda x: math.sin(x) / math.sin(phi)
    check("polar", cmath.exp(1j * th), C(th) + S(th) * u)
    check("pyth", C(th)**2 + 2*math.cos(phi)*C(th)*S(th) + S(th)**2, 1)
    check("addC", C(al+be), C(al)*C(be) - S(al)*S(be))
    check("addS", S(al+be), C(al)*S(be) + S(al)*C(be) + 2*math.cos(phi)*S(al)*S(be))
    for n in range(-4, 7):
        check("powers", G**n, -L**n * math.sin((n-1)*phi)/math.sin(phi) + L**(n-1)*math.sin(n*phi)/math.sin(phi)*G)
    check("step", G**4, 2*L*math.cos(phi)*G**3 - L**2*G**2)
    # Cauchy-Riemann: f = z^3 written as U + V G as a function of (a, b)
    f = lambda a_, b_: co((a_ + b_ * G)**3)
    h = 1e-6
    Ua = (f(a+h,b)[0]-f(a-h,b)[0])/(2*h); Va = (f(a+h,b)[1]-f(a-h,b)[1])/(2*h)
    Ub = (f(a,b+h)[0]-f(a,b-h)[0])/(2*h); Vb = (f(a,b+h)[1]-f(a,b-h)[1])/(2*h)
    check("CR1", Ub, -L**2 * Va, 1e-5); check("CR2", Vb, Ua + 2*L*math.cos(phi)*Va, 1e-5)

# elliptic pair converse
L, phi = 1.3, 1.1
p, q = -L**2, 2*L*math.cos(phi)
check("pair len", math.sqrt(-p), L); check("pair ang", math.acos(q/(2*math.sqrt(-p))), phi)
# Two angles of one length
A = 1.16043
p1 = 0.5 * math.asin(A*A - 1)
print("angles", round(math.degrees(p1), 2), round(90 - math.degrees(p1), 2))
# pi/4 and pi/5
G = math.sqrt(2) * cmath.exp(1j * math.pi / 4)
check("G^8", G**8, 16)
tau = (1 + 5**0.5) / 2
check("tau", 2*math.cos(math.pi/5), tau)
print("pi/5 coeffs", [round(math.sin(n*math.pi/5)/math.sin(math.pi/5), 6) for n in range(1, 6)], round(tau, 6))
# Example table
phi = math.radians(79.86); L = math.cos(phi) + math.sin(phi)
G, co = frame(phi, L); u = cmath.exp(1j*phi); z = 1 + G
J = 1j
r = lambda x: round(x, 5)
print("cos sin L", r(math.cos(phi)), r(math.sin(phi)), r(L))
for name, val in [("G^2", G*G), ("J", J), ("zbar", z.conjugate()), ("z^2", z*z), ("1/z", 1/z), ("(1+G)(1-G)", (1+G)*(1-G))]:
    print(name, [r(x) for x in co(val)])
print("Sh Ht", r(z.real), r(z.imag), "|z|", r(abs(z)), "angle", round(math.degrees(cmath.phase(z)), 2))
print("Com(J) length", r(abs(co(J)[1]) * L), "1/sin", r(1/math.sin(phi)))
print("Com(z^2) coeff", r(co(z*z)[1]), "2Sh", r(2*z.real))
print("pi/3 length", r(math.cos(math.pi/3)+math.sin(math.pi/3)))
print("checks passed:", ok)
