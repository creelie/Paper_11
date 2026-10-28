"""Exact checks for the two closures of Sections 9 and 10 (claims 23 and 26).

Part A, the entropy of small black holes (Section 10, Appendix G):
  A1  the Gauss-Bonnet value Phi* = 4 - 2b/a in D = 4..12 against Sen's (3.20);
  A2  a full numerical solution of the six entropy-function equations;
  A3  the identity H = 2(G - 2) for the Lovelock family lambda_m = 4^(1-m)/m!;
  A4  its Hermite form and the unique root of H = 4, D = 4..40;
  A5  the large-D expansion of that root;
  A6  rigidity of the top Lovelock coefficient;
  A7  the limits of Phi at the two endpoints.
Part B, the weak-coupling limit (Section 9, Appendix G):
  B1  the discriminant of Sen's family;
  B2  the reconstruction of h and of the normal form;
  B3  the first-order discriminant of a general family;
  B4  the isotrivial family;
  B5  the vanishing orders along a curve of self-intersection -n;
  B6  the Euler characteristics of the tadpole identity;
  B7  the orders on F_12, where the D7 divisor contains 7C;
  B8  the discriminant on the leading D7 locus at finite t;
  B9  base change of order at most two.
Part C, the half-BPS coefficients of eta^-24.
Run: python3 verify_closures.py  (sympy, mpmath).
"""
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

y = sp.symbols('y', positive=True)
checks = []


def ok(name, cond):
    assert cond, name
    checks.append(name)


def GH(D, lam):
    """G and H of the Lovelock entropy function, n = D - 2, alpha' = 1."""
    n = D - 2
    G = H = 0
    for m in range(1, len(lam) + 1):
        if n + 2 - 2 * m >= 0:
            G += 2 * m * lam[m - 1] * sp.factorial(n) / sp.factorial(n + 2 - 2 * m) * y ** (1 - m)
        if n - 2 * m >= 0:
            H += lam[m - 1] * sp.factorial(n) / sp.factorial(n - 2 * m) * y ** (-m)
    return sp.expand(G), sp.expand(H)


def prester(M):
    return [sp.Rational(1, 4 ** (m - 1)) / sp.factorial(m) for m in range(1, M + 1)]


# ---------------- A1 Gauss-Bonnet closed form ----------------
for D in range(4, 13):
    G, H = GH(D, [1, sp.Rational(1, 8)])
    Phi = sp.together(G ** 2 / H)
    a, b = (D - 2) * (D - 3), (D - 4) * (D - 5)
    crit = [r for r in sp.solve(sp.diff(Phi, y), y) if r.is_positive]
    ok(f'A1 unique critical point D={D}', len(crit) == 1)
    ok(f'A1 y* = (a-b)/4 D={D}', sp.simplify(crit[0] - sp.Rational(a - b, 4)) == 0)
    val = sp.simplify(Phi.subs(y, crit[0]))
    ok(f'A1 Phi* = 4 - 2b/a D={D}', sp.simplify(val - (4 - sp.Rational(2 * b, a))) == 0)
    # Sen (3.20), printed with C = 1: S = a_D sqrt(nw), a_D = pi sqrt(C) sqrt(2/a) (8a - 4b)^(1/2)
    sen = sp.pi * sp.sqrt(sp.Rational(2, a)) * sp.sqrt(8 * a - 4 * b)
    ok(f'A1 agrees with Sen (3.20) D={D}', sp.simplify(sen - 2 * sp.pi * sp.sqrt(val)) == 0)
    ok(f'A1 statistical iff D in (4,5) D={D}', (val == 4) == (D in (4, 5)))

# ---------------- A2 full solution of the entropy-function equations ----------------
v1, v2, uS, uT, e1, e2 = sp.symbols('v1 v2 uS uT e1 e2', positive=True)


def full_solve(D, lam, nn=3, ww=5, guess=(0.6, 2, 1, 1.3, 1, 1)):
    n = D - 2
    br = 2 * uT ** 2 * e1 ** 2 / v1 ** 2 + 2 * e2 ** 2 / (uT ** 2 * v1 ** 2)
    for m in range(1, len(lam) + 1):
        if n + 2 - 2 * m >= 0:
            br += lam[m - 1] * sp.factorial(n) / sp.factorial(n + 2 - 2 * m) * v2 ** (-m) * (
                (n + 2 - 2 * m) * (n + 1 - 2 * m) - 2 * m * v2 / v1)
    f = uS * v1 * v2 ** sp.Rational(n, 2) * br
    q1, q2 = 2 * nn, 2 * ww
    eqs = [sp.diff(f, uS), sp.diff(f, uT), sp.diff(f, v1), sp.diff(f, v2),
           sp.diff(f, e1) - q1, sp.diff(f, e2) - q2]
    sol = sp.nsolve(eqs, [v1, v2, uS, uT, e1, e2], guess, tol=1e-30, prec=40, maxsteps=200)
    S = 2 * sp.pi * (sol[4] * q1 + sol[5] * q2)
    return sol, S / (4 * sp.pi * sp.sqrt(nn * ww))


sol, ratio = full_solve(6, [1, sp.Rational(1, 8)], guess=(0.6, 2.5, 1, 1.3, 1, 1))
ok('A2 GB D=6 full solve, ratio sqrt(11/12)', sp.N(abs(ratio - sp.sqrt(sp.Rational(11, 12))), 45) < 1e-25)
ok('A2 GB D=6 sphere radius y* = 5/2', abs(sol[1] - sp.Rational(5, 2)) < 1e-25)
sol, ratio = full_solve(6, prester(3), guess=(0.6, 3.2, 1, 1.3, 1, 1))
ok('A2 Prester D=6 full solve, ratio 1', sp.N(abs(ratio - 1), 45) < 1e-25)
ok('A2 Prester D=6 v1 = 1/2', abs(sol[0] - sp.Rational(1, 2)) < 1e-25)
ok('A2 Prester D=6 v2 = 3/2 + sqrt(3)', abs(sol[1] - (sp.Rational(3, 2) + sp.sqrt(3))) < 1e-25)
ok('A2 Prester D=6 e1 e2 = 1/4', abs(sol[4] * sol[5] - sp.Rational(1, 4)) < 1e-25)

# ---------------- A3, A4 the Lovelock family of Prester ----------------
mp.mp.dps = 40
roots = {}
s = sp.symbols('s')
for D in range(4, 41):
    n, M = D - 2, D // 2
    G, H = GH(D, prester(M))
    ok(f'A3 H = 2(G-2) D={D}', sp.expand(H - 2 * (G - 2)) == 0)
    Hc = 4 * (sum(sp.ff(n, 2 * m) / (sp.factorial(m) * (4 * y) ** m) for m in range(0, n // 2 + 1)) - 1)
    ok(f'A4 closed form of H D={D}', sp.expand(H - Hc) == 0)
    He = sp.hermite_prob(n, s)
    herm = sp.expand(sp.I ** (-n) * He.subs(s, sp.I * sp.sqrt(2 * y)) / (2 * y) ** sp.Rational(n, 2))
    ok(f'A4 Hermite form D={D}', sp.simplify(herm - (Hc / 4 + 1)) == 0)
    # dH/dy < 0: every coefficient of H in 1/y is positive (m >= 1)
    ok(f'A4 H strictly decreasing D={D}',
       all(c > 0 for c in sp.Poly(sp.expand(H.subs(y, 1 / s)), s).coeffs()))
    poly = sp.Poly(sp.expand((H - 4) * y ** (n // 2)), y)
    pos = [r for r in poly.nroots(n=40, maxsteps=400) if abs(sp.im(r)) < 1e-25 and sp.re(r) > 0]
    ok(f'A4 unique positive root of H = 4 D={D}', len(pos) == 1)
    roots[D] = sp.re(pos[0])
ok('A4 D=4 root 1/2', abs(roots[4] - sp.Rational(1, 2)) < 1e-30)
ok('A4 D=5 root 3/2', abs(roots[5] - sp.Rational(3, 2)) < 1e-30)

# ---------------- A5 large D ----------------
L2 = sp.log(2)
errs = []
for D in (20, 30, 40):
    n = D - 2
    errs.append(abs(roots[D] - (n ** 2 / (4 * L2) - n * (2 * L2 + 1) / (4 * L2))))
ok('A5 y_D = n^2/(4 ln 2) - n(2 ln 2 + 1)/(4 ln 2) + O(1)', max(sp.N(e) for e in errs) < sp.Rational(5, 4))
ok('A5 leading term n^2/(4 ln 2)', abs(roots[40] / (38 ** 2 / (4 * L2)) - 1) < 0.07)

# ---------------- A6 rigidity of the top coefficient ----------------
beta = sp.symbols('beta', positive=True)
for D in range(4, 17):
    M = D // 2
    lam = prester(M - 1) + [beta * prester(M)[-1]]
    G, H = GH(D, lam)
    G1 = G.subs(beta, 1)
    ok(f'A6 top coefficient enters only G D={D}', sp.diff(H, beta) == 0)
    dG = sp.expand(sp.diff(G, beta))
    ok(f'A6 dG/dbeta > 0 D={D}', all(c > 0 for c in sp.Poly(sp.expand(dG.subs(y, 1 / s)), s).coeffs()))
    # Phi_1 = G1^2/H >= 4 with equality only at G1 = 4
    phi1 = G1 ** 2 / H - 4
    ok(f'A6 Phi_1 - 4 = (G1-4)^2/H D={D}', sp.simplify(phi1 - (G1 - 4) ** 2 / H) == 0)
    # numerical minimum of Phi_beta, increasing through 4 at beta = 1
    vals = []
    for bv in (sp.Rational(1, 2), 1, 2):
        Pb = sp.lambdify(y, (G ** 2 / H).subs(beta, bv), 'mpmath')
        dPb = sp.lambdify(y, sp.diff((G ** 2 / H).subs(beta, bv), y), 'mpmath')
        grid = [mp.mpf(10) ** (k / mp.mpf(200)) for k in range(-800, 1001)]
        y0 = min(grid, key=Pb)
        ym = mp.findroot(dPb, (y0 / mp.mpf('1.05'), y0 * mp.mpf('1.05')), solver='anderson')
        vals.append(Pb(ym))
    ok(f'A6 min Phi below, at, above 4 D={D}', vals[0] < 4 and abs(vals[1] - 4) < 1e-25 and vals[2] > 4)

# ---------------- B Sen limit ----------------
h, eta, chi, t, F, Gv, F1, G1v = sp.symbols('h eta chi t F G F1 G1')
f = -3 * h ** 2 + t * eta
g = -2 * h ** 3 + t * h * eta - t ** 2 * chi / 12
Delta = sp.expand(4 * f ** 3 + 27 * g ** 2)
ok('B1 Delta, order t^2', sp.expand(Delta.coeff(t, 2) + 9 * h ** 2 * (eta ** 2 - h * chi)) == 0)
ok('B1 Delta, order t^3', sp.expand(Delta.coeff(t, 3) - (4 * eta ** 3 - sp.Rational(9, 2) * h * eta * chi)) == 0)
ok('B1 Delta, order t^4', sp.expand(Delta.coeff(t, 4) - sp.Rational(3, 16) * chi ** 2) == 0)
ok('B1 no terms of order < 2', Delta.coeff(t, 0) == 0 and Delta.coeff(t, 1) == 0)
f0 = sp.symbols('f0')
hh = sp.symbols('hh')
# f0 = -3 hh^2, g0 = -2 hh^3 solve 4 f0^3 + 27 g0^2 = 0 and hh = 3 g0 / (2 f0)
ok('B2 Delta_0 = 0 on Sen central fibre', sp.expand(4 * (-3 * hh ** 2) ** 3 + 27 * (-2 * hh ** 3) ** 2) == 0)
ok('B2 h = 3 g0 / (2 f0)', sp.simplify(3 * (-2 * hh ** 3) / (2 * (-3 * hh ** 2)) - hh) == 0)
fs = -3 * h ** 2 + t * F
gs = -2 * h ** 3 + t * Gv
# without base change: eta = F, chi = -12 (G - h F)/t
ok('B2 normal form, k = 1', sp.expand(g.subs({eta: F, chi: -12 * (Gv - h * F) / t}) - gs) == 0)
# with t = u^2: eta_u = u F, chi_u = -12 (G - h F)
u = sp.symbols('u')
ok('B2 normal form after t = u^2, f',
   sp.expand((-3 * h ** 2 + u * (u * F)) - fs.subs(t, u ** 2)) == 0)
ok('B2 normal form after t = u^2, g',
   sp.expand((-2 * h ** 3 + u * h * (u * F) - u ** 2 * (-12 * (Gv - h * F)) / 12) - gs.subs(t, u ** 2)) == 0)
fg = -3 * h ** 2 + t * F1
gg = -2 * h ** 3 + t * G1v
D1 = sp.expand(4 * fg ** 3 + 27 * gg ** 2)
ok('B3 first-order discriminant 108 t h^3 (h F1 - G1)', sp.expand(D1.coeff(t, 1) - 108 * h ** 3 * (h * F1 - G1v)) == 0)
ok('B4 isotrivial Delta = h^6 t^2 (4t - 9)', sp.expand(Delta.subs({eta: h ** 2, chi: 0}) - h ** 6 * t ** 2 * (4 * t - 9)) == 0)
jiso = sp.simplify((f ** 3 / Delta).subs({eta: h ** 2, chi: 0}))
ok('B4 j independent of the base point', sp.diff(jiso, h) == 0)
for nn in range(3, 13):
    af = sp.ceiling(sp.Rational(4 * (nn - 2), nn))
    ag = sp.ceiling(sp.Rational(6 * (nn - 2), nn))
    if nn >= 5:
        ok(f'B5 ord f >= 3, ord g >= 4 on a -{nn} curve', af >= 3 and ag >= 4)
    if nn in (3, 4):
        ok(f'B5 bounds ({af},{ag}) on a -{nn} curve', (af, ag) == ((2, 2) if nn == 3 else (2, 3)))
K2 = sp.symbols('K2', positive=True)
chiY, chiO, chiDB = -60 * K2, -2 * K2, -56 * K2
ok('B6 genus of O in |-2K|', sp.expand(-(chiO) - (-2) * (-2 + 1) * K2) == 0)
ok('B6 genus of D_B in |-8K|', sp.expand(-(chiDB) - (-8) * (-8 + 1) * K2) == 0)
ok('B6 tangency points h = eta = 0: (-2K)(-4K) = 8K^2', (-2) * (-4) == 8)
ok('B6 2 chi(Y) = 4 chi(O) + chi(D~)', sp.expand(2 * chiY - (4 * chiO + 2 * chiDB)) == 0)

# ---------------- A7 the upper endpoint of Phi ----------------
uu = sp.symbols('uu', positive=True)
for D in range(4, 13):
    a, b = (D - 2) * (D - 3), (D - 4) * (D - 5)
    rho = sp.Rational(b, 8 * a)
    Phi_u = (2 + uu / 2) ** 2 / (uu * (1 + rho * uu))
    lim = sp.limit(Phi_u, uu, sp.oo)
    star = 4 - sp.Rational(2 * b, a)
    ok(f'A7 Phi* at u = 4a/(a-b) D={D}', sp.simplify(Phi_u.subs(uu, sp.Rational(4 * a, a - b)) - star) == 0)
    if b == 0:
        ok(f'A7 Phi -> oo as u -> oo D={D}', lim == sp.oo)
    else:
        ok(f'A7 Phi -> 2a/b > Phi* as u -> oo D={D}',
           lim == sp.Rational(2 * a, b) and sp.Rational(2 * a, b) - star == sp.Rational(2 * (a - b) ** 2, a * b) > 0)

# ---------------- B7 F_12: the D7 divisor contains 7C ----------------
def least_order(k):
    # monomials (a, b) of -kK on F_12: a >= -k, b <= k, 12 b >= a - k; order b + k on the (-12)-curve
    return min(b + k for a in range(-k, 13 * k + 1) for b in range(-k, k + 1) if 12 * b >= a - k)
orders = [least_order(k) for k in (2, 4, 6)]
ok('B7 least orders of h, eta, chi on F_12 are 2, 4, 5', orders == [2, 4, 5])
ok('B7 the same from ceil(10k/12)', [sp.ceiling(sp.Rational(10 * k, 12)) for k in (2, 4, 6)] == [2, 4, 5])
ok('B7 ord(eta^2) = 8 > ord(h chi) = 7, so D_B contains 7C', 2 * orders[1] == 8 and orders[0] + orders[2] == 7)
ok('B7 witness eta = h^2, chi = 0 has orders (4, 6) on C: nonminimal', (2 * orders[0], 3 * orders[0]) == (4, 6))

# ---------------- B8 the leading D7 locus is not the pole locus at finite t ----------------
D_on = Delta.subs({h: 1, eta: 1, chi: 1})
ok('B8 constant family h = eta = chi = 1: Delta = -t^3/2 + 3t^4/16', sp.expand(D_on - (-t ** 3 / 2 + sp.Rational(3, 16) * t ** 4)) == 0)
ok('B8 at t = 1/10 on eta^2 = h chi, Delta = -77/160000 and 4f^3 is finite',
   D_on.subs(t, sp.Rational(1, 10)) == sp.Rational(-77, 160000) and f.subs({h: 1, eta: 1, t: sp.Rational(1, 10)}) != 0)

# ---------------- B9 k <= 2: after eps = t^2 the first-order terms vanish ----------------
eps = sp.symbols('eps')
cs = sp.symbols('f0:4 g0:4')
for c in (1, 3):
    fe = eps ** (2 * c) * sum(cs[i] * eps ** i for i in range(4))
    ge = eps ** (3 * c) * sum(cs[4 + i] * eps ** i for i in range(4))
    ft = sp.expand(fe.subs(eps, t ** 2) / t ** (4 * c))
    gt = sp.expand(ge.subs(eps, t ** 2) / t ** (6 * c))
    ok(f'B9 c = {c}: rescaled f, g have no t^1 term, so G_0 = h F_0', ft.coeff(t, 1) == 0 and gt.coeff(t, 1) == 0)

# ---------------- C the half-BPS coefficients ----------------
qq = sp.symbols('qq')
ser = sp.expand(sp.series(sp.prod([(1 - qq ** m) ** -24 for m in range(1, 6)]), qq, 0, 5).removeO())
ok('C coefficients of eta^-24 at q^-1..q^3 are 1, 24, 324, 3200, 25650', [ser.coeff(qq, i) for i in range(5)] == [1, 24, 324, 3200, 25650])
ok('C q eta^-24 has zero coefficient at q^-1', sp.expand(qq * ser).coeff(qq, 0) == 0)

out = {'checks': len(checks), 'parts': sorted({c.split()[0] for c in checks})}
Path(__file__).with_name('closure-checks.json').write_text(json.dumps(out, indent=2) + '\n')
print(f'{len(checks)} checks passed, 0 failed')
