#!/usr/bin/env python3
"""Exact verification of the identities behind the positivity theorem for the
Z_2 CHL model at electric charge N Q^2/2 = 1.

Notation (sigma = the magnetic period-2 variable of 1/Phi~_6, x = e^{pi i sigma},
so p = e^{2 pi i sigma} = x^2; y = e^{2 pi i v}).  All series are exact
(integers / Fractions), keyed by (exponent of x, exponent of y).

  A(sigma,v) = phi_{-2,1}(sigma,v),  B(sigma,v) = phi_{0,1}(sigma,v)
  D(tau)     = eta(tau)^8 eta(2tau)^8,  tau = sigma/2  (so D is a series in x)
  E2^(2)(tau) = 2E_2(2tau) - E_2(tau),  E_4(tau), E_4(2tau)=E_4(sigma)
  phi_1      = coefficient of e^{2 pi i rho} in 1/Phi~_6 (chamber R expansion)
  L1, L2     = first two coefficients of the exponential of the product, built
               directly from the twisted elliptic genera
  Y          = (1/2)[theta_3(tau,v)theta_3(tau)^3 + theta_4(tau,v)theta_4(tau)^3]
  calH       = sum H(4n-r^2) p^n y^r (Hurwitz class numbers, H(0) = -1/12)

Checks:
 (a) exponents of the rho^{-1} factor of the product: phi_{-1} = 1/(D A)
 (b) D * A * phi_1 = (1/2) L1^2 + L2        (product formula vs explicit genera)
 (c) G1 := (1/2)L1^2 + L2 - [13/36 B^2 - 5/9 E2^(2) A B + (4/45 E4(tau) + 19/180 E4(sigma)) A^2] = 0
 (d) G3 := E2^(2) B - (1/5)(6 E4(sigma) - E4(tau)) A - 12 Y = 0, and Y = sum r_3(2e - l^2) x^e y^l
 (e) Bringmann-Murthy, Section 3: A*(B^2/A - 144 Av[y/(1-y)^2] - E4 A + 288 calH) = 0   (level one, variable q)
Each identity is checked to order x^EMAX (default 40), far beyond the Sturm bounds
(x^4 for (c) on Gamma(4), x^1 for (d) on Gamma^0(2)).
"""
import sys
from fractions import Fraction as Fr
from collections import defaultdict
from math import isqrt
import siegel_dyons as SD

EMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40

def mul(a, b, emax):
    out = defaultdict(int)
    for (ea, pa), ca in a.items():
        for (eb, pb), cb in b.items():
            e = ea + eb
            if e <= emax:
                out[(e, pa + pb)] += ca * cb
    return {k: v for k, v in out.items() if v}

def add(*terms):
    out = defaultdict(int)
    for f, a in terms:
        for k, v in a.items():
            out[k] += f * v
    return {k: v for k, v in out.items() if v}

def trunc(a, emax):
    return {k: v for k, v in a.items() if k[0] <= emax}

def sigma_k(n, k):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)

def eis(k, emax, scale):
    """E_k(scale*tau) as series in x = q_tau, k in {2,4}"""
    c = {2: -24, 4: 240}[k]
    s = {(0, 0): 1}
    for n in range(1, emax // scale + 1):
        s[(n * scale, 0)] = c * sigma_k(n, k - 1)
    return s

def qprod(factors, emax):
    """prod_n (1 - x^{step n})^{pw} for (step, pw) in factors; 1-variable list"""
    s = [0] * (emax + 1); s[0] = 1
    for step, pw in factors:
        n = 1
        while step * n <= emax:
            d = step * n
            for _ in range(abs(pw)):
                if pw > 0:
                    for e in range(emax, d - 1, -1): s[e] -= s[e - d]
                else:
                    for e in range(d, emax + 1): s[e] += s[e - d]
            n += 1
    return s

def lift(s, shift=0):
    return {(e + shift, 0): c for e, c in enumerate(s) if c}

def phim21(emax, scale):
    """phi_{-2,1}(scale*tau', v) in x where q_{tau'} = x (scale=2: argument sigma)"""
    s = {(0, 1): 1, (0, 0): -2, (0, -1): 1}
    n = 1
    while scale * n <= emax:
        d = scale * n
        for f in ({(0, 0): 1, (d, 1): -1},) * 2 + ({(0, 0): 1, (d, -1): -1},) * 2:
            s = mul(s, f, emax)
        n += 1
    return mul(s, lift(qprod([(scale, -4)], emax)), emax)

def phi01_from_genus(emax, scale):
    """phi_{0,1}(scale*tau', v) = 4 * A_code, A_code = sum theta_i^2(z)/theta_i^2(0) (from siegel_dyons)"""
    A, _ = SD.AB(8 * (emax // scale) + 16)
    out = {}
    for (e8, y), v in A.items():
        if e8 % 8 == 0 and (e8 // 8) * scale <= emax:
            v4 = 4 * v
            assert Fr(v4).denominator == 1
            out[((e8 // 8) * scale, y)] = int(v4)
    return out

def hurwitz(nmax):
    """Hurwitz class numbers H(0..nmax) exactly"""
    H = [Fr(0)] * (nmax + 1)
    H[0] = Fr(-1, 12)
    for D in range(1, nmax + 1):
        if D % 4 not in (0, 3):
            continue
        tot = Fr(0)
        a = 1
        while 3 * a * a <= D:
            for b in range(-a, a + 1):
                if (b * b + D) % (4 * a): continue
                c = (b * b + D) // (4 * a)
                if c < a: continue
                if b < 0 and (a == c or -b == a): continue   # reduced: |b|<=a<=c, b>=0 if |b|=a or a=c
                w = Fr(1)
                if a == b == c: w = Fr(1, 3)                  # x^2+xy+y^2 multiples
                elif b == 0 and a == c: w = Fr(1, 2)          # x^2+y^2 multiples
                tot += w
            a += 1
        H[D] = tot
    return H

def r3(n):
    if n < 0: return 0
    cnt = 0
    for a in range(-isqrt(n), isqrt(n) + 1):
        for b in range(-isqrt(n), isqrt(n) + 1):
            c2 = n - a * a - b * b
            if c2 < 0: continue
            c = isqrt(c2)
            if c * c == c2: cnt += 1 if c == 0 else 2
    return cnt

def main():
    E = EMAX
    report = {}
    # ---------- twisted genera (paper eq. for the 2A CHL model), in q^{1/8} units ----------
    Kq = E // 2 + 2                      # need q-powers up to E/2 for L2 at x^E (q^n -> x^n), up to E/2 for L1
    F = SD.twisted_genera('chl2', 8 * (E + 4) + 16)
    def coeff(r, s, n_half, j):
        """coefficient of q^{n_half/2} y^j in F^{(r,s)} (n_half integer)"""
        return F.get((r, s), {}).get((4 * n_half, j), 0)
    # (a) exponents of the l=0 factors: c_r(j) = sum_s [F^{(r,s)}]_{q^0 y^j}
    ca = {(r, j): sum(coeff(r, s, 0, j) for s in (0, 1)) for r in (0, 1) for j in range(-3, 4)}
    ok_a = (ca[(0, 0)] == 12 and ca[(0, 1)] == ca[(0, -1)] == 2 and ca[(1, 0)] == 8
            and ca[(1, 1)] == ca[(1, -1)] == 0 and all(ca[(r, j)] == 0 for r in (0, 1) for j in (2, 3, -2, -3)))
    report['a_phi_minus1_equals_1_over_DA'] = ok_a
    # ---------- L1, L2 from the genera ----------
    # L1 = sum_{r,s} (-1)^s sum_{k' in Z+r/2} c^{(r,s)}(q^{k'} y^j) p^{k'} y^j,  p^{k'} = x^{2k'}
    L1 = defaultdict(Fr)
    for r in (0, 1):
        for s in (0, 1):
            for (e8, j), v in F.get((r, s), {}).items():
                if e8 % 4: continue
                nh = e8 // 4              # 2k'
                if nh % 2 != r % 2: continue   # k' in Z + r/2
                if nh <= E:
                    L1[(nh, j)] += (-1) ** s * v
    L1 = {k: v for k, v in L1.items() if v}
    # L2a = sum_{r,s} sum_{k' in Z+r/2} c^{(r,s)}(q^{2k'} y^j) p^{k'} y^j: q-power n=2k' integer, n = r mod 2
    L2 = defaultdict(Fr)
    for r in (0, 1):
        for s in (0, 1):
            for (e8, j), v in F.get((r, s), {}).items():
                if e8 % 8: continue
                n = e8 // 8
                if n % 2 != r % 2: continue
                if n <= E:
                    L2[(n, j)] += v
    for (e, j), v in L1.items():       # + (1/2) L1(2 sigma, 2 v)
        if 2 * e <= E:
            L2[(2 * e, 2 * j)] += Fr(1, 2) * v
    L2 = {k: v for k, v in L2.items() if v}
    P = add((Fr(1, 2), mul(L1, L1, E)), (1, L2))
    assert all(Fr(v).denominator == 1 for v in P.values())
    # (b) compare with the product formula: D * A * phi_1 (chamber R), phi_1 from siegel_dyons
    Kb = min(E, 24)
    tab = SD_coeffs(Kb + 2, 3 * Kb + 12)
    phi1 = {(k2, p): v for (l, k2, p), v in tab.items() if l == 1}
    Dser = lift(qprod([(1, 8), (2, 8)], Kb + 2), 1)
    A_s = phim21(Kb + 2, 2)
    lhs = mul(mul(Dser, phi1, Kb), A_s, Kb)
    W = 2 * Kb + 4
    diffb = {k: v for k, v in add((1, lhs), (-1, trunc(P, Kb))).items() if abs(k[1]) <= W}
    report['b_DAphi1_equals_half_L1sq_plus_L2_to_x^%d' % Kb] = (len(diffb) == 0)
    # (c) G1 = 0
    B_s = phi01_from_genus(E, 2)
    A_s = phim21(E, 2)
    E22 = add((2, eis(2, E, 2)), (-1, eis(2, E, 1)))
    gam = add((Fr(4, 45), eis(4, E, 1)), (Fr(19, 180), eis(4, E, 2)))
    rhs = add((Fr(13, 36), mul(B_s, B_s, E)), (Fr(-5, 9), mul(mul(E22, A_s, E), B_s, E)), (1, mul(gam, mul(A_s, A_s, E), E)))
    G1 = add((1, P), (-1, rhs))
    report['c_G1_zero_to_x^%d' % E] = (len(G1) == 0)
    # (d) G3 = 0 and Y = r_3
    # theta_3(tau,v) = sum x^{n^2/2} y^n : use X = x^{1/2}
    th3v = {(n * n, n): 1 for n in range(-isqrt(2 * E) - 1, isqrt(2 * E) + 2) if n * n <= 2 * E}
    th4v = {(n * n, n): (1 if n % 2 == 0 else -1) for n in range(-isqrt(2 * E) - 1, isqrt(2 * E) + 2) if n * n <= 2 * E}
    def at0(s):
        o = defaultdict(int)
        for (e, p), v in s.items(): o[(e, 0)] += v
        return {k: v for k, v in o.items() if v}
    t3 = at0(th3v); t4 = at0(th4v)
    E2X = 2 * E
    Y2 = add((Fr(1, 2), mul(th3v, mul(t3, mul(t3, t3, E2X), E2X), E2X)),
             (Fr(1, 2), mul(th4v, mul(t4, mul(t4, t4, E2X), E2X), E2X)))
    Y = {}
    okY = True
    for (eX, p), v in Y2.items():
        if eX % 2: okY = False; continue
        Y[(eX // 2, p)] = v
    okY = okY and all(v == r3(2 * e - p * p) for (e, p), v in Y.items())
    okY = okY and all(Y.get((e, p), 0) == r3(2 * e - p * p) for e in range(E + 1) for p in range(-isqrt(2 * e), isqrt(2 * e) + 1))
    report['d_Y_coefficients_are_r3(2e-l^2)'] = okY
    G3 = add((1, mul(E22, B_s, E)), (Fr(-6, 5), mul(eis(4, E, 2), A_s, E)), (Fr(1, 5), mul(eis(4, E, 1), A_s, E)), (-12, Y))
    report['d_G3_zero_to_x^%d' % E] = (len(G3) == 0)
    # (e) level one identity of Bringmann-Murthy, Section 3, in q (scale 1): A*(B^2/A) = B^2 ; Av term times A computed in the annulus |q|<|y|<|q|^-1
    Eq = E // 2
    A1 = phim21(Eq, 1); B1 = phi01_from_genus(Eq, 1)
    H = hurwitz(4 * Eq + 8)
    calH = {}
    for n in range(0, Eq + 1):
        for r in range(-isqrt(4 * n + 1) - 1, isqrt(4 * n + 1) + 2):
            D = 4 * n - r * r
            if D >= 0 and H[D] != 0: calH[(n, r)] = H[D]
    # Av[y/(1-y)^2] * A :  s = 0 term gives A * y/(1-y)^2 = prod part of A
    prodA = mul(A1, {}, Eq)  # placeholder
    pa = {(0, 0): 1}
    n = 1
    while n <= Eq:
        for f in ({(0, 0): 1, (n, 1): -1},) * 2 + ({(0, 0): 1, (n, -1): -1},) * 2:
            pa = mul(pa, f, Eq)
        n += 1
    pa = mul(pa, lift(qprod([(1, -4)], Eq)), Eq)     # A = (y-2+1/y) * pa, and A*y/(1-y)^2 = pa
    AvA = dict(pa)
    for s in range(1, isqrt(Eq) + 2):
        # s>0: q^{s^2+s} y^{2s+1} sum_{k>=1} k (q^s y)^{k-1}
        t = {}
        k = 1
        while s * s + s + s * (k - 1) <= Eq:
            t[(s * s + s + s * (k - 1), 2 * s + 1 + (k - 1))] = k; k += 1
        AvA = add((1, AvA), (1, mul(t, A1, Eq)))
        # -s<0: q^{s^2-s} y^{-2s+1} /(1-q^{-s} y)^2 = q^{s^2-s} y^{-2s+1} q^{2s} y^{-2} sum k (q^{s} y^{-1})^{k-1}
        t = {}
        k = 1
        while s * s + s + s * (k - 1) <= Eq:
            t[(s * s + s + s * (k - 1), -2 * s - 1 - (k - 1))] = k; k += 1
        AvA = add((1, AvA), (1, mul(t, A1, Eq)))
    lhs_e = mul(B1, B1, Eq)
    rhs_e = add((144, AvA), (1, mul(mul(eis(4, Eq, 1), A1, Eq), A1, Eq)), (-288, mul(calH, A1, Eq)))
    De = add((1, lhs_e), (-1, rhs_e))
    report['e_BM_B2_over_A_finite_part_to_q^%d' % Eq] = (len(De) == 0)
    for k, v in report.items():
        print(k, v)
    return all(report.values())

_cache = {}
def SD_coeffs(K, P):
    import io, contextlib
    if (K, P) not in _cache:
        src = SD.coefficient_tables
        # allow K2 = -1 and L = -1 in the output (the filter in siegel_dyons drops them)
        import inspect, textwrap
        code = inspect.getsource(src).replace("if 0 <= L <= Lmax and 0 <= K2 <= Kmax2 and abs(p) <= Pmax:",
                                              "if -1 <= L <= Lmax and -1 <= K2 <= Kmax2 and abs(p) <= Pmax:")
        ns = dict(SD.__dict__)
        exec(code, ns)
        _cache[(K, P)] = ns['coefficient_tables']('chl2', 2, 1, K, P, 'R')
    return _cache[(K, P)]

if __name__ == '__main__':
    ok = main()
    print('ALL IDENTITIES VERIFIED' if ok else 'FAILURE')
    sys.exit(0 if ok else 1)
