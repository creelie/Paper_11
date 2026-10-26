#!/usr/bin/env python3
"""Fourier coefficients of 1/Phi_k for quarter-BPS dyon counting in N=4 string theory.

Exact integer arithmetic.  Reproduces (i) the heterotic/K3xT^2 table of Sen (1/Phi_10),
(ii) the wall-crossing jump between the two expansions in e^{+-2 pi i v}, and
(iii) the toroidal Z_2 orbifold table of Chattopadhyaya-David (JHEP 05 (2019) 005, Table 17-19),
whose known negative full-index entries require the hair and horizon distinction.

Conventions (Sen, JHEP 0705:039; David-Jatkar-Sen, JHEP 0606:064):
  Phi(rho,sigma,v) = e^{2 pi i (rho + sigma/N + v)} prod_{b,r} prod_{k' in Z + r/N, l in Z, j in 2Z+b}
                     (1 - e^{2 pi i (k' sigma + l rho + j v)})^{ sum_s e^{2 pi i s l/N} c_b^{(r,s)}(4 k' l - j^2) }
  with k', l >= 0 and j < 0 when k' = l = 0,
  F^{(r,s)}(tau,z) = sum_{b} sum_{j in 2Z+b, n in Z/N} c_b^{(r,s)}(4n - j^2) q^n y^j,
  d(Q,P) = (1/N)(-1)^{Q.P+1} * [coefficient of e^{2 pi i (rho N Q^2/2 + sigma P^2/(2N) + v Q.P)} in 1/Phi],
  the coefficient being read off in the chamber R (expansion of the (k'=l=0) factor in powers of e^{-2 pi i v}).
"""
from fractions import Fraction
from collections import defaultdict
import sys

# ---------- series in x = q^{1/8} (integer exponents) and y (Laurent) ----------
def smul(a, b, xmax):
    out = defaultdict(int)
    for (ea, ya), ca in a.items():
        for (eb, yb), cb in b.items():
            e = ea + eb
            if e <= xmax:
                out[(e, ya + yb)] += ca * cb
    return {k: v for k, v in out.items() if v != 0}

def sadd(a, b, fa=1, fb=1):
    out = defaultdict(Fraction)
    for k, v in a.items(): out[k] += Fraction(fa) * v
    for k, v in b.items(): out[k] += Fraction(fb) * v
    return {k: v for k, v in out.items() if v != 0}

def sinv(a, xmax):
    """inverse of a series whose lowest x-power term is a monomial c*x^e0*y^y0."""
    e0 = min(e for (e, y) in a)
    lead = [(e, y) for (e, y) in a if e == e0]
    assert len(lead) == 1, "leading term must be a monomial"
    (e0, y0) = lead[0]
    c0 = a[(e0, y0)]
    # a = c0 x^e0 y^y0 (1 + r),  1/a = (1/c0) x^-e0 y^-y0 sum (-r)^n
    r = {(e - e0, y - y0): Fraction(v, 1) / c0 for (e, y), v in a.items() if (e, y) != (e0, y0)}
    inv = {(0, 0): Fraction(1)}
    term = {(0, 0): Fraction(1)}
    for _ in range(xmax + 2):
        term = smul(term, {k: -v for k, v in r.items()}, xmax + e0)
        if not term: break
        inv = sadd(inv, term)
    return {(e - e0, y - y0): v / c0 for (e, y), v in inv.items() if e - e0 <= xmax}

def theta(kind, xmax):
    """theta_i(tau,z) as series in x=q^{1/8}, y^{1/2} units: we return exponents of y in HALF units (ky)."""
    s = defaultdict(Fraction)
    n = -40
    while n <= 40:
        if kind in (1, 2):
            e8 = (2 * n + 1) ** 2  # q^{(n+1/2)^2/2} = x^{(2n+1)^2}
            ky = 2 * n + 1
            if kind == 1: c = Fraction((-1) ** n) ; # -i factor dropped: theta_1^2 gets an overall -1 (handled below)
            else: c = Fraction(1)
        else:
            e8 = 4 * n * n  # q^{n^2/2} = x^{4 n^2}
            ky = 2 * n
            c = Fraction(1) if kind == 3 else Fraction((-1) ** n)
        if e8 <= xmax: s[(e8, ky)] += c
        n += 1
    return dict(s)

def sq_half(s, xmax):
    """square of a series in half-integer y powers; result in integer y powers."""
    out = smul(s, s, xmax)
    res = {}
    for (e, ky), v in out.items():
        assert ky % 2 == 0
        res[(e, ky // 2)] = v
    return res

def eta_pow(p, xmax):
    """eta(tau)^p / q^{p/24} as series in x = q^{1/8}: prod (1 - q^n)^p, with q = x^8."""
    s = {(0, 0): Fraction(1)}
    n = 1
    while 8 * n <= xmax:
        f = {(0, 0): Fraction(1), (8 * n, 0): Fraction(-1)}
        for _ in range(abs(p)):
            s = smul(s, f, xmax) if p > 0 else s
        if p < 0:
            finv = sinv(f, xmax)
            for _ in range(-p): s = smul(s, finv, xmax)
        n += 1
    return s

def E2_scaled(scale_num, scale_den, shift, xmax):
    """E_2((scale_num/scale_den) tau + shift/ scale_den... ) : we implement E2(tau), E2(2tau), E2(tau/2), E2((tau+1)/2).
    E2(tau) = 1 - 24 sum sigma(n) q^n.  E2(tau/2): q^{n/2} = x^{4n}; E2((tau+1)/2): (-1)^n x^{4n}."""
    s = {(0, 0): Fraction(1)}
    n = 1
    while True:
        e = 8 * n * scale_num // scale_den  # exponent of x for q^{n*scale}
        if e > xmax: break
        sig = sum(d for d in range(1, n + 1) if n % d == 0)
        s[(e, 0)] = s.get((e, 0), 0) + Fraction(-24 * sig * ((-1) ** (n * shift)))
        n += 1
    return s

def AB(xmax):
    th1 = theta(1, xmax); th2 = theta(2, xmax); th3 = theta(3, xmax); th4 = theta(4, xmax)
    th1_0 = None
    def at_zero(s):
        out = defaultdict(Fraction)
        for (e, ky), v in s.items(): out[(e, 0)] += v
        return dict(out)
    A = {}
    for th in (th2, th3, th4):
        num = sq_half(th, xmax)
        den = sq_half(at_zero(th), xmax)
        A = sadd(A, smul(num, sinv(den, xmax), xmax))
    # B = eta^{-6} theta_1^2 ; theta_1 = -i sum (-1)^n q^{(n+1/2)^2/2} y^{n+1/2}  => theta_1^2 = -(sum ...)^2
    th1sq = {k: -v for k, v in sq_half(th1, xmax).items()}
    # eta^{-6} = q^{-1/4} prod(1-q^n)^{-6} ; q^{-1/4} = x^{-2}
    etam6 = eta_pow(-6, xmax)
    B = {(e - 2, y): v for (e, y), v in smul(th1sq, etam6, xmax + 2).items()}
    return A, B

def twisted_genera(model, xmax):
    """return dict (r,s) -> F^{(r,s)} series in (x,y); x = q^{1/8}."""
    A, B = AB(xmax)
    # In David-Jatkar-Sen and Chattopadhyaya-David the symbol E_2 in the twisted genera is the
    # weight-2 form of Gamma_0(2),  E_N(tau) = N E(N tau) - E(tau)  (N = 2), E the Eisenstein series.
    Eis = lambda a, b, sh: E2_scaled(a, b, sh, xmax)
    E2 = sadd({k: 2 * v for k, v in Eis(2, 1, 0).items()}, Eis(1, 1, 0), 1, Fraction(-1))          # E_2(tau)
    E2_half = sadd({k: 2 * v for k, v in Eis(1, 1, 0).items()}, Eis(1, 2, 0), 1, Fraction(-1))     # E_2(tau/2)
    E2_half_shift = sadd({k: 2 * v for k, v in Eis(1, 1, 0).items()}, Eis(1, 2, 1), 1, Fraction(-1))  # E_2((tau+1)/2)
    F = {}
    if model == 'het':   # N = 1, K3 elliptic genus
        F[(0, 0)] = {k: 8 * v for k, v in A.items()}
    elif model == 'chl2':  # 2A CHL orbifold (David-Jatkar-Sen)
        F[(0, 0)] = {k: 4 * v for k, v in A.items()}
        F[(0, 1)] = sadd({k: Fraction(4, 3) * v for k, v in A.items()}, smul(B, E2, xmax), 1, Fraction(-2, 3))
        F[(1, 0)] = sadd({k: Fraction(4, 3) * v for k, v in A.items()}, smul(B, E2_half, xmax), 1, Fraction(1, 3))
        F[(1, 1)] = sadd({k: Fraction(4, 3) * v for k, v in A.items()}, smul(B, E2_half_shift, xmax), 1, Fraction(1, 3))
    elif model == 'tor2':  # Z_2 toroidal orbifold of type II on T^6 (Chattopadhyaya-David eq. 4.2)
        F[(0, 0)] = {}
        F[(0, 1)] = sadd({k: Fraction(8, 3) * v for k, v in A.items()}, smul(B, E2, xmax), 1, Fraction(-4, 3))
        F[(1, 0)] = sadd({k: Fraction(8, 3) * v for k, v in A.items()}, smul(B, E2_half, xmax), 1, Fraction(2, 3))
        F[(1, 1)] = sadd({k: Fraction(8, 3) * v for k, v in A.items()}, smul(B, E2_half_shift, xmax), 1, Fraction(2, 3))
    return F

def coefficient_tables(model, N, Lmax, Kmax2, Pmax, chamber='R'):
    """Return dict (l, k2, p) -> integer coefficient of e^{2 pi i (l rho + (k2/N) sigma + p v)} in 1/Phi,
    for 0<= l <= Lmax (rho exponent), 0 <= k2 <= Kmax2 (sigma exponent times N), |p| <= Pmax.
    The expansion of the (k'=l=0) factor is in powers of e^{-2 pi i v} (chamber R) or e^{+2 pi i v} (chamber L).
    The prefactor e^{-2 pi i (rho + sigma/N + v)} is included (exponents shifted)."""
    # we need F coefficients c^{(r,s)}(4 n - j^2) for n up to Kmax2/N * Lmax ... build genera to sufficient order
    # prefactor: e^{2 pi i(rho + sigma/N + v)} for K3-based models (shift k' by 1/N, i.e. k2 by 1);
    # the toroidal orbifolds have e^{2 pi i(rho + v)} only (Chattopadhyaya-David eq. 4.1).
    kshift = 0 if model == 'tor2' else 1
    nmax = Fraction(Kmax2 + kshift, N) * (Lmax + 1)
    xmax = int(8 * nmax) + 16
    F = twisted_genera(model, xmax)
    # c[(r,s)][(n8, j)] with n8 = exponent of x (q^{n} = x^{8n})
    # exponent for factor (k', l, j): sum_s e^{2 pi i s l/N} c^{(r,s)}(4 k' l - j^2), r = N k' mod N
    import cmath
    def expo(k2, l, j):
        r = k2 % N
        kp = Fraction(k2, N)
        n = kp * l           # q^n with n = k' l
        n8 = int(8 * n)
        assert 8 * n == n8
        assert n8 <= xmax - 8, (k2, l, n8, xmax)
        tot = Fraction(0)
        for s in range(N):
            Frs = F.get((r, s), {})
            c = Frs.get((n8, j), 0)
            if c == 0: continue
            ph = s * l % N
            # e^{2 pi i ph / N}: for N = 1,2 this is +-1 ; general N would need cyclotomic arithmetic
            assert N in (1, 2)
            tot += c * (1 if ph == 0 else -1)
        assert tot.denominator == 1, (k2, l, j, tot)
        return int(tot)
    # series keyed by (l, k2, p)
    ser = {(0, 0, 0): 1}
    def mul(a, b):
        out = defaultdict(int)
        for (la, ka, pa), ca in a.items():
            for (lb, kb, pb), cb in b.items():
                l, k, p = la + lb, ka + kb, pa + pb
                if l <= Lmax + 1 and k <= Kmax2 + kshift:
                    out[(l, k, p)] += ca * cb
        return {k: v for k, v in out.items() if v}
    def binom_series(l, k2, j, c):
        """(1 - X)^{-c} with X = e^{2 pi i (k' sigma + l rho + j v)}, truncated."""
        # number of powers allowed
        if l > 0: nmaxp = (Lmax + 1) // l
        else: nmaxp = (Kmax2 + kshift) // k2
        if k2 > 0: nmaxp = min(nmaxp, (Kmax2 + kshift) // k2)
        out = {(0, 0, 0): 1}
        if c == 0: return out
        # generalized binomial: coefficient of X^n in (1-X)^{-c} = binom(c+n-1, n) (for c>0), for c<0 it is (1-X)^{|c|}
        from math import comb
        for n in range(1, nmaxp + 1):
            if c > 0: coef = comb(c + n - 1, n)
            else:
                if n > -c: break
                coef = (-1) ** n * comb(-c, n)
            out[(l * n, k2 * n, j * n)] = coef
        return out
    # factors with (k', l) != (0, 0)
    for l in range(0, Lmax + 2):
        for k2 in range(0, Kmax2 + kshift + 1):
            if l == 0 and k2 == 0: continue
            n = Fraction(k2, N) * l
            jmax = int((4 * n + 1) ** 0.5) + 2
            for j in range(-jmax, jmax + 1):
                c = expo(k2, l, j)
                if c == 0: continue
                ser = mul(ser, binom_series(l, k2, j, c))
    # the (k'=l=0) factors: j < 0 : (1 - e^{2 pi i j v})^{-c(-j^2)}, c(-j^2) nonzero only for j = -1 (index 1)
    c1 = expo(0, 0, -1)
    assert expo(0, 0, -2) == 0
    # prefactor e^{-2 pi i (rho + sigma/N + v)} -> shift (l,k2,p) by (-1,-N,-1)
    # chamber R: e^{-2 pi i v}(1 - e^{-2 pi i v})^{-c1} = sum_{n>=0} binom(c1+n-1,n) e^{-2 pi i (n+1) v}
    # chamber L: e^{-2 pi i v}(1 - e^{-2 pi i v})^{-c1} = (-1)^{c1} e^{2 pi i (c1-1) v} (1 - e^{2 pi i v})^{-c1}
    from math import comb
    prange = max(abs(p) for (_, _, p) in ser) + Pmax + 2
    vser = {}
    for n in range(0, prange + 1):
        coef = comb(c1 + n - 1, n)
        if chamber == 'R': vser[(0, 0, -(n + 1))] = coef
        else: vser[(0, 0, c1 - 1 + n)] = (-1) ** c1 * coef
    full = mul(ser, vser)
    out = {}
    for (l, k2, p), v in full.items():
        L, K2 = l - 1, k2 - kshift
        if 0 <= L <= Lmax and 0 <= K2 <= Kmax2 and abs(p) <= Pmax:
            out[(L, K2, p)] = v
    return out

def d_index(coef, N, l, k2, p):
    """d(Q,P) with N Q^2/2 = l, P^2/(2N) = k2/N (i.e. P^2 = 2 k2), Q.P = p."""
    v = coef.get((l, k2, p), 0)
    # the 1/N in the integral is cancelled by the sigma-period N of the contour, so d is the
    # signed coefficient itself (Chattopadhyaya-David eq. 2.1-2.2)
    sgn = -1 if (p + 1) % 2 else 1
    return sgn * v

if __name__ == '__main__':
    import json
    results = {}
    # ---- heterotic on T^6 / IIB on K3 x T^2 : 1/Phi_10 ----
    cR = coefficient_tables('het', 1, 3, 3, 6, 'R')
    cL = coefficient_tables('het', 1, 3, 3, 6, 'L')
    print("1/Phi_10, chamber R (attractor chamber for Q.P>=0 etc.): rows Q^2/2, cols Q.P, blocks P^2/2")
    tabR = {}
    for n in range(1, 4):      # P^2/2
        for m in range(1, 4):  # Q^2/2
            row = [int(d_index(cR, 1, m, n, p)) for p in range(-2, 5)]
            tabR[(2 * m, 2 * n)] = row
            print(f"Q^2={2*m}, P^2={2*n}:", dict(zip(range(-2, 5), row)))
    # wall-crossing check: d_R - d_L = (-1)^{Q.P+1} * |Q.P| * d_el(Q) d_mag(P) in magnitude: Sen's formula
    # d_h(n) = coefficient of q^{n} in 1/eta^24 (q^{-1} prefactor shifted): 1/eta^24 = q^{-1} sum d(n) q^n
    def eta24inv(nmax):
        s = eta_pow(-24, 8 * nmax + 8)
        return {e // 8: int(v) for (e, y), v in s.items() if e % 8 == 0}
    dh = eta24inv(5)
    print("1/eta^24 coefficients d(n), n = Q^2/2 + 1:", {k: dh[k] for k in range(0, 5)})
    wc_ok = True
    for n in range(1, 4):
        for m in range(1, 4):
            for p in range(-2, 5):
                dR = d_index(cR, 1, m, n, p); dL = d_index(cL, 1, m, n, p)
                jump = dR - dL
                pred = (-1 if (p + 1) % 2 else 1) * abs(p) * dh[m + 1] * dh[n + 1]
                # Sen: d_> - d_< = -Q.P d_el d_mag after the (-1)^{Q.P+1}, i.e. jump = (-1)^{Q.P+1}|Q.P| d_el d_mag up to the orientation
                ok = abs(jump) == abs(pred)
                wc_ok &= ok
                if not ok: print("WALL-CROSSING MISMATCH", m, n, p, jump, pred)
    print("wall-crossing jump |d_R - d_L| = |Q.P| d_el(Q) d_mag(P) for all Q^2,P^2 in {2,4,6}, |Q.P|<=4:", wc_ok)
    results['Phi10_R'] = {f"{k[0]},{k[1]}": v for k, v in tabR.items()}
    results['Phi10_wall_crossing_check'] = wc_ok
    # ---- Z_2 toroidal orbifold of type II on T^6 ----
    cT = coefficient_tables('tor2', 2, 12, 4, 4, 'R')
    print("\nZ_2 toroidal orbifold (Chattopadhyaya-David Tables 17-19): d(Q,P), rows Q^2, cols P^2 in {0,2,4}; Q.P = 0,1,2")
    tabT = {}
    for p in range(0, 3):
        print(f"Q.P = {p}")
        for Q2 in range(0, 13):
            row = []
            for P2 in (0, 2, 4):
                # N Q^2/2 = Q^2 (N=2) -> l = Q2 ; P^2/(2N) = P2/4 -> k2 = P2/2
                v = cT.get((Q2, P2 // 2, p), 0)
                row.append((-1 if (p + 1) % 2 else 1) * v)   # without the 1/N factor; see text for the normalisation
            tabT[(p, Q2)] = row
            print(f"  Q^2={Q2}: ", row)
    results['tor2_R_no_1_over_N'] = {f"{k[0]},{k[1]}": v for k, v in tabT.items()}
    # ---- 2A CHL orbifold: a positivity spot check in the chamber R ----
    cC = coefficient_tables('chl2', 2, 6, 3, 3, 'R')
    print("\n2A CHL orbifold (heterotic on T^6 / Z_2): d(Q,P), rows Q^2 (in units where N Q^2/2 = Q^2), cols P^2 in {0,2,4}")
    tabC = {}
    for p in range(0, 3):
        print(f"Q.P = {p}")
        for Q2 in range(0, 7):
            row = [str(d_index(cC, 2, Q2, P2 // 2, p)) for P2 in (0, 2, 4)]
            tabC[(p, Q2)] = row
            print(f"  Q^2={Q2}: ", row)
    results['chl2_R'] = {f"{k[0]},{k[1]}": v for k, v in tabC.items()}
    # single-centred region R for Gamma_0(2) (Sen 1008.4209 eq. 3.9 = Chattopadhyaya-David eq. 4.7)
    def single_centred_Z2(q2, p2, qp):
        return (q2 > 0 and p2 > 0 and q2 * p2 > qp * qp and 2 * q2 + p2 - 3 * qp >= 0
                and 2 * q2 >= qp and p2 >= qp and qp >= 0)
    viol = []
    for l in range(0, 13):
        for k2 in range(0, 5):
            for p in range(0, 5):
                q2, p2 = l, 2 * k2          # N Q^2/2 = l with N = 2, P^2 = 2 k2
                if single_centred_Z2(q2, p2, p):
                    d = d_index(cT, 2, l, k2, p)
                    if d < 0: viol.append((q2, p2, p, d))
    print("\nZ_2 toroidal orbifold: single-centred charges (conditions 4.7) with negative index, Q^2<=12, P^2<=8, Q.P<=4:")
    for v in viol: print("  Q^2=%d P^2=%d Q.P=%d  d=%d" % v)
    results['tor2_negative_single_centred'] = viol
    # the same scan for heterotic and 2A CHL: no violation expected
    def single_centred_het(q2, p2, qp):
        return q2 > 0 and p2 > 0 and q2 * p2 > qp * qp and 0 <= qp <= q2 and qp <= p2
    for name, model, N, cond in (('het', 'het', 1, single_centred_het), ('chl2', 'chl2', 2, single_centred_Z2)):
        c = coefficient_tables(model, N, 6, 3, 4, 'R')
        bad = []
        for l in range(0, 7):
            for k2 in range(0, 4):
                for p in range(0, 5):
                    q2, p2 = Fraction(2 * l, N), 2 * k2
                    if cond(q2, p2, p) and d_index(c, N, l, k2, p) <= 0: bad.append((str(q2), p2, p))
        print("%s: single-centred charges with non-positive index in the scanned range: %s" % (name, bad))
        results[name + '_nonpositive_single_centred'] = bad
    json.dump(results, open(sys.argv[1] if len(sys.argv) > 1 else 'siegel_dyons_out.json', 'w'), indent=1)
