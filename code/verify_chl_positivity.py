"""Exact checks for Theorem 5.2: positivity of the single-centred index of the
Z_2 CHL model at P^2 = 2.

Series are in x = q^{1/2}, q = e^{2 pi i tau}, and y = e^{2 pi i z}; the Fourier-Jacobi
coefficient psi_1 of 1/Phi~_6 has q-exponent Q^2 and y-exponent Q.P.
  C1  psi_1 from the product formula (siegel_dyons.py) minus the polar part 52 AL/g
      equals the closed form (10), for Q^2 <= 28 and |Q.P| <= 7;
  C2  the theta identity A(2tau,z) = 3 T(tau,z) - B(2tau,z) E(tau)/4 and the identity
      E(tau)^2 = (E_4(tau) + 4 E_4(2tau))/5, which give the form (11);
  C3  G = E_4(x)/36 - 49 E_4(x^2)/180 + 11 E_4(x^4)/45, to order x^K, and the
      Sturm bound (coefficients x^0, x^1, x^2) that makes the identity exact;
  C4  S = -j_0 + j_1 = theta_3 G / eta(2tau)^6 on the theta components of (11);
  C5  the closed form of the coefficients g(j) of G;
  C6  the coefficients c(k) of theta_3 G are positive for k <= K, and the bound
      (14) is positive for k >= 8;
  C7  the inequalities used for n_0 and n_1, with H(N) <= N;
  C8  the resulting signs (-1)^{Q.P+1} c^F(Q^2, Q.P) > 0 for 2 Q^2 >= (Q.P)^2, Q^2 <= 28.
Run: python3 verify_chl_positivity.py  (standard library only).
"""
import json
import math
from collections import defaultdict
from fractions import Fraction as Fr
from pathlib import Path

import siegel_dyons as SD

checks = []


def ok(name, cond):
    assert cond, name
    checks.append(name)


L = 30          # largest Q^2
P = 9           # largest |Q.P|
K = 120         # x-order of the one-variable checks


# ---------- one-variable series in x ----------
def m1(a, b, n):
    c = [Fr(0)] * (n + 1)
    for i, u in enumerate(a[:n + 1]):
        if u:
            for j in range(n + 1 - i):
                if b[j]:
                    c[i + j] += u * b[j]
    return c


def pprod(factors, n):
    """prod_{r>=1} prod_{(s,e)} (1 - x^{s r})^e, to order x^n"""
    c = [Fr(0)] * (n + 1)
    c[0] = Fr(1)
    for s, e in factors:
        for r in range(1, n // s + 1):
            st = s * r
            for _ in range(abs(e)):
                if e > 0:
                    for i in range(n, st - 1, -1):
                        c[i] -= c[i - st]
                else:
                    for i in range(st, n + 1):
                        c[i] += c[i - st]
    return c


def sig(k, n):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)


def sig_odd(n):
    return sum(d for d in range(1, n + 1) if n % d == 0 and d % 2)


def hurwitz(N):
    """Hurwitz class number H(N), H(0) = -1/12"""
    if N < 0:
        return Fr(0)
    if N == 0:
        return Fr(-1, 12)
    if N % 4 in (1, 2):
        return Fr(0)
    tot = Fr(0)
    b = N % 2
    while 3 * b * b <= N:
        ac = (b * b + N) // 4
        a = max(b, 1)
        while a * a <= ac:
            if ac % a == 0:
                c = ac // a
                if a <= c:
                    w = Fr(1, 3) if a == b == c else (Fr(1, 2) if b == 0 and a == c else Fr(1))
                    tot += w if (b == 0 or a == b or a == c) else 2 * w
            a += 1
        b += 2
    return tot


def E4x(s, n):
    c = [Fr(0)] * (n + 1)
    c[0] = Fr(1)
    for r in range(1, n // s + 1):
        c[s * r] = Fr(240 * sig(3, r))
    return c


# ---------- two-variable series: keys (x-exponent, y-exponent) ----------
def m2(a, b, n):
    o = defaultdict(Fr)
    for (ea, ya), ca in a.items():
        for (eb, yb), cb in b.items():
            if ea + eb <= n:
                o[(ea + eb, ya + yb)] += ca * cb
    return {k: v for k, v in o.items() if v}


def lift(c):
    return {(i, 0): v for i, v in enumerate(c) if v}


def at_z0(a):
    o = defaultdict(Fr)
    for (e, y), v in a.items():
        o[e] += v
    return [o.get(i, Fr(0)) for i in range(max(o) + 1)]


def inv1(c, n):
    """inverse of a one-variable series with c[0] != 0"""
    b = [Fr(0)] * (n + 1)
    b[0] = 1 / c[0]
    for i in range(1, n + 1):
        b[i] = -sum(c[k] * b[i - k] for k in range(1, min(i, len(c) - 1) + 1)) / c[0]
    return b


N2 = 2 * (L + 2)            # x-order for the two-variable series (q-order L + 2)


def theta_w(kind, n):
    """theta_kind(tau, z) in w = q^{1/8}, keys (w-exponent, 2 * y-exponent); theta_1 without -i"""
    s = defaultdict(Fr)
    for k in range(-60, 61):
        if kind in (1, 2):
            e, ky, c = (2 * k + 1) ** 2, 2 * k + 1, ((1 if k % 2 == 0 else -1) if kind == 1 else 1)
        else:
            e, ky, c = 4 * k * k, 2 * k, (1 if kind == 3 or k % 2 == 0 else -1)
        if e <= n:
            s[(e, ky)] += c
    return dict(s)


def ratio_w(num, den, n):
    """num/den for two-variable w-series, den depending on w only through its leading term"""
    d = at_z0(den)
    lead = min(i for i, v in enumerate(d) if v)
    dd = d[lead:] + [Fr(0)] * (n + 1)
    return {(e - lead, y): v for (e, y), v in m2(num, lift(inv1(dd[:n + 1], n)), n + lead).items()
            if e - lead <= n}


def at_two_tau(S, n):
    """a w-series with even w-exponents and even y keys, read at 2 tau, as an x-series (x = q^{1/2})"""
    o = {}
    for (e, ky), v in S.items():
        assert e % 2 == 0 and ky % 2 == 0, (e, ky)
        if e // 2 <= n:
            o[(e // 2, ky // 2)] = v                 # q^{e/8} at 2 tau = q^{e/4} = x^{e/2}
    return o


NW = 2 * N2 + 16                                      # w-order needed: x^{N2} at 2 tau is w^{2 N2}
A1 = {}
for kind in (2, 3, 4):
    th = theta_w(kind, NW)
    th0 = lift(at_z0(th))
    term = ratio_w(m2(th, th, NW), m2(th0, th0, NW), NW)
    for k, v in term.items():
        A1[k] = A1.get(k, Fr(0)) + v
A1 = {k: v for k, v in A1.items() if v}
th4 = theta_w(4, NW)
th40 = lift(at_z0(th4))
T1 = ratio_w(m2(th4, th4, NW), m2(th40, th40, NW), NW)
th1 = theta_w(1, NW)
eta6 = [Fr(0)] * 2 + pprod([(8, 6)], NW)[:NW - 1]     # eta(tau)^6 = w^2 prod (1 - w^{8r})^6
B1 = ratio_w({k: -v for k, v in m2(th1, th1, NW).items()}, lift(eta6), NW)
A2 = at_two_tau(A1, N2)                               # A(2 tau, z)
T = at_two_tau(T1, N2)                                # T = theta_4(2tau,z)^2/theta_4(2tau)^2
B2 = at_two_tau(B1, N2)                               # B(2 tau, z)


def qser(coeffs_in_q):
    """a series in q given as dict n -> c, as a one-variable x-series"""
    c = [Fr(0)] * (N2 + 1)
    for nq, v in coeffs_in_q.items():
        if 2 * nq <= N2:
            c[2 * nq] = Fr(v)
    return c


calE = [Fr(0)] * (N2 + 1)
calE[0] = Fr(1)
for r in range(1, N2 // 2 + 1):
    calE[2 * r] = Fr(24 * sig_odd(r))                      # E(tau) = 2E_2(2tau) - E_2(tau)
E4t = qser({r: (240 * sig(3, r) if r else 1) for r in range(0, L + 3)})
E4tt = qser({2 * r: (240 * sig(3, r) if r else 1) for r in range(0, L + 3)})
Mser = [(a - b) / 240 for a, b in zip(E4t, E4tt)]          # M = (E_4(tau) - E_4(2tau))/240
H2 = {}
for nq in range(0, L + 3):
    for j in range(-2 * math.isqrt(nq) - 2, 2 * math.isqrt(nq) + 3):
        h = hurwitz(4 * nq - j * j)
        if h and 4 * nq <= N2:
            H2[(4 * nq, j)] = h                             # H(2tau,z): q^{2n} = x^{4n}

# closed form (10): N = g psi_1^F
Nform = {}
for c, S in ((Fr(-20, 9), m2(A2, lift(calE), N2)),
             (Fr(-4, 45), m2(B2, lift(E4t), N2)),
             (Fr(-7, 15), m2(B2, lift(E4tt), N2)),
             (Fr(-104), H2)):
    for k, v in S.items():
        Nform[k] = Nform.get(k, Fr(0)) + c * v
Nform = {k: v for k, v in Nform.items() if v}

# ---------- C1: product formula minus polar part ----------
coef = SD.coefficient_tables('chl2', 2, L, 1, P, 'R')
g_x = [Fr(0)] + pprod([(2, 8), (4, 8)], N2)[:-1]          # g = x^2 prod (1-x^{2r})^8 (1-x^{4r})^8 ... shifted below
g_x = [Fr(0)] * 2 + pprod([(2, 8), (4, 8)], N2)[:N2 - 1]   # g = q prod(1-q^r)^8(1-q^{2r})^8 = x^2 (...)
ginv_shift = inv1(pprod([(2, 8), (4, 8)], N2), N2)         # 1/g = x^{-2} ginv_shift
AL = defaultdict(Fr)                                       # sum_s q^{2(s^2+s)} y^{2s+1}/(1-q^{2s} y)^2, chamber R
for s in range(-8, 9):
    for k in range(0, 80):
        if s == 0:
            AL[(0, -(k + 1))] += k + 1
            if k > P + 2:
                break
        elif s > 0:
            e = 2 * (s * s + s) + 2 * s * k
            if 2 * e > N2 + 2:
                break
            AL[(2 * e, 2 * s + 1 + k)] += k + 1
        else:
            e = 2 * (s * s + s) - 2 * s * (k + 2)
            if 2 * e > N2 + 2:
                break
            AL[(2 * e, 2 * s + 1 - (k + 2))] += k + 1
AL = {k: Fr(v) for k, v in AL.items() if v}
polar = m2(AL, lift(ginv_shift), N2)                        # g psi^P / 52 * (1/g) = x^{-2} * polar
closed = m2(Nform, lift(ginv_shift), N2)                    # psi^F = x^{-2} * closed
c1 = True
for Q2 in range(0, L - 1):
    for p in range(-P + 2, P - 1):
        raw = coef.get((Q2, 1, p), 0)
        lhs = raw - 52 * polar.get((2 * Q2 + 2, p), 0)
        if lhs != closed.get((2 * Q2 + 2, p), 0):
            c1 = False
ok('C1 product formula minus the polar part 52 AL/g equals the closed form, Q^2 <= 28, |Q.P| <= 7', c1)

# ---------- C2: theta identities ----------
lhs = A2
rhs = {}
for c, S in ((Fr(3), T), (Fr(-1, 4), m2(B2, lift(calE), N2))):
    for k, v in S.items():
        rhs[k] = rhs.get(k, Fr(0)) + c * v
rhs = {k: v for k, v in rhs.items() if v}
ok('C2a A(2tau,z) = 3 T - B(2tau,z) E(tau)/4', {k: v for k, v in lhs.items() if k[0] <= N2 - 4} ==
   {k: v for k, v in rhs.items() if k[0] <= N2 - 4})
EE = m1(calE, calE, N2)
ok('C2b E(tau)^2 = (E_4(tau) + 4 E_4(2tau))/5', all(EE[i] == (E4t[i] + 4 * E4tt[i]) / 5 for i in range(N2 + 1)))
Nform2 = {}
for c, S in ((Fr(-20, 3), m2(T, lift(calE), N2)), (Fr(16, 3), m2(B2, lift(Mser), N2)), (Fr(-104), H2)):
    for k, v in S.items():
        Nform2[k] = Nform2.get(k, Fr(0)) + c * v
Nform2 = {k: v for k, v in Nform2.items() if v}
ok('C2c the closed form equals -20/3 E T + 16/3 M B(2tau,z) - 104 H(2tau,z)',
   {k: v for k, v in Nform.items() if k[0] <= N2 - 4} == {k: v for k, v in Nform2.items() if k[0] <= N2 - 4})

# ---------- C3: G is a combination of E_4(x), E_4(x^2), E_4(x^4) ----------
R = [Fr(0)] + pprod([(4, 8), (2, -4)], K)[:K]               # eta(2tau)^8/eta(tau)^4 = x prod(...)
calEK = [Fr(0)] * (K + 1)
calEK[0] = Fr(1)
for r in range(1, K // 2 + 1):
    calEK[2 * r] = Fr(24 * sig_odd(r))
MK = [Fr(0)] * (K + 1)
for r in range(1, K // 2 + 1):
    MK[2 * r] = Fr(sig(3, r) - (sig(3, r // 2) if r % 2 == 0 else 0))
G = [Fr(20, 3) * a - Fr(16, 3) * b for a, b in zip(m1(calEK, R, K), MK)]
Gform = [a / 36 - 49 * b / 180 + 11 * c / 45 for a, b, c in zip(E4x(1, K), E4x(2, K), E4x(4, K))]
ok('C3a G = E4(x)/36 - 49 E4(x^2)/180 + 11 E4(x^4)/45 on x^0, x^1, x^2 (Sturm bound for M_4(Gamma_0(4)))',
   G[:3] == Gform[:3])
ok('C3b the same identity to order x^K', G == Gform)

# ---------- C4: S = -j_0 + j_1 = theta_3 G / eta(2tau)^6 ----------
comp = {}
for (e, yv), v in Nform2.items():
    D = e - yv * yv
    if abs(yv) <= 6 and e <= N2 - 2:
        if (D, yv % 2) in comp:
            assert comp[(D, yv % 2)] == v
        comp[(D, yv % 2)] = v
j = {}
for (D, par), v in comp.items():
    j[(D, par)] = v + 104 * hurwitz(D)                    # remove the H(2tau,z) component
th3 = [Fr(0)] * (K + 1)
for n in range(-12, 13):
    if n * n <= K:
        th3[n * n] += 1
cK = m1(th3, G, K)
S = m1(cK, pprod([(4, -6)], K), K)                         # theta_3 G / eta(2tau)^6 = x^{-1} S
c4 = all((-j.get((t, 0), 0) if t % 2 == 0 else j.get((t, 1), 0)) == S[t + 1]
         for t in range(0, N2 - 4))
ok('C4 -j_0 + j_1 = theta_3 G / eta(2tau)^6 on the theta components', c4)

# ---------- C5: the coefficients of G ----------
def gform(jj):
    if jj <= 0:
        return Fr(0)
    if jj % 2:
        return Fr(20, 3) * sig(3, jj)
    a, m = 0, jj
    while m % 2 == 0:
        a, m = a + 1, m // 2
    return -Fr(16, 3) * 8 ** (a - 1) * sig(3, m)


ok('C5 g(j) = 20/3 sigma_3(j) for odd j and -16/3 8^(a-1) sigma_3(m) for j = 2^a m', all(G[i] == gform(i) for i in range(K + 1)))

# ---------- C6: positivity of theta_3 G ----------
ok('C6a c(k) > 0 for 1 <= k <= K', all(cK[k] > 0 for k in range(1, K + 1)))
z3 = Fr(12021, 10000)                                    # zeta(3) < 1.2021
lo = (Fr(20, 3) - Fr(2, 3) * z3) * Fr(16, 35)
hi = Fr(20, 3) + Fr(2, 3) * z3
ok('C6b the bound (14) is positive for k >= 8: lo^2 * 8 > hi^2', lo * lo * 8 > hi * hi)

# ---------- C7: n_0 and n_1 ----------
ok('C7a H(N) <= N for 3 <= N <= 400', all(hurwitz(n) <= n for n in range(3, 401)))
ok('C7b H(N) = 0 for N = 1, 2', hurwitz(1) == 0 and hurwitz(2) == 0)
s_of = lambda k: 8 * sig(1, k) + (16 * sig(1, k // 2) if k % 2 == 0 else 0)
a_k = {k: -(j[(2 * k, 0)] - 104 * hurwitz(2 * k)) for k in range(1, 30) if (2 * k, 0) in j}
ok('C7c a_k = -n_0(2k) > 2 s(k)/k for k <= 29', len(a_k) == 29 and all(a > Fr(2 * s_of(k), k) for k, a in a_k.items()))
# for k >= 4: a_k >= c(2k+1) >= (2k+1)^3 (3 lo - hi) and 2 s(k)/k <= 32 (1 + log k) <= 32 (1 + k)
ok('C7d 3 lo - hi > 23/40 and 23/40 (2k+1)^3 > 32 (1 + k) for k >= 4',
   3 * lo - hi > Fr(23, 40) and all(Fr(23, 40) * (2 * k + 1) ** 3 > 32 * (1 + k) for k in range(4, 400)))
ok('C7e s(k)/k <= 16 (1 + log k) for k <= 2000',
   all(s_of(k) <= 16 * k * (1 + math.log(k)) for k in range(1, 2001)))
n1 = {t: j[(t, 1)] - 104 * hurwitz(t) for t in range(1, 60, 2) if (t, 1) in j}
ok('C7f n_1(t) > 0 for odd t <= 59', len(n1) == 30 and all(v > 0 for v in n1.values()))
# for k = t + 1 >= 11: c(k) >= k^3 (lo sqrt k - hi) > 104 k > 104 H(t)
ok('C7g 11 lo^2 > (hi + 104/121)^2, so k^3 (lo sqrt k - hi) > 104 k for k >= 11',
   11 * lo * lo > (hi + Fr(104, 121)) ** 2)

# ---------- C8: the signs of the single-centred indices ----------
c8 = True
count = 0
for Q2 in range(0, L - 1):
    for p in range(-P + 2, P - 1):
        if 2 * Q2 - p * p >= 0:
            v = closed.get((2 * Q2 + 2, p), 0)
            if (-1) ** abs(p + 1) * v <= 0:
                c8 = False
            count += 1
ok(f'C8 (-1)^(Q.P+1) c^F(Q^2, Q.P) > 0 for the {count} charges with P^2 = 2, 2 Q^2 >= (Q.P)^2, Q^2 <= 28', c8)

for name in checks:
    print('ok ', name)
print(f'{len(checks)} checks passed')
Path(__file__).with_name('chl-positivity-checks.json').write_text(json.dumps({'checks': checks}, indent=1) + '\n')
