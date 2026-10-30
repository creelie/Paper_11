"""Z_N CHL, N = 3,5,7, N Q^2/2 = 1: the numerator calN(D), the finite coefficients c^F and their signs."""
import sys, os
from fractions import Fraction as Fr
from math import isqrt, log
import chl_unit_identities as I
from chl_odd_identities import DATA, EN, E4, cusp4, s3N, sigma, divisors

def series_div(a, b, n):
    """a/b to order n (lists of Fractions, b[0] != 0)"""
    out = [Fr(0)] * (n + 1)
    for k in range(n + 1):
        s = Fr(a[k]) if k < len(a) else Fr(0)
        for i in range(1, k + 1):
            if i < len(b) and b[i]: s -= b[i] * out[k - i]
        out[k] = s / b[0]
    return out

def disc_tables(Dmax):
    """alpha(Delta), C01(Delta) for -1 <= Delta <= Dmax: coefficients of phi_{-2,1}, phi_{0,1} by discriminant."""
    T = Dmax // 4 + 3                    # q-degree
    U = 4 * T + 4                        # units of q^{1/4}
    th0 = [0] * (U + 1)                  # sum q^{n^2}            (y^0 of theta_3(z)^2)  in units 1/4
    th1 = [0] * (U + 1)                  # sum q^{n^2-n+1/2}      (y^1 of theta_3(z)^2)
    t2a = [0] * (U + 1)                  # sum q^{(n+1/2)^2}      (y^0 of theta_2(z)^2)
    t2b = [0] * (U + 1)                  # sum q^{n^2+1/4}        (y^1 of theta_2(z)^2)
    R = isqrt(U) + 3
    for n in range(-R, R + 1):
        for arr, e in ((th0, 4 * n * n), (th1, 4 * n * n - 4 * n + 2), (t2a, (2 * n + 1) ** 2), (t2b, 4 * n * n + 1)):
            if 0 <= e <= U: arr[e] += 1
    th3sq = [0] * (U + 1); th4sq = [0] * (U + 1); th2sq = [0] * (U + 1)
    for n in range(-R, R + 1):
        for m in range(-R, R + 1):
            e = 2 * (n * n + m * m)
            if e <= U:
                th3sq[e] += 1; th4sq[e] += (-1) ** ((n + m) % 2)
            e2 = (2 * n + 1) ** 2 + (2 * m + 1) ** 2      # units 1/8 -> /2 for units 1/4
            if e2 // 2 <= U and e2 % 2 == 0: th2sq[e2 // 2] += 1
    # theta_2^2 starts at q^{1/4} (unit 1); shift both numerator and denominator
    sh = lambda a, k: a[k:] + [0] * k
    r0 = [Fr(0)] * (U + 1); r1 = [Fr(0)] * (U + 1)
    for num0, num1, den, s in ((th0, th1, th3sq, 0), (th0, [-v for v in th1], th4sq, 0), (t2a, t2b, th2sq, 1)):
        d = sh(den, s); a0 = sh(num0, s); a1 = sh(num1, s)
        q0 = series_div(a0, d, U - 2); q1 = series_div(a1, d, U - 2)
        for i in range(U - 1): r0[i] += 4 * q0[i]; r1[i] += 4 * q1[i]
    # r0: y^0 coefficients, units q^{1/4}, only multiples of 4 nonzero ; r1: y^1, likewise
    C01 = {}
    for t in range(T):
        C01[4 * t] = r0[4 * t]; C01[4 * t - 1] = r1[4 * t]
    al = I_alpha(Dmax + 8)
    return al, C01

def I_alpha(Dm):
    import chl_unit_positivity as P
    return P.alpha_table(Dm)

class Num:
    def __init__(self, N, Dmax):
        self.N = N; D = DATA[N]; self.D = D
        self.Dmax = Dmax
        tmax = Dmax // 4 + 2
        self.al, self.C01 = disc_tables(Dmax // N + 8)
        e4 = E4(tmax); e4N = E4(tmax, N); c4 = cusp4(N, tmax); en = EN(N, tmax)
        u, w, z = D['F4']
        self.f = [D['a'] * e4N[t] + u * e4[t] + w * e4N[t] + z * c4[t] for t in range(tmax + 1)]
        self.en = en
        uW, wW, zW = D['W']
        self.mu = [self.f[t] - D['c'] * (uW * e4[t] + wW * e4N[t] + zW * c4[t]) for t in range(tmax + 1)]
        self.H = I.hurwitz(Dmax // N + 8)
    def A(self, d): return self.al.get(d, 0) if d >= -1 else 0
    def B(self, d): return self.C01.get(d, 0) if d >= -1 else 0
    def Hd(self, D):
        if D % self.N: return Fr(0)
        m = D // self.N
        return self.H[m] if 0 <= m < len(self.H) else Fr(0)
    def calN(self, D):
        N, a, c = self.N, self.D['a'], self.D['c']
        s = Fr(0); t = 0
        while 4 * t <= D + N:
            r = D - 4 * t
            if r % N == 0:
                s += self.f[t] * self.A(r // N) - c * self.en[t] * self.B(r // N)
            t += 1
        return s - 288 * a * self.Hd(D)
    def calN_mu(self, D, Theta):
        """same quantity via mu A - 12 c Theta_L - 288 a calH, Theta = Theta_N(D)"""
        N, a, c = self.N, self.D['a'], self.D['c']
        s = Fr(0); t = 0
        while 4 * t <= D + N:
            r = D - 4 * t
            if r % N == 0: s += self.mu[t] * self.A(r // N)
            t += 1
        return s - 12 * c * Theta - 288 * a * self.Hd(D)

def ginv(N, n):
    """coefficients b_j of x/g(sigma/N) = prod (1-x^m)^{-(k+2)} (1-x^{Nm})^{-(k+2)}"""
    k = 24 // (N + 1)
    g = I.qprod([(1, k), (N, k)], n)
    b = [1] + [0] * n
    for m in range(1, n + 1): b[m] = -sum(g[i] * b[m - i] for i in range(1, m + 1))
    return b

def cF(Nm, b, Delta):
    return sum(b[j] * Nm(Delta + 4 - 4 * j) for j in range(0, (Delta + 4 + Nm.__self__.N) // 4 + 1))
