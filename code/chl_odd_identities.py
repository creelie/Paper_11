"""Z_N CHL models, N = 3, 5, 7, first electric charge N Q^2/2 = 1: the identities and inequalities."""
import sys, time
from fractions import Fraction as Fr
from math import isqrt
from chl_odd_series import S2
import os
JEXTRA = int(os.environ.get("JEXTRA", "8"))

def divisors(n): return [d for d in range(1, n + 1) if n % d == 0]
def sigma(n, k=1): return sum(d ** k for d in divisors(n))

def weak_forms(H, J):
    """phi_{-2,1}, phi_{0,1} in (q, y), q-degree <= H, |j| <= J (integers)."""
    X = H
    one = S2(X, J, {(0, 0): 1})
    p = S2(X, J, {(0, 1): 1, (0, 0): -2, (0, -1): 1})
    prod = one
    for n in range(1, X + 1):
        f1 = S2(X, 1, {(0, 0): 1, (n, 1): -1}); f2 = S2(X, 1, {(0, 0): 1, (n, -1): -1})
        prod = prod * f1; prod = prod * f1; prod = prod * f2; prod = prod * f2
        prod = S2(X, J, {k: v for k, v in prod.d.items() if abs(k[1]) <= J}, prod.den)
    # 1/prod(1-q^n)^4 : q-only series
    inv = [0] * (X + 1); inv[0] = 1
    for n in range(1, X + 1):
        for _ in range(4):
            for k in range(n, X + 1): inv[k] += inv[k - n]
    invS = S2(X, 0, {(n, 0): inv[n] for n in range(X + 1) if inv[n]})
    Pi = prod * invS                  # prod (1-q^n y)^2 (1-q^n/y)^2 / (1-q^n)^4
    Pi = S2(X, J, {k: v for k, v in Pi.d.items() if abs(k[1]) <= J}, Pi.den)
    phim = p * Pi
    phim = S2(X, J, {k: v for k, v in phim.d.items() if abs(k[1]) <= J}, phim.den)
    # phi01 = phi_{-2,1} (1 + 12 sum_n q^n sum_{d|n} d (y^d - 2 + y^-d)) + 12 Pi
    P = {(0, 0): 1}
    for n in range(1, X + 1):
        for d in divisors(n):
            if d <= J:
                P[(n, d)] = P.get((n, d), 0) + 12 * d
                P[(n, -d)] = P.get((n, -d), 0) + 12 * d
            P[(n, 0)] = P.get((n, 0), 0) - 24 * d
    phi01 = phim * S2(X, J, P) + Pi.scale(12)
    phi01 = S2(X, J, {k: v for k, v in phi01.d.items() if abs(k[1]) <= J}, phi01.den)
    return phim, phi01

def EN(N, T):
    c = Fr(24, N - 1)
    return [Fr(1)] + [c * (sigma(n) - (N * sigma(n // N) if n % N == 0 else 0)) for n in range(1, T + 1)]

class Model:
    def __init__(self, N, X):
        self.N, self.X = N, X
        self.H = (2 * X) // N + 3
        self.J = isqrt(8 * X // N + 4) + 2 + JEXTRA
        self.phim, self.phi01 = weak_forms(self.H, self.J)
        self.e = EN(N, N * self.H)
    def a(self, n, j): return self.phi01.get(n, j) / 4
    def b(self, n, j): return -self.phim.get(n, j)
    def BE(self, n, j):      # coefficient of q^n y^j in B(tau) E_N(tau)
        return sum(self.e[t] * self.b(n - t, j) for t in range(0, n + 1))
    def BEN(self, nN, j):    # coefficient of q^{nN/N} y^j in B(tau) E_N(tau/N)
        N = self.N; s = Fr(0); t = nN % N
        while t <= nN:
            s += self.e[t] * self.b((nN - t) // N, j); t += N
        return s
    def expo(self, k2, l, j):
        """sum_s e^{-2 pi i s l/N} c^{(r,s)}(4 k' l - j^2), k' = k2/N, r = k2 mod N."""
        N = self.N; r = k2 % N; nN = k2 * l; dl = (l % N == 0)
        if r == 0:
            n = nN // N
            F00 = Fr(8, N) * self.a(n, j)
            F01 = Fr(8, N * (N + 1)) * self.a(n, j) - Fr(2, N + 1) * self.BE(n, j)
            tot = F00 + ((N if dl else 0) - 1) * F01
        else:
            tot = Fr(2, N + 1) * self.BEN(nN, j)
            if dl: tot += Fr(8, N + 1) * self.a(nN // N, j)
        assert tot.denominator == 1, (N, k2, l, j, tot)
        return int(tot)
    def Psi(self):
        X, J = self.X, self.J
        L1 = {}; L2 = {}
        for k2 in range(0, X + 1):
            for l in (1, 2):
                jm = isqrt(4 * k2 * l // self.N + 1) + 1
                for j in range(-jm, jm + 1):
                    c = self.expo(k2, l, j)
                    if not c: continue
                    if l == 1:
                        L1[(k2, j)] = c
                        if 2 * k2 <= X: L2[(2 * k2, 2 * j)] = L2.get((2 * k2, 2 * j), 0) + Fr(c, 2)
                    else:
                        L2[(k2, j)] = L2.get((k2, j), 0) + c
        S1 = S2.from_dict(X, J, L1)
        return S2.from_dict(X, J, L2) + (S1 * S1).scale(Fr(1, 2))
    def in_x(self, s, Xmax=None):
        """phi(sigma, v) given in q = e(sigma) -> series in x = e(sigma/N)."""
        X = Xmax or self.X
        return S2(X, s.J, {(self.N * n, j): v for (n, j), v in s.d.items() if self.N * n <= X}, s.den)

def s3N(N, t): return sum(d ** 3 for d in divisors(t) if (t // d) % N)

def eta_series(exps, X):
    """prod_{(d,e)} prod_n (1 - q^{dn})^e as a list (q-prefactor omitted)."""
    s = [0] * (X + 1); s[0] = 1
    for d, e in exps:
        for n in range(1, X // d + 1):
            for _ in range(e):
                for k in range(X, d * n - 1, -1): s[k] -= s[k - d * n]
    return s

def cusp4(N, X):
    """the normalised cusp form of weight 4 for Gamma_0(N), N = 5, 7, as a list in q."""
    if N == 5:
        s = eta_series([(1, 4), (5, 4)], X); return [0] + s[:X]
    if N == 7:
        s = eta_series([(1, 3), (7, 3)], X); s = [0] + s[:X]
        th = [0] * (X + 1)
        R = isqrt(4 * X) + 2
        for u in range(-R, R + 1):
            for w in range(-R, R + 1):
                m = u * u + u * w + 2 * w * w
                if m <= X: th[m] += 1
        return [sum(s[i] * th[n - i] for i in range(n + 1)) for n in range(X + 1)]
    return [0] * (X + 1)

def E4(X, scale=1):
    return [Fr(1) if n == 0 else (240 * sigma(n // scale, 3) if n % scale == 0 else 0) for n in range(X + 1)]

DATA = {
    # a, c_N, (F4 coefficients of E4(sigma/N), E4(sigma), C4(sigma/N)), (W coefficients), kappa
    3: dict(a=Fr(3, 16), c=Fr(3, 8), F4=(Fr(3, 160), Fr(27, 160), 0), W=(Fr(-1, 20), Fr(21, 20), 0), kappa=Fr(9)),
    5: dict(a=Fr(7, 72), c=Fr(7, 36), F4=(Fr(11, 3744), Fr(353, 3744), Fr(6, 13)), W=(Fr(-1, 104), Fr(105, 104), Fr(108, 13)), kappa=Fr(15, 13)),
    7: dict(a=Fr(1, 16), c=Fr(1, 8), F4=(Fr(1, 1200), Fr(37, 600), Fr(3, 10)), W=(Fr(-1, 300), Fr(301, 300), Fr(24, 5)), kappa=Fr(3, 10)),
}

def form_x(N, coeffs, X):
    """u E4(sigma/N) + w E4(sigma) + z C4(sigma/N) as a q-only S2 in x."""
    u, w, z = coeffs
    A = E4(X); B = E4(X, N); C = cusp4(N, X)
    return S2.from_dict(X, 0, {(n, 0): u * A[n] + w * B[n] + z * C[n] for n in range(X + 1)})

def lattice_theta(N, X):
    """Theta_L(sigma, v) = sum_{lambda in L} x^{Q(lambda)} y^{B(lambda, pibar)/N}"""
    G = {3: [[2,-1,0,0],[-1,2,0,0],[0,0,2,-1],[0,0,-1,2]],
         5: [[2,-1,-1,-1],[-1,2,0,0],[-1,0,4,-1],[-1,0,-1,4]],
         7: [[2,0,-1,0],[0,2,0,-1],[-1,0,4,0],[0,-1,0,4]]}[N]
    pb = {3: (-2, -1, 0, 0), 5: (-2, -1, -2, -1), 7: (-1, 0, -2, 0)}[N]
    import itertools
    # Q(v) >= lam_min |v|^2 / 2 ; lam_min >= 1/2 for these Gram matrices (checked below)
    import numpy as np
    lmin = min(np.linalg.eigvalsh(np.array(G, dtype=float)))
    R = int((2 * X / lmin) ** 0.5) + 1
    Gp = [sum(G[i][k] * pb[k] for k in range(4)) for i in range(4)]
    assert all(g % N == 0 for g in Gp) and sum(pb[i] * Gp[i] for i in range(4)) == 2 * N
    T = {}
    for v in itertools.product(range(-R, R + 1), repeat=4):
        q2 = sum(v[i] * G[i][k] * v[k] for i in range(4) for k in range(4))
        if q2 <= 2 * X:
            f = sum(v[i] * Gp[i] for i in range(4)) // N
            T[(q2 // 2, f)] = T.get((q2 // 2, f), 0) + 1
    return T

def check_identities(N, X):
    t0 = time.time()
    M = Model(N, X); D = DATA[N]
    Psi = M.Psi()
    pm = M.in_x(M.phim); p01 = M.in_x(M.phi01)
    ENx = S2.from_dict(X, 0, {(n, 0): M.e[n] for n in range(X + 1)})
    F4 = form_x(N, D['F4'], X)
    rhs = (p01 * p01).scale(D['a']) - (ENx * pm * p01).scale(D['c']) + F4 * pm * pm
    G1 = Psi - rhs
    ok1 = all(n > X for (n, j) in G1.d) or G1.iszero()
    # second identity: E_N(sigma/N) phi01 = W phi_{-2,1} + 12 Theta_L
    X2 = min(X, 40)
    W = form_x(N, D['W'], X2)
    TL = lattice_theta(N, X2)
    TLs = S2.from_dict(X2, M.J, {k: 12 * v for k, v in TL.items()})
    G3 = S2(X2, M.J, {k: v for k, v in (ENx * p01).d.items() if k[0] <= X2}, (ENx * p01).den) - W * S2(X2, M.J, {k: v for k, v in pm.d.items() if k[0] <= X2}, pm.den) - TLs
    ok3 = G3.iszero()
    return ok1, ok3, time.time() - t0, M, Psi

if __name__ == '__main__':
    N = int(sys.argv[1]); X = int(sys.argv[2])
    ok1, ok3, dt, M, Psi = check_identities(N, X)
    print(f"N={N}: Psi identity through x^{X}: {ok1};  theta identity through x^{min(X,40)}: {ok3};  {dt:.1f}s")
