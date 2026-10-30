"""phi_1 in chamber R from the product formula: phi_1 = Psi / (g(sigma/N) phi_{-2,1}(sigma,v))."""
import sys, os
from fractions import Fraction as Fr
from chl_odd_identities import *

def chamberR(N, X, W):
    M = Model(N, X)
    Psi = M.Psi()                       # = g(sigma/N) phi_{-2,1} phi_1, series in x
    # 1/g(sigma/N) = x^{-1} prod (1-x^n)^{-(k+2)} (1-x^{Nn})^{-(k+2)}
    k2 = 24 // (N + 1)
    g0 = eta_series([(1, k2), (N, k2)], X + 2)
    b = [Fr(1)] + [0] * (X + 2)
    for n in range(1, X + 3): b[n] = -sum(g0[k] * b[n - k] for k in range(1, n + 1))
    # phi_{-2,1}(sigma,v) = (y-2+1/y) Pi ; 1/Pi as series in x
    Pi = {}
    phim = M.in_x(M.phim)
    # Pi = phim / (y-2+1/y): divide each x-coefficient Laurent poly exactly by (y-1)^2/y
    # easier: invert phim/(y-2+1/y) by Newton on x-series: Pi = 1 + O(x)
    from collections import defaultdict
    byx = defaultdict(dict)
    for (n, j), v in phim.items(): byx[n][j] = v
    def divp(p):   # exact division of Laurent poly p by y-2+1/y
        if not p: return {}
        lo = min(p); hi = max(p); q = {}; r = dict(p)
        for e in range(hi, lo + 1, -1):
            c = r.get(e, 0)
            if c:
                q[e - 1] = c
                r[e] = r.get(e, 0) - c; r[e - 1] = r.get(e - 1, 0) + 2 * c; r[e - 2] = r.get(e - 2, 0) - c
        assert all(v == 0 for v in r.values()), "not divisible"
        return {k: v for k, v in q.items() if v}
    PiX = {n: divp(byx[n]) for n in range(X + 1)}
    # inverse of PiX (PiX[0] == {0:1})
    assert PiX[0] == {0: 1}
    inv = {0: {0: Fr(1)}}
    for n in range(1, X + 1):
        acc = defaultdict(Fr)
        for k in range(1, n + 1):
            for e1, c1 in PiX[k].items():
                for e2, c2 in inv[n - k].items():
                    acc[e1 + e2] -= c1 * c2
        inv[n] = {e: v for e, v in acc.items() if v}
    # 1/(y-2+1/y) = sum_{m>=1} m y^{-m}  (|y| > 1)
    out = defaultdict(Fr)
    PsiX = defaultdict(dict)
    for (n, j), v in Psi.items(): PsiX[n][j] = v
    # phi_1 coefficient at x^e: sum over x-degrees  (1/g: x^{-1} b_i) * Psi(n1) * inv(n2), n1+n2+i-1 = e
    for e in range(-1, X - 1):
        acc = defaultdict(Fr)
        for i in range(0, e + 2):
            if not b[i]: continue
            for n1 in range(0, e + 2 - i):
                n2 = e + 1 - i - n1
                if n2 < 0 or n2 > X: continue
                for j1, c1 in PsiX.get(n1, {}).items():
                    for j2, c2 in inv[n2].items():
                        acc[j1 + j2] += b[i] * c1 * c2
        for jj, c in acc.items():
            for m in range(1, W + 1):
                out[(e, jj - m)] += m * c
    return {k: v for k, v in out.items() if v}

if __name__ == '__main__':
    N = int(sys.argv[1]); X = int(sys.argv[2])
    T = chamberR(N, X, 40)
    for P2 in (2, 4, 6):
        e = N * P2 // 2
        print((2, N), P2, [T.get((e, l), 0) for l in (-2, 0, 1, 2, 3, 4)])
    for e in range(-1, 7):
        print('e', e, [str(T.get((e, l), 0)) for l in (-2, 0, 1, 2, 3, 4)])
