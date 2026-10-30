"""Two-variable series sum c(n,j) x^n y^j, n = 0..X, |j| <= J, exact integers or rationals
(stored as integer numerators over a common denominator), products by Kronecker substitution."""
from fractions import Fraction as Fr
from math import lcm
try:
    import flint
    HAVE_FLINT = True
except ImportError:
    HAVE_FLINT = False

class S2:
    def __init__(self, X, J, d=None, den=1):
        self.X, self.J = X, J
        self.d = d if d is not None else {}
        self.den = den
    @staticmethod
    def from_dict(X, J, dd):
        den = 1
        for v in dd.values():
            if isinstance(v, Fr): den = lcm(den, v.denominator)
        d = {}
        for (n, j), v in dd.items():
            if n <= X and abs(j) <= J and v:
                w = Fr(v) * den
                assert w.denominator == 1
                d[(n, j)] = int(w)
        return S2(X, J, d, den)
    def get(self, n, j):
        return Fr(self.d.get((n, j), 0), self.den)
    def items(self):
        for k, v in self.d.items(): yield k, Fr(v, self.den)
    def __mul__(self, o):
        X = min(self.X, o.X); J = self.J + o.J
        M = 2 * J + 1
        def enc(s):
            coeffs = {}
            for (n, j), v in s.d.items():
                if n <= X: coeffs[n * M + j + s.J * (1 if s is self else 0)] = v
            return coeffs
        # encode with offsets: self offset J1, other offset J2; product offset J1+J2=J
        def poly(s, off):
            c = [0] * ((X + 1) * M)
            for (n, j), v in s.d.items():
                if n <= X: c[n * M + j + off] += v
            return c
        a = poly(self, self.J); b = poly(o, o.J)
        if HAVE_FLINT:
            pr = (flint.fmpz_poly(a) * flint.fmpz_poly(b)).coeffs()
            pr = [int(c) for c in pr]
        else:
            raise RuntimeError("python-flint is required")
        out = {}
        for i, v in enumerate(pr):
            if v:
                n, r = divmod(i, M)
                if n <= X:
                    j = r - J
                    out[(n, j)] = out.get((n, j), 0) + v
        # note: carry-free since n*M+r with 0<=r<M and product offsets in [0, 2J]
        return S2(X, J, {k: v for k, v in out.items() if v}, self.den * o.den)
    def scale(self, c):
        c = Fr(c)
        return S2(self.X, self.J, {k: v * c.numerator for k, v in self.d.items()}, self.den * c.denominator)
    def __add__(self, o):
        den = lcm(self.den, o.den)
        d = {}
        for k, v in self.d.items(): d[k] = d.get(k, 0) + v * (den // self.den)
        for k, v in o.d.items(): d[k] = d.get(k, 0) + v * (den // o.den)
        return S2(min(self.X, o.X), max(self.J, o.J), {k: v for k, v in d.items() if v and k[0] <= min(self.X, o.X)}, den)
    def __sub__(self, o): return self + o.scale(-1)
    def iszero(self): return not self.d
