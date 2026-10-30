#!/usr/bin/env python3
"""Positivity of the single-centred index of the Z_N CHL models, N = 3, 5, 7, at N Q^2/2 = 1.

Notation of the subsection ssec:chlunitodd: x = e(sigma/N), p = x^N, y = e(v),
Delta = 4n - N l^2 for the coefficient x^n y^l of phi_1, and
  calN_N(D) = sum_t mu_N(t) alpha((D-4t)/N) - 12 c_N r_N(D) - 288 a_N H(D/N),
  c^F(Delta) = sum_j beta(j) calN_N(Delta + 4 - 4j),   d* = (-1)^(l+1) c^F(Delta).

Checks (exact rational arithmetic; python-flint for the polynomial products):
 O1  the identity eq:chlunioddG1 through x^14, x^62, x^170 and eq:chlunioddG3 through x^40;
 O2  mu_N(t) = kappa_N (sigma_3^(N)(t) - c_N(t)) through x^40 (Step 4);
 O3  the Gram matrices have determinant N^2 and level N, Q(pibar) = N, <L, pibar> = N Z;
 O4  c^F(Delta) equals the coefficient of the product formula in chamber R for 0 <= l <= 2,
     where the polar part 144 a_N calB / g contributes nothing;
 O5  the bound of Lemma lem:chlunioddsigns(ii) is positive for e >= 8, 17, 25 and the exact
     values below it have the stated signs; mu_N(t) >= 0, with mu_N(t) > 0 for t >= 2;
 O6  d* > 0 for the listed charges, and the first values printed in the paper.
"""
import sys
from math import sqrt
from fractions import Fraction as Fr
import numpy as np
from chl_odd_identities import check_identities, cusp4, s3N, DATA
from chl_odd_chamber import chamberR
from chl_odd_positivity import Num, ginv

STURM = {3: 14, 5: 62, 7: 170}
EBOUND = {3: 8, 5: 17, 7: 25}
FIRST = {3: [27, 540, 378, 3294, 2646, 16200, 14094, 68418],
         5: [20, 100, 125, 460, 588, 1720, 2191, 5544],
         7: [3, 36, 18, 138, 81, 444, 273, 1242]}
GRAM = {3: [[2, -1, 0, 0], [-1, 2, 0, 0], [0, 0, 2, -1], [0, 0, -1, 2]],
        5: [[2, -1, -1, -1], [-1, 2, 0, 0], [-1, 0, 4, -1], [-1, 0, -1, 4]],
        7: [[2, 0, -1, 0], [0, 2, 0, -1], [-1, 0, 4, 0], [0, -1, 0, 4]]}
PIBAR = {3: (-2, -1, 0, 0), 5: (-2, -1, -2, -1), 7: (-1, 0, -2, 0)}
fails = []

def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok: fails.append(name)

for N in (3, 5, 7):
    print(f"--- N = {N}")
    ok1, ok3, dt, M, Psi = check_identities(N, STURM[N])
    check(f"O1 eq:chlunioddG1 through x^{STURM[N]}", ok1)
    check("O1 eq:chlunioddG3 through x^40", ok3)
    P = Num(N, 4 * 40 + 8)
    c4 = cusp4(N, 40)
    check("O2 mu_N = kappa_N (sigma_3^(N) - c_N) through x^40",
          all(P.mu[t] == DATA[N]['kappa'] * (s3N(N, t) - c4[t]) for t in range(1, 41)) and P.mu[0] == 0)
    G = np.array(GRAM[N], dtype=float); Gi = np.linalg.inv(G); pb = np.array(PIBAR[N])
    level = next(l for l in range(1, 100) if np.allclose(l * Gi, np.round(l * Gi))
                 and all(round(l * Gi[i][i]) % 2 == 0 for i in range(4)))
    Gp = [int(v) for v in np.array(GRAM[N]) @ pb]
    check("O3 det = N^2, level N, Q(pibar) = N, <L,pibar> = N Z",
          round(np.linalg.det(G)) == N * N and level == N and int(pb @ np.array(GRAM[N]) @ pb) == 2 * N
          and all(g % N == 0 for g in Gp))
    X = {3: 24, 5: 30, 7: 42}[N]
    T = chamberR(N, X, 40)
    b = ginv(N, X + 10)
    cache = {}
    def calN(D):
        if D not in cache: cache[D] = P.calN(D)
        return cache[D]
    def cF(Dl): return sum(b[j] * calN(Dl + 4 - 4 * j) for j in range(0, (Dl + 4 + N) // 4 + 1))
    check(f"O4 c^F = product formula for n <= {X - 3}, 0 <= l <= 2",
          all(cF(4 * e - N * l * l) == T.get((e, l), 0)
              for e in range(0, X - 2) for l in (0, 1, 2) if 4 * e - N * l * l >= -N))
    D = DATA[N]; a, c, k = float(D['a']), float(D['c']), float(D['kappa'])
    def bound(e):
        lead = 9 * e ** 3 if N == 3 else k * (e ** 3 - 2 * e ** 2)
        return lead - 288 * c / (N - 1) * e ** 2 - 288 * a * (4 / (3 * N) + 2 * sqrt(4 / (3 * N))) * e
    e0 = EBOUND[N]
    check(f"O5 bound positive for e >= {e0} (to 10^4), not at e = {e0 - 1}",
          all(bound(e) > 0 for e in range(e0, 10 ** 4)) and bound(e0 - 1) <= 0)
    check("O5 mu_N(t) >= 0, > 0 for t >= 2 (t <= 40)",
          all(P.mu[t] >= 0 for t in range(1, 41)) and all(P.mu[t] > 0 for t in range(2, 41)))
    check("O5 calN = 0 for D <= 0; < 0 at D = 4e; > 0 at D = 4e - N >= 5 (e < bound)",
          all(calN(Dd) == 0 for Dd in range(-12, 1))
          and all(calN(4 * e) < 0 for e in range(1, e0 + 1))
          and all(calN(4 * e - N) > 0 for e in range(1, e0 + 1) if 4 * e - N >= 5)
          and all(calN(4 * e - N) >= 0 for e in range(1, e0 + 1)))
    vals = []
    for Dl in range(1, 17):
        if Dl % 4 == 0: vals.append(-cF(Dl))
        elif Dl % 4 == (-N) % 4: vals.append(cF(Dl))
    check("O6 first values of d*: " + ", ".join(str(v) for v in vals), vals == FIRST[N])
    check("O6 d* > 0 for 1 <= Delta <= 120",
          all((cF(Dl) if Dl % 4 else -cF(Dl)) > 0 for Dl in range(1, 121) if Dl % 4 in (0, (-N) % 4)))

print(f"{'all checks passed' if not fails else str(len(fails)) + ' failed'}")
sys.exit(1 if fails else 0)
