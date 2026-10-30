#!/usr/bin/env python3
"""Independent check of the closed form on the full Laurent expansion:
   D(sigma/2) phi_1(sigma,v) = 52 Av[y/(1-y)^2] + 48 S(sigma/2) A - (20/3) Y - 104 calH,
with both sides expanded in chamber R (|x| small, |y| slightly above 1), phi_1 taken from
the product formula of 1/Phi~_6 and Av[y/(1-y)^2](sigma,v) = sum_s p^{s^2+s} y^{2s+1}/(1-p^s y)^2, p = x^2."""
import sys
from fractions import Fraction as Fr
from math import isqrt
import chl_unit_identities as I

K = int(sys.argv[1]) if len(sys.argv) > 1 else 20
W = 2 * K + 4
tab = I.SD_coeffs(K + 2, 3 * K + 12)
phi1 = {(k2, p): v for (l, k2, p), v in tab.items() if l == 1}
Dser = I.lift(I.qprod([(1, 8), (2, 8)], K + 2), 1)
Psi = {k: v for k, v in I.mul(Dser, phi1, K).items() if abs(k[1]) <= W}
# Appell-Lerch in chamber R
Av = {}
for n in range(1, W + 1): Av[(0, -n)] = n                       # s = 0: y/(1-y)^2 = sum n y^{-n}
for s in range(1, isqrt(K) + 2):
    k = 1
    while 2 * (s * s + s) + 2 * s * (k - 1) <= K:                 # s > 0
        Av[(2 * (s * s + s) + 2 * s * (k - 1), 2 * s + 1 + (k - 1))] = k; k += 1
    k = 1
    while 2 * (s * s + s) + 2 * s * (k - 1) <= K:                 # -s < 0
        Av[(2 * (s * s + s) + 2 * s * (k - 1), -2 * s - 1 - (k - 1))] = k; k += 1
A = I.phim21(K, 2)
S = {(n, 0): sum(d ** 3 for d in range(1, n + 1) if n % d == 0 and (n // d) % 2) for n in range(1, K + 1)}
H = I.hurwitz(2 * K + 4)
calH = {(e, l): H[2 * e - l * l] for e in range(K + 1) for l in range(-isqrt(2 * e), isqrt(2 * e) + 1)
        if 2 * e - l * l >= 0 and H[2 * e - l * l] != 0}
Y = {(e, l): I.r3(2 * e - l * l) for e in range(K + 1) for l in range(-isqrt(2 * e), isqrt(2 * e) + 1)
     if I.r3(2 * e - l * l)}
rhs = I.add((52, Av), (48, I.mul(S, A, K)), (Fr(-20, 3), Y), (-104, calH))
rhs = {k: v for k, v in rhs.items() if abs(k[1]) <= W and k[0] <= K}
diff = I.add((1, Psi), (-1, rhs))
diff = {k: v for k, v in diff.items() if abs(k[1]) <= W - 2}
print('closed form = D*phi_1 on chamber R, x^e with e<=%d, |l|<=%d:' % (K, W - 2), len(diff) == 0)
sys.exit(0 if not diff else 1)
