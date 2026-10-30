#!/usr/bin/env python3
"""Positivity of the single-centred index of the Z_2 CHL model at N Q^2/2 = 1.

With Delta' = 2e - l^2 (x = e^{pi i sigma}),
  D(sigma/2) phi_1^F(sigma,v) = 48 S(sigma/2) phi_{-2,1}(sigma,v) - (20/3) Y(sigma,v) - 104 calH(sigma,v),
so its coefficient at x^e y^l is
  N(Delta') = 48 sum_{n>=1} s3odd(n) alpha(Delta'-2n) - (20/3) r_3(Delta') - 104 H(Delta'),
  s3odd(n) = sum_{d|n, n/d odd} d^3,  phi_{-2,1} = sum alpha(4a-r^2) q^a y^r.
The coefficient of x^n y^l of phi_1^F is
  Cf(Delta) = sum_{j>=0} g(j) N(Delta+2-2j),  Delta = 2n - l^2,  1/D(tau) = sum_j g(j) q^{j-1},
and the single-centred index is d* = (-1)^{l+1} Cf(Delta).

The script (i) evaluates N and Cf exactly, (ii) checks them against the product
formula (chamber R) for 0 <= l <= 2, (iii) checks the hypotheses used in the proof
(sign pattern of alpha, N(odd) > 0 for Delta' in {1,3,5,7}, the bound for Delta' >= 9,
-N(2t) >= 96 s3odd(t)), and (iv) prints d* for small charges.
"""
import sys
from fractions import Fraction as Fr
from math import isqrt
import chl_unit_identities as I

DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000

def s3odd(n):
    return sum(d ** 3 for d in range(1, n + 1) if n % d == 0 and (n // d) % 2 == 1)

def eta_inv_pow(k, nmax):
    """coefficients of prod (1-q^n)^{-k}"""
    return I.qprod([(1, -k)], nmax)

def alpha_table(Dmax):
    """alpha(Delta) for -1 <= Delta <= Dmax from phi_{-2,1} = (theta_0 vt_{1,1} - theta_1 vt_{1,0})/eta^6"""
    nmax = Dmax // 4 + 3
    p6 = eta_inv_pow(6, nmax)
    al = {}
    # h_0 = -theta_1/eta^6 = -sum_{r odd} q^{(r^2-1)/4} * sum p6(n) q^n   (integer powers, Delta = 4*power)
    # h_1 = theta_0/eta^6 = q^{-1/4} sum_{r even} q^{r^2/4} sum p6(n) q^n  (Delta = 4*power - 1)
    h0 = [0] * (nmax + 1); h1 = [0] * (nmax + 1)
    for r in range(-2 * isqrt(4 * nmax) - 3, 2 * isqrt(4 * nmax) + 4):
        if r % 2:
            e = (r * r - 1) // 4
            for n in range(0, nmax + 1 - e): h0[e + n] -= p6[n]
        else:
            e = r * r // 4
            for n in range(0, nmax + 1 - e): h1[e + n] += p6[n] if e + n <= nmax else 0
    for t in range(nmax + 1):
        if 4 * t <= Dmax: al[4 * t] = h0[t]
        if 4 * t - 1 <= Dmax: al[4 * t - 1] = h1[t]
    return al

def main():
    Dm = DMAX
    out = {}
    al = alpha_table(Dm + 4)
    A = lambda D: al.get(D, 0) if D >= -1 else 0
    # check alpha against the bivariate phi_{-2,1}
    A2 = I.phim21(30, 1)
    out['alpha_matches_phi_-2,1'] = all(v == A(4 * e - p * p) for (e, p), v in A2.items())
    out['alpha_sign_pattern'] = all((A(D) > 0 if D % 4 == 3 else A(D) < 0 if D % 4 == 0 else A(D) == 0) for D in range(-1, Dm + 1))
    H = I.hurwitz(Dm + 4)
    Hf = lambda D: H[D] if 0 <= D <= Dm + 4 else Fr(0)
    r3c = [I.r3(n) for n in range(Dm + 5)]
    r3f = lambda D: r3c[D] if 0 <= D <= Dm + 4 else 0
    s3 = [0] + [s3odd(n) for n in range(1, Dm // 2 + 4)]
    def N(D):
        if D < -1: pass
        tot = Fr(0)
        n = 1
        while D - 2 * n >= -1:
            tot += 48 * s3[n] * A(D - 2 * n)
            n += 1
        return tot - Fr(20, 3) * r3f(D) - 104 * Hf(D)
    Nt = {D: N(D) for D in range(-6, Dm + 3)}
    out['N_integral'] = all(v.denominator == 1 for v in Nt.values())
    out['N_negative_args_vanish'] = all(Nt[D] == 0 for D in range(-6, 0))
    out['N(0)=2'] = (Nt[0] == 2)
    out['N(odd)>0 for 1<=D<=%d' % Dm] = all(Nt[D] > 0 for D in range(1, Dm + 1, 2))
    out['N(1),N(3),N(5),N(7)'] = [int(Nt[D]) for D in (1, 3, 5, 7)]
    out['-N(2t) >= 96 s3odd(t) for 1<=t<=%d' % (Dm // 2)] = all(-Nt[2 * t] >= 96 * s3[t] for t in range(1, Dm // 2 + 1))
    # analytic bound used for odd D >= 9: 6(D+1)^3 > 369 D, and the ingredient bounds
    out['bound_6(D+1)^3>369D_for_D>=9 (checked to Dm)'] = all(6 * (D + 1) ** 3 > 369 * D for D in range(9, Dm + 1, 2))
    out['r3(D)<=24D'] = all(r3f(D) <= 24 * D for D in range(1, Dm + 1))
    out['H(D)<=D/3+2sqrt(D/3)'] = all(float(Hf(D)) <= D / 3 + 2 * (D / 3) ** 0.5 + 1e-12 for D in range(1, Dm + 1))
    # N against the bivariate closed form 48 S A - 20/3 Y - 104 calH (x-grading), small range
    # Cf(Delta)
    jmax = Dm // 2 + 3
    g = I.qprod([(1, 8), (2, 8)], jmax)
    # invert prod(1-x^n)^8(1-x^2n)^8
    g = I.qprod([(1, -8), (2, -8)], jmax)
    def Cf(D):
        return sum(g[j] * Nt[D + 2 - 2 * j] for j in range(0, (D + 2 + 6) // 2 + 1) if D + 2 - 2 * j >= -6)
    Ct = {D: Cf(D) for D in range(-2, Dm - 4)}
    out['Cf integral'] = all(Fr(v).denominator == 1 for v in Ct.values())
    out['(-1)^(D+1) Cf(D) > 0 for -1<=D<=%d' % (Dm - 5)] = all(((-1) ** (D + 1)) * Ct[D] > 0 for D in range(-1, Dm - 4))
    # compare with product formula, chamber R, 0<=l<=2
    K = 24
    tab = I.SD_coeffs(K, 3 * K + 12)
    mism = []
    for k2 in range(-1, K + 1):
        for l in range(0, 3):
            D = 2 * k2 - l * l
            if D < -2: continue
            if tab.get((1, k2, l), 0) != Ct[D]: mism.append((k2, l, tab.get((1, k2, l), 0), Ct[D]))
    out['chamberR(1,k2,l)=Cf(2k2-l^2), 0<=l<=2, k2<=%d' % K] = (len(mism) == 0)
    if mism: print(mism[:5])
    for k, v in out.items(): print(k, v)
    print('d* (N Q^2/2 = 1) = (-1)^{l+1} Cf(Delta), Delta = Q^2 P^2 - (Q.P)^2 = 2n - l^2:')
    print(' '.join('%d:%d' % (D, (-1) ** (D + 1) * Ct[D]) for D in range(-1, 16)))
    good = all(v for k, v in out.items() if isinstance(v, bool))
    return good, Ct, Nt

if __name__ == '__main__':
    good, Ct, Nt = main()
    print('POSITIVITY HYPOTHESES VERIFIED' if good else 'FAILURE')
    sys.exit(0 if good else 1)
    sys.exit(0 if good else 1)
