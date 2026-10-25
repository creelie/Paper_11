# Lean 4 certificate

`SenConjectures.lean` carries out, inside the Lean 4 kernel, the finite
arithmetic on which the paper *The Conjectures of Ashoke Sen* rests, and proves
the corresponding statements by `decide +kernel`. It needs Lean 4.35.0
(`lean-toolchain`) and nothing else: no Mathlib, no `lake`, no native
evaluation.

```
lean SenConjectures.lean
```

prints one line per theorem, `'Sen.<name>' depends on axioms: [propext]` or
`'Sen.<name>' does not depend on any axioms`. The run takes about eleven
minutes on four cores and peaks at ten gigabytes of memory. The GitHub workflow `verify.yml` runs it on
every push and fails if a theorem depends on `sorryAx`, `Classical.choice` or
`Lean.ofReduceBool`, or if the number of theorems changes.

## What is proved

Everything is computed from first principles: the theta functions, the Dedekind
eta function and the Eisenstein series are defined as truncated series in
`x = q^{1/8}` and `y = e^{2 pi i z}`, the twisted elliptic genera of the
heterotic string and of the Z_2 toroidal orbifold are built from them, and the
Fourier coefficients of 1/Phi_10 and of the Siegel form of the orbifold are
obtained by expanding the Borcherds-type product, exactly as in
`code/siegel_dyons.py`.

| Theorem | Statement of the paper |
| --- | --- |
| `dh_values` | d_h(n) = 1, 24, 324, 3200, 25650 for Q^2 = -2, 0, 2, 4, 6 (coefficients of q eta^{-24}) |
| `het_table` | Table 1: the indices of 1/Phi_10 for Q^2, P^2 in {2,4,6}, Q.P = -2, ..., 4 |
| `wall_crossing` | the wall-crossing formula d_R - d_L = (-1)^{Q.P} Q.P d_h(Q) d_h(P) for Q^2, P^2 in {2,4,6}, abs(Q.P) <= 4 |
| `tor2_rows_P0`, `tor2_rows_P2`, `tor2_rows_P4` | Table 3: the indices of the Z_2 orbifold for Q^2 <= 12, P^2 in {0,2,4}, Q.P = 0, 1, 2 |
| `violations_in_table` | Table 2: the twenty-one listed charges have the listed negative indices |
| `violations_single_centred` | the twenty-one charges satisfy the single-centred conditions (4.7) of Chattopadhyaya and David |

The remaining theorems (`*_genera`, `*_step_*`, `*_upTo_*`, `*_finish`,
`*_table`, and the three unfolding lemmas) are the stages of the two
expansions. The expansion is too long for one kernel computation, so the file
stores the elliptic genera and the partial products after each power of
e^{2 pi i rho} as data blocks, encoded as natural numbers, and the kernel
checks that each block is what the definitions produce from the previous one;
the statements about the tables are then assembled from these steps by
`Eq.trans`. The blocks are certificates, not inputs: a wrong block makes its
step fail.

## What is not formalised

The verdicts of the paper are statements about string theory, and the kernel
checks none of them. In particular the following are outside Lean:

- the refutation of the positivity conjecture rests on the computed indices
  above together with the identification of the single-centred region and the
  argument of the paper on the hair of the black hole; only the arithmetic is
  checked;
- the caustic theorem for the dust of the rolling tachyon is an argument in
  ordinary differential equations and is proved by hand in the paper;
- the positivity spot checks of the heterotic string and of the 2A CHL model
  up to Q^2 = 12 and the Z_3 orbifold, which `code/siegel_dyons.py` performs,
  are not repeated in Lean;
- every verdict of the form "proved", "confirmed" or "incomplete" is a
  judgement about the literature and has no formal counterpart.
