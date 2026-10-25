# The Conjectures of Ashoke Sen

Source, data and verification code of the paper *The Conjectures of Ashoke Sen*
by Deep Bhattacharjee.

The paper goes through thirty-six statements that Ashoke Sen made as
conjectures or as results extended beyond their derivation, from the S-duality
of the heterotic string in 1994 to the soft graviton theorems and the
D-instanton normalisation, and gives each one a verdict: proved (twenty-two),
confirmed without a proof (ten), incomplete (two) or false as stated (two).
The two failures are the positivity of the single-centred dyon index, which
the paper refutes by an exact computation in the Z_2 toroidal orbifold of type
II string theory, and tachyon matter as the final state of an inhomogeneous
brane decay, which the paper refutes by a caustic theorem. For every claim
that is not proved the paper names the missing step and records that no
proof exists in the literature that cites it.

## Layout

| Path | Content |
| --- | --- |
| `main.pdf` | the compiled paper |
| `tex/` | the LaTeX source: `main.tex`, `sections/`, `tables/`, `figures/` |
| `code/siegel_dyons.py` | exact Fourier coefficients of 1/Phi_10 and of the Siegel form of the Z_2 orbifold, the wall-crossing check and the list of negative single-centred indices; writes `siegel_dyons_out.json` |
| `lean/SenConjectures.lean` | the same expansions carried out inside the Lean 4 kernel, with the theorems of the paper that are finite arithmetic proved by `decide` |
| `figures/make_figures.py` | the generator of the three figures |

## Reproducing

```
cd code && python3 siegel_dyons.py          # Python 3, standard library only
cd lean && lean SenConjectures.lean         # Lean 4.35.0, no Mathlib; prints the axioms of every theorem
cd tex && pdflatex main && pdflatex main && pdflatex main
```

The Lean file compiles against a bare toolchain (`lean-toolchain` pins it) and
needs no `lake`; `lake build` also works. Every theorem depends on no axiom or on
`propext` alone; none uses `sorry` or native evaluation. The verdicts of the
paper are statements about string theory and are not formalised; `lean/README.md`
says exactly which statements are.

## Releases

| Release | Zenodo DOI |
| --- | --- |
| concept DOI of all versions, resolving to the latest | [10.5281/zenodo.23652472](https://doi.org/10.5281/zenodo.23652472) |
| v1.0.0 | [10.5281/zenodo.23652473](https://doi.org/10.5281/zenodo.23652473) |
| v1.0.1 (the paper citing the archive by its concept DOI) | [10.5281/zenodo.23652641](https://doi.org/10.5281/zenodo.23652641) |

## Author and contact

Deep Bhattacharjee, ORCID
[0000-0003-0466-750X](https://orcid.org/0000-0003-0466-750X). Formerly,
Electro-Gravitational Space Propulsion Laboratory (EGSPL), Bhubaneswar, Odisha
751030, India; formerly, University of Calcutta (Heramba Chandra College),
23/49 Gariahat Road, Golpark, Kolkata 700029, West Bengal, India.
Emails: d.bhattacharjee@erl-forschung.de, itsdeep@live.com, egspl@proton.me.
For anything related to this repository write to egspl@proton.me.

## Licence

MIT, see `LICENSE`.
