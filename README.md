# Proving Conjectures of Ashoke Sen: F-Theory Limits, Small Black Holes and Dyon Indices

Deep Bhattacharjee. The paper examines thirty-six statements of Ashoke Sen on duality, dyon indices, tachyon condensation, F-theory and black hole entropy, and decides which are proved, which are supported but unproved, and where each stops. Two of them are proved in general form, and two new cases of a third: in the `Z_2` CHL model the single-centred dyon index is positive for every charge with `P^2 = 2` (Theorem 5.2) and for every charge with `N Q^2/2 = 1` (Theorem 5.7, claim 9). Every weak coupling degeneration of F-theory on an elliptic Calabi-Yau variety is, after base change and rescaling, Sen's family, and one exists exactly when `-2K_B` has a section (claim 23). The entropy of a two-charge small black hole is fixed by one critical value for every curvature-corrected action; the Gauss-Bonnet truncation matches the statistical entropy only in four and five dimensions, and Prester's Lovelock completion is the unique one that matches it in every dimension (claim 26). For the complete heterotic action in `D >= 6`, which is not a claim of Sen, the paper proves that the central-charge route does not give the statistical value under its own hypotheses (Theorem 10.10, in a highlighted box); the statement itself remains open. Three highlights mark what Sen got wrong (rose), what is fatal (burnt orange) and what is left unproved (yellow); no statement of Sen examined is wrong as printed or fatal, and a table lists what is left unproved, claim by claim. Exact product expansions reproduce the negative full dyon indices of the toroidal orbifolds; with the charge cone of Sen's equation (3.9) transcribed correctly, the scan over `Q^2 <= 12`, `P^2 <= 8`, `Q.P <= 4` contains thirty negative-index charge triples. The dust equations of tachyon matter are given their characteristic flow, the Hessian sign criterion and the attainment hypothesis under which the first focusing time is a minimum. Two appendices construct a finite-energy soft-emission model and a subtracted tree-emission bound; the full four-dimensional string-field-theory infrared problem is not proved. Every result the paper uses is taken from refereed publications.

## Layout

- `tex/`: the LaTeX source (`main.tex`, `sections/`, `tables/`, `figures/`) and the compiled `main.pdf`.
- `dist/`: the submission files, `ConjecturesOfAshokeSen_final2.pdf`, `ConjecturesOfAshokeSen-tex.zip` and `ConjecturesOfAshokeSen-arxiv.tar.gz`.
- `code/siegel_dyons.py`: exact integer arithmetic for the Fourier coefficients of `1/Phi_10` and of the Siegel form of the `Z_2` toroidal orbifold; it writes `siegel_dyons_out.json`, which the tables of the paper are read from. `verify_revision.py` rechecks the thirty-triple scan, `verify_infrared_estimates.py` the explicit constants of Appendices E and F, `verify_closures.py` (sympy, mpmath) the identities and numbers of Appendix G for claims 23 and 26, and `verify_chl_positivity.py` (standard library) the closed form, the theta identities, the Eisenstein identity and the values used in the proof of Theorem 5.2; their outputs are the four `*-checks.json` files. `chl_unit_identities.py`, `chl_unit_positivity.py` and `chl_unit_fullcheck.py` (standard library) check the Jacobi form identities, the closed form and the sign lemmas of Theorem 5.7.
- `lean/`: the Lean 4 certificates. `SenConjectures.lean` (Lean 4.34.1, no library) recomputes the coefficient tables inside the kernel from the definitions of the theta functions and checks the wall-crossing identities and the twenty-one negative triples with `P^2 <= 4`, `Q.P <= 2`; `AnalyticScopeChecks.lean` and `InfraredBounds.lean` check elementary identities of the analytic sections against Mathlib. See `lean/README.md` for what is and is not formalised.
- `figures/make_figures.py`: the figure generator.

## Checks

```
cd code && python3 siegel_dyons.py /tmp/out.json && python3 verify_revision.py && python3 verify_infrared_estimates.py && python3 verify_closures.py && python3 verify_chl_positivity.py && python3 chl_unit_identities.py 48 && python3 chl_unit_positivity.py && python3 chl_unit_fullcheck.py
cd lean && lean SenConjectures.lean
```

The Lean run prints the axioms of each of the fifty-seven theorems; every one is axiom-free or depends on `propext` only. The `verify` workflow runs both on every push.

## Releases

| Release | Zenodo DOI |
| --- | --- |
| v1.0.0 (25 October 2026) | [10.5281/zenodo.23652473](https://doi.org/10.5281/zenodo.23652473) |
| v1.0.1 (26 October 2026) | [10.5281/zenodo.23652641](https://doi.org/10.5281/zenodo.23652641) |
| v1.1.0 (26 October 2026) | [10.5281/zenodo.23656435](https://doi.org/10.5281/zenodo.23656435) |
| v2.0.0 (26 October 2026) | [10.5281/zenodo.23669248](https://doi.org/10.5281/zenodo.23669248) |
| v2.0.1 (26 October 2026) | [10.5281/zenodo.23670581](https://doi.org/10.5281/zenodo.23670581) |
| v2.1.0 (28 October 2026) | [10.5281/zenodo.23746042](https://doi.org/10.5281/zenodo.23746042) |
| v2.2.0 (28 October 2026) | [10.5281/zenodo.23747112](https://doi.org/10.5281/zenodo.23747112) |
| v2.2.1 (28 October 2026) | [10.5281/zenodo.23749469](https://doi.org/10.5281/zenodo.23749469) |
| v3.0.0 (30 October 2026) | [10.5281/zenodo.23789251](https://doi.org/10.5281/zenodo.23789251) |
| v3.1.0 (30 October 2026, this revision) | [10.5281/zenodo.23790654](https://doi.org/10.5281/zenodo.23790654) |

The concept DOI of all versions is [10.5281/zenodo.23652472](https://doi.org/10.5281/zenodo.23652472).

## Contact

Deep Bhattacharjee, d.bhattacharjee@erl-forschung.de, itsdeep@live.com, egspl@proton.me; formerly EGSPL, Bhubaneswar, and formerly University of Calcutta (Heramba Chandra College), Golpark, Kolkata; GitHub [creelie](https://github.com/creelie). For anything related to this repository (GitHub-related contact) write to egspl@proton.me.

## Licence

MIT, see `LICENSE`.
