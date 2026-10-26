# Formal verification scope

The core certificate and the two new supplements use separate toolchains. A successful check is not a formal proof of the full manuscript or four-dimensional string-field-theory infrared finiteness.

- SenConjectures.lean: bare Lean 4.34.1; 57 named statements in the axiom audit. It expands truncated products and checks the original 21 coefficient triples using the corrected linear cone predicate. It is not an exhaustive check of the 30-entry Python scan. Only propext, or no axioms, occur in its printed audit. No Mathlib is needed for this file.
- AnalyticScopeChecks.lean: Lean 4.34.0-rc2 with Mathlib commit 85e3a25e006c35636f0e53b0e9296caca2685bc0; 12 elementary sign, polynomial, index-logic and exponential identities.
- InfraredBounds.lean: the same Mathlib environment; six complex-norm, corner-algebra and homogeneous-degree identities. These are not a formalisation of an amplitude or a phase-space integral.

There are 75 named audited statements. The Mathlib supplements use propext, Classical.choice and Quot.sound. No admitted proof or custom axiom appears in the successful audits. Formal receipts record hashes and completed runs.

The bare-core command is `lean SenConjectures.lean` under the pinned lean-toolchain. It does not build the Mathlib supplements. The supplements require the separately pinned Mathlib build and its dependencies. Local run scripts in the verification directory record the exact invocation.

Not formalised: dust PDE and global geometric arguments, probability construction, tree amplitude factorisation, phase-space integrals, uniform loop estimates, flat-function analysis, physical interpretations and literature verdicts. The full string-field-theory infrared problem remains open in this manuscript.
