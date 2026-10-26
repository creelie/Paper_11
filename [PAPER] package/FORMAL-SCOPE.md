# Exactly what is checked

All three current Lean source modules compile successfully. Lean's own kernel replay also succeeds for each compiled module, relative to its imports. The receipt lists 75 named audited statements, their source hashes and axioms. Replay is not a separate proof assistant or a verification of the physical assumptions.

| File | Toolchain | Named statements | Scope |
| --- | --- | ---: | --- |
| SenConjectures.lean | Lean 4.34.1, bare core | 57 | Finite truncated-series calculations, the original 21 stored negative triples, and their corrected linear cone. |
| AnalyticScopeChecks.lean | Lean 4.34.0-rc2, Mathlib 85e3a25e006c35636f0e53b0e9296caca2685bc0 | 12 | Elementary sign, index logic, discriminant, logarithmic coefficient and scalar exponential identities. |
| InfraredBounds.lean | Same pinned Mathlib environment | 6 | Squared complex-norm bound, weighted bound, corner Jacobian/density algebra, homogeneous degree identity/positivity. |

The core module's printed axioms are empty or propext. The Mathlib supplements use propext, Classical.choice and Quot.sound. The completed audits contain no sorryAx, native_decide or custom theorem assumption. Standard imported library dependencies remain part of the trust base.

The exact Python scan checks 30 negative full-index triples with a larger-truncation control. Only the original 21 are covered by the finite coefficient certificate. Finite moment-polynomial recurrence checks through degree 40 are transcription controls.

Not formalised: the entire dust PDE theorem, global geometry, probability construction, string amplitude factorisation, phase-space integration and logarithmic moments, uniform loop remainders, flat-function analysis or the catalogue's physical interpretations. There is no Lean theorem proving full four-dimensional infrared finiteness in this package. No such theorem is claimed in the manuscript.

See verification/all-lean-audit.json and the per-module compile/replay logs. For the core file run `lean SenConjectures.lean` under lean/lean-toolchain. Build the supplements in the separately pinned Mathlib environment. The local scripts record the actual machine paths used; those paths must be adapted on another computer.
