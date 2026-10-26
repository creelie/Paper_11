import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Analysis.SpecialFunctions.Exp

/-! Algebra and elementary sign logic only. This module does not formalise
the dust PDE, global geometry, string spectra, or the literature catalogue. -/
namespace SenScope

theorem compression_sign (h d : ℝ) (hd : 0 < d) : -h / d < 0 ↔ 0 < h := by
  constructor
  · intro H
    rcases div_neg_iff.mp H with H | H
    · exact False.elim (lt_asymm hd H.2)
    · exact neg_neg_iff_pos.mp H.1
  · intro H
    exact div_neg_iff.mpr (Or.inr ⟨neg_neg_iff_pos.mpr H, hd⟩)

theorem no_compression_for_nonpositive_hessian (h d : ℝ) (hd : 0 < d)
    (hh : h ≤ 0) : 0 ≤ -h / d :=
  div_nonneg (neg_nonneg.mpr hh) (le_of_lt hd)

theorem focal_time_positive (l : ℝ) (hl : l < 0) : 0 < -1 / l :=
  div_pos_of_neg_of_neg (by norm_num) hl

theorem focal_jacobian_zero (l : ℝ) (hl : l < 0) : 1 + (-1 / l) * l = 0 := by
  have hn : l ≠ 0 := ne_of_lt hl
  field_simp
  norm_num

theorem discriminant_expansion (h e a c : ℚ) :
    4 * (-3*h^2+e*a)^3 + 27 * (-2*h^3+e*h*a-e^2*c/12)^2 =
    -9*e^2*h^2*(a^2-h*c) + e^3*(4*a^3-(9/2:ℚ)*h*a*c) + (3/16:ℚ)*e^4*c^2 := by
  ring

theorem cover_pullback_factorisation (a x : ℚ) : a^2-x^2=(a-x)*(a+x) := by ring

theorem base_transverse_derivative_nonzero : (-1 : ℚ) ≠ 0 := by norm_num

theorem positive_index_does_not_exclude_fermions :
    (3 : ℤ)-2 > 0 ∧ (2 : ℤ)>0 := by norm_num

theorem negative_index_rejects_conjunction (A B C : Prop) (i : ℤ)
    (h : A ∧ B ∧ C → 0 < i) (hi : i < 0) : ¬ (A ∧ B ∧ C) := by
  intro habc
  exact (lt_asymm (h habc) hi)

theorem cone_linear_example : (3 : ℤ)*2 ≤ 2*3+2 ∧ ¬ ((3 : ℤ)*2^2 ≤ 2*3+2) := by norm_num

theorem leading_soft_exponent_cancel (L E U b : ℝ) :
    Real.exp (-b*(L-U)) * Real.exp (b*(E-U)) = Real.exp (-b*(L-E)) := by
  rw [← Real.exp_add]
  congr 1
  ring

theorem cutoff_log_coefficients (a n ns : ℚ) :
    a+3+4-n = a+7-n ∧ -3-2-2+3+ns = ns-4 := by constructor <;> ring

end SenScope

#print axioms SenScope.compression_sign
#print axioms SenScope.no_compression_for_nonpositive_hessian
#print axioms SenScope.focal_time_positive
#print axioms SenScope.focal_jacobian_zero
#print axioms SenScope.discriminant_expansion
#print axioms SenScope.cover_pullback_factorisation
#print axioms SenScope.base_transverse_derivative_nonzero
#print axioms SenScope.positive_index_does_not_exclude_fermions
#print axioms SenScope.negative_index_rejects_conjunction
#print axioms SenScope.cone_linear_example

#print axioms SenScope.leading_soft_exponent_cancel
#print axioms SenScope.cutoff_log_coefficients
