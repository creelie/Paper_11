import Mathlib.Analysis.Complex.Basic
import Mathlib.Tactic.Abel
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-! Elementary bounds used in Appendix F.  These theorems do not formalise
string amplitudes, pole factorisation, phase-space integrals or infrared
finiteness.  All complex amplitudes here are arbitrary complex numbers. -/
namespace SenInfrared

theorem squared_norm_error (x y : ℂ) :
    |‖x + y‖ ^ 2 - ‖x‖ ^ 2| ≤ 2 * ‖x‖ * ‖y‖ + ‖y‖ ^ 2 := by
  have hu := norm_add_le x y
  have hd := abs_norm_sub_norm_le (x + y) x
  have hc : x + y - x = y := by abel
  rw [hc] at hd
  have hf : ‖x + y‖ ^ 2 - ‖x‖ ^ 2 =
      (‖x + y‖ - ‖x‖) * (‖x + y‖ + ‖x‖) := by ring
  rw [hf, abs_mul, abs_of_nonneg (add_nonneg (norm_nonneg _) (norm_nonneg _))]
  calc
    _ ≤ ‖y‖ * (‖x + y‖ + ‖x‖) :=
      mul_le_mul_of_nonneg_right hd (add_nonneg (norm_nonneg _) (norm_nonneg _))
    _ ≤ ‖y‖ * (2 * ‖x‖ + ‖y‖) :=
      mul_le_mul_of_nonneg_left (by linarith) (norm_nonneg _)
    _ = _ := by ring

theorem weighted_difference_bound (x y : ℂ) (f f0 B D : ℝ)
    (hB : |f| ≤ B) (hD : |f - f0| ≤ D) :
    |f * ‖x + y‖ ^ 2 - f0 * ‖x‖ ^ 2| ≤
      D * ‖x‖ ^ 2 + B * (2 * ‖x‖ * ‖y‖ + ‖y‖ ^ 2) := by
  have hb : 0 ≤ B := le_trans (abs_nonneg _) hB
  have he : f * ‖x + y‖ ^ 2 - f0 * ‖x‖ ^ 2 =
      (f - f0) * ‖x‖ ^ 2 + f * (‖x + y‖ ^ 2 - ‖x‖ ^ 2) := by ring
  rw [he]
  calc
    _ ≤ |(f - f0) * ‖x‖ ^ 2| + |f * (‖x + y‖ ^ 2 - ‖x‖ ^ 2)| :=
      abs_add_le _ _
    _ = |f - f0| * ‖x‖ ^ 2 + |f| * |‖x + y‖ ^ 2 - ‖x‖ ^ 2| := by
      rw [abs_mul, abs_mul, abs_of_nonneg (sq_nonneg ‖x‖)]
    _ ≤ D * ‖x‖ ^ 2 + B * (2 * ‖x‖ * ‖y‖ + ‖y‖ ^ 2) :=
      add_le_add (mul_le_mul_of_nonneg_right hD (sq_nonneg _))
        (mul_le_mul hB (squared_norm_error x y) (abs_nonneg _) hb)

theorem corner_jacobian (r t : ℝ) : t * (-r) - r * (1 - t) = -r := by
  ring

theorem blown_up_corner_density (r t : ℝ) (hr : r ≠ 0) :
    r / (r * t + r * (1 - t)) ^ 2 = 1 / r := by
  have hh : r * t + r * (1 - t) = r := by ring
  rw [hh]
  field_simp

theorem soft_degree_identity (l i v e d : ℤ)
    (hl : l = i - v) (he : e = 2 * i - d) :
    4 * l - 2 * i - e + 2 * v = d - 2 * v := by
  linarith

theorem positive_soft_degree (l i v e d : ℤ)
    (hl : l = i - v) (he : e = 2 * i - d)
    (hv : 0 < v) (hd : 3 * v ≤ d) :
    0 < 4 * l - 2 * i - e + 2 * v := by
  linarith

end SenInfrared

#print axioms SenInfrared.squared_norm_error
#print axioms SenInfrared.weighted_difference_bound
#print axioms SenInfrared.corner_jacobian
#print axioms SenInfrared.blown_up_corner_density
#print axioms SenInfrared.soft_degree_identity
#print axioms SenInfrared.positive_soft_degree
