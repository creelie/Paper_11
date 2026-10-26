"""Exact finite algebra checks for Appendix F, not an amplitude calculation.

The general moment identities are proved analytically in the manuscript.
This script checks their polynomial recurrence through degree 40, as a
transcription control. No numerical check is treated as an infrared theorem.
"""
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
import json


def moment_polynomial(m, rate):
    coefficients = [Fraction(0) for _ in range(m+1)]
    for j in range(m+1):
        coefficients[m-j] = Fraction(comb(m, j)*factorial(j), rate**j)
    return coefficients


def verify():
    checked = []
    for m in range(41):
        for rate in (1, 2):
            coefficients = moment_polynomial(m, rate)
            for degree in range(m+1):
                derivative_coefficient = ((degree+1)*coefficients[degree+1]
                                          if degree < m else Fraction(0))
                target = Fraction(rate) if degree == m else Fraction(0)
                assert rate*coefficients[degree]-derivative_coefficient == target
            checked.append({'degree': m, 'rate': rate})
    # Pure-graviton homogeneous degree: non-exhaustive arithmetic control.
    for vertices in range(9):
        for extra_valence in range(9):
            valence_sum = 3*vertices+extra_valence
            for internal_lines in range((valence_sum+1)//2, valence_sum+5):
                attachments = 2*internal_lines-valence_sum
                loops = internal_lines-vertices
                degree = 4*loops-2*internal_lines-attachments+2*vertices
                assert degree == valence_sum-2*vertices
                if vertices:
                    assert degree >= vertices > 0
    report = {
        'status': 'passed',
        'polynomial_recurrences_checked': len(checked),
        'maximum_degree': 40,
        'exact_arithmetic': 'fractions.Fraction and integer factorials',
        'scope': 'Finite transcription checks of moment-polynomial and homogeneous-degree identities.',
        'not_established': ['String amplitude factorisation', 'Uniform loop-remainder bound',
                            'Phase-space integration theorem in Lean',
                            'Full four-dimensional infrared finiteness']
    }
    target = Path(__file__).with_name('infrared-estimate-checks.json')
    target.write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    verify()
