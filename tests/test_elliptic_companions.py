"""Independent representations, explicit tails, and exact exclusion witnesses."""
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import unittest

RESEARCH = Path(__file__).resolve().parents[1]/'research'/'elliptic-companions'
sys.path.insert(0, str(RESEARCH))
from intervals import I, SCALE, pi_interval, exp_positive, gauss_log_jet
from verify_bounded_reduction import build
from verify_fourier import finite_moments


class EllipticCompanionTests(unittest.TestCase):
    def test_certified_fourier_against_gauss_and_published_snapshot(self):
        data = json.loads(subprocess.check_output(
            [sys.executable, str(RESEARCH/'verify_fourier.py')], text=True))
        stored = json.loads((RESEARCH/'certificates'/'fourier.json').read_text())
        self.assertEqual(data, stored)
        self.assertTrue(data['frequency_balance_certificate']['independent_gauss_intersection'])
        self.assertEqual(data['frequency_balance_certificate']['positive_odd_harmonics'], 61)

    def test_exact_exclusion_witness_and_negative_control(self):
        data = json.loads((RESEARCH/'certificates'/'bounded_reduction.json').read_text())
        labels, columns, target = build(data['maximum_q_degree'])
        self.assertEqual(labels, data['basis'])
        w = [F(data['witness'].get(str(n),'0')) for n in range(len(target))]
        dot = lambda v: sum((a*b for a,b in zip(w,v)), F(0))
        self.assertTrue(all(dot(col)==0 for col in columns))
        self.assertEqual(dot(target), 1)
        # One changed target coefficient must invalidate this exact witness.
        index = next(i for i,a in enumerate(w) if a)
        target[index] += 1
        self.assertNotEqual(dot(target), 1)

    def test_modular_primitive_matches_gauss_with_tail(self):
        # Independent q recurrence; |C_n| <= 35 n^3, proof in BOUNDED_REDUCTION.md.
        N = 150
        _, _, target = build(N-1)
        q = exp_positive(pi_interval()).reciprocal()
        total = I.exact(0)
        for value in reversed(target):
            total = total*q + I.exact(value.numerator)/value.denominator
        tail = 35*N**3*q**N/(1-q*I.exact(N+1)**3/N**3)
        result = I(total.lo-tail.hi, total.hi)
        reference = gauss_log_jet()[4]
        self.assertLessEqual(max(result.lo,reference.lo),min(result.hi,reference.hi))
        self.assertLess(result.hi-result.lo,SCALE//10**150)

    def test_independent_phase_model_misses_the_resonance(self):
        # Use exact rational amplitudes, unrelated to the saved numerical target.
        a, b = I.exact(2), I.exact(1)/4
        moments = finite_moments({1:a,3:b},4)
        joint = (moments[4]-3*moments[2]**2)/24
        independent = -(a**4+b**4)/4
        difference = joint-independent
        resonance = a**3*b/3
        self.assertGreater(difference.lo,0)
        self.assertLessEqual(max(difference.lo,resonance.lo),min(difference.hi,resonance.hi))


if __name__ == '__main__':
    unittest.main()
