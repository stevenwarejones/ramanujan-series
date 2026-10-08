"""Independent pi checks and a negative control for the certified summation."""
from pathlib import Path
import json
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))

from flint import arb, ctx, fmpz_poly
from benchmark import coefficients, certified_sum
from identity32 import ROOT, isolated_root, refined_root


class IdentityTests(unittest.TestCase):
    def setUp(self):
        self.original_dps = ctx.dps
        ctx.dps = 640

    def tearDown(self):
        ctx.dps = self.original_dps

    def test_classical_identities_enclose_pi(self):
        for method in ('ramanujan-1103', 'chudnovsky'):
            with self.subTest(method=method):
                A, B, x, s = coefficients(method, 640)
                _, error, approximation, enclosure = certified_sum(A, B, x, s, 500)
                reference = arb.pi()
                self.assertTrue(enclosure.contains(reference))
                self.assertLess(error, arb(10)**(-500))
                self.assertLess(abs(approximation-reference), arb(10)**(-500))

    def test_degree32_algebraic_construction_encloses_pi(self):
        A, B, x, s = coefficients('cm-degree-32', 640)
        terms, error, approximation, enclosure = certified_sum(A, B, x, s, 500)
        reference = arb.pi()
        self.assertEqual(terms, 2)
        self.assertTrue(enclosure.contains(reference))
        self.assertLess(error, arb(10)**(-500))
        self.assertLess(abs(approximation-reference), arb(10)**(-500))

    def test_perturbed_coefficient_fails_independent_pi_check(self):
        A, B, x, s = coefficients('cm-degree-32', 640)
        # A small error, invisible to a low-precision spot check, must be detected.
        A += arb(10)**(-400)
        _, _, approximation, enclosure = certified_sum(A, B, x, s, 500)
        reference = arb.pi()
        self.assertFalse(enclosure.contains(reference))
        self.assertGreater(abs(approximation-reference), arb(10)**(-410))

    def test_refined_root_contains_independent_high_precision_isolation(self):
        poly=fmpz_poly(json.loads((ROOT/'results/candidate32.json').read_text())['p'])
        # Two targets exercise several precision-doubling steps. The reference
        # uses FLINT's separate all-roots algorithm at still higher precision.
        for dps in (640, 2128):
            with self.subTest(dps=dps):
                ctx.dps=dps
                target=ctx.prec
                root=refined_root(poly)
                self.assertEqual(ctx.prec,target)
                self.assertGreaterEqual(root.rel_accuracy_bits(),target-12)
                ctx.dps=dps+100
                self.assertTrue(root.contains(isolated_root(poly)))


if __name__ == '__main__':
    unittest.main()
