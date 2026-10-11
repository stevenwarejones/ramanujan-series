# Elliptic companions: a circle, a torus, and interacting frequencies

A circle can hide surprisingly organized arithmetic. This study follows one
hypergeometric family behind elliptic integrals and asks what happens when we
change its parameters. Its fourth logarithmic coefficient is

```math
c=[s^{4}]\log {}_2F_1(\tfrac12+s,\tfrac12+s;1;\tfrac12)
=-0.1399920769176353288990254098142611333152\ldots.
```

The notation means “the coefficient of the fourth power of s.” It measures a
response beyond the simplest quadratic one. It is not another name for π.

The geometric path lives on a square torus. In suitable coordinates its
logarithmic weight becomes a sum of odd cosine waves. Averaging around the
circle keeps only combinations of frequencies that add to zero. At quadratic
order opposite frequencies cancel. At fourth order, different frequencies can
interact: **1 + 1 + 1 − 3 = 0** is the first example.

This one interaction accounts for about **99.89% of the correction to the
single-wave approximation of c**. The percentage concerns that correction,
not the whole value of c.

| Consecutive odd harmonics retained | Absolute error in c, approximately |
|---|---:|
| 1 | 2.62 × 10⁻³ |
| 1, 3 | 2.85 × 10⁻⁶ |
| 1, 3, 5 | 3.51 × 10⁻⁹ |
| Through 9 | 6.65 × 10⁻¹⁵ |
| Through 19 | 5.89 × 10⁻²⁹ |

The [research note](NOTE.md) derives the path, exact coefficient rule, and
monotone approximations with rigorous error bounds. A certificate using 61
positive odd harmonics encloses c in an interval narrower than 10⁻¹⁵⁰.
The amplitudes already contain π; this is **not a new algorithm for π or a
convergence record for Ramanujan series**.

## The bounded follow-up and its stopping point

The [reduction audit](BOUNDED_REDUCTION.md) fixes an explicit alphabet of
fourth-order torsion points, checks the square-torus symmetry without dropping
path information, and tests a precisely stated 28-function family of simpler
expressions. An exact rational witness rules out that functional ansatz.
The one-part logarithmic alphabet is not preserved by the quarter-turn: a
second differential appears.

**No finite evaluation of c in familiar constants was obtained.** The functional
exclusion says nothing about an extra relation only at this special numerical
point. The work does not claim irrationality, a minimal basis of all elliptic
periods, or novelty established by a literature search. Further open-ended
search is paused at this checkpoint.

## Reproduce

From the repository root, with Python 3.11 or newer:

```sh
python research/elliptic-companions/verify_fourier.py
python research/elliptic-companions/verify_bounded_reduction.py
python -m unittest discover -s tests -p test_elliptic_companions.py -v
```

These checks use only the standard library. They write JSON to stdout and
leave tracked files unchanged. Do not run Python with `-O`, which disables
assertions. The repository's existing CI runs the regression tests on Python
3.11 and 3.12, including against the extracted research bundle.

Optional exact symbolic geometry and floating cross-checks:

```sh
python -m pip install -r research/elliptic-companions/requirements-audit.txt
python research/elliptic-companions/audit_geometry.py
python research/elliptic-companions/check_torus_path.py
```

The geometry audit uses exact SymPy algebra. The path checks use mpmath and
are explicitly **not** rigorous interval quadrature. Published outputs are in
[certificates](certificates/). The integer interval calculations are independent
of those optional packages and of stored decimal digits of special constants.

The mathematical framework is classical; the specific derivations, bounds,
and implementation here have not been independently refereed or formally
verified in a proof assistant. Sources and the scope of each claim appear in
the two notes. This directory is a focused companion study, separate from the
repository's exhaustive algebraic-degree frontier claims.
