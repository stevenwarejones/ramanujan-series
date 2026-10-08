# Ramanujan-Type Series for 1/π

This project began with a circle: imagine admiring one on a museum wall, or
tracing your way around a circular room. How does something so simple lead to
Ramanujan's extraordinary formulas for π, full of factorials, square roots, and
seemingly mysterious whole numbers? We started by exploring those series with
that sense of curiosity about the circle's beauty and the arithmetic behind it.

The exploration led to a concrete question we could investigate with mathematics
and code:

**How much convergence can simple algebraic constants buy?**

Ramanujan found that π, the number connecting a circle's circumference to its
diameter, can be recovered from a remarkable sum:

```math
\frac{1}{\pi}
=\frac{2\sqrt{2}}{9801}
\sum_{n=0}^{\infty}
\frac{(4n)!\,(1103+26390n)}{(n!)^{4}\,396^{4n}}.
```

Start at n = 0, then add the terms for n = 1, 2, and so on; the exclamation mark
means factorial, a product of successive positive integers, with 0! = 1.
Multiply the sum by the factor in front to approach 1/π, then take the reciprocal
to recover π. Each additional term contributes roughly
eight decimal digits of accuracy. This is Ramanujan's classic identity, not a
formula discovered by this project; see the [literature and sources](REPORT.md#sources).

This repository asks: **can similar formulas gain many more digits per term
without making their constants arbitrarily complicated?** It studies four
classical families of these identities, explains some of their arithmetic
structure, and provides computations that readers can reproduce.

For the geometric motivation—from circles and torus paths to higher-dimensional
questions—read [Geometric Intuition and Open Questions](EXPLORATIONS.md). It
separates intuition, established mathematics, repository results, and open directions.

## What “complicated” means here

We measure **algebraic degree**, one way to describe the arithmetic complexity
of a number. Ordinary fractions have degree 1. The square root of 2 has degree
2 because its simplest polynomial equation is x² − 2 = 0. When several constants
appear in a formula, we count the number system they need **together**: the
square roots of 2 and 3 together require degree 4, even though each has degree 2.

Degree does not count how many digits a constant needs to describe it, or how
long a computer takes to calculate it. We therefore report convergence and
actual evaluation time separately.

## What the calculations establish

Within the **four precisely specified families**, exhaustive computations and
coefficient-field arguments give these optimal convergence rates:

| Allowed joint degree, at most | Decimal digits gained per term, asymptotically |
| ---: | ---: |
| 2 | 19.293494 |
| 4 | 26.373311 |
| 8 | 70.976946 |
| 16 | 149.503901 |
| 32 | 366.521245 |

The degree-32 formula has a construction that does not use π as an input.
Four terms give an approximation to π with certified absolute error below
10⁻¹⁴⁶³. The [report](REPORT.md#6-results-and-what-each-row-establishes) also lists
higher-degree candidates, whose optimality remains open in this project.

The arithmetic explanation matters as much as the table. Symmetries can make
apparently different constants share the same number system. Certain changes
to the underlying geometry increase convergence much more than they increase
the degree. [Section 4](REPORT.md#4-hidden-arithmetic-organization-why-the-winners-change)
explains that mechanism.

These are optima in a defined mathematical class, not among every conceivable
formula for π. The construction methods are classical. Novelty of the exact
optimization table has not been established, and the work has not been
independently refereed or formally verified in a proof assistant.

## Fewer terms does not necessarily mean less time

At a target absolute error below 10⁻¹⁰⁰⁰⁰, the default implementation gives:

| Identity | Terms | Coefficient setup | Sum and error bound | Total |
| --- | ---: | ---: | ---: | ---: |
| Ramanujan (1103) | 1,253 | 0.00007 s | 0.67815 s | 0.67824 s |
| Chudnovsky | 706 | 0.00014 s | 0.41704 s | 0.41719 s |
| CM degree 32 | 28 | 0.08693 s | 0.01498 s | 0.10271 s |

These are medians of five trials on one shared Linux machine, with one FLINT
thread. Every trial rebuilds its coefficients. All three use the same simple
summation and rigorous error-bound routine; none uses optimized binary splitting.

Finding all 32 polynomial roots at full precision costs much more than refining
the one we need. Isolating it at low precision and then using certified Newton
refinement reduces the degree-32 total from about 3.20 seconds to 0.103 seconds
in a comparison run. That beats this simple Chudnovsky implementation at 10,000
digits. At 1,000 digits, Chudnovsky is still faster: about 1.56 ms versus 11.4 ms.
Digits per term, setup cost, and target accuracy all matter. This is a comparison
of these implementations, not a speed record or a benchmark against optimized
π software.

[Full method, environment, before/after comparison, and 1,000-digit results](REPORT.md#practical-evaluation-cost)
are in the report; [raw samples](results/benchmark_newton.json) are included.

## Technical scope

For readers familiar with hypergeometric series, the fixed normalized kernels
and identities are:

```math
c_{s,n}=\frac{(s)_n(1/2)_n(1-s)_n}{(n!)^{3}},
\qquad s\in\{1/6,1/4,1/3,1/2\},
\qquad \frac{1}{\pi}=\sum_{n=0}^{\infty}(A+Bn)c_{s,n}x^{n}.
```

The search ranges over canonical complex multiplication (CM) constructions,
including conjugate parameters and nonmaximal quadratic orders. CM means the
modular parameter is an imaginary quadratic point. The cost is
$`[\mathbb{Q}(A,B,x):\mathbb{Q}]`$, after absorbing every outside prefactor;
terms are not regrouped. The objective is $`-\log_{10}|x|`$ for $`0<|x|<1`$.

[REPORT.md](REPORT.md) gives the finite reduction, coefficient-field descent,
exact construction, completeness ledger, benchmarks, literature comparison,
and open tasks. The current completeness result reaches degree 32; the larger
candidate rows do not extend that completeness claim.

## Reproduce

Use Python 3.11 or later. From the repository root:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python code/verify.py
```

On Windows, activate the environment with `.venv\Scripts\activate` instead.
For an exhaustive replay or a timing comparison:

```sh
python code/verify.py --full
python code/benchmark.py --digits 1000 10000 --repeats 5
```

The quick verifier checks the degree-32 parameter/radical certificate, the
degree-16 ranking, genus-character examples, the modular polynomial through an
exact Sturm-bound test, and the algebraic construction and error bounds for the
degree-32 identity. The full command reconstructs the candidate universe and
replays the exhaustive degree-32 exclusion; it can take tens of minutes.

GitHub Actions runs regression tests, quick verification, and a benchmark smoke
test on Python 3.11 and 3.12. It checks that reproduction leaves the checkout
clean, then builds and verifies a downloadable research bundle. To request an
exhaustive replay, use **Actions → Verify → Run workflow** and enable `full`.
Version tags publish verified bundles as GitHub Releases; see
[verification and releases](RELEASING.md). Timings are not CI performance thresholds.

Run without `python -O`: the programs intentionally use assertions. No network
access is needed after dependencies are installed. The benchmark and verifier
write fresh results under ignored `results/latest/`, leaving the tracked
snapshots in `results/` unchanged. The verifier prints its own fresh run directory;
the benchmark defaults to `results/latest/benchmark.json`. Use `--output PATH`
to select another benchmark destination. Published snapshots are updated only
by deliberately reviewing and copying a run's output.

## Repository guide and trust assumptions

| Location | Contents |
| --- | --- |
| [REPORT.md](REPORT.md) | Mathematical arguments, results, timings, sources, and limits |
| [EXPLORATIONS.md](EXPLORATIONS.md) | Geometric motivation, research questions, and open directions |
| [RELEASING.md](RELEASING.md) | CI tests, research bundles, and tagged releases |
| [tests/](tests/) | Arithmetic, root-refinement, and output-isolation regression tests |
| [code/](code/) | Construction, enumeration, exact checks, and benchmark |
| [data/](data/) | Fixed inputs and search records |
| [results/](results/) | Certificates, exclusion ledgers, timing samples, and logs |
| [MANIFEST.json](MANIFEST.json) | SHA-256 hashes of source and fixed inputs |

The computations rely on the cited classification and CM theorems, FLINT's
exact polynomial arithmetic and certified Hilbert class polynomials, and Arb
ball arithmetic. The manifest detects changed inputs; it is not itself a proof.
Generated results are excluded from the manifest because timings and interval
representations can vary. This is not a cross-CAS replication.

[CONTRIBUTING.md](CONTRIBUTING.md) sets the evidence required for additional
results. Research snapshot: 8 October 2026.
