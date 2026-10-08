# Ramanujan-Type Series for 1/π

**Exact arithmetic, convergence, and algebraic complexity.**

How quickly can a Ramanujan-type series converge when the algebraic complexity
of **all** its coefficients is bounded? This repository studies that question
in four classical families arising from complex multiplication (CM).

## Current results

The current computation gives the following optimal rates **within the four
specified canonical CM families**. A degree budget means *at most* that degree.

| Full coefficient degree budget | Decimal digits per term (asymptotic) |
| --- | ---: |
| 2 | 19.293494 |
| 4 | 26.373311 |
| 8 | 70.976946 |
| 16 | 149.503901 |
| 32 | 366.521245 |

The degree-32 identity has an algebraic construction that does not use π as an
input. Four terms give an approximation to π with certified absolute error below
10⁻¹⁴⁶³. The report also includes higher-degree candidates; their optimality has
**not** been established. These are convergence rates, not runtime benchmarks.

See [the research report](REPORT.md) for definitions, derivations, literature,
the completeness argument, and the limits of the claims. General construction
methods are classical; novelty of the exact optimization table remains unverified.

## Scope

Research snapshot, 8 October 2026. Read [REPORT.md](REPORT.md) for scope, proofs, results,
literature comparison, and remaining gaps. This is a reproducible research package,
not a claim of a new general construction of pi series.

The fixed kernels are `(s)_n (1/2)_n (1-s)_n/(n!)^3` for
`s=1/6,1/4,1/3,1/2`. Cost is the **joint** degree `[Q(A,B,x):Q]` after
writing the right side as `sum (A+Bn)c_n x^n = 1/pi`. Terms are not grouped.

## Reproduce

Use Python 3.11 or later. From the repository root:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python code/verify.py
```

On Windows, activate the environment with `.venv\Scripts\activate` instead.
For the exhaustive replay:

```sh
python code/verify.py --full
```

The quick command checks the degree-32 parameter/radical certificate, the
previous degree-16 ranking, the genus-character examples, the modular polynomial
by an exact Sturm-bound test, and the pi-free algebraic construction and error
bounds for the new degree-32 identity. The full command also reconstructs the
entire low-degree candidate universe and reruns the exhaustive degree-32
exclusion. Full verification can take tens of minutes depending on hardware.

The manifest covers source code and fixed input data, not regenerated output
files whose timings and ball precision strings can vary between runs. Run each
command without `python -O`: the programs intentionally use assertions.

GitHub Actions runs the quick verification on pushes and pull requests. An
exhaustive replay can be requested through **Actions → Verify → Run workflow**
by enabling `full`. The local commands work without GitHub.

## Repository layout

- [`code/`](code/): exact arithmetic, construction, enumeration, and verification.
- [`data/`](data/): fixed inputs and earlier search records.
- [`results/`](results/): saved certificates, exclusion ledgers, and run logs.
- [`REPORT.md`](REPORT.md): mathematical arguments, sources, and research status.
- [`MANIFEST.json`](MANIFEST.json): SHA-256 hashes of source and fixed inputs.

See [CONTRIBUTING.md](CONTRIBUTING.md) for the evidence required for new results.

## What is trusted

The classification theorem of Watkins, standard complex multiplication and
modular-function identities described in the report, and FLINT's exact integer
factorization, certified Hilbert class polynomial computation, and Arb ball
arithmetic. This is not a Lean/Coq formalization and is not a cross-CAS replication.

`data/` contains the earlier search records and the downloaded classical modular
polynomial. `results/` contains machine-readable checks and exclusion ledgers.
The modular polynomial is independently checked from its integer coefficients.
No external network access is required to run the package after dependencies
have been installed.
