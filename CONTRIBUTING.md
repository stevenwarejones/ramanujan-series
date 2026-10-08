# Contributing

The research question is how convergence depends on the full algebraic degree
of a Ramanujan-type identity. Contributions should make that comparison more
precise, reproducible, or understandable.

For a new result, include:

1. **The search universe.** Specify the series kernel, allowed constructions,
   normalization, embeddings, and whether the search is exhaustive or bounded.
2. **Every coefficient.** Measure the joint field degree after absorbing any
   outside prefactor. Do not replace it with the largest individual degree.
3. **An explanation.** Identify the modular or arithmetic mechanism behind the
   result and distinguish established theory from new observations.
4. **Checkable evidence.** Supply exact defining equations, branch isolation,
   and certificates or rigorous error bounds. Numerical agreement alone does
   not establish an identity, a minimal polynomial, or optimality.

Separate candidate lower bounds from proved optima, and digits per term from
end-to-end runtime. Cite primary literature and avoid novelty claims based only
on not finding a matching formula in a search.

Run `python -m unittest discover -s tests -v` and `python code/verify.py` for
routine changes. Changes to enumeration,
exclusion logic, or completeness claims also require
`python code/verify.py --full`. Do not use Python's `-O` option: the verifier
uses assertions.

`MANIFEST.json` records SHA-256 hashes of source files, repository documentation,
tests, workflows, and fixed inputs. Update the affected hashes when
intentionally changing those files. It excludes itself and the regenerated
`results/` snapshots. A matching hash detects a changed input; it does not prove
that the input or the mathematics is correct.

Fresh output goes under ignored `results/latest/`. Each verifier invocation
creates and prints a fresh subdirectory shared by its child scripts, so a full
replay consumes its own newly generated catalogue and selections without picking
up a previous run. Standalone scripts write directly in `results/latest/`; a
pipeline consumer prefers a predecessor there, falling back to the published
snapshot if absent. Set `RAMANUJAN_RESULTS_DIR` to give a standalone pipeline its
own directory. The defining `results/candidate32.json` certificate remains a
fixed input; regenerating it does not silently replace the published identity.
Review and deliberately copy generated files when updating a tracked snapshot.

See [RELEASING.md](RELEASING.md) for the CI checks and release process.

For documentation figures, use the optional rendering dependencies and workflow
in [assets/figures/README.md](assets/figures/README.md). Preview generation writes
to ignored `results/latest/figures/`; publishing revised SVGs is an explicit step.

Record any substantive mathematical change in `REPORT.md`, including its
scope, supporting sources, verification performed, and remaining limitations.

For GitHub math rendering, use fenced `math` blocks and GitHub's dollar/backtick
inline delimiters. Use explicit braces for superscripts and avoid Markdown
emphasis characters in TeX; for a superscript star, use `^{\ast}`. Use
`\mathrm{Im}` for the imaginary-part label. This keeps Markdown preprocessing
from changing the expressions before the math renderer receives them.
See [GitHub's math formatting guide](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).
