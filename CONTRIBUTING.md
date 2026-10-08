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

Run `python code/verify.py` for routine changes. Changes to enumeration,
exclusion logic, or completeness claims also require
`python code/verify.py --full`. Do not use Python's `-O` option: the verifier
uses assertions.

`MANIFEST.json` records SHA-256 hashes of source files, repository documentation,
the verification workflow, and fixed inputs. Update the affected hashes when
intentionally changing those files. It excludes itself and the regenerated
`results/` directory. A matching hash detects a changed input; it does not prove
that the input or the mathematics is correct.

Record any substantive mathematical change in `REPORT.md`, including its
scope, supporting sources, verification performed, and remaining limitations.
