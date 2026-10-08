# Verification and releases

The repository delivers a research bundle: source, documentation, fixed inputs,
and saved certificates. It does not deploy a website or publish a Python package.

## Every branch push and pull request

The **Verify** workflow runs on Python 3.11 and 3.12. It executes the arithmetic
regression tests and the quick certificate verifier. The tests check classical
identities and the degree-32 construction against an independent Arb value of π
at 500-digit absolute accuracy. A negative control perturbs a coefficient by
10⁻⁴⁰⁰ and checks that the independent π comparison detects the error.
The suite also compares certified Newton refinement with independent root
isolation and checks that fresh verification runs cannot consume stale outputs.
A small benchmark exercises the CLI, and a clean-checkout check catches any
accidental writes to tracked files or unignored output directories.

After both versions pass, a separate job builds a ZIP from the checked-out Git
commit. It extracts the ZIP into a fresh directory and reruns the tests and quick
verifier there, ensuring the downloadable bundle contains the required files.
The `research-bundle` Actions artifact is retained for 14 days and contains:

- `ramanujan-series.zip`, with the committed files under `ramanujan-series/`;
- `COMMIT.txt`, identifying the source commit;
- `SHA256SUMS`, covering both files.

On a pull request, the commit may be GitHub's synthetic merge commit. The
provenance file records the actual checked-out commit, rather than assuming it
is the contributor's branch head.

## Exhaustive replay

Use **Actions → Verify → Run workflow**, choose a branch, and enable `full`.
This adds the exhaustive replay on Python 3.12. Ordinary pull requests retain
the quick checks. Timings and benchmark comparisons are not CI pass/fail gates.

## Publish a version

After the workflow changes are merged and the intended source is on `main`, a
maintainer can publish by pushing a version tag. For example, after selecting
the appropriate unused version number:

```sh
git switch main
git pull --ff-only
git tag -a v0.1.0 -m "Research snapshot v0.1.0"
git push origin v0.1.0
```

A pushed `v*` tag triggers **Release**. It reuses the same tests, verification,
and packaging workflow with **exhaustive replay enabled**. Only after all those
jobs pass does it publish a GitHub Release with the ZIP, checksums, source commit,
and generated release notes. A release can therefore take tens of minutes.

The publish job checks that the downloaded bundle belongs to the tagged commit.
Only that job has `contents: write`; it uses GitHub's built-in token, so no
additional publishing secret is needed. The release command requires an
existing tag and does not overwrite an existing release. If a run fails,
inspect its logs before rerunning it; do not treat a failed gate as verification.

Publishing requires a separate version-tag push. A release packages the report's stated claims and
limitations; it does not confer peer review or mathematical novelty.

## Local checks

```sh
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python code/verify.py
```

When changing tracked source or documentation, refresh the corresponding
`MANIFEST.json` hashes. Do not use `python -O`: the certificate programs use
assertions. Full replay remains necessary for changes to completeness logic,
as described in [CONTRIBUTING.md](CONTRIBUTING.md).
