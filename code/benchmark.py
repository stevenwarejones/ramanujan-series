"""Setup-inclusive, certified direct-summation timings for three pi identities.

This measures the implementations here, not optimized pi algorithms. Every trial
reconstructs its coefficients. Imports, theorem replay, and the independent pi
reference check are outside the timer. No timed path uses pi as an input.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import statistics
import time

import flint
from flint import arb, ctx
from identity32 import construct
from run_paths import output_path

ROOT = Path(__file__).resolve().parents[1]
METHODS = ("ramanujan-1103", "chudnovsky", "cm-degree-32")


def coefficients(method, dps, root_method='newton'):
    ctx.dps = dps
    if method == "ramanujan-1103":
        factor = 2 * arb(2).sqrt() / 9801
        return factor * 1103, factor * 26390, arb(1) / 99**4, arb(1) / 4
    if method == "chudnovsky":
        factor = 12 / (arb(640320) * arb(640320).sqrt())
        return factor * 13591409, factor * 545140134, -arb(1728) / 640320**3, arb(1) / 6
    if method == "cm-degree-32":
        A, B, x, _ = construct(dps, root_method=root_method)
        return A, B, x, arb(1) / 4
    raise ValueError(method)


def certified_sum(A, B, x, s, digits):
    """Stop on an absolute-error proof, including coefficient rounding and tail.

    c[n+1]/c[n] = (n+s)(n+1/2)(n+1-s)/(n+1)^3 < 1.
    For n >= m, the absolute summand ratio is at most
    |x| (1 + B/(A+B*m)). Hence first_omitted/(1-ratio) bounds the tail.
    """
    assert A > 0 and B > 0 and abs(x) < 1 and 0 < s < 1
    threshold = arb(10) ** (-digits)
    total, c, power = arb(0), arb(1), arb(1)
    for n in range(2 * digits + 100):
        total += (A + B * n) * c * power
        m = n + 1
        c *= (n + s) * (n + arb(1)/2) * (n + 1 - s) / (n + 1)**3
        power *= x
        first = abs((A + B * m) * c * power)
        ratio = abs(x) * (1 + B / (A + B * m))
        assert ratio < 1
        tail = first / (1 - ratio)
        enclosure = total + arb(0, tail.upper())
        if enclosure.contains(0):
            continue
        pi_interval = 1 / enclosure
        approximation = (1 / total).mid()
        error = abs(pi_interval.mid() - approximation) + pi_interval.rad()
        if error < threshold:
            return m, error, approximation, pi_interval
    raise ArithmeticError("Precision budget insufficient; target error not certified")


def trial(method, digits, guard, root_method='newton'):
    start = time.perf_counter()
    A, B, x, s = coefficients(method, digits + guard, root_method)
    setup_done = time.perf_counter()
    terms, error, approximation, pi_interval = certified_sum(A, B, x, s, digits)
    end = time.perf_counter()
    # Independent check only; excluded from coefficient construction and timing.
    reference = arb.pi()
    assert pi_interval.contains(reference)
    assert abs(approximation - reference) < arb(10) ** (-digits)
    return dict(terms=terms, setup_seconds=setup_done-start,
                summation_and_bound_seconds=end-setup_done, total_seconds=end-start,
                absolute_error_upper_bound=str(error.upper()))


def cpu_model():
    path = Path('/proc/cpuinfo')
    if path.exists():
        for line in path.read_text().splitlines():
            if line.startswith('model name'):
                return line.split(':', 1)[1].strip()
    return platform.processor() or 'unreported'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--digits', nargs='+', type=int, default=[1000, 10000])
    ap.add_argument('--repeats', type=int, default=5)
    ap.add_argument('--guard', type=int, default=128)
    ap.add_argument('--output', type=Path, default=output_path('benchmark.json'))
    ap.add_argument('--root-method', choices=['newton', 'all-roots'], default='newton',
                    help='degree-32 root construction; all-roots retains the baseline')
    args = ap.parse_args()
    if args.repeats < 1 or args.guard < 32 or any(d < 1 for d in args.digits):
        ap.error('positive digits/repeats and at least 32 guard digits are required')
    ctx.threads = 1
    rows = []
    for digits in args.digits:
        # Warm imports, library paths, and filesystem caches. Rebuild coefficients
        # for every measured trial; no coefficient value is carried between trials.
        for method in METHODS:
            trial(method, digits, args.guard, args.root_method)
        samples = {method: [] for method in METHODS}
        for repeat in range(args.repeats):
            # Rotate order to reduce consistent first/last-position bias.
            order = METHODS[repeat % len(METHODS):] + METHODS[:repeat % len(METHODS)]
            for method in order:
                samples[method].append(trial(method, digits, args.guard, args.root_method))
        for method in METHODS:
            trials = samples[method]
            assert len({t['terms'] for t in trials}) == 1
            row = dict(method=method, target_absolute_error_digits=digits,
                       working_decimal_digits=digits+args.guard,
                       terms=trials[0]['terms'], trials=trials)
            for field in ('setup_seconds', 'summation_and_bound_seconds', 'total_seconds'):
                row['median_'+field] = statistics.median(t[field] for t in trials)
            row['min_total_seconds'] = min(t['total_seconds'] for t in trials)
            row['max_total_seconds'] = max(t['total_seconds'] for t in trials)
            rows.append(row)
            print(json.dumps({k:v for k,v in row.items() if k != 'trials'}), flush=True)
    inputs = ['code/benchmark.py', 'code/identity32.py', 'code/run_paths.py', 'requirements.txt',
              'results/candidate32.json', 'data/phi_j_17.txt']
    result = dict(status='PASS', timestamp_utc=datetime.now(timezone.utc).isoformat(),
                  python=platform.python_version(), python_flint=flint.__version__,
                  platform=platform.platform(), cpu=cpu_model(), flint_threads=ctx.threads,
                  repeats=args.repeats, warmups_per_method_and_target=1,
                  guard_digits=args.guard, degree32_root_method=args.root_method,
                  timed='coefficient reconstruction plus direct summation, tail and rounding bounds',
                  excluded='interpreter/import startup, theorem replay, independent pi check, JSON output',
                  binary_splitting=False, cached_coefficient_values=False,
                  note='Operating-system/library caches are warm; shared-machine load can vary.',
                  sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
                  rows=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
