"""Keep published snapshots separate from freshly generated run outputs."""
import os
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[1]
RUN_ENV = 'RAMANUJAN_RESULTS_DIR'


def output_path(name):
    directory = Path(os.environ.get(RUN_ENV, ROOT / 'results' / 'latest'))
    directory.mkdir(parents=True, exist_ok=True)
    return directory / name


def generated_input(name):
    """Read a pipeline predecessor, falling back to its published snapshot.

    verify.py gives every invocation a fresh directory shared by its children,
    so an older run cannot supply a predecessor to a new verification run.
    Standalone scripts can chain outputs in results/latest/ (or RUN_ENV).
    Fixed certificates, such as candidate32.json, are read explicitly instead.
    """
    generated = output_path(name)
    return generated if generated.exists() else ROOT / 'results' / name


def fresh_run_directory():
    parent = ROOT / 'results' / 'latest'
    parent.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix='verify-', dir=parent))
