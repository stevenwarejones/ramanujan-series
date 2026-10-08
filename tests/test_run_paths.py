"""Ensure verification cannot pick up stale results from an earlier run."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import run_paths


class RunPathTests(unittest.TestCase):
    def test_fresh_verification_uses_snapshot_then_current_predecessor(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            snapshots=root/'results'
            snapshots.mkdir()
            (snapshots/'catalogue64.json').write_text('published')
            with patch.object(run_paths,'ROOT',root), patch.dict(os.environ,{},clear=True):
                run_paths.output_path('catalogue64.json').write_text('stale standalone run')
                first=run_paths.fresh_run_directory()
                with patch.dict(os.environ,{run_paths.RUN_ENV:str(first)}):
                    self.assertEqual(run_paths.generated_input('catalogue64.json').read_text(),'published')
                    run_paths.output_path('catalogue64.json').write_text('current full replay')
                    self.assertEqual(run_paths.generated_input('catalogue64.json').read_text(),'current full replay')
                second=run_paths.fresh_run_directory()
                self.assertNotEqual(first,second)
                with patch.dict(os.environ,{run_paths.RUN_ENV:str(second)}):
                    self.assertEqual(run_paths.generated_input('catalogue64.json').read_text(),'published')
                self.assertEqual((snapshots/'catalogue64.json').read_text(),'published')


if __name__=='__main__':
    unittest.main()
