import pathlib as pl
import os
import shutil
import subprocess
import textwrap

import openmc
import pytest

from tests.testing_harness import TestHarness

def model():
  subprocess.run(['gcc','-o gen_dummt_mcpl.out','-lm','-lmcpl','gen_dummy_mcpl.c'] )
  subprocess.run(['./gen_dummy_mcpl.out'])
    
class SourceMCPLFileTestHarness(TestHarness):
  def execute_test(self):
    """Run OpenMC with the appropriate arguments and check the outputs."""
    try:
      self._create_input()
      self._run_openmc()
      self._test_output_created()
      results = self._get_results()
      self._write_results(results)
      self._compare_results()
    finally:
      self._cleanup()

  def update_results(self):
    """Update the results_true using the current version of OpenMC."""
    try:
      self._create_input()
      self._run_openmc()
      self._test_output_created()
      results = self._get_results()
      self._write_results(results)
      self._overwrite_results()
    finally:
      self._cleanup()

  def _test_output_created(self):
    """Check that the output files were created"""
    stat        epoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))
    assert len(statepoint) == 1, 'Either multiple or no statepoint files ' \
      'exist.'
    assert statepoint[0].endswith('h5'), \
      'statepoint file does not end with h5.'

def test_mcpl_source_file():
  harness = SourceMCPLFileTestHarness('statepoint.10.h5')
  harness.main()
  
