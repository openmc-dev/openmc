import pathlib as pl
import os
import shutil
import subprocess
import textwrap
import glob
import openmc
import pytest

from tests.testing_harness import TestHarness


class SourceMCPLFileTestHarness(TestHarness):
  def execute_test(self):
    """Run OpenMC with the appropriate arguments and check the outputs."""
    try:
      self._create_input()
      self._test_input_created()
      self._run_openmc()
      self._test_output_created()
      results = self._get_results()
      self._write_results(results)
      self._compare_results()
    finally:
      self._cleanup()

  def _create_input(self):
    subprocess.run(['gcc','-o','gen_dummy_mcpl.out','gen_dummy_mcpl.c','-lm','-lmcpl'])
    subprocess.run(['./gen_dummy_mcpl.out'])

  def update_results(self):
    """Update the results_true using the current version of OpenMC."""
    try:
      self._create_input()
      self._test_input_created()
      self._run_openmc()
      self._test_output_created()
      results = self._get_results()
      self._write_results(results)
      self._overwrite_results()
    finally:
      self._cleanup()

  def _test_input_created(self):
    """Check that the input mcpl.file was generated as it should"""
    mcplfile=glob.glob(os.path.join(os.getcwd(),'source.10.mcpl'))
    assert len(mcplfile) == 1, 'Either multiple or no mcpl files ' \
      'exist.'
    assert mcplfile[0].endswith('mcpl'), \
      'output file does not end with mcpl.'

  def _test_output_created(self):
    """Check that the output files were created"""
    statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))
    assert len(statepoint) == 1, 'Either multiple or no statepoint files ' \
      'exist.'
    assert statepoint[0].endswith('h5'), \
      'statepoint file does not end with h5.'

  def _cleanup(self):
    super()._cleanup()
    source_mcpl=glob.glob(os.path.join(os.getcwd(),'source*.mcpl'))
    for f in source_mcpl:
      if (os.path.exists(f)):
        os.remove(f)

def test_mcpl_source_file():
  harness = SourceMCPLFileTestHarness('source.10.mcpl')
  harness.main()
