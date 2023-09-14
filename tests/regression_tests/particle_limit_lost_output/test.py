from tests.testing_harness import TestHarness
import functools
import glob
import os
from pytest import approx

def _cleanup_wrapper(func):
  @functools.wraps(func)
  def _cleanup_particles():
    func()
  return _cleanup_particles

class TestHarness_limit_particle_output(TestHarness):
    #redefine the check output method to do something else
    def _check_output_limit_lost(self):
      glb=glob.glob('particle_[0-9]_*.h5')
      assert (5== approx(len(glb),abs=2))

    def _cleanup(self):
      super()._cleanup()
      output=glob.glob('particle_[0-9]_*.h5')
      for f in output:
        if os.path.exists(f):
          os.remove(f)

    def _test_output_created(self):
      super()._test_output_created()
      self._check_output_limit_lost()


def test_limit_particle_output():
    harness = TestHarness_limit_particle_output('statepoint.10.h5')
    harness.main()
