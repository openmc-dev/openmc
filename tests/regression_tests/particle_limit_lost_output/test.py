from tests.testing_harness import TestHarness
import glob
import os
from tests.regression_tests import config

class TestHarness_limit_particle_output(TestHarness):

    def mpi_np(self):
      #returns the number of mpi-tasks for this run.
      if config['mpi']:
        try:
          mpi_np=int(config['mpi_np'])
        except:
          mpi_np=1
        return mpi_np
      else:
        return 1

    def omp_nt(self):
      #get the number of concurrent omp-threads - as this
      #is the maximum number of extra particle-files that may
      #be written per task.
      try:
        omp_nt=int(os.environ['OMP_NUM_THREADS'])
      except:
        omp_nt=len(os.sched_getaffinity(0))
      return omp_nt

    #redefine the check output method to do something else
    def _check_output_limit_lost(self):
      glb=glob.glob('particle_[0-9]_*.h5')
      assert (len(glb) <= 5*self.mpi_np() + self.omp_nt())

    def _cleanup(self):
      super()._cleanup()
      output=glob.glob('particle_[0-9]*_*.h5')
      for f in output:
        if os.path.exists(f):
          os.remove(f)

    def _test_output_created(self):
      super()._test_output_created()
      self._check_output_limit_lost()

def test_limit_particle_output():
    harness = TestHarness_limit_particle_output('statepoint.10.h5')
    harness.main()
