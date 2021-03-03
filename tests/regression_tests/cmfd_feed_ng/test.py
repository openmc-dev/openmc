from tests.testing_harness import CMFDTestHarness
from openmc import cmfd
import numpy as np


def test_cmfd_feed_ng():
    """Test n group CMFD solver with CMFD feedback"""
    # Initialize and set CMFD mesh
    cmfd_mesh = cmfd.CMFDMesh()
    cmfd_mesh.lower_left = (-1.25984, -1.25984, -1.0)
    cmfd_mesh.upper_right = (1.25984, 1.25984, 1.0)
    cmfd_mesh.dimension = (2, 2, 1)
    cmfd_mesh.energy = (0.0, 0.625, 5.53080, 20000000)
    cmfd_mesh.albedo = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

    # Initialize and run CMFDRun object
    cmfd_run = cmfd.CMFDRun()
    cmfd_run.mesh = cmfd_mesh
    cmfd_run.reset = [5]
    cmfd_run.tally_begin = 10
    cmfd_run.solver_begin = 10
    cmfd_run.display = {'dominance': True}
    cmfd_run.feedback = True
    cmfd_run.downscatter = True
    cmfd_run.gauss_seidel_tolerance = [1.e-15, 1.e-20]
    cmfd_run.run()

    # Initialize and run CMFD test harness
    harness = CMFDTestHarness('statepoint.20.h5', cmfd_run)
    harness.main()
