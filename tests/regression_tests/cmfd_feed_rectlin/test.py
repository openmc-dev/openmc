from tests.testing_harness import CMFDTestHarness
from openmc import cmfd
import numpy as np
import scipy.sparse


def test_cmfd_feed_rectlin():
    """Test 1 group CMFD solver with CMFD feedback"""
    # Initialize and set CMFD mesh
    cmfd_mesh = cmfd.CMFDMesh()
    cmfd_mesh.mesh_type = 'rectilinear'
    x_grid = [-10., -9., -7., -6., -4., -3., -1., 0., 1., 3., 4., 6., 7., 9., 10.]
    y_grid = [-1., 1.]
    z_grid = [-1., 1.]
    cmfd_mesh.grid = [x_grid, y_grid, z_grid]
    cmfd_mesh.albedo = (0.0, 0.0, 1.0, 1.0, 1.0, 1.0)

    # Initialize and run CMFDRun object
    cmfd_run = cmfd.CMFDRun()
    cmfd_run.mesh = cmfd_mesh
    cmfd_run.tally_begin = 5
    cmfd_run.solver_begin = 5
    cmfd_run.display = {'dominance': True}
    cmfd_run.feedback = True
    cmfd_run.gauss_seidel_tolerance = [1.e-15, 1.e-20]
    cmfd_run.run()

    # Initialize and run CMFD test harness
    harness = CMFDTestHarness('statepoint.20.h5', cmfd_run)
    harness.main()
