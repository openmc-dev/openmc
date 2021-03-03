from tests.testing_harness import CMFDTestHarness
from openmc import cmfd


def test_cmfd_feed_rolling_window():
    """Test 1 group CMFD solver with CMFD feedback"""
    # Initialize and set CMFD mesh
    cmfd_mesh = cmfd.CMFDMesh()
    cmfd_mesh.lower_left = (-10.0, -1.0, -1.0)
    cmfd_mesh.upper_right = (10.0, 1.0, 1.0)
    cmfd_mesh.dimension = (10, 1, 1)
    cmfd_mesh.albedo = (0.0, 0.0, 1.0, 1.0, 1.0, 1.0)

    # Initialize and run CMFDRun object
    cmfd_run = cmfd.CMFDRun()
    cmfd_run.mesh = cmfd_mesh
    cmfd_run.tally_begin = 5
    cmfd_run.solver_begin = 10
    cmfd_run.feedback = True
    cmfd_run.gauss_seidel_tolerance = [1.e-15, 1.e-20]
    cmfd_run.window_type = 'rolling'
    cmfd_run.window_size = 5
    cmfd_run.run()

    # Initialize and run CMFD test harness
    harness = CMFDTestHarness('statepoint.20.h5', cmfd_run)
    harness.main()
