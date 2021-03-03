from tests.testing_harness import CMFDTestHarness
from openmc import cmfd
import numpy as np
import scipy.sparse


def test_cmfd_physical_adjoint():
    """Test physical adjoint functionality of CMFD

    This test runs CMFD with a physical adjoint calculation and asserts that
    the adjoint k-effective and flux vector are equal to the non-adjoint
    k-effective and flux vector at the last batch (equivalent for 1 group
    problems).

    """
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
    cmfd_run.solver_begin = 5
    cmfd_run.feedback = True
    cmfd_run.gauss_seidel_tolerance = [1.e-15, 1.e-20]
    cmfd_run.run_adjoint = True
    cmfd_run.adjoint_type = 'physical'
    cmfd_run.run()
    assert(np.all(cmfd_run._phi == cmfd_run._adj_phi))
    assert(cmfd_run._adj_keff == cmfd_run._keff)


def test_cmfd_math_adjoint():
    """Test mathematical adjoint functionality of CMFD

    This test runs CMFD with a mathematical adjoint calculation and asserts
    that the adjoint k-effective and flux vector are equal to the non-adjoint
    k-effective and flux vector at the last batch (equivalent for 1 group
    problems).

    """
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
    cmfd_run.solver_begin = 5
    cmfd_run.feedback = True
    cmfd_run.gauss_seidel_tolerance = [1.e-15, 1.e-20]
    cmfd_run.run_adjoint = True
    cmfd_run.adjoint_type = 'math'
    cmfd_run.run()
    assert(np.all(cmfd_run._phi == cmfd_run._adj_phi))
    assert(cmfd_run._adj_keff == cmfd_run._keff)


def test_cmfd_write_matrices():
    """Test write matrices functionality of CMFD

    This test runs CMFD with feedback and loads the loss/production matrices
    and flux vector that are saved to disk, and checks to make sure these
    values are consistent with each other and simulation results.

    """
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
    cmfd_run.solver_begin = 5
    cmfd_run.display = {'dominance': True}
    cmfd_run.feedback = True
    cmfd_run.gauss_seidel_tolerance = [1.e-15, 1.e-20]
    cmfd_run.write_matrices = True
    cmfd_run.run()

    # Load loss matrix from numpy output file
    loss_np = scipy.sparse.load_npz('loss.npz').todense()
    # Load loss matrix from data file
    loss_dat = np.loadtxt("loss.dat", delimiter=',')

    # Go through each element of loss_dat and compare to loss_np
    for elem in loss_dat:
        assert(np.isclose(loss_np[int(elem[0]), int(elem[1])], elem[2]))

    # Load production matrix from numpy output file
    prod_np = scipy.sparse.load_npz('prod.npz').todense()
    # Load production matrix from data file
    prod_dat = np.loadtxt("prod.dat", delimiter=',')

    # Go through each element of prod_dat and compare to prod_np
    for elem in prod_dat:
        assert(np.isclose(prod_np[int(elem[0]), int(elem[1])], elem[2]))

    # Load flux vector from numpy output file
    flux_np = np.load('fluxvec.npy')
    # Load flux from data file
    flux_dat = np.loadtxt("fluxvec.dat", delimiter='\n')

    # Compare flux from numpy file, .dat file, and from simulation
    assert(np.all(np.isclose(flux_np, cmfd_run._phi)))
    assert(np.all(np.isclose(flux_np, flux_dat)))


def test_cmfd_feed():
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
    cmfd_run.solver_begin = 5
    cmfd_run.display = {'dominance': True}
    cmfd_run.feedback = True
    cmfd_run.gauss_seidel_tolerance = [1.e-15, 1.e-20]
    cmfd_run.run()

    # Initialize and run CMFD test harness
    harness = CMFDTestHarness('statepoint.20.h5', cmfd_run)
    harness.main()

def test_cmfd_multithread():
    """Test 1 group CMFD solver with all available threads"""
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
    cmfd_run.solver_begin = 5
    cmfd_run.display = {'dominance': True}
    cmfd_run.feedback = True
    cmfd_run.gauss_seidel_tolerance = [1.e-15, 1.e-20]
    cmfd_run.use_all_threads = True
    cmfd_run.run()

    # Initialize and run CMFD test harness
    harness = CMFDTestHarness('statepoint.20.h5', cmfd_run)
    harness.main()
