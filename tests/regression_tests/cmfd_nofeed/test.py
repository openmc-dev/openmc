from tests.testing_harness import CMFDTestHarness
from openmc import cmfd
import numpy as np


def test_cmfd_nofeed():
    """ Test 1 group CMFD solver without CMFD feedback"""
    # Initialize and set CMFD mesh
    cmfd_mesh = cmfd.CMFDMesh()
    cmfd_mesh.lower_left = -10.0, -1.0, -1.0
    cmfd_mesh.upper_right = 10.0, 1.0, 1.0
    cmfd_mesh.dimension = 10, 1, 1
    cmfd_mesh.albedo = 0.0, 0.0, 1.0, 1.0, 1.0, 1.0

    # Initialize and run CMFDRun object
    cmfd_run = cmfd.CMFDRun()
    cmfd_run.mesh = cmfd_mesh
    cmfd_run.begin = 5
    cmfd_run.display = {'dominance': True}
    cmfd_run.feedback = False
    cmfd_run.gauss_seidel_tolerance = [1.e-15, 1.e-20]
    cmfd_run.run()

    # Create output string of all CMFD results to pass into testing harness
    outstr = 'cmfd indices\n'
    outstr += '\n'.join(['{0:12.6E}'.format(x) for x in cmfd_run.indices])
    outstr += '\nk cmfd\n'
    outstr += '\n'.join(['{0:12.6E}'.format(x) for x in cmfd_run.k_cmfd])
    outstr += '\ncmfd entropy\n'
    outstr += '\n'.join(['{0:12.6E}'.format(x) for x in cmfd_run.entropy])
    outstr += '\ncmfd balance\n'
    outstr += '\n'.join(['{0:12.6E}'.format(x) for x in cmfd_run.balance])
    outstr += '\ncmfd dominance ratio\n'
    outstr += '\n'.join(['{0:10.3E}'.format(x) for x in cmfd_run.dom])
    outstr += '\ncmfd openmc source comparison\n'
    outstr += '\n'.join(['{0:12.6E}'.format(x) for x in cmfd_run.src_cmp])
    outstr += '\ncmfd source\n'
    cmfdsrc = np.reshape(cmfd_run.cmfd_src, np.product(cmfd_run.indices),
                         order='F')
    outstr += '\n'.join(['{0:12.6E}'.format(x) for x in cmfdsrc])
    outstr += '\n'

    # Initialize and run CMFD test harness
    harness = CMFDTestHarness('statepoint.20.h5', cmfd_results=outstr)
    harness.main()
