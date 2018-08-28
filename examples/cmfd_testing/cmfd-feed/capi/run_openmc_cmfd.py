import openmc
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD

# Initialize CMFD Mesh
cmfd_mesh = openmc.CMFDMesh()
cmfd_mesh.lower_left = [-10, -1, -1]
cmfd_mesh.upper_right = [10, 1, 1]
cmfd_mesh.dimension = [10, 1, 1]
cmfd_mesh.albedo = [0., 0., 1., 1., 1., 1.]
cmfd_mesh.energy = [0., 1000000., 10000000000.]
cmfd_mesh.map = [0,1,1,1,1,1,1,1,1,0]

# Initialize CMFDRun object
cmfd_run = openmc.CMFDRun()

# Set all runtime parameters (cmfd_mesh, tolerances, tally_resets, etc)
# All error checking done under the hood when setter function called
cmfd_run.cmfd_mesh = cmfd_mesh
cmfd_run.cmfd_begin = 5
cmfd_run.cmfd_display = 'dominance'
cmfd_run.cmfd_feedback = True

# Run CMFD
cmfd_run.run()
