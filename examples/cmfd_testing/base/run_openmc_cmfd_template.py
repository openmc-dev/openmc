import openmc
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD

# Initialize CMFD Mesh
cmfd_mesh = openmc.CMFDMesh()
cmfd_mesh.lower_left = {cmfd_mesh_ll}
cmfd_mesh.upper_right = {cmfd_mesh_ur}
cmfd_mesh.dimension = {cmfd_mesh_dim}
cmfd_mesh.albedo = {albedo}

# Initialize CMFDRun object
cmfd_run = openmc.CMFDRun()

# Set all runtime parameters (cmfd_mesh, tolerances, tally_resets, etc)
# All error checking done under the hood when setter function called
cmfd_run.cmfd_mesh = cmfd_mesh
cmfd_run.cmfd_begin = {cmfd_begin}
cmfd_run.cmfd_display = 'dominance'
cmfd_run.cmfd_feedback = True

# Run CMFD
cmfd_run.run(vectorized={vectorized}, omp_num_threads={omp_num_threads})
