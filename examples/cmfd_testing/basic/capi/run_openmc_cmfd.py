import openmc
import numpy as np

# Initialize CMFD Mesh
cmfd_mesh = openmc.CMFDMesh()
cmfd_mesh.lower_left = [-10, -1, -1]
cmfd_mesh.upper_right = [10, 1, 1]
cmfd_mesh.dimension = [10, 1, 1]
cmfd_mesh.albedo = [0., 0., 1., 1., 1., 1.]
cmfd_mesh.energy = [0.001, 0.01, 0.1]

# Initialize CMFDRun object
cmfd_run = openmc.CMFDRun()

# Set all runtime parameters (cmfd_mesh, tolerances, tally_resets, etc)

# All error checking done under the hood when setter function called
cmfd_run.cmfd_mesh = cmfd_mesh

# Run CMFD
cmfd_run.run()
