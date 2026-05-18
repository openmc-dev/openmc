"""Runs OpenMC on a 2x2x2 cube with different temperature values for each cell.
The cube has reflective boundary conditions.

The entire geometry is defined as a DAGMC geometry which includes the
boundary conditions.

A "results.dat" file is created and can be compared to the "results_true.dat"
file from the corresponding case.

"""

import openmc
import numpy as np

model = openmc.Model()

# Material
mat = openmc.Material(name="fuel")
mat.add_nuclide("U235", 0.2)
mat.add_nuclide("U238", 0.8)
mat.add_element("O", 2.0)
mat.add_element("H", 4.0)
mat.set_density("g/cm3", 5.0)
mat.add_s_alpha_beta('c_H_in_H2O')
model.materials = openmc.Materials([mat])

# Create mesh for tallying
dim = 2
lower_left = (0.0, 0.0, 0.0)
upper_right = (5.0, 5.0, 5.0)
mesh = openmc.RegularMesh()
mesh.lower_left = lower_left
mesh.upper_right = upper_right
mesh.dimension = (dim, dim, dim)

# Read DAGMC universe
dag_univ = openmc.DAGMCUniverse(filename="single_cube_multicell_reflecting.h5m", auto_geom_ids=True)
model.geometry = openmc.Geometry(dag_univ)

# Settings
settings = openmc.Settings()
settings.batches = 20
settings.particles = 200
spatial_dist = openmc.stats.Box(lower_left, upper_right)
settings.source = openmc.IndependentSource(
    space=spatial_dist, constraints={"fissionable": True})
settings.temperature = {'tolerance': 1000, 'multipole': True}
model.settings = settings

# Add tallies
mesh_filter = openmc.MeshFilter(mesh)
mesh_tally = openmc.Tally(name="total reaction rate")
mesh_tally.filters = [mesh_filter]
mesh_tally.scores = ["total"]
tallies = openmc.Tallies([mesh_tally])
model.tallies = tallies

model.run()

with openmc.StatePoint(f"statepoint.20.h5") as sp:

    outstr = 'k-combined:\n'
    form = '{0:12.6E} {1:12.6E}\n'
    outstr += form.format(sp.keff.n, sp.keff.s)

    for i, tally_ind in enumerate(sp.tallies):
        tally = sp.tallies[tally_ind]
        results = np.zeros((tally.sum.size * 2, ))
        results[0::2] = tally.sum.ravel()
        results[1::2] = tally.sum_sq.ravel()
        results = ['{0:12.6E}'.format(x) for x in results]

        outstr += 'tally {}:\n'.format(i + 1)
        outstr += '\n'.join(results) + '\n'

    with open("results.dat", "w") as res:
        res.write(outstr)
