from math import log10

import random
import numpy as np
import matplotlib.pyplot as plt
import openmc
import sys

###############################################################################
# Create materials for the problem

uo2 = openmc.Material(name='UO2 fuel at 2.4% wt enrichment')
uo2.set_density('g/cm3', 10.29769)
uo2.add_element('U', 1., enrichment=2.4)
uo2.add_element('O', 2.)

iron = openmc.Material(name="Iron")
iron.set_density('g/cm3', 7.874)
iron.add_element('Fe', 1.0, 'wo')

air = openmc.Material(name="Air")
air.set_density('g/cm3', 0.001225)
air.add_element('O', 0.23, 'wo')
air.add_element('N', 0.75, 'wo')
air.add_element('Ar', 0.02, 'wo')

helium = openmc.Material(name='Helium for gap')
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)

zircaloy = openmc.Material(name='Zircaloy 4')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_element('Sn', 0.014  , 'wo')
zircaloy.add_element('Fe', 0.00165, 'wo')
zircaloy.add_element('Cr', 0.001  , 'wo')
zircaloy.add_element('Zr', 0.98335, 'wo')

borated_water = openmc.Material(name='Borated water')
borated_water.set_density('g/cm3', 0.740582)
borated_water.add_element('B', 4.0e-5)
borated_water.add_element('H', 5.0e-2)
borated_water.add_element('O', 2.4e-2)
borated_water.add_s_alpha_beta('c_H_in_H2O')

# Collect the materials together and export to XML
materials = openmc.Materials([uo2, helium, zircaloy, borated_water, iron, air])
materials.export_to_xml()

###############################################################################
# Define problem geometry

objects = []

box_size=100.0
sphere_radius=0.5

box_i = openmc.openmc.model.RectangularParallelepiped(0.0, box_size, 0.0, box_size, 0.0, box_size, boundary_type='reflective')
box_o = openmc.openmc.model.RectangularParallelepiped(-0.1, box_size + 0.1, -0.1, box_size + 0.1, -0.1, box_size + 0.1, boundary_type='reflective')

water_region=-box_i

for i in range(0, 1000, 1):
    pos = [random.uniform(sphere_radius, box_size-sphere_radius), random.uniform(sphere_radius, box_size-sphere_radius), random.uniform(sphere_radius, box_size-sphere_radius)]
    sphere = openmc.Sphere(x0=pos[0], y0=pos[1], z0=pos[2], r=sphere_radius, name='Iron Sphere')
    water_region = water_region & +sphere
    iron_sphere = openmc.Cell(fill=uo2, region=-sphere)
    objects.extend([iron_sphere])

boundary=openmc.Cell(fill=zircaloy, region=+box_i & -box_o)
water = openmc.Cell(fill=air, region=water_region)
objects.extend([boundary, water])

# Create a geometry and export to XML
geometry = openmc.Geometry(objects)
geometry.export_to_xml()

###############################################################################
# Define problem settings

# Indicate how many particles to run
settings = openmc.Settings()
settings.batches = 11 # I am only using 11 batches to speed up debugging
settings.inactive = 10
settings.particles = int(sys.argv[1])

# Create an initial uniform spatial source distribution over fissionable zones
lower_left = (box_size / 2 - 0.001, box_size / 2 - 0.001, box_size / 2 - 0.001)
upper_right = (box_size / 2 + 0.001, box_size / 2 + 0.001, box_size / 2 + 0.001)

lower_left = (0, 0, 0)
upper_right = (box_size, box_size, box_size)

source = openmc.Source()
source.space = openmc.stats.Point(xyz=(box_size / 2, box_size / 2, box_size / 2))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([10.0e6], [1.0])
source.time = openmc.stats.Uniform(0, 1e-6)
settings.source = source

# For source convergence checks, add a mesh that can be used to calculate the
# Shannon entropy
entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = lower_left
entropy_mesh.upper_right = upper_right
entropy_mesh.dimension = (10, 10, 10)
settings.entropy_mesh = entropy_mesh
settings.export_to_xml()

###############################################################################
# Define tallies

# Create a mesh that will be used for tallying
mesh = openmc.RegularMesh()
mesh.dimension = (100, 100)
mesh.lower_left = (0.0, 0.0)
mesh.upper_right = (box_size, box_size)

# Create a mesh filter that can be used in a tally
mesh_filter = openmc.MeshFilter(mesh)

# Now use the mesh filter in a tally and indicate what scores are desired
mesh_tally = openmc.Tally(name="Mesh tally")
mesh_tally.filters = [mesh_filter]
mesh_tally.scores = ['flux', 'fission', 'nu-fission']

# Let's also create a tally to get the flux energy spectrum. We start by
# creating an energy filter
e_min, e_max = 1e-5, 20.0e6
groups = 500
energies = np.logspace(log10(e_min), log10(e_max), groups + 1)
energy_filter = openmc.EnergyFilter(energies)

spectrum_tally = openmc.Tally(name="Flux spectrum")
spectrum_tally.filters = [energy_filter]
spectrum_tally.scores = ['flux']

# Instantiate a Tallies collection and export to XML
tallies = openmc.Tallies([mesh_tally, spectrum_tally])
tallies.export_to_xml()

openmc.run()
universe = openmc.Universe(cells=objects)
image = universe.plot(width=(box_size + 0.3, box_size + 0.3), origin=(box_size/2, box_size/2, 0.1))
image.write_png(fname='plot.png')