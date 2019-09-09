import numpy as np
import openmc

###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 500
inactive = 100
particles = 10000


###############################################################################
#                 Exporting to OpenMC materials.xml file
###############################################################################

# Instantiate a Material and register the Nuclide
fuel = openmc.Material(material_id=1, name='fuel')
fuel.set_density('g/cc', 4.5)
fuel.add_nuclide('U235', 1.)

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([fuel])
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml file
###############################################################################

# Instantiate Surfaces
surf1 = openmc.ZPlane(surface_id=1, z0=0.0, name='surf 1')
surf2 = openmc.Sphere(surface_id=2, x0=0.0, y0=0.0, z0=0.0, r=0.39218, name='surf 2')

surf1.boundary_type = 'vacuum'
surf2.boundary_type = 'white'

# Instantiate Cell
cell = openmc.Cell(cell_id=1, name='cell 1')

# Use surface half-spaces to define region
cell.region = -surf1 & -surf2 

# Register Material with Cell
cell.fill = fuel

# Instantiate Universes
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cell with Universe
root.add_cell(cell)

# Instantiate a Geometry, register the root Universe, and export to XML
geometry = openmc.Geometry(root)
geometry.export_to_xml()


###############################################################################
#                   Exporting to OpenMC settings.xml file
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles

# Create an initial uniform spatial source distribution over fissionable zones
uniform_dist = openmc.stats.Point(xyz=(0, 0, -0.1))
settings_file.source = openmc.source.Source(space=uniform_dist)

settings_file.export_to_xml()
