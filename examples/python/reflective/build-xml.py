import numpy as np
import openmc

###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 500
inactive = 10
particles = 10000


###############################################################################
#                 Exporting to OpenMC materials.xml file
###############################################################################

# Instantiate a Nuclides
u235 = openmc.Nuclide('U-235')

# Instantiate a Material and register the Nuclide
fuel = openmc.Material(material_id=1, name='fuel')
fuel.set_density('g/cc', 4.5)
fuel.add_nuclide(u235, 1.)

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([fuel])
materials_file.default_xs = '71c'
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml file
###############################################################################

# Instantiate Surfaces
surf1 = openmc.XPlane(surface_id=1, x0=-1, name='surf 1')
surf2 = openmc.XPlane(surface_id=2, x0=+1, name='surf 2')
surf3 = openmc.YPlane(surface_id=3, y0=-1, name='surf 3')
surf4 = openmc.YPlane(surface_id=4, y0=+1, name='surf 4')
surf5 = openmc.ZPlane(surface_id=5, z0=-1, name='surf 5')
surf6 = openmc.ZPlane(surface_id=6, z0=+1, name='surf 6')

surf1.boundary_type = 'vacuum'
surf2.boundary_type = 'vacuum'
surf3.boundary_type = 'reflective'
surf4.boundary_type = 'reflective'
surf5.boundary_type = 'reflective'
surf6.boundary_type = 'reflective'

# Instantiate Cell
cell = openmc.Cell(cell_id=1, name='cell 1')

# Use surface half-spaces to define region
cell.region = +surf1 & -surf2 & +surf3 & -surf4 & +surf5 & -surf6

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
uniform_dist = openmc.stats.Box(*cell.region.bounding_box,
                                only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

settings_file.export_to_xml()
