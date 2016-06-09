import numpy as np
import openmc

###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 15
inactive = 5
particles = 10000


###############################################################################
#                 Exporting to OpenMC materials.xml File
###############################################################################

# Instantiate some Nuclides
h1 = openmc.Nuclide('H-1')
o16 = openmc.Nuclide('O-16')
u235 = openmc.Nuclide('U-235')
u238 = openmc.Nuclide('U-238')

# Instantiate some Materials and register the appropriate Nuclides
fuel1 = openmc.Material(material_id=1, name='fuel')
fuel1.set_density('g/cc', 4.5)
fuel1.add_nuclide(u235, 1.)

fuel2 = openmc.Material(material_id=2, name='depleted fuel')
fuel2.set_density('g/cc', 4.5)
fuel2.add_nuclide(u238, 1.)

moderator = openmc.Material(material_id=3, name='moderator')
moderator.set_density('g/cc', 1.0)
moderator.add_nuclide(h1, 2.)
moderator.add_nuclide(o16, 1.)
moderator.add_s_alpha_beta('HH2O', '71t')

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([fuel1, fuel2, moderator])
materials_file.default_xs = '71c'
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml file
###############################################################################

# Instantiate planar surfaces
x1 = openmc.XPlane(surface_id=1, x0=-10)
x2 = openmc.XPlane(surface_id=2, x0=-7)
x3 = openmc.XPlane(surface_id=3, x0=-4)
x4 = openmc.XPlane(surface_id=4, x0=4)
x5 = openmc.XPlane(surface_id=5, x0=7)
x6 = openmc.XPlane(surface_id=6, x0=10)
y1 = openmc.YPlane(surface_id=11, y0=-10)
y2 = openmc.YPlane(surface_id=12, y0=-7)
y3 = openmc.YPlane(surface_id=13, y0=-4)
y4 = openmc.YPlane(surface_id=14, y0=4)
y5 = openmc.YPlane(surface_id=15, y0=7)
y6 = openmc.YPlane(surface_id=16, y0=10)
z1 = openmc.ZPlane(surface_id=21, z0=-10)
z2 = openmc.ZPlane(surface_id=22, z0=-7)
z3 = openmc.ZPlane(surface_id=23, z0=-4)
z4 = openmc.ZPlane(surface_id=24, z0=4)
z5 = openmc.ZPlane(surface_id=25, z0=7)
z6 = openmc.ZPlane(surface_id=26, z0=10)

# Set vacuum boundary conditions on outside
for surface in [x1, x6, y1, y6, z1, z6]:
    surface.boundary_type = 'vacuum'

# Instantiate Cells
inner_box = openmc.Cell(cell_id=1, name='inner box')
middle_box = openmc.Cell(cell_id=2, name='middle box')
outer_box = openmc.Cell(cell_id=3, name='outer box')

# Use each set of six planes to create solid cube regions. We can then use these
# to create cubic shells.
inner_cube = +x3 & -x4 & +y3 & -y4 & +z3 & -z4
middle_cube = +x2 & -x5 & +y2 & -y5 & +z2 & -z5
outer_cube = +x1 & -x6 & +y1 & -y6 & +z1 & -z6
outside_inner_cube = -x3 | +x4 | -y3 | +y4 | -z3 | +z4

# Use surface half-spaces to define regions
inner_box.region = inner_cube
middle_box.region = middle_cube & outside_inner_cube
outer_box.region = outer_cube & ~middle_cube

# Register Materials with Cells
inner_box.fill = fuel1
middle_box.fill = fuel2
outer_box.fill = moderator

# Instantiate root universe
root = openmc.Universe(universe_id=0, name='root universe')
root.add_cells([inner_box, middle_box, outer_box])

# Instantiate a Geometry, register the root Universe, and export to XML
geometry = openmc.Geometry(root)
geometry.export_to_xml()


###############################################################################
#                   Exporting to OpenMC settings.xml File
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles

# Create an initial uniform spatial source distribution over fissionable zones
uniform_dist = openmc.stats.Box(*outer_cube.bounding_box, only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

settings_file.export_to_xml()

###############################################################################
#                   Exporting to OpenMC plots.xml File
###############################################################################

plot = openmc.Plot(plot_id=1)
plot.origin = [0, 0, 0]
plot.width = [20, 20]
plot.pixels = [200, 200]
plot.color = 'cell'

# Instantiate a Plots collection and export to XML
plot_file = openmc.Plots([plot])
plot_file.export_to_xml()
