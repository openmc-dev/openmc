import openmc


###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 20
inactive = 10
particles = 10000


###############################################################################
#                 Exporting to OpenMC materials.xml File
###############################################################################

# Instantiate some Nuclides
h1 = openmc.Nuclide('H-1')
o16 = openmc.Nuclide('O-16')
u235 = openmc.Nuclide('U-235')

# Instantiate some Materials and register the appropriate Nuclides
fuel = openmc.Material(material_id=1, name='fuel')
fuel.set_density('g/cc', 4.5)
fuel.add_nuclide(u235, 1.)

moderator = openmc.Material(material_id=2, name='moderator')
moderator.set_density('g/cc', 1.0)
moderator.add_nuclide(h1, 2.)
moderator.add_nuclide(o16, 1.)
moderator.add_s_alpha_beta('HH2O', '71t')

# Instantiate a MaterialsFile, register all Materials, and export to XML
materials_file = openmc.MaterialsFile()
materials_file.set_default_xs('71c')
materials_file.add_materials([moderator, fuel])
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# Instantiate Surfaces
left = openmc.XPlane(surface_id=1, x0=-2, name='left')
right = openmc.XPlane(surface_id=2, x0=2, name='right')
bottom = openmc.YPlane(surface_id=3, y0=-2, name='bottom')
top = openmc.YPlane(surface_id=4, y0=2, name='top')
fuel1 = openmc.ZCylinder(surface_id=5, x0=0, y0=0, R=0.4)
fuel2 = openmc.ZCylinder(surface_id=6, x0=0, y0=0, R=0.3)
fuel3 = openmc.ZCylinder(surface_id=7, x0=0, y0=0, R=0.2)

left.set_boundary_type('vacuum')
right.set_boundary_type('vacuum')
top.set_boundary_type('vacuum')
bottom.set_boundary_type('vacuum')

# Instantiate Cells
cell1 = openmc.Cell(cell_id=1, name='Cell 1')
cell2 = openmc.Cell(cell_id=2, name='Cell 2')
cell3 = openmc.Cell(cell_id=101, name='cell 3')
cell4 = openmc.Cell(cell_id=102, name='cell 4')
cell5 = openmc.Cell(cell_id=201, name='cell 5')
cell6 = openmc.Cell(cell_id=202, name='cell 6')
cell7 = openmc.Cell(cell_id=301, name='cell 7')
cell8 = openmc.Cell(cell_id=302, name='cell 8')

# Register Surfaces with Cells
cell1.add_surface(left, halfspace=+1)
cell1.add_surface(right, halfspace=-1)
cell1.add_surface(bottom, halfspace=+1)
cell1.add_surface(top, halfspace=-1)
cell2.add_surface(left, halfspace=+1)
cell2.add_surface(right, halfspace=-1)
cell2.add_surface(bottom, halfspace=+1)
cell2.add_surface(top, halfspace=-1)
cell3.add_surface(fuel1, halfspace=-1)
cell4.add_surface(fuel1, halfspace=+1)
cell5.add_surface(fuel2, halfspace=-1)
cell6.add_surface(fuel2, halfspace=+1)
cell7.add_surface(fuel3, halfspace=-1)
cell8.add_surface(fuel3, halfspace=+1)

# Register Materials with Cells
cell3.set_fill(fuel)
cell4.set_fill(moderator)
cell5.set_fill(fuel)
cell6.set_fill(moderator)
cell7.set_fill(fuel)
cell8.set_fill(moderator)

# Instantiate Universe
univ1 = openmc.Universe(universe_id=1)
univ2 = openmc.Universe(universe_id=2)
univ3 = openmc.Universe(universe_id=3)
univ4 = openmc.Universe(universe_id=5)
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cells with Universe
univ1.add_cells([cell3, cell4])
univ2.add_cells([cell5, cell6])
univ3.add_cells([cell7, cell8])
root.add_cell(cell1)
univ4.add_cell(cell2)

# Instantiate nested Lattices
lattice1 = openmc.Lattice(lattice_id=4, name='4x4 assembly')
lattice1.set_dimension([2, 2])
lattice1.set_lower_left([-1., -1.])
lattice1.set_width([1., 1.])
lattice1.set_universes([[univ1, univ2],
                       [univ2, univ3]])

lattice2 = openmc.Lattice(lattice_id=6, name='4x4 core')
lattice2.set_dimension([2, 2])
lattice2.set_lower_left([-2., -2.])
lattice2.set_width([2., 2.])
lattice2.set_universes([[univ4, univ4],
                       [univ4, univ4]])

# Fill Cell with the Lattice
cell1.set_fill(lattice2)
cell2.set_fill(lattice1)

# Instantiate a Geometry and register the root Universe
geometry = openmc.Geometry()
geometry.set_root_universe(root)

# Instantiate a GeometryFile, register Geometry, and export to XML
geometry_file = openmc.GeometryFile()
geometry_file.set_geometry(geometry)
geometry_file.export_to_xml()


###############################################################################
#                   Exporting to OpenMC settings.xml File
###############################################################################

# Instantiate a SettingsFile, set all runtime parameters, and export to XML
settings_file = openmc.SettingsFile()
settings_file.set_batches(batches)
settings_file.set_inactive(inactive)
settings_file.set_particles(particles)
settings_file.set_source_space('box', [-1, -1, -1, 1, 1, 1])
settings_file.export_to_xml()


###############################################################################
#                   Exporting to OpenMC plots.xml File
###############################################################################

plot = openmc.Plot(plot_id=1)
plot.set_origin([0, 0, 0])
plot.set_width([4, 4])
plot.set_pixels([400, 400])
plot.set_color('mat')

# Instantiate a PlotsFile, add Plot, and export to XML
plot_file = openmc.PlotsFile()
plot_file.add_plot(plot)
plot_file.export_to_xml()


###############################################################################
#                   Exporting to OpenMC tallies.xml File
###############################################################################

# Instantiate a tally mesh
mesh = openmc.Mesh(mesh_id=1)
mesh.set_type('rectangular')
mesh.set_dimension([4, 4])
mesh.set_lower_left([-2, -2])
mesh.set_width([1, 1])

# Instantiate tally Filter
mesh_filter = openmc.Filter()
mesh_filter.set_mesh(mesh)

# Instantiate the Tally
tally = openmc.Tally(tally_id=1)
tally.add_filter(mesh_filter)
tally.add_score('total')

# Instantiate a TalliesFile, register Tally/Mesh, and export to XML
tallies_file = openmc.TalliesFile()
tallies_file.add_mesh(mesh)
tallies_file.add_tally(tally)
tallies_file.export_to_xml()
