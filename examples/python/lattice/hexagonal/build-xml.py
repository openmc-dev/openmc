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
fe56 = openmc.Nuclide('Fe-56')

# Instantiate some Materials and register the appropriate Nuclides
fuel = openmc.Material(material_id=1, name='fuel')
fuel.set_density('g/cc', 4.5)
fuel.add_nuclide(u235, 1.)

moderator = openmc.Material(material_id=2, name='moderator')
moderator.set_density('g/cc', 1.0)
moderator.add_nuclide(h1, 2.)
moderator.add_nuclide(o16, 1.)
moderator.add_s_alpha_beta('HH2O', '71t')

iron = openmc.Material(material_id=3, name='iron')
iron.set_density('g/cc', 7.9)
iron.add_nuclide(fe56, 1.)

# Instantiate a MaterialsFile, register all Materials, and export to XML
materials_file = openmc.MaterialsFile()
materials_file.set_default_xs('71c')
materials_file.add_materials([moderator, fuel, iron])
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# Instantiate Surfaces
left = openmc.XPlane(surface_id=1, x0=-3, name='left')
right = openmc.XPlane(surface_id=2, x0=3, name='right')
bottom = openmc.YPlane(surface_id=3, y0=-4, name='bottom')
top = openmc.YPlane(surface_id=4, y0=4, name='top')
fuel_surf = openmc.ZCylinder(surface_id=5, x0=0, y0=0, R=0.4)

left.set_boundary_type('vacuum')
right.set_boundary_type('vacuum')
top.set_boundary_type('vacuum')
bottom.set_boundary_type('vacuum')

# Instantiate Cells
cell1 = openmc.Cell(cell_id=1, name='Cell 1')
cell2 = openmc.Cell(cell_id=101, name='cell 2')
cell3 = openmc.Cell(cell_id=102, name='cell 3')
cell4 = openmc.Cell(cell_id=500, name='cell 4')
cell5 = openmc.Cell(cell_id=600, name='cell 5')
cell6 = openmc.Cell(cell_id=601, name='cell 6')

# Register Surfaces with Cells
cell1.add_surface(left, halfspace=+1)
cell1.add_surface(right, halfspace=-1)
cell1.add_surface(bottom, halfspace=+1)
cell1.add_surface(top, halfspace=-1)
cell2.add_surface(fuel_surf, halfspace=-1)
cell3.add_surface(fuel_surf, halfspace=+1)
cell5.add_surface(fuel_surf, halfspace=-1)
cell6.add_surface(fuel_surf, halfspace=+1)

# Register Materials with Cells
cell2.set_fill(fuel)
cell3.set_fill(moderator)
cell4.set_fill(moderator)
cell5.set_fill(iron)
cell6.set_fill(moderator)

# Instantiate Universe
univ1 = openmc.Universe(universe_id=1)
univ2 = openmc.Universe(universe_id=3)
univ3 = openmc.Universe(universe_id=4)
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cells with Universe
univ1.add_cells([cell2, cell3])
univ2.add_cells([cell4])
univ3.add_cells([cell5, cell6])
root.add_cell(cell1)

# Instantiate a Lattice
lattice = openmc.HexLattice(lattice_id=5)
lattice.set_center([0., 0., 0.])
lattice.set_pitch([1., 2.])
lattice.set_universes([
     [ [univ2] + [univ3]*11, [univ2] + [univ3]*5, [univ3] ],
     [ [univ2] + [univ1]*11, [univ2] + [univ1]*5, [univ1] ],
     [ [univ2] + [univ3]*11, [univ2] + [univ3]*5, [univ3] ]])
lattice.set_outer(univ2)

# Fill Cell with the Lattice
cell1.set_fill(lattice)

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

plot_xy = openmc.Plot(plot_id=1)
plot_xy.set_filename('plot_xy')
plot_xy.set_origin([0, 0, 0])
plot_xy.set_width([6, 6])
plot_xy.set_pixels([400, 400])
plot_xy.set_color('mat')

plot_yz = openmc.Plot(plot_id=2)
plot_yz.set_filename('plot_yz')
plot_yz.set_basis('yz')
plot_yz.set_origin([0, 0, 0])
plot_yz.set_width([8, 8])
plot_yz.set_pixels([400, 400])
plot_yz.set_color('mat')

# Instantiate a PlotsFile, add Plot, and export to XML
plot_file = openmc.PlotsFile()
plot_file.add_plot(plot_xy)
plot_file.add_plot(plot_yz)
plot_file.export_to_xml()
