import numpy as np
import openmc
import openmc.mgxs

###############################################################################
# Create multigroup data

# Instantiate the energy group data
group_edges = [1e-5, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6]
groups = openmc.mgxs.EnergyGroups(group_edges)

# Instantiate the 7-group (C5G7) cross section data
uo2_xsdata = openmc.XSdata('UO2', groups)
uo2_xsdata.order = 0
uo2_xsdata.set_total(
    [0.1779492, 0.3298048, 0.4803882, 0.5543674, 0.3118013, 0.3951678,
     0.5644058])
uo2_xsdata.set_absorption([8.0248e-03, 3.7174e-03, 2.6769e-02, 9.6236e-02,
                           3.0020e-02, 1.1126e-01, 2.8278e-01])
scatter_matrix = np.array(
    [[[0.1275370, 0.0423780, 0.0000094, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
      [0.0000000, 0.3244560, 0.0016314, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
      [0.0000000, 0.0000000, 0.4509400, 0.0026792, 0.0000000, 0.0000000, 0.0000000],
      [0.0000000, 0.0000000, 0.0000000, 0.4525650, 0.0055664, 0.0000000, 0.0000000],
      [0.0000000, 0.0000000, 0.0000000, 0.0001253, 0.2714010, 0.0102550, 0.0000000],
      [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0012968, 0.2658020, 0.0168090],
      [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0085458, 0.2730800]]])
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
uo2_xsdata.set_scatter_matrix(scatter_matrix)
uo2_xsdata.set_fission([7.21206e-03, 8.19301e-04, 6.45320e-03,
                        1.85648e-02, 1.78084e-02, 8.30348e-02,
                        2.16004e-01])
uo2_xsdata.set_nu_fission([2.005998e-02, 2.027303e-03, 1.570599e-02,
                           4.518301e-02, 4.334208e-02, 2.020901e-01,
                           5.257105e-01])
uo2_xsdata.set_chi([5.8791e-01, 4.1176e-01, 3.3906e-04, 1.1761e-07, 0.0000e+00,
                    0.0000e+00, 0.0000e+00])

h2o_xsdata = openmc.XSdata('LWTR', groups)
h2o_xsdata.order = 0
h2o_xsdata.set_total([0.15920605, 0.412969593, 0.59030986, 0.58435,
                      0.718, 1.2544497, 2.650379])
h2o_xsdata.set_absorption([6.0105e-04, 1.5793e-05, 3.3716e-04,
                           1.9406e-03, 5.7416e-03, 1.5001e-02,
                           3.7239e-02])
scatter_matrix = np.array(
    [[[0.0444777, 0.1134000, 0.0007235, 0.0000037, 0.0000001, 0.0000000, 0.0000000],
      [0.0000000, 0.2823340, 0.1299400, 0.0006234, 0.0000480, 0.0000074, 0.0000010],
      [0.0000000, 0.0000000, 0.3452560, 0.2245700, 0.0169990, 0.0026443, 0.0005034],
      [0.0000000, 0.0000000, 0.0000000, 0.0910284, 0.4155100, 0.0637320, 0.0121390],
      [0.0000000, 0.0000000, 0.0000000, 0.0000714, 0.1391380, 0.5118200, 0.0612290],
      [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0022157, 0.6999130, 0.5373200],
      [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1324400, 2.4807000]]])
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
h2o_xsdata.set_scatter_matrix(scatter_matrix)

mg_cross_sections_file = openmc.MGXSLibrary(groups)
mg_cross_sections_file.add_xsdatas([uo2_xsdata, h2o_xsdata])
mg_cross_sections_file.export_to_hdf5()

###############################################################################
# Create materials for the problem

# Instantiate some Materials and register the appropriate macroscopic data
uo2 = openmc.Material(name='UO2 fuel')
uo2.set_density('macro', 1.0)
uo2.add_macroscopic('UO2')

water = openmc.Material(name='Water')
water.set_density('macro', 1.0)
water.add_macroscopic('LWTR')

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([uo2, water])
materials_file.cross_sections = "mgxs.h5"
materials_file.export_to_xml()

###############################################################################
# Define problem geometry

# The geometry we will define a simplified pincell with fuel radius 0.54 cm
# surrounded by moderator (same as in the multigroup example).
# In random ray, we typically want several radial regions and azimuthal
# sectors in both the fuel and moderator areas of the pincell. This is
# due to the flat source approximation requiring that source regions are
# small compared to the typical mean free path of a neutron. Below we
# sudivide the basic pincell into 8 aziumthal sectors (pizza slices) and
# 5 concentric rings in both the fuel and moderator.

# TODO: When available in OpenMC, use cylindrical lattice instead to
# simplify definition and improve runtime performance.

pincell_base = openmc.Universe()

# These are the subdivided radii (creating 5 concentric regions in the
# fuel and moderator)
ring_radii = [0.241, 0.341, 0.418, 0.482, 0.54, 0.572, 0.612, 0.694, 0.786]
fills = [uo2, uo2, uo2, uo2, uo2, water, water, water, water, water]

# We then create cells representing the bounded rings, with special
# treatment for both the innermost and outermost cells
cells = []
for r in range(10):
    cell = []
    if r == 0:
        outer_bound = openmc.ZCylinder(r=ring_radii[r])
        cell = openmc.Cell(fill=fills[r], region=-outer_bound)
    elif r == 9:
        inner_bound = openmc.ZCylinder(r=ring_radii[r-1])
        cell = openmc.Cell(fill=fills[r], region=+inner_bound)
    else:
        inner_bound = openmc.ZCylinder(r=ring_radii[r-1])
        outer_bound = openmc.ZCylinder(r=ring_radii[r])
        cell = openmc.Cell(fill=fills[r], region=+inner_bound & -outer_bound)
    pincell_base.add_cell(cell)

# We then generate 8 planes to bound 8 azimuthal sectors
azimuthal_planes = []
for i in range(8):
    angle = 2 * i * openmc.pi / 8
    normal_vector = (-openmc.sin(angle), openmc.cos(angle), 0)
    azimuthal_planes.append(openmc.Plane(a=normal_vector[0], b=normal_vector[1], c=normal_vector[2], d=0))

# Create a cell for each azimuthal sector using the pincell base class
azimuthal_cells = []
for i in range(8):
    azimuthal_cell = openmc.Cell(name=f'azimuthal_cell_{i}')
    azimuthal_cell.fill = pincell_base
    azimuthal_cell.region = +azimuthal_planes[i] & -azimuthal_planes[(i+1) % 8]
    azimuthal_cells.append(azimuthal_cell)

# Create the (subdivided) geometry with the azimuthal universes
pincell = openmc.Universe(cells=azimuthal_cells)

# Create a region represented as the inside of a rectangular prism
pitch = 1.26
box = openmc.model.RectangularPrism(pitch, pitch, boundary_type='reflective')
pincell_bounded = openmc.Cell(fill=pincell, region=-box, name='pincell')

# Create a geometry (specifying merge surfaces option to remove
# all the redundant cylinder/plane surfaces) and export to XML
geometry = openmc.Geometry([pincell_bounded], merge_surfaces=True)
geometry.export_to_xml()

###############################################################################
# Define problem settings

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings = openmc.Settings()
settings.energy_mode = "multi-group"
settings.batches = 600
settings.inactive = 300
settings.particles = 50

# Create an initial uniform spatial source distribution for sampling rays.
# Note that this must be uniform in space and angle.
lower_left = (-pitch/2, -pitch/2, -1)
upper_right = (pitch/2, pitch/2, 1)
uniform_dist = openmc.stats.Box(lower_left, upper_right)
settings.random_ray['ray_source'] = openmc.IndependentSource(space=uniform_dist)
settings.random_ray['distance_inactive'] = 40.0
settings.random_ray['distance_active'] = 400.0

settings.export_to_xml()

###############################################################################
# Define tallies

# Create a mesh that will be used for tallying
mesh = openmc.RegularMesh()
mesh.dimension = (2, 2)
mesh.lower_left = (-pitch/2, -pitch/2)
mesh.upper_right = (pitch/2, pitch/2)

# Create a mesh filter that can be used in a tally
mesh_filter = openmc.MeshFilter(mesh)

# Let's also create a filter to measure each group
# indepdendently
energy_filter = openmc.EnergyFilter(group_edges)

# Now use the mesh filter in a tally and indicate what scores are desired
tally = openmc.Tally(name="Mesh and Energy tally")
tally.filters = [mesh_filter, energy_filter]
tally.scores = ['flux', 'fission', 'nu-fission']

# Instantiate a Tallies collection and export to XML
tallies = openmc.Tallies([tally])
tallies.export_to_xml()

###############################################################################
#                   Exporting to OpenMC plots.xml file
###############################################################################

plot = openmc.VoxelPlot()
plot.origin = [0, 0, 0]
plot.width = [pitch, pitch, pitch]
plot.pixels = [1000, 1000, 1]

# Instantiate a Plots collection and export to XML
plots = openmc.Plots([plot])
plots.export_to_xml()
