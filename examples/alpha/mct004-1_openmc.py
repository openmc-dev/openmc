"""
MIX-COMP-THERM-004-1: 23x23 MOX pins with 1.825 cm pitch; Hc=59.55
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# MOX
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("Pu238", 2.000300e-06)
mat1.add_nuclide("Pu239", 2.749100e-04)
mat1.add_nuclide("Pu240", 8.841700e-05)
mat1.add_nuclide("Pu241", 2.792300e-05)
mat1.add_nuclide("Pu242", 8.123400e-06)
mat1.add_nuclide("Am241", 1.353100e-06)
mat1.add_nuclide("U234", 7.174900e-07)
mat1.add_nuclide("U235", 9.392600e-05)
mat1.add_nuclide("U238", 1.295100e-02)
mat1.add_nuclide("O16", 2.783700e-02)
mat1.add_nuclide("B10", 6.041800e-08)
mat1.add_nuclide("B11", 2.431900e-07)

# Cladding (with air gap)
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_element("Zr", 3.777200e-02)
mat2.add_element("Sn", 4.373700e-04)
mat2.add_element("Fe", 8.857000e-05)
mat2.add_element("Cr", 6.611900e-05)
mat2.add_element("Ni", 3.586400e-05)

# Water
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_nuclide("H1", 6.673500e-02)
mat3.add_nuclide("O16", 3.336800e-02)
mat3.add_s_alpha_beta("c_H_in_H2O")

# Aluminum
mat4 = openmc.Material(material_id=4)
mat4.set_density("sum")
mat4.add_element("Al", 6.022400e-02)

# SS304L
mat5 = openmc.Material(material_id=5)
mat5.set_density("sum")
mat5.add_element("C", 1.192800e-04)
mat5.add_element("Si", 1.700300e-03)
mat5.add_element("Mn", 1.738500e-03)
mat5.add_element("P", 6.938100e-05)
mat5.add_element("S", 4.467300e-05)
mat5.add_element("Ni", 8.950600e-03)
mat5.add_element("Cr", 1.745000e-02)
mat5.add_element("Fe", 5.720200e-02)

# Ordinary
mat6 = openmc.Material(material_id=6)
mat6.set_density("sum")
mat6.add_nuclide("H1", 1.374200e-02)
mat6.add_nuclide("O16", 4.591900e-02)
mat6.add_element("C", 1.153200e-04)
mat6.add_element("Na", 9.639500e-04)
mat6.add_element("Mg", 1.238800e-04)
mat6.add_element("Al", 1.740900e-03)
mat6.add_element("Si", 1.661700e-02)
mat6.add_element("K", 4.605200e-04)
mat6.add_element("Ca", 1.502500e-03)
mat6.add_element("Fe", 3.449200e-04)
mat6.add_s_alpha_beta("c_H_in_H2O")

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6])

# ==============================================================================
# Geometry
# ==============================================================================

# Fuel
surf1 = openmc.ZCylinder(surface_id=1, r=0.5325)
# Aluminum
surf2 = openmc.ZCylinder(surface_id=2, r=0.5325)
# Void
surf3 = openmc.ZCylinder(surface_id=3, r=0.5325)
# Cladding
surf4 = openmc.ZCylinder(surface_id=4, r=0.6115)
# Lattice plane boundaries
surf5 = openmc.model.RectangularParallelepiped(-20.9875, 20.9875, -20.9875, 20.9875, -499.995, 499.995)
# BCD
surf6 = openmc.model.RectangularParallelepiped(-50.9875, 50.9875, -50.9875, 50.9875, -71.60000000000001, 81.17, boundary_type="vacuum")
# Middle grid, bottom
surf10 = openmc.ZPlane(surface_id=10, z0=80.57)
# Water height
surf11 = openmc.ZPlane(surface_id=11, z0=59.55)
# Al lower grid, top
surf12 = openmc.ZPlane(surface_id=12, z0=-11.785)
# Al lower grid, bottom
surf13 = openmc.ZPlane(surface_id=13, z0=-12.385)
# Al support plate, top
surf14 = openmc.ZPlane(surface_id=14, z0=-16.83)
# Al support plate, bottom; SS304L, top
surf15 = openmc.ZPlane(surface_id=15, z0=-18.10)
# SS304L support plate, bottom
surf16 = openmc.ZPlane(surface_id=16, z0=-20.30)
# SS304L liner, top
surf17 = openmc.ZPlane(surface_id=17, z0=-34.10)
# SS304L liner, bottom
surf18 = openmc.ZPlane(surface_id=18, z0=-34.60)

# ------------------------------------------------------------------------------
# Universes
# ------------------------------------------------------------------------------

u1_cell0 = openmc.Cell(fill=mat1)
u1_cell0.region = -surf1 & -surf4
u1_cell1 = openmc.Cell(fill=mat4)
u1_cell1.region = +surf1 & -surf2 & -surf4
u1_cell2 = openmc.Cell(fill=mat2)
u1_cell2.region = +surf1 & +surf2 & +surf3 & -surf4
universe1 = openmc.Universe(universe_id=1, cells=[u1_cell0, u1_cell1, u1_cell2])

u2_cell0 = openmc.Cell(fill=mat4)
u2_cell0.region = -surf6 & +surf10
u2_cell1 = openmc.Cell(fill=mat3)
u2_cell1.region = -surf6 & -surf11 & +surf12
u2_cell2 = openmc.Cell(fill=mat4)
u2_cell2.region = -surf6 & -surf12 & +surf13
u2_cell3 = openmc.Cell(fill=mat3)
u2_cell3.region = -surf6 & -surf13 & +surf14
u2_cell4 = openmc.Cell(fill=mat4)
u2_cell4.region = -surf6 & -surf14 & +surf15
u2_cell5 = openmc.Cell(fill=mat5)
u2_cell5.region = -surf6 & -surf15 & +surf16
u2_cell6 = openmc.Cell(fill=mat3)
u2_cell6.region = -surf6 & -surf16 & +surf17
u2_cell7 = openmc.Cell(fill=mat5)
u2_cell7.region = -surf6 & -surf17 & +surf18
u2_cell8 = openmc.Cell(fill=mat6)
u2_cell8.region = -surf6 & -surf18
universe2 = openmc.Universe(universe_id=2, cells=[u2_cell0, u2_cell1, u2_cell2, u2_cell3, u2_cell4, u2_cell5, u2_cell6, u2_cell7, u2_cell8])

universe3 = openmc.Universe(universe_id=3, cells=[])

# Lattice 4: 23x23 array
lattice4 = openmc.RectLattice(lattice_id=4)
lattice4.lower_left = [-20.9875, -20.9875]
lattice4.pitch = [1.825000, 1.825000]
lattice4.universes = [
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
    [universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3, universe3],
]
universe4 = openmc.Universe(universe_id=4)
universe4.add_cell(openmc.Cell(fill=lattice4))

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# lttc
cell1 = openmc.Cell(cell_id=1, fill=universe4)
cell1.region = -surf5 & -surf6

# stuff
cell2 = openmc.Cell(cell_id=2, fill=universe2)
cell2.region = +surf5 & -surf6

# Clad
cell6 = openmc.Cell(cell_id=6, fill=mat2)
cell6.region = +surf1 & +surf2 & +surf3 & -surf4

# Conc
cell16 = openmc.Cell(cell_id=16, fill=mat6)
cell16.region = -surf6 & -surf18

# Alles
cell17 = openmc.Cell(cell_id=17, fill=universe2)
cell17.region = +surf4 & -surf6

root_universe = openmc.Universe(cells=[cell1, cell2, cell6, cell16, cell17])
geometry = openmc.Geometry(root_universe)

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 15000
settings.batches = 4400
settings.inactive = 100
settings.run_mode = "eigenvalue"

source = openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 29.775))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
