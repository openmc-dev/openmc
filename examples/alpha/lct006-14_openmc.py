"""
LCT006-14: Tank-Type Critical Assembly (TCA) 15x15 rods with 2.293 cm pitch; Hc=90.75 cm
Converted from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# U(2.596)O2
mat1 = openmc.Material(material_id=1)
mat1.set_density("sum")
mat1.add_nuclide("U234", 4.887200e-06)
mat1.add_nuclide("U235", 6.083000e-04)
mat1.add_nuclide("U238", 2.253100e-02)
mat1.add_nuclide("O16", 4.721400e-02)

# Aluminum clad
mat2 = openmc.Material(material_id=2)
mat2.set_density("sum")
mat2.add_element("Al", 5.513700e-02)

# Water
mat3 = openmc.Material(material_id=3)
mat3.set_density("sum")
mat3.add_nuclide("H1", 6.673500e-02)
mat3.add_nuclide("O16", 3.336800e-02)
mat3.add_s_alpha_beta("c_H_in_H2O")

materials = openmc.Materials([mat1, mat2, mat3])

# ==============================================================================
# Geometry
# ==============================================================================

# Fuel
surf1 = openmc.ZCylinder(surface_id=1, r=0.625)
# Clad
surf2 = openmc.ZCylinder(surface_id=2, r=0.7085)
# Water critical height, Hc
surf3 = openmc.ZPlane(surface_id=3, z0=90.75)
# DX=DY=15*2.293        = 34.395 (core planar bdy)
surf4 = openmc.model.RectangularParallelepiped(-17.1975, 17.1975, -17.1975, 17.1975, -30.0, 144.15)
# DX=DY=15*2.293 + 2*30 = 94.395 (refl planar bdy)
surf5 = openmc.model.RectangularParallelepiped(-47.1975, 47.1975, -47.1975, 47.1975, -30.0, 144.15, boundary_type="vacuum")

# ------------------------------------------------------------------------------
# Universes
# ------------------------------------------------------------------------------

u1_cell0 = openmc.Cell(fill=mat1)
u1_cell0.region = -surf1 & -surf2
u1_cell1 = openmc.Cell(fill=mat2)
u1_cell1.region = +surf1 & -surf2
u1_cell2 = openmc.Cell(fill=mat3)
u1_cell2.region = +surf2 & -surf3 & -surf5
universe1 = openmc.Universe(universe_id=1, cells=[u1_cell0, u1_cell1, u1_cell2])

# Lattice 2: 15x15 array
lattice2 = openmc.RectLattice(lattice_id=2)
lattice2.lower_left = [-17.1975, -17.1975]
lattice2.pitch = [2.293000, 2.293000]
lattice2.universes = [
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
    [universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1, universe1],
]
universe2 = openmc.Universe(universe_id=2)
universe2.add_cell(openmc.Cell(fill=lattice2))

# ------------------------------------------------------------------------------
# Root Cells
# ------------------------------------------------------------------------------

# core
cell1 = openmc.Cell(cell_id=1, fill=universe2)
cell1.region = -surf4 & -surf5

# refl
cell2 = openmc.Cell(cell_id=2, fill=mat3)
cell2.region = -surf3 & +surf4 & -surf5

# H2O
cell6 = openmc.Cell(cell_id=6, fill=mat3)
cell6.region = +surf2 & -surf3 & -surf5

root_universe = openmc.Universe(cells=[cell1, cell2, cell6])
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
source.space = openmc.stats.Point((0.0, 0.0, 45.375))
settings.source = source

# ==============================================================================
# Export and Run
# ==============================================================================

materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
