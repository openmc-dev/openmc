"""
LCT006-14: Tank-Type Critical Assembly (TCA) 15x15 rods with 2.293 cm pitch; Hc=90.75 cm
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# U(2.596)O2
mat1 = openmc.Material(material_id=1, name="U(2.596)O2")
mat1.set_density("atom/b-cm", 7.035819e-02)
mat1.add_nuclide("U234", 4.8872e-6)
mat1.add_nuclide("U235", 6.0830e-4)
mat1.add_nuclide("U238", 2.2531e-2)
mat1.add_nuclide("O16", 4.7214e-2)

# Aluminum clad
mat2 = openmc.Material(material_id=2, name="Aluminum clad")
mat2.set_density("atom/b-cm", 5.513700e-02)
mat2.add_element("Al", 5.5137e-2)

# Water
mat3 = openmc.Material(material_id=3, name="Water")
mat3.set_density("atom/b-cm", 1.001030e-01)
mat3.add_element("H", 6.6735e-2)
mat3.add_nuclide("O16", 3.3368e-2)
mat3.add_s_alpha_beta("c_H_in_H2O")

materials = openmc.Materials([mat1, mat2, mat3])
materials.export_to_xml()

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
surf4 = openmc.model.RectangularParallelepiped(
    -17.1975, 17.1975, -17.1975, 17.1975, -30.0, 144.15, surface_id=4)

# DX=DY=15*2.293 + 2*30 = 94.395 (refl planar bdy)
surf5 = openmc.model.RectangularParallelepiped(
    -47.1975, 47.1975, -47.1975, 47.1975, -30.0, 144.15, surface_id=5)


# ==============================================================================
# Universes (from COG define unit blocks)
# ==============================================================================

# Unit 1: Unit cell with one fuel rod
u1_cell0 = openmc.Cell(fill=mat1, name="UO2")
u1_cell0.region = -surf1 & -surf2
u1_cell1 = openmc.Cell(fill=mat2, name="Al")
u1_cell1.region = +surf1 & -surf2
u1_cell2 = openmc.Cell(fill=mat3, name="H2O")
u1_cell2.region = +surf2 & -surf3 & -surf5
universe1 = openmc.Universe(universe_id=1, cells=[u1_cell0, u1_cell1, u1_cell2])

# Unit 2: Array of 15x15 fuel rods with 2.293 cm pitch
# Lattice for unit 2: 15x15 array
# Pitch: (2.293000, 2.293000) cm
# Lower left: (-17.1975, -17.1975)

# Lattice filled with universe1
lattice2 = openmc.RectLattice(lattice_id=2)
lattice2.lower_left = [-17.1975, -17.1975]
lattice2.pitch = [2.293000, 2.293000]
lattice2.universes = [[universe1]*15]*15
universe2 = lattice2  # Lattice can be used as universe

# Cell using unit 2: core
cell0 = openmc.Cell(cell_id=0, fill=universe2, name="core")
cell0.region = -surf4 & -surf5

# Cell: refl
cell1 = openmc.Cell(cell_id=1, fill=mat3, name="refl")
cell1.region = -surf3 & +surf4 & -surf5

# Cell: H2O
cell2 = openmc.Cell(cell_id=2, fill=mat3, name="H2O")
cell2.region = +surf2 & -surf3 & -surf5

# ==============================================================================
# Boundary Conditions
# ==============================================================================

# Create outer bounding box with vacuum boundary
# TODO: Adjust dimensions to encompass your entire geometry
boundary_box = openmc.model.RectangularParallelepiped(
    -200, 200, -200, 200, -200, 200,  # xmin, xmax, ymin, ymax, zmin, zmax
    boundary_type="vacuum")

# Create outer void cell (everything outside geometry but inside boundary)
# Particles are killed at the vacuum boundary
outer_region = -boundary_box
outer_region = outer_region & ~cell0.region
outer_region = outer_region & ~cell1.region
outer_region = outer_region & ~cell2.region
outer_cell = openmc.Cell(cell_id=3, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, outer_cell])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# ==============================================================================
# Settings
# ==============================================================================

settings = openmc.Settings()
settings.particles = 15000
settings.batches = 4400
settings.inactive = 100
settings.run_mode = "eigenvalue"

# Source definition
source = openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 45.375))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Watt(a=0.988e6, b=2.249e-6)
settings.source = source

# Enable delayed neutron kinetics and alpha eigenvalue calculations
settings.calculate_prompt_k = True
settings.calculate_alpha = True

settings.export_to_xml()

# ==============================================================================
# Tallies
# ==============================================================================

tallies = openmc.Tallies()
tallies.export_to_xml()

# ==============================================================================
# Run OpenMC
# ==============================================================================

openmc.run()
