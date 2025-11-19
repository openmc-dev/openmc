"""
LCT006-4: Tank-Type Critical Assembly (TCA) 17x17 rods with 1.956 cm pitch; Hc=114.59 cm
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
surf3 = openmc.ZPlane(surface_id=3, z0=114.59)

# DX=DY=17*1.956        = 33.252 (core planar bdy)
# Box (6 planes): xmin=-16.626, xmax=16.626, ymin=-16.626, ymax=16.626, zmin=-30.0, zmax=144.15
surf4_xmin = openmc.XPlane(surface_id=10000, x0=-16.626)
surf4_xmax = openmc.XPlane(surface_id=10001, x0=16.626)
surf4_ymin = openmc.YPlane(surface_id=10002, y0=-16.626)
surf4_ymax = openmc.YPlane(surface_id=10003, y0=16.626)
surf4_zmin = openmc.ZPlane(surface_id=10004, z0=-30.0)
surf4_zmax = openmc.ZPlane(surface_id=10005, z0=144.15)

# DX=DY=17*1.956 + 2*30 = 93.252 (refl planar bdy)
# Box (6 planes): xmin=-46.626, xmax=46.626, ymin=-46.626, ymax=46.626, zmin=-30.0, zmax=144.15
surf5_xmin = openmc.XPlane(surface_id=10006, x0=-46.626)
surf5_xmax = openmc.XPlane(surface_id=10007, x0=46.626)
surf5_ymin = openmc.YPlane(surface_id=10008, y0=-46.626)
surf5_ymax = openmc.YPlane(surface_id=10009, y0=46.626)
surf5_zmin = openmc.ZPlane(surface_id=10010, z0=-30.0)
surf5_zmax = openmc.ZPlane(surface_id=10011, z0=144.15)


# ==============================================================================
# Universes (from COG define unit blocks)
# ==============================================================================

# Unit 1: Unit cell with one fuel rod
u1_cell0 = openmc.Cell(fill=mat1, name="UO2")
u1_cell0.region = -surf1 & -surf2
u1_cell1 = openmc.Cell(fill=mat2, name="Al")
u1_cell1.region = +surf1 & -surf2
u1_cell2 = openmc.Cell(fill=mat3, name="H2O")
u1_cell2.region = +surf2 & -surf3 & (+surf5_xmin & -surf5_xmax & +surf5_ymin & -surf5_ymax & +surf5_zmin & -surf5_zmax)
universe1 = openmc.Universe(universe_id=1, cells=[u1_cell0, u1_cell1, u1_cell2])

# Unit 2: Array of 17x17 fuel rods with 1.956 cm pitch
# Lattice for unit 2: 17x17 array
# Pitch: (1.956000, 1.956000) cm
# Lower left: (-16.626, -16.626)

# Lattice filled with universe1
lattice2 = openmc.RectLattice(lattice_id=2)
lattice2.lower_left = [-16.626, -16.626]
lattice2.pitch = [1.956000, 1.956000]
lattice2.universes = [[universe1]*17]*17
universe2 = lattice2  # Lattice can be used as universe

# Cell using unit 2: core
cell0 = openmc.Cell(cell_id=0, fill=universe2, name="core")
cell0.region = (+surf4_xmin & -surf4_xmax & +surf4_ymin & -surf4_ymax & +surf4_zmin & -surf4_zmax) & (+surf5_xmin & -surf5_xmax & +surf5_ymin & -surf5_ymax & +surf5_zmin & -surf5_zmax)

# Cell: refl
cell1 = openmc.Cell(cell_id=1, fill=mat3, name="refl")
cell1.region = -surf3 & (-surf4_xmin | +surf4_xmax | -surf4_ymin | +surf4_ymax | -surf4_zmin | +surf4_zmax) & (+surf5_xmin & -surf5_xmax & +surf5_ymin & -surf5_ymax & +surf5_zmin & -surf5_zmax)

# Cell: H2O
cell2 = openmc.Cell(cell_id=2, fill=mat3, name="H2O")
cell2.region = +surf2 & -surf3 & (+surf5_xmin & -surf5_xmax & +surf5_ymin & -surf5_ymax & +surf5_zmin & -surf5_zmax)

# ==============================================================================
# Boundary Conditions
# ==============================================================================

# Create outer bounding box with vacuum boundary (6 planes)
# TODO: Adjust dimensions to encompass your entire geometry
boundary_xmin = openmc.XPlane(surface_id=10012, x0=-200, boundary_type="vacuum")
boundary_xmax = openmc.XPlane(surface_id=10013, x0=200, boundary_type="vacuum")
boundary_ymin = openmc.YPlane(surface_id=10014, y0=-200, boundary_type="vacuum")
boundary_ymax = openmc.YPlane(surface_id=10015, y0=200, boundary_type="vacuum")
boundary_zmin = openmc.ZPlane(surface_id=10016, z0=-200, boundary_type="vacuum")
boundary_zmax = openmc.ZPlane(surface_id=10017, z0=200, boundary_type="vacuum")

# Create outer void cell (everything outside geometry but inside boundary)
# Particles are killed at the vacuum boundary
outer_region = +boundary_xmin & -boundary_xmax & +boundary_ymin & -boundary_ymax & +boundary_zmin & -boundary_zmax
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
source.space = openmc.stats.Point((0.0, 0.0, 57.295))
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
