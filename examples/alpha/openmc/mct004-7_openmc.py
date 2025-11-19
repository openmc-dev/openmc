"""
MIX-COMP-THERM-004-7: 20x20 MOX pins with 2.225 cm pitch; Hc=60.32
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# MOX
mat1 = openmc.Material(material_id=1, name="MOX")
mat1.set_density("atom/b-cm", 4.027335e-04)
mat1.add_nuclide("Pu238", 2.0003e-6)
mat1.add_nuclide("Pu239", 2.7491e-4)
mat1.add_nuclide("Pu240", 8.8417e-5)
mat1.add_nuclide("Pu241", 2.8133e-5)
mat1.add_nuclide("Pu242", 8.1234e-6)
mat1.add_nuclide("Am241", 1.1498e-6)

# Cladding (with air gap)
mat2 = openmc.Material(material_id=2, name="Cladding (with air gap)")
mat2.set_density("atom/b-cm", 3.839992e-02)
mat2.add_element("Zr", 3.7772e-2)
mat2.add_element("Sn", 4.3737e-4)
mat2.add_element("Fe", 8.8570e-5)
mat2.add_element("Cr", 6.6119e-5)
mat2.add_element("Ni", 3.5864e-5)

# Water
mat3 = openmc.Material(material_id=3, name="Water")
mat3.set_density("atom/b-cm", 1.001030e-01)
mat3.add_element("H", 6.6735e-2)
mat3.add_nuclide("O16", 3.3368e-2)
mat3.add_s_alpha_beta("c_H_in_H2O")

# Aluminum
mat4 = openmc.Material(material_id=4, name="Aluminum")
mat4.set_density("atom/b-cm", 6.022400e-02)
mat4.add_element("Al", 6.0224e-2)

# SS304L
mat5 = openmc.Material(material_id=5, name="SS304L")
mat5.set_density("atom/b-cm", 1.262273e-02)
mat5.add_element("C", 1.1928e-4)
mat5.add_element("Si", 1.7003e-3)
mat5.add_element("Mn", 1.7385e-3)
mat5.add_element("P", 6.9381e-5)
mat5.add_element("S", 4.4673e-5)
mat5.add_element("Ni", 8.9506e-3)
mat5.add_s_alpha_beta("c_Graphite")

# Ordinary
mat6 = openmc.Material(material_id=6, name="Ordinary")
mat6.set_density("atom/b-cm", 6.260505e-02)
mat6.add_element("H", 1.3742e-2)
mat6.add_nuclide("O16", 4.5919e-2)
mat6.add_element("C", 1.1532e-4)
mat6.add_element("Na", 9.6395e-4)
mat6.add_element("Mg", 1.2388e-4)
mat6.add_element("Al", 1.7409e-3)
mat6.add_s_alpha_beta("c_H_in_H2O")
mat6.add_s_alpha_beta("c_Graphite")

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6])
materials.export_to_xml()

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
surf5 = openmc.model.RectangularParallelepiped(
    -22.25, 22.25, -22.25, 22.25, -499.995, 499.995, surface_id=5)

# BCD
surf6 = openmc.model.RectangularParallelepiped(
    -52.25, 52.25, -52.25, 52.25, -71.60000000000001, 81.17, surface_id=6)

# Middle grid, bottom
surf10 = openmc.ZPlane(surface_id=10, z0=80.57)

# Water height
surf11 = openmc.ZPlane(surface_id=11, z0=60.32)

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


# ==============================================================================
# Universes (from COG define unit blocks)
# ==============================================================================

# Unit 1: Mox fuel pin
u1_cell0 = openmc.Cell(fill=mat1, name="Mox")
u1_cell0.region = -surf1 & -surf4
u1_cell1 = openmc.Cell(fill=mat4, name="Al")
u1_cell1.region = +surf1 & -surf2 & -surf4
u1_cell2 = openmc.Cell(fill=mat2, name="Clad")
u1_cell2.region = +surf1 & +surf2 & +surf3 & -surf4
universe1 = openmc.Universe(universe_id=1, cells=[u1_cell0, u1_cell1, u1_cell2])

# Unit 2: Everything else
u2_cell0 = openmc.Cell(fill=mat4, name="Al")
u2_cell0.region = -surf6 & +surf10
u2_cell1 = openmc.Cell(fill=mat3, name="H2O")
u2_cell1.region = -surf6 & -surf11 & +surf12
u2_cell2 = openmc.Cell(fill=mat4, name="Al")
u2_cell2.region = -surf6 & -surf12 & +surf13
u2_cell3 = openmc.Cell(fill=mat3, name="H2O")
u2_cell3.region = -surf6 & -surf13 & +surf14
u2_cell4 = openmc.Cell(fill=mat4, name="Al")
u2_cell4.region = -surf6 & -surf14 & +surf15
u2_cell5 = openmc.Cell(fill=mat5, name="SS304L")
u2_cell5.region = -surf6 & -surf15 & +surf16
u2_cell6 = openmc.Cell(fill=mat3, name="H2O")
u2_cell6.region = -surf6 & -surf16 & +surf17
u2_cell7 = openmc.Cell(fill=mat5, name="SS304L")
u2_cell7.region = -surf6 & -surf17 & +surf18
u2_cell8 = openmc.Cell(fill=mat6, name="Conc")
u2_cell8.region = -surf6 & -surf18
universe2 = openmc.Universe(universe_id=2, cells=[u2_cell0, u2_cell1, u2_cell2, u2_cell3, u2_cell4, u2_cell5, u2_cell6, u2_cell7, u2_cell8])

# Unit 3: Unit cell with a fuel rod
universe3 = openmc.Universe(universe_id=3, cells=[])

# Unit 4: 20x20 lattice
# Lattice for unit 4: 20x20 array
# Pitch: (2.225000, 2.225000) cm
# Lower left: (-22.25, -22.25)

# Lattice filled with universe3
lattice4 = openmc.RectLattice(lattice_id=4)
lattice4.lower_left = [-22.25, -22.25]
lattice4.pitch = [2.225000, 2.225000]
lattice4.universes = [[universe3]*20]*20
universe4 = lattice4  # Lattice can be used as universe

# Cell using unit 2: stuff
cell0 = openmc.Cell(cell_id=0, fill=universe2, name="stuff")
cell0.region = +surf5 & -surf6

# Cell: Clad
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="Clad")
cell1.region = +surf1 & +surf2 & +surf3 & -surf4

# Cell using unit 2: Alles
cell2 = openmc.Cell(cell_id=2, fill=universe2, name="Alles")
cell2.region = +surf4 & -surf6

# Cell using unit 4: lttc
cell3 = openmc.Cell(cell_id=3, fill=universe4, name="lttc")
cell3.region = -surf5 & -surf6

# Cell: Conc
cell4 = openmc.Cell(cell_id=4, fill=mat6, name="Conc")
cell4.region = -surf6 & -surf18

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
outer_region = outer_region & ~cell3.region
outer_region = outer_region & ~cell4.region
outer_cell = openmc.Cell(cell_id=5, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, outer_cell])
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
source.space = openmc.stats.Point((0.0, 0.0, 0.0))
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
