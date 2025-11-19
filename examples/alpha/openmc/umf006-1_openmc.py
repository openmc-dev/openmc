"""
U233-MET-FAST-006: 5.740 kg 233U(98.1) in FLATTOP
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# U-233 per Table 2
mat1 = openmc.Material(material_id=1, name="U-233 per Table 2")
mat1.set_density("atom/b-cm", 1.204200e+02)
mat1.add_element("P", 18.42)
mat1.add_nuclide("U233", 98.13)
mat1.add_nuclide("U234", 1.24)
mat1.add_nuclide("U235", 0.03)
mat1.add_nuclide("U238", 0.60)
mat1.add_element("Table", 2)

# Nat-U
mat2 = openmc.Material(material_id=2, name="Nat-U")
mat2.set_density("atom/b-cm", 1.190000e+02)
mat2.add_element("F", 19.00)
mat2.add_nuclide("U234", 0.005)
mat2.add_nuclide("U235", 0.72)
mat2.add_nuclide("U238", 99.275)

materials = openmc.Materials([mat1, mat2])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# Dimensions per
surf1 = openmc.Sphere(surface_id=1, r=4.2058)

# Section 3.2
surf2 = openmc.Sphere(surface_id=2, r=24.1194)


# Cell: U233
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="U233")
cell0.region = -surf1

# Cell: Tu
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="Tu")
cell1.region = +surf1 & -surf2

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
outer_cell = openmc.Cell(cell_id=2, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, outer_cell])
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
