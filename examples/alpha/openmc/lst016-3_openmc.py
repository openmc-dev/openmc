"""
LST016-3: Water reflected 28-cm-thick slab tank with 10% enriched uranyl nitrate solution @ H/X=608 (Run 125)
Translated from COG to OpenMC
"""

import openmc
import numpy as np

# ==============================================================================
# Materials
# ==============================================================================

# Solution: Run 125
mat1 = openmc.Material(material_id=1, name="Solution: Run 125")
mat1.set_density("atom/b-cm", 1.250409e+02)
mat1.add_nuclide("U234", 7.6555e-7)
mat1.add_nuclide("U235", 9.4999e-5)
mat1.add_nuclide("U236", 9.4881e-8)
mat1.add_nuclide("U238", 8.4617e-4)
mat1.add_nuclide("N", 2.3658e-3)
mat1.add_nuclide("O16", 3.7641e-2)
mat1.add_nuclide("Run", 125)
mat1.add_s_alpha_beta("c_H_in_H2O")

# Stainless Steel
mat2 = openmc.Material(material_id=2, name="Stainless Steel")
mat2.set_density("atom/b-cm", 8.668297e-02)
mat2.add_nuclide("C0", 7.1567e-5)
mat2.add_nuclide("Si", 7.1415e-4)
mat2.add_nuclide("Mn", 9.9095e-4)
mat2.add_nuclide("P", 5.0879e-5)
mat2.add_nuclide("S", 1.0424e-5)
mat2.add_nuclide("Ni", 8.5600e-3)
mat2.add_nuclide("Cr", 1.6725e-2)
mat2.add_nuclide("Fe", 5.9560e-2)

# Water
mat3 = openmc.Material(material_id=3, name="Water")
mat3.set_density("atom/b-cm", 3.332900e-02)
mat3.add_nuclide("O16", 3.3329e-2)
mat3.add_s_alpha_beta("c_H_in_H2O")

# Air
mat4 = openmc.Material(material_id=4, name="Air")
mat4.set_density("atom/b-cm", 4.942500e-05)
mat4.add_nuclide("N", 3.9016e-5)
mat4.add_nuclide("O16", 1.0409e-5)

materials = openmc.Materials([mat1, mat2, mat3, mat4])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# Slab tank/inner
surf1 = openmc.model.RectangularParallelepiped(
    -14.04, 14.04, -34.515, 34.515, 0.0, 149.75, surface_id=1)

# Slab tank/outer
surf2 = openmc.model.RectangularParallelepiped(
    -16.57, 16.57, -37.045, 37.045, -2.039999999999992, 152.63, surface_id=2)

# Water reflector
surf3 = openmc.model.RectangularParallelepiped(
    -46.57, 46.57, -67.045, 67.045, -32.03999999999999, 172.63, surface_id=3)

# Hc
surf4 = openmc.ZPlane(surface_id=4, z0=51.37)


# Cell: Soln
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="Soln")
cell0.region = -surf1 & -surf4

# Cell: SST
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="SST")
cell1.region = +surf1 & -surf2

# Cell: Water
cell2 = openmc.Cell(cell_id=2, fill=mat3, name="Water")
cell2.region = +surf2 & -surf3

# Cell: Air
cell3 = openmc.Cell(cell_id=3, fill=mat4, name="Air")
cell3.region = -surf1 & +surf4

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
outer_cell = openmc.Cell(cell_id=4, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, outer_cell])
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
source.space = openmc.stats.Point((0.0, 0.0, 25.685))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Watt(a=0.988e6, b=2.249e-6)
settings.source = source

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
