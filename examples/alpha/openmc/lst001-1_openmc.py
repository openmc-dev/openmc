"""
LST001-1: SHEBA-II
Translated from COG to OpenMC
"""

import openmc
import numpy as np

# ==============================================================================
# Materials
# ==============================================================================

# Fuel
mat1 = openmc.Material(material_id=1, name="Fuel")
mat1.set_density("atom/b-cm", 2.476457e-03)
mat1.add_nuclide("U234", 6.7855e-7)
mat1.add_nuclide("U235", 1.2377e-4)
mat1.add_nuclide("U236", 1.2085e-6)
mat1.add_nuclide("U238", 2.3508e-3)

# Air
mat2 = openmc.Material(material_id=2, name="Air")
mat2.set_density("atom/b-cm", 5.030600e-05)
mat2.add_nuclide("N", 3.5214e-5)
mat2.add_nuclide("O16", 1.5092e-5)

# SS304L
mat3 = openmc.Material(material_id=3, name="SS304L")
mat3.set_density("atom/b-cm", 8.534700e-02)
mat3.add_nuclide("Cr", 1.6348e-2)
mat3.add_nuclide("Mn", 1.7192e-3)
mat3.add_nuclide("Fe", 6.0038e-2)
mat3.add_nuclide("Ni", 7.2418e-3)

materials = openmc.Materials([mat1, mat2, mat3])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

surf1 = openmc.ZCylinder(surface_id=1, r=2.54)

surf2 = openmc.ZCylinder(surface_id=2, r=3.175)

surf3 = openmc.ZCylinder(surface_id=3, r=24.4475)

surf4 = openmc.ZCylinder(surface_id=4, r=25.4)

surf5 = openmc.ZPlane(surface_id=5, z0=44.8)


# Cell: Soln
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="Soln")
cell0.region = +surf2 & -surf3 & -surf5

# Cell: Air
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="Air")
cell1.region = +surf2 & -surf3 & +surf5

# Cell: Air
cell2 = openmc.Cell(cell_id=2, fill=mat2, name="Air")
cell2.region = -surf1 & -surf4

# Cell: SS304L
cell3 = openmc.Cell(cell_id=3, fill=mat3, name="SS304L")
cell3.region = +surf1 & +surf3 & -surf4

# Cell: SS304L
cell4 = openmc.Cell(cell_id=4, fill=mat3, name="SS304L")
cell4.region = +surf1 & -surf2 & -surf3

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
source.space = openmc.stats.Point((-5.0, 0.0, 22.4))
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
