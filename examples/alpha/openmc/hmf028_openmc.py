"""
HEU-MET-FAST-028-1:  17.840 kg Oy(93.24) sphere in 7.09" Nat-U (Flattop)
Translated from COG to OpenMC
"""

import openmc
import numpy as np

# ==============================================================================
# Materials
# ==============================================================================

# Oy(93.24) @ 18.62 g/cc
mat1 = openmc.Material(material_id=1, name="Oy(93.24) @ 18.62 g/cc")
mat1.set_density("atom/b-cm", 4.767449e-02)
mat1.add_nuclide("U234", 4.8869e-4)
mat1.add_nuclide("U235", 4.4482e-2)
mat1.add_nuclide("U238", 2.7038e-3)

# Natural-U @ 19.00 g/cc
mat2 = openmc.Material(material_id=2, name="Natural-U @ 19.00 g/cc")
mat2.set_density("atom/b-cm", 4.806974e-02)
mat2.add_nuclide("U234", 2.6438e-6)
mat2.add_nuclide("U235", 3.4610e-4)
mat2.add_nuclide("U238", 4.7721e-2)

materials = openmc.Materials([mat1, mat2])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================


# Cell: Oy
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="Oy")
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
