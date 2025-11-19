"""
HEU-MET-FAST-001:  GODIVA (52.420 kg U(93.71) @ 18.74 g/cc)
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1, name="")
mat1.set_density("atom/b-cm", 4.798424e-02)
mat1.add_nuclide("U234", 4.9184e-4)
mat1.add_nuclide("U235", 4.4994e-2)
mat1.add_nuclide("U238", 2.4984e-3)

materials = openmc.Materials([mat1])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# per Section 3.2
surf1 = openmc.Sphere(surface_id=1, r=8.7407)


# Cell: CORE
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="CORE")
cell0.region = -surf1

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
outer_cell = openmc.Cell(cell_id=1, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, outer_cell])
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
