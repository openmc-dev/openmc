"""
SHE-8 (JAERI-1257) Homogenized model
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Core
mat1 = openmc.Material(material_id=1, name="Core")
mat1.set_density("atom/b-cm", 7.974066e-02)
mat1.add_nuclide("U235", 0.3415e-4)
mat1.add_nuclide("U238", 1.379e-4)
mat1.add_element("C", 0.7911e-1)
mat1.add_element("H", 0.7191e-4)
mat1.add_nuclide("O16", 0.3867e-3)
mat1.add_s_alpha_beta("c_H_in_H2O")
mat1.add_s_alpha_beta("c_Graphite")

# Reflector
mat2 = openmc.Material(material_id=2, name="Reflector")
mat2.set_density("atom/b-cm", 7.742810e-02)
mat2.add_element("C", 0.7732e-1)
mat2.add_element("H", 0.7137e-4)
mat2.add_nuclide("O16", 0.3673e-4)
mat2.add_s_alpha_beta("c_H_in_H2O")
mat2.add_s_alpha_beta("c_Graphite")

materials = openmc.Materials([mat1, mat2])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

surf1 = openmc.ZCylinder(surface_id=1, r=28.73)

# Hexagonal prism (COG "pri 6" surface with edge_length=150, z from -120 to 120)
surf2 = openmc.model.HexagonalPrism(edge_length=150.0,
                                     origin=(0.0, 0.0), orientation='x')
surf_zmin = openmc.ZPlane(surface_id=3, z0=-120.0)
surf_zmax = openmc.ZPlane(surface_id=4, z0=120.0)


# Cell: Core
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="Core")
cell0.region = -surf1 & -surf2 & +surf_zmin & -surf_zmax

# Cell: Refl
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="Refl")
cell1.region = +surf1 & -surf2 & +surf_zmin & -surf_zmax

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
settings.particles = 100
settings.batches = 520
settings.inactive = 21
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

# Add a simple flux tally to avoid warning
tally = openmc.Tally(tally_id=1, name='flux')
tally.scores = ['flux']
tallies.append(tally)

tallies.export_to_xml()

# ==============================================================================
# Run OpenMC
# ==============================================================================

openmc.run()
