"""
LST004-3: Water-reflected 60-cm-diameter cylindrical tank with 10% enriched uranyl nitrate solution @ H/X=842 (Run 33)
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Solution: Run 33
mat1 = openmc.Material(material_id=1, name="Solution: Run 33")
mat1.set_density("atom/b-cm", 3.309929e+01)
mat1.add_nuclide("U234", 5.5579e-7)
mat1.add_nuclide("U235", 6.8970e-5)
mat1.add_nuclide("U236", 6.8884e-8)
mat1.add_nuclide("U238", 6.1432e-4)
mat1.add_element("H", 5.8085e-2)
mat1.add_element("N", 2.6927e-3)
mat1.add_nuclide("O16", 3.7826e-2)
mat1.add_element("Run", 33)
mat1.add_s_alpha_beta("c_H_in_H2O")

# Stainless Steel
mat2 = openmc.Material(material_id=2, name="Stainless Steel")
mat2.set_density("atom/b-cm", 8.684498e-02)
mat2.add_element("C", 4.3736e-5)
mat2.add_element("Si", 1.0627e-3)
mat2.add_element("Mn", 1.1561e-3)
mat2.add_element("P", 4.3170e-5)
mat2.add_element("S", 2.9782e-6)
mat2.add_element("Ni", 8.3403e-3)
mat2.add_element("Cr", 1.6775e-2)
mat2.add_element("Fe", 5.9421e-2)
mat2.add_s_alpha_beta("c_Graphite")

# Water
mat3 = openmc.Material(material_id=3, name="Water")
mat3.set_density("atom/b-cm", 9.998700e-02)
mat3.add_element("H", 6.6658e-2)
mat3.add_nuclide("O16", 3.3329e-2)
mat3.add_s_alpha_beta("c_H_in_H2O")

# Air
mat4 = openmc.Material(material_id=4, name="Air")
mat4.set_density("atom/b-cm", 4.942500e-05)
mat4.add_element("N", 3.9016e-5)
mat4.add_nuclide("O16", 1.0409e-5)

materials = openmc.Materials([mat1, mat2, mat3, mat4])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

surf1 = openmc.ZCylinder(surface_id=1, r=29.5)

surf2 = openmc.ZCylinder(surface_id=2, r=29.8)

surf3 = openmc.ZCylinder(surface_id=3, r=59.8)

# Hc
surf4 = openmc.ZPlane(surface_id=4, z0=52.93)


# Cell: Soln
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="Soln")
cell0.region = -surf1 & -surf4

# Cell: SST
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="SST")
cell1.region = +surf1 & -surf2

# Cell: H2O
cell2 = openmc.Cell(cell_id=2, fill=mat3, name="H2O")
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
source.space = openmc.stats.Point((0.0, 0.0, 26.465))
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
