"""
LST021-1: Unreflected 80-cm diameter cylindrical tank with 10% enriched uranyl nitrate solution @ H/X=971 (Run 215)
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Solution: Run 215
mat1 = openmc.Material(material_id=1, name="Solution: Run 215")
mat1.set_density("atom/b-cm", 2.150993e+02)
mat1.add_nuclide("U234", 5.0042e-7)
mat1.add_nuclide("U235", 6.2098e-5)
mat1.add_nuclide("U236", 6.2021e-8)
mat1.add_nuclide("U238", 5.5312e-4)
mat1.add_element("H", 6.0297e-2)
mat1.add_element("N", 1.8157e-3)
mat1.add_nuclide("O16", 3.6535e-2)
mat1.add_element("Run", 215)
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

# Air
mat3 = openmc.Material(material_id=3, name="Air")
mat3.set_density("atom/b-cm", 4.942400e-05)
mat3.add_element("N", 3.9014e-5)
mat3.add_nuclide("O16", 1.0410e-5)

materials = openmc.Materials([mat1, mat2, mat3])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================

# Tank/inner
surf1 = openmc.ZCylinder(surface_id=1, r=39.505)

# Tank/outer
surf2 = openmc.ZCylinder(surface_id=2, r=39.815)

# Base plate
# Box (6 planes): xmin=-60.2, xmax=39.8, ymin=-50.0, ymax=50.0, zmin=-19.0, zmax=-16.0
surf3_xmin = openmc.XPlane(surface_id=10000, x0=-60.2)
surf3_xmax = openmc.XPlane(surface_id=10001, x0=39.8)
surf3_ymin = openmc.YPlane(surface_id=10002, y0=-50.0)
surf3_ymax = openmc.YPlane(surface_id=10003, y0=50.0)
surf3_zmin = openmc.ZPlane(surface_id=10004, z0=-19.0)
surf3_zmax = openmc.ZPlane(surface_id=10005, z0=-16.0)

# Hole in base plate
surf4 = openmc.ZCylinder(surface_id=4, r=7.76)

# boundary condition
surf5 = openmc.ZCylinder(surface_id=5, r=79.815)

# Hc
surf6 = openmc.ZPlane(surface_id=6, z0=43.98)


# Cell: Soln
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="Soln")
cell0.region = -surf1 & -surf6

# Cell: SST
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="SST")
cell1.region = +surf1 & -surf2

# Cell: SST
cell2 = openmc.Cell(cell_id=2, fill=mat2, name="SST")
cell2.region = (+surf3_xmin & -surf3_xmax & +surf3_ymin & -surf3_ymax & +surf3_zmin & -surf3_zmax) & +surf4

# ==============================================================================
# Boundary Conditions
# ==============================================================================

# Create outer bounding box with vacuum boundary (6 planes)
# TODO: Adjust dimensions to encompass your entire geometry
boundary_xmin = openmc.XPlane(surface_id=10006, x0=-200, boundary_type="vacuum")
boundary_xmax = openmc.XPlane(surface_id=10007, x0=200, boundary_type="vacuum")
boundary_ymin = openmc.YPlane(surface_id=10008, y0=-200, boundary_type="vacuum")
boundary_ymax = openmc.YPlane(surface_id=10009, y0=200, boundary_type="vacuum")
boundary_zmin = openmc.ZPlane(surface_id=10010, z0=-200, boundary_type="vacuum")
boundary_zmax = openmc.ZPlane(surface_id=10011, z0=200, boundary_type="vacuum")

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
source.space = openmc.stats.Point((0.0, 0.0, 21.99))
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
