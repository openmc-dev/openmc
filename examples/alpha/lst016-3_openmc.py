"""
LST016-3: Water reflected 28-cm-thick slab tank with 10% enriched uranyl nitrate solution @ H/X=608 (Run 125)
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

# Solution: Run 125
mat1 = openmc.Material(material_id=1, name="Solution: Run 125")
mat1.set_density("atom/b-cm", 1.250987e+02)
mat1.add_nuclide("U234", 7.6555e-7)
mat1.add_nuclide("U235", 9.4999e-5)
mat1.add_nuclide("U236", 9.4881e-8)
mat1.add_nuclide("U238", 8.4617e-4)
mat1.add_element("H", 5.7800e-2)
mat1.add_element("N", 2.3658e-3)
mat1.add_nuclide("O16", 3.7641e-2)
mat1.add_element("Run", 125)
mat1.add_s_alpha_beta("c_H_in_H2O")

# Stainless Steel
mat2 = openmc.Material(material_id=2, name="Stainless Steel")
mat2.set_density("atom/b-cm", 8.668297e-02)
mat2.add_element("C", 7.1567e-5)
mat2.add_element("Si", 7.1415e-4)
mat2.add_element("Mn", 9.9095e-4)
mat2.add_element("P", 5.0879e-5)
mat2.add_element("S", 1.0424e-5)
mat2.add_element("Ni", 8.5600e-3)
mat2.add_element("Cr", 1.6725e-2)
mat2.add_element("Fe", 5.9560e-2)
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

# Slab tank/inner
# Box (6 planes): xmin=-14.04, xmax=14.04, ymin=-34.515, ymax=34.515, zmin=0.0, zmax=149.75
surf1_xmin = openmc.XPlane(surface_id=10000, x0=-14.04)
surf1_xmax = openmc.XPlane(surface_id=10001, x0=14.04)
surf1_ymin = openmc.YPlane(surface_id=10002, y0=-34.515)
surf1_ymax = openmc.YPlane(surface_id=10003, y0=34.515)
surf1_zmin = openmc.ZPlane(surface_id=10004, z0=0.0)
surf1_zmax = openmc.ZPlane(surface_id=10005, z0=149.75)

# Slab tank/outer
# Box (6 planes): xmin=-16.57, xmax=16.57, ymin=-37.045, ymax=37.045, zmin=-2.039999999999992, zmax=152.63
surf2_xmin = openmc.XPlane(surface_id=10006, x0=-16.57)
surf2_xmax = openmc.XPlane(surface_id=10007, x0=16.57)
surf2_ymin = openmc.YPlane(surface_id=10008, y0=-37.045)
surf2_ymax = openmc.YPlane(surface_id=10009, y0=37.045)
surf2_zmin = openmc.ZPlane(surface_id=10010, z0=-2.039999999999992)
surf2_zmax = openmc.ZPlane(surface_id=10011, z0=152.63)

# Water reflector
# Box (6 planes): xmin=-46.57, xmax=46.57, ymin=-67.045, ymax=67.045, zmin=-32.03999999999999, zmax=172.63
surf3_xmin = openmc.XPlane(surface_id=10012, x0=-46.57)
surf3_xmax = openmc.XPlane(surface_id=10013, x0=46.57)
surf3_ymin = openmc.YPlane(surface_id=10014, y0=-67.045)
surf3_ymax = openmc.YPlane(surface_id=10015, y0=67.045)
surf3_zmin = openmc.ZPlane(surface_id=10016, z0=-32.03999999999999)
surf3_zmax = openmc.ZPlane(surface_id=10017, z0=172.63)

# Hc
surf4 = openmc.ZPlane(surface_id=4, z0=51.37)


# Cell: Soln
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="Soln")
cell0.region = (+surf1_xmin & -surf1_xmax & +surf1_ymin & -surf1_ymax & +surf1_zmin & -surf1_zmax) & -surf4

# Cell: SST
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="SST")
cell1.region = (-surf1_xmin | +surf1_xmax | -surf1_ymin | +surf1_ymax | -surf1_zmin | +surf1_zmax) & (+surf2_xmin & -surf2_xmax & +surf2_ymin & -surf2_ymax & +surf2_zmin & -surf2_zmax)

# Cell: Water
cell2 = openmc.Cell(cell_id=2, fill=mat3, name="Water")
cell2.region = (-surf2_xmin | +surf2_xmax | -surf2_ymin | +surf2_ymax | -surf2_zmin | +surf2_zmax) & (+surf3_xmin & -surf3_xmax & +surf3_ymin & -surf3_ymax & +surf3_zmin & -surf3_zmax)

# Cell: Air
cell3 = openmc.Cell(cell_id=3, fill=mat4, name="Air")
cell3.region = (+surf1_xmin & -surf1_xmax & +surf1_ymin & -surf1_ymax & +surf1_zmin & -surf1_zmax) & +surf4

# ==============================================================================
# Boundary Conditions
# ==============================================================================

# Create outer bounding box with vacuum boundary (6 planes)
# TODO: Adjust dimensions to encompass your entire geometry
boundary_xmin = openmc.XPlane(surface_id=10018, x0=-200, boundary_type="vacuum")
boundary_xmax = openmc.XPlane(surface_id=10019, x0=200, boundary_type="vacuum")
boundary_ymin = openmc.YPlane(surface_id=10020, y0=-200, boundary_type="vacuum")
boundary_ymax = openmc.YPlane(surface_id=10021, y0=200, boundary_type="vacuum")
boundary_zmin = openmc.ZPlane(surface_id=10022, z0=-200, boundary_type="vacuum")
boundary_zmax = openmc.ZPlane(surface_id=10023, z0=200, boundary_type="vacuum")

# Create outer void cell (everything outside geometry but inside boundary)
# Particles are killed at the vacuum boundary
outer_region = +boundary_xmin & -boundary_xmax & +boundary_ymin & -boundary_ymax & +boundary_zmin & -boundary_zmax
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
