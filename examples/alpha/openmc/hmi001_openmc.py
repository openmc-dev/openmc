"""
HEU-MET-INTER-001: ZPR-9/34 Loading 303; Benchmark Model
Translated from COG to OpenMC
"""

import openmc

# ==============================================================================
# Materials
# ==============================================================================

mat1 = openmc.Material(material_id=1, name="")
mat1.set_density("atom/b-cm", 4.800706e-03)
mat1.add_nuclide("U235", 1.13762e-3)
mat1.add_nuclide("U238", 6.59390e-5)
mat1.add_nuclide("U234", 1.11576e-5)
mat1.add_nuclide("U236", 5.34915e-6)
mat1.add_element("Cr", 2.48397e-3)
mat1.add_element("Ni", 1.09667e-3)

mat2 = openmc.Material(material_id=2, name="")
mat2.set_density("atom/b-cm", 4.830156e-03)
mat2.add_nuclide("U235", 1.13252e-3)
mat2.add_nuclide("U238", 6.56432e-5)
mat2.add_nuclide("U234", 1.11077e-5)
mat2.add_nuclide("U236", 5.32515e-6)
mat2.add_element("Cr", 2.50498e-3)
mat2.add_element("Ni", 1.11058e-3)

mat3 = openmc.Material(material_id=3, name="")
mat3.set_density("atom/b-cm", 4.509931e-03)
mat3.add_nuclide("U235", 1.02854e-3)
mat3.add_nuclide("U238", 5.96166e-5)
mat3.add_nuclide("U234", 1.00880e-5)
mat3.add_nuclide("U236", 4.83609e-6)
mat3.add_element("Cr", 2.32285e-3)
mat3.add_element("Ni", 1.08400e-3)

mat4 = openmc.Material(material_id=4, name="")
mat4.set_density("atom/b-cm", 4.575898e-03)
mat4.add_nuclide("U235", 1.03816e-3)
mat4.add_nuclide("U238", 6.01743e-5)
mat4.add_nuclide("U234", 1.01821e-5)
mat4.add_nuclide("U236", 4.88142e-6)
mat4.add_element("Cr", 2.35918e-3)
mat4.add_element("Ni", 1.10332e-3)

mat5 = openmc.Material(material_id=5, name="")
mat5.set_density("atom/b-cm", 7.606109e-02)
mat5.add_element("Cr", 1.49836e-2)
mat5.add_element("Ni", 6.61034e-3)
mat5.add_element("Fe", 5.32325e-2)
mat5.add_element("Al", 1.15377e-3)
mat5.add_element("C", 7.23704e-5)
mat5.add_element("Mo", 8.51061e-6)
mat5.add_s_alpha_beta("c_Graphite")

mat6 = openmc.Material(material_id=6, name="")
mat6.set_density("atom/b-cm", 7.666411e-02)
mat6.add_element("Cr", 1.49334e-2)
mat6.add_element("Ni", 6.45653e-3)
mat6.add_element("Fe", 5.27897e-2)
mat6.add_element("Al", 2.35245e-3)
mat6.add_element("C", 8.21637e-5)
mat6.add_element("Mo", 4.98640e-5)
mat6.add_s_alpha_beta("c_Graphite")

mat7 = openmc.Material(material_id=7, name="")
mat7.set_density("atom/b-cm", 4.814991e-02)
mat7.add_nuclide("U235", 8.59476e-5)
mat7.add_nuclide("U238", 3.84988e-2)
mat7.add_element("Cr", 1.86022e-3)
mat7.add_element("Ni", 8.08917e-4)
mat7.add_element("Fe", 6.86642e-3)
mat7.add_element("C", 2.96046e-5)
mat7.add_s_alpha_beta("c_Graphite")

materials = openmc.Materials([mat1, mat2, mat3, mat4, mat5, mat6, mat7])
materials.export_to_xml()

# ==============================================================================
# Geometry
# ==============================================================================


# Cell: CORE1
cell0 = openmc.Cell(cell_id=0, fill=mat1, name="CORE1")
cell0.region = -surf1 & +surf6 & -surf7

# Cell: CORE2
cell1 = openmc.Cell(cell_id=1, fill=mat2, name="CORE2")
cell1.region = -surf1 & +surf5 & -surf6

# Cell: CORE3
cell2 = openmc.Cell(cell_id=2, fill=mat3, name="CORE3")
cell2.region = -surf1 & +surf4 & -surf5

# Cell: CORE4
cell3 = openmc.Cell(cell_id=3, fill=mat4, name="CORE4")
cell3.region = -surf1 & +surf7 & -surf8

# Cell: AXSST
cell4 = openmc.Cell(cell_id=4, fill=mat5, name="AXSST")
cell4.region = -surf1 & -surf3 & +surf8

# Cell: AXSST
cell5 = openmc.Cell(cell_id=5, fill=mat5, name="AXSST")
cell5.region = -surf1 & -surf3 & -surf4

# Cell: RDSST
cell6 = openmc.Cell(cell_id=6, fill=mat6, name="RDSST")
cell6.region = +surf1 & -surf2 & -surf3

# Cell: DURRS
cell7 = openmc.Cell(cell_id=7, fill=mat7, name="DURRS")
cell7.region = +surf2 & -surf3

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
outer_region = outer_region & ~cell5.region
outer_region = outer_region & ~cell6.region
outer_region = outer_region & ~cell7.region
outer_cell = openmc.Cell(cell_id=8, name="outer_void")
outer_cell.region = outer_region
outer_cell.fill = None  # Void

# Create root universe and geometry
root_universe = openmc.Universe(cells=[cell0, cell1, cell2, cell3, cell4, cell5, cell6, cell7, outer_cell])
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
