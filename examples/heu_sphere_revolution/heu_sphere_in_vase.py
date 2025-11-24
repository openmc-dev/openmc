"""
HEU Sphere in Vase-like Revolution Container

This example demonstrates the Revolution surface feature in OpenMC by modeling
a subcritical HEU (93.3% U-235) sphere inside a vase-like container created
using the Revolution surface (a solid of revolution).

The geometry consists of:
1. A subcritical HEU sphere (3 cm radius, ~2 kg)
2. A vase-like aluminum container modeled using Revolution
3. Air filling the space between the sphere and container

The model generates a 3D voxel plot that can be converted to VTK format
for visualization in VisIt or ParaView.

To run this example:
    python heu_sphere_in_vase.py

To visualize in VisIt:
    1. Run the script (generates plot_1.h5 and geometry.vti)
    2. Open geometry.vti in VisIt
    3. Add a Pseudocolor plot of the "id" variable
    4. Apply a threshold or subset to visualize different materials
"""

import openmc


def create_materials():
    """Create materials for the model."""

    # ===========================================================================
    # HEU Material (93.3% U-235 enriched uranium)
    # ===========================================================================
    # Isotopic composition for 93.3% enriched uranium
    # Based on typical HEU isotopic vector:
    # - U-234: ~1.0 wt%
    # - U-235: 93.3 wt%
    # - U-238: ~5.7 wt%
    # Metal density: 18.74 g/cc

    heu = openmc.Material(name='HEU')
    heu.set_density('g/cc', 18.74)
    heu.add_nuclide('U234', 1.0, percent_type='wo')
    heu.add_nuclide('U235', 93.3, percent_type='wo')
    heu.add_nuclide('U238', 5.7, percent_type='wo')

    # ===========================================================================
    # Aluminum (container material)
    # ===========================================================================
    aluminum = openmc.Material(name='Aluminum')
    aluminum.set_density('g/cc', 2.7)
    aluminum.add_element('Al', 1.0)

    # ===========================================================================
    # Air (gap filler)
    # ===========================================================================
    air = openmc.Material(name='Air')
    air.set_density('g/cc', 0.001205)
    air.add_element('N', 0.7809, percent_type='wo')
    air.add_element('O', 0.2095, percent_type='wo')
    air.add_element('Ar', 0.0093, percent_type='wo')
    air.add_element('C', 0.0003, percent_type='wo')

    return openmc.Materials([heu, aluminum, air])


def create_geometry(materials):
    """Create the geometry with HEU sphere inside a vase-like Revolution container."""

    heu = materials[0]
    aluminum = materials[1]
    air = materials[2]

    # ===========================================================================
    # Surfaces
    # ===========================================================================

    # HEU sphere (3 cm radius - subcritical)
    # Critical mass for bare HEU sphere is ~52 kg at r=8.74 cm
    # This 3 cm sphere has mass ~2.1 kg, well below critical
    sphere_radius = 3.0
    heu_sphere = openmc.Sphere(r=sphere_radius, name='HEU Sphere')

    # Vase-like container using Revolution surface
    # The vase profile is defined as (r, z) coordinate pairs
    # The profile is revolved around the z-axis
    # Container dimensions chosen to enclose the sphere with clearance

    container_thickness = 0.5  # cm wall thickness

    # Outer vase profile (r, z) - a classic vase shape
    # Starting from bottom and going up
    vase_outer_profile = [
        (2.0, -6.0),    # Narrow base
        (4.0, -4.0),    # Widening
        (6.0, -1.0),    # Bulge starts
        (7.0, 2.0),     # Maximum bulge (encloses sphere)
        (6.5, 5.0),     # Narrowing toward top
        (5.0, 7.0),     # Neck
        (4.5, 9.0),     # Flared rim
        (5.5, 10.0),    # Rim top
    ]

    # Inner vase profile (offset inward by wall thickness)
    vase_inner_profile = [
        (1.5, -5.5),    # Narrow base (inner)
        (3.5, -3.5),    # Widening
        (5.5, -0.5),    # Bulge starts
        (6.5, 2.0),     # Maximum bulge (sphere fits here)
        (6.0, 5.0),     # Narrowing
        (4.5, 7.0),     # Neck
        (4.0, 9.0),     # Near rim
        (5.0, 10.0),    # Rim top (inner)
    ]

    # Create Revolution surfaces
    vase_outer = openmc.Revolution(
        rz=vase_outer_profile,
        axis='z',
        origin=(0., 0., 0.),
        name='Vase Outer'
    )

    vase_inner = openmc.Revolution(
        rz=vase_inner_profile,
        axis='z',
        origin=(0., 0., 0.),
        name='Vase Inner'
    )

    # Bounding planes for the vase (top and bottom caps)
    bottom_plane = openmc.ZPlane(z0=-6.0, name='Bottom')
    top_plane = openmc.ZPlane(z0=10.0, name='Top')

    # Outer boundary sphere (vacuum boundary)
    outer_boundary = openmc.Sphere(
        r=15.0,
        boundary_type='vacuum',
        name='Outer Boundary'
    )

    # ===========================================================================
    # Cells
    # ===========================================================================

    # HEU sphere cell
    cell_heu = openmc.Cell(name='HEU Sphere', fill=heu)
    cell_heu.region = -heu_sphere

    # Air gap between sphere and inner vase wall
    cell_air = openmc.Cell(name='Air Gap', fill=air)
    cell_air.region = +heu_sphere & -vase_inner & +bottom_plane & -top_plane

    # Aluminum vase wall (between inner and outer Revolution surfaces)
    cell_vase = openmc.Cell(name='Vase Wall', fill=aluminum)
    cell_vase.region = +vase_inner & -vase_outer & +bottom_plane & -top_plane

    # Void outside vase but inside boundary
    cell_void = openmc.Cell(name='External Void')
    cell_void.region = (+vase_outer | -bottom_plane | +top_plane) & -outer_boundary

    # ===========================================================================
    # Universe and Geometry
    # ===========================================================================

    root_universe = openmc.Universe(cells=[cell_heu, cell_air, cell_vase, cell_void])
    geometry = openmc.Geometry(root_universe)

    return geometry


def create_settings():
    """Create settings for eigenvalue calculation."""

    settings = openmc.Settings()

    # Use fewer particles for demonstration (increase for production)
    settings.particles = 5000
    settings.batches = 50
    settings.inactive = 10
    settings.run_mode = 'eigenvalue'

    # Source in the HEU sphere
    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((0.0, 0.0, 0.0))
    settings.source = source

    return settings


def create_plots():
    """Create plots for geometry visualization."""

    plots = openmc.Plots()

    # ===========================================================================
    # 2D Slice plots for quick visualization
    # ===========================================================================

    # XZ slice (side view)
    plot_xz = openmc.Plot(name='XZ Slice')
    plot_xz.basis = 'xz'
    plot_xz.origin = (0., 0., 2.)
    plot_xz.width = (20., 20.)
    plot_xz.pixels = (800, 800)
    plot_xz.color_by = 'material'
    plots.append(plot_xz)

    # XY slice (top view through sphere center)
    plot_xy = openmc.Plot(name='XY Slice')
    plot_xy.basis = 'xy'
    plot_xy.origin = (0., 0., 0.)
    plot_xy.width = (20., 20.)
    plot_xy.pixels = (800, 800)
    plot_xy.color_by = 'material'
    plots.append(plot_xy)

    # ===========================================================================
    # 3D Voxel plot for VisIt visualization
    # ===========================================================================

    voxel_plot = openmc.Plot(name='3D Voxel')
    voxel_plot.type = 'voxel'
    voxel_plot.origin = (0., 0., 2.)  # Center on the geometry
    voxel_plot.width = (18., 18., 20.)  # Cover full geometry
    voxel_plot.pixels = (180, 180, 200)  # Resolution (adjust for file size)
    voxel_plot.color_by = 'material'
    plots.append(voxel_plot)

    return plots


def main():
    """Main function to build and export the model."""

    print("=" * 70)
    print("HEU Sphere in Vase-like Revolution Container")
    print("=" * 70)
    print()

    # Create all components
    print("Creating materials...")
    materials = create_materials()

    print("Creating geometry with Revolution surface...")
    geometry = create_geometry(materials)

    print("Creating settings...")
    settings = create_settings()

    print("Creating plots (including 3D voxel for VisIt)...")
    plots = create_plots()

    # Export to XML
    print()
    print("Exporting to XML files...")
    materials.export_to_xml()
    geometry.export_to_xml()
    settings.export_to_xml()
    plots.export_to_xml()

    print("  - materials.xml")
    print("  - geometry.xml")
    print("  - settings.xml")
    print("  - plots.xml")

    # Generate plots
    print()
    print("Generating geometry plots...")
    openmc.plot_geometry()

    # Convert voxel plot to VTK for VisIt
    print()
    print("Converting voxel plot to VTK format for VisIt...")
    try:
        vtk_file = openmc.voxel_to_vtk('plot_3.h5', 'geometry.vti')
        print(f"  Created: {vtk_file}")
        print()
        print("To visualize in VisIt:")
        print("  1. Open geometry.vti in VisIt")
        print("  2. Add a Pseudocolor plot of 'id'")
        print("  3. Use Subset or Threshold to isolate materials")
    except ImportError:
        print("  Note: VTK package not installed. Install with: pip install vtk")
        print("  The voxel HDF5 file (plot_3.h5) was still created.")
    except Exception as e:
        print(f"  Warning: Could not convert to VTK: {e}")
        print("  The voxel HDF5 file was still created.")

    print()
    print("=" * 70)
    print("Model files created successfully!")
    print()
    print("To run the eigenvalue calculation:")
    print("  openmc")
    print()
    print("Note: This model uses a SUBCRITICAL HEU sphere (k-eff < 1)")
    print("      The 3 cm radius sphere has mass ~2.1 kg (critical mass ~52 kg)")
    print("=" * 70)


if __name__ == '__main__':
    main()
