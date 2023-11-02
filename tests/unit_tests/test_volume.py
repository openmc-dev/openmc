import numpy as np
import pytest

import openmc


def test_infinity_handling():
    surf1 = openmc.Sphere(boundary_type="vacuum")
    cell1 = openmc.Cell(region=-surf1)

    lower_left = (-2, -np.inf, -2)
    upper_right = (np.inf, 2, 2)

    with pytest.raises(ValueError, match="must be finite"):
        openmc.VolumeCalculation([cell1], 100, lower_left, upper_right)


def test_volume_no_cross_section(run_in_tmpdir):
    # removing path to cross section
    del openmc.config['cross_sections']

    # setting the simulation
    uo2 = openmc.Material(1, "uo2")
    mat = openmc.Material()
    # Add nuclides to uo2
    uo2.add_nuclide('U235', 0.03)
    uo2.add_nuclide('U238', 0.97)
    uo2.add_nuclide('O16', 2.0)
    uo2.set_density('g/cm3', 10.0)
    zirconium = openmc.Material(name="zirconium")
    zirconium.add_element('Zr', 1.0)
    zirconium.set_density('g/cm3', 6.6)

    water = openmc.Material(name="h2o")
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')

    materials = openmc.Materials([uo2, zirconium, water])
    materials.export_to_xml()

    fuel_outer_radius = openmc.ZCylinder(r=0.39)
    clad_inner_radius = openmc.ZCylinder(r=0.40)
    clad_outer_radius = openmc.ZCylinder(r=0.46)
    fuel_region = -fuel_outer_radius
    gap_region = +fuel_outer_radius & -clad_inner_radius
    clad_region = +clad_inner_radius & -clad_outer_radius

    fuel = openmc.Cell(name='fuel')
    fuel.fill = uo2
    fuel.region = fuel_region

    gap = openmc.Cell(name='air gap')
    gap.region = gap_region

    clad = openmc.Cell(name='clad')
    clad.fill = zirconium
    clad.region = clad_region

    pitch = 1.26
    left = openmc.XPlane(-pitch/2, boundary_type='reflective')
    right = openmc.XPlane(pitch/2, boundary_type='reflective')
    bottom = openmc.YPlane(-pitch/2, boundary_type='reflective')
    top = openmc.YPlane(pitch/2, boundary_type='reflective')

    lower_left = (-pitch/2, -pitch/2, -pitch/2)
    upper_right = (pitch/2, pitch/2, pitch/2)

    water_region = +left & -right & +bottom & -top & +clad_outer_radius

    moderator = openmc.Cell(name='moderator')
    moderator.fill = water
    moderator.region = water_region
    box = openmc.rectangular_prism(width=pitch, height=pitch,
                               boundary_type='reflective')

    water_region = box & +clad_outer_radius

    root_universe = openmc.Universe(cells=(fuel, gap, clad, moderator))
    geometry = openmc.Geometry(root_universe)
    geometry.export_to_xml()
    volumes_calc = openmc.VolumeCalculation([fuel, gap, clad, moderator], 100000, lower_left, upper_right)
    settings = openmc.Settings()
    settings.volume_calculations = [volumes_calc]
    settings.export_to_xml()

    # This should not fail
    openmc.calculate_volumes()
