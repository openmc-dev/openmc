from math import pi, cos

import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


class PinRadialAzimuthalTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        ###############################################################################
        # Create materials for the problem

        uo2 = openmc.Material(name='UO2 fuel at 2.4% wt enrichment')
        uo2.set_density('g/cm3', 10.29769)
        uo2.add_element('U', 1., enrichment=2.4)
        uo2.add_element('O', 2.)

        helium = openmc.Material(name='Helium for gap')
        helium.set_density('g/cm3', 0.001598)
        helium.add_element('He', 2.4044e-4)

        zircaloy = openmc.Material(name='Zircaloy 4')
        zircaloy.set_density('g/cm3', 6.55)
        zircaloy.add_element('Sn', 0.014  , 'wo')
        zircaloy.add_element('Fe', 0.00165, 'wo')
        zircaloy.add_element('Cr', 0.001  , 'wo')
        zircaloy.add_element('Zr', 0.98335, 'wo')

        borated_water = openmc.Material(name='Borated water')
        borated_water.set_density('g/cm3', 0.740582)
        borated_water.add_element('B', 4.0e-5)
        borated_water.add_element('H', 5.0e-2)
        borated_water.add_element('O', 2.4e-2)
        borated_water.add_s_alpha_beta('c_H_in_H2O')


        ###############################################################################
        # Define problem geometry

        # Create a region represented as the inside of a rectangular prism
        pitch = 1.25984
        box = openmc.rectangular_prism(pitch, pitch, boundary_type='reflective')

        cell1 = openmc.Cell(name='Cell 1')
        cell1.region = box

        # Create cylindrical surfaces
        fuel_or = openmc.ZCylinder(r=0.39218, name='Fuel OR')
        clad_ir = openmc.ZCylinder(r=0.40005, name='Clad IR')
        clad_or = openmc.ZCylinder(r=0.45720, name='Clad OR')
        corner =  openmc.ZCylinder(r=pitch/2.0, name='Clad OR')

        surfs = [fuel_or, clad_ir, clad_or, corner]
        mats = [uo2, helium, zircaloy, borated_water, borated_water.clone()]
        subdivs_r = {
                0 : 3,
                2 : 2,
                }
        subdivs_a = {
                0 : 4,
                2 : 4,
                3 : 8
                }
        rad_div_types = {
                0 : "area",
                2 : "area",
                3 : "radius"
        }

        pin_universe = openmc.model.pin_radial_azimuthal(surfs, mats, subdivisions_r=subdivs_r, subdivisions_a=subdivs_a, rad_div_types=rad_div_types, implicit_azi_div=8)
    
        cell1.fill = pin_universe

        root = openmc.Universe(name='root universe')
        root.add_cell(cell1)

        # Create a geometry
        self._model.geometry = openmc.Geometry(root)
        self._model.geometry.remove_redundant_surfaces()

        # Collect the materials together
        self._model.materials = openmc.Materials(mats)

        ###############################################################################
        # Define problem settings

        # Indicate how many particles to run
        settings = openmc.Settings()
        settings.batches = 100
        settings.inactive = 10
        settings.particles = 1000

        # Create an initial uniform spatial source distribution over fissionable zones
        lower_left = (-pitch/2, -pitch/2, -1)
        upper_right = (pitch/2, pitch/2, 1)
        uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=False)
        settings.source = openmc.source.Source(space=uniform_dist)

        self._model.settings = settings


def test_source():
    harness = PinRadialAzimuthalTestHarness('statepoint.100.h5', model=openmc.Model())
    harness.main()


test_source()
