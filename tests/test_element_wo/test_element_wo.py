#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
from input_set import PinCellInputSet
import openmc
import openmc.mgxs


class ElementWOTestHarness(PyAPITestHarness):
    def _build_inputs(self):

        # Set the input set to use the pincell model
        self._input_set = PinCellInputSet()

        # Define materials.
        fuel = openmc.Material(name='Fuel')
        fuel.set_density('g/cm3', 10.29769)
        fuel.add_element("U", 0.88, 'wo')
        fuel.add_element("O", 0.12, 'wo')

        # Add some natural elements to the fuel
        for element in [ 'Ir', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Th', 'Pa']:
            fuel.add_element(element, 1.e-8, 'wo')

        clad = openmc.Material(name='Cladding')
        clad.set_density('g/cm3', 6.55)
        clad.add_element("Zr", 1.0, 'wo')

        hot_water = openmc.Material(name='Hot borated water')
        hot_water.set_density('g/cm3', 0.740582)
        hot_water.add_element("H",  2./18., 'wo')
        hot_water.add_element("O", 16./18., 'wo')
        hot_water.add_s_alpha_beta('c_H_in_H2O')

        # Define the materials file.
        self._input_set.materials += (fuel, clad, hot_water)

        # Instantiate ZCylinder surfaces
        fuel_or = openmc.ZCylinder(x0=0, y0=0, R=0.39218, name='Fuel OR')
        clad_or = openmc.ZCylinder(x0=0, y0=0, R=0.45720, name='Clad OR')
        left = openmc.XPlane(x0=-0.63, name='left')
        right = openmc.XPlane(x0=0.63, name='right')
        bottom = openmc.YPlane(y0=-0.63, name='bottom')
        top = openmc.YPlane(y0=0.63, name='top')

        left.boundary_type = 'reflective'
        right.boundary_type = 'reflective'
        top.boundary_type = 'reflective'
        bottom.boundary_type = 'reflective'

        # Instantiate Cells
        fuel_pin = openmc.Cell(name='cell 1')
        cladding = openmc.Cell(name='cell 3')
        water = openmc.Cell(name='cell 2')

        # Use surface half-spaces to define regions
        fuel_pin.region = -fuel_or
        cladding.region = +fuel_or & -clad_or
        water.region = +clad_or & +left & -right & +bottom & -top

        # Register Materials with Cells
        fuel_pin.fill = fuel
        cladding.fill = clad
        water.fill = hot_water

        # Instantiate Universe
        root = openmc.Universe(universe_id=0, name='root universe')

        # Register Cells with Universe
        root.add_cells([fuel_pin, cladding, water])

        # Instantiate a Geometry, register the root Universe, and export to XML
        self._input_set.geometry.root_universe = root


        mat_filter = openmc.MaterialFilter((fuel.id,))
        flux_tally = openmc.Tally()
        flux_tally.filters = [mat_filter]
        flux_tally.scores = ['flux']

        self._input_set.tallies = openmc.Tallies()
        self._input_set.tallies += [flux_tally]
        self._input_set.build_default_settings()
        self._input_set.export()

    def _cleanup(self):
        super(ElementWOTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'tallies.xml')
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    harness = ElementWOTestHarness('statepoint.10.*', True)
    harness.main()
