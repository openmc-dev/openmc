#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
from testing_harness import TestHarness, PyAPITestHarness
from input_set import MGInputSet
import openmc

class MGNuclideInputSet(MGInputSet):
    def build_default_materials_and_geometry(self):
        # Define materials needed for 1D/1G slab problem
        uo2_data = openmc.Macroscopic('uo2_iso', '71c')
        uo2 = openmc.Material(name='UO2', material_id=1)
        uo2.set_density('macro', 1.0)
        uo2.add_macroscopic(uo2_data)

        clad_data = openmc.Macroscopic('clad_iso', '71c')
        clad = openmc.Material(name='Clad', material_id=2)
        clad.set_density('macro', 1.0)
        clad.add_macroscopic(clad_data)

        water_data = openmc.Macroscopic('lwtr_iso', '71c')
        water = openmc.Material(name='LWTR', material_id=3)
        water.set_density('macro', 1.0)
        water.add_macroscopic(water_data)

        # Define the materials file.
        self.materials.default_xs = '71c'
        self.materials += (uo2, clad, water)

        # Define surfaces.

        # Assembly/Problem Boundary
        left   = openmc.XPlane(x0=0.0, surface_id=200,
                               boundary_type='reflective')
        right  = openmc.XPlane(x0=10.0, surface_id=201,
                               boundary_type='reflective')
        bottom = openmc.YPlane(y0=0.0, surface_id=300,
                               boundary_type='reflective')
        top    = openmc.YPlane(y0=10.0, surface_id=301,
                               boundary_type='reflective')

        down   = openmc.ZPlane(z0=0.0, surface_id=0,
                               boundary_type='reflective')
        fuel_clad_intfc = openmc.ZPlane(z0=2.0, surface_id=1)
        clad_lwtr_intfc = openmc.ZPlane(z0=2.4, surface_id=2)
        up     = openmc.ZPlane(z0=5.0, surface_id=3,
                               boundary_type='reflective')

        # Define cells
        c1 = openmc.Cell(cell_id=1)
        c1.region = +left & -right & +bottom & -top & +down & -fuel_clad_intfc
        c1.fill = uo2
        c2 = openmc.Cell(cell_id=2)
        c2.region = +left & -right & +bottom & -top & +fuel_clad_intfc & -clad_lwtr_intfc
        c2.fill = clad
        c3 = openmc.Cell(cell_id=3)
        c3.region = +left & -right & +bottom & -top & +clad_lwtr_intfc & -up
        c3.fill = water

        # Define root universe.
        root = openmc.Universe(universe_id=0, name='root universe')

        root.add_cells((c1,c2,c3))

        # Define the geometry file.
        geometry = openmc.Geometry()
        geometry.root_universe = root

        self.geometry = geometry

class MGMaxOrderTestHarness(PyAPITestHarness):
    def __init__(self, statepoint_name, tallies_present, mg=False):
        PyAPITestHarness.__init__(self, statepoint_name, tallies_present)
        self._input_set = MGNuclideInputSet()

    def _build_inputs(self):
        # Set P1 scattering
        self._input_set.settings.max_order = 1
        # Call standard input build
        super(MGMaxOrderTestHarness, self)._build_inputs()

if __name__ == '__main__':
    harness = MGMaxOrderTestHarness('statepoint.10.*', False, mg=True)
    harness.main()
