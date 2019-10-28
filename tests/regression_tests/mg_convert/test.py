import os
import hashlib

import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config

# OpenMC simulation parameters
batches = 10
inactive = 5
particles = 100


def build_mgxs_library(convert):
    # Instantiate the energy group data
    groups = openmc.mgxs.EnergyGroups(group_edges=[1e-5, 0.625, 20.0e6])

    # Instantiate the 2-group (C5G7) cross section data
    uo2_xsdata = openmc.XSdata('UO2', groups)
    uo2_xsdata.order = 2
    uo2_xsdata.set_total([2., 2.])
    uo2_xsdata.set_absorption([1., 1.])
    scatter_matrix = np.array([[[0.75, 0.25],
                                [0.00, 1.00]],
                               [[0.75 / 3., 0.25 / 3.],
                                [0.00 / 3., 1.00 / 3.]],
                               [[0.75 / 4., 0.25 / 4.],
                                [0.00 / 4., 1.00 / 4.]]])
    scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
    uo2_xsdata.set_scatter_matrix(scatter_matrix)
    uo2_xsdata.set_fission([0.5, 0.5])
    uo2_xsdata.set_nu_fission([1., 1.])
    uo2_xsdata.set_chi([1., 0.])

    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    mg_cross_sections_file.add_xsdatas([uo2_xsdata])

    if convert is not None:
        if isinstance(convert[0], list):
            for conv in convert:
                if conv[0] in ['legendre', 'tabular', 'histogram']:
                    mg_cross_sections_file = \
                        mg_cross_sections_file.convert_scatter_format(
                            conv[0], conv[1])
                elif conv[0] in ['angle', 'isotropic']:
                    mg_cross_sections_file = \
                        mg_cross_sections_file.convert_representation(
                            conv[0], conv[1], conv[1])
        elif convert[0] in ['legendre', 'tabular', 'histogram']:
            mg_cross_sections_file = \
                mg_cross_sections_file.convert_scatter_format(
                    convert[0], convert[1])
        elif convert[0] in ['angle', 'isotropic']:
            mg_cross_sections_file = \
                mg_cross_sections_file.convert_representation(
                    convert[0], convert[1], convert[1])

    mg_cross_sections_file.export_to_hdf5()


class MGXSTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Instantiate some Macroscopic Data
        uo2_data = openmc.Macroscopic('UO2')

        # Instantiate some Materials and register the appropriate objects
        mat = openmc.Material(material_id=1, name='UO2 fuel')
        mat.set_density('macro', 1.0)
        mat.add_macroscopic(uo2_data)

        # Instantiate a Materials collection and export to XML
        materials_file = openmc.Materials([mat])
        materials_file.cross_sections = "./mgxs.h5"
        materials_file.export_to_xml()

        # Instantiate ZCylinder surfaces
        left = openmc.XPlane(surface_id=4, x0=-5., name='left')
        right = openmc.XPlane(surface_id=5, x0=5., name='right')
        bottom = openmc.YPlane(surface_id=6, y0=-5., name='bottom')
        top = openmc.YPlane(surface_id=7, y0=5., name='top')

        left.boundary_type = 'reflective'
        right.boundary_type = 'vacuum'
        top.boundary_type = 'reflective'
        bottom.boundary_type = 'reflective'

        # Instantiate Cells
        fuel = openmc.Cell(cell_id=1, name='cell 1')

        # Use surface half-spaces to define regions
        fuel.region = +left & -right & +bottom & -top

        # Register Materials with Cells
        fuel.fill = mat

        # Instantiate Universe
        root = openmc.Universe(universe_id=0, name='root universe')

        # Register Cells with Universe
        root.add_cells([fuel])

        # Instantiate a Geometry, register the root Universe, and export to XML
        geometry = openmc.Geometry(root)
        geometry.export_to_xml()

        settings_file = openmc.Settings()
        settings_file.energy_mode = "multi-group"
        settings_file.batches = batches
        settings_file.inactive = inactive
        settings_file.particles = particles

        # Create an initial uniform spatial source distribution
        bounds = [-5, -5, -5, 5, 5, 5]
        uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
        settings_file.source = openmc.source.Source(space=uniform_dist)

        settings_file.export_to_xml()

    def _run_openmc(self):
        # Run multiple conversions to compare results
        cases = [['legendre', 2], ['legendre', 0],
                 ['tabular', 33], ['histogram', 32],
                 [['tabular', 33], ['legendre', 1]],
                 [['tabular', 33], ['tabular', 3]],
                 [['tabular', 33], ['histogram', 32]],
                 [['histogram', 32], ['legendre', 1]],
                 [['histogram', 32], ['tabular', 3]],
                 [['histogram', 32], ['histogram', 16]],
                 ['angle', 2], [['angle', 2], ['isotropic', None]]]

        outstr = ''
        for case in cases:
            build_mgxs_library(case)

            if config['mpi']:
                mpi_args = [config['mpiexec'], '-n', config['mpi_np']]
                openmc.run(openmc_exec=config['exe'], mpi_args=mpi_args)

            else:
                openmc.run(openmc_exec=config['exe'])

            with openmc.StatePoint('statepoint.{}.h5'.format(batches)) as sp:
                # Write out k-combined.
                outstr += 'k-combined:\n'
                form = '{:12.6E} {:12.6E}\n'
                outstr += form.format(sp.k_combined.n, sp.k_combined.s)

        return outstr

    def _get_results(self, outstr, hash_output=False):
        # Hash the results if necessary.
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr

    def _cleanup(self):
        super()._cleanup()
        f = os.path.join(os.getcwd(), 'mgxs.h5')
        if os.path.exists(f):
            os.remove(f)

    def execute_test(self):
        """Build input XMLs, run OpenMC, and verify correct results."""
        try:
            self._build_inputs()
            inputs = self._get_inputs()
            self._write_inputs(inputs)
            self._compare_inputs()
            outstr = self._run_openmc()
            results = self._get_results(outstr)
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def update_results(self):
        """Update results_true.dat and inputs_true.dat"""
        try:
            self._build_inputs()
            inputs = self._get_inputs()
            self._write_inputs(inputs)
            self._overwrite_inputs()
            outstr = self._run_openmc()
            results = self._get_results(outstr)
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()


def test_mg_convert():
    harness = MGXSTestHarness('statepoint.10.h5')
    harness.main()
