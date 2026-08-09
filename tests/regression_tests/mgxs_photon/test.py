import hashlib

import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


class PhotonMGXSTestHarness(PyAPITestHarness):
    """Run photon transport and verify photon MGXS post-processing."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.material = openmc.Material()
        self.material.set_density('g/cm3', 11.35)
        self.material.add_element('Pb', 1.0)
        self._model.materials = openmc.Materials([self.material])

        sphere = openmc.Sphere(r=1.0, boundary_type='vacuum')
        cell = openmc.Cell(fill=self.material, region=-sphere)
        self._model.geometry = openmc.Geometry([cell])

        self._model.settings.run_mode = 'fixed source'
        self._model.settings.particles = 2000
        self._model.settings.batches = 2
        self._model.settings.photon_transport = True
        self._model.settings.atomic_relaxation = True
        self._model.settings.electron_treatment = 'ttb'
        self._model.settings.cutoff = {'energy_photon': 1000.0}
        self._model.settings.source = openmc.IndependentSource(
            particle='photon',
            space=openmc.stats.Point((0.0, 0.0, 0.0)),
            energy=openmc.stats.Discrete([1.0e6], [1.0]))

        groups = openmc.mgxs.EnergyGroups(
            group_edges=[1.0e3, 1.0e5, 5.0e5, 1.1e6])
        self.mgxs_lib = openmc.mgxs.Library(
            self._model.geometry,
            mgxs_types=['total', 'absorption', 'photon-production matrix'],
            particle_type='photon')
        self.mgxs_lib.energy_groups = groups
        self.mgxs_lib.correction = None
        self.mgxs_lib.domain_type = 'material'
        self.mgxs_lib.build_library()
        self.mgxs_lib.add_to_tallies(self._model.tallies, merge=False)

    def _get_results(self, hash_output=False):
        with openmc.StatePoint(self._sp_name) as statepoint:
            self.mgxs_lib.load_from_statepoint(statepoint)

            production = self.mgxs_lib.get_mgxs(
                self.material, 'photon-production matrix')
            absorption = self.mgxs_lib.get_mgxs(self.material, 'absorption')

            primary = production.tallies[
                'primary photon production'].mean.sum()
            secondary = production.tallies[
                'secondary photon production'].mean.sum()
            if secondary <= 0.0:
                raise AssertionError(
                    'Photon transport did not score any secondary photons')

            production_xs = production.get_xs()
            absorption_xs = absorption.get_xs()

            # Exercise the XSdata conversion used by subsequent MG/RR
            # workflows.
            mg_library = self.mgxs_lib.create_mg_library()
            xsdata = mg_library.xsdatas[0]
            if xsdata.multiplicity_matrix[0] is not None:
                raise AssertionError(
                    'Photon production should not create a multiplicity '
                    'matrix')
            if not np.allclose(
                    xsdata.scatter_matrix[0][:, :, 0], production_xs):
                raise AssertionError('Photon production matrix was not '
                                     'exported')
            if not np.allclose(xsdata.absorption[0], absorption_xs):
                raise AssertionError('Photon absorption was changed during '
                                     'export')

        output = [
            f'primary photons: {primary:.8e}',
            f'secondary photons: {secondary:.8e}',
            'production matrix:',
            np.array2string(production_xs, precision=8),
            'absorption:',
            np.array2string(absorption_xs, precision=8),
        ]
        output = '\n'.join(output) + '\n'

        if hash_output:
            digest = hashlib.sha512(output.encode('utf-8'))
            output = digest.hexdigest()
        return output


def test_photon_mgxs():
    harness = PhotonMGXSTestHarness(
        'statepoint.2.h5', model=openmc.Model())
    harness.main()
