import hashlib
from pathlib import Path

import numpy as np
import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.Model()

    material = openmc.Material()
    material.set_density('g/cm3', 11.35)
    material.add_element('Pb', 1.0)
    model.materials.append(material)

    sphere = openmc.Sphere(r=1.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=material, region=-sphere)
    model.geometry = openmc.Geometry([cell])

    model.settings.run_mode = 'fixed source'
    model.settings.particles = 2000
    model.settings.batches = 2
    model.settings.photon_transport = True
    # Survival biasing applies only to neutrons and must not alter photon
    # analog scatter or absorption scores.
    model.settings.survival_biasing = True
    model.settings.source = openmc.IndependentSource(
        particle='photon',
        energy=openmc.stats.delta_function(1.0e6))

    return model


class PhotonMGXSTestHarness(PyAPITestHarness):
    """Run photon transport and verify photon MGXS post-processing."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.material = self._model.materials[0]

        groups = openmc.mgxs.EnergyGroups(
            group_edges=[1.0e3, 1.0e5, 5.0e5, 1.1e6])
        self.mgxs_lib = openmc.mgxs.Library(
            self._model.geometry,
            mgxs_types=['total', 'absorption', 'nu-scatter matrix'],
            particle_type='photon')
        self.mgxs_lib.energy_groups = groups
        self.mgxs_lib.correction = None
        self.mgxs_lib.domain_type = 'material'
        self.mgxs_lib.build_library()
        self.mgxs_lib.add_to_tallies(self._model.tallies, merge=False)

    def _run_openmc(self):
        super()._run_openmc()

        with openmc.StatePoint(self._sp_name) as statepoint:
            self.mgxs_lib.load_from_statepoint(statepoint)

            production = self.mgxs_lib.get_mgxs(
                self.material, 'nu-scatter matrix')
            absorption = self.mgxs_lib.get_mgxs(self.material, 'absorption')

            self.primary = production.tallies[
                'primary photon production'].mean.sum()
            self.secondary = production.tallies[
                'secondary photon production'].mean.sum()
            if not all(np.all(np.isfinite(tally.mean))
                       for tally in production.tallies.values()):
                raise AssertionError(
                    'Photon production tallies contain non-finite values')
            if self.secondary <= 0.0:
                raise AssertionError(
                    'Photon transport did not score any secondary photons')

            self.production_xs = production.get_xs()
            self.absorption_xs = absorption.get_xs()
            self.ce_flux = production.tallies['flux'].mean[:, 0, 0]

            # Exercise the XSdata conversion used by subsequent MG/RR
            # workflows.
            mg_library = self.mgxs_lib.create_mg_library()
            xsdata = mg_library.xsdatas[0]
            multiplicity = xsdata.multiplicity_matrix[0]
            if multiplicity is None:
                raise AssertionError(
                    'Photon production should create a multiplicity matrix')
            if not np.allclose(
                    xsdata.scatter_matrix[0][:, :, 0], self.production_xs):
                raise AssertionError('Photon production matrix was not '
                                     'exported')
            if not np.allclose(xsdata.absorption[0], self.absorption_xs):
                raise AssertionError('Photon absorption was changed during '
                                     'export')
            physical_scatter = xsdata.total[0] - xsdata.absorption[0]
            scatter = xsdata.scatter_matrix[0][:, :, 0] / multiplicity
            if not np.allclose(np.sum(scatter, axis=1), physical_scatter):
                raise AssertionError(
                    'Photon multiplicity does not recover physical scatter')

        self._model.mgxs_file, self._model.materials, \
            self._model.geometry = self.mgxs_lib.create_mg_mode()
        self._model.materials.cross_sections = './mgxs.h5'
        self._model.settings.energy_mode = 'multi-group'

        flux = openmc.Tally(name='multigroup photon flux')
        flux.filters = [
            openmc.EnergyFilter(self.mgxs_lib.energy_groups.group_edges),
            openmc.ParticleFilter('photon'),
        ]
        flux.scores = ['flux']
        self._model.tallies = openmc.Tallies([flux])
        self._model.export_to_model_xml()
        self._model.mgxs_file.export_to_hdf5()
        super()._run_openmc()

    def _get_results(self, hash_output=False):
        with openmc.StatePoint(self._sp_name) as statepoint:
            mg_flux = statepoint.get_tally(
                name='multigroup photon flux').mean[:, 0, 0]

        ce_spectrum = self.ce_flux / self.ce_flux.sum()
        mg_spectrum = mg_flux / mg_flux.sum()
        if not np.allclose(mg_spectrum, ce_spectrum, rtol=0.2, atol=0.02):
            raise AssertionError(
                'Multigroup photon spectrum does not reproduce the '
                'continuous-energy spectrum')

        output = [
            f'primary photons: {self.primary:.8e}',
            f'secondary photons: {self.secondary:.8e}',
            'production matrix:',
            np.array2string(self.production_xs, precision=8),
            'absorption:',
            np.array2string(self.absorption_xs, precision=8),
            'continuous-energy flux spectrum:',
            np.array2string(ce_spectrum, precision=8),
            'multigroup flux spectrum:',
            np.array2string(mg_spectrum, precision=8),
        ]
        output = '\n'.join(output) + '\n'

        if hash_output:
            digest = hashlib.sha512(output.encode('utf-8'))
            output = digest.hexdigest()
        return output

    def _cleanup(self):
        super()._cleanup()
        Path('mgxs.h5').unlink(missing_ok=True)

def test_photon_mgxs(model):
    harness = PhotonMGXSTestHarness('statepoint.2.h5', model)
    harness.main()
