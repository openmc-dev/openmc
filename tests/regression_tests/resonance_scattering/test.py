import openmc

from tests.testing_harness import PyAPITestHarness


class ResonanceScatteringTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Materials
        mat = openmc.Material(material_id=1)
        mat.set_density('g/cc', 1.0)
        mat.add_nuclide('U238', 1.0)
        mat.add_nuclide('U235', 0.02)
        mat.add_nuclide('Pu239', 0.02)
        mat.add_nuclide('H1', 20.0)

        self._model.materials = openmc.Materials([mat])

        # Geometry
        dumb_surface = openmc.XPlane(100, boundary_type='reflective')
        c1 = openmc.Cell(cell_id=1, fill=mat, region=-dumb_surface)
        root_univ = openmc.Universe(universe_id=0, cells=[c1])
        self._model.geometry = openmc.Geometry(root_univ)

        # Resonance elastic scattering settings
        res_scat_settings = {
            'enable': True,
            'energy_min': 1.0,
            'energy_max': 210.0,
            'method': 'rvs',
            'nuclides': ['U238', 'U235', 'Pu239']
        }

        settings = openmc.Settings()
        settings.batches = 10
        settings.inactive = 5
        settings.particles = 1000
        settings.source = openmc.IndependentSource(
             space=openmc.stats.Box([-4, -4, -4], [4, 4, 4]))
        settings.resonance_scattering = res_scat_settings
        self._model.settings = settings


def test_resonance_scattering():
    harness = ResonanceScatteringTestHarness('statepoint.10.h5',
                                             model=openmc.Model())
    harness.main()
