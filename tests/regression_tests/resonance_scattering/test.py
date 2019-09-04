import openmc

from tests.testing_harness import PyAPITestHarness


class ResonanceScatteringTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Materials
        mat = openmc.Material(material_id=1)
        mat.set_density('g/cc', 1.0)
        mat.add_nuclide('U238', 1.0)
        mat.add_nuclide('U235', 0.02)
        mat.add_nuclide('Pu239', 0.02)
        mat.add_nuclide('H1', 20.0)

        mats_file = openmc.Materials([mat])
        mats_file.export_to_xml()

        # Geometry
        dumb_surface = openmc.XPlane(100, 'reflective')
        c1 = openmc.Cell(cell_id=1, fill=mat, region=-dumb_surface)
        root_univ = openmc.Universe(universe_id=0, cells=[c1])
        geometry = openmc.Geometry(root_univ)
        geometry.export_to_xml()

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
        settings.source = openmc.source.Source(
             space=openmc.stats.Box([-4, -4, -4], [4, 4, 4]))
        settings.resonance_scattering = res_scat_settings
        settings.export_to_xml()


def test_resonance_scattering():
    harness = ResonanceScatteringTestHarness('statepoint.10.h5')
    harness.main()
