import openmc

from tests.testing_harness import PyAPITestHarness


class EnergyCutoffTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set energy cutoff
        energy_cutoff = 4.0

        # Material is composed of H-1
        mat = openmc.Material(material_id=1, name='mat')
        mat.set_density('atom/b-cm', 0.069335)
        mat.add_nuclide('H1', 40.0)
        self._model.materials = openmc.Materials([mat])

        # Cell is box with reflective boundary
        x1 = openmc.XPlane(surface_id=1, x0=-1)
        x2 = openmc.XPlane(surface_id=2, x0=1)
        y1 = openmc.YPlane(surface_id=3, y0=-1)
        y2 = openmc.YPlane(surface_id=4, y0=1)
        z1 = openmc.ZPlane(surface_id=5, z0=-1)
        z2 = openmc.ZPlane(surface_id=6, z0=1)
        for surface in [x1, x2, y1, y2, z1, z2]:
            surface.boundary_type = 'reflective'
        box = openmc.Cell(cell_id=1, name='box')
        box.region = +x1 & -x2 & +y1 & -y2 & +z1 & -z2
        box.fill = mat
        root = openmc.Universe(universe_id=0, name='root universe')
        root.add_cell(box)
        self._model.geometry = openmc.Geometry(root)

        # Set the running parameters
        settings_file = openmc.Settings()
        settings_file.run_mode = 'fixed source'
        settings_file.batches = 10
        settings_file.particles = 100
        settings_file.cutoff = {'energy_neutron': energy_cutoff}
        bounds = [-1, -1, -1, 1, 1, 1]
        uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
        watt_dist = openmc.stats.Watt()
        settings_file.source = openmc.IndependentSource(space=uniform_dist,
                                                        energy=watt_dist)
        self._model.settings = settings_file

        # Tally flux under energy cutoff
        tallies = openmc.Tallies()
        tally = openmc.Tally(1)
        tally.scores = ['flux']
        energy_filter = openmc.filter.EnergyFilter((0.0, energy_cutoff))
        tally.filters = [energy_filter]
        tallies.append(tally)
        self._model.tallies = tallies

    def _get_results(self):
        """Digest info in the statepoint and return as a string."""
        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        # Write out tally data.
        outstr = ''
        t = sp.get_tally()
        outstr += 'tally {}:\n'.format(t.id)
        outstr += 'sum = {:12.6E}\n'.format(t.sum[0, 0, 0])
        outstr += 'sum_sq = {:12.6E}\n'.format(t.sum_sq[0, 0, 0])

        return outstr


def test_energy_cutoff():
    harness = EnergyCutoffTestHarness('statepoint.10.h5', model=openmc.Model())
    harness.main()
