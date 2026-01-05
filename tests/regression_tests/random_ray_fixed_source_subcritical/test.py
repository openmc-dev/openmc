import os

import openmc
from openmc.examples import random_ray_lattice
from openmc.utility_funcs import change_directory
import pytest

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


@pytest.mark.parametrize("shape", ["flat", "linear_xy"])
def test_random_ray_fixed_source_subcritical(shape):
    with change_directory(shape):
        openmc.reset_auto_ids()

        # The general strategy is to reuse the random_ray_lattice model,
        # but redfine some of the geometry to make it a good
        # subcritical multiplication problem. We then also add in
        # a fixed source term.

        model = random_ray_lattice()

        # Begin by updating the random ray settings for fixed source
        settings = model.settings
        settings.random_ray['source_shape'] = shape
        settings.run_mode = 'fixed source'
        settings.particles = 30
        settings.random_ray['distance_active'] = 40.0
        settings.random_ray['distance_inactive'] = 40.0
        settings.random_ray['volume_normalized_flux_tallies'] = False

        # This problem needs about 2k iterations to converge,
        # but for regression testing we only need a few hundred
        # to ensure things are working as expected. With
        # only 100 inactive batches, tallies will still be off
        # by 3x or more. For validation against MGMC, be sure
        # to increase the batch counts.
        settings.batches = 125
        settings.inactive = 100

        ########################################
        # Define the alternative geometry

        pitch = 1.26

        for material in model.materials:
            if material.name == 'Water':
                water = material

        # The new geometry replaces two of the fuel pins with
        # moderator, reducing k-eff to around 0.84. We also
        # add a special universe in the corner of one of the moderator
        # regions to use as a domain constraint for the source
        moderator_infinite = openmc.Cell(fill=water, name='moderator infinite')
        mu = openmc.Universe(cells=[moderator_infinite])

        moderator_infinite2 = openmc.Cell(fill=water, name='moderator infinite 2')
        mu2 = openmc.Universe(cells=[moderator_infinite2])

        n_sub = 10

        lattice = openmc.RectLattice()
        lattice.lower_left = [-pitch/2.0, -pitch/2.0]
        lattice.pitch = [pitch/n_sub, pitch/n_sub]
        lattice.universes = [[mu] * n_sub for _ in range(n_sub)]

        lattice2 = openmc.RectLattice()
        lattice2.lower_left = [-pitch/2.0, -pitch/2.0]
        lattice2.pitch = [pitch/n_sub, pitch/n_sub]
        lattice2.universes = [[mu] * n_sub for _ in range(n_sub)]
        lattice2.universes[n_sub-1][n_sub-1] = mu2

        mod_lattice_cell = openmc.Cell(fill=lattice)
        mod_lattice_uni = openmc.Universe(cells=[mod_lattice_cell])

        mod_lattice_cell2 = openmc.Cell(fill=lattice2)
        mod_lattice_uni2 = openmc.Universe(cells=[mod_lattice_cell2])

        lattice2x2 = openmc.RectLattice()
        lattice2x2.lower_left = [-pitch, -pitch]
        lattice2x2.pitch = [pitch, pitch]

        universes = model.geometry.get_all_universes()
        for universe in universes.values():
            if universe.name == 'pincell':
                pincell = universe

        lattice2x2.universes = [
            [pincell, mod_lattice_uni],
            [mod_lattice_uni, mod_lattice_uni2]
        ]

        box = openmc.model.RectangularPrism(
            pitch*2, pitch*2, boundary_type='reflective')

        assembly = openmc.Cell(fill=lattice2x2, region=-box, name='assembly')

        root = openmc.Universe(name='root universe', cells=[assembly])
        model.geometry = openmc.Geometry(root)

        ########################################
        # Define the fixed source term

        s = 1.0 / 7.0
        strengths = [s, s, s, s, s, s, s]
        midpoints = [2.0e-5, 0.0735, 20.0, 2.0e2, 2.0e3, 0.75e6, 2.0e6]
        energy_distribution = openmc.stats.Discrete(x=midpoints, p=strengths)

        lower_left_src = [pitch - pitch/10.0, -pitch, -1.0]
        upper_right_src = [pitch, -pitch + pitch/10.0, 1.0]
        spatial_distribution = openmc.stats.Box(
            lower_left_src, upper_right_src, only_fissionable=False)

        settings.source = openmc.IndependentSource(
            space=spatial_distribution,
            energy=energy_distribution,
            constraints={'domains': [mu2]},
            strength=1.0
        )

        ########################################
        # Run test

        harness = MGXSTestHarness('statepoint.125.h5', model)
        harness.main()
