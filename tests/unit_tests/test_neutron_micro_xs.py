"""Neutron micro cross sections of a nuclide absent from the scoring region.

A tally binned by nuclide with multiply_density off scores nuclides the
scoring region need not contain, and the cross sections of those are evaluated
at scoring time rather than by transport. This checks that they are evaluated
at the energy the score belongs to, by comparing against a region that does
contain the nuclide, which makes the check independent of the data library.
"""

import openmc
import pytest


@pytest.mark.parametrize('estimator', ['tracklength', 'collision'])
def test_nuclide_absent_from_the_scoring_region(run_in_tmpdir, estimator):
    """The cache of a nuclide absent from the scoring region is refreshed.

    The inner sphere holds no silicon at all, so nothing but the tally
    evaluates silicon there. The energy bin is narrow enough that only
    neutrons at the source energy fall in it, so the microscopic total over
    the flux in that bin is the silicon cross section at the source energy,
    and must match the same ratio in the shell, which does contain silicon.

    The collision estimator scores after the collision, when elastic
    scattering off iron has already changed the neutron energy, so the
    refresh has to use the pre-collision energy.
    """
    openmc.reset_auto_ids()

    iron = openmc.Material()
    iron.add_nuclide('Fe56', 1.0)
    iron.set_density('g/cm3', 7.87)

    silicon = openmc.Material()
    silicon.add_nuclide('Si28', 1.0)
    silicon.set_density('g/cm3', 2.33)

    inner = openmc.Sphere(r=5.0)
    outer = openmc.Sphere(r=10.0, boundary_type='vacuum')
    iron_cell = openmc.Cell(fill=iron, region=-inner)
    silicon_cell = openmc.Cell(fill=silicon, region=+inner & -outer)

    model = openmc.Model()
    model.geometry = openmc.Geometry([iron_cell, silicon_cell])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 200
    model.settings.batches = 2
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1.0e6], [1.0]))

    filters = [
        openmc.CellFilter([iron_cell, silicon_cell]),
        openmc.EnergyFilter([1.0e6*(1 - 1e-9), 1.0e6*(1 + 1e-9)]),
    ]
    flux = openmc.Tally()
    flux.filters = filters
    flux.scores = ['flux']
    flux.estimator = estimator

    total = openmc.Tally()
    total.filters = filters
    total.nuclides = ['Si28']
    total.scores = ['total']
    total.multiply_density = False
    total.estimator = estimator
    model.tallies = openmc.Tallies([flux, total])

    model.run(apply_tally_results=True)

    iron_flux, silicon_flux = flux.mean.ravel()
    iron_total, silicon_total = total.mean.ravel()
    assert iron_flux > 0.0
    assert silicon_total > 0.0
    assert iron_total / iron_flux == pytest.approx(
        silicon_total / silicon_flux, rel=1e-10)
