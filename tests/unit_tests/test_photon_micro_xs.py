"""Photon micro cross section tallies are looked up per element.

Photon data does not distinguish isotopes, so the cache holding photon micro
cross sections is indexed by element while a tally's nuclide bins are indexed by
nuclide. Two isotopes of the same element present at the same atom density must
therefore score the same microscopic cross section, whatever that cross section
happens to be, which makes this check independent of the data library.
"""

import openmc
import pytest


def test_isotopes_of_one_element_score_alike(run_in_tmpdir):
    openmc.reset_auto_ids()

    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.add_nuclide('U238', 1.0)
    mat.set_density('g/cm3', 10.0)

    sphere = openmc.Sphere(r=5.0, boundary_type='vacuum')
    model = openmc.Model()
    model.geometry = openmc.Geometry([openmc.Cell(fill=mat, region=-sphere)])
    model.settings.run_mode = 'fixed source'
    model.settings.photon_transport = True
    model.settings.particles = 200
    model.settings.batches = 2
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1.0e6], [1.0]),
        particle='photon')

    tally = openmc.Tally()
    tally.nuclides = ['U235', 'U238']
    tally.scores = ['total']
    model.tallies = openmc.Tallies([tally])

    model.run(apply_tally_results=True)

    u235, u238 = tally.mean.ravel()
    assert u235 > 0.0
    assert u238 == pytest.approx(u235, rel=1e-12)


def test_nuclide_absent_from_the_material(run_in_tmpdir):
    """A tally may name a nuclide no material contains.

    Doing so adds it to the global nuclide list without adding an element, so
    its index reaches past the end of the element-indexed photon cache unless
    the lookup goes through the nuclide-to-element mapping.

    multiply_density has to be off for the score to say anything: an absent
    nuclide has an atom density of zero, which multiplies away whatever came
    out of the cache. With it off the score is the microscopic cross section
    itself. The element is still present in the material, so transport keeps
    its cache current and this does not exercise the refresh for an element
    the scoring region lacks. test_element_absent_from_the_scoring_region does.
    """
    openmc.reset_auto_ids()

    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.set_density('g/cm3', 10.0)

    sphere = openmc.Sphere(r=5.0, boundary_type='vacuum')
    model = openmc.Model()
    model.geometry = openmc.Geometry([openmc.Cell(fill=mat, region=-sphere)])
    model.settings.run_mode = 'fixed source'
    model.settings.photon_transport = True
    model.settings.particles = 200
    model.settings.batches = 2
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1.0e6], [1.0]),
        particle='photon')

    tally = openmc.Tally()
    tally.nuclides = ['U235', 'U238']  # U238 is not in the material
    tally.scores = ['total']
    tally.multiply_density = False
    model.tallies = openmc.Tallies([tally])

    model.run(apply_tally_results=True)

    # Same element as U235, so the same microscopic cross section
    u235, u238 = tally.mean.ravel()
    assert u235 > 0.0
    assert u238 == pytest.approx(u235, rel=1e-12)


@pytest.mark.parametrize('estimator', ['tracklength', 'collision'])
def test_element_absent_from_the_scoring_region(run_in_tmpdir, estimator):
    """The cache of an element absent from the scoring region is refreshed.

    In the test above the tallied nuclide shares its element with the material,
    so transport keeps that element's cache current and the refresh is never
    needed. Here the inner sphere holds no silicon at all. Uncollided photons
    crossing it are all at the source energy and have not yet reached the
    silicon shell, so nothing but the tally evaluates silicon for them. The
    microscopic total over the flux in the uncollided energy bin is then the
    silicon cross section at the source energy, and must match the same ratio
    in the shell, which does contain silicon.

    The collision estimator scores after the collision, when photoelectric
    absorption has already set the photon energy to zero, so the refresh has to
    use the pre-collision energy. The source is at 100 keV so that photoelectric
    absorption in the iron is common enough to catch that.
    """
    openmc.reset_auto_ids()

    iron = openmc.Material()
    iron.add_nuclide('Fe56', 1.0)
    iron.set_density('g/cm3', 0.1)

    silicon = openmc.Material()
    silicon.add_nuclide('Si28', 1.0)
    silicon.set_density('g/cm3', 0.1)

    inner = openmc.Sphere(r=5.0)
    outer = openmc.Sphere(r=10.0, boundary_type='vacuum')
    iron_cell = openmc.Cell(fill=iron, region=-inner)
    silicon_cell = openmc.Cell(fill=silicon, region=+inner & -outer)

    model = openmc.Model()
    model.geometry = openmc.Geometry([iron_cell, silicon_cell])
    model.settings.run_mode = 'fixed source'
    model.settings.photon_transport = True
    model.settings.particles = 200
    model.settings.batches = 2
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1.0e5], [1.0]),
        particle='photon')

    filters = [
        openmc.CellFilter([iron_cell, silicon_cell]),
        openmc.EnergyFilter([0.999e5, 1.001e5]),
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
    assert silicon_total > 0.0
    assert iron_total / iron_flux == pytest.approx(
        silicon_total / silicon_flux, rel=1e-10)
