"""Photon micro cross section tallies are looked up per element.

Photon data does not distinguish isotopes, so the cache holding photon micro
cross sections is indexed by element while a tally's nuclide bins are indexed by
nuclide. Two isotopes of the same element present at the same atom density must
therefore score the same microscopic cross section, whatever that cross section
happens to be -- which makes this check independent of the data library.
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
    itself, and it is also the branch that refreshes the cache for a nuclide
    the scoring region does not contain.
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
