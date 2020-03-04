import openmc
from pytest import approx, raises

from openmc.data import NATURAL_ABUNDANCE, atomic_mass


def test_expand_no_enrichment():
    """ Expand Li in natural compositions"""
    lithium = openmc.Element('Li')

    # Verify the expansion into ATOMIC fraction against natural composition
    for isotope in lithium.expand(100.0, 'ao'):
        assert isotope[1] == approx(NATURAL_ABUNDANCE[isotope[0]] * 100.0)

    # Verify the expansion into WEIGHT fraction against natural composition
    natural = {'Li6': NATURAL_ABUNDANCE['Li6'] * atomic_mass('Li6'),
               'Li7': NATURAL_ABUNDANCE['Li7'] * atomic_mass('Li7')}
    li_am = sum(natural.values())
    for key in natural:
        natural[key] /= li_am

    for isotope in lithium.expand(100.0, 'wo'):
        assert isotope[1] == approx(natural[isotope[0]] * 100.0)


def test_expand_enrichment():
    """ Expand and verify enrichment of Li """
    lithium = openmc.Element('Li')

    # Verify the enrichment by atoms
    ref = {'Li6': 75.0, 'Li7': 25.0}
    for isotope in lithium.expand(100.0, 'ao', 25.0, 'Li7', 'ao'):
        assert isotope[1] == approx(ref[isotope[0]])

    # Verify the enrichment by weight
    for isotope in lithium.expand(100.0, 'wo', 25.0, 'Li7', 'wo'):
        assert isotope[1] == approx(ref[isotope[0]])


def test_expand_exceptions():
    """ Test that correct exceptions are raised for invalid input """

    # 1 Isotope Element
    with raises(ValueError):
        element = openmc.Element('Be')
        element.expand(70.0, 'ao', 4.0, 'Be9')

    # 3 Isotope Element
    with raises(ValueError):
        element = openmc.Element('Cr')
        element.expand(70.0, 'ao', 4.0, 'Cr52')

    # Non-present Enrichment Target
    with raises(ValueError):
        element = openmc.Element('H')
        element.expand(70.0, 'ao', 4.0, 'H4')

    # Enrichment Procedure for Uranium if not Uranium
    with raises(ValueError):
        element = openmc.Element('Li')
        element.expand(70.0, 'ao', 4.0)

    # Missing Enrichment Target
    with raises(ValueError):
        element = openmc.Element('Li')
        element.expand(70.0, 'ao', 4.0, enrichment_type='ao')

    # Invalid Enrichment Type Entry
    with raises(ValueError):
        element = openmc.Element('Li')
        element.expand(70.0, 'ao', 4.0, 'Li7', 'Grand Moff Tarkin')

    # Trying to enrich Uranium
    with raises(ValueError):
        element = openmc.Element('U')
        element.expand(70.0, 'ao', 4.0, 'U235', 'wo')

    # Trying to enrich Uranium with wrong enrichment_target
    with raises(ValueError):
        element = openmc.Element('U')
        element.expand(70.0, 'ao', 4.0, enrichment_type='ao')
