import random

import openmc
import pytest


@pytest.mark.parametrize("metal", [False, True])
def test_waste_classification_long(metal):
    """Test classification when determined by long-lived radionuclides"""
    f = 10.0 if metal else 1.0
    limit = 8.0*f
    mat = openmc.Material()
    mat.add_nuclide('C14', 1e-9*f)
    assert mat.get_activity('Ci/m3') < 0.1 * limit
    assert mat.waste_classification(metal=metal) == 'Class A'

    mat = openmc.Material()
    mat.add_nuclide('C14', 1e-8*f)
    assert 0.1 * limit < mat.get_activity('Ci/m3') < limit
    assert mat.waste_classification(metal=metal) == 'Class C'

    mat = openmc.Material()
    mat.add_nuclide('C14', 1e-7*f)
    assert mat.get_activity('Ci/m3') > limit
    assert mat.waste_classification(metal=metal) == 'GTCC'


@pytest.mark.parametrize("metal", [False, True])
def test_waste_classification_short(metal):
    """Test classification when determined by short-lived radionuclides"""
    f = 10.0 if metal else 1.0
    col1, col2, col3 = 3.5*f, 70.0*f, 700.0*f

    mat = openmc.Material()
    mat.add_nuclide('Ni63', 1e-10*f)
    assert mat.get_activity('Ci/m3') < col1
    assert mat.waste_classification(metal=metal) == 'Class A'

    mat = openmc.Material()
    mat.add_nuclide('Ni63', 1e-10*10*f)
    assert col1 < mat.get_activity('Ci/m3') < col2
    assert mat.waste_classification(metal=metal) == 'Class B'

    mat = openmc.Material()
    mat.add_nuclide('Ni63', 1e-10*200*f)
    assert col2 < mat.get_activity('Ci/m3') < col3
    assert mat.waste_classification(metal=metal) == 'Class C'

    mat = openmc.Material()
    mat.add_nuclide('Ni63', 1e-10*2000*f)
    assert mat.get_activity('Ci/m3') > col3
    assert mat.waste_classification(metal=metal) == 'GTCC'


def test_waste_classification_mix():
    """Test classification when determined by a mix of radionuclides"""
    # Check example from 10 CFR 61.55 with mix of Sr90 and Cs137
    mat = openmc.Material()
    mat.add_nuclide('Sr90', 2.425e-9)
    mat.add_nuclide('Cs137', 1.115e-9)

    # In example, activity of Sr90 is 50.0 Ci/m3 and Cs137 is 22.0 Ci/m3
    activity = mat.get_activity(units='Ci/m3', by_nuclide=True)
    assert activity['Sr90'] == pytest.approx(50.0, 0.01)
    assert activity['Cs137'] == pytest.approx(22.0, 0.01)

    # According to example, the waste should be class B
    assert mat.waste_classification() == 'Class B'


def test_waste_rating_fetter():
    """Test waste classification using the Fetter limits"""
    # For Tc99, Fetter has a more strict limit. Here, we create a material with
    # Tc99 at 1 Ci/m3 which exceeds Fetter but not NRC
    density = 3.5561e-7
    mat = openmc.Material()
    mat.add_nuclide('Tc99', density)
    assert mat.get_activity('Ci/m3') == pytest.approx(1.0, 1e-3)
    assert mat.waste_disposal_rating(limits='NRC_short_C') < 1.0
    assert mat.waste_disposal_rating(limits='Fetter') > 1.0

    # With a lower density, it should be Class C under Fetter limits and Class A
    # under NRC limits
    mat = openmc.Material()
    mat.add_nuclide('Tc99', 5.0e-2*density)
    assert mat.waste_disposal_rating(limits='NRC_short_A') < 1.0
    assert mat.waste_disposal_rating(limits='Fetter') < 1.0


def test_waste_disposal_rating():
    """Test waste_disposal_rating method"""
    mat = openmc.Material()
    mat.add_nuclide('K40', random.random())

    # Check for correct classification based on actual activity
    ci_m3 = mat.get_activity('Ci/m3')
    assert mat.waste_disposal_rating(limits={'K40': 2*ci_m3}) < 1.0
    assert mat.waste_disposal_rating(limits={'K40': 0.5*ci_m3}) > 1.0

    wdr = mat.waste_disposal_rating(limits={'K40': 4*ci_m3}, by_nuclide=True)
    assert isinstance(wdr, dict)
    assert wdr['K40'] == pytest.approx(1/4)


def test_strlschv_unrestricted():
    """Test German StrlSchV unrestricted clearance rating"""
    # Co-60 unrestricted clearance limit is 0.1 Bq/g.
    # Create a material with Co60 and check the rating against its activity.
    mat = openmc.Material()
    mat.add_nuclide('Co60', 1e-12)
    bq_g = mat.get_activity('Bq/g', by_nuclide=True)['Co60']

    # Rating should equal activity / limit
    expected = bq_g / 0.1
    rating = mat.waste_disposal_rating(limits='StrlSchV_unrestricted')
    assert rating == pytest.approx(expected, rel=1e-3)

    # by_nuclide should return dict with Co60 entry
    wdr = mat.waste_disposal_rating(
        limits='StrlSchV_unrestricted', by_nuclide=True
    )
    assert isinstance(wdr, dict)
    assert 'Co60' in wdr
    assert wdr['Co60'] == pytest.approx(expected, rel=1e-3)


def test_strlschv_metal_recycling():
    """Test German StrlSchV metal scrap recycling clearance rating"""
    # Co-60 metal recycling limit is 0.6 Bq/g, which is less restrictive
    # than the unrestricted clearance limit of 0.1 Bq/g.
    mat = openmc.Material()
    mat.add_nuclide('Co60', 1e-12)

    unrestricted = mat.waste_disposal_rating(limits='StrlSchV_unrestricted')
    metal = mat.waste_disposal_rating(limits='StrlSchV_metal_recycling')

    # Metal recycling limit is 6x higher, so rating should be 6x lower
    assert unrestricted == pytest.approx(6.0 * metal, rel=1e-3)


def test_strlschv_sum_rule():
    """Test that the sum-of-fractions rule works for StrlSchV limits"""
    # Create material with two nuclides, each at half their unrestricted limit
    # Co-60 limit = 0.1 Bq/g, Cs-137 limit = 0.1 Bq/g
    mat = openmc.Material()
    mat.add_nuclide('Co60', 1e-12)
    mat.add_nuclide('Cs137', 1e-12)

    wdr = mat.waste_disposal_rating(
        limits='StrlSchV_unrestricted', by_nuclide=True
    )
    assert 'Co60' in wdr
    assert 'Cs137' in wdr

    # Total rating should be the sum of individual contributions
    total = mat.waste_disposal_rating(limits='StrlSchV_unrestricted')
    assert total == pytest.approx(wdr['Co60'] + wdr['Cs137'], rel=1e-6)


@pytest.mark.parametrize("pathway", [
    'StrlSchV_unrestricted',
    'StrlSchV_metal_recycling',
    'StrlSchV_landfill_100',
    'StrlSchV_landfill_1000',
    'StrlSchV_incineration_100',
    'StrlSchV_incineration_1000',
    'StrlSchV_soil',
    'StrlSchV_rubble',
])
def test_strlschv_all_pathways(pathway):
    """Test that all StrlSchV pathways return valid ratings"""
    mat = openmc.Material()
    mat.add_nuclide('Co60', 1e-12)

    rating = mat.waste_disposal_rating(limits=pathway)
    assert isinstance(rating, float)
    assert rating > 0.0

    wdr = mat.waste_disposal_rating(limits=pathway, by_nuclide=True)
    assert isinstance(wdr, dict)
    assert rating == pytest.approx(sum(wdr.values()), rel=1e-6)
