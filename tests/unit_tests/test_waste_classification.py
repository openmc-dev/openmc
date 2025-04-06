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
    assert mat.waste_classification(method='NRC') == 'Class B'


def test_waste_classification_fetter():
    """Test waste classification using the Fetter limits"""
    # For Tc99, Fetter has a more strict limit. Here, we create a material with
    # Tc99 at 1 Ci/m3 which exceeds Fetter but not NRC
    density = 3.5561e-7
    mat = openmc.Material()
    mat.add_nuclide('Tc99', density)
    assert mat.get_activity('Ci/m3') == pytest.approx(1.0, 1e-3)
    assert mat.waste_classification(method='NRC') == 'Class C'
    assert mat.waste_classification(method='Fetter') == 'GTCC'

    # With a lower density, it should be Class C under Fetter limits and Class A
    # under NRC limits
    mat = openmc.Material()
    mat.add_nuclide('Tc99', 5.0e-2*density)
    assert mat.waste_classification(method='NRC') == 'Class A'
    assert mat.waste_classification(method='Fetter') == 'Class C'


def test_waste_classification_limits():
    """Test override of specific activity limits"""
    mat = openmc.Material()
    mat.add_nuclide('K40', random.random())

    # Check for correct classification based on actual activity
    ci_m3 = mat.get_activity('Ci/m3')
    assert mat.waste_classification('Fetter', limits={'K40': 2*ci_m3}) == 'Class C'
    assert mat.waste_classification('Fetter', limits={'K40': 0.5*ci_m3}) == 'GTCC'
