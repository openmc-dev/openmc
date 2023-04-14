import numpy as np
import openmc
import pytest
from matplotlib.figure import Figure


@pytest.fixture(scope='module')
def test_mat():
    mat_1 = openmc.Material()
    mat_1.add_element("H", 4.0, "ao")
    mat_1.add_element("O", 4.0, "ao")
    mat_1.add_element("C", 4.0, "ao")
    return mat_1


def test_calculate_cexs_elem_mat_sab(test_mat):
    """Checks that sab cross sections are included in the
    _calculate_cexs_elem_mat method and have the correct shape"""

    test_mat.add_s_alpha_beta("c_C6H6")
    test_mat.set_density("g/cm3", 0.865)

    energy_grid, data = openmc.plotter._calculate_cexs_elem_mat(
        test_mat,
        ["inelastic"],
        sab_name="c_C6H6",
    )

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 1
    assert len(data[0]) == len(energy_grid)


@pytest.mark.parametrize("this", ["Li", "Li6"])
def test_calculate_cexs_with_nuclide_and_element(this):
    # single type (reaction)
    energy_grid, data = openmc.plotter.calculate_cexs(
        this=this, types=[205]
    )

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 1
    assert len(data[0]) == len(energy_grid)

    # two types (reactions)
    energy_grid, data = openmc.plotter.calculate_cexs(
        this=this, types=[2, "elastic"]
    )

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 2
    assert len(data[0]) == len(energy_grid)
    assert len(data[0]) == len(energy_grid)
    # reactions are both the same MT number 2 is elastic
    assert np.array_equal(data[0], data[1])


def test_calculate_cexs_with_materials(test_mat):
    energy_grid, data = openmc.plotter.calculate_cexs(
        this=test_mat, types=[205]
    )

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 1
    assert len(data[0]) == len(energy_grid)


@pytest.mark.parametrize("this", ["Be", "Be9"])
def test_plot_xs(this):
    fig = openmc.plotter.plot_xs({this: ['total', 'elastic']})
    assert isinstance(fig, Figure)


def test_plot_xs_mat(test_mat):
    assert isinstance(openmc.plotter.plot_xs({test_mat: ['total']}), Figure)


def test_plot_axes_labels():

    # just nuclides
    axis_label = openmc.plotter._get_yaxis_label(
        this_and_types = {
            'Li6': [205],
            'Li7': [205],
        }, divisor_types=False
    )
    assert axis_label == 'Microscopic Cross Section [b]'

    # just elements
    axis_label = openmc.plotter._get_yaxis_label(
        this_and_types = {
            'Li': [205],
            'Be': [16],
        }, divisor_types=False
    )
    assert axis_label == 'Microscopic Cross Section [b]'

    # mixed nuclide and element
    axis_label = openmc.plotter._get_yaxis_label(
        this_and_types = {
            'Li': [205],
            'Li7': [205],
        }, divisor_types=False
    )
    assert axis_label == 'Microscopic Cross Section [b]'

    # just materials
    mat1 = openmc.Material()
    mat1.add_nuclide('Fe56', 1)
    mat1.set_density('g/cm3', 1)
    mat2 = openmc.Material()
    mat2.add_element('Fe', 1)
    mat2.add_nuclide('Fe55', 1)
    mat2.set_density('g/cm3', 1)
    axis_label = openmc.plotter._get_yaxis_label(
        this_and_types = {
            mat1: [205],
            mat2: [16],
        }, divisor_types=False
    )
    assert axis_label == 'Macroscopic Cross Section [1/cm]'

    # mixed materials and nuclides
    with pytest.raises(TypeError):
        openmc.plotter._get_yaxis_label(
                this_and_types = {
                    'Li6': [205],
                    mat2: [16],
                }, divisor_types=False
            )

    # mixed materials and elements
    with pytest.raises(TypeError):
        openmc.plotter._get_yaxis_label(
                this_and_types = {
                    'Li': [205],
                    mat2: [16],
                }, divisor_types=False
            )


def test_get_title():
    title = openmc.plotter._get_title(this_and_types = {'Li': [205]})
    assert title == 'Cross Section Plot For Li'
    title = openmc.plotter._get_title(this_and_types = {'Li6': [205]})
    assert title == 'Cross Section Plot For Li6'
    title = openmc.plotter._get_title(this_and_types = {
        'Li6': [205],
        'Li7': [205]
    })
    assert title == 'Cross Section Plot'

    mat1 = openmc.Material()
    mat1.add_nuclide('Fe56', 1)
    mat1.set_density('g/cm3', 1)
    mat1.name = 'my_mat'
    title = openmc.plotter._get_title(this_and_types = {mat1: [205]})
    assert title == 'Cross Section Plot For my_mat'