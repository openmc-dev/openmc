import numpy as np
import openmc
import pytest
from matplotlib.figure import Figure


@pytest.fixture(scope="module")
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


@pytest.mark.parametrize("this,data_type", [("Li", "element"), ("Li6", "nuclide")])
def test_calculate_cexs_with_element(this, data_type):

    # single type (reaction)
    energy_grid, data = openmc.plotter.calculate_cexs(
        this=this, data_type=data_type, types=[205]
    )

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 1
    assert len(data[0]) == len(energy_grid)

    # two types (reaction)
    energy_grid, data = openmc.plotter.calculate_cexs(
        this=this, data_type=data_type, types=[2, "elastic"]
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
        this=test_mat, types=[205], data_type="material"
    )

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 1
    assert len(data[0]) == len(energy_grid)


@pytest.mark.parametrize(("this,data_type"), [("Be", "element"), ("Be9", "nuclide")])
def test_plot_xs(this, data_type):
    assert isinstance(
        openmc.plotter.plot_xs(this, data_type=data_type, types=["total"]), Figure
    )


def test_plot_xs_mat(test_mat):
    assert isinstance(
        openmc.plotter.plot_xs(test_mat, data_type="material", types=["total"]), Figure
    )
