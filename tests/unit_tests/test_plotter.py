import numpy as np
import pytest

import openmc


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


@pytest.mark.parametrize("this", ["Li", "Li6"])
def test_calculate_cexs_with_nuclide_and_element(this):
    # single type (reaction)
    energy_grid, data = openmc.plotter.calculate_cexs(this=this, types=[205])

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 1
    assert len(data[0]) == len(energy_grid)

    # two types (reactions)
    energy_grid, data = openmc.plotter.calculate_cexs(this=this, types=[2, "elastic"])

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 2
    assert len(data[0]) == len(energy_grid)
    assert len(data[0]) == len(energy_grid)
    # reactions are both the same MT number 2 is elastic
    assert np.array_equal(data[0], data[1])


def test_calculate_cexs_with_materials(test_mat):
    energy_grid, data = openmc.plotter.calculate_cexs(this=test_mat, types=[205])

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 1
    assert len(data[0]) == len(energy_grid)


@pytest.mark.parametrize("this", ["Be", "Be9"])
def test_plot_xs(this):
    from matplotlib.figure import Figure

    assert isinstance(
        openmc.plot_xs({this: ["total", "elastic", 16, "(n,2n)"]}), Figure
    )


def test_plot_xs_mat(test_mat):
    from matplotlib.figure import Figure

    assert isinstance(openmc.plot_xs({test_mat: ["total"]}), Figure)


@pytest.mark.parametrize("units", ["eV", "keV", "MeV"])
def test_plot_xs_energy_axis(units):
    plot = openmc.plot_xs({"Be9": ["(n,2n)"]}, energy_axis_units=units)
    axis_text = plot.get_axes()[0].get_xaxis().get_label().get_text()
    assert axis_text == f"Energy [{units}]"


def test_plot_axes_labels():
    # just nuclides
    axis_label = openmc.plotter._get_yaxis_label(
        reactions={
            "Li6": [205],
            "Li7": [205],
        },
        divisor_types=False,
    )
    assert axis_label == "Microscopic Cross Section [b]"

    # just elements
    axis_label = openmc.plotter._get_yaxis_label(
        reactions={
            "Li": [205],
            "Be": [16],
        },
        divisor_types=False,
    )
    assert axis_label == "Microscopic Cross Section [b]"

    # mixed nuclide and element
    axis_label = openmc.plotter._get_yaxis_label(
        reactions={
            "Li": [205],
            "Li7": [205],
        },
        divisor_types=False,
    )
    assert axis_label == "Microscopic Cross Section [b]"

    axis_label = openmc.plotter._get_yaxis_label(
        reactions={
            "Li": ["heating", "heating-local"],
            "Li7": ["heating"],
            "Be": ["damage-energy"],
        },
        divisor_types=False,
    )
    assert axis_label == "Heating Cross Section [eV-barn]"

    with pytest.raises(TypeError):
        axis_label = openmc.plotter.plot_xs(
            reactions={"Li": ["heating", "heating-local"], "Be9": ["(n,2n)"]}
        )

    # just materials
    mat1 = openmc.Material()
    mat1.add_nuclide("Fe56", 1)
    mat1.set_density("g/cm3", 1)
    mat2 = openmc.Material()
    mat2.add_element("Fe", 1)
    mat2.add_nuclide("Fe55", 1)
    mat2.set_density("g/cm3", 1)
    axis_label = openmc.plotter._get_yaxis_label(
        reactions={
            mat1: [205],
            mat2: [16],
        },
        divisor_types=False,
    )
    assert axis_label == "Macroscopic Cross Section [1/cm]"

    # mixed materials and nuclides
    with pytest.raises(TypeError):
        openmc.plotter._get_yaxis_label(
            reactions={"Li6": [205], mat2: [16]}, divisor_types=False
        )

    # mixed materials and elements
    with pytest.raises(TypeError):
        openmc.plotter._get_yaxis_label(
            reactions={"Li": [205], mat2: [16]}, divisor_types=False
        )


def test_get_title():
    title = openmc.plotter._get_title(reactions={"Li": [205]})
    assert title == "Cross Section Plot For Li"
    title = openmc.plotter._get_title(reactions={"Li6": [205]})
    assert title == "Cross Section Plot For Li6"
    title = openmc.plotter._get_title(reactions={"Li6": [205], "Li7": [205]})
    assert title == "Cross Section Plot"

    mat1 = openmc.Material()
    mat1.add_nuclide("Fe56", 1)
    mat1.set_density("g/cm3", 1)
    mat1.name = "my_mat"
    title = openmc.plotter._get_title(reactions={mat1: [205]})
    assert title == "Cross Section Plot For my_mat"


def _any_photon_mt(element_symbol, cross_sections=None):
    """Return a photon MT that is guaranteed to exist for the given element
    in the configured cross sections library.
    """
    if cross_sections is None:
        cross_sections = openmc.config.get("cross_sections")

    library = openmc.data.DataLibrary.from_xml(cross_sections)
    lib = library.get_by_material(element_symbol, data_type="photon")
    if lib is None:
        raise RuntimeError(f"No photon library entry found for {element_symbol}")

    inc = openmc.data.IncidentPhoton.from_hdf5(lib["path"])
    # `reactions` is a dict keyed by MT
    return next(iter(inc.reactions.keys()))


@pytest.mark.parametrize("this", ["Be", "Be9"])
def test_calculate_cexs_photon_with_element_and_nuclide(this):
    mt = _any_photon_mt("Be")

    # Use a common photoatomic MT (total) and verify basic shape/types
    energy_grid, data = openmc.plotter.calculate_cexs(
        this=this, types=[mt], incident_particle="photon"
    )

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 1
    assert len(data[0]) == len(energy_grid)


def test_calculate_cexs_photon_requires_integer_mts():
    # Photon cross sections can only be requested with integer MT numbers
    with pytest.raises(TypeError):
        openmc.plotter.calculate_cexs(
            this="Be", types=["total"], incident_particle="photon"
        )

    with pytest.raises(TypeError):
        openmc.plotter.calculate_cexs(
            this="Be", types=[502, "elastic"], incident_particle="photon"
        )


def test_calculate_cexs_photon_with_material():
    mat = openmc.Material()
    mat.add_element("Be", 1.0, "ao")
    mat.set_density("g/cm3", 1.85)

    mt = _any_photon_mt("Be")

    energy_grid, data = openmc.plotter.calculate_cexs(
        this=mat, types=[mt], incident_particle="photon"
    )

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert len(energy_grid) > 1
    assert len(data) == 1
    assert len(data[0]) == len(energy_grid)


def _any_photon_mt(element_symbol="C", cross_sections=None):
    """Pick an MT that actually exists in the configured photon library."""
    if cross_sections is None:
        cross_sections = openmc.config.get("cross_sections")

    library = openmc.data.DataLibrary.from_xml(cross_sections)
    lib = library.get_by_material(element_symbol, data_type="photon")
    if lib is None:
        raise RuntimeError(f"No photon library entry found for {element_symbol}")

    inc = openmc.data.IncidentPhoton.from_hdf5(lib["path"])
    return next(iter(inc.reactions.keys()))


def test_calculate_cexs_photon_material_element_vs_explicit_natural_abundance():
    mt = _any_photon_mt("C")

    # Material 1: defined as a single element (uses natural abundance implicitly)
    mat_elem = openmc.Material()
    mat_elem.add_element("C", 1.0, "ao")
    mat_elem.set_density("g/cm3", 1.0)

    # Material 2: defined by explicitly specifying natural isotopic abundance
    # (values are standard natural abundances for carbon)
    mat_iso = openmc.Material()
    mat_iso.add_nuclide("C12", 0.988922, "ao")
    mat_iso.add_nuclide("C13", 0.011078, "ao")
    mat_iso.set_density("g/cm3", 1.0)

    E1, xs1 = openmc.plotter.calculate_cexs(
        this=mat_elem, types=[mt], incident_particle="photon"
    )
    E2, xs2 = openmc.plotter.calculate_cexs(
        this=mat_iso, types=[mt], incident_particle="photon"
    )

    assert isinstance(E1, np.ndarray)
    assert isinstance(E2, np.ndarray)
    assert isinstance(xs1, np.ndarray)
    assert isinstance(xs2, np.ndarray)

    assert len(E1) > 1
    assert len(E2) > 1
    assert len(xs1) == 1
    assert len(xs2) == 1
    assert len(xs1[0]) == len(E1)
    assert len(xs2[0]) == len(E2)

    # For photon data, isotopes map to the same element library, so these should match.
    assert np.array_equal(E1, E2)
    assert np.allclose(xs1[0], xs2[0], rtol=1e-12, atol=0.0)


def test_calculate_cexs_photon_missing_mt_fallback():
    # Use an MT that should never exist in photon data
    energy_grid, data = openmc.plotter.calculate_cexs(
        this="Be", types=[9999], incident_particle="photon"
    )

    assert isinstance(energy_grid, np.ndarray)
    assert isinstance(data, np.ndarray)
    assert np.allclose(energy_grid, [openmc.plotter._MIN_E, openmc.plotter._MAX_E])
    assert data.shape == (1, 2)
    assert np.allclose(data, 0.0)


def test_calculate_cexs_photon_total_attenuation_reference_values():
    """Check total photon interaction XS for Pb and V at two reference energies.

    Total interaction is approximated by summing MTs: 502, 504, 515, 517, 522.
    Reference mass attenuation data from NIST.
    """
    openmc.reset_auto_ids()

    # Total interaction channels (library must contain these to run the test)
    types = [501]
    energies = np.array([1.0e5, 1.0e6])  # eV
    # data from https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z23.html
    v_density = 6.11  # g/cm3
    v_ref = np.array(
        [
            2.877e-01,
            5.794e-02,
        ]
    )  # cm2/g
    v_expected = v_ref * v_density

    # data from https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z82.html
    pb_density = 11.35  # g/cm3
    pb_ref = np.array(
        [
            5.549e00,
            7.102e-02,
        ]
    )  # cm2/g
    pb_expected = pb_ref * pb_density

    def _run_element(symbol: str, density: float):
        mat = openmc.Material()
        mat.add_element(symbol, 1.0)
        mat.set_density("g/cm3", density)

        # Compute macroscopic total XS for the material
        e_grid, data = openmc.plotter.calculate_cexs(
            this=mat, types=types, incident_particle="photon"
        )
        xs_grid = data.sum(axis=0)
        tab = openmc.data.Tabulated1D(e_grid, xs_grid, [len(e_grid)], [5])
        xs_mat_eval = tab(energies)

        return xs_mat_eval

    try:
        pb_vals = _run_element("Pb", pb_density)
        v_vals = _run_element("V", v_density)
    except Exception:
        pytest.skip(
            "Pb or V photon data / required MTs not available in cross section library."
        )

    assert np.allclose(pb_vals, pb_expected, rtol=5e-3, atol=1e-8)
    assert np.allclose(v_vals, v_expected, rtol=5e-3, atol=1e-8)
