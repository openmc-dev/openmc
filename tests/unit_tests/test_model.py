from math import pi
from pathlib import Path
import os

import numpy as np
import pytest

import openmc
import openmc.lib


@pytest.fixture(scope='function')
def pin_model_attributes():
    uo2 = openmc.Material(material_id=1, name='UO2')
    uo2.set_density('g/cm3', 10.29769)
    uo2.add_element('U', 1., enrichment=2.4)
    uo2.add_element('O', 2.)
    uo2.depletable = True

    zirc = openmc.Material(material_id=2, name='Zirc')
    zirc.set_density('g/cm3', 6.55)
    zirc.add_element('Zr', 1.)
    zirc.depletable = False

    borated_water = openmc.Material(material_id=3, name='Borated water')
    borated_water.set_density('g/cm3', 0.740582)
    borated_water.add_element('B', 4.0e-5)
    borated_water.add_element('H', 5.0e-2)
    borated_water.add_element('O', 2.4e-2)
    borated_water.add_s_alpha_beta('c_H_in_H2O')
    borated_water.depletable = False

    mats = openmc.Materials([uo2, zirc, borated_water])

    pitch = 1.25984
    fuel_or = openmc.ZCylinder(r=0.39218, name='Fuel OR')
    clad_or = openmc.ZCylinder(r=0.45720, name='Clad OR')
    box = openmc.model.RectangularPrism(pitch, pitch,
                                        boundary_type='reflective')

    # Define cells
    fuel_inf_cell = openmc.Cell(cell_id=1, name='inf fuel', fill=uo2)
    fuel_inf_univ = openmc.Universe(universe_id=1, cells=[fuel_inf_cell])
    fuel = openmc.Cell(cell_id=2, name='fuel',
                       fill=fuel_inf_univ, region=-fuel_or)
    clad = openmc.Cell(cell_id=3, fill=zirc, region=+fuel_or & -clad_or)
    water = openmc.Cell(cell_id=4, fill=borated_water, region=+clad_or & -box)

    # Define overall geometry
    geom = openmc.Geometry([fuel, clad, water])
    uo2.volume = pi * fuel_or.r**2

    settings = openmc.Settings()
    settings.batches = 100
    settings.inactive = 10
    settings.particles = 1000

    # Create a uniform spatial source distribution over fissionable zones
    bounds = [-0.62992, -0.62992, -1, 0.62992, 0.62992, 1]
    uniform_dist = openmc.stats.Box(
        bounds[:3], bounds[3:], only_fissionable=True)
    settings.source = openmc.IndependentSource(space=uniform_dist)

    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.lower_left = [-0.39218, -0.39218, -1.e50]
    entropy_mesh.upper_right = [0.39218, 0.39218, 1.e50]
    entropy_mesh.dimension = [10, 10, 1]
    settings.entropy_mesh = entropy_mesh

    tals = openmc.Tallies()
    tal = openmc.Tally(tally_id=1, name='test')
    tal.filters = [openmc.MaterialFilter(bins=[uo2])]
    tal.scores = ['flux', 'fission']
    tals.append(tal)

    plot1 = openmc.Plot(plot_id=1)
    plot1.origin = (0., 0., 0.)
    plot1.width = (pitch, pitch)
    plot1.pixels = (300, 300)
    plot1.color_by = 'material'
    plot1.filename = 'test'
    plot2 = openmc.Plot(plot_id=2)
    plot2.origin = (0., 0., 0.)
    plot2.width = (pitch, pitch)
    plot2.pixels = (300, 300)
    plot2.color_by = 'cell'
    plots = openmc.Plots((plot1, plot2))

    chain = './test_chain.xml'

    chain_file_xml = """<?xml version="1.0"?>
<depletion_chain>
  <nuclide name="Xe136" decay_modes="0" reactions="0" />
  <nuclide name="U235" decay_modes="0" reactions="1">
    <reaction type="fission" Q="200000000."/>
    <neutron_fission_yields>
      <energies>2.53000e-02</energies>
      <fission_yields energy="2.53000e-02">
        <products>Xe136</products>
        <data>1.0</data>
      </fission_yields>
    </neutron_fission_yields>
  </nuclide>
</depletion_chain>
"""
    operator_kwargs = {'chain_file': chain}

    return (mats, geom, settings, tals, plots, operator_kwargs, chain_file_xml)


def test_init(run_in_tmpdir, pin_model_attributes, mpi_intracomm):
    mats, geom, settings, tals, plots, _, _ = \
        pin_model_attributes

    openmc.reset_auto_ids()
    # Check blank initialization of a model
    test_model = openmc.Model()
    assert test_model.geometry.root_universe is None
    assert len(test_model.materials) == 0
    ref_settings = openmc.Settings()
    assert sorted(test_model.settings.__dict__.keys()) == \
        sorted(ref_settings.__dict__.keys())
    for ref_k, ref_v in ref_settings.__dict__.items():
        assert test_model.settings.__dict__[ref_k] == ref_v
    assert len(test_model.tallies) == 0
    assert len(test_model.plots) == 0
    assert test_model._materials_by_id == {}
    assert test_model._materials_by_name == {}
    assert test_model._cells_by_id == {}
    assert test_model._cells_by_name == {}
    assert test_model.is_initialized is False

    # Now check proper init of an actual model. Assume no interference between
    # parameters and so we can apply them all at once instead of testing one
    # parameter initialization at a time
    test_model = openmc.Model(geom, mats, settings, tals, plots)
    assert test_model.geometry is geom
    assert test_model.materials is mats
    assert test_model.settings is settings
    assert test_model.tallies is tals
    assert test_model.plots is plots
    assert test_model._materials_by_id == {1: mats[0], 2: mats[1], 3: mats[2]}
    assert test_model._materials_by_name == {
        'UO2': {mats[0]}, 'Zirc': {mats[1]}, 'Borated water': {mats[2]}}
    # The last cell is the one that contains the infinite fuel
    assert test_model._cells_by_id == \
        {2: geom.root_universe.cells[2], 3: geom.root_universe.cells[3],
         4: geom.root_universe.cells[4],
         1: geom.root_universe.cells[2].fill.cells[1]}
    # No cell name for 2 and 3, so we expect a blank name to be assigned to
    # cell 3 due to overwriting
    assert test_model._cells_by_name == {
        'fuel': {geom.root_universe.cells[2]},
        '': {geom.root_universe.cells[3], geom.root_universe.cells[4]},
        'inf fuel': {geom.root_universe.cells[2].fill.cells[1]}}
    assert test_model.is_initialized is False

    # Finally test the parameter type checking by passing bad types and
    # obtaining the right exception types
    def_params = [geom, mats, settings, tals, plots]
    for i in range(len(def_params)):
        args = def_params.copy()
        # Try an integer, as that is a bad type for all arguments
        args[i] = i
        with pytest.raises(TypeError):
            test_model = openmc.Model(*args)


def test_from_xml(run_in_tmpdir, pin_model_attributes):
    mats, geom, settings, tals, plots, _, _ = pin_model_attributes

    # This test will write the individual files to xml and then init that way
    # and run the same sort of test as in test_init
    mats.export_to_xml()
    geom.export_to_xml()
    settings.export_to_xml()
    tals.export_to_xml()
    plots.export_to_xml()

    # This from_xml method cannot load chain and fission_q
    test_model = openmc.Model.from_xml()
    assert test_model.geometry.root_universe.cells.keys() == \
        geom.root_universe.cells.keys()
    assert [c.fill.name for c in
            test_model.geometry.root_universe.cells.values()] == \
        [c.fill.name for c in geom.root_universe.cells.values()]
    assert [mat.name for mat in test_model.materials] == \
        [mat.name for mat in mats]
    # We will assume the attributes of settings that are custom objects are
    # OK if the others are so we dotn need to implement explicit comparisons
    no_test = ['_source', '_entropy_mesh']
    assert sorted(k for k in test_model.settings.__dict__.keys()
                  if k not in no_test) == \
        sorted(k for k in settings.__dict__.keys() if k not in no_test)
    keys = sorted(k for k in settings.__dict__.keys() if k not in no_test)
    for ref_k in keys:
        assert test_model.settings.__dict__[ref_k] == settings.__dict__[ref_k]
    assert len(test_model.tallies) == 1
    assert len(test_model.plots) == 2
    assert test_model._materials_by_id == \
        {1: test_model.materials[0], 2: test_model.materials[1],
         3: test_model.materials[2]}
    assert test_model._materials_by_name == {
        'UO2': {test_model.materials[0]}, 'Zirc': {test_model.materials[1]},
        'Borated water': {test_model.materials[2]}}
    assert test_model._cells_by_id == {
        2: test_model.geometry.root_universe.cells[2],
        3: test_model.geometry.root_universe.cells[3],
        4: test_model.geometry.root_universe.cells[4],
        1: test_model.geometry.root_universe.cells[2].fill.cells[1]}
    # No cell name for 2 and 3, so we expect a blank name to be assigned to
    # cell 3 due to overwriting
    assert test_model._cells_by_name == {
        'fuel': {test_model.geometry.root_universe.cells[2]},
        '': {test_model.geometry.root_universe.cells[3],
             test_model.geometry.root_universe.cells[4]},
        'inf fuel': {test_model.geometry.root_universe.cells[2].fill.cells[1]}}
    assert test_model.is_initialized is False


def test_init_finalize_lib(run_in_tmpdir, pin_model_attributes, mpi_intracomm):
    # We are going to init and then make sure data is loaded
    mats, geom, settings, tals, plots, _, _ = pin_model_attributes
    test_model = openmc.Model(geom, mats, settings, tals, plots)
    test_model.init_lib(output=False, intracomm=mpi_intracomm)

    # First check that the API is advertised as initialized
    assert openmc.lib.is_initialized is True
    assert test_model.is_initialized is True
    # Now make sure it actually is initialized by making a call to the lib
    c_mat = openmc.lib.find_material((0.6, 0., 0.))
    # This should be Borated water
    assert c_mat.name == 'Borated water'
    assert c_mat.id == 3

    # Ok, now lets test that we can clear the data and check that it is cleared
    test_model.finalize_lib()

    # First check that the API is advertised as initialized
    assert openmc.lib.is_initialized is False
    assert test_model.is_initialized is False
    # Note we cant actually test that a sys call fails because we should get a
    # seg fault


def test_import_properties(run_in_tmpdir, mpi_intracomm):
    """Test importing properties on the Model class """

    # Create PWR pin cell model and write XML files
    openmc.reset_auto_ids()
    model = openmc.examples.pwr_pin_cell()
    model.init_lib(output=False, intracomm=mpi_intracomm)

    # Change fuel temperature and density and export properties
    cell = openmc.lib.cells[1]
    cell.set_temperature(600.0)
    cell.fill.set_density(5.0, 'g/cm3')
    openmc.lib.export_properties(output=False)

    # Import properties to existing model
    model.import_properties("properties.h5")

    # Check to see that values are assigned to the C and python representations
    # First python
    cell = model.geometry.get_all_cells()[1]
    assert cell.temperature == [600.0]
    assert cell.fill.get_mass_density() == pytest.approx(5.0)
    # Now C
    assert openmc.lib.cells[1].get_temperature() == 600.
    assert openmc.lib.materials[1].get_density('g/cm3') == pytest.approx(5.0)

    # Clear the C API
    openmc.lib.finalize()

    # Verify the attributes survived by exporting to XML and re-creating
    model.export_to_xml("with_properties")

    # Load model with properties and confirm temperature/density changed
    model_with_properties = openmc.Model.from_xml(
        'with_properties/geometry.xml',
        'with_properties/materials.xml',
        'with_properties/settings.xml'
    )
    cell = model_with_properties.geometry.get_all_cells()[1]
    assert cell.temperature == [600.0]
    assert cell.fill.get_mass_density() == pytest.approx(5.0)


def test_run(run_in_tmpdir, pin_model_attributes, mpi_intracomm):
    mats, geom, settings, tals, plots, _, _ = pin_model_attributes
    test_model = openmc.Model(geom, mats, settings, tals, plots)

    # This case will run by getting the k-eff and tallies for command-line and
    # C API execution modes and ensuring they give the same result.
    sp_path = test_model.run(output=False)
    with openmc.StatePoint(sp_path) as sp:
        cli_keff = sp.keff
        cli_flux = sp.get_tally(id=1).get_values(scores=['flux'])[0, 0, 0]
        cli_fiss = sp.get_tally(id=1).get_values(scores=['fission'])[0, 0, 0]

    test_model.init_lib(output=False, intracomm=mpi_intracomm)
    sp_path = test_model.run(output=False)
    with openmc.StatePoint(sp_path) as sp:
        lib_keff = sp.keff
        lib_flux = sp.get_tally(id=1).get_values(scores=['flux'])[0, 0, 0]
        lib_fiss = sp.get_tally(id=1).get_values(scores=['fission'])[0, 0, 0]

    # and lets compare results
    assert lib_keff.n == pytest.approx(cli_keff.n, abs=1e-13)
    assert lib_flux == pytest.approx(cli_flux, abs=1e-13)
    assert lib_fiss == pytest.approx(cli_fiss, abs=1e-13)

    # Now we should make sure that the flags for items which should be handled
    # by init are properly set
    with pytest.raises(ValueError):
        test_model.run(threads=1)
    with pytest.raises(ValueError):
        test_model.run(geometry_debug=True)
    with pytest.raises(ValueError):
        test_model.run(restart_file='1.h5')
    with pytest.raises(ValueError):
        test_model.run(tracks=True)

    test_model.finalize_lib()


def test_plots(run_in_tmpdir, pin_model_attributes, mpi_intracomm):
    mats, geom, settings, tals, plots, _, _ = pin_model_attributes
    test_model = openmc.Model(geom, mats, settings, tals, plots)

    # This test cannot check the correctness of the plot, but it can
    # check that a plot was made and that the expected png files are there

    # We will run the test twice, the first time without C API, the second with
    for i in range(2):
        if i == 1:
            test_model.init_lib(output=False, intracomm=mpi_intracomm)
        test_model.plot_geometry(output=False)

        # Now look for the files
        for fname in ('test.png', 'plot_2.png'):
            test_file = Path(fname)
            assert test_file.exists()
            test_file.unlink()

    test_model.finalize_lib()


def test_py_lib_attributes(run_in_tmpdir, pin_model_attributes, mpi_intracomm):
    mats, geom, settings, tals, plots, _, _ = pin_model_attributes
    test_model = openmc.Model(geom, mats, settings, tals, plots)

    test_model.init_lib(output=False, intracomm=mpi_intracomm)

    # Now we can call rotate_cells, translate_cells, update_densities,
    # and update_cell_temperatures and make sure the changes have taken hold.
    # For each we will first try bad inputs to make sure we get the right
    # errors and then we do a good one which calls the material by name and
    # then id to make sure it worked

    # The rotate_cells and translate_cells will work on the cell named fill, as
    # it is filled with a universe and thus the operation will be valid

    # First rotate_cells
    with pytest.raises(TypeError):
        # Make sure it tells us we have a bad names_or_ids type
        test_model.rotate_cells(None, (0, 0, 90))
    with pytest.raises(TypeError):
        test_model.rotate_cells([None], (0, 0, 90))
    with pytest.raises(openmc.exceptions.InvalidIDError):
        # Make sure it tells us we had a bad id
        test_model.rotate_cells([7200], (0, 0, 90))
    with pytest.raises(openmc.exceptions.InvalidIDError):
        # Make sure it tells us we had a bad id
        test_model.rotate_cells(['bad_name'], (0, 0, 90))
    # Now a good one
    assert np.all(openmc.lib.cells[2].rotation == (0., 0., 0.))
    test_model.rotate_cells([2], (0, 0, 90))
    assert np.all(openmc.lib.cells[2].rotation == (0., 0., 90.))

    # And same thing by name
    test_model.rotate_cells(['fuel'], (0, 0, 180))

    # Now translate_cells. We dont need to re-check the TypeErrors/bad ids,
    # because the other functions use the same hidden method as rotate_cells
    assert np.all(openmc.lib.cells[2].translation == (0., 0., 0.))
    test_model.translate_cells([2], (0, 0, 10))
    assert np.all(openmc.lib.cells[2].translation == (0., 0., 10.))

    # Now lets do the density updates.
    # Check initial conditions
    assert openmc.lib.materials[1].get_density('atom/b-cm') == \
        pytest.approx(0.06891296988603757, abs=1e-13)
    mat_a_dens = np.sum(
        list(test_model.materials[0].get_nuclide_atom_densities().values()))
    assert mat_a_dens == pytest.approx(0.06891296988603757, abs=1e-8)
    # Change the density
    test_model.update_densities(['UO2'], 2.)
    assert openmc.lib.materials[1].get_density('atom/b-cm') == \
        pytest.approx(2., abs=1e-13)
    mat_a_dens = np.sum(
        list(test_model.materials[0].get_nuclide_atom_densities().values()))
    assert mat_a_dens == pytest.approx(2., abs=1e-8)

    # Now lets do the cell temperature updates.
    # Check initial conditions
    assert test_model._cells_by_id == \
        {2: geom.root_universe.cells[2], 3: geom.root_universe.cells[3],
         4: geom.root_universe.cells[4],
         1: geom.root_universe.cells[2].fill.cells[1]}
    assert openmc.lib.cells[3].get_temperature() == \
        pytest.approx(293.6, abs=1e-13)
    assert test_model.geometry.root_universe.cells[3].temperature is None
    # Change the temperature
    test_model.update_cell_temperatures([3], 600.)
    assert openmc.lib.cells[3].get_temperature() == \
        pytest.approx(600., abs=1e-13)
    assert test_model.geometry.root_universe.cells[3].temperature == \
        pytest.approx(600., abs=1e-13)

    # And finally material volume
    assert openmc.lib.materials[1].volume == \
        pytest.approx(0.4831931368640985, abs=1e-13)
    # The temperature on the material will be None because its just the default
    assert test_model.materials[0].volume == \
        pytest.approx(0.4831931368640985, abs=1e-13)
    # Change the temperature
    test_model.update_material_volumes(['UO2'], 2.)
    assert openmc.lib.materials[1].volume == pytest.approx(2., abs=1e-13)
    assert test_model.materials[0].volume == pytest.approx(2., abs=1e-13)

    test_model.finalize_lib()


def test_deplete(run_in_tmpdir, pin_model_attributes, mpi_intracomm):
    mats, geom, settings, tals, plots, op_kwargs, chain_file_xml = \
        pin_model_attributes
    with open('test_chain.xml', 'w') as f:
        f.write(chain_file_xml)
    test_model = openmc.Model(geom, mats, settings, tals, plots)

    initial_mat = mats[0].clone()
    initial_u = initial_mat.get_nuclide_atom_densities()['U235']

    # Note that the chain file includes only U-235 fission to a stable Xe136 w/
    # a yield of 100%. Thus all the U235 we lose becomes Xe136

    # In this test we first run without pre-initializing the shared library
    # data and then compare. Then we repeat with the C API already initialized
    # and make sure we get the same answer
    test_model.deplete([1e6], 'predictor', final_step=False,
                       operator_kwargs=op_kwargs,
                       power=1., output=False)
    # Get the new Xe136 and U235 atom densities
    after_xe = mats[0].get_nuclide_atom_densities()['Xe136']
    after_u = mats[0].get_nuclide_atom_densities()['U235']
    assert after_xe + after_u == pytest.approx(initial_u, abs=1e-15)
    assert test_model.is_initialized is False

    # check the tally output
    def check_tally_output():
        with openmc.StatePoint('openmc_simulation_n0.h5') as sp:
            flux = sp.get_tally(id=1).get_values(scores=['flux'])[0, 0, 0]
            fission = sp.get_tally(id=1).get_values(
                scores=['fission'])[0, 0, 0]

            # we're mainly just checking that the result was produced,
            # so a rough numerical comparison doesn't hurt to have.
            assert flux == pytest.approx(13.1, abs=0.2)
            assert fission == pytest.approx(0.47, abs=0.2)

    check_tally_output()

    # Reset the initial material densities
    mats[0].nuclides.clear()
    densities = initial_mat.get_nuclide_atom_densities()
    tot_density = 0.
    for nuc, density in densities.items():
        mats[0].add_nuclide(nuc, density)
        tot_density += density
    mats[0].set_density('atom/b-cm', tot_density)

    # Now we can re-run with the pre-initialized API
    test_model.init_lib(output=False, intracomm=mpi_intracomm)
    test_model.deplete([1e6], 'predictor', final_step=False,
                       operator_kwargs=op_kwargs,
                       power=1., output=False)
    # Get the new Xe136 and U235 atom densities
    after_lib_xe = mats[0].get_nuclide_atom_densities()['Xe136']
    after_lib_u = mats[0].get_nuclide_atom_densities()['U235']
    assert after_lib_xe + after_lib_u == pytest.approx(initial_u, abs=1e-15)
    assert test_model.is_initialized is True

    # And end by comparing to the previous case
    assert after_xe == pytest.approx(after_lib_xe, abs=1e-15)
    assert after_u == pytest.approx(after_lib_u, abs=1e-15)

    check_tally_output()

    test_model.finalize_lib()


def test_calc_volumes(run_in_tmpdir, pin_model_attributes, mpi_intracomm):
    mats, geom, settings, tals, plots, _, _ = pin_model_attributes

    test_model = openmc.Model(geom, mats, settings, tals, plots)

    # With no vol calcs, it should fail
    with pytest.raises(ValueError):
        test_model.calculate_volumes(output=False)

    # Add a cell and mat volume calc
    material_vol_calc = openmc.VolumeCalculation(
        [mats[2]], samples=1000, lower_left=(-.63, -.63, -100.),
        upper_right=(.63, .63, 100.))
    cell_vol_calc = openmc.VolumeCalculation(
        [geom.root_universe.cells[3]], samples=1000,
        lower_left=(-.63, -.63, -100.), upper_right=(.63, .63, 100.))
    test_model.settings.volume_calculations = \
        [material_vol_calc, cell_vol_calc]

    # Now lets compute the volumes and check to see if it was applied
    # First lets do without using the C-API
    # Make sure the volumes are unassigned first
    assert mats[2].volume is None
    assert geom.root_universe.cells[3].volume is None
    test_model.calculate_volumes(output=False, apply_volumes=True)

    # Now let's test that we have volumes assigned; we arent checking the
    # value, just that the value was changed
    assert mats[2].volume > 0.
    assert geom.root_universe.cells[3].volume > 0.

    # Now reset the values
    mats[2].volume = None
    geom.root_universe.cells[3].volume = None

    # And do again with an initialized library
    for file in ['volume_1.h5', 'volume_2.h5']:
        file = Path(file)
        file.unlink()
    test_model.init_lib(output=False, intracomm=mpi_intracomm)
    test_model.calculate_volumes(output=False, apply_volumes=True)
    assert mats[2].volume > 0.
    assert geom.root_universe.cells[3].volume > 0.
    assert openmc.lib.materials[3].volume == mats[2].volume

    test_model.finalize_lib()


def test_model_xml(run_in_tmpdir):

    # load a model from examples
    pwr_model = openmc.examples.pwr_core()

    # export to separate XMLs manually
    pwr_model.settings.export_to_xml('settings_ref.xml')
    pwr_model.materials.export_to_xml('materials_ref.xml')
    pwr_model.geometry.export_to_xml('geometry_ref.xml')

    # now write and read a model.xml file
    pwr_model.export_to_model_xml()
    new_model = openmc.Model.from_model_xml()

    # make sure we can also export this again to separate
    # XML files
    new_model.export_to_xml()


def test_single_xml_exec(run_in_tmpdir):

    pincell_model = openmc.examples.pwr_pin_cell()

    pincell_model.export_to_model_xml('pwr_pincell.xml')

    openmc.run(path_input='pwr_pincell.xml')

    with pytest.raises(RuntimeError, match='ex-em-ell.xml'):
        openmc.run(path_input='ex-em-ell.xml')

    # test that a file in a different directory can be used
    os.mkdir('inputs')
    pincell_model.export_to_model_xml('./inputs/pincell.xml')
    openmc.run(path_input='./inputs/pincell.xml')

    with pytest.raises(RuntimeError, match='input_dir'):
        openmc.run(path_input='input_dir/pincell.xml')
