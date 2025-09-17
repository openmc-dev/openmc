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
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
    settings.source = openmc.IndependentSource(
        space=uniform_dist, constraints={'fissionable': True})

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
    test_model.deplete(timesteps=[1e6], method='predictor', final_step=False,
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
    test_model.deplete(timesteps=[1e6], method='predictor', final_step=False,
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

    # Make sure path can be specified with run
    pincell_model.run(path='my_model.xml')

    os.mkdir('subdir')
    pincell_model.run(path='subdir')


def test_nuclides_to_ignore(run_in_tmpdir, pin_model_attributes):
    """Test nuclides_to_ignore when exporting a model XML"""
    materials, geometry, settings = pin_model_attributes[:3]
    model = openmc.Model(geometry=geometry, settings=settings)

    # grab one of the nuclides present in this model as a test
    test_nuclide = list(materials[0].get_nuclides())[0]

    # exclude the test nuclide from the XML file during export
    model.export_to_model_xml(nuclides_to_ignore=[test_nuclide])

    # ensure that the nuclide doesn't appear after reading in
    # the resulting XML model
    xml_model = openmc.Model.from_model_xml()
    for material in xml_model.materials:
        assert test_nuclide not in material.get_nuclides()


def test_model_plot():
    # plots the geometry with source location and checks the resulting
    # matplotlib includes the correct coordinates for the scatter plot for all
    # basis.

    surface = openmc.Sphere(r=600, boundary_type="vacuum")
    cell = openmc.Cell(region=-surface)
    geometry = openmc.Geometry([cell])
    source = openmc.IndependentSource(space=openmc.stats.Point((1, 2, 3)))
    settings = openmc.Settings(particles=1, batches=1, source=source)
    model = openmc.Model(geometry, settings=settings)

    plot = model.plot(n_samples=1, plane_tolerance=4.0, basis="xy")
    coords = plot.axes.collections[0].get_offsets().data.flatten()
    assert (coords == np.array([1.0, 2.0])).all()

    plot = model.plot(n_samples=1, plane_tolerance=4.0, basis="xz")
    coords = plot.axes.collections[0].get_offsets().data.flatten()
    assert (coords == np.array([1.0, 3.0])).all()

    plot = model.plot(n_samples=1, plane_tolerance=4.0, basis="yz")
    coords = plot.axes.collections[0].get_offsets().data.flatten()
    assert (coords == np.array([2.0, 3.0])).all()

    plot = model.plot(n_samples=1, plane_tolerance=0.1, basis="xy")
    coords = plot.axes.collections[0].get_offsets().data.flatten()
    assert (coords == np.array([])).all()

    # modify model to include another cell that overlaps the original cell entirely
    model.geometry.root_universe.add_cell(openmc.Cell(region=-surface))
    axes = model.plot(show_overlaps=True)
    white = np.array((1.0, 1.0, 1.0))
    red = np.array((1.0, 0.0, 0.0))
    axes_image = axes.get_images()[0]
    image_data = axes_image.get_array()
    # ensure that all of the data in the image data is either white or red
    test_mask = (image_data == white) | (image_data == red)
    assert np.all(test_mask), "Colors other than white or red found in overlap plot image"

    # Close plots to avoid warning
    import matplotlib.pyplot as plt
    plt.close('all')


def test_model_id_map_initialization(run_in_tmpdir):
    model = openmc.examples.pwr_assembly()
    model.init_lib(output=False)

    id_map = model.id_map(
        pixels=(100, 100),
        basis='xy',
        origin=(0, 0, 0),
        width=(10, 10),
    )

    assert id_map.shape == (100, 100, 3)
    assert id_map.dtype == np.int32

    max_cell_id = max(model.geometry.get_all_cells().keys())
    max_material_id = max(model.geometry.get_all_materials().keys())

    # add some spot checks for the id_map
    # Check that the array contains valid cell/material IDs (not all -2)
    # The -2 values indicate outside the geometry
    assert not np.all(id_map == -2), "All values are -2, indicating no valid geometry found"

    # Check that we have valid cell IDs (first dimension)
    valid_cell_ids = id_map[:, :, 0]
    assert np.any(valid_cell_ids >= 0), "No valid cell IDs found in the id_map"

    # Check that we have valid material IDs (third dimension)
    valid_material_ids = id_map[:, :, 2]
    assert np.any(valid_material_ids >= 0), "No valid material IDs found in the id_map"

    # Check that the middle dimension (cell instances) is consistent
    # Cell instances should be >= 0 when cell IDs are valid
    cell_instances = id_map[:, :, 1]
    valid_cells = valid_cell_ids >= 0
    if np.any(valid_cells):
        assert np.all(cell_instances[valid_cells] >= 0), "Invalid cell instances found for valid cells"

    # Check that the array contains reasonable ranges of values
    # Cell IDs should be within the expected range for the assembly
    if np.any(valid_cell_ids >= 0):
        max_map_cell_id = np.max(valid_cell_ids)
        assert max_map_cell_id <= max_cell_id, \
            f"Cell ID {max_map_cell_id} in the map is greater than the maximum cell ID {max_cell_id}"

    # Material IDs should be within the expected range
    if np.any(valid_material_ids >= 0):
        max_map_material_id = np.max(valid_material_ids)
        assert max_map_material_id <= max_material_id, \
            f"Material ID {max_map_material_id} in the map is greater than the maximum material ID {max_material_id}"

    # Test id_map with pixels outside the model geometry
    # Use a plot that's far from the model center to ensure we get -2 values
    outside_id_map = model.id_map(
        pixels=(50, 50),
        basis='xy',
        origin=(1000, 1000, 0),  # Far from the model center
        width=(10, 10),
    )

    assert outside_id_map.shape == (50, 50, 3)
    assert outside_id_map.dtype == np.int32

    # All values should be -2 (outside geometry) for this plot
    assert np.all(outside_id_map == -2), "Expected all values to be -2 for plot outside model geometry"

    # Verify that the outside plot has the correct structure
    assert np.all(outside_id_map[:, :, 0] == -2), "Cell IDs should all be -2 outside geometry"
    assert np.all(outside_id_map[:, :, 1] == -2), "Cell instances should all be -2 outside geometry"
    assert np.all(outside_id_map[:, :, 2] == -2), "Material IDs should all be -2 outside geometry"

    # if the model is already initialized, it should not be finalized
    # after calling this method
    model.id_map(
        pixels=(100, 100),
        basis='xy',
        origin=(0, 0, 0),
        width=(10, 10),
    )
    assert model.is_initialized

    # if the model is not initialized, it should be finalized
    # before exiting this method
    model.finalize_lib()
    model.id_map(
        pixels=(100, 100),
        basis='xy',
        origin=(0, 0, 0),
        width=(10, 10),
    )
    assert not model.is_initialized


def test_id_map_aligned_model():
    """Test id_map with a 2x2 lattice where pixel boundaries align to cell boundaries"""
    # Create materials -- identical compositions, different IDs
    mat1 = openmc.Material(material_id=1, name='Material 1')
    mat1.set_density('g/cm3', 1.0)
    mat1.add_element('H', 1.0)

    mat2 = openmc.Material(material_id=2, name='Material 2')
    mat2.set_density('g/cm3', 1.0)
    mat2.add_element('H', 1.0)

    mat3 = openmc.Material(material_id=3, name='Material 3')
    mat3.set_density('g/cm3', 1.0)
    mat3.add_element('H', 1.0)

    mat4 = openmc.Material(material_id=4, name='Material 4')
    mat4.set_density('g/cm3', 1.0)
    mat4.add_element('H', 1.0)

    outer_mat = openmc.Material(material_id=5, name='Material 5')
    outer_mat.set_density('g/cm3', 1.0)
    outer_mat.add_element('H', 1.0)

    inner_materials = [mat1, mat2, mat3, mat4]

    # Create square surface that fits inside the lattice cell
    # Lattice cell is 1 cm x 1 cm, so square will be 0.6 cm x 0.6 cm centered on the origin
    square = openmc.model.RectangularPrism(0.6, 0.6, boundary_type='transmission')

    # Create cells for this universe
    inner_cell = openmc.Cell(cell_id=10, region=-square, name='inner_cell')
    inner_cell.fill = inner_materials

    outer_cell = openmc.Cell(cell_id=20, region=+square, name='outer_cell')
    outer_cell.fill = outer_mat

    # Create universe
    universe = openmc.Universe(universe_id=100, cells=[inner_cell, outer_cell])

    # Create 2x2 lattice
    lattice = openmc.RectLattice(lattice_id=1)
    lattice.lower_left = [-1.0, -1.0]
    lattice.pitch = [1.0, 1.0]
    lattice.universes = [[universe, universe], [universe, universe]]

    # Create outer boundary
    outer_boundary = openmc.model.RectangularPrism(2.0, 2.0, boundary_type='vacuum')

    # Create root cell
    root_cell = openmc.Cell(cell_id=1, name='root', fill=lattice, region=-outer_boundary)

    # Create geometry
    geometry = openmc.Geometry([root_cell])

    # Create settings
    settings = openmc.Settings()
    settings.particles = 1000
    settings.batches = 10

    # Create model
    model = openmc.Model(settings=settings, geometry=geometry)

    # Generate id_map with pixel boundaries aligned to cell boundaries
    # The model is 2 cm x 2 cm, so we'll use 200x200 pixels to get 0.01 cm resolution
    # This allows us to align pixels with the squares inside each lattice cell
    id_map = model.id_map(
        pixels=(200, 200),
        basis='xy',
        origin=(0.0, 0.0, 0.0),  # Align with lattice lower_left
        width=(2.0, 2.0),     # Align with lattice size
    )

    # Verify id_map properties
    assert id_map.shape == (200, 200, 3)
    assert id_map.dtype == np.int32

    cell_id_map = id_map[:, :, 0]
    material_ids_map = id_map[:, :, 2]

    # Check that we have valid cell IDs (not all -2)
    assert np.any(cell_id_map >= 0), "No valid cell IDs found in the id_map"

    # Check that we have valid material IDs
    assert np.any(material_ids_map >= 0), "No valid material IDs found in the id_map"

    # Check that the expected cell IDs are present
    expected_cell_ids = [10, 20]  # Root cell, inner cell, outer cell
    found_cell_ids = np.unique(cell_id_map[cell_id_map >= 0])
    for cell_id in expected_cell_ids:
        assert cell_id in found_cell_ids, f"Expected cell ID {cell_id} not found in id_map"

    # Check that the expected material IDs are present
    expected_material_ids = [1, 2, 3, 4, 5]  # All materials defined above
    found_material_ids = np.unique(material_ids_map[material_ids_map >= 0])
    for mat_id in expected_material_ids:
        assert mat_id in found_material_ids, f"Expected material ID {mat_id} not found in id_map"

    # Test specific regions to verify lattice structure
    # Check center of each lattice cell (should be inner cells)
    # Lattice cell centers are at (-0.5, -0.5), (0.5, -0.5), (-0.5, 0.5), (0.5, 0.5)
    # With 200x200 pixels over 2x2 units, each pixel is 0.01 units

    # Bottom-left lattice cell center (should be inner cell 10)
    bl_cell, bl_instance, bl_material = id_map[-50, 50]
    assert bl_cell == 10, f"Expected cell ID 10 at bottom-left center, got {bl_cell}"
    assert bl_instance == 0, f"Expected cell instance 0 at bottom-left center, got {bl_instance}"
    assert bl_material == 1, f"Expected material ID 1 at bottom-left center, got {bl_material}"

    # Bottom-right lattice cell center (should be inner cell 10)
    br_cell, br_instance, br_material = id_map[-50, 150]
    assert br_cell == 10, f"Expected cell ID 10 at bottom-right center, got {br_cell}"
    assert br_instance == 1, f"Expected cell instance 1 at bottom-right center, got {br_instance}"
    assert br_material == 2, f"Expected material ID 2 at bottom-right center, got {br_material}"

    # Top-left lattice cell center (should be inner cell 10)
    tl_cell, tl_instance, tl_material = id_map[-150, 50]
    assert tl_cell == 10, f"Expected cell ID 10 at top-left center, got {tl_cell}"
    assert tl_instance == 2, f"Expected cell instance 2 at top-left center, got {tl_instance}"
    assert tl_material == 3, f"Expected material ID 3 at top-left center, got {tl_material}"

    # Top-right lattice cell center (should be inner cell 10)
    tr_cell, tr_instance, tr_material = id_map[-150, 150]
    assert tr_cell == 10, f"Expected cell ID 10 at top-right center, got {tr_cell}"
    assert tr_instance == 3, f"Expected cell instance 3 at top-right center, got {tr_instance}"
    assert tr_material == 4, f"Expected material ID 4 at top-right center, got {tr_material}"

    # Check that the model is properly finalized after id_map call
    assert not model.is_initialized, "Model should be finalized after id_map call"

    # Check that the values at the corners are correctly set as the outer cell and material
    bl_cell, bl_instance, bl_material = id_map[-1, 0]
    assert bl_cell == 20, f"Expected cell ID 20 at bottom-left corner, got {bl_cell}"
    assert bl_instance == 0, f"Expected cell instance 0 at bottom-left corner, got {bl_instance}"
    assert bl_material == 5, f"Expected material ID 5 at bottom-left corner, got {bl_material}"

    br_cell, br_instance, br_material = id_map[-1, -1]
    assert br_cell == 20, f"Expected cell ID 20 at bottom-right corner, got {br_cell}"
    assert br_instance == 1, f"Expected cell instance 1 at bottom-right corner, got {br_instance}"
    assert br_material == 5, f"Expected material ID 5 at bottom-right corner, got {br_material}"

    tl_cell, tl_instance, tl_material = id_map[0, 0]
    assert tl_cell == 20, f"Expected cell ID 20 at top-left corner, got {tl_cell}"
    assert tl_instance == 2, f"Expected cell instance 2 at top-left corner, got {tl_instance}"
    assert tl_material == 5, f"Expected material ID 5 at top-left corner, got {tl_material}"

    tr_cell, tr_instance, tr_material = id_map[0, -1]
    assert tr_cell == 20, f"Expected cell ID 20 at top-right corner, got {tr_cell}"
    assert tr_instance == 3, f"Expected cell instance 3 at top-right corner, got {tr_instance}"
    assert tr_material == 5, f"Expected material ID 5 at top-right corner, got {tr_material}"

def test_setter_from_list():
    mat = openmc.Material()
    model = openmc.Model(materials=[mat])
    assert isinstance(model.materials, openmc.Materials)

    tally = openmc.Tally()
    model = openmc.Model(tallies=[tally])
    assert isinstance(model.tallies, openmc.Tallies)

    plot = openmc.Plot()
    model = openmc.Model(plots=[plot])
    assert isinstance(model.plots, openmc.Plots)
