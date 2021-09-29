import numpy as np
import pytest
from pathlib import Path
from shutil import which
import openmc
import openmc.lib


def test_init(run_in_tmpdir, pin_model_attributes):
    mats, geom, settings, tals, plots, chain, fission_q, _ = \
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
    assert test_model.chain_file is None
    assert test_model.fission_q is None
    assert test_model._materials_by_id == {}
    assert test_model._materials_by_name == {}
    assert test_model._cells_by_id == {}
    assert test_model._cells_by_name == {}
    assert test_model.C_init is False

    # Now check proper init of an actual model. Assume no interference between
    # parameters and so we can apply them all at once instead of testing one
    # parameter initialization at a time
    test_model = openmc.Model(geom, mats, settings, tals, plots, chain,
                              fission_q)
    assert test_model.geometry is geom
    assert test_model.materials is mats
    assert test_model.settings is settings
    assert test_model.tallies is tals
    assert test_model.plots is plots
    assert test_model.chain_file == Path(chain).resolve()
    assert test_model.fission_q is fission_q
    assert test_model._materials_by_id == {1: mats[0], 2: mats[1], 3: mats[2]}
    assert test_model._materials_by_name == {'UO2': mats[0], 'Zirc': mats[1],
                                             'Borated water': mats[2]}
    # The last cell is the one that contains the infinite fuel
    assert test_model._cells_by_id == \
        {2: geom.root_universe.cells[2], 3: geom.root_universe.cells[3],
         4: geom.root_universe.cells[4],
         1: geom.root_universe.cells[2].fill.cells[1]}
    # No cell name for 2 and 3, so we expect a blank name to be assigned to
    # cell 3 due to overwriting
    assert test_model._cells_by_name == {
        'fuel': geom.root_universe.cells[2], '': geom.root_universe.cells[4],
        'inf fuel': geom.root_universe.cells[2].fill.cells[1]}
    assert test_model.C_init is False

    # Finally test the parameter type checking by passing bad types and
    # obtaining the right exception types
    def_params = [geom, mats, settings, tals, plots, chain, fission_q]
    for i in range(len(def_params)):
        args = def_params.copy()
        # Try an integer, as that is a bad type for all arguments
        args[i] = i
        with pytest.raises(TypeError):
            test_model = openmc.Model(*args)


def test_from_xml(run_in_tmpdir, pin_model_attributes):
    mats, geom, settings, tals, plots, _, _, _ = pin_model_attributes

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
    assert len(test_model.tallies) == 0
    assert len(test_model.plots) == 0
    assert test_model.chain_file is None
    assert test_model.fission_q is None
    assert test_model._materials_by_id == \
        {1: test_model.materials[0], 2: test_model.materials[1],
         3: test_model.materials[2]}
    assert test_model._materials_by_name == {
        'UO2': test_model.materials[0], 'Zirc': test_model.materials[1],
        'Borated water': test_model.materials[2]}
    assert test_model._cells_by_id == {
        2: test_model.geometry.root_universe.cells[2],
        3: test_model.geometry.root_universe.cells[3],
        4: test_model.geometry.root_universe.cells[4],
        1: test_model.geometry.root_universe.cells[2].fill.cells[1]}
    # No cell name for 2 and 3, so we expect a blank name to be assigned to
    # cell 3 due to overwriting
    assert test_model._cells_by_name == {
        'fuel': test_model.geometry.root_universe.cells[2],
        '': test_model.geometry.root_universe.cells[4],
        'inf fuel': test_model.geometry.root_universe.cells[2].fill.cells[1]}
    assert test_model.C_init is False


def test_init_clear_C_api(run_in_tmpdir, pin_model_attributes, mpi_intracomm):
    # We are going to init and then make sure data is loaded
    mats, geom, settings, tals, plots, _, _, _ = pin_model_attributes
    test_model = openmc.Model(geom, mats, settings, tals, plots)
    # Set the depletion module to use the test intracomm as the depletion mods
    # intracomm is the one that init uses. We will not perturb system state
    # outside of this test by re-setting the openmc.deplete.comm to what it was
    # before when we exist
    if mpi_intracomm is not None:
        orig_comm = openmc.deplete.comm
        openmc.deplete.comm = mpi_intracomm
    test_model.init_C_api()

    # First check that the API is advertised as initialized
    assert openmc.lib.LIB_INIT is True
    assert test_model.C_init is True
    # Now make sure it actually is initialized by making a call to the lib
    c_mat = openmc.lib.find_material((0.6, 0., 0.))
    # This should be Borated water
    assert c_mat.name == 'Borated water'
    assert c_mat.id == 3

    # Ok, now lets test that we can clear the data and check that it is cleared
    test_model.clear_C_api()

    # First check that the API is advertised as initialized
    assert openmc.lib.LIB_INIT is False
    assert test_model.C_init is False
    # Note we cant actually test that a sys call fails because we should get a
    # seg fault

    # And before done, reset the deplete communicator
    if mpi_intracomm is not None:
        openmc.deplete.comm = orig_comm

    # TODO: in above test, tests are necessary for all the arguments of init


def test_import_properties(run_in_tmpdir, mpi_intracomm):
    """Test importing properties on the Model class """

    # Create PWR pin cell model and write XML files
    openmc.reset_auto_ids()
    model = openmc.examples.pwr_pin_cell()
    model.init_C_api()

    # Change fuel temperature and density and export properties
    cell = openmc.lib.cells[1]
    cell.set_temperature(600.0)
    cell.fill.set_density(5.0, 'g/cm3')
    openmc.lib.export_properties()

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

    # Clear the C-API
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
    mats, geom, settings, tals, plots, _, _, _ = pin_model_attributes
    test_model = openmc.Model(geom, mats, settings, tals, plots)

    # This case will run by getting the k-eff and tallies for command-line and
    # C-API execution modes and ensuring they give the same result.
    sp_path = test_model.run(output=False)
    with openmc.StatePoint(sp_path) as sp:
        cli_keff = sp.k_combined
        cli_flux = sp.get_tally(id=1).get_values()[0, 0, 0]

    if mpi_intracomm is not None:
        # Set the depletion module to use the test intracomm as the depletion
        # module's intracomm is the one that init uses. We will not perturb
        # system state outside of this test by re-setting the
        # openmc.deplete.comm to what it was before when we exist
        orig_comm = openmc.deplete.comm
        openmc.deplete.comm = mpi_intracomm
    test_model.init_C_api()
    sp_path = test_model.run(output=False)
    with openmc.StatePoint(sp_path) as sp:
        C_keff = sp.k_combined
        C_flux = sp.get_tally(id=1).get_values()[0, 0, 0]

    # and lets compare results
    assert abs(C_keff - cli_keff) < 1e-13
    assert abs(C_flux - cli_flux) < 1e-13

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

    test_model.clear_C_api()

    # And before done, reset the deplete communicator
    if mpi_intracomm is not None:
        openmc.deplete.comm = orig_comm


def test_plots(run_in_tmpdir, pin_model_attributes, mpi_intracomm):
    mats, geom, settings, tals, plots, _, _, _ = pin_model_attributes
    test_model = openmc.Model(geom, mats, settings, tals, plots)

    if mpi_intracomm is not None:
        # Set the depletion module to use the test intracomm as the depletion
        # module's intracomm is the one that init uses. We will not perturb
        # system state outside of this test by re-setting the
        # openmc.deplete.comm to what it was before when we exist
        orig_comm = openmc.deplete.comm
        openmc.deplete.comm = mpi_intracomm

    # This test cannot check the correctness of the plot, but it can
    # check that a plot was made and that the expected ppm and png files are
    # there

    # We will only test convert if it is on the system, so as not to add an
    # extra dependency just for tests
    convert = which('convert') is not None
    if convert:
        exts = ['ppm', 'png']
    else:
        exts = ['ppm']

    # We will run the test twice, the first time without C-API, the second with
    for i in range(2):
        if i == 1:
            test_model.init_C_api()
        test_model.plot_geometry(output=True, convert=convert)

        # Now look for the files, expect to find test.ppm, plot_2.ppm, and if
        # convert is True, test.png, plot_2.png
        for fname in ['test.', 'plot_2.']:
            for ext in exts:
                test_file = Path('./{}{}'.format(fname, ext))
                assert test_file.exists()
                test_file.unlink()

    test_model.clear_C_api()

    # And before done, reset the deplete communicator
    if mpi_intracomm is not None:
        openmc.deplete.comm = orig_comm


def test_py_C_attributes(run_in_tmpdir, pin_model_attributes, mpi_intracomm):
    mats, geom, settings, tals, plots, _, _, _ = pin_model_attributes
    test_model = openmc.Model(geom, mats, settings, tals, plots)

    if mpi_intracomm is not None:
        # Set the depletion module to use the test intracomm as the depletion
        # module's intracomm is the one that init uses. We will not perturb
        # system state outside of this test by re-setting the
        # openmc.deplete.comm to what it was before when we exist
        orig_comm = openmc.deplete.comm
        openmc.deplete.comm = mpi_intracomm
    test_model.init_C_api()

    # Now we can call rotate_cells, translate_cells, update_densities,
    # update_cell_temperatures, and update_material_temperatures and make sure
    # the changes have taken hold.
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
    assert abs(openmc.lib.materials[1].get_density(
        'atom/b-cm') - 0.06891296988603757) < 1e-13
    mat_a_dens = np.sum(
        [v[1] for v in test_model.materials[0].
            get_nuclide_atom_densities().values()])
    assert abs(mat_a_dens - 0.06891296988603757) < 1e-8
    # Change the density
    test_model.update_densities(['UO2'], 2.)
    assert abs(openmc.lib.materials[1].get_density('atom/b-cm') - 2.) < 1e-13
    mat_a_dens = np.sum(
        [v[1] for v in test_model.materials[0].
            get_nuclide_atom_densities().values()])
    assert abs(mat_a_dens - 2.) < 1e-8

    # Now lets do the cell temperature updates.
    # Check initial conditions
    assert test_model._cells_by_id == \
        {2: geom.root_universe.cells[2], 3: geom.root_universe.cells[3],
         4: geom.root_universe.cells[4],
         1: geom.root_universe.cells[2].fill.cells[1]}
    assert abs(openmc.lib.cells[3].get_temperature() - 293.6) < 1e-13
    assert test_model.geometry.root_universe.cells[3].temperature is None
    # Change the temperature
    test_model.update_cell_temperatures([3], 600.)
    assert abs(openmc.lib.cells[3].get_temperature() - 600.) < 1e-13
    assert abs(test_model.geometry.root_universe.cells[3].temperature -
               600.) < 1e-13

    # Now lets do the material temperature updates.
    # Check initial conditions
    with pytest.raises(NotImplementedError):
        test_model.update_material_temperatures(['UO2'], 600.)
    # TODO: When C-API material.set_temperature is implemented, uncomment below
    # assert abs(openmc.lib.materials[1].temperature - 293.6) < 1e-13
    # # The temperature on the material will be None because its just the
    # # default assert test_model.materials[0].temperature is None
    # # Change the temperature
    # test_model.update_material_temperatures(['UO2'], 600.)
    # assert abs(openmc.lib.materials[1].temperature - 600.) < 1e-13
    # assert abs(test_model.materials[0].temperature - 600.) < 1e-13

    # And finally material volume
    # import pdb; pdb.set_trace()
    assert abs(openmc.lib.materials[1].volume - 0.4831931368640985) < 1e-13
    # The temperature on the material will be None because its just the default
    assert abs(test_model.materials[0].volume - 0.4831931368640985) < 1e-13
    # Change the temperature
    test_model.update_material_volumes(['UO2'], 2.)
    assert abs(openmc.lib.materials[1].volume - 2.) < 1e-13
    assert abs(test_model.materials[0].volume - 2.) < 1e-13

    test_model.clear_C_api()

    # And before done, reset the deplete communicator
    if mpi_intracomm is not None:
        openmc.deplete.comm = orig_comm
