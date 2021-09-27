from typing import Type
import pytest
from pathlib import Path
import openmc
import openmc.lib


def test_init(run_in_tmpdir, pin_model_attributes):
    mats, geom, settings, tals, plots, chain, fission_q, _  = \
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
    assert test_model._cells_by_id == {1: geom.root_universe.cells[1],
                                       2: geom.root_universe.cells[2],
                                       3: geom.root_universe.cells[3]}
    # No cell name for 2 and 3, so we expect a blank name to be assigned to
    # cell 3 due to overwriting
    assert test_model._cells_by_name == {
        'fuel': geom.root_universe.cells[1], '': geom.root_universe.cells[3]}
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
    assert [c.fill.name for c in test_model.geometry.root_universe.cells.values()] == \
        [c.fill.name for c in geom.root_universe.cells.values()]
    assert [mat.name for mat in test_model.materials] == [mat.name for mat in mats]
    # We will assume the attributes of settings that are custom objects are
    # OK if the others are so we dotn need to implement explicit comparisons
    no_test = ['_source', '_entropy_mesh']
    assert sorted(k for k in test_model.settings.__dict__.keys() if k not in no_test) == \
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
    assert test_model._cells_by_id == {\
        1: test_model.geometry.root_universe.cells[1],
        2: test_model.geometry.root_universe.cells[2],
        3: test_model.geometry.root_universe.cells[3]}
    # No cell name for 2 and 3, so we expect a blank name to be assigned to
    # cell 3 due to overwriting
    assert test_model._cells_by_name == {
        'fuel': test_model.geometry.root_universe.cells[1],
        '': test_model.geometry.root_universe.cells[3]}
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
    openmc.lib.finalize()

    # Import properties to existing model and re-export to new directory
    model.import_properties("properties.h5")
    model.export_to_xml("with_properties")

    # Load model with properties and confirm temperature/density has been changed
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
    # Set the depletion module to use the test intracomm as the depletion mods
    # intracomm is the one that init uses. We will not perturb system state
    # outside of this test by re-setting the openmc.deplete.comm to what it was
    # before when we exist

    # This case will run by getting the k-eff and tallies for command-line and
    # C-API execution modes and ensuring they give the same result.
    sp_path = test_model.run(output=False)
    with openmc.StatePoint(sp_path) as sp:
        cli_keff = sp.k_combined
        cli_flux = sp.get_tally(id=1).get_values()[0, 0, 0]

    if mpi_intracomm is not None:
        orig_comm = openmc.deplete.comm
        openmc.deplete.comm = mpi_intracomm
    test_model.settings.verbosity = 1
    test_model.init_C_api()
    sp_path = test_model.run()
    with openmc.StatePoint(sp_path) as sp:
        C_keff = sp.k_combined
        C_flux = sp.get_tally(id=1).get_values()[0, 0, 0]
    test_model.clear_C_api()

    # and lets compare results
    assert (C_keff - cli_keff) < 1e-15
    assert (C_flux - cli_flux) < 1e-15

    # And before done, reset the deplete communicator
    if mpi_intracomm is not None:
        openmc.deplete.comm = orig_comm


