import pytest
import openmc
import openmc.lib


def test_import_properties(run_in_tmpdir, mpi_intracomm):
    """Test importing properties on the Model class """

    # Create PWR pin cell model and write XML files
    openmc.reset_auto_ids()
    model = openmc.examples.pwr_pin_cell()
    model.export_to_xml()

    # Change fuel temperature and density and export properties
    openmc.lib.init(intracomm=mpi_intracomm)
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
