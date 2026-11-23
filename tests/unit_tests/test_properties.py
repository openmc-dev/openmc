import openmc.lib
import pytest


def test_properties_volume_mode(sphere_model, run_in_tmpdir, mpi_intracomm, monkeypatch):
    monkeypatch.setattr(
        sphere_model.settings, 'temperature',
        {'method': 'interpolation', 'default': 600.0},
    )
    sphere_model.export_to_model_xml()
    (cell_id,) = sphere_model.geometry.get_all_cells()

    with openmc.lib.run_in_memory(args=['-c'], intracomm=mpi_intracomm):
        openmc.lib.export_properties('properties.h5')
        openmc.lib.import_properties('properties.h5')
        assert openmc.lib.cells[cell_id].get_temperature() == pytest.approx(600.0)
