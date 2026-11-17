from pathlib import Path
from dataclasses import replace
import vtk
import numpy as np
from sandbox.tokomak_analysis import Config,build_model
from sandbox.tokomak_analysis.visualization.export_vtk import visualize


def test_native_vtk_roundtrip(tmp_path):
    a=build_model(replace(Config(),vtk_dimension=(12,12,10),enable_mesh=False),tallies=False)
    metadata=visualize(a,tmp_path,plots=False,vtk=True)
    assert metadata['vtk_readback_passed']
    reader=vtk.vtkXMLImageDataReader();reader.SetFileName(str(tmp_path/'visualization/reactor.vti'));reader.Update()
    image=reader.GetOutput();assert image.GetNumberOfCells()==12*12*10
    for name in ('cell_id','material_id','component_id','module_id','sector_id','inboard_outboard'):
        assert image.GetCellData().GetArray(name) is not None
