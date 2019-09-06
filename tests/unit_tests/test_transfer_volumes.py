"""Regression tests for openmc.deplete.Results.transfer_volumes method.


"""

from pytest import approx
import openmc
from openmc.deplete import PredictorIntegrator, ResultsList

from tests import dummy_operator


def test_transfer_volumes(run_in_tmpdir):
    """Unit test of volume transfer in restart calculations."""

    op = dummy_operator.DummyOperator()
    op.output_dir = "test_transfer_volumes"

    # Perform simulation using the predictor algorithm
    dt = [0.75]
    power = 1.0
    PredictorIntegrator(op, dt, power).integrate()

    # Load the files
    res = openmc.deplete.ResultsList.from_hdf5(op.output_dir / "depletion_results.h5")

    # Create a dictionary of volumes to transfer
    res[0].volume['1'] = 1.5
    res[0].volume['2'] = 2.5

    # Create dummy geometry
    mat1 = openmc.Material(material_id=1)
    mat1.depletable = True
    mat2 = openmc.Material(material_id=2)

    cell = openmc.Cell()
    cell.fill = [mat1, mat2]
    root = openmc.Universe()
    root.add_cell(cell)
    geometry = openmc.Geometry(root)

    # Transfer volumes
    res[0].transfer_volumes(geometry)

    assert mat1.volume == 1.5
    assert mat2.volume is None
