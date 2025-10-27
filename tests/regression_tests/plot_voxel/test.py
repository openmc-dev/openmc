import pytest

from tests.testing_harness import PlotTestHarness
from tests.regression_tests import config


vtk = pytest.importorskip('vtk')

def test_plot_voxel():
    harness = PlotTestHarness(('plot_4.h5', 'plot.vti'), voxel_convert_checks=['plot_4.h5'])
    harness.main()
