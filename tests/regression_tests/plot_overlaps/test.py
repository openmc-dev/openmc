from tests.testing_harness import PlotTestHarness
from tests.regression_tests import config


def test_plot_overlap():
    harness = PlotTestHarness(('plot_1.png', 'plot_2.png', 'plot_3.png',
                               'plot_4.h5'))
    harness.main()
