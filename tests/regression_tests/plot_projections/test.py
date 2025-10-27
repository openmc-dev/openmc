from tests.testing_harness import PlotTestHarness
from tests.regression_tests import config

def test_plot():
    harness = PlotTestHarness(('plot_1.png', 'example1.png', 'example2.png',
                               'example3.png', 'orthographic_example1.png',
                               'phong.png', 'phong_diffuse.png',
                               'phong_move_light.png'))
    harness.main()
