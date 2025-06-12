from difflib import unified_diff
import glob
import filecmp
import os
from pathlib import Path

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness, colorize

# use a few models from other tests to make sure the same results are
# produced when using a single model.xml file as input
from ..adj_cell_rotation.test import model as adj_cell_rotation_model
from ..lattice_multiple.test import model as lattice_multiple_model
from ..energy_laws.test import model as energy_laws_model
from ..photon_production.test import model as photon_production_model


class ModelXMLTestHarness(PyAPITestHarness):
    """Accept a results file to check against and assume inputs_true is the contents of a model.xml file.
    """
    def __init__(self, model=None, inputs_true=None, results_true=None):
        statepoint_name = f'statepoint.{model.settings.batches}.h5'
        super().__init__(statepoint_name, model, inputs_true)

        self.results_true = 'results_true.dat' if results_true is None else results_true

    def _build_inputs(self):
        self._model.export_to_model_xml()

    def _get_inputs(self):
        return open('model.xml').read()

    def _compare_results(self):
        """Make sure the current results agree with the reference."""
        compare = filecmp.cmp('results_test.dat', self.results_true)
        if not compare:
            expected = open(self.results_true).readlines()
            actual = open('results_test.dat').readlines()
            diff = list(unified_diff(expected, actual, self.results_true,
                                'results_test.dat'))
            if diff:
                print('Result differences:')
                print(''.join(colorize(diff)))
                os.rename('results_test.dat', 'results_error.dat')
            else:
                compare = True
        assert compare, 'Results do not agree'

    def _cleanup(self):
        super()._cleanup()
        if os.path.exists('model.xml'):
            os.remove('model.xml')


test_names = [
    'adj_cell_rotation',
    'lattice_multiple',
    'energy_laws',
    'photon_production'
]


@pytest.mark.parametrize("test_name", test_names, ids=lambda test: test)
def test_model_xml(test_name, request):
    openmc.reset_auto_ids()

    test_path = '../' + test_name
    results = test_path + "/results_true.dat"
    inputs = test_name + "_inputs_true.dat"
    model_name = test_name + "_model"
    harness = ModelXMLTestHarness(request.getfixturevalue(model_name), inputs, results)
    harness.main()

def test_input_arg(run_in_tmpdir):

    pincell = openmc.examples.pwr_pin_cell()

    pincell.settings.particles = 100

    # export to separate XML files and run
    pincell.export_to_xml()
    openmc.run()

    # make sure the executable isn't falling back on the separate XMLs
    for f in glob.glob('*.xml'):
        os.remove(f)
    # now export to a single XML file with a custom name
    pincell.export_to_model_xml('pincell.xml')
    assert Path('pincell.xml').exists()

    # run by specifying that single file
    openmc.run(path_input='pincell.xml')

    # check that this works for plotting too
    openmc.plot_geometry(path_input='pincell.xml')

    # now ensure we get an error for an incorrect filename,
    # even in the presence of other, valid XML files
    pincell.export_to_model_xml()
    with pytest.raises(RuntimeError, match='ex-em-ell.xml'):
        openmc.run(path_input='ex-em-ell.xml')
