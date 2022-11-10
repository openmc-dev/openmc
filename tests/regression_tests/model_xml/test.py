from difflib import unified_diff
import filecmp
import os

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness, colorize

# use a few models from other tests to make sure the same results are
# produced when using a single model.xml file as input
from ..adj_cell_rotation import model as adj_cell_rotation_model
from ..lattice_multiple import model as lattice_mult_model
from ..energy_laws import model as energy_laws_model
from ..photon_production import model as photon_prod_model


class ModelXMLTestHarness(PyAPITestHarness):
    """Accept a results file to check against and assume inputs_true is the contents of a model.xml file.
    """
    def __init__(self, model=None, inputs_true=None, results_true=None):
        statepoint_name = f'statepoint.{model.settings.batches}.h5'
        super().__init__(statepoint_name, model, inputs_true)

        self.results_true = 'results_true.dat' if results_true is None else results_true

    def _build_inputs(self):
        self._model.export_to_xml(separate_xmls=False)

    def _get_inputs(self):
        return open('model.xml').read()

    def _compare_inputs(self):
        """Skip input comparisons for now
        """
        pass

    def _compare_results(self):
        """Make sure the current results agree with the reference."""
        compare = filecmp.cmp('results_test.dat', self.results_true)
        if not compare:
            expected = open(self.results_true).readlines()
            actual = open('results_test.dat').readlines()
            diff = unified_diff(expected, actual, self.results_true,
                                'results_test.dat')
            print('Result differences:')
            print(''.join(colorize(diff)))
            os.rename('results_test.dat', 'results_error.dat')
        assert compare, 'Results do not agree'

    def _cleanup(self):
        super()._cleanup()
        if os.path.exists('model.xml'):
            os.remove('model.xml')


models = [
'adj_cell_rotation_model',
'lattice_mult_model',
'energy_laws_model',
'photon_prod_model'
]
paths = [
'../adj_cell_rotation',
'../lattice_multiple',
'../energy_laws',
'../photon_production'
]


@pytest.mark.parametrize("model, path", zip(models, paths))
def test_model_xml(model, path, request):
    inputs = path + "/inputs_true.dat"
    results = path + "/results_true.dat"
    harness = ModelXMLTestHarness(request.getfixturevalue(model), inputs, results)
    harness.main()