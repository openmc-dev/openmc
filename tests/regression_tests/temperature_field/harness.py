from openmc.utility_funcs import change_directory

from tests.testing_harness import PyAPITestHarness


class TemperatureFieldTestHarness(PyAPITestHarness):
    """Verify equivalence to cell temperatures before accepting a regression."""

    def __init__(self, model, reference_model):
        super().__init__('statepoint.20.h5', model)
        self.reference_model = reference_model

    def _run_openmc(self):
        super()._run_openmc()
        actual = self._get_results()
        reference = PyAPITestHarness(self._sp_name, self.reference_model)
        with change_directory(tmpdir=True):
            reference._build_inputs()
            reference._run_openmc()
            expected = reference._get_results()
        assert actual == expected, (
            'Temperature field differs from explicit cell temperatures:\n'
            f'Field:\n{actual}\nReference:\n{expected}')

    def _overwrite_inputs(self):
        # PyAPITestHarness normally writes inputs before running the model.
        # Defer both baseline writes until the equivalence assertion succeeds.
        pass

    def _overwrite_results(self):
        super()._overwrite_inputs()
        super()._overwrite_results()
