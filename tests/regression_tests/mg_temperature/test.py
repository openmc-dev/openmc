import os
from tests.regression_tests.mg_temperature.build_2g import *
from tests.testing_harness import *


class MgTemperatureTestHarness(TestHarness):

    def execute_test(self):
        """Run OpenMC with the appropriate arguments and check the outputs."""
        base_dir = os.getcwd()
        print("Base dir is {}".format(base_dir))
        macro_xs = create_openmc_2mg_libs(names)
        dirs = ('micro/nearest/case1', 'micro/nearest/case2',
                'micro/nearest/case3', 'micro/interpolation/case1',
                'micro/interpolation/case2',
                'macro/nearest/case1', 'macro/nearest/case2',
                'macro/nearest/case3', 'macro/interpolation/case1',
                'macro/interpolation/case2')
        temperatures = (300., 600., 900.,
                        520., 600.,
                        300., 600., 900.,
                        520., 600)
        methods = 2 * (3 * ('nearest',) + 2 * ('interpolation',))
        analyt_interp = 10 * [None]
        analyt_interp[3] = (600. - 520.) / 300.
        analyt_interp[8] = (600. - 520.) / 300.
        try:
            for d, t, m, ai in zip(dirs, temperatures, methods, analyt_interp):
                os.chdir(os.path.join(base_dir, d))
                if (d[:5] == 'macro'):
                    build_inf_model(['macro'], '../../../macro_2g.h5', t, m)
                else:
                    build_inf_model(names, '../../../micro_2g.h5', t, m)
                if not ai:
                    kanalyt = analytical_solution_2g_therm(macro_xs[t])
                else:
                    kanalyt = analytical_solution_2g_therm(macro_xs[300],
                                                           macro_xs[600], ai)
                self._run_openmc()
                self._test_output_created()
                results = self._get_results()
                results += "k-analytical:\n"
                results += "{:12.6E}".format(kanalyt)
                self._write_results(results)
                self._compare_results()
        finally:
            for d in dirs:
                os.chdir(os.path.join(base_dir, d))
                self._cleanup()
            os.chdir(base_dir)
            for f in ['micro_2g.h5', 'macro_2g.h5']:
                if os.path.exists(f):
                    os.remove(f)

    def update_results(self):
        """Update the results_true using the current version of OpenMC."""
        base_dir = os.getcwd()
        print("Base dir is {}".format(base_dir))
        macro_xs = create_openmc_2mg_libs(names)
        dirs = ('micro/nearest/case1', 'micro/nearest/case2',
                'micro/nearest/case3', 'micro/interpolation/case1',
                'micro/interpolation/case2',
                'macro/nearest/case1', 'macro/nearest/case2',
                'macro/nearest/case3', 'macro/interpolation/case1',
                'macro/interpolation/case2')
        temperatures = (300., 600., 900.,
                        520., 600.,
                        300., 600., 900.,
                        520., 600)
        methods = 2 * (3 * ('nearest',) + 2 * ('interpolation',))
        analyt_interp = 10 * [None]
        analyt_interp[3] = (600. - 520.) / 300.
        analyt_interp[8] = (600. - 520.) / 300.
        try:
            for d, t, m, ai in zip(dirs, temperatures, methods, analyt_interp):
                os.chdir(os.path.join(base_dir, d))
                if (d[:5] == 'macro'):
                    build_inf_model(['macro'], '../../../macro_2g.h5', t)
                else:
                    build_inf_model(names, '../../../micro_2g.h5', t)
                if not ai:
                    kanalyt = analytical_solution_2g_therm(macro_xs[t])
                else:
                    kanalyt = analytical_solution_2g_therm(macro_xs[300],
                                                           macro_xs[600], ai)
                self._run_openmc()
                self._test_output_created()
                results = self._get_results()
                results += "k-analytical:\n"
                results += "{:12.6E}".format(kanalyt)
                self._write_results(results)
                self._overwrite_results(results)
                self._compare_results()
        finally:
            for d in dirs:
                os.chdir(os.path.join(base_dir, d))
                self._cleanup()
            os.chdir(base_dir)
            for f in ['micro_2g.h5', 'macro_2g.h5']:
                if os.path.exists(f):
                    os.remove(f)


def test_mg_temperature():
    harness = MgTemperatureTestHarness('statepoint.15.h5')
    harness.main()
