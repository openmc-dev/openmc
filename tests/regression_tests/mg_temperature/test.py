import os
from tests.regression_tests.mg_temperature.build_2g import *
from tests.testing_harness import *
import shutil


class MgTemperatureTestHarness(TestHarness):

    def execute_test(self):
        """Run OpenMC with the appropriate arguments and check the outputs."""
        base_dir = os.getcwd()
        overall_results = []
        macro_xs = create_openmc_2mg_libs(names)
        types = ('micro', 'micro',
                 'micro', 'micro',
                 'micro',
                 'macro', 'macro',
                 'macro', 'macro',
                 'macro')
        temperatures = (300., 600., 900.,
                        520., 600.,
                        300., 600., 900.,
                        520., 600)
        methods = 2 * (3 * ('nearest',) + 2 * ('interpolation',))
        analyt_interp = 10 * [None]
        analyt_interp[3] = (600. - 520.) / 300.
        analyt_interp[8] = (600. - 520.) / 300.
        try:
            if (os.path.isdir("./temp")):
                shutil.rmtree("./temp")
            os.mkdir("temp")
            os.chdir(os.path.join(base_dir, "temp"))
            for cs, t, m, ai in zip(types, temperatures, methods, analyt_interp):
                if (cs == 'macro'):
                    build_inf_model(['macro'], '../macro_2g.h5', t, m)
                else:
                    build_inf_model(names, '../micro_2g.h5', t, m)
                if not ai:
                    kanalyt = analytical_solution_2g_therm(macro_xs[t])
                else:
                    kanalyt = analytical_solution_2g_therm(macro_xs[300],
                                                           macro_xs[600], ai)
                self._run_openmc()
                self._test_output_created()
                string = "{}, method: {}, t: {}, {}kanalyt\n{:12.6E}\n"
                results = string.format(cs, m, t, self._get_results(), kanalyt)
                overall_results.append(results)
            os.chdir(base_dir)
            self._write_results("".join(overall_results))
            self._compare_results()
        finally:
            os.chdir(base_dir)
            shutil.copyfile("results_test.dat", "results_true.dat")
            if (os.path.isdir("./temp")):
                shutil.rmtree("./temp")
            self._cleanup()
            for f in ['micro_2g.h5', 'macro_2g.h5']:
                if os.path.exists(f):
                    os.remove(f)

    def update_results(self):
        """Update the results_true using the current version of OpenMC."""
        base_dir = os.getcwd()
        overall_results = []
        macro_xs = create_openmc_2mg_libs(names)
        types = ('micro', 'micro',
                 'micro', 'micro',
                 'micro',
                 'macro', 'macro',
                 'macro', 'macro',
                 'macro')
        temperatures = (300., 600., 900.,
                        520., 600.,
                        300., 600., 900.,
                        520., 600)
        methods = 2 * (3 * ('nearest',) + 2 * ('interpolation',))
        analyt_interp = 10 * [None]
        analyt_interp[3] = (600. - 520.) / 300.
        analyt_interp[8] = (600. - 520.) / 300.
        try:
            if (os.path.isdir("./temp")):
                shutil.rmtree("./temp")
            os.mkdir("temp")
            os.chdir(os.path.join(base_dir, "temp"))
            for cs, t, m, ai in zip(types, temperatures, methods, analyt_interp):
                if (cs == 'macro'):
                    build_inf_model(['macro'], '../macro_2g.h5', t, m)
                else:
                    build_inf_model(names, '../micro_2g.h5', t, m)
                if not ai:
                    kanalyt = analytical_solution_2g_therm(macro_xs[t])
                else:
                    kanalyt = analytical_solution_2g_therm(macro_xs[300],
                                                           macro_xs[600], ai)
                self._run_openmc()
                self._test_output_created()
                string = "{}, method: {}, t: {}, {}kanalyt\n{:12.6E}\n"
                results = string.format(cs, m, t, self._get_results(), kanalyt)
                overall_results.append(results)
            os.chdir(base_dir)
            self._write_results("".join(overall_results))
            self._compare_results()
        finally:
            os.chdir(base_dir)
            shutil.copyfile("results_test.dat", "results_true.dat")
            if (os.path.isdir("./temp")):
                shutil.rmtree("./temp")
            self._cleanup()
            for f in ['micro_2g.h5', 'macro_2g.h5']:
                if os.path.exists(f):
                    os.remove(f)


def test_mg_temperature():
    harness = MgTemperatureTestHarness('statepoint.200.h5')
    harness.main()
