import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    m = openmc.Material()
    m.set_density('g/cm3', 10.0)
    m.add_nuclide('Am241', 1.0)
    model.materials.append(m)

    s = openmc.Sphere(r=100.0, boundary_type='vacuum')
    c = openmc.Cell(fill=m, region=-s)
    model.geometry = openmc.Geometry([c])

    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000

    # Define Am242m / Am242 branching ratio from ENDF/B-VII.1 data.
    x = [1e-5, 3.69e-1, 1e3, 1e5, 6e5, 1e6, 2e6, 4e6, 3e7]
    y = [0.1, 0.1, 0.1333, 0.158, 0.18467, 0.25618, 0.4297, 0.48, 0.48]

    # Make an EnergyFunctionFilter directly from the x and y lists.
    filt1 = openmc.EnergyFunctionFilter(x, y)

    # check interpolatoin property setter
    assert filt1.interpolation == 'linear-linear'

    with pytest.raises(ValueError):
        filt1.interpolation = '🥏'

    # Also make a filter with the .from_tabulated1d constructor.  Make sure
    # the filters are identical.
    tab1d = openmc.data.Tabulated1D(x, y)
    filt2 = openmc.EnergyFunctionFilter.from_tabulated1d(tab1d)
    assert filt1 == filt2, 'Error with the .from_tabulated1d constructor'

    filt3 = openmc.EnergyFunctionFilter(x, y)
    filt3.interpolation = 'log-log'

    filt4 = openmc.EnergyFunctionFilter(x, y)
    filt4.interpolation = 'linear-log'

    filt5 = openmc.EnergyFunctionFilter(x, y)
    filt5.interpolation = 'log-linear'

    # define a trapezoidal function for comparison
    x = [0.0, 5e6, 1e7, 1.5e7]
    y = [0.2, 0.7, 0.7, 0.2]

    filt6 = openmc.EnergyFunctionFilter(x, y)

    filt7 = openmc.EnergyFunctionFilter(x, y)
    filt7.interpolation = 'quadratic'

    filt8 = openmc.EnergyFunctionFilter(x, y)
    filt8.interpolation = 'cubic'

    filt9 = openmc.EnergyFunctionFilter(x, y)
    filt9.interpolation = 'histogram'

    filters = [filt1, filt3, filt4, filt5, filt6, filt7, filt8, filt9]
    # Make tallies
    tallies = [openmc.Tally() for _ in range(len(filters) + 1)]
    for t in tallies:
        t.scores = ['(n,gamma)']
        t.nuclides = ['Am241']

    for t, f in zip(tallies[1:], filters):
        t.filters = [f]

    model.tallies.extend(tallies)

    interpolation_vals = \
    list(openmc.EnergyFunctionFilter.INTERPOLATION_SCHEMES.keys())
    for i_val in interpolation_vals:
        # breakpoint here is fake and unused
        t1d = openmc.data.Tabulated1D(x,
                                      y,
                                      breakpoints=[1],
                                      interpolation=[i_val])
        f = openmc.EnergyFunctionFilter.from_tabulated1d(t1d)
        assert f.interpolation == \
        openmc.EnergyFunctionFilter.INTERPOLATION_SCHEMES[i_val]

    return model


class FilterEnergyFunHarness(PyAPITestHarness):
    def _get_results(self):
        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        dataframes_string = ""
        # Use tally arithmetic to compute the branching ratio.
        br_tally = sp.tallies[2] / sp.tallies[1]
        dataframes_string += br_tally.get_pandas_dataframe().to_string() + '\n'

        for t_id in (3, 4, 5, 6):
            ef_tally = sp.tallies[t_id]
            dataframes_string += ef_tally.get_pandas_dataframe().to_string() + '\n'

        # Output the tally in a Pandas DataFrame.
        return dataframes_string

    def _compare_results(self):
        super()._compare_results()

        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        # statepoint file round-trip checks

        # linear-linear interpolation tally
        sp_lin_lin_tally = sp.get_tally(id=2)
        sp_lin_lin_filt = sp_lin_lin_tally.find_filter(openmc.EnergyFunctionFilter)

        model_lin_lin_tally = self._model.tallies[1]
        model_lin_lin_filt = model_lin_lin_tally.find_filter(openmc.EnergyFunctionFilter)

        assert sp_lin_lin_filt.interpolation == 'linear-linear'
        assert all(sp_lin_lin_filt.energy == model_lin_lin_filt.energy)
        assert all(sp_lin_lin_filt.y == model_lin_lin_filt.y)

        # log-log interpolation tally
        sp_log_log_tally = sp.get_tally(id=3)
        sp_log_log_filt = sp_log_log_tally.find_filter(openmc.EnergyFunctionFilter)

        model_log_log_tally = self._model.tallies[2]
        model_log_log_filt = model_log_log_tally.find_filter(openmc.EnergyFunctionFilter)

        assert sp_log_log_filt.interpolation == 'log-log'
        assert all(sp_log_log_filt.energy == model_log_log_filt.energy)
        assert all(sp_log_log_filt.y == model_log_log_filt.y)

        # because the values of y are monotonically increasing,
        # we expect the log-log tally to have a higher value
        assert all(sp_lin_lin_tally.mean < sp_log_log_tally.mean)

        sp_lin_log_tally = self._model.tallies[3]
        sp_lin_log_filt = sp_lin_log_tally.find_filter(openmc.EnergyFunctionFilter)
        assert sp_lin_log_filt.interpolation == 'linear-log'

        sp_log_lin_tally = self._model.tallies[4]
        sp_log_lin_filt = sp_log_lin_tally.find_filter(openmc.EnergyFunctionFilter)
        assert sp_log_lin_filt.interpolation == 'log-linear'

        # check that the cubic interpolation provides a higher value
        # than linear-linear
        contrived_lin_lin_tally = sp.get_tally(id=6)
        contrived_quadratic_tally = sp.get_tally(id=7)
        contrived_cubic_tally = sp.get_tally(id=8)

        assert all(contrived_lin_lin_tally.mean < contrived_quadratic_tally.mean)
        assert all(contrived_lin_lin_tally.mean < contrived_cubic_tally.mean)

        # check that the histogram tally is less than the quadratic/cubic interpolations
        histogram_tally = sp.get_tally(id=9)
        assert all(histogram_tally.mean < contrived_quadratic_tally.mean)
        assert all(histogram_tally.mean < contrived_cubic_tally.mean)

def test_filter_energyfun(model):
    harness = FilterEnergyFunHarness('statepoint.5.h5', model)
    harness.main()

if __name__ == '__main__':
    pytest.main(['-v', __file__])