from ctypes import c_int

import h5py
import numpy as np
import openmc
import openmc.lib
from openmc.exceptions import InvalidArgumentError, InvalidTypeError
import pytest


def test_lifetime_moment_filter():
    n = 3
    f = openmc.LifetimeMomentFilter(n)
    assert f.order == n
    assert f.bins == ['tau^0', 'tau^1', 'tau^2', 'tau^3']
    assert f.num_bins == n + 1
    assert f.shape == (n + 1,)

    # Changing the order updates the bins
    f.order = 1
    assert f.bins == ['tau^0', 'tau^1']
    f.order = n

    # Make sure __repr__ works
    repr(f)

    # Equality and hashing depend on the order only
    assert f == openmc.LifetimeMomentFilter(n)
    assert f != openmc.LifetimeMomentFilter(n + 1)
    assert f != openmc.LegendreFilter(n)
    assert hash(f) == hash(openmc.LifetimeMomentFilter(n))

    # to_xml_element()
    elem = f.to_xml_element()
    assert elem.tag == 'filter'
    assert elem.attrib['type'] == 'lifetimemoment'
    assert elem.attrib['id'] == str(f.id)
    assert elem.find('order').text == str(n)
    assert elem.find('bins') is None

    # from_xml_element()
    new_f = openmc.Filter.from_xml_element(elem)
    assert isinstance(new_f, openmc.LifetimeMomentFilter)
    assert new_f.id == f.id
    assert new_f.order == n
    assert new_f.bins == f.bins

    # Merging keeps the higher order
    merged = f.merge(openmc.LifetimeMomentFilter(n + 2))
    assert isinstance(merged, openmc.LifetimeMomentFilter)
    assert merged.order == n + 2


def test_lifetime_moment_filter_invalid_order():
    with pytest.raises(ValueError):
        openmc.LifetimeMomentFilter(-1)
    with pytest.raises(TypeError):
        openmc.LifetimeMomentFilter(1.5)
    with pytest.raises(TypeError):
        openmc.LifetimeMomentFilter('2')

    f = openmc.LifetimeMomentFilter(0)
    assert f.bins == ['tau^0']
    with pytest.raises(ValueError):
        f.order = -2
    assert f.order == 0


def test_lifetime_moment_filter_from_hdf5(run_in_tmpdir):
    # Mimic the group that is written to a statepoint for this filter
    with h5py.File('filter.h5', 'w') as fh:
        group = fh.create_group('tallies/filters/filter 17')
        group.create_dataset('type', data=np.bytes_('lifetimemoment'))
        group.create_dataset('n_bins', data=3)
        group.create_dataset('order', data=2)

    with h5py.File('filter.h5', 'r') as fh:
        f = openmc.Filter.from_hdf5(fh['tallies/filters/filter 17'])
    assert isinstance(f, openmc.LifetimeMomentFilter)
    assert f.id == 17
    assert f.order == 2
    assert isinstance(f.order, int)
    assert f.bins == ['tau^0', 'tau^1', 'tau^2']


def test_lifetime_moment_filter_dataframe():
    f = openmc.LifetimeMomentFilter(2)

    # Filter with a stride of two (one more filter with two bins after it) in a
    # tally with another filter of two bins before it
    df = f.get_pandas_dataframe(12, 2)
    assert list(df.columns) == ['lifetimemoment']
    expected = ['tau^0', 'tau^0', 'tau^1', 'tau^1', 'tau^2', 'tau^2'] * 2
    assert list(df['lifetimemoment']) == expected


@pytest.fixture
def one_group_model(run_in_tmpdir):
    groups = openmc.mgxs.EnergyGroups([0.0, 20.0e6])
    xsdata = openmc.XSdata('medium', groups)
    xsdata.order = 0
    xsdata.set_total([1.0])
    xsdata.set_absorption([0.25])
    xsdata.set_scatter_matrix([[[0.75]]])
    xsdata.set_inverse_velocity([1.0e-5])
    library = openmc.MGXSLibrary(groups)
    library.add_xsdata(xsdata)
    library.export_to_hdf5('one_group.h5')

    medium = openmc.Material()
    medium.set_density('macro', 1.0)
    medium.add_macroscopic('medium')

    model = openmc.Model()
    model.materials = openmc.Materials([medium])
    model.materials.cross_sections = 'one_group.h5'
    box = openmc.model.RectangularParallelepiped(
        -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, boundary_type='reflective')
    model.geometry = openmc.Geometry([openmc.Cell(fill=medium, region=-box)])
    model.settings.energy_mode = 'multi-group'
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 2
    model.settings.particles = 100
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1.0e6], [1.0]))

    tally = openmc.Tally(tally_id=7)
    tally.filters = [openmc.LifetimeMomentFilter(2, filter_id=11)]
    tally.scores = ['absorption']
    model.tallies = [tally]
    return model


def test_lifetime_moment_filter_lib(one_group_model):
    one_group_model.init_lib(output=False)
    try:
        # The filter is read from XML and the collision estimator is chosen
        tally = openmc.lib.tallies[7]
        assert tally.estimator == 'collision'
        lib_filter = tally.filters[0]
        assert isinstance(lib_filter, openmc.lib.LifetimeMomentFilter)
        assert lib_filter.id == 11
        assert lib_filter.order == 2
        assert lib_filter.n_bins == 3
        assert isinstance(openmc.lib.filters[11],
                          openmc.lib.LifetimeMomentFilter)

        # A filter created through the C API can change its order
        new_filter = openmc.lib.LifetimeMomentFilter(1)
        assert new_filter.order == 1
        assert new_filter.n_bins == 2
        new_filter.order = 4
        assert new_filter.order == 4
        assert new_filter.n_bins == 5

        # Invalid orders and filter types are reported as errors
        with pytest.raises(InvalidArgumentError):
            new_filter.order = -1
        assert new_filter.order == 4
        energy_filter = openmc.lib.EnergyFilter([0.0, 20.0e6])
        with pytest.raises(InvalidTypeError):
            openmc.lib.filter._dll.openmc_lifetime_moment_filter_get_order(
                energy_filter._index, c_int())
    finally:
        one_group_model.finalize_lib()


def test_lifetime_moment_filter_lib_runtime_tally(one_group_model):
    one_group_model.init_lib(output=False)
    try:
        # Tallies created at run time keep the estimator they are given until
        # the simulation is initialized
        runtime = openmc.lib.Tally()
        runtime.filters = [openmc.lib.LifetimeMomentFilter(2)]
        runtime.scores = ['absorption']
        assert runtime.estimator == 'tracklength'
        analog = openmc.lib.Tally()
        analog.filters = [openmc.lib.LifetimeMomentFilter(2)]
        analog.scores = ['absorption']
        analog.estimator = 'analog'

        # The track-length estimator is replaced by a collision estimator, and
        # an analog estimator is kept
        openmc.lib.run(output=False)
        assert runtime.estimator == 'collision'
        assert analog.estimator == 'analog'

        # The run-time tally scores the same events as the tally read from
        # tallies.xml, which also uses a collision estimator. Scores from
        # several OpenMP threads are summed in a run-dependent order, so the
        # two means agree to rounding (observed ~1e-16 relative), not bitwise
        xml_tally = openmc.lib.tallies[7]
        assert xml_tally.estimator == 'collision'
        assert runtime.num_realizations == xml_tally.num_realizations > 0
        np.testing.assert_allclose(runtime.mean, xml_tally.mean,
                                   rtol=1e-12, atol=0.0)
        assert np.all(runtime.mean > 0.0)
    finally:
        one_group_model.finalize_lib()
