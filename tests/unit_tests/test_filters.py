from math import sqrt, pi

import openmc
from pytest import fixture, approx


@fixture(scope='module')
def box_model():
    model = openmc.model.Model()
    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    m.set_density('g/cm3', 1.0)

    box = openmc.model.rectangular_prism(10., 10., boundary_type='vacuum')
    c = openmc.Cell(fill=m, region=box)
    model.geometry.root_universe = openmc.Universe(cells=[c])

    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.inactive = 0
    model.settings.source = openmc.Source(space=openmc.stats.Point())
    return model


def test_cell_instance():
    c1 = openmc.Cell()
    c2 = openmc.Cell()
    f = openmc.CellInstanceFilter([(c1, 0), (c1, 1), (c1, 2), (c2, 0), (c2, 1)])

    # Make sure __repr__ works
    repr(f)

    # to_xml_element()
    elem = f.to_xml_element()
    assert elem.tag == 'filter'
    assert elem.attrib['type'] == 'cellinstance'
    bins = [int(x) for x in elem.find('bins').text.split()]
    assert all(x == c1.id for x in bins[:6:2])
    assert all(x == c2.id for x in bins[6::2])

    # get_pandas_dataframe()
    df = f.get_pandas_dataframe(f.num_bins, 1)
    cells = df['cellinstance', 'cell']
    instances = df['cellinstance', 'instance']
    assert cells.apply(lambda x: x in (c1.id, c2.id)).all()
    assert instances.apply(lambda x: x in (0, 1, 2)).all()


def test_legendre():
    n = 5
    f = openmc.LegendreFilter(n)
    assert f.order == n
    assert f.bins[0] == 'P0'
    assert f.bins[-1] == 'P5'
    assert len(f.bins) == n + 1

    # Make sure __repr__ works
    repr(f)

    # to_xml_element()
    elem = f.to_xml_element()
    assert elem.tag == 'filter'
    assert elem.attrib['type'] == 'legendre'
    assert elem.find('order').text == str(n)


def test_spatial_legendre():
    n = 5
    axis = 'x'
    f = openmc.SpatialLegendreFilter(n, axis, -10., 10.)
    assert f.order == n
    assert f.axis == axis
    assert f.minimum == -10.
    assert f.maximum == 10.
    assert f.bins[0] == 'P0'
    assert f.bins[-1] == 'P5'
    assert len(f.bins) == n + 1

    # Make sure __repr__ works
    repr(f)

    # to_xml_element()
    elem = f.to_xml_element()
    assert elem.tag == 'filter'
    assert elem.attrib['type'] == 'spatiallegendre'
    assert elem.find('order').text == str(n)
    assert elem.find('axis').text == str(axis)


def test_spherical_harmonics():
    n = 3
    f = openmc.SphericalHarmonicsFilter(n)
    f.cosine = 'particle'
    assert f.order == n
    assert f.bins[0] == 'Y0,0'
    assert f.bins[-1] == 'Y{0},{0}'.format(n)
    assert len(f.bins) == (n + 1)**2

    # Make sure __repr__ works
    repr(f)

    # to_xml_element()
    elem = f.to_xml_element()
    assert elem.tag == 'filter'
    assert elem.attrib['type'] == 'sphericalharmonics'
    assert elem.attrib['cosine'] == f.cosine
    assert elem.find('order').text == str(n)


def test_zernike():
    n = 4
    f = openmc.ZernikeFilter(n, 0., 0., 1.)
    assert f.order == n
    assert f.bins[0] == 'Z0,0'
    assert f.bins[-1] == 'Z{0},{0}'.format(n)
    assert len(f.bins) == (n + 1)*(n + 2)//2

    # Make sure __repr__ works
    repr(f)

    # to_xml_element()
    elem = f.to_xml_element()
    assert elem.tag == 'filter'
    assert elem.attrib['type'] == 'zernike'
    assert elem.find('order').text == str(n)

def test_zernike_radial():
    n = 4
    f = openmc.ZernikeRadialFilter(n, 0., 0., 1.)
    assert f.order == n
    assert f.bins[0] == 'Z0,0'
    assert f.bins[-1] == 'Z{},0'.format(n)
    assert len(f.bins) == n//2 + 1

    # Make sure __repr__ works
    repr(f)

    # to_xml_element()
    elem = f.to_xml_element()
    assert elem.tag == 'filter'
    assert elem.attrib['type'] == 'zernikeradial'
    assert elem.find('order').text == str(n)


def test_first_moment(run_in_tmpdir, box_model):
    plain_tally = openmc.Tally()
    plain_tally.scores = ['flux', 'scatter']

    # Create tallies with expansion filters
    leg_tally = openmc.Tally()
    leg_tally.filters = [openmc.LegendreFilter(3)]
    leg_tally.scores = ['scatter']
    leg_sptl_tally = openmc.Tally()
    leg_sptl_tally.filters = [openmc.SpatialLegendreFilter(3, 'x', -5., 5.)]
    leg_sptl_tally.scores = ['scatter']
    sph_scat_filter = openmc.SphericalHarmonicsFilter(5)
    sph_scat_filter.cosine = 'scatter'
    sph_scat_tally = openmc.Tally()
    sph_scat_tally.filters = [sph_scat_filter]
    sph_scat_tally.scores = ['scatter']
    sph_flux_filter = openmc.SphericalHarmonicsFilter(5)
    sph_flux_filter.cosine = 'particle'
    sph_flux_tally = openmc.Tally()
    sph_flux_tally.filters = [sph_flux_filter]
    sph_flux_tally.scores = ['flux']
    zernike_tally = openmc.Tally()
    zernike_tally.filters = [openmc.ZernikeFilter(3, r=10.)]
    zernike_tally.scores = ['scatter']

    # Add tallies to model and ensure they all use the same estimator
    box_model.tallies = [plain_tally, leg_tally, leg_sptl_tally,
                         sph_scat_tally, sph_flux_tally, zernike_tally]
    for t in box_model.tallies:
        t.estimator = 'analog'

    sp_name = box_model.run()

    # Check that first moment matches the score from the plain tally
    with openmc.StatePoint(sp_name) as sp:
        # Get scores from tally without expansion filters
        flux, scatter = sp.tallies[plain_tally.id].mean.ravel()

        # Check that first moment matches
        first_score = lambda t: sp.tallies[t.id].mean.ravel()[0]
        assert first_score(leg_tally) == scatter
        assert first_score(leg_sptl_tally) == scatter
        assert first_score(sph_scat_tally) == scatter
        assert first_score(sph_flux_tally) == approx(flux)
        assert first_score(zernike_tally) == approx(scatter)
