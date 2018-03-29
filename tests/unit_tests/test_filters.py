import openmc


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
    assert f.bins[-1] == 'Y{0},{0}'.format(n, n)
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
