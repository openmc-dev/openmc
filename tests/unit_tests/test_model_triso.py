#!/usr/bin/env python

from math import pi

import numpy as np
from numpy.linalg import norm
import openmc
import openmc.model
import pytest
import scipy.spatial


_PACKING_FRACTION = 0.35
_RADIUS = 1.
_PARAMS = [{'shape' : 'cube', 'length' : 20., 'radius' : 0.},
           {'shape' : 'cylinder', 'length' : 10., 'radius' : 10.},
           {'shape' : 'sphere', 'length' : 0., 'radius' : 10.}]


@pytest.fixture(scope='module', params=_PARAMS)
def domain(request):
    return request.param


@pytest.fixture(scope='module')
def trisos(domain):
    return openmc.model.pack_trisos(
        radius=_RADIUS,
        fill=openmc.Universe(),
        domain_shape=domain['shape'],
        domain_length=domain['length'],
        domain_radius=domain['radius'],
        domain_center=[0., 0., 0.],
        initial_packing_fraction=0.2,
        packing_fraction=_PACKING_FRACTION
    )


def test_overlap(trisos):
    """Check that no TRISO particles overlap."""
    tree = scipy.spatial.cKDTree([t.center for t in trisos])
    r = min(tree.query([t.center for t in trisos], k=2)[0][:,1])/2
    assert r > _RADIUS or r == pytest.approx(_RADIUS)


def test_contained(trisos, domain):
    """Make sure all particles are entirely contained within the domain."""
    if domain['shape'] == 'cube':
        x = 2*(max(np.hstack([abs(t.center) for t in trisos])) + _RADIUS)
        assert x < domain['length'] or x == pytest.approx(domain['length'])

    elif domain['shape'] == 'cylinder':
        r = max([norm(t.center[0:2]) for t in trisos]) + _RADIUS
        z = 2*(max([abs(t.center[2]) for t in trisos]) + _RADIUS)
        assert r < domain['radius'] or r == pytest.approx(domain['radius'])
        assert z < domain['length'] or z == pytest.approx(domain['length'])

    elif domain['shape'] == 'sphere':
        r = max([norm(t.center) for t in trisos]) + _RADIUS
        assert r < domain['radius'] or r == pytest.approx(domain['radius'])


def test_packing_fraction(trisos, domain):
    """Check that the actual PF is close to the requested PF."""
    if domain['shape'] == 'cube':
        volume = domain['length']**3

    elif domain['shape'] == 'cylinder':
        volume = domain['length']*pi*domain['radius']**2

    elif domain['shape'] == 'sphere':
        volume = 4/3*pi*domain['radius']**3

    pf = len(trisos)*4/3*pi*_RADIUS**3/volume
    assert pf == pytest.approx(_PACKING_FRACTION, rel=1e-2)
