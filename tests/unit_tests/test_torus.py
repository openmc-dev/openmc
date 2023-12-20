from itertools import combinations
from random import uniform
import openmc
import pytest


def get_torus_keff(cls, R, r, center=(0, 0, 0)):
    model = openmc.Model()
    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.set_density('g/cm3', 10.0)

    x, y, z = center
    torus = cls(x0=x, y0=y, z0=z, a=R, b=r, c=r)
    sphere = openmc.Sphere(x0=x, y0=y, z0=z, r=R + r + 1, boundary_type="vacuum")
    torus_cell = openmc.Cell(fill=mat, region=-torus)
    outer_cell = openmc.Cell(region=+torus & -sphere)
    model.geometry = openmc.Geometry([torus_cell, outer_cell])

    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point(center))
    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000

    sp_path = model.run()
    with openmc.StatePoint(sp_path) as sp:
        return sp.keff


@pytest.mark.parametrize("R,r", [(2.1, 2.0), (3.0, 1.0)])
def test_torus_keff(R, r, run_in_tmpdir):
    random_point = lambda: (uniform(-5, 5), uniform(-5, 5), uniform(-5, 5))
    keffs = [
        get_torus_keff(openmc.XTorus, R, r),
        get_torus_keff(openmc.XTorus, R, r, random_point()),
        get_torus_keff(openmc.YTorus, R, r),
        get_torus_keff(openmc.YTorus, R, r, random_point()),
        get_torus_keff(openmc.ZTorus, R, r),
        get_torus_keff(openmc.ZTorus, R, r, random_point())
    ]

    # For each combination of keff values, their difference should be within
    # uncertainty (3 std dev)
    for k1, k2 in combinations(keffs, 2):
        print(k1, k2)
        diff = k1 - k2
        assert abs(diff.n) < 3*diff.s
