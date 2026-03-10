from itertools import product, permutations

import openmc
import numpy as np

import pytest

pitch = 1.25    

@pytest.fixture(params=['x','y'])
def model(request):
    openmc.reset_auto_ids()

    orientation = request.param

    water = openmc.Material(name='water')
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.set_density('g/cc', 1.0)

    outer = openmc.Cell(fill=water, cell_id=100)
    cell10 = openmc.Cell(fill=water, cell_id=10)
    cell00 = openmc.Cell(fill=water, cell_id=0)
    cell01 = openmc.Cell(fill=water, cell_id=1)
    cell02 = openmc.Cell(fill=water, cell_id=2)
    cell03 = openmc.Cell(fill=water, cell_id=3)
    cell04 = openmc.Cell(fill=water, cell_id=4)
    cell05 = openmc.Cell(fill=water, cell_id=5)

    univ10 = openmc.Universe(cells=[cell10])
    univ00 = openmc.Universe(cells=[cell00])
    univ01 = openmc.Universe(cells=[cell01])
    univ02 = openmc.Universe(cells=[cell02])
    univ03 = openmc.Universe(cells=[cell03])
    univ04 = openmc.Universe(cells=[cell04])
    univ05 = openmc.Universe(cells=[cell05])

    plane1 = openmc.ZPlane(-10.0, boundary_type='vacuum')
    plane2 = openmc.ZPlane(10.0, boundary_type='vacuum')

    lat = openmc.HexLattice()
    lat.center = (0., 0.)
    lat.pitch = (pitch,)
    lat.universes = [[univ00, univ01, univ02, univ03, univ04, univ05], [univ10]]
    lat.outer = openmc.Universe(cells=[outer])
    lat.orientation = orientation

    hex_prism = openmc.model.HexagonalPrism(
                    edge_length=2*pitch,
                    orientation=orientation,
                    boundary_type='vacuum'
                )
    cell = openmc.Cell(region=-hex_prism & +plane1 & -plane2, fill=lat)

    geom = openmc.Geometry([cell])

    source = openmc.IndependentSource()
    source.space = openmc.stats.Point()
    source.energy = openmc.stats.Discrete([10000], [1.0])

    settings = openmc.Settings()
    settings.particles = 2000
    settings.batches = 10
    settings.run_mode = 'fixed source'

    # build
    mesh = openmc.HexagonalMesh(
        z_grid=[-10.0, 10.0],
        num_rings = 2,
        pitch = pitch,
        orientation = orientation,
    )
    tally = openmc.Tally()

    mesh_filter = openmc.MeshFilter(mesh)
    cell_filter = openmc.CellFilter([cell10, cell00, cell01, cell02, cell03, cell04, cell05])
    tally.filters.append(mesh_filter)
    tally.filters.append(cell_filter)

    tally.scores.append("total")

    tallies = openmc.Tallies([tally])

    return openmc.Model(geometry=geom, settings=settings, tallies=tallies)


@pytest.mark.parametrize("estimator", ["collision", "tracklength"])
def test_correct_locations(model, run_in_tmpdir, estimator):
    tally, = model.tallies
    tally.estimator = estimator
    model.run(apply_tally_results=True)
    df = tally.get_pandas_dataframe()
    df = df[df['mean']>0]
    df = df.set_index("cell")["mesh 1"]
    assert len(df) == 7

    if tally.filters[0].mesh.orientation == "y":
        assert tuple(df.loc[10]) == (0,0,0,0)
        assert tuple(df.loc[0]) == (0,-1,1,0)
        assert tuple(df.loc[1]) == (1,-1,0,0)
        assert tuple(df.loc[2]) == (1,0,-1,0)
        assert tuple(df.loc[3]) == (0,1,-1,0)
        assert tuple(df.loc[4]) == (-1,1,0,0)
        assert tuple(df.loc[5]) == (-1,0,1,0)
    else:
        assert tuple(df.loc[10]) == (0,0,0,0)
        assert tuple(df.loc[0]) == (1,0,-1,0)
        assert tuple(df.loc[1]) == (0,1,-1,0)
        assert tuple(df.loc[2]) == (-1,1,0,0)
        assert tuple(df.loc[3]) == (-1,0,1,0)
        assert tuple(df.loc[4]) == (0,-1,1,0)
        assert tuple(df.loc[5]) == (1,-1,0,0)

def test_meshsurface(model):
    tally, = model.tallies
    orientation = tally.filters[0].mesh.orientation
    
    tally = openmc.Tally()
    mesh = openmc.HexagonalMesh(
    z_grid=[-10.0, 10.0],
    num_rings = 2,
    pitch = pitch,
    orientation = orientation,
    )

    tally.filters = [openmc.MeshSurfaceFilter(mesh)]
    tally.scores = ['current']

    model.tallies = [tally]

    model.run(apply_tally_results=True)

    df = tally.get_pandas_dataframe()

    z_max_in = df['mesh 2', 'surf']=='z-max in'
    z_min_in = df['mesh 2', 'surf']=='z-min in'
    
    assert np.all(df[z_max_in & z_min_in]['mean']==0.0)
