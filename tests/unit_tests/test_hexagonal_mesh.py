from itertools import product, permutations

import openmc
import numpy as np

import pytest

@pytest.fixture()
def model():
    openmc.reset_auto_ids()

    water = openmc.Material(name='water')
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.set_density('g/cc', 1.0)

    outer = openmc.Cell(fill=water)
    cell10 = openmc.Cell(fill=water)
    cell00 = openmc.Cell(fill=water)
    cell01 = openmc.Cell(fill=water)
    cell02 = openmc.Cell(fill=water)
    cell03 = openmc.Cell(fill=water)
    cell04 = openmc.Cell(fill=water)
    cell05 = openmc.Cell(fill=water)

    univ10 = openmc.Universe(cells=[cell10])
    univ00 = openmc.Universe(cells=[cell00])
    univ01 = openmc.Universe(cells=[cell01])
    univ02 = openmc.Universe(cells=[cell02])
    univ03 = openmc.Universe(cells=[cell03])
    univ04 = openmc.Universe(cells=[cell04])
    univ05 = openmc.Universe(cells=[cell05])

    plane1 = openmc.ZPlane(-10.0, boundary_type='vacuum')
    plane2 = openmc.ZPlane(10.0, boundary_type='vacuum')

    pitch = 1.25

    lat = openmc.HexLattice()
    lat.center = (0., 0.)
    lat.pitch = (pitch,)
    lat.universes = [[univ00, univ01, univ02, univ03, univ04, univ05], [univ10]]
    lat.outer = openmc.Universe(cells=[outer])
    lat.orientation = 'y'

    hex_prism = openmc.model.HexagonalPrism(
                    edge_length=2*pitch,
                    orientation='y',
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
    )
    tally = openmc.Tally()

    mesh_filter = openmc.MeshFilter(mesh)
    cell_filter = openmc.CellFilter([cell10, cell00, cell01, cell02, cell03, cell04, cell05])
    tally.filters.append(mesh_filter)
    tally.filters.append(cell_filter)

    tally.scores.append("total")
    tally.estimator = "collision"

    tallies = openmc.Tallies([tally])

    return openmc.Model(geometry=geom, settings=settings, tallies=tallies)

def test_correct_locations(model, run_in_tmpdir):
    tally, = model.tallies
    model.run(apply_tally_results=True)
    df = tally.get_pandas_dataframe()
    df = df[df['mean']>0]

    assert len(df) == 7
    
