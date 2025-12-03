import openmc
import pandas as pd


def test_distribcell_filter_apply_tally_results(run_in_tmpdir):
    mat = openmc.Material()
    mat.add_nuclide("U235", 1.0)
    mat.set_density("g/cm3", 1.0)

    # Define 2x2 lattice with a cylinder in each universe
    cyl = openmc.ZCylinder(r=1.0)
    cell1 = openmc.Cell(fill=mat, region=-cyl)
    cell2 = openmc.Cell(fill=None, region=+cyl)
    univ = openmc.Universe(cells=[cell1, cell2])
    lattice = openmc.RectLattice()
    lattice.lower_left = (-3.0, -3.0)
    lattice.pitch = (3.0, 3.0)
    lattice.universes = [[univ, univ], [univ, univ]]
    box = openmc.model.RectangularPrism(6., 6., boundary_type='reflective')
    root_cell = openmc.Cell(region=-box, fill=lattice)
    geometry = openmc.Geometry([root_cell])

    # Create model and add tally with distribcell filter
    model = openmc.Model(geometry)
    model.settings.batches = 10
    model.settings.particles = 1000
    tally = openmc.Tally()
    distribcell_filter = openmc.DistribcellFilter(cell1)
    tally.filters = [distribcell_filter]
    tally.scores = ['flux']
    model.tallies = [tally]

    # Run OpenMC and apply tally results
    model.run(apply_tally_results=True)

    # Check that mean and standard deviation are available on tally
    assert tally.mean.shape == (4, 1, 1)
    assert tally.std_dev.shape == (4, 1, 1)

    # Make sure paths attribute on filter is correct
    assert distribcell_filter.paths == [
        'u3->c3->l2(0,0)->u1->c1',
        'u3->c3->l2(1,0)->u1->c1',
        'u3->c3->l2(0,1)->u1->c1',
        'u3->c3->l2(1,1)->u1->c1',
    ]

    # Check that we can get a DataFrame from the tally
    df = tally.get_pandas_dataframe()
    assert isinstance(df, pd.DataFrame)
