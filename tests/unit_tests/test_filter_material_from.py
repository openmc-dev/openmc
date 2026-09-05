"""MaterialFromFilter must decompose scores the same way CellFromFilter does.

The geometry below gives each cell its own material, so the two decompositions
are one-to-one and must agree bin for bin. Both must also sum to the
undecomposed total. These are per-history identities within a single run, not
statistical statements, so they hold to round-off.

Note the volumetric tests below score `total`, not `flux`. Scoring `flux` here
would silently turn them into surface tallies measuring a different quantity.
"""

import numpy as np
import openmc
import pytest


@pytest.fixture
def model():
    openmc.reset_auto_ids()
    model = openmc.Model()

    water = openmc.Material(name='water')
    water.set_density('g/cm3', 1.0)
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)

    iron = openmc.Material(name='iron')
    iron.set_density('g/cm3', 7.87)
    iron.add_nuclide('Fe56', 1.0)

    inner_surf = openmc.Sphere(r=5.0)
    outer_surf = openmc.Sphere(r=15.0, boundary_type='vacuum')
    inner = openmc.Cell(name='inner', fill=water, region=-inner_surf)
    outer = openmc.Cell(name='outer', fill=iron,
                        region=+inner_surf & -outer_surf)
    model.geometry = openmc.Geometry([inner, outer])

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 10
    model.settings.particles = 2000
    # Source at the centre, so every particle scoring in the outer cell got
    # there by crossing the surface out of the inner cell. The shell is thick
    # enough in mean free paths that most of them scatter there repeatedly,
    # which is what exercises the post-collision reset.
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0.0, 0.0, 0.0)),
        energy=openmc.stats.delta_function(2.0e6),
    )

    in_outer = openmc.CellFilter([outer])

    total = openmc.Tally(name='total')
    total.filters = [in_outer]
    total.scores = ['total']

    by_cell = openmc.Tally(name='by cell')
    by_cell.filters = [in_outer, openmc.CellFromFilter([inner, outer])]
    by_cell.scores = ['total']

    by_material = openmc.Tally(name='by material')
    by_material.filters = [in_outer, openmc.MaterialFromFilter([water, iron])]
    by_material.scores = ['total']

    model.tallies = openmc.Tallies([total, by_cell, by_material])
    return model


def test_material_from_matches_cell_from(model, run_in_tmpdir):
    with openmc.StatePoint(model.run()) as sp:
        total = sp.get_tally(name='total').mean.ravel()
        by_cell = sp.get_tally(name='by cell').mean.ravel()
        by_material = sp.get_tally(name='by material').mean.ravel()

    # Neither decomposition may lose scores. Before the fix the material
    # decomposition fell short of the total by everything scored after a
    # particle's first collision in the cell.
    assert by_cell.sum() == pytest.approx(total.sum(), rel=1e-9)
    assert by_material.sum() == pytest.approx(total.sum(), rel=1e-9)

    # One material per cell, so the two decompositions are the same partition
    np.testing.assert_allclose(by_material, by_cell, rtol=1e-9)


def test_from_bin_survives_repeated_collisions(model, run_in_tmpdir):
    """The 'from inner' bin holds all of the outer cell's reaction rate.

    A particle that has crossed from the inner cell keeps `cell_last` 
    pointing at that cell for as long as it stays in the outer one,
    however many times it collides. The 'from outer' bin can only be
    populated by a particle that leaves the outer cell and comes back, which
    this geometry does not allow since the outer boundary is vacuum and the
    inner cell is convex. So the 'from inner' bin should carry the whole total.
    """
    with openmc.StatePoint(model.run()) as sp:
        total = sp.get_tally(name='total').mean.ravel()
        by_material = sp.get_tally(name='by material').mean.ravel()

    from_water, from_iron = by_material
    assert from_water == pytest.approx(total.sum(), rel=1e-9)
    assert from_iron == pytest.approx(0.0, abs=1e-12)


@pytest.fixture
def distribcell_model():
    """A 2x2 lattice of one repeated cell carrying four distributed materials.

    Every instance has the same composition, so the four "from" bins are
    physically equivalent and should come out equal within statistics. The
    point is that they are non-zero at all: the material fill is resolved per
    instance, so deriving the previous material from cell_last() alone would
    attribute every score to the first material.
    """
    openmc.reset_auto_ids()
    model = openmc.Model()

    def fuel():
        m = openmc.Material()
        m.set_density('g/cm3', 10.0)
        m.add_nuclide('U238', 1.0)
        m.add_nuclide('O16', 2.0)
        return m

    fuels = [fuel() for _ in range(4)]

    reflector = openmc.Material(name='reflector')
    reflector.set_density('g/cm3', 1.0)
    reflector.add_nuclide('H1', 2.0)
    reflector.add_nuclide('O16', 1.0)

    # One cell, filled with a list of materials, repeated at every lattice
    # position, so it has four distribcell instances
    pin = openmc.Cell(name='pin', fill=fuels)
    pin_universe = openmc.Universe(cells=[pin])

    lattice = openmc.RectLattice()
    lattice.lower_left = (-2.0, -2.0)
    lattice.pitch = (2.0, 2.0)
    lattice.universes = [[pin_universe, pin_universe],
                         [pin_universe, pin_universe]]

    box = openmc.model.RectangularPrism(4.0, 4.0)
    sphere = openmc.Sphere(r=20.0, boundary_type='vacuum')
    lattice_cell = openmc.Cell(name='lattice', fill=lattice,
                               region=-box & -sphere)
    outside = openmc.Cell(name='outside', fill=reflector,
                          region=+box & -sphere)
    model.geometry = openmc.Geometry([lattice_cell, outside])

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 10
    model.settings.particles = 4000
    # Born throughout the lattice, so all four instances are populated
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-2.0, -2.0, -1.0), (2.0, 2.0, 1.0)),
        energy=openmc.stats.delta_function(2.0e6),
    )

    in_outside = openmc.CellFilter([outside])

    total = openmc.Tally(name='total')
    total.filters = [in_outside]
    total.scores = ['total']

    by_material = openmc.Tally(name='by material')
    by_material.filters = [in_outside, openmc.MaterialFromFilter(fuels)]
    by_material.scores = ['total']

    model.tallies = openmc.Tallies([total, by_material])
    return model


def test_distributed_materials_resolve_per_instance(distribcell_model,
                                                    run_in_tmpdir):
    with openmc.StatePoint(distribcell_model.run()) as sp:
        total = sp.get_tally(name='total').mean.ravel().sum()
        by_material = sp.get_tally(name='by material').mean.ravel()

    # Everything scoring in the reflector arrived from some lattice instance,
    # so no score may go unattributed
    assert by_material.sum() == pytest.approx(total, rel=1e-9)

    # All four instances contribute. Resolving the fill without the instance
    # would put the whole total in the first bin and leave the rest at zero.
    assert (by_material > 0.0).all()

    # The instances are geometrically and compositionally equivalent, so the
    # four bins should be comparable; the bound is loose enough to be
    # insensitive to counting statistics
    assert by_material.max() < 1.5 * by_material.min()


def test_partial_current_decompositions_agree(model, run_in_tmpdir):
    """The surface-tally use of these filters is unaffected.

    A `flux` score with a "from" filter is a cell-to-cell partial current, not a
    volumetric flux: it is scored in event_cross_surface() right after the
    crossing.
    """
    inner, outer = model.geometry.get_all_cells().values()
    water = inner.fill
    in_outer = openmc.CellFilter([outer])

    from_cell = openmc.Tally(name='current from cell')
    from_cell.filters = [in_outer, openmc.CellFromFilter([inner])]
    from_cell.scores = ['flux']

    from_mat = openmc.Tally(name='current from material')
    from_mat.filters = [in_outer, openmc.MaterialFromFilter([water])]
    from_mat.scores = ['flux']

    model.tallies = openmc.Tallies([from_cell, from_mat])

    with openmc.StatePoint(model.run()) as sp:
        by_cell = sp.get_tally(name='current from cell').mean.ravel()
        by_material = sp.get_tally(name='current from material').mean.ravel()

    # Every particle from a point source at the centre crosses into the outer
    # cell, so this is not a vacuous comparison of two zeros
    assert by_cell.sum() > 0.0
    np.testing.assert_allclose(by_material, by_cell, rtol=1e-9)
