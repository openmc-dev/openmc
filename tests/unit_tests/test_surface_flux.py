"""Tests for surface flux tallying via flux score + SurfaceFilter."""

import math
import pytest

import openmc


@pytest.fixture
def two_cell_model():
    """Simple two-cell slab model with a monodirectional fixed source.

    Cell1 occupies x in [-10, 0], cell2 x in [0, 10].  The source fires all
    particles from (-5, 0, 0) in the +x direction with weight 1.  Every
    particle therefore crosses the surface at x=0 from cell1 into cell2 at
    mu = 1 (normal incidence).
    """
    openmc.reset_auto_ids()
    model = openmc.Model()

    xmin = openmc.XPlane(-10.0, boundary_type="vacuum")
    xmid = openmc.XPlane(0.0)
    xmax = openmc.XPlane(10.0, boundary_type="vacuum")
    ymin = openmc.YPlane(-10.0, boundary_type="vacuum")
    ymax = openmc.YPlane(10.0, boundary_type="vacuum")
    zmin = openmc.ZPlane(-10.0, boundary_type="vacuum")
    zmax = openmc.ZPlane(10.0, boundary_type="vacuum")

    cell1 = openmc.Cell(region=+xmin & -xmid & +ymin & -ymax & +zmin & -zmax)
    cell2 = openmc.Cell(region=+xmid & -xmax & +ymin & -ymax & +zmin & -zmax)
    model.geometry = openmc.Geometry([cell1, cell2])

    src = openmc.IndependentSource()
    src.space = openmc.stats.Point((-5.0, 0.0, 0.0))
    src.angle = openmc.stats.Monodirectional((1.0, 0.0, 0.0))

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 100
    model.settings.source = src

    return model, xmid, cell1, cell2


def test_surface_filter_flux_normal_incidence(two_cell_model, run_in_tmpdir):
    """SurfaceFilter + flux at mu=1 gives w/|mu| = 1.0 per source particle."""
    model, xmid, *_ = two_cell_model

    surf_filter = openmc.SurfaceFilter([xmid])
    flux_tally = openmc.Tally()
    flux_tally.filters = [surf_filter]
    flux_tally.scores = ['flux']
    model.tallies = [flux_tally]

    model.run(apply_tally_results=True)
    flux_mean = flux_tally.mean.flat[0]

    # Every particle crosses at mu=1 with weight 1, so flux = 1.0
    assert flux_mean == pytest.approx(1.0, rel=1e-8)


def test_surface_filter_current_outward(two_cell_model, run_in_tmpdir):
    """SurfaceFilter + current gives +1.0 for purely outward crossings."""
    model, xmid, *_ = two_cell_model

    surf_filter = openmc.SurfaceFilter([xmid])
    current_tally = openmc.Tally()
    current_tally.filters = [surf_filter]
    current_tally.scores = ['current']
    model.tallies = [current_tally]

    model.run(apply_tally_results=True)
    current_mean = current_tally.mean.flat[0]

    # All crossings are outward → net current = +1.0
    assert current_mean == pytest.approx(1.0)


def test_surface_filter_flux_angled(two_cell_model, run_in_tmpdir):
    """Surface flux at 60-degree incidence gives w/|mu| = 2.0."""
    model, xmid, *_ = two_cell_model

    # Modify source to use 60-degree angle from normal: mu = cos(60°) = 0.5
    mu = ux = 0.5
    uy = math.sqrt(1.0 - ux**2)
    model.settings.source[0].angle = openmc.stats.Monodirectional((ux, uy, 0.0))

    surf_filter = openmc.SurfaceFilter([xmid])
    flux_tally = openmc.Tally()
    flux_tally.filters = [surf_filter]
    flux_tally.scores = ['flux']
    model.tallies = [flux_tally]

    model.run(apply_tally_results=True)
    flux_mean = flux_tally.mean.flat[0]

    # flux = w/|mu| = 1/0.5 = 2.0
    assert flux_mean == pytest.approx(1.0 / mu)


def test_cellfrom_filter_flux_directional(two_cell_model, run_in_tmpdir):
    """SurfaceFilter + CellFromFilter + flux scores only the correct direction."""
    model, xmid, cell1, cell2 = two_cell_model

    surf_filter = openmc.SurfaceFilter([xmid])
    cellfrom_filter = openmc.CellFromFilter([cell1, cell2])

    tally = openmc.Tally()
    tally.filters = [surf_filter, cellfrom_filter]
    tally.scores = ['flux']

    model.tallies = [tally]
    model.run(apply_tally_results=True)
    mean_from1 = tally.mean.flat[0]
    mean_from2 = tally.mean.flat[1]

    # All particles cross xmid from cell1 at mu=1 → flux = 1.0
    assert mean_from1 == pytest.approx(1.0)
    # No particles cross xmid from cell2 → flux = 0
    assert mean_from2 == pytest.approx(0.0)
