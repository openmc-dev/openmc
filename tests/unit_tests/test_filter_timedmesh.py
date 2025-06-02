import math

import numpy as np
import pytest
from uncertainties import unumpy

import openmc


def test_filter_timed_Mesh(run_in_tmpdir):
    """
    Test that TimeFilter+MeshFilter with collision estimator agree
    with TimedMeshFilter with tracklength estimator.
    """

    # ===========================================================================
    # Set Material
    # ===========================================================================

    mat = openmc.Material()
    mat.add_nuclide('Fe56', 1.0)
    mat.set_density('g/cm3', 7.8)

    # ===========================================================================
    # Set geometry
    # ===========================================================================

    # Instantiate ZCylinder surfaces
    surf_Z1 = openmc.XPlane(surface_id=1, x0=-1e10, boundary_type="reflective")
    surf_Z2 = openmc.XPlane(surface_id=2, x0=1e10, boundary_type="reflective")

    # Instantiate Cells
    cell_F = openmc.Cell(cell_id=1, name="F")

    # Use surface half-spaces to define regions
    cell_F.region = +surf_Z1 & -surf_Z2

    # Register Materials with Cells
    cell_F.fill = mat

    # Instantiate Universes
    root = openmc.Universe(universe_id=0, name="root universe", cells=[cell_F])

    # Instantiate a Geometry, register the root Universe, and export to XML
    geometry = openmc.Geometry(root)

    # ===========================================================================
    # Settings
    # ===========================================================================

    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings_file = openmc.Settings()
    settings_file.run_mode = "fixed source"
    settings_file.particles = 100000
    settings_file.batches = 10
    settings_file.output = {"tallies": False}
    settings_file.cutoff = {"time_neutron": 1E-7}

    # Create an initial uniform spatial source distribution over fissionable zones
    delta_dist = openmc.stats.Point()
    isotropic = openmc.stats.Isotropic()
    settings_file.source = openmc.IndependentSource(space=delta_dist, angle=isotropic)

    # ===========================================================================
    # Set tallies
    # ===========================================================================

    # Create a mesh filter that can be used in a tally
    mesh = openmc.RegularMesh()
    mesh.dimension = (21, 1, 1)
    mesh.lower_left = (-20.5, -1e10, -1e10)
    mesh.upper_right = (20.5, 1e10, 1e10)
    time_grid = np.linspace(0.0, 1E-7, 21)

    mesh_filter = openmc.MeshFilter(mesh)
    time_filter = openmc.TimeFilter(time_grid)
    timed_mesh_filter = openmc.TimedMeshFilter(mesh, time_grid)

    # Now use the mesh filter in a tally and indicate what scores are desired
    tally1 = openmc.Tally(name="collision")
    tally1.estimator = "collision"
    tally1.filters = [time_filter, mesh_filter]
    tally1.scores = ["flux"]

    tally2 = openmc.Tally(name="tracklength")
    tally2.estimator = "tracklength"
    tally2.filters = [timed_mesh_filter]
    tally2.scores = ["flux"]

    # Instantiate a tallies collection
    tallies = openmc.Tallies([tally1, tally2])
    
    # ===========================================================================
    # Set the model
    # ===========================================================================

    model = openmc.Model()
    model.geometry = geometry
    model.settings = settings_file
    model.tallies = tallies

    # ===========================================================================
    # Run and post-process
    # ===========================================================================

    sp_filename = model.run()

    # Get radial flux distribution
    with openmc.StatePoint(sp_filename) as sp:
        flux_collision = sp.tallies[tally1.id].mean.ravel()
        flux_collision_unc = sp.tallies[tally1.id].std_dev.ravel()
        flux_tracklength = sp.tallies[tally2.id].mean.ravel()
        flux_tracklength_unc = sp.tallies[tally2.id].std_dev.ravel()

    # Construct arrays with uncertainties
    collision = unumpy.uarray(flux_collision, flux_collision_unc)
    tracklength = unumpy.uarray(flux_tracklength, flux_tracklength_unc)
    delta = collision - tracklength

    # Check that difference is within uncertainty
    diff = unumpy.nominal_values(delta)
    std_dev = unumpy.std_devs(delta)
    assert np.all(diff < 3*std_dev)
