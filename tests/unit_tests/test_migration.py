import openmc
import pytest
import numpy as np


def water():
    material = openmc.Material(name="Hydrogen")
    material.add_nuclide("H1", 1.0)
    material.set_density('g/cm3', 1.0)
    material.add_s_alpha_beta('c_H_in_H2O')
    return material


def migration_model(region, source_point=None, particles=10000, batches=10):
    """Fixed source model over `region` with a single migration-area tally."""
    model = openmc.Model()
    cell = openmc.Cell(region=region, fill=water())
    model.geometry = openmc.Geometry([cell])

    model.settings.particles = particles
    model.settings.batches = batches
    model.settings.run_mode = 'fixed source'
    if source_point is None:
        model.settings.source = openmc.IndependentSource()
    else:
        model.settings.source = openmc.IndependentSource(
            space=openmc.stats.Point(source_point))

    tally = openmc.Tally()
    tally.scores = ["migration-area"]
    model.tallies = [tally]

    return model, tally


def migration_area(model, tally):
    """Run `model` and return the tallied (mean, std_dev)."""
    model.run(apply_tally_results=True)
    return tally.mean.squeeze(), tally.std_dev.squeeze()


def assert_agrees(result1, result2):
    """Assert two migration areas agree to within three combined sigma."""
    mean1, std1 = result1
    mean2, std2 = result2
    std = (std1**2 + std2**2)**0.5
    assert np.abs(mean2 - mean1) <= 3*std


def sphere_model(radius, boundary_type):
    sphere = openmc.Sphere(r=radius, boundary_type=boundary_type)
    return migration_model(-sphere)


def test_reflective_is_equivalent_to_large_model(run_in_tmpdir):
    """A reflected sphere unfolds to the infinite medium a large sphere models."""
    openmc.reset_auto_ids()
    refl = migration_area(*sphere_model(2.5, "reflective"))
    large = migration_area(*sphere_model(100, "vacuum"))
    assert_agrees(refl, large)


def test_translational_periodic_is_equivalent_to_large_model(run_in_tmpdir):
    """A periodic cube tiles into the same infinite medium as a large sphere."""
    openmc.reset_auto_ids()

    d = 2.5
    xmin = openmc.XPlane(-d, boundary_type='periodic')
    xmax = openmc.XPlane(d, boundary_type='periodic')
    xmin.periodic_surface = xmax
    ymin = openmc.YPlane(-d, boundary_type='periodic')
    ymax = openmc.YPlane(d, boundary_type='periodic')
    ymin.periodic_surface = ymax
    zmin = openmc.ZPlane(-d, boundary_type='periodic')
    zmax = openmc.ZPlane(d, boundary_type='periodic')
    zmin.periodic_surface = zmax

    cube = +xmin & -xmax & +ymin & -ymax & +zmin & -zmax
    periodic = migration_area(*migration_model(cube))
    large = migration_area(*sphere_model(100, "vacuum"))
    assert_agrees(periodic, large)


def test_rotational_periodic_is_equivalent_to_full_cylinder(run_in_tmpdir):
    """A rotationally periodic wedge is the full cylinder folded about z."""
    openmc.reset_auto_ids()

    # 90 degree wedge about the z axis, so the wedge unfolds to the cylinder
    p1 = openmc.Plane.from_points(
        (0., 0., 0.), (1., 0., 0.), (0., 0., 1.), boundary_type='periodic')
    p2 = openmc.Plane.from_points(
        (0., 0., 0.), (0., 1., 0.), (0., 0., 1.), boundary_type='periodic')
    p1.periodic_surface = p2

    zcyl = openmc.ZCylinder(r=5., boundary_type='vacuum')
    zmin = openmc.ZPlane(-5., boundary_type='vacuum')
    zmax = openmc.ZPlane(5., boundary_type='vacuum')

    # Pick the sides of the planes that contain a point inside the wedge
    source_point = (1., 1., 0.)
    r1 = -p1 if source_point in -p1 else +p1
    r2 = -p2 if source_point in -p2 else +p2

    axial = +zmin & -zmax
    wedge = migration_area(*migration_model(
        r1 & r2 & -zcyl & axial, source_point=source_point))
    full = migration_area(*migration_model(
        -zcyl & axial, source_point=source_point))
    assert_agrees(wedge, full)


def test_white_boundary_rejected(run_in_tmpdir):
    """No isometry unfolds a white boundary, so it cannot be scored."""
    openmc.reset_auto_ids()
    model, _ = sphere_model(2.5, "white")
    with pytest.raises(RuntimeError, match='white boundary conditions'):
        model.run()


def test_albedo_rejected(run_in_tmpdir):
    """An albedo attenuates the longest displacements and biases the score."""
    openmc.reset_auto_ids()
    sphere = openmc.Sphere(r=2.5, boundary_type='reflective', albedo=0.5)
    model, _ = migration_model(-sphere)
    with pytest.raises(RuntimeError, match='albedo'):
        model.run()


@pytest.mark.parametrize("setup", [
    lambda t: setattr(t, 'estimator', 'analog'),
    lambda t: setattr(t, 'estimator', 'collision'),
    # another score in the same tally must not change the estimator undetected
    lambda t: setattr(t, 'scores', ['migration-area', 'nu-scatter']),
], ids=['analog', 'collision', 'analog-forced-by-other-score'])
def test_non_tracklength_estimator_rejected(setup, run_in_tmpdir):
    """The estimator is only final once every score has been processed."""
    openmc.reset_auto_ids()
    model, tally = sphere_model(2.5, 'vacuum')
    setup(tally)
    with pytest.raises(RuntimeError, match='tracklength estimator'):
        model.run()


def test_disallowed_filter_rejected(run_in_tmpdir):
    openmc.reset_auto_ids()
    model, tally = sphere_model(2.5, 'vacuum')
    tally.filters = [openmc.CellFilter(model.geometry.get_all_cells()[1])]
    with pytest.raises(RuntimeError, match='filters other than'):
        model.run()


def test_unit_albedo_accepted(run_in_tmpdir):
    """An albedo of one removes no weight and must not be rejected."""
    openmc.reset_auto_ids()
    sphere = openmc.Sphere(r=2.5, boundary_type='reflective', albedo=1.0)
    model, tally = migration_model(-sphere, particles=100, batches=2)
    migration_area(model, tally)


def test_meshborn_with_nonvacuum_boundary_rejected(run_in_tmpdir):
    """A MeshBorn filter cannot follow the birth point across a boundary."""
    openmc.reset_auto_ids()
    model, _ = sphere_model(2.5, "reflective")

    mesh = openmc.RegularMesh()
    mesh.lower_left = (-2.5, -2.5, -2.5)
    mesh.upper_right = (2.5, 2.5, 2.5)
    mesh.dimension = (2, 2, 2)

    flux_tally = openmc.Tally()
    flux_tally.filters = [openmc.MeshBornFilter(mesh)]
    flux_tally.scores = ['flux']
    model.tallies.append(flux_tally)

    with pytest.raises(RuntimeError, match='MeshBorn'):
        model.run()
