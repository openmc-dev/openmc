import numpy as np
import pytest
import openmc
import openmc.lib


@pytest.fixture(params=['spheres and rods', 'spheres'])
def many_objects_model(request):
    """Objects in a box, with the background defined as the space outside of
    all of them. The background region, and the region of each object outside
    of the preceding ones, are searched with a tree over the boxes of the
    objects. With only spheres, the background is an intersection of
    half-spaces."""
    rng = np.random.default_rng(1)
    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    mat.set_density('g/cm3', 0.1)

    box = openmc.model.RectangularParallelepiped(
        -10, 10, -10, 10, -10, 10, boundary_type='vacuum')
    objects = []
    # Spheres, as for TRISO particles
    n_spheres = 40 if request.param == 'spheres and rods' else 60
    for _ in range(n_spheres):
        center = rng.uniform(-8, 8, 3)
        objects.append(-openmc.Sphere(*center, r=rng.uniform(0.3, 1.0)))
    # Finite rods, whose complements are unions
    n_rods = 20 if request.param == 'spheres and rods' else 0
    for _ in range(n_rods):
        x0, y0 = rng.uniform(-8, 8, 2)
        z0 = rng.uniform(-8, 6)
        objects.append(-openmc.ZCylinder(x0, y0, r=rng.uniform(0.2, 0.8)) &
                       +openmc.ZPlane(z0) & -openmc.ZPlane(z0 + 2.0))

    # Objects may overlap, so each cell is the part of an object outside of
    # the preceding objects
    cells = []
    for i, obj in enumerate(objects):
        region = openmc.Intersection([obj])
        for other in objects[:i]:
            region &= ~other
        cells.append(openmc.Cell(fill=mat, region=region))
    background = openmc.Intersection([-box])
    for obj in objects:
        background &= ~obj
    cells.append(openmc.Cell(fill=mat, region=background))

    model = openmc.Model()
    model.geometry = openmc.Geometry(cells)
    model.materials = openmc.Materials([mat])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1000
    model.settings.batches = 2
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-9, -9, -9), (9, 9, 9)))
    return model


def test_find_cell_outside_many_objects(run_in_tmpdir, many_objects_model):
    model = many_objects_model
    model.export_to_model_xml()
    cells = model.geometry.root_universe.cells
    points = np.random.default_rng(2).uniform(-9.9, 9.9, size=(2000, 3))
    openmc.lib.init()
    try:
        for p in points:
            cell, _ = openmc.lib.find_cell(p)
            expected = [c.id for c in cells.values() if tuple(p) in c.region]
            assert [cell.id] == expected
    finally:
        openmc.lib.finalize()


def test_transport_outside_many_objects(run_in_tmpdir, many_objects_model):
    model = many_objects_model
    model.settings.max_lost_particles = 1
    model.run()
