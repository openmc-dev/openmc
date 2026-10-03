import numpy as np
import pytest
import openmc
import openmc.lib


@pytest.fixture
def many_cells_model():
    """Spheres, finite rods and their overlaps in a box, with a background
    cell outside of all of them, so that most cells are skipped based on their
    bounding boxes when searching for the cell containing a point."""
    rng = np.random.default_rng(1)
    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    mat.set_density('g/cm3', 0.1)

    box = openmc.model.RectangularParallelepiped(
        -10, 10, -10, 10, -10, 10, boundary_type='vacuum')
    objects = []
    for _ in range(40):
        center = rng.uniform(-8, 8, 3)
        objects.append(-openmc.Sphere(*center, r=rng.uniform(0.3, 1.0)))
    for _ in range(20):
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


def test_find_cell_many_cells(run_in_tmpdir, many_cells_model):
    model = many_cells_model
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


def test_transport_many_cells(run_in_tmpdir, many_cells_model):
    model = many_cells_model
    model.settings.max_lost_particles = 1
    model.run()


def random_region(rng, surfaces, depth):
    """Random region expression with unions, intersections and complements"""
    if depth == 0 or rng.random() < 0.25:
        s = surfaces[rng.integers(len(surfaces))]
        return -s if rng.random() < 0.5 else +s
    n = rng.integers(2, 4)
    terms = [random_region(rng, surfaces, depth - 1) for _ in range(n)]
    region = openmc.Union(terms) if rng.random() < 0.5 else openmc.Intersection(terms)
    return ~region if rng.random() < 0.3 else region


@pytest.mark.parametrize('seed', range(5))
def test_find_cell_random_regions(run_in_tmpdir, seed):
    """Cells found for random points agree with Python's evaluation of random
    region expressions, and the points lie within the bounding boxes of the
    cells, which the search relies on."""
    rng = np.random.default_rng(seed)
    surfaces = [
        openmc.XPlane(-1.0), openmc.XPlane(2.0), openmc.YPlane(0.5),
        openmc.ZPlane(1.5), openmc.Sphere(x0=1.0, r=3.0),
        openmc.ZCylinder(y0=-1.0, r=2.0), openmc.XCylinder(z0=1.0, r=2.5),
        openmc.Plane(a=1.0, b=1.0, c=0.0, d=0.5),
    ]
    outer = openmc.Sphere(r=8.0, boundary_type='vacuum')

    # Split space into cells: each cell is a random region outside of the
    # preceding ones
    regions = [random_region(rng, surfaces, 3) for _ in range(6)]
    cells = []
    remaining = openmc.Intersection([-outer])
    for region in regions:
        cells.append(openmc.Cell(region=remaining & region))
        remaining = remaining & ~region
    cells.append(openmc.Cell(region=remaining))

    model = openmc.Model()
    model.geometry = openmc.Geometry(cells)
    model.settings.run_mode = 'fixed source'
    model.settings.particles = 1
    model.settings.batches = 1
    model.export_to_model_xml()

    points = rng.uniform(-5.6, 5.6, size=(500, 3))
    openmc.lib.init()
    try:
        for p in points:
            expected = [c.id for c in cells if tuple(p) in c.region]
            if len(expected) != 1:
                continue
            cell, _ = openmc.lib.find_cell(p)
            assert cell.id == expected[0]
            lower_left, upper_right = cell.bounding_box
            assert np.all(lower_left <= p) and np.all(p <= upper_right)
    finally:
        openmc.lib.finalize()
