import numpy as np
import pytest

import openmc


@pytest.fixture(params=[
    (openmc.XPlane, 'x0', (1., 0., 0.)),
    (openmc.YPlane, 'y0', (0., 1., 0.)),
    (openmc.ZPlane, 'z0', (0., 0., 1.)),
    (openmc.XCylinder, 'r', (0., 1., 0.)),
    (openmc.YCylinder, 'r', (1., 0., 0.)),
    (openmc.ZCylinder, 'r', (1., 0., 0.)),
    (openmc.Sphere, 'r', (1., 0., 0.)),
])
def ordered_surfaces(request):
    surface_type, coordinate, direction = request.param
    surfaces = [surface_type(**{coordinate: value}) for value in (1., 2., 4.)]
    return surfaces, np.asarray(direction)


def test_subdivide_defaults(ordered_surfaces):
    surfaces, direction = ordered_surfaces
    regions = openmc.model.subdivide(surfaces)
    expected = [-surfaces[0], +surfaces[0] & -surfaces[1],
                +surfaces[1] & -surfaces[2], +surfaces[2]]
    assert [str(region) for region in regions] == [str(r) for r in expected]


@pytest.mark.parametrize('include_first', [False, True])
@pytest.mark.parametrize('include_last', [False, True])
@pytest.mark.parametrize('as_tuple', [False, True])
def test_subdivide_end_regions(ordered_surfaces, include_first, include_last,
                               as_tuple):
    surfaces, direction = ordered_surfaces
    if as_tuple:
        surfaces = tuple(surfaces)
    regions = openmc.model.subdivide(
        surfaces, include_first=include_first, include_last=include_last)
    assert len(regions) == 2 + include_first + include_last

    # Each remaining region contains its sample point and no others.
    points = [distance * direction for distance in (0., 1.5, 3., 5.)]
    indices = [1, 2]
    if include_first:
        indices.insert(0, 0)
    if include_last:
        indices.append(3)
    for region, index in zip(regions, indices):
        assert [point in region for point in points] == [
            i == index for i in range(len(points))]

    # Subdivision reuses the input surfaces without creating or modifying them.
    used_surfaces = {}
    for region in regions:
        used_surfaces.update(region.get_surfaces())
    assert all(surface is used_surfaces[surface.id] for surface in surfaces)


@pytest.mark.parametrize('include_first', [False, True])
@pytest.mark.parametrize('include_last', [False, True])
def test_subdivide_single_surface(ordered_surfaces, include_first,
                                  include_last):
    surfaces, direction = ordered_surfaces
    surface = surfaces[0]
    regions = openmc.model.subdivide(
        [surface], include_first=include_first, include_last=include_last)
    expected = []
    if include_first:
        expected.append(-surface)
    if include_last:
        expected.append(+surface)
    assert [str(region) for region in regions] == [str(r) for r in expected]


@pytest.mark.parametrize('option', ['include_first', 'include_last'])
@pytest.mark.parametrize('value', [None, 1, 'False'])
def test_subdivide_invalid_option(option, value):
    with pytest.raises(TypeError, match=option):
        openmc.model.subdivide([openmc.ZPlane()], **{option: value})
