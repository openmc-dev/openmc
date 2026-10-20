from contextlib import nullcontext

import numpy as np
import pytest
import openmc


@pytest.mark.parametrize('samples', [216, (6, 6, 6)])
@pytest.mark.parametrize('defect', ['none', 'overlap', 'undefined'])
def test_geometry_debug(samples, defect, run_in_tmpdir, capfd):
    """Exercise actual slice generation with a single-voxel defect."""
    outer = -openmc.model.RectangularParallelepiped(0, 6, 0, 6, 0, 6)
    inner = -openmc.model.RectangularParallelepiped(1, 2, 3, 4, 2, 3)
    cell = openmc.Cell(region=outer & ~inner if defect == 'undefined' else outer)
    cells = [cell]
    if defect == 'overlap':
        cells.append(openmc.Cell(region=inner))
    model = openmc.Model(geometry=openmc.Geometry(cells))

    context = (pytest.warns(UserWarning, match='Sampling resolution')
               if defect == 'undefined' else nullcontext())
    with context:
        result = model.geometry_debug((0, 0, 0), (6, 6, 6), samples,
                                      print_summary=True)
    output = capfd.readouterr().out
    assert bool(result['overlap_boxes']) == (defect == 'overlap')
    assert bool(result['undefined_boxes']) == (defect == 'undefined')
    assert result['under_resolved'] == (defect == 'undefined')
    if defect == 'none':
        assert 'Overlap bounding boxes: None' in output
        assert 'Undefined bounding boxes: None' in output
    else:
        boxes = result[f'{defect}_boxes']
        assert len(boxes) == 1
        np.testing.assert_allclose(boxes[0]['bbox'].lower_left, (1, 3, 2))
        np.testing.assert_allclose(boxes[0]['bbox'].upper_right, (2, 4, 3))
        if defect == 'overlap':
            assert boxes[0]['key'] == (
                model.geometry.root_universe.id, *sorted(c.id for c in cells))
            assert 'Overlaps found: 1' in output
        else:
            assert 'Undefined regions found: 1' in output


def test_undefined_regions(run_in_tmpdir):
    """Separate holes get separate boxes and resolution flags."""
    outer = -openmc.model.RectangularParallelepiped(0, 10, 0, 10, 0, 10)
    thick = -openmc.model.RectangularParallelepiped(1, 4, 1, 4, 1, 4)
    thin = -openmc.model.RectangularParallelepiped(6, 7, 6, 7, 6, 7)
    cell = openmc.Cell(region=outer & ~thick & ~thin)
    model = openmc.Model(geometry=openmc.Geometry([cell]))
    with pytest.warns(UserWarning, match='1 undefined region'):
        result = model.geometry_debug((0, 0, 0), (10, 10, 10), 1000)
    boxes = sorted(result['undefined_boxes'], key=lambda b: b['bbox'].lower_left[0])
    assert len(boxes) == 2
    np.testing.assert_allclose(boxes[0]['bbox'].lower_left, (1, 1, 1))
    np.testing.assert_allclose(boxes[0]['bbox'].upper_right, (4, 4, 4))
    assert not boxes[0]['under_resolved']
    np.testing.assert_allclose(boxes[1]['bbox'].lower_left, (6, 6, 6))
    np.testing.assert_allclose(boxes[1]['bbox'].upper_right, (7, 7, 7))
    assert boxes[1]['under_resolved']
