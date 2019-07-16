from tests.testing_harness import TestHarness

import sys

import openmc.capi

def test_complex_cell():
    harness = TestHarness('statepoint.10.h5')
    harness.main()

def test_complex_cell_capi():
    # initialize
    openmc.capi.init([])

    inf = sys.float_info.max

    expected_boxes = { 1 : (( -4.,  -4., -inf), ( 4.,  4., inf)),
                       2 : (( -7.,  -7., -inf), ( 7.,  7., inf)),
                       3 : ((-10., -10., -inf), (10., 10., inf)),
                       4 : ((-10., -10., -inf), (10., 10., inf)) }

    for cell_id, cell in openmc.capi.cells.items():
        cell_box = cell.bounding_box

        assert tuple(cell_box[0]) == expected_boxes[cell_id][0]
        assert tuple(cell_box[1]) == expected_boxes[cell_id][1]

        cell_box = openmc.capi.bounding_box("Cell", cell_id)

        assert tuple(cell_box[0]) == expected_boxes[cell_id][0]
        assert tuple(cell_box[1]) == expected_boxes[cell_id][1]

    openmc.capi.finalize()
