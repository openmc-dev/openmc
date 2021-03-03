import numpy as np
import openmc.lib
import pytest

@pytest.fixture(autouse=True)
def complex_cell(run_in_tmpdir, mpi_intracomm):

    openmc.reset_auto_ids()

    model = openmc.model.Model()

    u235 = openmc.Material()
    u235.set_density('g/cc', 4.5)
    u235.add_nuclide("U235", 1.0)

    u238 = openmc.Material()
    u238.set_density('g/cc', 4.5)
    u238.add_nuclide("U238", 1.0)

    zr90 = openmc.Material()
    zr90.set_density('g/cc', 2.0)
    zr90.add_nuclide("Zr90", 1.0)

    n14 = openmc.Material()
    n14.set_density('g/cc', 0.1)
    n14.add_nuclide("N14", 1.0)

    model.materials = (u235, u238, zr90, n14)

    s1 = openmc.XPlane(-10.0, boundary_type='vacuum')
    s2 = openmc.XPlane(-7.0)
    s3 = openmc.XPlane(-4.0)
    s4 = openmc.XPlane(4.0)
    s5 = openmc.XPlane(7.0)
    s6 = openmc.XPlane(10.0, boundary_type='vacuum')
    s7 = openmc.XPlane(0.0)

    s11 = openmc.YPlane(-10.0, boundary_type='vacuum')
    s12 = openmc.YPlane(-7.0)
    s13 = openmc.YPlane(-4.0)
    s14 = openmc.YPlane(4.0)
    s15 = openmc.YPlane(7.0)
    s16 = openmc.YPlane(10.0, boundary_type='vacuum')
    s17 = openmc.YPlane(0.0)

    c1 = openmc.Cell(fill=u235)
    c1.region = ~(-s3 | +s4 | ~(+s13 & -s14))

    c2 = openmc.Cell(fill=u238)
    c2.region =  ~(+s3 & -s4 & +s13 & -s14) & +s2 & -s5 & +s12 & -s15

    c3 = openmc.Cell(fill=zr90)
    c3.region = ((+s1 & -s7 & +s17 & -s16) | (+s7 & -s6 & +s11 & -s17)) \
                & (-s2 | +s5 | -s12 | +s15)

    c4 = openmc.Cell(fill=n14)
    c4.region = ((+s1 & -s7 & +s11 & -s17) | (+s7 & -s6 & +s17 & -s16)) & \
                ~(+s2 & -s5 & +s12 & -s15)

    c5 = openmc.Cell(fill=n14)
    c5.region = ~(+s1 & -s6 & +s11 & -s16)

    model.geometry.root_universe = openmc.Universe()
    model.geometry.root_universe.add_cells([c1, c2, c3, c4, c5])

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 100
    model.settings.source = openmc.Source(space=openmc.stats.Box(
        [-10., -10., -1.], [10., 10., 1.]))

    model.settings.verbosity = 1

    model.export_to_xml()

    openmc.lib.finalize()
    openmc.lib.init(intracomm=mpi_intracomm)

    yield

    openmc.lib.finalize()


expected_results = ( (1, (( -4.,  -4., -np.inf),
                          (  4.,   4.,  np.inf))),
                     (2, (( -7.,  -7., -np.inf),
                          (  7.,   7.,  np.inf))),
                     (3, ((-10., -10., -np.inf),
                          ( 10.,  10.,  np.inf))),
                     (4, ((-10., -10., -np.inf),
                          ( 10.,  10.,  np.inf))),
                     (5, ((-np.inf, -np.inf, -np.inf),
                          ( np.inf,  np.inf,  np.inf))) )
@pytest.mark.parametrize("cell_id,expected_box", expected_results)
def test_cell_box(cell_id, expected_box):
    cell_box = openmc.lib.cells[cell_id].bounding_box
    assert tuple(cell_box[0]) == expected_box[0]
    assert tuple(cell_box[1]) == expected_box[1]
