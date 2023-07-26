from pathlib import Path

import openmc
import openmc.examples
import pytest
from matplotlib.colors import LogNorm


@pytest.fixture(scope='module')
def myplot():
    plot = openmc.Plot(name='myplot')
    plot.width = (100., 100.)
    plot.origin = (2., 3., -10.)
    plot.pixels = (500, 500)
    plot.filename = 'myplot'
    plot.type = 'slice'
    plot.basis = 'yz'
    plot.background = 'black'
    plot.background = (0, 0, 0)

    plot.color_by = 'material'
    m1, m2 = openmc.Material(), openmc.Material()
    plot.colors = {m1: (0, 255, 0), m2: (0, 0, 255)}
    plot.colors = {m1: 'green', m2: 'blue'}

    plot.mask_components = [openmc.Material()]
    plot.mask_background = 'white'
    plot.mask_background = (255, 255, 255)

    plot.overlap_color = (255, 211, 0)
    plot.overlap_color = 'yellow'
    plot.show_overlaps = True

    plot.level = 1
    plot.meshlines = {
        'type': 'tally',
        'id': 1,
        'linewidth': 2,
        'color': (40, 30, 20)
    }
    return plot


@pytest.fixture(scope='module')
def myprojectionplot():
    plot = openmc.ProjectionPlot(name='myprojectionplot')
    plot.look_at = (0.0, 0.0, 0.0)
    plot.camera_position = (4.0, 3.0, 0.0)
    plot.pixels = (500, 500)
    plot.filename = 'myprojectionplot'
    plot.background = (0, 0, 0)
    plot.background = 'black'

    plot.color_by = 'material'
    m1, m2 = openmc.Material(), openmc.Material()
    plot.colors = {m1: (0, 255, 0), m2: (0, 0, 255)}
    plot.colors = {m1: 'green', m2: 'blue'}
    plot.xs = {m1: 1.0, m2: 0.01}

    plot.mask_components = [openmc.Material()]
    plot.mask_background = (255, 255, 255)
    plot.mask_background = 'white'

    plot.overlap_color = (255, 211, 0)
    plot.overlap_color = 'yellow'

    plot.wireframe_thickness = 2

    plot.level = 1
    return plot


def test_voxel_plot(run_in_tmpdir):
    # attempt to preload VTK and skip this test if unavailable
    vtk = pytest.importorskip('vtk')
    surf1 = openmc.Sphere(r=500, boundary_type='vacuum')
    cell1 = openmc.Cell(region=-surf1)
    geometry = openmc.Geometry([cell1])
    geometry.export_to_xml()
    materials = openmc.Materials()
    materials.export_to_xml()
    vox_plot = openmc.Plot()
    vox_plot.type = 'voxel'
    vox_plot.id = 12
    vox_plot.width = (1500., 1500., 1500.)
    vox_plot.pixels = (200, 200, 200)
    vox_plot.color_by = 'cell'
    vox_plot.to_vtk('test_voxel_plot.vti')

    assert Path('plot_12.h5').is_file()
    assert Path('test_voxel_plot.vti').is_file()

    vox_plot.filename = 'h5_voxel_plot'
    vox_plot.to_vtk('another_test_voxel_plot.vti')

    assert Path('h5_voxel_plot.h5').is_file()
    assert Path('another_test_voxel_plot.vti').is_file()

    slice_plot = openmc.Plot()
    with pytest.raises(ValueError):
        slice_plot.to_vtk('shimmy.vti')


def test_attributes(myplot):
    assert myplot.name == 'myplot'


def test_attributes_proj(myprojectionplot):
    assert myprojectionplot.name == 'myprojectionplot'


def test_repr(myplot):
    r = repr(myplot)
    assert isinstance(r, str)


def test_repr_proj(myprojectionplot):
    r = repr(myprojectionplot)
    assert isinstance(r, str)


def test_from_geometry():
    width = 25.
    s = openmc.Sphere(r=width/2, boundary_type='vacuum')
    c = openmc.Cell(region=-s)
    univ = openmc.Universe(cells=[c])
    geom = openmc.Geometry(univ)

    for basis in ('xy', 'yz', 'xz'):
        plot = openmc.Plot.from_geometry(geom, basis)
        assert plot.origin == pytest.approx((0., 0., 0.))
        assert plot.width == pytest.approx((width, width))
        assert plot.basis == basis


def test_highlight_domains():
    plot = openmc.Plot()
    plot.color_by = 'material'
    plots = openmc.Plots([plot])

    model = openmc.examples.pwr_pin_cell()
    mats = {m for m in model.materials if 'UO2' in m.name}
    plots.highlight_domains(model.geometry, mats)


def test_xml_element(myplot):
    elem = myplot.to_xml_element()
    assert 'id' in elem.attrib
    assert 'color_by' in elem.attrib
    assert 'type' in elem.attrib
    assert elem.find('origin') is not None
    assert elem.find('width') is not None
    assert elem.find('pixels') is not None
    assert elem.find('background').text == '0 0 0'

    newplot = openmc.Plot.from_xml_element(elem)
    attributes = ('id', 'color_by', 'filename', 'type', 'basis', 'level',
                  'meshlines', 'show_overlaps', 'origin', 'width', 'pixels',
                  'background', 'mask_background')
    for attr in attributes:
        assert getattr(newplot, attr) == getattr(myplot, attr), attr


def test_to_xml_element_proj(myprojectionplot):
    elem = myprojectionplot.to_xml_element()
    assert 'id' in elem.attrib
    assert 'color_by' in elem.attrib
    assert 'type' in elem.attrib
    assert elem.find('camera_position') is not None
    assert elem.find('wireframe_thickness') is not None
    assert elem.find('look_at') is not None
    assert elem.find('pixels') is not None
    assert elem.find('background').text == '0 0 0'


def test_plots(run_in_tmpdir):
    p1 = openmc.Plot(name='plot1')
    p1.origin = (5., 5., 5.)
    p1.colors = {10: (255, 100, 0)}
    p1.mask_components = [2, 4, 6]
    p2 = openmc.Plot(name='plot2')
    p2.origin = (-3., -3., -3.)
    plots = openmc.Plots([p1, p2])
    assert len(plots) == 2

    p3 = openmc.ProjectionPlot(name='plot3')
    plots = openmc.Plots([p1, p2, p3])
    assert len(plots) == 3

    p4 = openmc.Plot(name='plot4')
    plots.append(p4)
    assert len(plots) == 4

    plots.export_to_xml()

    # from_xml
    new_plots = openmc.Plots.from_xml()
    assert len(new_plots)
    assert new_plots[0].origin == p1.origin
    assert new_plots[0].colors == p1.colors
    assert new_plots[0].mask_components == p1.mask_components
    assert new_plots[1].origin == p2.origin


def test_voxel_plot_roundtrip():
    # Define a voxel plot and create XML element
    plot = openmc.Plot(name='my voxel plot')
    plot.type = 'voxel'
    plot.filename = 'voxel1'
    plot.pixels = (50, 50, 50)
    plot.origin = (0., 0., 0.)
    plot.width = (75., 75., 75.)
    plot.color_by = 'material'
    elem = plot.to_xml_element()

    # Read back from XML and make sure it hasn't changed
    new_plot = plot.from_xml_element(elem)
    assert new_plot.name == plot.name
    assert new_plot.filename == plot.filename
    assert new_plot.type == plot.type
    assert new_plot.pixels == plot.pixels
    assert new_plot.origin == plot.origin
    assert new_plot.width == plot.width
    assert new_plot.color_by == plot.color_by


# def test_plot_mesh_tally(run_in_tmpdir):

mat1 = openmc.Material()
mat1.add_nuclide('Li6', 1, percent_type='ao')
mats = openmc.Materials([mat1])

# this shape is chose to create axis with testable ranges
surface = openmc.model.RectangularParallelepiped(
    -100, 50, -200, 250, -300, 350, boundary_type='vacuum'
)
cell = openmc.Cell(region=-surface)
cell.fill = mat1
geom = openmc.Geometry([cell])

source = openmc.IndependentSource()
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14e6], [1])
# puts the source in the center of the RectangularParallelepiped
source.space = openmc.stats.Point(cell.bounding_box.center)

sett = openmc.Settings()
sett.batches = 2
sett.inactive = 0
sett.particles = 5000
sett.run_mode = 'fixed source'
sett.source = source

mesh = openmc.RegularMesh().from_domain(geom, dimension=[10, 20, 30])
mesh_filter = openmc.MeshFilter(mesh)
mesh_tally = openmc.Tally(name='mesh-tal')
mesh_tally.filters = [mesh_filter]
mesh_tally.scores = ['flux']
tallies = openmc.Tallies([mesh_tally])

model = openmc.Model(geom, mats, sett, tallies)
sp_filename = model.run()
statepoint = openmc.StatePoint(sp_filename)
tally_result = statepoint.get_tally(name='mesh-tal')

plot = openmc.plot_mesh_tally(
    tally=tally_result,
    basis='xy',
    slice_index=29  # max value of slice selected
)
# axis_units defaults to cm
assert plot.xaxis.get_label().get_text() == 'x [cm]'
assert plot.yaxis.get_label().get_text() == 'y [cm]'
assert plot.get_xlim() == (-100., 50)
assert plot.get_ylim() == (-200., 250.)

plot = openmc.plot_mesh_tally(
    tally=tally_result,
    basis='yz',
    axis_units='m',
    slice_index=9,  # max value of slice selected
    value= 'std_dev'
)
assert plot.xaxis.get_label().get_text() == 'y [m]'
assert plot.yaxis.get_label().get_text() == 'z [m]'
assert plot.get_xlim() == (-2., 2.5)  # note that units are in m
assert plot.get_ylim() == (-3., 3.5)

plot = openmc.plot_mesh_tally(
    tally=tally_result,
    basis='xz',
    slice_index=19,  # max value of slice selected
    axis_units='mm',
    score='flux',
    value= 'mean',
    norm=LogNorm(vmin=0.1, vmax=100)
)
plot.colorbar()
assert plot.xaxis.get_label().get_text() == 'x [mm]'
assert plot.yaxis.get_label().get_text() == 'z [mm]'
assert plot.get_xlim() == (-1000., 500)  # note that units are in mm
assert plot.get_ylim() == (-3000.0, 3500.0)
plot.figure.savefig('z.png')