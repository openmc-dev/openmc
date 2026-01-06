from pathlib import Path

import openmc
import openmc.examples
import pytest

from openmc.plots import _SVG_COLORS


@pytest.fixture(scope='module')
def myplot():
    plot = openmc.SlicePlot(name='myplot')
    plot.width = (100., 100.)
    plot.origin = (2., 3., -10.)
    plot.pixels = (500, 500)
    plot.filename = './not-a-dir/myplot'
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
    plot = openmc.WireframeRayTracePlot(name='myprojectionplot')
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
    vox_plot = openmc.VoxelPlot()
    vox_plot.id = 12
    vox_plot.width = (1500., 1500., 1500.)
    vox_plot.pixels = (200, 200, 200)
    vox_plot.color_by = 'cell'
    vox_plot.to_vtk('test_voxel_plot.vti')

    assert Path('plot_12.h5').is_file()
    assert Path('test_voxel_plot.vti').is_file()

    vox_plot.filename = 'h5_voxel_plot'
    vox_plot.to_vtk(Path('another_test_voxel_plot.vti'))

    assert Path('h5_voxel_plot.h5').is_file()
    assert Path('another_test_voxel_plot.vti').is_file()

    # SlicePlot should not have to_vtk method
    slice_plot = openmc.SlicePlot()
    with pytest.raises(AttributeError):
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


def test_projection_plot_roundtrip(myprojectionplot):

    elem = myprojectionplot.to_xml_element()

    xml_plot = openmc.WireframeRayTracePlot.from_xml_element(elem)

    svg_colors = _SVG_COLORS

    assert xml_plot.name == myprojectionplot.name
    assert xml_plot.look_at == myprojectionplot.look_at
    assert xml_plot.camera_position == myprojectionplot.camera_position
    assert xml_plot.pixels == myprojectionplot.pixels
    assert xml_plot.filename == myprojectionplot.filename
    assert xml_plot.background == svg_colors[myprojectionplot.background]
    assert xml_plot.color_by == myprojectionplot.color_by
    expected_colors = {m.id: svg_colors[c] for m, c in myprojectionplot.colors.items()}
    assert xml_plot.colors == expected_colors
    # TODO: needs geometry information
    # assert xml_plot.mask_components == myprojectionplot.mask_components
    assert xml_plot.mask_background == svg_colors[myprojectionplot.mask_background]
    # assert xml_plot.overlap_color == svg_colors[myprojectionplot.overlap_color]
    assert xml_plot.wireframe_thickness == myprojectionplot.wireframe_thickness
    assert xml_plot.level == myprojectionplot.level


def test_from_geometry():
    width = 25.
    s = openmc.Sphere(r=width/2, boundary_type='vacuum')
    c = openmc.Cell(region=-s)
    univ = openmc.Universe(cells=[c])
    geom = openmc.Geometry(univ)

    for basis in ('xy', 'yz', 'xz'):
        plot = openmc.SlicePlot.from_geometry(geom, basis)
        assert plot.origin == pytest.approx((0., 0., 0.))
        assert plot.width == pytest.approx((width, width))
        assert plot.basis == basis


def test_highlight_domains():
    plot = openmc.SlicePlot()
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

    newplot = openmc.SlicePlot.from_xml_element(elem)
    attributes = ('id', 'color_by', 'filename', 'basis', 'level',
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
    p1 = openmc.SlicePlot(name='plot1')
    p1.origin = (5., 5., 5.)
    p1.colors = {10: (255, 100, 0)}
    p1.mask_components = [2, 4, 6]
    p2 = openmc.SlicePlot(name='plot2')
    p2.origin = (-3., -3., -3.)
    plots = openmc.Plots([p1, p2])
    assert len(plots) == 2

    p3 = openmc.WireframeRayTracePlot(name='plot3')
    plots = openmc.Plots([p1, p2, p3])
    assert len(plots) == 3

    p4 = openmc.VoxelPlot(name='plot4')
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
    plot = openmc.VoxelPlot(name='my voxel plot')
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
    assert new_plot.pixels == plot.pixels
    assert new_plot.origin == plot.origin
    assert new_plot.width == plot.width
    assert new_plot.color_by == plot.color_by


def test_phong_plot_roundtrip():
    plot = openmc.SolidRayTracePlot(name='my phong plot')
    plot.id = 2300
    plot.filename = 'phong1'
    plot.pixels = (50, 50)
    plot.look_at = (11., 12., 13.)
    plot.camera_position = (22., 23., 24.)
    plot.diffuse_fraction = 0.5
    plot.horizontal_field_of_view = 90.0
    plot.color_by = 'material'
    plot.light_position = (8., 9., 10.)
    plot.opaque_domains = [6, 7, 8]

    elem = plot.to_xml_element()

    repr(plot)

    new_plot = openmc.SolidRayTracePlot.from_xml_element(elem)

    assert new_plot.name == plot.name
    assert new_plot.id == plot.id
    assert new_plot.filename == plot.filename
    assert new_plot.pixels == plot.pixels
    assert new_plot.look_at == plot.look_at
    assert new_plot.camera_position == plot.camera_position
    assert new_plot.diffuse_fraction == plot.diffuse_fraction
    assert new_plot.horizontal_field_of_view == plot.horizontal_field_of_view
    assert new_plot.color_by == plot.color_by
    assert new_plot.light_position == plot.light_position
    assert new_plot.opaque_domains == plot.opaque_domains

    # ensure the new object is valid to re-write to XML
    new_elem = new_plot.to_xml_element()


def test_plot_directory(run_in_tmpdir):
    pwr_pin = openmc.examples.pwr_pin_cell()

    # create a standard slice plot, expected to work
    plot = openmc.SlicePlot()
    plot.filename = 'plot_1'
    plot.pixels = (10, 10)
    plot.color_by = 'material'
    plot.width = (100., 100.)
    pwr_pin.plots = [plot]
    pwr_pin.plot_geometry()

    # use current directory, also expected to work
    plot.filename = './plot_1'
    pwr_pin.plot_geometry()

    # use a non-existent directory, should raise an error
    plot.filename = './not-a-dir/plot_1'
    with pytest.raises(RuntimeError, match='does not exist'):
        pwr_pin.plot_geometry()
