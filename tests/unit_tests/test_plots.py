import openmc
import openmc.examples
import pytest


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
    assert len(plots)
    assert plots[0].origin == p1.origin
    assert plots[0].colors == p1.colors
    assert plots[0].mask_components == p1.mask_components
    assert plots[1].origin == p2.origin
