from math import pi
import openmc
import pytest

from tests.regression_tests import config


@pytest.fixture(scope='module')
def pin_model_attributes():
    uo2 = openmc.Material(name='UO2')
    uo2.set_density('g/cm3', 10.29769)
    uo2.add_element('U', 1., enrichment=2.4)
    uo2.add_element('O', 2.)

    zirc = openmc.Material(name='Zirc')
    zirc.set_density('g/cm3', 6.55)
    zirc.add_element('Zr', 1.)

    borated_water = openmc.Material(name='Borated water')
    borated_water.set_density('g/cm3', 0.740582)
    borated_water.add_element('B', 4.0e-5)
    borated_water.add_element('H', 5.0e-2)
    borated_water.add_element('O', 2.4e-2)
    borated_water.add_s_alpha_beta('c_H_in_H2O')

    mats = openmc.Materials([uo2, zirc, borated_water])

    pitch = 1.25984
    fuel_or = openmc.ZCylinder(r=0.39218, name='Fuel OR')
    clad_or = openmc.ZCylinder(r=0.45720, name='Clad OR')
    box = openmc.model.rectangular_prism(pitch, pitch, boundary_type='reflective')

    # Define cells
    fuel = openmc.Cell(name='fuel', fill=uo2, region=-fuel_or)
    clad = openmc.Cell(fill=zirc, region=+fuel_or & -clad_or)
    water = openmc.Cell(fill=borated_water, region=+clad_or & box)

    # Define overall geometry
    geom = openmc.Geometry([fuel, clad, water])
    uo2.volume = pi * fuel_or.r**2

    settings = openmc.Settings()
    settings.batches = 100
    settings.inactive = 10
    settings.particles = 1000

    # Create an initial uniform spatial source distribution over fissionable zones
    bounds = [-0.62992, -0.62992, -1, 0.62992, 0.62992, 1]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    settings.source = openmc.source.Source(space=uniform_dist)

    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.lower_left = [-0.39218, -0.39218, -1.e50]
    entropy_mesh.upper_right = [0.39218, 0.39218, 1.e50]
    entropy_mesh.dimension = [10, 10, 1]
    settings.entropy_mesh = entropy_mesh

    tals = openmc.Tallies()
    tal = openmc.Tally(name='test')
    tal.scores = ['flux']
    tals.append(tal)

    plot = openmc.Plot()
    plot.origin = (0., 0., 0.)
    plot.width = (pitch, pitch)
    plot.pixels = (300, 300)
    plot.color_by = 'material'
    plots = openmc.Plots((plot,))

    chain = './chain_simple.xml'
    fission_q = {'U235': 200e6}

    chain_file_xml = """<?xml version="1.0"?>
<depletion_chain>
  <nuclide name="I135" decay_modes="1" reactions="1" half_life="2.36520E+04">
    <decay type="beta" target="Xe135" branching_ratio="1.0" />
    <reaction type="(n,gamma)" Q="0.0" target="Xe136" /> <!-- Not precisely true, but whatever -->
  </nuclide>
  <nuclide name="Xe135" decay_modes="1" reactions="1" half_life="3.29040E+04">
    <decay type=" beta" target="Cs135" branching_ratio="1.0" />
    <reaction type="(n,gamma)" Q="0.0" target="Xe136" />
  </nuclide>
  <nuclide name="Xe136" decay_modes="0" reactions="0" />
  <nuclide name="Cs135" decay_modes="0" reactions="0" />
  <nuclide name="Gd157" decay_modes="0" reactions="1"  >
    <reaction type="(n,gamma)" Q="0.0" target="Nothing" />
  </nuclide>
  <nuclide name="Gd156" decay_modes="0" reactions="1">
    <reaction type="(n,gamma)" Q="0.0" target="Gd157" />
  </nuclide>
  <nuclide name="U234" decay_modes="0" reactions="1">
    <reaction type="fission" Q="191840000."/>
    <neutron_fission_yields>
      <energies>2.53000e-02</energies>
      <fission_yields energy="2.53000e-02">
        <products>Gd157 Gd156 I135 Xe135 Xe136 Cs135</products>
        <data>1.093250e-04 2.087260e-04 2.780820e-02 6.759540e-03 2.392300e-02 4.356330e-05</data>
      </fission_yields>
    </neutron_fission_yields>
  </nuclide>
  <nuclide name="U235" decay_modes="0" reactions="1">
    <reaction type="fission" Q="193410000."/>
    <neutron_fission_yields>
      <energies>2.53000e-02</energies>
      <fission_yields energy="2.53000e-02">
        <products>Gd157 Gd156 I135 Xe135 Xe136 Cs135</products>
        <data>6.142710e-5 1.483250e-04 0.0292737 0.002566345 0.0219242 4.9097e-6</data>
      </fission_yields>
    </neutron_fission_yields>
  </nuclide>
  <nuclide name="U238" decay_modes="0" reactions="1">
    <reaction type="fission" Q="197790000."/>
    <neutron_fission_yields>
      <energies>2.53000e-02</energies>
      <fission_yields energy="2.53000e-02">
        <products>Gd157 Gd156 I135 Xe135 Xe136 Cs135</products>
        <data>4.141120e-04 7.605360e-04 0.0135457 0.00026864 0.0024432 3.7100E-07</data>
      </fission_yields>
    </neutron_fission_yields>
  </nuclide>
</depletion_chain>
"""

    return (mats, geom, settings, tals, plots, chain, fission_q,
            chain_file_xml)

@pytest.fixture(scope='module')
def mpi_intracomm():
    if config['mpi']:
        from mpi4py import MPI
        return MPI.COMM_WORLD
    else:
        return None


@pytest.fixture(scope='module')
def uo2():
    m = openmc.Material(material_id=100, name='UO2')
    m.add_nuclide('U235', 1.0)
    m.add_nuclide('O16', 2.0)
    m.set_density('g/cm3', 10.0)
    m.depletable = True
    return m


@pytest.fixture(scope='module')
def water():
    m = openmc.Material(name='light water')
    m.add_nuclide('H1', 2.0)
    m.add_nuclide('O16', 1.0)
    m.set_density('g/cm3', 1.0)
    m.add_s_alpha_beta('c_H_in_H2O')
    return m


@pytest.fixture(scope='module')
def sphere_model():
    model = openmc.model.Model()
    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    m.set_density('g/cm3', 1.0)
    model.materials.append(m)

    sph = openmc.Sphere(boundary_type='vacuum')
    c = openmc.Cell(fill=m, region=-sph)
    model.geometry.root_universe = openmc.Universe(cells=[c])

    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.Source(space=openmc.stats.Point())
    return model


@pytest.fixture
def cell_with_lattice():
    m_inside = [openmc.Material(), openmc.Material(), None, openmc.Material()]
    m_outside = openmc.Material()

    cyl = openmc.ZCylinder(r=1.0)
    inside_cyl = openmc.Cell(fill=m_inside, region=-cyl)
    outside_cyl = openmc.Cell(fill=m_outside, region=+cyl)
    univ = openmc.Universe(cells=[inside_cyl, outside_cyl])

    lattice = openmc.RectLattice(name='My Lattice')
    lattice.lower_left = (-4.0, -4.0)
    lattice.pitch = (4.0, 4.0)
    lattice.universes = [[univ, univ], [univ, univ]]
    main_cell = openmc.Cell(fill=lattice)

    return ([inside_cyl, outside_cyl, main_cell],
            [m_inside[0], m_inside[1], m_inside[3], m_outside],
            univ, lattice)

@pytest.fixture
def mixed_lattice_model(uo2, water):
    cyl = openmc.ZCylinder(r=0.4)
    c1 = openmc.Cell(fill=uo2, region=-cyl)
    c1.temperature = 600.0
    c2 = openmc.Cell(fill=water, region=+cyl)
    pin = openmc.Universe(cells=[c1, c2])

    empty = openmc.Cell()
    empty_univ = openmc.Universe(cells=[empty])

    hex_lattice = openmc.HexLattice()
    hex_lattice.center = (0.0, 0.0)
    hex_lattice.pitch = (1.2, 10.0)
    outer_ring = [pin]*6
    inner_ring = [empty_univ]
    axial_level = [outer_ring, inner_ring]
    hex_lattice.universes = [axial_level]*3
    hex_lattice.outer = empty_univ

    cell_hex = openmc.Cell(fill=hex_lattice)
    u = openmc.Universe(cells=[cell_hex])
    rotated_cell_hex = openmc.Cell(fill=u)
    rotated_cell_hex.rotation = (0., 0., 30.)
    ur = openmc.Universe(cells=[rotated_cell_hex])

    d = 6.0
    rect_lattice = openmc.RectLattice()
    rect_lattice.lower_left = (-d, -d)
    rect_lattice.pitch = (d, d)
    rect_lattice.outer = empty_univ
    rect_lattice.universes = [
        [ur, empty_univ],
        [empty_univ, u]
    ]

    xmin = openmc.XPlane(-d, boundary_type='periodic')
    xmax = openmc.XPlane(d, boundary_type='periodic')
    xmin.periodic_surface = xmax
    ymin = openmc.YPlane(-d, boundary_type='periodic')
    ymax = openmc.YPlane(d, boundary_type='periodic')
    main_cell = openmc.Cell(fill=rect_lattice,
                            region=+xmin & -xmax & +ymin & -ymax)

    # Create geometry and use unique material in each fuel cell
    geometry = openmc.Geometry([main_cell])
    geometry.determine_paths()
    c1.fill = [water.clone() for i in range(c1.num_instances)]

    return openmc.model.Model(geometry)
