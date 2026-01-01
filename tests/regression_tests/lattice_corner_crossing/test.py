import openmc
import pytest

import numpy as np
from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()
    
    
    # length of lattice on each side
    lat_size = 20.0

    # angle that we're crossing the corner at
    phi = np.pi/4.0

    # number of lattice elements across
    numlat = 2

    air = openmc.Material(name='air')
    air.set_density('g/cc', 0.001)
    air.add_nuclide('N14', 1.)

    metal = openmc.Material(name='metal')
    metal.set_density('g/cc', 7)
    metal.add_nuclide('Fe56', 1.)

    model.materials.extend([air, metal])


    left = openmc.XPlane(x0=-lat_size/2.0, name='left')
    right = openmc.XPlane(x0=lat_size/2.0, name='right')
    bottom = openmc.YPlane(y0=-lat_size/2.0, name='bottom')
    top = openmc.YPlane(y0=lat_size/2.0, name='top')
    cyl = openmc.ZCylinder(x0=0, y0=0, r=lat_size)
    cyl.boundary_type = 'vacuum'

    in_lattice_region = +left & -right & +bottom & -top

    outside_lattice = openmc.Cell(name='outside_lattice')
    outside_lattice.region = -cyl & ~in_lattice_region
    outside_lattice.fill = air

    inside_lattice = openmc.Cell(name='inside_lattice')
    inside_lattice.region = in_lattice_region

    root = openmc.Universe(name='root universe')

    # can we straightforwardly define infinite universes?
    metal_uni = openmc.Universe()
    metal_cell = openmc.Cell()
    metal_cell.fill = metal
    metal_cell.region = -openmc.Sphere(r=1e90)
    metal_uni.add_cell(metal_cell)

    air_uni   = openmc.Universe()
    air_cell = openmc.Cell()
    air_cell.fill = air
    air_cell.region = -openmc.Sphere(r=1e90)
    air_uni.add_cell(air_cell)

    # define a checkerboard lattice
    lattice = openmc.RectLattice()
    lattice.lower_left = [-lat_size/2.0, -lat_size/2.0]
    lattice.pitch = [lat_size/numlat, lat_size/numlat]
    lattice.universes = [[metal_uni if (i%2)==(j%2) else air_uni for i in range(numlat)]
                            for j in range(numlat)]
    inside_lattice.fill = lattice
    root.add_cells([outside_lattice, inside_lattice])
    
    model.geometry = openmc.Geometry(root)


    # Instantiate a Settings object, set all runtime parameters
    model.settings = openmc.Settings() 
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 10
    model.settings.particles = int(1e3)
    model.settings.max_lost_particles = int(1e8)

    # define a source located outside the lattice and pointing straight into its
    # corner at 45 degrees
    angular_distr = openmc.stats.Monodirectional(reference_uvw=[np.cos(phi), np.sin(phi), 0.0])
    offset = 0.0 #1e-6
    spatial_distr = openmc.stats.Point(xyz=[-np.cos(phi)+offset, -np.sin(phi), 0.0])
    model.settings.source = openmc.source.Source(space=spatial_distr, angle=angular_distr)

    # Instantiate a tally mesh
    mesh = openmc.RegularMesh()
    mesh.dimension = [10, 10, 1]
    mesh.lower_left = [-lat_size, -lat_size, -1e9]
    mesh.upper_right = [lat_size, lat_size, 1e9]
    mesh_filter = openmc.MeshFilter(mesh)

    # Instantiate the Tally
    tally = openmc.Tally(tally_id=1)
    tally.filters = [mesh_filter]
    tally.scores = ['flux']
    tally.estimator = 'tracklength'
    model.tallies = openmc.Tallies([tally])
    return model
    
    

def test_lattice_corner_crossing(model):
    # Ensure we account for potential corner crossings
    # in floating point precision.
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
