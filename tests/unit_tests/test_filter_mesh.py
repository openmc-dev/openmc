import math

import numpy as np
import pytest
from uncertainties import unumpy

import openmc


def test_spherical_mesh_estimators(run_in_tmpdir):
    """Test that collision/tracklength estimators agree for SphericalMesh"""

    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.set_density('g/cm3', 10.0)

    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 1_000
    model.settings.inactive = 10
    model.settings.batches = 20

    sph_mesh = openmc.SphericalMesh(
        r_grid=np.linspace(0.0, 5.0**3, 20)**(1/3)
    )
    tally1 = openmc.Tally()
    tally1.filters = [openmc.MeshFilter(sph_mesh)]
    tally1.scores = ['flux']
    tally1.estimator = 'collision'

    sph_mesh = openmc.SphericalMesh(
        r_grid=np.linspace(0.0, 5.0**3, 20)**(1/3)
    )
    tally2 = openmc.Tally()
    tally2.filters = [openmc.MeshFilter(sph_mesh)]
    tally2.scores = ['flux']
    tally2.estimator = 'tracklength'

    model.tallies = openmc.Tallies([tally1, tally2])

    # Run OpenMC
    sp_filename = model.run()

    # Get radial flux distribution
    with openmc.StatePoint(sp_filename) as sp:
        flux_collision = sp.tallies[tally1.id].mean.ravel()
        flux_collision_unc = sp.tallies[tally1.id].std_dev.ravel()
        flux_tracklength = sp.tallies[tally2.id].mean.ravel()
        flux_tracklength_unc = sp.tallies[tally2.id].std_dev.ravel()

    # Construct arrays with uncertainties
    collision = unumpy.uarray(flux_collision, flux_collision_unc)
    tracklength = unumpy.uarray(flux_tracklength, flux_tracklength_unc)
    delta = collision - tracklength

    # Check that difference is within uncertainty
    diff = unumpy.nominal_values(delta)
    std_dev = unumpy.std_devs(delta)
    assert np.all(diff < 3*std_dev)


def test_cylindrical_mesh_estimators(run_in_tmpdir):
    """Test that collision/tracklength estimators agree for CylindricalMesh"""

    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.set_density('g/cm3', 10.0)

    cyl = openmc.model.RightCircularCylinder((0., 0., -5.), 10., 10.0,
                                             boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-cyl)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 1_000
    model.settings.inactive = 10
    model.settings.batches = 20

    cyl_mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(0.0, 5.0**3, 20)**(1/3),
        z_grid=[-5., 5.]
    )
    tally1 = openmc.Tally()
    tally1.filters = [openmc.MeshFilter(cyl_mesh)]
    tally1.scores = ['flux']
    tally1.estimator = 'collision'

    cyl_mesh = openmc.CylindricalMesh(
        r_grid=np.linspace(0.0, 5.0**3, 20)**(1/3),
        z_grid=[-5., 5.]
    )
    tally2 = openmc.Tally()
    tally2.filters = [openmc.MeshFilter(cyl_mesh)]
    tally2.scores = ['flux']
    tally2.estimator = 'tracklength'

    model.tallies = openmc.Tallies([tally1, tally2])

    # Run OpenMC
    sp_filename = model.run()

    # Get radial flux distribution
    with openmc.StatePoint(sp_filename) as sp:
        flux_collision = sp.tallies[tally1.id].mean.ravel()
        flux_collision_unc = sp.tallies[tally1.id].std_dev.ravel()
        flux_tracklength = sp.tallies[tally2.id].mean.ravel()
        flux_tracklength_unc = sp.tallies[tally2.id].std_dev.ravel()

    # Construct arrays with uncertainties
    collision = unumpy.uarray(flux_collision, flux_collision_unc)
    tracklength = unumpy.uarray(flux_tracklength, flux_tracklength_unc)
    delta = collision - tracklength

    # Check that difference is within uncertainty
    diff = unumpy.nominal_values(delta)
    std_dev = unumpy.std_devs(delta)
    assert np.all(diff < 3*std_dev)


@pytest.mark.parametrize("scale", [0.1, 1.0, 1e2, 1e4, 1e5])
def test_cylindrical_mesh_coincident(scale, run_in_tmpdir):
    """Test for cylindrical mesh boundary being coincident with a cell boundary"""

    fuel = openmc.Material()
    fuel.add_nuclide('U235', 1.)
    fuel.set_density('g/cm3', 4.5)

    zcyl = openmc.ZCylinder(r=1.25*scale)
    box = openmc.model.RectangularPrism(4*scale, 4*scale, boundary_type='reflective')
    cell1 = openmc.Cell(fill=fuel, region=-zcyl)
    cell2 = openmc.Cell(fill=None, region=+zcyl & -box)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.inactive = 0

    cyl_mesh = openmc.CylindricalMesh(
        r_grid=[0., 1.25*scale],
        phi_grid=[0., 2*math.pi],
        z_grid=[-1e10, 1e10]
    )
    cyl_mesh_filter = openmc.MeshFilter(cyl_mesh)
    cell_filter = openmc.CellFilter([cell1])

    tally1 = openmc.Tally()
    tally1.filters = [cyl_mesh_filter]
    tally1.scores = ['flux']
    tally2 = openmc.Tally()
    tally2.filters = [cell_filter]
    tally2.scores = ['flux']
    model.tallies = openmc.Tallies([tally1, tally2])

    # Run OpenMC
    sp_filename = model.run()

    # Get flux for each of the two tallies
    with openmc.StatePoint(sp_filename) as sp:
        t1 = sp.tallies[tally1.id]
        t2 = sp.tallies[tally2.id]
        mean1 = t1.mean.ravel()[0]
        mean2 = t2.mean.ravel()[0]

    # The two tallies should be exactly the same
    assert mean1 == pytest.approx(mean2)


@pytest.mark.parametrize("scale", [0.1, 1.0, 1e2, 1e4, 1e5])
def test_spherical_mesh_coincident(scale, run_in_tmpdir):
    """Test for spherical mesh boundary being coincident with a cell boundary"""

    fuel = openmc.Material()
    fuel.add_nuclide('U235', 1.)
    fuel.set_density('g/cm3', 4.5)

    sph = openmc.Sphere(r=1.25*scale)
    rcc = openmc.model.RectangularParallelepiped(
        -2*scale, 2*scale, -2*scale, 2*scale, -2*scale, 2*scale,
        boundary_type='reflective')
    cell1 = openmc.Cell(fill=fuel, region=-sph)
    cell2 = openmc.Cell(fill=None, region=+sph & -rcc)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.inactive = 0

    sph_mesh = openmc.SphericalMesh(
        r_grid=[0., 1.25*scale],
        phi_grid=[0., 2*math.pi],
        theta_grid=[0., math.pi],
    )

    sph_mesh_filter = openmc.MeshFilter(sph_mesh)
    cell_filter = openmc.CellFilter([cell1])

    tally1 = openmc.Tally()
    tally1.filters = [sph_mesh_filter]
    tally1.scores = ['flux']
    tally2 = openmc.Tally()
    tally2.filters = [cell_filter]
    tally2.scores = ['flux']
    model.tallies = openmc.Tallies([tally1, tally2])

    # Run OpenMC
    sp_filename = model.run()

    # Get flux for each of the two tallies
    with openmc.StatePoint(sp_filename) as sp:
        t1 = sp.tallies[tally1.id]
        t2 = sp.tallies[tally2.id]
        mean1 = t1.mean.ravel()[0]
        mean2 = t2.mean.ravel()[0]

    # The two tallies should be exactly the same
    assert mean1 == pytest.approx(mean2)


def test_get_reshaped_data(run_in_tmpdir):
    """Test that expanding MeshFilter dimensions works as expected"""

    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.set_density('g/cm3', 10.0)

    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 1_000
    model.settings.inactive = 10
    model.settings.batches = 20

    sph_mesh = openmc.SphericalMesh(
        r_grid=np.linspace(0.0, 5.0**3, 20)**(1/3),
        theta_grid=np.linspace(0, math.pi, 4),
        phi_grid=np.linspace(0, 2*math.pi, 3)
    )
    tally1 = openmc.Tally()
    efilter = openmc.EnergyFilter([0, 1e5, 1e8])
    meshfilter = openmc.MeshFilter(sph_mesh)
    assert meshfilter.shape == (19, 3, 2)
    tally1.filters = [efilter, meshfilter]
    tally1.scores = ['flux']

    model.tallies = openmc.Tallies([tally1])

    # Run OpenMC
    sp_filename = model.run()

    # Get flux tally as reshaped data
    with openmc.StatePoint(sp_filename) as sp:
        t1 = sp.tallies[tally1.id]
        data1 = t1.get_reshaped_data()
        data2 = t1.get_reshaped_data(expand_dims=True)

    assert data1.shape == (2, 19*3*2, 1, 1)
    assert data2.shape == (2, 19, 3, 2, 1, 1)

def test_mesh_filter_rotation_roundtrip(run_in_tmpdir):
    """Test that MeshFilter rotation works as expected"""


    mesh = openmc.RegularMesh()
    mesh.lower_left = [-10, -10, -10]
    mesh.upper_right = [10, 10, 10]
    mesh.dimension = [2, 3, 4]

    # check that rotatoin is round-tripped correctly for a set of angles
    mesh_filter = openmc.MeshFilter(mesh)
    mesh_filter.rotation = [0, 0, 90]  # Rotate around z-axis by 90 degrees

    elem = mesh_filter.to_xml_element()
    mesh_filter_xml = openmc.MeshFilter.from_xml_element(elem, meshes={mesh.id: mesh})
    assert all(mesh_filter_xml.rotation == mesh_filter.rotation)

    # check that rotation matrix is round-tripped correctly for a rotation matrix
    mesh_filter.rotation = np.array([[0.7071, 0, 0.7071],
                                     [0, 1, 0],
                                     [-0.7071, 0, 0.7071]])

    elem = mesh_filter.to_xml_element()
    mesh_filter_xml = openmc.MeshFilter.from_xml_element(elem, meshes={mesh.id: mesh})
    assert np.allclose(mesh_filter_xml.rotation, mesh_filter.rotation)
