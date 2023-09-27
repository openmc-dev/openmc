
import numpy as np
import pytest

import openmc
import openmc.model

@pytest.mark.parametrize('mesh_type', ('rectangular', 'cylindrical'))
def test_source_mesh(mesh_type):
    """
    A void model containing a single box
    """
    min, max = -10, 10
    box = openmc.model.RectangularParallelepiped(min, max, min, max, min, max, boundary_type='vacuum')

    geometry = openmc.Geometry([openmc.Cell(region=-box)])

    settings = openmc.Settings()
    settings.particles = 100
    settings.batches = 10
    settings.run_mode = 'fixed source'

    model = openmc.Model(geometry=geometry, settings=settings)

    # define a 2 x 2 x 2 mesh
    if mesh_type == 'rectangular':
        mesh = openmc.RegularMesh.from_domain(model.geometry, (2, 2, 2))
    elif mesh_type == 'cylindrical':
        mesh = openmc.CylindricalMesh.from_domain(model.geometry, (1, 4, 2))

    energy = openmc.stats.Discrete([1.e6], [1.0])

    # create sources with only one non-zero strength for the source in the
    # mesh voxel occupyting the lowest octant. Direct source particles straight
    # out of the problem from there. This demonstrates that
    # 1) particles are only being sourced within the intented mesh voxel based on source strength
    # 2) particles are respecting the angle distributions assigned to each voxel
    x, y, z = mesh.centroids
    if mesh_type == 'cylindrical':
        x, y, z = mesh._convert_to_cartesian(mesh.centroids, mesh.origin)

    sources = np.ndarray(mesh.dimension, dtype=openmc.SourceBase)
    for i, j, k in mesh.indices:
        # mesh.indices is currently one-indexed, adjust for Python arrays
        idx = (i-1, j-1, k-1)

        # get the centroid of the ijk mesh element, set a particle source
        # vector based on the
        centroid = np.array((x[idx], y[idx], z[idx]))
        vec = np.sign(centroid, dtype=float)
        vec /= np.linalg.norm(vec)
        angle = openmc.stats.Monodirectional(vec)
        sources[idx] = openmc.Source(energy=energy, angle=angle, strength=0.0)

    # create and apply the mesh source
    mesh_source = openmc.MeshSource(mesh, sources)
    model.settings.source = mesh_source

    # tally the flux on the mesh
    mesh_filter = openmc.MeshFilter(mesh)
    tally = openmc.Tally()
    tally.filters = [mesh_filter]
    tally.scores = ['flux']

    model.tallies = openmc.Tallies([tally])

    # for each element, set a single-non zero source with particles
    # traveling out of the mesh (and geometry) w/o crossing any other
    # mesh elements
    for i, j, k in mesh.indices:
        ijk = (i-1, j-1, k-1)
        print(ijk)
        # zero-out all source strengths and set the strength
        # on the element of interest
        mesh_source.strength = 0.0
        mesh_source.sources[ijk].strength = 1.0

        sp_file = model.run()

        with openmc.StatePoint(sp_file) as sp:
            tally_out = sp.get_tally(id=tally.id)
            mean = tally_out.get_reshaped_data(expand_dims=True)

        # remove nuclides and scores axes
        mean = mean[..., 0, 0]

        assert mean[ijk] != 0
        mean[ijk] = 0
        assert np.all(mean == 0)

    # check strength adjustment methods
    assert mesh_source.strength == 1.0
    mesh_source.strength = 100.0
    assert mesh_source.strength == 100.0

