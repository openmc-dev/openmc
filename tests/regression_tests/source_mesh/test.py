
import numpy as np
import pytest

import openmc
import openmc.model


def test_source_mesh():
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
    mesh = openmc.RegularMesh.from_domain(model.geometry, (2, 2, 2))

    strengths = np.ones((2, 2, 2))

    energy = openmc.stats.Discrete([1.e6], [1.0])

    # create sources with only one non-zero strength for the source in the
    # mesh voxel occupyting the lowest octant. Direct source particles straight
    # out of the problem from there. This demonstrates that
    # 1) particles are only being sourced within the intented mesh voxel based on source strength
    # 2) particles are respecting the angle distributions assigned to each voxel

    x, y, z = mesh.centroids

    sources = []
    for idx in mesh.indices:

        idx = tuple(i-1 for i in idx)
        centroid = np.array((x[idx], y[idx], z[idx]))
        print(centroid)
        vec = np.sign(centroid, dtype=float)
        print(vec)
        vec /= np.linalg.norm(vec)
        angle = openmc.stats.Monodirectional(vec)
        sources.append(openmc.Source(energy=energy, angle=angle, strength=0.0))

    for idx in range(mesh.num_mesh_cells):

        for s in sources:
            s.strength = 0.0

        sources[idx].strength = 1.0

        mesh_source = openmc.MeshSource(mesh, sources=sources)

        model.settings.source = mesh_source


        # tally the flux on the mesh
        mesh_filter = openmc.MeshFilter(mesh)
        tally = openmc.Tally()
        tally.filters = [mesh_filter]
        tally.scores = ['flux']

        model.tallies = openmc.Tallies([tally])

        sp_file = model.run()

        with openmc.StatePoint(sp_file) as sp:
            tally_out = sp.get_tally(id=tally.id)
            mean = tally_out.mean

        assert mean[idx] != 0
        mean[idx] = 0
        assert np.all(mean == 0)

