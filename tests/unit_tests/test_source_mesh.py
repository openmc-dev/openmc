
import numpy as np
import pytest

import openmc
import openmc.lib
import openmc.model

@pytest.mark.parametrize('mesh_type', ('rectangular', 'cylindrical'))
def test_source_mesh(run_in_tmpdir, mesh_type):
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
    sources = np.ndarray(mesh.dimension, dtype=openmc.SourceBase)
    centroids = mesh.centroids
    x, y, z = np.swapaxes(mesh.centroids, -1, 0)
    for i, j, k in mesh.indices:
        # mesh.indices is currently one-indexed, adjust for Python arrays
        ijk = (i-1, j-1, k-1)

        # get the centroid of the ijk mesh element, set a particle source
        # vector based on the
        centroid = centroids[ijk]
        vec = np.sign(centroid, dtype=float)
        vec /= np.linalg.norm(vec)
        angle = openmc.stats.Monodirectional(vec)
        sources[ijk] = openmc.IndependentSource(energy=energy, angle=angle, strength=0.0)

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
        # the mesh elment with a non-zero source strength should have a value
        assert mean[ijk] != 0
        # all other values should be zero
        mean[ijk] = 0
        assert np.all(mean == 0), f'Failed on index {ijk} with centroid {mesh.centroids[ijk]}'

        # test roundtrip
        xml_model = openmc.Model.from_model_xml()
        xml_source = xml_model.settings.source[0]
        assert isinstance(xml_source, openmc.MeshSource)
        assert xml_source.strength == 1.0
        assert isinstance(xml_source.mesh, type(mesh_source.mesh))
        assert xml_source.mesh.dimension == mesh_source.mesh.dimension
        assert xml_source.mesh.id == mesh_source.mesh.id
        assert len(xml_source.sources) == len(mesh_source.sources)

    # check strength adjustment methods
    assert mesh_source.strength == 1.0
    mesh_source.strength = 100.0
    assert mesh_source.strength == 100.0

    mesh_source.normalize_source_strengths()
    assert mesh_source.strength == 1.0


def test_file_source(run_in_tmpdir):
    source_particle = openmc.SourceParticle(r=(0.0, 0.0, 0.0),
                                            u=(0.0, 0.0, 1.0),
                                            E=1e6,
                                            time=10.0)

    openmc.write_source_file([source_particle], 'source.h5')

    file_source = openmc.FileSource('source.h5')

    model = openmc.Model()

    rect_prism = openmc.model.RectangularParallelepiped(-5.0, 5.0,
                                                        -5.0, 5.0,
                                                        -5.0, 5.0,
                                                        boundary_type='vacuum')
    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)

    model.geometry = openmc.Geometry([openmc.Cell(fill=mat, region=-rect_prism)])
    model.settings.particles = 1000
    model.settings.batches = 10
    model.settings.inactive = 1
    model.settings.run_mode = 'fixed source'

    mesh = openmc.RegularMesh()
    mesh.lower_left = (-1, -2, -3)
    mesh.upper_right = (2, 3, 4)
    mesh.dimension = (1, 1, 1)

    mesh_source_arr = np.asarray([file_source]).reshape(mesh.dimension)
    source = openmc.MeshSource(mesh, mesh_source_arr)

    model.settings.source = source

    model.export_to_model_xml()

    openmc.lib.init()
    openmc.lib.simulation_init()
    sites = openmc.lib.sample_external_source(10)
    openmc.lib.simulation_finalize()
    openmc.lib.finalize()


    # the mesh bounds do not contain the point of the lone
    # source site in the file source, so it should not appear
    # in the set of source sites produced from the mesh source.
    # Additionally, the source should be located within the mesh
    bbox = mesh.bounding_box
    for site in sites:
        assert site.r != (0, 0, 0)
        assert site.E == source_particle.E
        assert site.u == source_particle.u
        assert site.time == source_particle.time
        assert site.r in bbox

