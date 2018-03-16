import openmc

from tests.testing_harness import HashedPyAPITestHarness


class FilterMeshTestHarness(HashedPyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Initialize Meshes
        mesh_1d = openmc.Mesh(mesh_id=1)
        mesh_1d.type = 'regular'
        mesh_1d.dimension = [17]
        mesh_1d.lower_left = [-182.07]
        mesh_1d.upper_right = [182.07]

        mesh_2d = openmc.Mesh(mesh_id=2)
        mesh_2d.type = 'regular'
        mesh_2d.dimension = [17, 17]
        mesh_2d.lower_left = [-182.07, -182.07]
        mesh_2d.upper_right = [182.07, 182.07]

        mesh_3d = openmc.Mesh(mesh_id=3)
        mesh_3d.type = 'regular'
        mesh_3d.dimension = [17, 17, 17]
        mesh_3d.lower_left = [-182.07, -182.07, -183.00]
        mesh_3d.upper_right = [182.07, 182.07, 183.00]

        # Initialize the filters
        mesh_1d_filter = openmc.MeshFilter(mesh_1d)
        mesh_2d_filter = openmc.MeshFilter(mesh_2d)
        mesh_3d_filter = openmc.MeshFilter(mesh_3d)
        meshsurf_1d_filter = openmc.MeshSurfaceFilter(mesh_1d)
        meshsurf_2d_filter = openmc.MeshSurfaceFilter(mesh_2d)
        meshsurf_3d_filter = openmc.MeshSurfaceFilter(mesh_3d)

        # Initialized the tallies
        tally = openmc.Tally(name='tally 1')
        tally.filters = [mesh_1d_filter]
        tally.scores = ['total']
        self._model.tallies.append(tally)

        tally = openmc.Tally(name='tally 2')
        tally.filters = [meshsurf_1d_filter]
        tally.scores = ['current']
        self._model.tallies.append(tally)

        tally = openmc.Tally(name='tally 3')
        tally.filters = [mesh_2d_filter]
        tally.scores = ['total']
        self._model.tallies.append(tally)

        tally = openmc.Tally(name='tally 4')
        tally.filters = [meshsurf_2d_filter]
        tally.scores = ['current']
        self._model.tallies.append(tally)

        tally = openmc.Tally(name='tally 5')
        tally.filters = [mesh_3d_filter]
        tally.scores = ['total']
        self._model.tallies.append(tally)

        tally = openmc.Tally(name='tally 6')
        tally.filters = [meshsurf_3d_filter]
        tally.scores = ['current']
        self._model.tallies.append(tally)


def test_filter_mesh():
    harness = FilterMeshTestHarness('statepoint.10.h5')
    harness.main()
