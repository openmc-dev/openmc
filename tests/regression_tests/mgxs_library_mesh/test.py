import hashlib

import openmc
import openmc.mgxs
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.set_density('g/cm3', 10.0)
    fuel.add_nuclide('U235', 1.0)
    zr = openmc.Material()
    zr.set_density('g/cm3', 1.0)
    zr.add_nuclide('Zr90', 1.0)
    model.materials.extend([fuel, zr])

    box1 = openmc.model.rectangular_prism(10.0, 10.0)
    box2 = openmc.model.rectangular_prism(20.0, 20.0, boundary_type='reflective')
    top = openmc.ZPlane(z0=10.0, boundary_type='vacuum')
    bottom = openmc.ZPlane(z0=-10.0, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=fuel, region=box1 & +bottom & -top)
    cell2 = openmc.Cell(fill=zr, region=~box1 & box2 & +bottom & -top)
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000

    # Initialize a one-group structure
    energy_groups = openmc.mgxs.EnergyGroups([0, 20.e6])

    # Initialize MGXS Library for a few cross section types
    # for one material-filled cell in the geometry
    model.mgxs_lib = openmc.mgxs.Library(model.geometry)
    model.mgxs_lib.by_nuclide = False

    # Test all MGXS types
    model.mgxs_lib.mgxs_types = openmc.mgxs.MGXS_TYPES + openmc.mgxs.MDGXS_TYPES
    model.mgxs_lib.energy_groups = energy_groups
    model.mgxs_lib.num_delayed_groups = 6
    model.mgxs_lib.correction = None  # Avoid warning about P0 correction
    model.mgxs_lib.legendre_order = 3
    model.mgxs_lib.domain_type = 'mesh'

    # Instantiate a tally mesh
    mesh = openmc.RegularMesh(mesh_id=1)
    mesh.dimension = [2, 2]
    mesh.lower_left = [-100., -100.]
    mesh.width = [100., 100.]

    model.mgxs_lib.domains = [mesh]
    model.mgxs_lib.build_library()

    # Add tallies
    model.mgxs_lib.add_to_tallies_file(model.tallies, merge=False)

    return model


class MGXSTestHarness(PyAPITestHarness):
    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        # Load the MGXS library from the statepoint
        mgxs_lib = self._model.mgxs_lib
        mgxs_lib.load_from_statepoint(sp)

        # Build a string from Pandas Dataframe for each 1-group MGXS
        outstr = ''
        for domain in mgxs_lib.domains:
            for mgxs_type in mgxs_lib.mgxs_types:
                mgxs = mgxs_lib.get_mgxs(domain, mgxs_type)
                df = mgxs.get_pandas_dataframe()
                outstr += df.to_string() + '\n'

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr


def test_mgxs_library_mesh(model):
    harness = MGXSTestHarness('statepoint.5.h5', model)
    harness.main()
