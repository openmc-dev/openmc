from collections.abc import Iterable
from numbers import Real, Integral
import warnings

from xml.etree import ElementTree as ET
import numpy as np

from openmc.filter import _PARTICLES
from openmc.mesh import MeshBase, RectilinearMesh, UnstructuredMesh
import openmc.checkvalue as cv

from ._xml import get_text
from .mixin import IDManagerMixin


class WeightWindows(IDManagerMixin):
    """Mesh-based weight windows

    This class enables you to specify weight window parameters that are used in
    a simulation. Multiple sets of weight windows can be defined for different
    meshes and different particles. An iterable of :class:`WeightWindows`
    instances can be assigned to the :attr:`openmc.Settings.weight_windows`
    attribute, which is then exported to XML.

    Weight window lower/upper bounds are to be specified for each combination of
    a mesh element and an energy bin. Thus the total number of bounds should be
    equal to the product of the number of mesh bins and the number of energy
    bins.

    .. versionadded:: 0.13

    Parameters
    ----------
    mesh : openmc.MeshBase
        Mesh for the weight windows
    lower_ww_bounds : Iterable of Real
        A list of values for which each value is the lower bound of a weight
        window
    upper_ww_bounds : Iterable of Real
        A list of values for which each value is the upper bound of a weight
        window
    upper_bound_ratio : float
        Ratio of the lower to upper weight window bounds
    energy_bounds : Iterable of Real
        A list of values for which each successive pair constitutes a range of
        energies in [eV] for a single bin
    particle_type : {'neutron', 'photon'}
        Particle type the weight windows apply to
    survival_ratio : float
        Ratio of the survival weight to the lower weight window bound for
        rouletting
    max_lower_bound_ratio : float
        Maximum allowed ratio of a particle's weight to the weight window's
        lower bound. A factor will be applied to raise the weight window to be
        lower than the particle's weight by a factor of max_lower_bound_ratio
        during transport if exceeded.
    max_split : int
        Maximum allowable number of particles when splitting
    weight_cutoff : float
        Threshold below which particles will be terminated
    id : int
       Unique identifier for the weight window settings. If not specified, an
       identifier will automatically be assigned.

    Attributes
    ----------
    id : int
       Unique identifier for the weight window settings.
    mesh : openmc.MeshBase
        Mesh for the weight windows with dimension (ni, nj, nk)
    particle_type : str
        Particle type the weight windows apply to
    energy_bounds : Iterable of Real
        A list of values for which each successive pair constitutes a range of
        energies in [eV] for a single bin
    num_energy_bins : int
        Number of energy bins
    lower_ww_bounds : numpy.ndarray of float
        An array of values for which each value is the lower bound of a weight
        window. Shape: (ni, nj, nk, num_energy_bins) for StructuredMesh;
        (num_elements, num_energy_bins) for UnstructuredMesh
    upper_ww_bounds : numpy.ndarray of float
        An array of values for which each value is the upper bound of a weight
        window. Shape: (ni, nj, nk, num_energy_bins) for StructuredMesh;
        (num_elements, num_energy_bins) for UnstructuredMesh
    survival_ratio : float
        Ratio of the survival weight to the lower weight window bound for
        rouletting
    max_lower_bound_ratio: float
        Maximum allowed ratio of a particle's weight to the weight window's
        lower bound. (Default: 1.0)
    max_split : int
        Maximum allowable number of particles when splitting
    weight_cutoff : float
        Threshold below which particles will be terminated

    See Also
    --------
    openmc.Settings

    """
    next_id = 1
    used_ids = set()

    def __init__(self, mesh, lower_ww_bounds,
                 upper_ww_bounds=None,
                 upper_bound_ratio=None,
                 energy_bounds=None,
                 particle_type='neutron',
                 survival_ratio=3,
                 max_lower_bound_ratio=None,
                 max_split=10,
                 weight_cutoff=1.e-38,
                 id=None):
        self.mesh = mesh
        self.id = id
        self.particle_type = particle_type
        self.energy_bounds = energy_bounds
        self.lower_ww_bounds = lower_ww_bounds

        cv.check_length('Lower window bounds',
                        self.lower_ww_bounds,
                        len(self.energy_bounds))

        if upper_ww_bounds is not None and upper_bound_ratio:
            raise ValueError("Exactly one of upper_ww_bounds and "
                             "upper_bound_ratio must be present.")

        if upper_ww_bounds is None and upper_bound_ratio is None:
            raise ValueError("Exactly one of upper_ww_bounds and "
                             "upper_bound_ratio must be present.")

        if upper_bound_ratio:
            self.upper_ww_bounds = [
                lb * upper_bound_ratio for lb in self.lower_ww_bounds
            ]

        if upper_ww_bounds is not None:
            self.upper_ww_bounds = upper_ww_bounds

        if len(self.lower_ww_bounds) != len(self.upper_ww_bounds):
            raise ValueError('Size of the lower and upper weight '
                             'window bounds do not match')

        self.survival_ratio = survival_ratio

        self._max_lower_bound_ratio = None
        if max_lower_bound_ratio is not None:
            self.max_lower_bound_ratio = max_lower_bound_ratio

        self.max_split = max_split
        self.weight_cutoff = weight_cutoff

    def __repr__(self):
        string = type(self).__name__ + '\n'
        string += '{: <16}=\t{}\n'.format('\tID', self._id)
        string += '{: <16}=\t{}\n'.format('\tMesh:', self.mesh)
        string += '{: <16}=\t{}\n'.format('\tParticle Type', self._particle_type)
        string += '{: <16}=\t{}\n'.format('\tEnergy Bounds', self._energy_bounds)
        string += '{: <16}=\t{}\n'.format('\tLower WW Bounds', self._lower_ww_bounds)
        string += '{: <16}=\t{}\n'.format('\tUpper WW Bounds', self._upper_ww_bounds)
        string += '{: <16}=\t{}\n'.format('\tSurvival Ratio', self._survival_ratio)
        string += '{: <16}=\t{}\n'.format('\tMax Split', self._max_split)
        string += '{: <16}=\t{}\n'.format('\tWeight Cutoff', self._weight_cutoff)
        return string

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        cv.check_type('Weight window mesh', mesh, MeshBase)
        self._mesh = mesh

    @property
    def particle_type(self):
        return self._particle_type

    @particle_type.setter
    def particle_type(self, pt):
        cv.check_value('Particle type', pt, _PARTICLES)
        self._particle_type = pt

    @property
    def energy_bounds(self):
        return self._energy_bounds

    @energy_bounds.setter
    def energy_bounds(self, bounds):
        cv.check_type('Energy bounds', bounds, Iterable, Real)
        self._energy_bounds = np.asarray(bounds)

    @property
    def num_energy_bins(self):
        if self.energy_bounds is None:
            raise ValueError('Energy bounds are not set')
        return self.energy_bounds.size - 1

    @property
    def lower_ww_bounds(self):
        return self._lower_ww_bounds

    @lower_ww_bounds.setter
    def lower_ww_bounds(self, bounds):
        cv.check_iterable_type('Lower WW bounds',
                               bounds,
                               Real,
                               min_depth=1,
                               max_depth=4)
        # reshape data according to mesh and energy bins
        bounds = np.asarray(bounds)
        if isinstance(self.mesh, UnstructuredMesh):
            bounds.reshape(-1, self.num_energy_bins)
        else:
            bounds.reshape(*self.mesh.dimension, self.num_energy_bins)
        self._lower_ww_bounds = bounds

    @property
    def upper_ww_bounds(self):
        return self._upper_ww_bounds

    @upper_ww_bounds.setter
    def upper_ww_bounds(self, bounds):
        cv.check_iterable_type('Upper WW bounds',
                               bounds,
                               Real,
                               min_depth=1,
                               max_depth=4)
        # reshape data according to mesh and energy bins
        bounds = np.asarray(bounds)
        if isinstance(self.mesh, UnstructuredMesh):
            bounds.reshape(-1, self.num_energy_bins)
        else:
            bounds.reshape(*self.mesh.dimension, self.num_energy_bins)
        self._upper_ww_bounds = bounds

    @property
    def survival_ratio(self):
        return self._survival_ratio

    @survival_ratio.setter
    def survival_ratio(self, val):
        cv.check_type('Survival ratio', val, Real)
        cv.check_greater_than('Survival ratio', val, 1.0, True)
        self._survival_ratio = val

    @property
    def max_lower_bound_ratio(self):
        return self._max_lower_bound_ratio

    @max_lower_bound_ratio.setter
    def max_lower_bound_ratio(self, val):
        cv.check_type('Maximum lower bound ratio', val, Real)
        cv.check_greater_than('Maximum lower bound ratio', val, 1.0)
        self._max_lower_bound_ratio = val

    @property
    def max_split(self):
        return self._max_split

    @max_split.setter
    def max_split(self, val):
        cv.check_type('Max split', val, Integral)
        self._max_split = val

    @property
    def weight_cutoff(self):
        return self._weight_cutoff

    @weight_cutoff.setter
    def weight_cutoff(self, cutoff):
        cv.check_type('Weight cutoff', cutoff, Real)
        cv.check_greater_than('Weight cutoff', cutoff, 0.0, True)
        self._weight_cutoff = cutoff

    def to_xml_element(self):
        """Return an XML representation of the weight window settings

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing the weight window information
        """
        element = ET.Element('weight_windows')

        element.set('id', str(self._id))

        subelement = ET.SubElement(element, 'mesh')
        subelement.text = str(self.mesh.id)

        subelement = ET.SubElement(element, 'particle_type')
        subelement.text = self.particle_type

        subelement = ET.SubElement(element, 'energy_bounds')
        subelement.text = ' '.join(str(e) for e in self.energy_bounds)

        subelement = ET.SubElement(element, 'lower_ww_bounds')
        subelement.text = ' '.join(str(b) for b in self.lower_ww_bounds)

        subelement = ET.SubElement(element, 'upper_ww_bounds')
        subelement.text = ' '.join(str(b) for b in self.upper_ww_bounds)

        subelement = ET.SubElement(element, 'survival_ratio')
        subelement.text = str(self.survival_ratio)

        if self.max_lower_bound_ratio is not None:
            subelement = ET.SubElement(element, 'max_lower_bound_ratio')
            subelement.text = str(self.max_lower_bound_ratio)

        subelement = ET.SubElement(element, 'max_split')
        subelement.text = str(self.max_split)

        subelement = ET.SubElement(element, 'weight_cutoff')
        subelement.text = str(self.weight_cutoff)

        return element

    @classmethod
    def from_xml_element(cls, elem, root):
        """Generate weight window settings from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element
        root : xml.etree.ElementTree.Element
            Root element for the file where meshes can be found

        Returns
        -------
        openmc.WeightWindows
            Weight windows object
        """
        # Get mesh for weight windows
        mesh_id = int(get_text(elem, 'mesh'))
        path = f"./mesh[@id='{mesh_id}']"
        mesh_elem = root.find(path)
        if mesh_elem is not None:
            mesh = MeshBase.from_xml_element(mesh_elem)

        # Read all other parameters
        lower_ww_bounds = [float(l) for l in get_text(elem, 'lower_ww_bounds').split()]
        upper_ww_bounds = [float(u) for u in get_text(elem, 'upper_ww_bounds').split()]
        e_bounds = [float(b) for b in get_text(elem, 'energy_bounds').split()]
        particle_type = get_text(elem, 'particle_type')
        survival_ratio = float(get_text(elem, 'survival_ratio'))

        max_lower_bound_ratio = None
        if get_text(elem, 'max_lower_bound_ratio'):
            max_lower_bound_ratio = float(get_text(elem, 'max_lower_bound_ratio'))

        max_split = int(get_text(elem, 'max_split'))
        weight_cutoff = float(get_text(elem, 'weight_cutoff'))
        id = int(get_text(elem, 'id'))

        return cls(
            mesh=mesh,
            lower_ww_bounds=lower_ww_bounds,
            upper_ww_bounds=upper_ww_bounds,
            energy_bounds=e_bounds,
            particle_type=particle_type,
            survival_ratio=survival_ratio,
            max_lower_bound_ratio=max_lower_bound_ratio,
            max_split=max_split,
            weight_cutoff=weight_cutoff,
            id=id
        )

    @classmethod
    def from_hdf5(cls, group, meshes):
        """Create weight windows from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            Group in HDF5 file
        meshes : dict
            Dictionary mapping IDs to mesh objects

        Returns
        -------
        openmc.WeightWindows
            A weight window object
        """

        id = int(group.name.split('/')[-1].lstrip('weight_windows'))
        mesh_id = group['mesh'][()]
        ptype = group['particle_type'][()].decode()
        e_bounds = group['energy_bounds'][()]
        lower_ww_bounds = group['lower_ww_bounds'][()]
        upper_ww_bounds = group['upper_ww_bounds'][()]
        survival_ratio = group['survival_ratio'][()]

        max_lower_bound_ratio = None
        if group.get('max_lower_bound_ratio') is not None:
            max_lower_bound_ratio = group['max_lower_bound_ratio'][()]

        max_split = group['max_split'][()]
        weight_cutoff = group['weight_cutoff'][()]

        return cls(
            mesh=meshes[mesh_id],
            lower_ww_bounds=lower_ww_bounds,
            upper_ww_bounds=upper_ww_bounds,
            energy_bounds=e_bounds,
            particle_type=ptype,
            survival_ratio=survival_ratio,
            max_lower_bound_ratio=max_lower_bound_ratio,
            max_split=max_split,
            weight_cutoff=weight_cutoff,
            id=id
        )

def __wwinp_reader(path):
    """
    Generator that returns the next value in a wwinp file.

    path : str or pathlib.Path
        Location of the wwinp file
    """
    fh = open(path, 'r')

    # read the first line of the file and
    # keep only the first four entries
    while(True):
        line = next(fh)
        if line and not line.startswith('c'):
            break

    values = line.strip().split()[:4]
    for value in values:
        yield value

    # the remainder of the file can be read as
    # sequential values
    while(True):
        line = next(fh)
        # skip empty or commented lines
        if not line or line.startswith('c'):
            continue
        values = line.strip().split()
        for value in values:
            yield value

def wwinp_to_wws(path):
    """Creates WeightWindows classes from a wwinp file

    Parameters
    ----------
    path : str
        Path to the wwinp file.

    Returns
    -------
    list of openmc.WeightWindows
    """
    # create generator for getting the next parameter from the file
    wwinp = __wwinp_reader(path)

    # first parameter, if, of wwinp file is unused
    next(wwinp)

    # check time parameter, iv
    if int(float(next(wwinp))) > 1:
        raise ValueError('Time-dependent weight windows are not yet supported.')

    # number of particle types, ni
    n_particle_types = int(float(next(wwinp)))

    # read an indicator of the mesh type.
    # this will be 10 if a rectilinear mesh
    # and 16 for cylindrical or spherical meshes
    mesh_chars = int(float(next(wwinp)))

    if mesh_chars != 10:
        # TODO: read the first entry by default and display a warning
        raise ValueError('Cylindrical and Spherical meshes are not currently supported')

    # read the number of energy groups for each particle, ne
    n_egroups = [int(next(wwinp)) for _ in range(n_particle_types)]

    # order that supported particle types will appear in the file
    particle_types = ['neutron', 'photon']

    # add particle to list if at least one energy group is present
    particles = [p for e, p in zip(n_egroups, particle_types) if e > 0]
    # truncate list of energy groups if needed
    n_egroups = [e for e in n_egroups if e > 0]

    if n_particle_types > 2:
        msg = ('More than two particle types are present. '
                'Only neutron and photon weight windows will be read.')
        warnings.warn(msg)

    # read total number of fine mesh elements in each coarse
    # element (nfx, nfy, nfz)
    n_fine_x = int(float(next(wwinp)))
    n_fine_y = int(float(next(wwinp)))
    n_fine_z = int(float(next(wwinp)))
    header_mesh_dims = (n_fine_x, n_fine_y, n_fine_z)

    # read the mesh origin: x0, y0, z0
    llc = tuple(float(next(wwinp)) for _ in range(3))

    # read the number of coarse mesh elements (ncx, ncy, ncz)
    n_coarse_x = int(float(next(wwinp)))
    n_coarse_y = int(float(next(wwinp)))
    n_coarse_z = int(float(next(wwinp)))

    # skip the value defining the geometry type, nwg, we already know this
    # 1 - rectilinear mesh
    # 2 - cylindrical mesh
    # 3 - spherical mesh
    mesh_type = int(float(next(wwinp)))

    if mesh_type != 1:
        # TODO: support additional mesh types
        raise ValueError('Cylindrical and Spherical meshes are not currently supported')

    # internal function for parsing mesh coordinates
    def _read_mesh_coords(wwinp, n_coarse_bins):
        coords = [float(next(wwinp))]

        for _ in range(n_coarse_bins):
            # number of fine mesh elements in this coarse element, sx
            sx = int(float(next(wwinp)))
            # value of next coordinate, px
            px = float(next(wwinp))
            # fine mesh ratio, qx, is currently unused
            qx = next(wwinp)
            # append the fine mesh coordinates for this coarse element
            coords += list(np.linspace(coords[-1], px, sx + 1))[1:]

        return np.asarray(coords)

    # read the coordinates for each dimension into a rectilinear mesh
    mesh = RectilinearMesh()
    mesh.x_grid = _read_mesh_coords(wwinp, n_coarse_x)
    mesh.y_grid = _read_mesh_coords(wwinp, n_coarse_y)
    mesh.z_grid = _read_mesh_coords(wwinp, n_coarse_z)

    dims = ('x', 'y', 'z')
    # check consistency of mesh coordinates
    mesh_llc = mesh_val = (mesh.x_grid[0], mesh.y_grid[0], mesh.z_grid[0])
    for dim, header_val, mesh_val in zip(dims, llc, mesh_llc):
        if header_val != mesh_val:
            msg = ('The {} corner of the mesh ({}) does not match '
                    'the value read in block 1 of the wwinp file ({})')
            raise ValueError(msg.format(dim, mesh_val, header_val))

    # check total number of mesh elements in each direction
    mesh_dims = mesh.dimension
    for dim, header_val, mesh_val in zip(dims, header_mesh_dims, mesh_dims):
        if header_val != mesh_val:
            msg = ('Total number of mesh elements read in the {} '
                    'direction ({}) is inconsistent with the '
                    'number read in block 1 of the wwinp file ({})')
            raise ValueError(msg.format(dim, mesh_val, header_val))

    # total number of fine mesh elements, nft
    n_elements = n_fine_x * n_fine_y * n_fine_z
    # read energy bins and weight window values for each particle
    wws = []
    for particle, ne in zip(particles, n_egroups):
        # read upper energy bounds
        # it is implied that zero is always the first bound in MCNP
        e_bounds = np.asarray([0.0] + [float(next(wwinp)) for _ in range(ne)])
        # adjust energy from MeV to eV
        e_bounds *= 1E6

        # create an array for weight window lower bounds
        ww_lb = np.zeros((*mesh.dimension, ne))
        for ijk in mesh.indices:
            idx = tuple([v - 1 for v in ijk] + [slice(None)])
            ww_lb[idx] = [float(next(wwinp)) for _ in range(ne)]
        # for e in range(ne):
        #     # MCNP ordering for weight windows matches that of OpenMC
        #     # ('xyz' with x changing fastest)
        #     ww_vals = [float(next(wwinp)) for _ in range(n_elements)]
        #     ww_lb[e, :] = np.asarray(ww_vals).reshape(mesh.dimension)
        settings = WeightWindows(id=None,
                                    mesh=mesh,
                                    lower_ww_bounds=ww_lb.flatten(),
                                    upper_bound_ratio=5.0,
                                    energy_bins=e_bounds,
                                    particle_type=particle)
        wws.append(settings)

    return wws