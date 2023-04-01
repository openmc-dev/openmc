from collections.abc import Iterable
from numbers import Real, Integral

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

    def __eq__(self, other):
        # ensure that `other` is a WeightWindows object
        if not isinstance(other, WeightWindows):
            return False

        # TODO: add ability to check mesh equality

        # check several attributes directly
        attrs = ('particle_type',
                 'survival_ratio',
                 'max_lower_bound_ratio',
                 'max_split',
                 'weight_cutoff')
        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False

        # save most expensive checks for last
        if not np.array_equal(self.energy_bounds, other.energy_bounds):
            return False

        if not np.array_equal(self.lower_ww_bounds, other.lower_ww_bounds):
            return False

        if not np.array_equal(self.upper_ww_bounds, other.upper_ww_bounds):
            return False

        return True

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
            bounds = bounds.reshape(-1, self.num_energy_bins)
        else:
            bounds = bounds.reshape(*self.mesh.dimension, self.num_energy_bins)
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
            bounds = bounds.reshape(-1, self.num_energy_bins)
        else:
            bounds = bounds.reshape(*self.mesh.dimension, self.num_energy_bins)
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
        subelement.text = ' '.join(str(b) for b in self.lower_ww_bounds.ravel('F'))

        subelement = ET.SubElement(element, 'upper_ww_bounds')
        subelement.text = ' '.join(str(b) for b in self.upper_ww_bounds.ravel('F'))

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

        ww_shape = (len(e_bounds) - 1,) + mesh.dimension[::-1]
        lower_ww_bounds = np.array(lower_ww_bounds).reshape(ww_shape).T
        upper_ww_bounds = np.array(upper_ww_bounds).reshape(ww_shape).T

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


def wwinp_to_wws(path):
    """Create WeightWindows instances from a wwinp file

    .. versionadded:: 0.13.1

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the wwinp file

    Returns
    -------
    list of openmc.WeightWindows
    """

    with open(path) as wwinp:
        # BLOCK 1
        header = wwinp.readline().split(None, 4)
        # read file type, time-dependence, number of
        # particles, mesh type and problem identifier
        _if, iv, ni, nr = [int(x) for x in header[:4]]
        probid = header[4] if len(header) > 4 else ""

        # header value checks
        if _if != 1:
            raise ValueError(f'Found incorrect file type, if: {_if}')

        if iv > 1:
            # read number of time bins for each particle, 'nt(1...ni)'
            nt = np.fromstring(wwinp.readline(), sep=' ', dtype=int)

            # raise error if time bins are present for now
            raise ValueError('Time-dependent weight windows '
                             'are not yet supported')
        else:
            nt = ni * [1]

        if nr == 16:
            raise NotImplementedError('Cylindrical and spherical mesh '
                                      'types are not yet supported')

        # read number of energy bins for each particle, 'ne(1...ni)'
        ne = np.fromstring(wwinp.readline(), sep=' ', dtype=int)

        # read coarse mesh dimensions and lower left corner
        mesh_description = np.fromstring(wwinp.readline(), sep=' ')
        nfx, nfy, nfz = mesh_description[:3].astype(int)
        xyz0 = mesh_description[3:]

        # read cylindrical and spherical mesh vectors if present
        if nr == 16:
            # read number of coarse bins
            line_arr = np.fromstring(wwinp.readline(), sep=' ')
            ncx, ncy, ncz = line_arr[:3].astype(int)
            # read polar vector (x1, y1, z1)
            xyz1 = line_arr[3:] - xyz0
            polar_vec = xyz1 / np.linalg.norm(xyz1)
            # read azimuthal vector (x2, y2, z2)
            line_arr = np.fromstring(wwinp.readline(), sep=' ')
            xyz2 = line_arr[:3] - xyz0
            azimuthal_vec = xyz2 / np.linalg.norm(xyz2)
            # read geometry type
            nwg = int(line_arr[-1])
        elif nr == 10:
            # read rectilinear data:
            # number of coarse mesh bins and mesh type
            ncx, ncy, ncz, nwg = \
                np.fromstring(wwinp.readline(), sep=' ').astype(int)
        else:
            raise RuntimeError(f'Invalid mesh description (nr) found: {nr}')

        # read BLOCK 2 and BLOCK 3 data into a single array
        ww_data = np.fromstring(wwinp.read(), sep=' ')

    # extract mesh data from the ww_data array
    start_idx = 0

    # first values in the mesh definition arrays are the first
    # coordinate of the grid
    end_idx = start_idx + 1 + 3 * ncx
    x0, x_vals = ww_data[start_idx], ww_data[start_idx+1:end_idx]
    start_idx = end_idx

    end_idx = start_idx + 1 + 3 * ncy
    y0, y_vals = ww_data[start_idx], ww_data[start_idx+1:end_idx]
    start_idx = end_idx

    end_idx = start_idx + 1 + 3 * ncz
    z0, z_vals = ww_data[start_idx], ww_data[start_idx+1:end_idx]
    start_idx = end_idx

    # mesh consistency checks
    if nr == 16 and nwg == 1:
        raise ValueError(f'Mesh description in header ({nr}) '
                         f'does not match the mesh type ({nwg})')

    if (xyz0 != (x0, y0, z0)).any():
        raise ValueError(f'Mesh origin in the header ({xyz0}) '
                         f' does not match the origin in the mesh '
                         f' description ({x0, y0, z0})')

    # create openmc mesh object
    grids = []
    mesh_definition = [(x0, x_vals, nfx), (y0, y_vals, nfy), (z0, z_vals, nfz)]
    for grid0, grid_vals, n_pnts in mesh_definition:
        # file spec checks for the mesh definition
        if (grid_vals[2::3] != 1.0).any():
            raise ValueError('One or more mesh ratio value, qx, '
                             'is not equal to one')

        if grid_vals[::3].sum() != n_pnts:
            raise ValueError('Sum of the fine bin entries, s, does '
                             'not match the number of fine bins')

        # extend the grid based on the next coarse bin endpoint, px
        # and the number of fine bins in the coarse bin, sx
        intervals = grid_vals.reshape(-1, 3)
        coords = [grid0]
        for sx, px, qx in intervals:
            coords += np.linspace(coords[-1], px, int(sx + 1)).tolist()[1:]

        grids.append(np.array(coords))

    mesh = RectilinearMesh()
    mesh.x_grid, mesh.y_grid, mesh.z_grid = grids

    # extract weight window values from array
    wws = []
    for ne_i, nt_i, particle_type in zip(ne, nt, ('neutron', 'photon')):
        # no information to read for this particle if
        # either the energy bins or time bins are empty
        if ne_i == 0 or nt_i == 0:
            continue

        if iv > 1:
            # time bins are parsed but unused for now
            end_idx = start_idx + nt_i
            time_bounds = ww_data[start_idx:end_idx]
            np.insert(time_bounds, (0,), (0.0,))
            start_idx = end_idx

        # read energy boundaries
        end_idx = start_idx + ne_i
        energy_bounds = np.insert(ww_data[start_idx:end_idx], (0,), (0.0,))
        # convert from MeV to eV
        energy_bounds *= 1e6
        start_idx = end_idx

        # read weight window values
        end_idx = start_idx + (nfx * nfy * nfz) * nt_i * ne_i

        # read values and reshape according to ordering
        # slowest to fastest: t, e, z, y, x
        # reorder with transpose since our ordering is x, y, z, e, t
        ww_shape = (nt_i, ne_i, nfz, nfy, nfx)
        ww_values = ww_data[start_idx:end_idx].reshape(ww_shape).T
        # Only use first time bin since we don't support time dependent weight
        # windows yet.
        ww_values = ww_values[:, :, :, :, 0]
        start_idx = end_idx

        # create a weight window object
        ww = WeightWindows(id=None,
                           mesh=mesh,
                           lower_ww_bounds=ww_values,
                           upper_bound_ratio=5.0,
                           energy_bounds=energy_bounds,
                           particle_type=particle_type)
        wws.append(ww)

    return wws
