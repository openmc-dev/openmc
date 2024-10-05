from __future__ import annotations
from abc import ABC, abstractmethod
from collections.abc import Iterable, Sequence
from math import cos, pi
from numbers import Real
from warnings import warn

import lxml.etree as ET
import numpy as np

import openmc
import openmc.checkvalue as cv
from .._xml import get_text
from ..mesh import MeshBase
from .univariate import PowerLaw, Uniform, Univariate


class UnitSphere(ABC):
    """Distribution of points on the unit sphere.

    This abstract class is used for angular distributions, since a direction is
    represented as a unit vector (i.e., vector on the unit sphere).

    Parameters
    ----------
    reference_uvw : Iterable of float
        Direction from which polar angle is measured

    Attributes
    ----------
    reference_uvw : Iterable of float
        Direction from which polar angle is measured

    """
    def __init__(self, reference_uvw=None):
        self._reference_uvw = None
        if reference_uvw is not None:
            self.reference_uvw = reference_uvw

    @property
    def reference_uvw(self):
        return self._reference_uvw

    @reference_uvw.setter
    def reference_uvw(self, uvw):
        cv.check_type('reference direction', uvw, Iterable, Real)
        uvw = np.asarray(uvw)
        self._reference_uvw = uvw/np.linalg.norm(uvw)

    @abstractmethod
    def to_xml_element(self):
        return ''

    @classmethod
    @abstractmethod
    def from_xml_element(cls, elem):
        distribution = get_text(elem, 'type')
        if distribution == 'mu-phi':
            return PolarAzimuthal.from_xml_element(elem)
        elif distribution == 'isotropic':
            return Isotropic.from_xml_element(elem)
        elif distribution == 'monodirectional':
            return Monodirectional.from_xml_element(elem)


class PolarAzimuthal(UnitSphere):
    """Angular distribution represented by polar and azimuthal angles

    This distribution allows one to specify the distribution of the cosine of
    the polar angle and the azimuthal angle independently of one another. The
    polar angle is measured relative to the reference angle.

    Parameters
    ----------
    mu : openmc.stats.Univariate
        Distribution of the cosine of the polar angle
    phi : openmc.stats.Univariate
        Distribution of the azimuthal angle in radians
    reference_uvw : Iterable of float
        Direction from which polar angle is measured. Defaults to the positive
        z-direction.

    Attributes
    ----------
    mu : openmc.stats.Univariate
        Distribution of the cosine of the polar angle
    phi : openmc.stats.Univariate
        Distribution of the azimuthal angle in radians

    """

    def __init__(self, mu=None, phi=None, reference_uvw=(0., 0., 1.)):
        super().__init__(reference_uvw)
        if mu is not None:
            self.mu = mu
        else:
            self.mu = Uniform(-1., 1.)

        if phi is not None:
            self.phi = phi
        else:
            self.phi = Uniform(0., 2*pi)

    @property
    def mu(self):
        return self._mu

    @mu.setter
    def mu(self, mu):
        cv.check_type('cosine of polar angle', mu, Univariate)
        self._mu = mu

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, phi):
        cv.check_type('azimuthal angle', phi, Univariate)
        self._phi = phi

    def to_xml_element(self):
        """Return XML representation of the angular distribution

        Returns
        -------
        element : lxml.etree._Element
            XML element containing angular distribution data

        """
        element = ET.Element('angle')
        element.set("type", "mu-phi")
        if self.reference_uvw is not None:
            element.set("reference_uvw", ' '.join(map(str, self.reference_uvw)))
        element.append(self.mu.to_xml_element('mu'))
        element.append(self.phi.to_xml_element('phi'))
        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate angular distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.PolarAzimuthal
            Angular distribution generated from XML element

        """
        mu_phi = cls()
        uvw = get_text(elem, 'reference_uvw')
        if uvw is not None:
            mu_phi.reference_uvw = [float(x) for x in uvw.split()]
        mu_phi.mu = Univariate.from_xml_element(elem.find('mu'))
        mu_phi.phi = Univariate.from_xml_element(elem.find('phi'))
        return mu_phi


class Isotropic(UnitSphere):
    """Isotropic angular distribution."""

    def __init__(self):
        super().__init__()

    def to_xml_element(self):
        """Return XML representation of the isotropic distribution

        Returns
        -------
        element : lxml.etree._Element
            XML element containing isotropic distribution data

        """
        element = ET.Element('angle')
        element.set("type", "isotropic")
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate isotropic distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.Isotropic
            Isotropic distribution generated from XML element

        """
        return cls()


class Monodirectional(UnitSphere):
    """Monodirectional angular distribution.

    A monodirectional angular distribution is one for which the polar and
    azimuthal angles are always the same. It is completely specified by the
    reference direction vector.

    Parameters
    ----------
    reference_uvw : Iterable of float
        Direction from which polar angle is measured. Defaults to the positive
        x-direction.

    """

    def __init__(self, reference_uvw: Sequence[float] = [1., 0., 0.]):
        super().__init__(reference_uvw)

    def to_xml_element(self):
        """Return XML representation of the monodirectional distribution

        Returns
        -------
        element : lxml.etree._Element
            XML element containing monodirectional distribution data

        """
        element = ET.Element('angle')
        element.set("type", "monodirectional")
        if self.reference_uvw is not None:
            element.set("reference_uvw", ' '.join(map(str, self.reference_uvw)))
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate monodirectional distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.Monodirectional
            Monodirectional distribution generated from XML element

        """
        monodirectional = cls()
        uvw = get_text(elem, 'reference_uvw')
        if uvw is not None:
            monodirectional.reference_uvw = [float(x) for x in uvw.split()]
        return monodirectional


class Spatial(ABC):
    """Distribution of locations in three-dimensional Euclidean space.

    Classes derived from this abstract class can be used for spatial
    distributions of source sites.

    """
    @abstractmethod
    def to_xml_element(self):
        return ''

    @classmethod
    @abstractmethod
    def from_xml_element(cls, elem, meshes=None):
        distribution = get_text(elem, 'type')
        if distribution == 'cartesian':
            return CartesianIndependent.from_xml_element(elem)
        elif distribution == 'cylindrical':
            return CylindricalIndependent.from_xml_element(elem)
        elif distribution == 'spherical':
            return SphericalIndependent.from_xml_element(elem)
        elif distribution == 'box' or distribution == 'fission':
            return Box.from_xml_element(elem)
        elif distribution == 'point':
            return Point.from_xml_element(elem)
        elif distribution == 'mesh':
            return MeshSpatial.from_xml_element(elem, meshes)
        elif distribution == 'cloud':
            return PointCloud.from_xml_element(elem)


class CartesianIndependent(Spatial):
    """Spatial distribution with independent x, y, and z distributions.

    This distribution allows one to specify coordinates whose x-, y-, and z-
    components are sampled independently from one another.

    Parameters
    ----------
    x : openmc.stats.Univariate
        Distribution of x-coordinates
    y : openmc.stats.Univariate
        Distribution of y-coordinates
    z : openmc.stats.Univariate
        Distribution of z-coordinates

    Attributes
    ----------
    x : openmc.stats.Univariate
        Distribution of x-coordinates
    y : openmc.stats.Univariate
        Distribution of y-coordinates
    z : openmc.stats.Univariate
        Distribution of z-coordinates

    """

    def __init__(
        self,
        x: openmc.stats.Univariate,
        y: openmc.stats.Univariate,
        z: openmc.stats.Univariate
    ):
        self.x = x
        self.y = y
        self.z = z

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        cv.check_type('x coordinate', x, Univariate)
        self._x = x

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        cv.check_type('y coordinate', y, Univariate)
        self._y = y

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, z):
        cv.check_type('z coordinate', z, Univariate)
        self._z = z

    def to_xml_element(self):
        """Return XML representation of the spatial distribution

        Returns
        -------
        element : lxml.etree._Element
            XML element containing spatial distribution data

        """
        element = ET.Element('space')
        element.set('type', 'cartesian')
        element.append(self.x.to_xml_element('x'))
        element.append(self.y.to_xml_element('y'))
        element.append(self.z.to_xml_element('z'))
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate spatial distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.CartesianIndependent
            Spatial distribution generated from XML element

        """
        x = Univariate.from_xml_element(elem.find('x'))
        y = Univariate.from_xml_element(elem.find('y'))
        z = Univariate.from_xml_element(elem.find('z'))
        return cls(x, y, z)


class SphericalIndependent(Spatial):
    r"""Spatial distribution represented in spherical coordinates.

    This distribution allows one to specify coordinates whose :math:`r`,
    :math:`\theta`, and :math:`\phi` components are sampled independently
    from one another and centered on the coordinates (x0, y0, z0).

    .. versionadded:: 0.12

    .. versionchanged:: 0.13.1
        Accepts ``cos_theta`` instead of ``theta``

    Parameters
    ----------
    r : openmc.stats.Univariate
        Distribution of r-coordinates in a reference frame specified by
        the origin parameter
    cos_theta : openmc.stats.Univariate
        Distribution of the cosine of the theta-coordinates (angle relative to
        the z-axis) in a reference frame specified by the origin parameter
    phi : openmc.stats.Univariate
        Distribution of phi-coordinates (azimuthal angle) in a reference frame
        specified by the origin parameter
    origin: Iterable of float, optional
        coordinates (x0, y0, z0) of the center of the spherical reference frame
        for the source. Defaults to (0.0, 0.0, 0.0)

    Attributes
    ----------
    r : openmc.stats.Univariate
        Distribution of r-coordinates in the local reference frame
    cos_theta : openmc.stats.Univariate
        Distribution of the cosine of the theta-coordinates (angle relative to
        the z-axis) in the local reference frame
    phi : openmc.stats.Univariate
        Distribution of phi-coordinates (azimuthal angle) in the local
        reference frame
    origin: Iterable of float, optional
        coordinates (x0, y0, z0) of the center of the spherical reference
        frame. Defaults to (0.0, 0.0, 0.0)

    """

    def __init__(self, r, cos_theta, phi, origin=(0.0, 0.0, 0.0)):
        self.r = r
        self.cos_theta = cos_theta
        self.phi = phi
        self.origin = origin

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, r):
        cv.check_type('r coordinate', r, Univariate)
        self._r = r

    @property
    def cos_theta(self):
        return self._cos_theta

    @cos_theta.setter
    def cos_theta(self, cos_theta):
        cv.check_type('cos_theta coordinate', cos_theta, Univariate)
        self._cos_theta = cos_theta

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, phi):
        cv.check_type('phi coordinate', phi, Univariate)
        self._phi = phi

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, origin):
        cv.check_type('origin coordinates', origin, Iterable, Real)
        origin = np.asarray(origin)
        self._origin = origin

    def to_xml_element(self):
        """Return XML representation of the spatial distribution

        Returns
        -------
        element : lxml.etree._Element
            XML element containing spatial distribution data

        """
        element = ET.Element('space')
        element.set('type', 'spherical')
        element.append(self.r.to_xml_element('r'))
        element.append(self.cos_theta.to_xml_element('cos_theta'))
        element.append(self.phi.to_xml_element('phi'))
        element.set("origin", ' '.join(map(str, self.origin)))
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate spatial distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.SphericalIndependent
            Spatial distribution generated from XML element

        """
        r = Univariate.from_xml_element(elem.find('r'))
        cos_theta = Univariate.from_xml_element(elem.find('cos_theta'))
        phi = Univariate.from_xml_element(elem.find('phi'))
        origin = [float(x) for x in elem.get('origin').split()]
        return cls(r, cos_theta, phi, origin=origin)


class CylindricalIndependent(Spatial):
    r"""Spatial distribution represented in cylindrical coordinates.

    This distribution allows one to specify coordinates whose :math:`r`,
    :math:`\phi`, and :math:`z` components are sampled independently from
    one another and in a reference frame whose origin is specified by the
    coordinates (x0, y0, z0).

    .. versionadded:: 0.12

    Parameters
    ----------
    r : openmc.stats.Univariate
        Distribution of r-coordinates in a reference frame specified by the
        origin parameter
    phi : openmc.stats.Univariate
        Distribution of phi-coordinates (azimuthal angle) in a reference frame
        specified by the origin parameter
    z : openmc.stats.Univariate
        Distribution of z-coordinates in a reference frame specified by the
        origin parameter
    origin: Iterable of float, optional
        coordinates (x0, y0, z0) of the center of the cylindrical reference
        frame. Defaults to (0.0, 0.0, 0.0)

    Attributes
    ----------
    r : openmc.stats.Univariate
        Distribution of r-coordinates in the local reference frame
    phi : openmc.stats.Univariate
        Distribution of phi-coordinates (azimuthal angle) in the local
        reference frame
    z : openmc.stats.Univariate
        Distribution of z-coordinates in the local reference frame
    origin: Iterable of float, optional
        coordinates (x0, y0, z0) of the center of the cylindrical reference
        frame. Defaults to (0.0, 0.0, 0.0)

    """

    def __init__(self, r, phi, z, origin=(0.0, 0.0, 0.0)):
        self.r = r
        self.phi = phi
        self.z = z
        self.origin = origin

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, r):
        cv.check_type('r coordinate', r, Univariate)
        self._r = r

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, phi):
        cv.check_type('phi coordinate', phi, Univariate)
        self._phi = phi

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, z):
        cv.check_type('z coordinate', z, Univariate)
        self._z = z

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, origin):
        cv.check_type('origin coordinates', origin, Iterable, Real)
        origin = np.asarray(origin)
        self._origin = origin

    def to_xml_element(self):
        """Return XML representation of the spatial distribution

        Returns
        -------
        element : lxml.etree._Element
            XML element containing spatial distribution data

        """
        element = ET.Element('space')
        element.set('type', 'cylindrical')
        element.append(self.r.to_xml_element('r'))
        element.append(self.phi.to_xml_element('phi'))
        element.append(self.z.to_xml_element('z'))
        element.set("origin", ' '.join(map(str, self.origin)))
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate spatial distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.CylindricalIndependent
            Spatial distribution generated from XML element

        """
        r = Univariate.from_xml_element(elem.find('r'))
        phi = Univariate.from_xml_element(elem.find('phi'))
        z = Univariate.from_xml_element(elem.find('z'))
        origin = [float(x) for x in elem.get('origin').split()]
        return cls(r, phi, z, origin=origin)


class MeshSpatial(Spatial):
    """Spatial distribution for a mesh.

    This distribution specifies a mesh to sample over with source strengths
    specified for each mesh element.

    .. versionadded:: 0.13.3

    Parameters
    ----------
    mesh : openmc.MeshBase
        The mesh instance used for sampling
    strengths : iterable of float, optional
        An iterable of values that represents the weights of each element. If no
        source strengths are specified, they will be equal for all mesh
        elements.
    volume_normalized : bool, optional
        Whether or not the strengths will be multiplied by element volumes at
        runtime. Default is True.

    Attributes
    ----------
    mesh : openmc.MeshBase
        The mesh instance used for sampling
    strengths : numpy.ndarray or None
        An array of source strengths for each mesh element
    volume_normalized : bool
        Whether or not the strengths will be multiplied by element volumes at
        runtime.
    """

    def __init__(self, mesh, strengths=None, volume_normalized=True):
        self.mesh = mesh
        self.strengths = strengths
        self.volume_normalized = volume_normalized

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        if mesh is not None:
            cv.check_type('mesh instance', mesh, MeshBase)
        self._mesh = mesh

    @property
    def volume_normalized(self):
        return self._volume_normalized

    @volume_normalized.setter
    def volume_normalized(self, volume_normalized):
        cv.check_type('Multiply strengths by element volumes', volume_normalized, bool)
        self._volume_normalized = volume_normalized

    @property
    def strengths(self):
        return self._strengths

    @strengths.setter
    def strengths(self, given_strengths):
        if given_strengths is not None:
            cv.check_type('strengths array passed in', given_strengths, Iterable, Real)
            self._strengths = np.asarray(given_strengths, dtype=float).flatten()
        else:
            self._strengths = None

    @property
    def num_strength_bins(self):
        if self.strengths is None:
            raise ValueError('Strengths are not set')
        return self.strengths.size

    def to_xml_element(self):
        """Return XML representation of the spatial distribution

        Returns
        -------
        element : lxml.etree._Element
            XML element containing spatial distribution data

        """
        element = ET.Element('space')

        element.set('type', 'mesh')
        element.set("mesh_id", str(self.mesh.id))
        element.set("volume_normalized", str(self.volume_normalized))

        if self.strengths is not None:
            subelement = ET.SubElement(element, 'strengths')
            subelement.text = ' '.join(str(e) for e in self.strengths)

        return element

    @classmethod
    def from_xml_element(cls, elem, meshes):
        """Generate spatial distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element
        meshes : dict
            A dictionary with mesh IDs as keys and openmc.MeshBase instances as
            values

        Returns
        -------
        openmc.stats.MeshSpatial
            Spatial distribution generated from XML element

        """

        mesh_id = int(elem.get('mesh_id'))

        # check if this mesh has been read in from another location already
        if mesh_id not in meshes:
            raise ValueError(f'Could not locate mesh with ID "{mesh_id}"')

        volume_normalized = elem.get("volume_normalized")
        volume_normalized = get_text(elem, 'volume_normalized').lower() == 'true'
        strengths = get_text(elem, 'strengths')
        if strengths is not None:
            strengths = [float(b) for b in get_text(elem, 'strengths').split()]

        return cls(meshes[mesh_id], strengths, volume_normalized)


class PointCloud(Spatial):
    """Spatial distribution from a point cloud.

    This distribution specifies a discrete list of points, 
    each with different relative probability.

    .. versionadded:: 0.15.x

    Parameters
    ----------
    positions : iterable of 3-tuples
        The points in space to be sampled
    strengths : iterable of float, optional
        An iterable of values that represents the relative probabilty of each point.

    Attributes
    ----------
    psoitions: numpy.ndarray (3xN)
        The points in space to be sampled
    strengths : numpy.ndarray or None
        An array of relative probabilities for each mesh point
    """

    def __init__(self, positions, strengths=None):
        self.positions = positions
        self.strengths = strengths

    @property
    def positions(self):
        return self._positions

    @positions.setter
    def positions(self, given_positions):
        if given_positions is None:
            raise ValueError('No positions were provided')
        cv.check_iterable_type('position list passed in', given_positions, Real, 2, 2)

        if isinstance(given_positions, list):
            cv.check_length('first position entry', given_positions[0], 3, 3)
            self._positions = np.asarray(given_positions)
        elif isinstance(given_positions, np.ndarray):
            self._positions = given_positions
        else:
            raise ValueError('Unable to interpret list of positions')

    @property
    def strengths(self):
        return self._strengths

    @strengths.setter
    def strengths(self, given_strengths):
        if given_strengths is not None:
            cv.check_type('strengths array passed in', given_strengths, Iterable, Real)
            self._strengths = np.asarray(given_strengths, dtype=float).flatten()
        else:
            self._strengths = None

    @property
    def num_strength_bins(self):
        if self.strengths is None:
            raise ValueError('Strengths are not set')
        return self.strengths.size

    def to_xml_element(self):
        """Return XML representation of the spatial distribution

        Returns
        -------
        element : lxml.etree._Element
            XML element containing spatial distribution data

        """
        element = ET.Element('space')

        element.set('type', 'cloud')

        for idx, axis in enumerate(('x','y','z')):
            subelement = ET.SubElement(element, axis)
            subelement.text = ' '.join(str(e) for e in self.positions[idx,:])

        if self.strengths is not None:
            subelement = ET.SubElement(element, 'strengths')
            subelement.text = ' '.join(str(e) for e in self.strengths)

        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate spatial distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.PointCloud
            Spatial distribution generated from XML element


        """
        coord = {}

        for axis in enumerate(('x','y','z')):
            coord_data = get_text(elem, axis)
            if coord_data is not None:
                coord[axis] = [float(b) for b in coord_data.split]
        
        positions = np.column_stack([coord[axis] for axis in ('x','y','z')])

        strengths = get_text(elem, 'strengths')
        if strengths is not None:
            strengths = [float(b) for b in strengths.split()]

        return cls(positions, strengths)



class Box(Spatial):
    """Uniform distribution of coordinates in a rectangular cuboid.

    Parameters
    ----------
    lower_left : Iterable of float
        Lower-left coordinates of cuboid
    upper_right : Iterable of float
        Upper-right coordinates of cuboid
    only_fissionable : bool, optional
        Whether spatial sites should only be accepted if they occur in
        fissionable materials

        .. deprecated:: 0.15.0
            Use the `constraints` argument when defining a source object instead.

    Attributes
    ----------
    lower_left : Iterable of float
        Lower-left coordinates of cuboid
    upper_right : Iterable of float
        Upper-right coordinates of cuboid
    only_fissionable : bool, optional
        Whether spatial sites should only be accepted if they occur in
        fissionable materials

        .. deprecated:: 0.15.0
            Use the `constraints` argument when defining a source object instead.

    """

    def __init__(
        self,
        lower_left: Sequence[float],
        upper_right: Sequence[float],
        only_fissionable: bool = False
    ):
        self.lower_left = lower_left
        self.upper_right = upper_right
        self.only_fissionable = only_fissionable

    @property
    def lower_left(self):
        return self._lower_left

    @lower_left.setter
    def lower_left(self, lower_left):
        cv.check_type('lower left coordinate', lower_left, Iterable, Real)
        cv.check_length('lower left coordinate', lower_left, 3)
        self._lower_left = lower_left

    @property
    def upper_right(self):
        return self._upper_right

    @upper_right.setter
    def upper_right(self, upper_right):
        cv.check_type('upper right coordinate', upper_right, Iterable, Real)
        cv.check_length('upper right coordinate', upper_right, 3)
        self._upper_right = upper_right

    @property
    def only_fissionable(self):
        return self._only_fissionable

    @only_fissionable.setter
    def only_fissionable(self, only_fissionable):
        cv.check_type('only fissionable', only_fissionable, bool)
        self._only_fissionable = only_fissionable
        if only_fissionable:
            warn("The 'only_fissionable' has been deprecated. Use the "
                 "'constraints' argument when defining a source instead.",
                 FutureWarning)

    def to_xml_element(self):
        """Return XML representation of the box distribution

        Returns
        -------
        element : lxml.etree._Element
            XML element containing box distribution data

        """
        element = ET.Element('space')
        if self.only_fissionable:
            element.set("type", "fission")
        else:
            element.set("type", "box")
        params = ET.SubElement(element, "parameters")
        params.text = ' '.join(map(str, self.lower_left)) + ' ' + \
                      ' '.join(map(str, self.upper_right))
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate box distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.Box
            Box distribution generated from XML element

        """
        only_fissionable = get_text(elem, 'type') == 'fission'
        params = [float(x) for x in get_text(elem, 'parameters').split()]
        lower_left = params[:len(params)//2]
        upper_right = params[len(params)//2:]
        return cls(lower_left, upper_right, only_fissionable)


class Point(Spatial):
    """Delta function in three dimensions.

    This spatial distribution can be used for a point source where sites are
    emitted at a specific location given by its Cartesian coordinates.

    Parameters
    ----------
    xyz : Iterable of float, optional
        Cartesian coordinates of location. Defaults to (0., 0., 0.).

    Attributes
    ----------
    xyz : Iterable of float
        Cartesian coordinates of location

    """

    def __init__(self, xyz: Sequence[float] = (0., 0., 0.)):
        self.xyz = xyz

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, xyz):
        cv.check_type('coordinate', xyz, Iterable, Real)
        cv.check_length('coordinate', xyz, 3)
        self._xyz = xyz

    def to_xml_element(self):
        """Return XML representation of the point distribution

        Returns
        -------
        element : lxml.etree._Element
            XML element containing point distribution location

        """
        element = ET.Element('space')
        element.set("type", "point")
        params = ET.SubElement(element, "parameters")
        params.text = ' '.join(map(str, self.xyz))
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate point distribution from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.stats.Point
            Point distribution generated from XML element

        """
        xyz = [float(x) for x in get_text(elem, 'parameters').split()]
        return cls(xyz)


def spherical_uniform(
        r_outer: float,
        r_inner: float = 0.0,
        thetas: Sequence[float] = (0., pi),
        phis: Sequence[float] = (0., 2*pi),
        origin: Sequence[float] = (0., 0., 0.)
    ):
    """Return a uniform spatial distribution over a spherical shell.

    This function provides a uniform spatial distribution over a spherical
    shell between `r_inner` and `r_outer`. Optionally, the range of angles
    can be restricted by the `thetas` and `phis` arguments.

    .. versionadded:: 0.13.1

    Parameters
    ----------
    r_outer : float
        Outer radius of the spherical shell in [cm]
    r_inner : float
        Inner radius of the spherical shell in [cm]
    thetas : iterable of float
        Starting and ending theta coordinates (angle relative to
        the z-axis) in radius in a reference frame centered at `origin`
    phis : iterable of float
        Starting and ending phi coordinates (azimuthal angle) in
        radians in a reference frame centered at `origin`
    origin: iterable of float
        Coordinates (x0, y0, z0) of the center of the spherical
        reference frame for the distribution.

    Returns
    -------
    openmc.stats.SphericalIndependent
        Uniform distribution over the spherical shell
    """

    r_dist = PowerLaw(r_inner, r_outer, 2)
    cos_thetas_dist = Uniform(cos(thetas[1]), cos(thetas[0]))
    phis_dist = Uniform(phis[0], phis[1])

    return SphericalIndependent(r_dist, cos_thetas_dist, phis_dist, origin)
