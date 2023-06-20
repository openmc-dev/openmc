from __future__ import annotations
from abc import ABC, abstractmethod
from collections.abc import Iterable
from enum import IntEnum
from numbers import Real
import warnings
import typing  # imported separately as py3.8 requires typing.Iterable
# also required to prevent typing.Union namespace overwriting Union
from typing import Optional, Sequence
import lxml.etree as ET

import numpy as np
import h5py

import openmc
import openmc.checkvalue as cv
from openmc.checkvalue import PathLike
from openmc.stats.multivariate import UnitSphere, Spatial
from openmc.stats.univariate import Univariate
from ._xml import get_text
from .mesh import MeshBase


class SourceBase(ABC):
    """Base class for external sources

    Parameters
    ----------
    strength : float
        Strength of the source

    Attributes
    ----------
    type : {'independent', 'file', 'compiled'}
        Indicator of source type.
    strength : float
        Strength of the source

    """

    def __init__(self, strength=1.0):
        self.strength = strength

    @property
    def strength(self):
        return self._strength

    @strength.setter
    def strength(self, strength):
        cv.check_type('source strength', strength, Real)
        cv.check_greater_than('source strength', strength, 0.0, True)
        self._strength = strength

    @abstractmethod
    def populate_xml_element(self, element):
        """Add necessary source information to an XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """
        pass

    def to_xml_element(self) -> ET.Element:
        """Return XML representation of the source

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing source data

        """
        element = ET.Element("source")
        element.set("type", self.type)
        element.set("strength", str(self.strength))
        self.populate_xml_element(element)
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element, meshes=None) -> openmc.SourceBase:
        """Generate source from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element
        meshes : dict
            Dictionary with mesh IDs as keys and openmc.MeshBase instances as
            values

        Returns
        -------
        openmc.SourceBase
            Source generated from XML element

        """
        source_type = get_text(elem, 'type')

        if source_type is None:
            # attempt to determine source type based on attributes
            # for backward compatibility
            if get_text(elem, 'file') is not None:
                return FileSource.from_xml_element(elem)
            elif get_text(elem, 'library') is not None:
                return CompiledSource.from_xml_element(elem)
            else:
                return IndependentSource.from_xml_element(elem)
        else:
            if source_type == 'independent':
                return IndependentSource.from_xml_element(elem, meshes)
            elif source_type == 'compiled':
                return CompiledSource.from_xml_element(elem)
            elif source_type == 'file':
                return FileSource.from_xml_element(elem)
            else:
                raise ValueError(f'Source type {source_type} is not recognized')


class IndependentSource(SourceBase):
    """Distribution of phase space coordinates for source sites.

    Parameters
    ----------
    space : openmc.stats.Spatial
        Spatial distribution of source sites
    angle : openmc.stats.UnitSphere
        Angular distribution of source sites
    energy : openmc.stats.Univariate
        Energy distribution of source sites
    time : openmc.stats.Univariate
        time distribution of source sites
    strength : float
        Strength of the source
    particle : {'neutron', 'photon'}
        Source particle type
    domains : iterable of openmc.Cell, openmc.Material, or openmc.Universe
        Domains to reject based on, i.e., if a sampled spatial location is not
        within one of these domains, it will be rejected.

    Attributes
    ----------
    space : openmc.stats.Spatial or None
        Spatial distribution of source sites
    angle : openmc.stats.UnitSphere or None
        Angular distribution of source sites
    energy : openmc.stats.Univariate or None
        Energy distribution of source sites
    time : openmc.stats.Univariate or None
        time distribution of source sites
    strength : float
        Strength of the source
    type : str
        Indicator of source type: 'independent'

    .. versionadded:: 0.13.4

    particle : {'neutron', 'photon'}
        Source particle type
    ids : Iterable of int
        IDs of domains to use for rejection
    domain_type : {'cell', 'material', 'universe'}
        Type of domain to use for rejection

    """

    def __init__(
        self,
        space: Optional[openmc.stats.Spatial] = None,
        angle: Optional[openmc.stats.UnitSphere] = None,
        energy: Optional[openmc.stats.Univariate] = None,
        time: Optional[openmc.stats.Univariate] = None,
        strength: float = 1.0,
        particle: str = 'neutron',
        domains: Optional[Sequence[typing.Union[openmc.Cell, openmc.Material, openmc.Universe]]] = None
    ):
        super().__init__(strength)

        self._space = None
        self._angle = None
        self._energy = None
        self._time = None

        if space is not None:
            self.space = space
        if angle is not None:
            self.angle = angle
        if energy is not None:
            self.energy = energy
        if time is not None:
            self.time = time
        self.strength = strength
        self.particle = particle

        self._domain_ids = []
        self._domain_type = None
        if domains is not None:
            if isinstance(domains[0], openmc.Cell):
                self.domain_type = 'cell'
            elif isinstance(domains[0], openmc.Material):
                self.domain_type = 'material'
            elif isinstance(domains[0], openmc.Universe):
                self.domain_type = 'universe'
            self.domain_ids = [d.id for d in domains]

    @property
    def type(self) -> str:
        return 'independent'

    def __getattr__(self, name):
        cls_names = {'file': 'FileSource', 'library': 'CompiledSource',
                     'parameters': 'CompiledSource'}
        if name in cls_names:
            raise AttributeError(
                f'The "{name}" attribute has been deprecated on the '
                f'IndependentSource class. Please use the {cls_names[name]} class.')
        else:
            super().__getattribute__(name)

    def __setattr__(self, name, value):
        if name in ('file', 'library', 'parameters'):
            # Ensure proper AttributeError is thrown
            getattr(self, name)
        else:
            super().__setattr__(name, value)

    @property
    def space(self):
        return self._space

    @space.setter
    def space(self, space):
        cv.check_type('spatial distribution', space, Spatial)
        self._space = space

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, angle):
        cv.check_type('angular distribution', angle, UnitSphere)
        self._angle = angle

    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self, energy):
        cv.check_type('energy distribution', energy, Univariate)
        self._energy = energy

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, time):
        cv.check_type('time distribution', time, Univariate)
        self._time = time

    @property
    def particle(self):
        return self._particle

    @particle.setter
    def particle(self, particle):
        cv.check_value('source particle', particle, ['neutron', 'photon'])
        self._particle = particle

    @property
    def domain_ids(self):
        return self._domain_ids

    @domain_ids.setter
    def domain_ids(self, ids):
        cv.check_type('domain IDs', ids, Iterable, Real)
        self._domain_ids = ids

    @property
    def domain_type(self):
        return self._domain_type

    @domain_type.setter
    def domain_type(self, domain_type):
        cv.check_value('domain type', domain_type, ('cell', 'material', 'universe'))
        self._domain_type = domain_type

    def populate_xml_element(self, element):
        """Add necessary source information to an XML element
        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """
        super().populate_xml_element(element)
        element.set("particle", self.particle)
        if self.space is not None:
            element.append(self.space.to_xml_element())
        if self.angle is not None:
            element.append(self.angle.to_xml_element())
        if self.energy is not None:
            element.append(self.energy.to_xml_element('energy'))
        if self.time is not None:
            element.append(self.time.to_xml_element('time'))
        if self.domain_ids:
            dt_elem = ET.SubElement(element, "domain_type")
            dt_elem.text = self.domain_type
            id_elem = ET.SubElement(element, "domain_ids")
            id_elem.text = ' '.join(str(uid) for uid in self.domain_ids)

    @classmethod
    def from_xml_element(cls, elem: ET.Element, meshes=None) -> 'openmc.SourceBase':
        """Generate source from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element
        meshes : dict
            Dictionary with mesh IDs as keys and openmc.MeshBase instaces as
            values

        Returns
        -------
        openmc.Source
            Source generated from XML element

        """
        domain_type = get_text(elem, "domain_type")
        if domain_type is not None:
            domain_ids = [int(x) for x in get_text(elem, "domain_ids").split()]

            # Instantiate some throw-away domains that are used by the
            # constructor to assign IDs
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', openmc.IDWarning)
                if domain_type == 'cell':
                    domains = [openmc.Cell(uid) for uid in domain_ids]
                elif domain_type == 'material':
                    domains = [openmc.Material(uid) for uid in domain_ids]
                elif domain_type == 'universe':
                    domains = [openmc.Universe(uid) for uid in domain_ids]
        else:
            domains = None

        source = cls(domains=domains)

        strength = get_text(elem, 'strength')
        if strength is not None:
            source.strength = float(strength)

        particle = get_text(elem, 'particle')
        if particle is not None:
            source.particle = particle

        filename = get_text(elem, 'file')
        if filename is not None:
            source.file = filename

        space = elem.find('space')
        if space is not None:
            source.space = Spatial.from_xml_element(space, meshes)

        angle = elem.find('angle')
        if angle is not None:
            source.angle = UnitSphere.from_xml_element(angle)

        energy = elem.find('energy')
        if energy is not None:
            source.energy = Univariate.from_xml_element(energy)

        time = elem.find('time')
        if time is not None:
            source.time = Univariate.from_xml_element(time)

        return source


def Source(*args, **kwargs):
    """
    A function for backward compatibility of sources. Will be removed in the
    future. Please update to IndependentSource.
    """
    warnings.warn("This class is deprecated in favor of 'IndependentSource'", FutureWarning)
    return openmc.IndependentSource(*args, **kwargs)


class CompiledSource(SourceBase):
    """A source based on a compiled shared library

    .. versionadded:: 0.13.4

    Parameters
    ----------
    library : str or None
        Path to a compiled shared library
    parameters : str
        Parameters to be provided to the compiled shared library function
    strength : float
        Strength of the source

    Attributes
    ----------
    library : str or None
        Path to a compiled shared library
    parameters : str
        Parameters to be provided to the compiled shared library function
    strength : float
        Strength of the source
    type : str
        Indicator of source type: 'compiled'

    """
    def __init__(self, library: Optional[str] = None, parameters: Optional[str] = None, strength=1.0) -> None:
        super().__init__(strength=strength)

        self._library = None
        if library is not None:
            self.library = library

        self._parameters = None
        if parameters is not None:
            self.parameters = parameters

    @property
    def type(self) -> str:
        return "compiled"

    @property
    def library(self) -> str:
        return self._library

    @library.setter
    def library(self, library_name):
        cv.check_type('library', library_name, str)
        self._library = library_name

    @property
    def parameters(self) -> str:
        return self._parameters

    @parameters.setter
    def parameters(self, parameters_path):
        cv.check_type('parameters', parameters_path, str)
        self._parameters = parameters_path

    def populate_xml_element(self, element):
        """Add necessary compiled source information to an XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """
        super().populate_xml_element(element)

        element.set("library", self.library)

        if self.parameters is not None:
            element.set("parameters", self.parameters)

    @classmethod
    def from_xml_element(cls, elem: ET.Element, meshes=None) -> openmc.CompiledSource:
        """Generate a compiled source from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element
        meshes : dict
            Dictionary with mesh IDs as keys and openmc.MeshBase instances as
            values

        Returns
        -------
        openmc.CompiledSource
            Source generated from XML element

        """
        library = get_text(elem, 'library')

        source = cls(library)

        strength = get_text(elem, 'strength')
        if strength is not None:
            source.strength = float(strength)

        parameters = get_text(elem, 'parameters')
        if parameters is not None:
            source.parameters = parameters

        return source


class FileSource(SourceBase):
    """A source based on particles stored in a file

    .. versionadded:: 0.13.4

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the source file from which sites should be sampled
    strength : float
        Strength of the source (default is 1.0)

    Attributes
    ----------
    path : Pathlike
        Source file from which sites should be sampled
    strength : float
        Strength of the source
    type : str
        Indicator of source type: 'file'

    """
    def __init__(self, path: Optional[PathLike] = None, strength=1.0) -> None:
        super().__init__(strength=strength)

        self._path = None

        if path is not None:
            self.path = path

    @property
    def type(self) -> str:
        return "file"

    @property
    def path(self) -> PathLike:
        return self._path

    @path.setter
    def path(self, p: PathLike):
        cv.check_type('source file', p, str)
        self._path = p

    def populate_xml_element(self, element):
        """Add necessary file source information to an XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """
        super().populate_xml_element(element)

        if self.path is not None:
            element.set("file", self.path)

    @classmethod
    def from_xml_element(cls, elem: ET.Element) -> openmc.FileSource:
        """Generate file source from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element
        meshes : dict
            Dictionary with mesh IDs as keys and openmc.MeshBase instances as
            values

        Returns
        -------
        openmc.FileSource
            Source generated from XML element

        """

        filename = get_text(elem, 'file')

        source = cls(filename)

        strength = get_text(elem, 'strength')
        if strength is not None:
            source.strength = float(strength)

        return source


class ParticleType(IntEnum):
    """
    IntEnum class representing a particle type. Type
    values mirror those found in the C++ class.
    """
    NEUTRON = 0
    PHOTON = 1
    ELECTRON = 2
    POSITRON = 3

    @classmethod
    def from_string(cls, value: str):
        """
        Constructs a ParticleType instance from a string.

        Parameters
        ----------
        value : str
            The string representation of the particle type.

        Returns
        -------
        The corresponding ParticleType instance.
        """
        try:
            return cls[value.upper()]
        except KeyError:
            raise ValueError(f"Invalid string for creation of {cls.__name__}: {value}")

    def __repr__(self) -> str:
        """
        Returns a string representation of the ParticleType instance.

        Returns:
            str: The lowercase name of the ParticleType instance.
        """
        return self.name.lower()

    # needed for < Python 3.11
    def __str__(self) -> str:
        return self.__repr__()

    # needed for <= 3.7, IntEnum will use the mixed-in type's `__format__` method otherwise
    # this forces it to default to the standard object format, relying on __str__ under the hood
    def __format__(self, spec):
        return object.__format__(self, spec)


class SourceParticle:
    """Source particle

    This class can be used to create source particles that can be written to a
    file and used by OpenMC

    Parameters
    ----------
    r : iterable of float
        Position of particle in Cartesian coordinates
    u : iterable of float
        Directional cosines
    E : float
        Energy of particle in [eV]
    time : float
        Time of particle in [s]
    wgt : float
        Weight of the particle
    delayed_group : int
        Delayed group particle was created in (neutrons only)
    surf_id : int
        Surface ID where particle is at, if any.
    particle : ParticleType
        Type of the particle

    """
    def __init__(
        self,
        r: typing.Iterable[float] = (0., 0., 0.),
        u: typing.Iterable[float] = (0., 0., 1.),
        E: float = 1.0e6,
        time: float = 0.0,
        wgt: float = 1.0,
        delayed_group: int = 0,
        surf_id: int = 0,
        particle: ParticleType = ParticleType.NEUTRON
    ):

        self.r = tuple(r)
        self.u = tuple(u)
        self.E = float(E)
        self.time = float(time)
        self.wgt = float(wgt)
        self.delayed_group = delayed_group
        self.surf_id = surf_id
        self.particle = particle

    def __repr__(self):
        name = self.particle.name.lower()
        return f'<SourceParticle: {name} at E={self.E:.6e} eV>'

    def to_tuple(self) -> tuple:
        """Return source particle attributes as a tuple

        Returns
        -------
        tuple
            Source particle attributes

        """
        return (self.r, self.u, self.E, self.time, self.wgt,
                self.delayed_group, self.surf_id, self.particle.value)


def write_source_file(
    source_particles: typing.Iterable[SourceParticle],
    filename: PathLike, **kwargs
):
    """Write a source file using a collection of source particles

    Parameters
    ----------
    source_particles : iterable of SourceParticle
        Source particles to write to file
    filename : str or path-like
        Path to source file to write
    **kwargs
        Keyword arguments to pass to :class:`h5py.File`

    See Also
    --------
    openmc.SourceParticle

    """
    # Create compound datatype for source particles
    pos_dtype = np.dtype([('x', '<f8'), ('y', '<f8'), ('z', '<f8')])
    source_dtype = np.dtype([
        ('r', pos_dtype),
        ('u', pos_dtype),
        ('E', '<f8'),
        ('time', '<f8'),
        ('wgt', '<f8'),
        ('delayed_group', '<i4'),
        ('surf_id', '<i4'),
        ('particle', '<i4'),
    ])

    # Create array of source particles
    cv.check_iterable_type("source particles", source_particles, SourceParticle)
    arr = np.array([s.to_tuple() for s in source_particles], dtype=source_dtype)

    # Write array to file
    kwargs.setdefault('mode', 'w')
    with h5py.File(filename, **kwargs) as fh:
        fh.attrs['filetype'] = np.string_("source")
        fh.create_dataset('source_bank', data=arr, dtype=source_dtype)
