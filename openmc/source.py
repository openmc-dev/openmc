from collections.abc import Iterable
from enum import Enum
from numbers import Real
import warnings
import typing  # imported separately as py3.8 requires typing.Iterable
# also required to prevent typing.Union namespace overwriting Union
from typing import Optional, Sequence
from xml.etree import ElementTree as ET

import numpy as np
import h5py

import openmc
import openmc.checkvalue as cv
from openmc.checkvalue import PathLike
from openmc.stats.multivariate import UnitSphere, Spatial
from openmc.stats.univariate import Univariate
from ._xml import get_text


class Source:
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
    filename : str
        Source file from which sites should be sampled
    library : str
        Path to a custom source library
    parameters : str
        Parameters to be provided to the custom source library

        .. versionadded:: 0.12
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
    file : str or None
        Source file from which sites should be sampled
    library : str or None
        Path to a custom source library
    parameters : str
        Parameters to be provided to the custom source library
    strength : float
        Strength of the source
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
        filename: Optional[str] = None,
        library: Optional[str] = None,
        parameters: Optional[str] = None,
        strength: float = 1.0,
        particle: str = 'neutron',
        domains: Optional[Sequence[typing.Union[openmc.Cell, openmc.Material, openmc.Universe]]] = None
    ):
        self._space = None
        self._angle = None
        self._energy = None
        self._time = None
        self._file = None
        self._library = None
        self._parameters = None

        if space is not None:
            self.space = space
        if angle is not None:
            self.angle = angle
        if energy is not None:
            self.energy = energy
        if time is not None:
            self.time = time
        if filename is not None:
            self.file = filename
        if library is not None:
            self.library = library
        if parameters is not None:
            self.parameters = parameters
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
    def file(self):
        return self._file

    @property
    def library(self):
        return self._library

    @property
    def parameters(self):
        return self._parameters

    @property
    def space(self):
        return self._space

    @property
    def angle(self):
        return self._angle

    @property
    def energy(self):
        return self._energy

    @property
    def time(self):
        return self._time

    @property
    def strength(self):
        return self._strength

    @property
    def particle(self):
        return self._particle

    @property
    def domain_ids(self):
        return self._domain_ids

    @property
    def domain_type(self):
        return self._domain_type

    @domain_ids.setter
    def domain_ids(self, ids):
        cv.check_type('domain IDs', ids, Iterable, Real)
        self._domain_ids = ids

    @domain_type.setter
    def domain_type(self, domain_type):
        cv.check_value('domain type', domain_type, ('cell', 'material', 'universe'))
        self._domain_type = domain_type

    @file.setter
    def file(self, filename):
        cv.check_type('source file', filename, str)
        self._file = filename

    @library.setter
    def library(self, library_name):
        cv.check_type('library', library_name, str)
        self._library = library_name

    @parameters.setter
    def parameters(self, parameters_path):
        cv.check_type('parameters', parameters_path, str)
        self._parameters = parameters_path

    @space.setter
    def space(self, space):
        cv.check_type('spatial distribution', space, Spatial)
        self._space = space

    @angle.setter
    def angle(self, angle):
        cv.check_type('angular distribution', angle, UnitSphere)
        self._angle = angle

    @energy.setter
    def energy(self, energy):
        cv.check_type('energy distribution', energy, Univariate)
        self._energy = energy

    @time.setter
    def time(self, time):
        cv.check_type('time distribution', time, Univariate)
        self._time = time

    @strength.setter
    def strength(self, strength):
        cv.check_type('source strength', strength, Real)
        cv.check_greater_than('source strength', strength, 0.0, True)
        self._strength = strength

    @particle.setter
    def particle(self, particle):
        cv.check_value('source particle', particle, ['neutron', 'photon'])
        self._particle = particle

    def to_xml_element(self) -> ET.Element:
        """Return XML representation of the source

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing source data

        """
        element = ET.Element("source")
        element.set("strength", str(self.strength))
        if self.particle != 'neutron':
            element.set("particle", self.particle)
        if self.file is not None:
            element.set("file", self.file)
        if self.library is not None:
            element.set("library", self.library)
        if self.parameters is not None:
            element.set("parameters", self.parameters)
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
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element) -> 'openmc.Source':
        """Generate source from an XML element

        Parameters
        ----------
        elem : xml.etree.ElementTree.Element
            XML element

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

        library = get_text(elem, 'library')
        if library is not None:
            source.library = library

        parameters = get_text(elem, 'parameters')
        if parameters is not None:
            source.parameters = parameters

        space = elem.find('space')
        if space is not None:
            source.space = Spatial.from_xml_element(space)

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


class ParticleType(Enum):
    NEUTRON = 0
    PHOTON = 1
    ELECTRON = 2
    POSITRON = 3


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
