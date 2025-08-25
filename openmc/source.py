from __future__ import annotations
from abc import ABC, abstractmethod
from collections.abc import Iterable, Sequence
from enum import IntEnum
from numbers import Real
from pathlib import Path
import warnings
from typing import Any
from pathlib import Path

import lxml.etree as ET
import numpy as np
import h5py
import pandas as pd

import openmc
import openmc.checkvalue as cv
from openmc.checkvalue import PathLike
from openmc.stats.multivariate import UnitSphere, Spatial
from openmc.stats.univariate import Univariate
from ._xml import get_elem_list, get_text
from .mesh import MeshBase, StructuredMesh, UnstructuredMesh
from .utility_funcs import input_path


class SourceBase(ABC):
    """Base class for external sources

    Parameters
    ----------
    strength : float
        Strength of the source
    constraints : dict
        Constraints on sampled source particles. Valid keys include 'domains',
        'time_bounds', 'energy_bounds', 'fissionable', and 'rejection_strategy'.
        For 'domains', the corresponding value is an iterable of
        :class:`openmc.Cell`, :class:`openmc.Material`, or
        :class:`openmc.Universe` for which sampled sites must be within. For
        'time_bounds' and 'energy_bounds', the corresponding value is a sequence
        of floats giving the lower and upper bounds on time in [s] or energy in
        [eV] that the sampled particle must be within. For 'fissionable', the
        value is a bool indicating that only sites in fissionable material
        should be accepted. The 'rejection_strategy' indicates what should
        happen when a source particle is rejected: either 'resample' (pick a new
        particle) or 'kill' (accept and terminate).

    Attributes
    ----------
    type : {'independent', 'file', 'compiled', 'mesh'}
        Indicator of source type.
    strength : float
        Strength of the source
    constraints : dict
        Constraints on sampled source particles. Valid keys include
        'domain_type', 'domain_ids', 'time_bounds', 'energy_bounds',
        'fissionable', and 'rejection_strategy'.

    """

    def __init__(
        self,
        strength: float | None = 1.0,
        constraints: dict[str, Any] | None = None
    ):
        self.strength = strength
        self.constraints = constraints

    @property
    def strength(self):
        return self._strength

    @strength.setter
    def strength(self, strength):
        cv.check_type('source strength', strength, Real, none_ok=True)
        if strength is not None:
            cv.check_greater_than('source strength', strength, 0.0, True)
        self._strength = strength

    @property
    def constraints(self) -> dict[str, Any]:
        return self._constraints

    @constraints.setter
    def constraints(self, constraints: dict[str, Any] | None):
        self._constraints = {}
        if constraints is None:
            return

        for key, value in constraints.items():
            if key == 'domains':
                cv.check_type('domains', value, Iterable,
                              (openmc.Cell, openmc.Material, openmc.Universe))
                if isinstance(value[0], openmc.Cell):
                    self._constraints['domain_type'] = 'cell'
                elif isinstance(value[0], openmc.Material):
                    self._constraints['domain_type'] = 'material'
                elif isinstance(value[0], openmc.Universe):
                    self._constraints['domain_type'] = 'universe'
                self._constraints['domain_ids'] = [d.id for d in value]
            elif key == 'time_bounds':
                cv.check_type('time bounds', value, Iterable, Real)
                self._constraints['time_bounds'] = tuple(value)
            elif key == 'energy_bounds':
                cv.check_type('energy bounds', value, Iterable, Real)
                self._constraints['energy_bounds'] = tuple(value)
            elif key == 'fissionable':
                cv.check_type('fissionable', value, bool)
                self._constraints['fissionable'] = value
            elif key == 'rejection_strategy':
                cv.check_value('rejection strategy', value, ('resample', 'kill'))
                self._constraints['rejection_strategy'] = value
            else:
                raise ValueError(f'Unknown key in constraints dictionary: {key}')

    @abstractmethod
    def populate_xml_element(self, element):
        """Add necessary source information to an XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """

    def to_xml_element(self) -> ET.Element:
        """Return XML representation of the source

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing source data

        """
        element = ET.Element("source")
        element.set("type", self.type)
        if self.strength is not None:
            element.set("strength", str(self.strength))
        self.populate_xml_element(element)
        constraints = self.constraints
        if constraints:
            constraints_elem = ET.SubElement(element, "constraints")
            if "domain_ids" in constraints:
                dt_elem = ET.SubElement(constraints_elem, "domain_type")
                dt_elem.text = constraints["domain_type"]
                id_elem = ET.SubElement(constraints_elem, "domain_ids")
                id_elem.text = ' '.join(str(uid) for uid in constraints["domain_ids"])
            if "time_bounds" in constraints:
                dt_elem = ET.SubElement(constraints_elem, "time_bounds")
                dt_elem.text = ' '.join(str(t) for t in constraints["time_bounds"])
            if "energy_bounds" in constraints:
                dt_elem = ET.SubElement(constraints_elem, "energy_bounds")
                dt_elem.text = ' '.join(str(E) for E in constraints["energy_bounds"])
            if "fissionable" in constraints:
                dt_elem = ET.SubElement(constraints_elem, "fissionable")
                dt_elem.text = str(constraints["fissionable"]).lower()
            if "rejection_strategy" in constraints:
                dt_elem = ET.SubElement(constraints_elem, "rejection_strategy")
                dt_elem.text = constraints["rejection_strategy"]

        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element, meshes=None) -> SourceBase:
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
            elif source_type == 'mesh':
                return MeshSource.from_xml_element(elem, meshes)
            else:
                raise ValueError(f'Source type {source_type} is not recognized')

    @staticmethod
    def _get_constraints(elem: ET.Element) -> dict[str, Any]:
        # Find element containing constraints
        constraints_elem = elem.find("constraints")
        elem = constraints_elem if constraints_elem is not None else elem

        constraints = {}
        domain_type = get_text(elem, "domain_type")
        if domain_type is not None:
            domain_ids = get_elem_list(elem, "domain_ids", int)

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
            constraints['domains'] = domains

        time_bounds = get_elem_list(elem, "time_bounds", float)
        if time_bounds is not None:
            constraints['time_bounds'] = time_bounds

        energy_bounds = get_elem_list(elem, "energy_bounds", float)
        if energy_bounds is not None:
            constraints['energy_bounds'] = energy_bounds

        fissionable = get_text(elem, "fissionable")
        if fissionable is not None:
            constraints['fissionable'] = fissionable in ('true', '1')

        rejection_strategy = get_text(elem, "rejection_strategy")
        if rejection_strategy is not None:
            constraints['rejection_strategy'] = rejection_strategy

        return constraints


class IndependentSource(SourceBase):
    """Distribution of phase space coordinates for source sites.

    .. versionadded:: 0.14.0

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

        .. deprecated:: 0.15.0
            Use the `constraints` argument instead.
    constraints : dict
        Constraints on sampled source particles. Valid keys include 'domains',
        'time_bounds', 'energy_bounds', 'fissionable', and 'rejection_strategy'.
        For 'domains', the corresponding value is an iterable of
        :class:`openmc.Cell`, :class:`openmc.Material`, or
        :class:`openmc.Universe` for which sampled sites must be within. For
        'time_bounds' and 'energy_bounds', the corresponding value is a sequence
        of floats giving the lower and upper bounds on time in [s] or energy in
        [eV] that the sampled particle must be within. For 'fissionable', the
        value is a bool indicating that only sites in fissionable material
        should be accepted. The 'rejection_strategy' indicates what should
        happen when a source particle is rejected: either 'resample' (pick a new
        particle) or 'kill' (accept and terminate).

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

    .. versionadded:: 0.14.0

    particle : {'neutron', 'photon'}
        Source particle type
    constraints : dict
        Constraints on sampled source particles. Valid keys include
        'domain_type', 'domain_ids', 'time_bounds', 'energy_bounds',
        'fissionable', and 'rejection_strategy'.

    """

    def __init__(
        self,
        space: openmc.stats.Spatial | None = None,
        angle: openmc.stats.UnitSphere | None = None,
        energy: openmc.stats.Univariate | None = None,
        time: openmc.stats.Univariate | None = None,
        strength: float = 1.0,
        particle: str = 'neutron',
        domains: Sequence[openmc.Cell | openmc.Material | openmc.Universe] | None = None,
        constraints: dict[str, Any] | None = None
    ):
        if domains is not None:
            warnings.warn("The 'domains' arguments has been replaced by the "
                          "'constraints' argument.", FutureWarning)
            constraints = {'domains': domains}

        super().__init__(strength=strength, constraints=constraints)

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
        self.particle = particle

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

    def populate_xml_element(self, element):
        """Add necessary source information to an XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """
        element.set("particle", self.particle)
        if self.space is not None:
            element.append(self.space.to_xml_element())
        if self.angle is not None:
            element.append(self.angle.to_xml_element())
        if self.energy is not None:
            element.append(self.energy.to_xml_element('energy'))
        if self.time is not None:
            element.append(self.time.to_xml_element('time'))

    @classmethod
    def from_xml_element(cls, elem: ET.Element, meshes=None) -> SourceBase:
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
        constraints = cls._get_constraints(elem)
        source = cls(constraints=constraints)

        strength = get_text(elem, 'strength')
        if strength is not None:
            source.strength = float(strength)

        particle = get_text(elem, 'particle')
        if particle is not None:
            source.particle = particle

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


class MeshSource(SourceBase):
    """A source with a spatial distribution over mesh elements

    This class represents a mesh-based source in which random positions are
    uniformly sampled within mesh elements and each element can have independent
    angle, energy, and time distributions. The element sampled is chosen based
    on the relative strengths of the sources applied to the elements. The
    strength of the mesh source as a whole is the sum of all source strengths
    applied to the elements.

    .. versionadded:: 0.15.0

    Parameters
    ----------
    mesh : openmc.MeshBase
        The mesh over which source sites will be generated.
    sources : sequence of openmc.SourceBase
        Sources for each element in the mesh. Sources must be specified as
        either a 1-D array in the order of the mesh indices or a
        multidimensional array whose shape matches the mesh shape. If spatial
        distributions are set on any of the source objects, they will be ignored
        during source site sampling.
    constraints : dict
        Constraints on sampled source particles. Valid keys include 'domains',
        'time_bounds', 'energy_bounds', 'fissionable', and 'rejection_strategy'.
        For 'domains', the corresponding value is an iterable of
        :class:`openmc.Cell`, :class:`openmc.Material`, or
        :class:`openmc.Universe` for which sampled sites must be within. For
        'time_bounds' and 'energy_bounds', the corresponding value is a sequence
        of floats giving the lower and upper bounds on time in [s] or energy in
        [eV] that the sampled particle must be within. For 'fissionable', the
        value is a bool indicating that only sites in fissionable material
        should be accepted. The 'rejection_strategy' indicates what should
        happen when a source particle is rejected: either 'resample' (pick a new
        particle) or 'kill' (accept and terminate).

    Attributes
    ----------
    mesh : openmc.MeshBase
        The mesh over which source sites will be generated.
    sources : numpy.ndarray of openmc.SourceBase
        Sources to apply to each element
    strength : float
        Strength of the source
    type : str
        Indicator of source type: 'mesh'
    constraints : dict
        Constraints on sampled source particles. Valid keys include
        'domain_type', 'domain_ids', 'time_bounds', 'energy_bounds',
        'fissionable', and 'rejection_strategy'.

    """
    def __init__(
            self,
            mesh: MeshBase,
            sources: Sequence[SourceBase],
            constraints: dict[str, Any] | None = None,
    ):
        super().__init__(strength=None, constraints=constraints)
        self.mesh = mesh
        self.sources = sources

    @property
    def type(self) -> str:
        return "mesh"

    @property
    def mesh(self) -> MeshBase:
        return self._mesh

    @property
    def strength(self) -> float:
        return sum(s.strength for s in self.sources)

    @property
    def sources(self) -> np.ndarray:
        return self._sources

    @mesh.setter
    def mesh(self, m):
        cv.check_type('source mesh', m, MeshBase)
        self._mesh = m

    @sources.setter
    def sources(self, s):
        cv.check_iterable_type('mesh sources', s, SourceBase, max_depth=3)

        s = np.asarray(s)

        if isinstance(self.mesh, StructuredMesh):
            if s.size != self.mesh.num_mesh_cells:
                raise ValueError(
                    f'The length of the source array ({s.size}) does not match '
                    f'the number of mesh elements ({self.mesh.num_mesh_cells}).')

            # If user gave a multidimensional array, flatten in the order
            # of the mesh indices
            if s.ndim > 1:
                s = s.ravel(order='F')

        elif isinstance(self.mesh, UnstructuredMesh):
            if s.ndim > 1:
                raise ValueError('Sources must be a 1-D array for unstructured mesh')

        self._sources = s
        for src in self._sources:
            if isinstance(src, IndependentSource) and src.space is not None:
                warnings.warn('Some sources on the mesh have spatial '
                              'distributions that will be ignored at runtime.')
                break

    @strength.setter
    def strength(self, val):
        if val is not None:
            cv.check_type('mesh source strength', val, Real)
            self.set_total_strength(val)

    def set_total_strength(self, strength: float):
        """Scales the element source strengths based on a desired total strength.

        Parameters
        ----------
        strength : float
            Total source strength

        """
        current_strength = self.strength if self.strength != 0.0 else 1.0

        for s in self.sources:
            s.strength *= strength / current_strength

    def normalize_source_strengths(self):
        """Update all element source strengths such that they sum to 1.0."""
        self.set_total_strength(1.0)

    def populate_xml_element(self, elem: ET.Element):
        """Add necessary source information to an XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """
        elem.set("mesh", str(self.mesh.id))

        # write in the order of mesh indices
        for s in self.sources:
            elem.append(s.to_xml_element())

    @classmethod
    def from_xml_element(cls, elem: ET.Element, meshes) -> openmc.MeshSource:
        """
        Generate MeshSource from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element
        meshes : dict
            A dictionary with mesh IDs as keys and openmc.MeshBase instances as
            values

        Returns
        -------
        openmc.MeshSource
            MeshSource generated from the XML element
        """
        mesh_id = int(get_text(elem, 'mesh'))
        mesh = meshes[mesh_id]

        sources = [SourceBase.from_xml_element(e) for e in elem.iterchildren('source')]
        constraints = cls._get_constraints(elem)
        return cls(mesh, sources, constraints=constraints)


def Source(*args, **kwargs):
    """
    A function for backward compatibility of sources. Will be removed in the
    future. Please update to IndependentSource.
    """
    warnings.warn("This class is deprecated in favor of 'IndependentSource'", FutureWarning)
    return openmc.IndependentSource(*args, **kwargs)


class CompiledSource(SourceBase):
    """A source based on a compiled shared library

    .. versionadded:: 0.14.0

    Parameters
    ----------
    library : path-like
        Path to a compiled shared library
    parameters : str
        Parameters to be provided to the compiled shared library function
    strength : float
        Strength of the source
    constraints : dict
        Constraints on sampled source particles. Valid keys include 'domains',
        'time_bounds', 'energy_bounds', 'fissionable', and 'rejection_strategy'.
        For 'domains', the corresponding value is an iterable of
        :class:`openmc.Cell`, :class:`openmc.Material`, or
        :class:`openmc.Universe` for which sampled sites must be within. For
        'time_bounds' and 'energy_bounds', the corresponding value is a sequence
        of floats giving the lower and upper bounds on time in [s] or energy in
        [eV] that the sampled particle must be within. For 'fissionable', the
        value is a bool indicating that only sites in fissionable material
        should be accepted. The 'rejection_strategy' indicates what should
        happen when a source particle is rejected: either 'resample' (pick a new
        particle) or 'kill' (accept and terminate).

    Attributes
    ----------
    library : pathlib.Path
        Path to a compiled shared library
    parameters : str
        Parameters to be provided to the compiled shared library function
    strength : float
        Strength of the source
    type : str
        Indicator of source type: 'compiled'
    constraints : dict
        Constraints on sampled source particles. Valid keys include
        'domain_type', 'domain_ids', 'time_bounds', 'energy_bounds',
        'fissionable', and 'rejection_strategy'.

    """
    def __init__(
        self,
        library: PathLike,
        parameters: str | None = None,
        strength: float = 1.0,
        constraints: dict[str, Any] | None = None
    ) -> None:
        super().__init__(strength=strength, constraints=constraints)
        self.library = library
        self._parameters = None
        if parameters is not None:
            self.parameters = parameters

    @property
    def type(self) -> str:
        return "compiled"

    @property
    def library(self) -> Path:
        return self._library

    @library.setter
    def library(self, library_name: PathLike):
        cv.check_type('library', library_name, PathLike)
        self._library = input_path(library_name)

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
        element.set("library", str(self.library))

        if self.parameters is not None:
            element.set("parameters", self.parameters)

    @classmethod
    def from_xml_element(cls, elem: ET.Element) -> openmc.CompiledSource:
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
        kwargs = {'constraints': cls._get_constraints(elem)}
        kwargs['library'] = get_text(elem, 'library')

        source = cls(**kwargs)

        strength = get_text(elem, 'strength')
        if strength is not None:
            source.strength = float(strength)

        parameters = get_text(elem, 'parameters')
        if parameters is not None:
            source.parameters = parameters

        return source


class FileSource(SourceBase):
    """A source based on particles stored in a file

    .. versionadded:: 0.14.0

    Parameters
    ----------
    path : path-like
        Path to the source file from which sites should be sampled
    strength : float
        Strength of the source (default is 1.0)
    constraints : dict
        Constraints on sampled source particles. Valid keys include 'domains',
        'time_bounds', 'energy_bounds', 'fissionable', and 'rejection_strategy'.
        For 'domains', the corresponding value is an iterable of
        :class:`openmc.Cell`, :class:`openmc.Material`, or
        :class:`openmc.Universe` for which sampled sites must be within. For
        'time_bounds' and 'energy_bounds', the corresponding value is a sequence
        of floats giving the lower and upper bounds on time in [s] or energy in
        [eV] that the sampled particle must be within. For 'fissionable', the
        value is a bool indicating that only sites in fissionable material
        should be accepted. The 'rejection_strategy' indicates what should
        happen when a source particle is rejected: either 'resample' (pick a new
        particle) or 'kill' (accept and terminate).

    Attributes
    ----------
    path : Pathlike
        Source file from which sites should be sampled
    strength : float
        Strength of the source
    type : str
        Indicator of source type: 'file'
    constraints : dict
        Constraints on sampled source particles. Valid keys include
        'domain_type', 'domain_ids', 'time_bounds', 'energy_bounds',
        'fissionable', and 'rejection_strategy'.

    """

    def __init__(
        self,
        path: PathLike,
        strength: float = 1.0,
        constraints: dict[str, Any] | None = None
    ):
        super().__init__(strength=strength, constraints=constraints)
        self.path = path

    @property
    def type(self) -> str:
        return "file"

    @property
    def path(self) -> PathLike:
        return self._path

    @path.setter
    def path(self, p: PathLike):
        cv.check_type('source file', p, PathLike)
        self._path = input_path(p)

    def populate_xml_element(self, element):
        """Add necessary file source information to an XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """
        if self.path is not None:
            element.set("file", str(self.path))

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
        kwargs = {'constraints': cls._get_constraints(elem)}
        kwargs['path'] = get_text(elem, 'file')
        strength = get_text(elem, 'strength')
        if strength is not None:
            kwargs['strength'] = float(strength)

        return cls(**kwargs)


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

    @classmethod
    def from_pdg_number(cls, pdg_number: int) -> ParticleType:
        """Constructs a ParticleType instance from a PDG number.

        The Particle Data Group at LBNL publishes a Monte Carlo particle
        numbering scheme as part of the `Review of Particle Physics
        <10.1103/PhysRevD.110.030001>`_. This method maps PDG numbers to the
        corresponding :class:`ParticleType`.

        Parameters
        ----------
        pdg_number : int
            The PDG number of the particle type.

        Returns
        -------
        The corresponding ParticleType instance.
        """
        try:
            return {
                2112: ParticleType.NEUTRON,
                22: ParticleType.PHOTON,
                11: ParticleType.ELECTRON,
                -11: ParticleType.POSITRON,
            }[pdg_number]
        except KeyError:
            raise ValueError(f"Unrecognized PDG number: {pdg_number}")

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
        r: Iterable[float] = (0., 0., 0.),
        u: Iterable[float] = (0., 0., 1.),
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
    source_particles: Iterable[SourceParticle],
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
    cv.check_iterable_type("source particles", source_particles, SourceParticle)
    pl = ParticleList(source_particles)
    pl.export_to_hdf5(filename, **kwargs)


class ParticleList(list):
    """A collection of SourceParticle objects.

    Parameters
    ----------
    particles : list of SourceParticle
        Particles to collect into the list

    """
    @classmethod
    def from_hdf5(cls, filename: PathLike) -> ParticleList:
        """Create particle list from an HDF5 file.

        Parameters
        ----------
        filename : path-like
            Path to source file to read.

        Returns
        -------
        ParticleList instance

        """
        with h5py.File(filename, 'r') as fh:
            filetype = fh.attrs['filetype']
            arr = fh['source_bank'][...]

        if filetype != b'source':
            raise ValueError(f'File {filename} is not a source file')

        source_particles = [
            SourceParticle(*params, ParticleType(particle))
            for *params, particle in arr
        ]
        return cls(source_particles)

    @classmethod
    def from_mcpl(cls, filename: PathLike) -> ParticleList:
        """Create particle list from an MCPL file.

        Parameters
        ----------
        filename : path-like
            Path to MCPL file to read.

        Returns
        -------
        ParticleList instance

        """
        import mcpl
        # Process .mcpl file
        particles = []
        with mcpl.MCPLFile(filename) as f:
            for particle in f.particles:
                # Determine particle type based on the PDG number
                try:
                    particle_type = ParticleType.from_pdg_number(particle.pdgcode)
                except ValueError:
                    particle_type = "UNKNOWN"

                # Create a source particle instance. Note that MCPL stores
                # energy in MeV and time in ms.
                source_particle = SourceParticle(
                    r=tuple(particle.position),
                    u=tuple(particle.direction),
                    E=1.0e6*particle.ekin,
                    time=1.0e-3*particle.time,
                    wgt=particle.weight,
                    particle=particle_type
                )
                particles.append(source_particle)

        return cls(particles)

    def __getitem__(self, index):
        """
        Return a new ParticleList object containing the particle(s)
        at the specified index or slice.

        Parameters
        ----------
        index : int, slice or list
            The index, slice or list to select from the list of particles

        Returns
        -------
        openmc.ParticleList or openmc.SourceParticle
            A new object with the selected particle(s)
        """
        if isinstance(index, int):
            # If it's a single integer, return the corresponding particle
            return super().__getitem__(index)
        elif isinstance(index, slice):
            # If it's a slice, return a new ParticleList object with the
            # sliced particles
            return ParticleList(super().__getitem__(index))
        elif isinstance(index, list):
            # If it's a list of integers, return a new ParticleList object with
            # the selected particles. Note that Python 3.10 gets confused if you
            # use super() here, so we call list.__getitem__ directly.
            return ParticleList([list.__getitem__(self, i) for i in index])
        else:
            raise TypeError(f"Invalid index type: {type(index)}. Must be int, "
                            "slice, or list of int.")

    def to_dataframe(self) -> pd.DataFrame:
        """A dataframe representing the source particles

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the source particles attributes.
        """
        # Extract the attributes of the source particles into a list of tuples
        data = [(sp.r[0], sp.r[1], sp.r[2], sp.u[0], sp.u[1], sp.u[2],
                 sp.E, sp.time, sp.wgt, sp.delayed_group, sp.surf_id,
                 sp.particle.name.lower()) for sp in self]

        # Define the column names for the DataFrame
        columns = ['x', 'y', 'z', 'u_x', 'u_y', 'u_z', 'E', 'time', 'wgt',
                   'delayed_group', 'surf_id', 'particle']

        # Create the pandas DataFrame from the data
        return pd.DataFrame(data, columns=columns)

    def export_to_hdf5(self, filename: PathLike, **kwargs):
        """Export particle list to an HDF5 file.

        This method write out an .h5 file that can be used as a source file in
        conjunction with the :class:`openmc.FileSource` class.

        Parameters
        ----------
        filename : path-like
            Path to source file to write
        **kwargs
            Keyword arguments to pass to :class:`h5py.File`

        See Also
        --------
        openmc.FileSource

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
        arr = np.array([s.to_tuple() for s in self], dtype=source_dtype)

        # Write array to file
        kwargs.setdefault('mode', 'w')
        with h5py.File(filename, **kwargs) as fh:
            fh.attrs['filetype'] = np.bytes_("source")
            fh.create_dataset('source_bank', data=arr, dtype=source_dtype)


def read_source_file(filename: PathLike) -> ParticleList:
    """Read a source file and return a list of source particles.

    .. versionadded:: 0.15.0

    Parameters
    ----------
    filename : str or path-like
        Path to source file to read

    Returns
    -------
    openmc.ParticleList

    See Also
    --------
    openmc.SourceParticle

    """
    filename = Path(filename)
    if filename.suffix not in ('.h5', '.mcpl'):
        raise ValueError('Source file must have a .h5 or .mcpl extension.')

    if filename.suffix == '.h5':
        return ParticleList.from_hdf5(filename)
    else:
        return ParticleList.from_mcpl(filename)


def read_collision_track_hdf5(filename):

    with h5py.File(filename, 'r') as file:
        data = file['collision_track_bank'][:]

    return data


def read_collision_track_mcpl(file_path):
    import mcpl
    myfile = mcpl.MCPLFile(file_path)
    data = {
        'r': [],  # for position (x, y, z)
        'u': [],  # for direction (ux, uy, uz)
        'E': [], 'dE': [], 'time': [],
        'wgt': [], 'event_mt': [], 'delayed_group': [],
        'cell_id': [], 'nuclide_id': [], 'material_id': [],
        'universe_id': [], 'n_collision': [], 'particle': [],
        'parent_id': [], 'progeny_id': []
    }

    # Read and collect data from the MCPL file
    for i, p in enumerate(myfile.particles):
        if 'blob_' + str(i) in myfile.blobs:
            blob_data = myfile.blobs['blob_' + str(i)]
            decoded_str = blob_data.decode('utf-8')
            pairs = decoded_str.split(';')
            values_dict = {k.strip(): v.strip()
                           for k, v in (pair.split(':') for pair in pairs if pair.strip())}

            data['r'].append((p.x, p.y, p.z))  # Append as tuple
            data['u'].append((p.ux, p.uy, p.uz))  # Append as tuple
            data['E'].append(p.ekin * 1e6)
            data['dE'].append(float(values_dict.get('dE', 0)))
            data['time'].append(p.time * 1e-3)
            data['wgt'].append(p.weight)
            data['event_mt'].append(int(values_dict.get('event_mt', 0)))
            data['delayed_group'].append(
                int(values_dict.get('delayed_group', 0)))
            data['cell_id'].append(int(values_dict.get('cell_id', 0)))
            data['nuclide_id'].append(int(values_dict.get('nuclide_id', 0)))
            data['material_id'].append(int(values_dict.get('material_id', 0)))
            data['universe_id'].append(int(values_dict.get('universe_id', 0)))
            data['n_collision'].append(int(values_dict.get('n_collision', 0)))
            data['particle'].append(p.pdgcode)
            data['parent_id'].append(int(values_dict.get('parent_id', 0)))
            data['progeny_id'].append(int(values_dict.get('progeny_id', 0)))

    dtypes = [
        ('r', [('x', 'f8'), ('y', 'f8'), ('z', 'f8')]),
        ('u', [('x', 'f8'), ('y', 'f8'), ('z', 'f8')]),
        ('E', 'f8'), ('dE', 'f8'), ('time', 'f8'), ('wgt', 'f8'),
        ('event_mt', 'f8'), ('delayed_group', 'i4'), ('cell_id', 'i4'),
        ('nuclide_id', 'i4'), ('material_id', 'i4'), ('universe_id', 'i4'),
        ('n_collision', 'i4'), ('particle', 'i4'),
        ('parent_id', 'i8'), ('progeny_id', 'i8')
    ]

    structured_array = np.zeros(len(data['r']), dtype=dtypes)
    for key in data:
        structured_array[key] = data[key]  # Assign data

    return structured_array


def read_collision_track_file(filename):

    filename = Path(filename)
    if filename.suffix not in ('.h5', '.mcpl'):
        raise ValueError('Retina file must have a .h5 or .mcpl extension.')

    if filename.suffix == '.h5':
        return read_collision_track_hdf5(filename)
    else:
        return read_collision_track_mcpl(filename)
