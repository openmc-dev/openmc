from __future__ import annotations
from abc import ABC, abstractmethod
from collections.abc import Iterable, Sequence
from numbers import Integral, Real
import os
from pathlib import Path
import warnings
from typing import Any, Self

import lxml.etree as ET
import numpy as np
import h5py
import pandas as pd
from scipy.io import netcdf_file

import openmc
import openmc.checkvalue as cv
from openmc.checkvalue import PathLike
from openmc.stats.multivariate import UnitSphere, Spatial
from openmc.stats.univariate import Univariate
from ._xml import get_elem_list, get_text
from .mesh import MeshBase, StructuredMesh, UnstructuredMesh
from .particle_type import ParticleType
from .statepoint import _VERSION_STATEPOINT
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
    type : {'independent', 'file', 'compiled', 'mesh', 'tokamak', 'stellarator'}
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
                cv.check_value('rejection strategy',
                               value, ('resample', 'kill'))
                self._constraints['rejection_strategy'] = value
            else:
                raise ValueError(
                    f'Unknown key in constraints dictionary: {key}')

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
                id_elem.text = ' '.join(str(uid)
                                        for uid in constraints["domain_ids"])
            if "time_bounds" in constraints:
                dt_elem = ET.SubElement(constraints_elem, "time_bounds")
                dt_elem.text = ' '.join(str(t)
                                        for t in constraints["time_bounds"])
            if "energy_bounds" in constraints:
                dt_elem = ET.SubElement(constraints_elem, "energy_bounds")
                dt_elem.text = ' '.join(str(E)
                                        for E in constraints["energy_bounds"])
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
            elif source_type == 'tokamak':
                return TokamakSource.from_xml_element(elem)
            elif source_type == 'stellarator':
                return StellaratorSource.from_xml_element(elem)
            else:
                raise ValueError(
                    f'Source type {source_type} is not recognized')

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
    particle : str or int or openmc.ParticleType
        Source particle type (name, PDG number, or type)
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
    particle : str or int or openmc.ParticleType
        Source particle type (alias, PDG number, or GNDS nuclide name)
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
        particle: str | int | ParticleType = 'neutron',
        domains: Sequence[openmc.Cell | openmc.Material |
                          openmc.Universe] | None = None,
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
    def particle(self) -> ParticleType:
        return self._particle

    @particle.setter
    def particle(self, particle):
        self._particle = ParticleType(particle)

    def populate_xml_element(self, element):
        """Add necessary source information to an XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """
        element.set("particle", str(self.particle))
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
            if s.size != self.mesh.n_elements:
                raise ValueError(
                    f'The length of the source array ({s.size}) does not match '
                    f'the number of mesh elements ({self.mesh.n_elements}).')

            # If user gave a multidimensional array, flatten in the order
            # of the mesh indices
            if s.ndim > 1:
                s = s.ravel(order='F')

        elif isinstance(self.mesh, UnstructuredMesh):
            if s.ndim > 1:
                raise ValueError(
                    'Sources must be a 1-D array for unstructured mesh')

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

        sources = [SourceBase.from_xml_element(
            e) for e in elem.iterchildren('source')]
        constraints = cls._get_constraints(elem)
        return cls(mesh, sources, constraints=constraints)


def Source(*args, **kwargs):
    """
    A function for backward compatibility of sources. Will be removed in the
    future. Please update to IndependentSource.
    """
    warnings.warn(
        "This class is deprecated in favor of 'IndependentSource'", FutureWarning)
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


class TokamakSource(SourceBase):
    r"""A source representing neutron emission from a tokamak plasma.

    This source samples neutron positions from a tokamak plasma geometry using
    Miller-style flux surface parameterization. The user provides an emission
    profile S(r/a) as a function of normalized minor radius, along with one or
    more energy distributions.

    The flux surface parameterization is

    .. math::

        \begin{aligned}
        R &= R_0 + r \cos\left(\alpha + \delta \sin\alpha\right)
             + \Delta \left[1 - \left(\frac{r}{a}\right)^2\right] \\
        Z &= Z_\mathrm{shift} + \kappa r \sin\alpha
        \end{aligned}

    where :math:`R_0` is major radius, :math:`a` is minor radius,
    :math:`\kappa` is elongation, :math:`\delta` is triangularity,
    :math:`\Delta` is the Shafranov shift, and :math:`Z_\mathrm{shift}` is
    the vertical shift.

    .. versionadded:: 0.16.0

    Parameters
    ----------
    major_radius : float
        Major radius R0 in [cm]
    minor_radius : float
        Minor radius a in [cm]
    elongation : float
        Plasma elongation κ (must be > 0)
    triangularity : float
        Plasma triangularity δ (must be in [-1, 1])
    shafranov_shift : float
        Shafranov shift Δ in [cm] (must be >= 0 and < a/2)
    r_over_a : numpy.ndarray
        Normalized minor radius grid points, must start at 0 and end at 1
    emission_density : numpy.ndarray
        Emission density S(r) at each r/a point (arbitrary units, must be >= 0).
        Values are linearly interpolated between grid points and refined on an
        internal grid for radial sampling. Must have the same length as
        ``r_over_a`` and contain at least one positive value.
    energy : openmc.stats.Univariate or Sequence[openmc.stats.Univariate]
        Energy distribution(s). Either a single distribution used at all radii,
        or one distribution per ``r_over_a`` grid point. When one distribution
        per grid point is given, the energy of a sampled particle is drawn from
        one of the two distributions bracketing its sampled radius, selected
        stochastically with probability proportional to the proximity of the
        radius to each grid point (stochastic interpolation).
    time : openmc.stats.Univariate, optional
        Time distribution of the source. If None, particles are born at
        :math:`t=0`, matching the default behavior of
        :class:`openmc.IndependentSource`.
    phi_start : float
        Starting toroidal angle in [rad] (default: 0)
    phi_extent : float
        Toroidal angle extent in [rad] (default: 2π)
    n_alpha : int
        Number of poloidal angle grid points for CDF sampling (default: 101)
    vertical_shift : float
        Vertical shift of the plasma center in [cm] (default: 0)
    strength : float
        Strength of the source (default: 1.0)
    constraints : dict
        Constraints on sampled source particles. See :class:`SourceBase` for
        valid keys and values.

    Attributes
    ----------
    major_radius : float
        Major radius R0 in [cm]
    minor_radius : float
        Minor radius a in [cm]
    elongation : float
        Plasma elongation κ
    triangularity : float
        Plasma triangularity δ
    shafranov_shift : float
        Shafranov shift Δ in [cm]
    r_over_a : numpy.ndarray
        Normalized minor radius grid points
    emission_density : numpy.ndarray
        Emission density S(r) at each r/a point
    energy : list of openmc.stats.Univariate
        Energy distribution(s)
    time : openmc.stats.Univariate or None
        Time distribution of the source
    phi_start : float
        Starting toroidal angle in [rad]
    phi_extent : float
        Toroidal angle extent in [rad]
    n_alpha : int
        Number of poloidal angle grid points
    vertical_shift : float
        Vertical shift of the plasma center in [cm]
    strength : float
        Strength of the source
    type : str
        Indicator of source type: 'tokamak'
    constraints : dict
        Constraints on sampled source particles

    """

    def __init__(
        self,
        major_radius: float,
        minor_radius: float,
        elongation: float,
        triangularity: float,
        shafranov_shift: float,
        r_over_a: Sequence[float],
        emission_density: Sequence[float],
        energy: Univariate | Sequence[Univariate],
        time: Univariate | None = None,
        phi_start: float = 0.0,
        phi_extent: float = 2.0 * np.pi,
        n_alpha: int = 101,
        vertical_shift: float = 0.0,
        strength: float = 1.0,
        constraints: dict[str, Any] | None = None
    ):
        super().__init__(strength=strength, constraints=constraints)
        self.major_radius = major_radius
        self.minor_radius = minor_radius
        self.elongation = elongation
        self.triangularity = triangularity
        self.shafranov_shift = shafranov_shift
        self.r_over_a = r_over_a
        self.emission_density = emission_density
        self.phi_start = phi_start
        self.phi_extent = phi_extent
        self.n_alpha = n_alpha
        self.vertical_shift = vertical_shift
        self.energy = energy
        self.time = time

        self._validate()

    def _validate(self):
        """Validate relationships between tokamak source parameters."""
        if self.minor_radius >= self.major_radius:
            raise ValueError(
                f"minor_radius ({self.minor_radius}) must be smaller than "
                f"major_radius ({self.major_radius})")
        if self.shafranov_shift >= 0.5 * self.minor_radius:
            raise ValueError(
                f"shafranov_shift ({self.shafranov_shift}) must be smaller "
                f"than half the minor_radius ({0.5 * self.minor_radius})")
        if len(self.emission_density) != len(self.r_over_a):
            raise ValueError(
                f"emission_density (length {len(self.emission_density)}) must "
                f"have the same length as r_over_a (length {len(self.r_over_a)})")
        if not np.any(self.emission_density > 0.0):
            raise ValueError("emission_density must contain a positive value")
        if len(self.energy) not in (1, len(self.r_over_a)):
            raise ValueError(
                f"Number of energy distributions ({len(self.energy)}) must be "
                f"either 1 or equal to the number of r_over_a grid points "
                f"({len(self.r_over_a)})")

    @property
    def type(self) -> str:
        return "tokamak"

    @property
    def major_radius(self) -> float:
        return self._major_radius

    @major_radius.setter
    def major_radius(self, value: float):
        cv.check_type('major radius', value, Real)
        cv.check_greater_than('major radius', value, 0.0)
        self._major_radius = value

    @property
    def minor_radius(self) -> float:
        return self._minor_radius

    @minor_radius.setter
    def minor_radius(self, value: float):
        cv.check_type('minor radius', value, Real)
        cv.check_greater_than('minor radius', value, 0.0)
        self._minor_radius = value

    @property
    def elongation(self) -> float:
        return self._elongation

    @elongation.setter
    def elongation(self, value: float):
        cv.check_type('elongation', value, Real)
        cv.check_greater_than('elongation', value, 0.0)
        self._elongation = value

    @property
    def triangularity(self) -> float:
        return self._triangularity

    @triangularity.setter
    def triangularity(self, value: float):
        cv.check_type('triangularity', value, Real)
        cv.check_greater_than('triangularity', value, -1.0, equality=True)
        cv.check_less_than('triangularity', value, 1.0, equality=True)
        self._triangularity = value

    @property
    def shafranov_shift(self) -> float:
        return self._shafranov_shift

    @shafranov_shift.setter
    def shafranov_shift(self, value: float):
        cv.check_type('Shafranov shift', value, Real)
        cv.check_greater_than('Shafranov shift', value, 0.0, equality=True)
        self._shafranov_shift = value

    @property
    def r_over_a(self) -> np.ndarray:
        return self._r_over_a

    @r_over_a.setter
    def r_over_a(self, value: Sequence[float]):
        value = np.asarray(value, dtype=float)
        if value.ndim != 1 or len(value) < 2:
            raise ValueError("r_over_a must be a 1-D array with at least 2 points")
        if value[0] != 0.0:
            raise ValueError("r_over_a must start at 0")
        if value[-1] != 1.0:
            raise ValueError("r_over_a must end at 1")
        if not np.all(np.diff(value) > 0):
            raise ValueError("r_over_a must be strictly increasing")
        self._r_over_a = value

    @property
    def emission_density(self) -> np.ndarray:
        return self._emission_density

    @emission_density.setter
    def emission_density(self, value: Sequence[float]):
        value = np.asarray(value, dtype=float)
        if value.ndim != 1:
            raise ValueError("emission_density must be a 1-D array")
        if np.any(value < 0):
            raise ValueError("emission_density values cannot be negative")
        self._emission_density = value

    @property
    def energy(self) -> list[Univariate]:
        return self._energy

    @energy.setter
    def energy(self, value: Univariate | Sequence[Univariate]):
        if isinstance(value, Univariate):
            self._energy = [value]
        else:
            cv.check_iterable_type('energy distributions', value, Univariate)
            self._energy = list(value)

    @property
    def time(self) -> Univariate | None:
        return self._time

    @time.setter
    def time(self, value: Univariate | None):
        if value is not None:
            cv.check_type('time distribution', value, Univariate)
        self._time = value

    @property
    def phi_start(self) -> float:
        return self._phi_start

    @phi_start.setter
    def phi_start(self, value: float):
        cv.check_type('phi_start', value, Real)
        self._phi_start = value

    @property
    def phi_extent(self) -> float:
        return self._phi_extent

    @phi_extent.setter
    def phi_extent(self, value: float):
        cv.check_type('phi_extent', value, Real)
        cv.check_greater_than('phi_extent', value, 0.0)
        cv.check_less_than('phi_extent', value, 2.0 * np.pi, equality=True)
        self._phi_extent = value

    @property
    def n_alpha(self) -> int:
        return self._n_alpha

    @n_alpha.setter
    def n_alpha(self, value: int):
        cv.check_type('n_alpha', value, Integral)
        cv.check_greater_than('n_alpha', value, 2)
        if value < 51:
            warnings.warn(
                "n_alpha values below 51 may introduce noticeable "
                "discretization bias in tokamak source sampling", stacklevel=2)
        self._n_alpha = value

    @property
    def vertical_shift(self) -> float:
        return self._vertical_shift

    @vertical_shift.setter
    def vertical_shift(self, value: float):
        cv.check_type('vertical shift', value, Real)
        self._vertical_shift = value

    def populate_xml_element(self, element):
        """Add necessary tokamak source information to an XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """
        self._validate()

        # Geometry parameters
        ET.SubElement(element, "major_radius").text = str(self.major_radius)
        ET.SubElement(element, "minor_radius").text = str(self.minor_radius)
        ET.SubElement(element, "elongation").text = str(self.elongation)
        ET.SubElement(element, "triangularity").text = str(self.triangularity)
        ET.SubElement(element, "shafranov_shift").text = str(self.shafranov_shift)

        # Toroidal angle bounds
        ET.SubElement(element, "phi_start").text = str(self.phi_start)
        ET.SubElement(element, "phi_extent").text = str(self.phi_extent)

        # Poloidal sampling resolution
        ET.SubElement(element, "n_alpha").text = str(self.n_alpha)

        # Vertical shift
        if self.vertical_shift != 0.0:
            ET.SubElement(element, "vertical_shift").text = str(self.vertical_shift)

        # Emission profile
        ET.SubElement(element, "r_over_a").text = ' '.join(str(r) for r in self.r_over_a)
        ET.SubElement(element, "emission_density").text = ' '.join(str(s) for s in self.emission_density)

        # Energy distribution(s)
        for dist in self.energy:
            element.append(dist.to_xml_element('energy'))

        # Time distribution
        if self.time is not None:
            element.append(self.time.to_xml_element('time'))

    @classmethod
    def from_xml_element(cls, elem: ET.Element) -> TokamakSource:
        """Generate tokamak source from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.TokamakSource
            Source generated from XML element

        """
        # Read geometry parameters
        major_radius = float(get_text(elem, 'major_radius'))
        minor_radius = float(get_text(elem, 'minor_radius'))
        elongation = float(get_text(elem, 'elongation'))
        triangularity = float(get_text(elem, 'triangularity'))
        shafranov_shift = float(get_text(elem, 'shafranov_shift'))

        # Read optional parameters
        phi_start_text = get_text(elem, 'phi_start')
        phi_start = float(phi_start_text) if phi_start_text else 0.0

        phi_extent_text = get_text(elem, 'phi_extent')
        phi_extent = float(phi_extent_text) if phi_extent_text else 2.0 * np.pi

        n_alpha_text = get_text(elem, 'n_alpha')
        n_alpha = int(n_alpha_text) if n_alpha_text else 101

        vertical_shift_text = get_text(elem, 'vertical_shift')
        vertical_shift = float(vertical_shift_text) if vertical_shift_text else 0.0

        # Read emission profile
        r_over_a = np.array([float(x) for x in get_text(elem, 'r_over_a').split()])
        emission_density = np.array([float(x) for x in get_text(elem, 'emission_density').split()])

        # Read energy distributions
        energy = [Univariate.from_xml_element(e) for e in elem.findall('energy')]
        if len(energy) == 1:
            energy = energy[0]

        # Read time distribution
        time_elem = elem.find('time')
        time = Univariate.from_xml_element(time_elem) if time_elem is not None else None

        # Read constraints and strength
        constraints = cls._get_constraints(elem)
        strength_text = get_text(elem, 'strength')
        strength = float(strength_text) if strength_text else 1.0

        return cls(
            major_radius=major_radius,
            minor_radius=minor_radius,
            elongation=elongation,
            triangularity=triangularity,
            shafranov_shift=shafranov_shift,
            r_over_a=r_over_a,
            emission_density=emission_density,
            energy=energy,
            time=time,
            phi_start=phi_start,
            phi_extent=phi_extent,
            n_alpha=n_alpha,
            vertical_shift=vertical_shift,
            strength=strength,
            constraints=constraints
        )


class StellaratorSource(SourceBase):
    r"""A source representing neutron emission from a stellarator plasma.

    This source samples neutron positions from a 3-D stellarator plasma
    described by the flux-surface Fourier representation shared by the VMEC
    and DESC equilibrium codes. The flux coordinates are
    :math:`(\rho, \theta, \zeta)` where :math:`\rho = \sqrt{s}` is the square
    root of the normalized toroidal flux (proportional to the average minor
    radius), :math:`\theta` is the poloidal angle, and :math:`\zeta` is the
    toroidal angle, which coincides with the cylindrical azimuthal angle
    :math:`\phi` in both codes. Flux surfaces are given by

    .. math::

        \begin{aligned}
        R(\rho,\theta,\zeta) &= \sum_k \left[ R^c_k(\rho)
            \cos(m_k\theta - n_k N_{fp} \zeta)
            + R^s_k(\rho) \sin(m_k\theta - n_k N_{fp}\zeta) \right] \\
        Z(\rho,\theta,\zeta) &= \sum_k \left[ Z^s_k(\rho)
            \sin(m_k\theta - n_k N_{fp}\zeta)
            + Z^c_k(\rho) \cos(m_k\theta - n_k N_{fp}\zeta) \right]
        \end{aligned}

    where :math:`N_{fp}` is the number of field periods. For
    stellarator-symmetric equilibria only :math:`R^c` and :math:`Z^s` are
    non-zero.

    Because :math:`\phi = \zeta`, the volume element is
    :math:`dV = R\,|\tau|\, d\rho\, d\theta\, d\zeta` with the poloidal-plane
    Jacobian :math:`\tau = \partial_\rho R\, \partial_\theta Z -
    \partial_\theta R \,\partial_\rho Z`, so for an emission density
    :math:`S(\rho)` constant on flux surfaces the joint density is
    :math:`p(\rho,\theta,\zeta) \propto S(\rho) R |\tau|`. The radial
    coordinate is sampled from the marginal
    :math:`p(\rho) \propto S(\rho) V'(\rho)` (with
    :math:`V'(\rho) = \oint\oint R|\tau| \,d\theta\, d\zeta` the differential
    volume for the interpolated geometry) by integrating a cubic PDF in each
    radial interval and using rejection sampling within the selected interval.
    The two angles are then sampled from the conditional
    :math:`p(\theta,\zeta|\rho) \propto R|\tau|` by rejection against a
    precomputed per-radial-bin majorant.

    .. versionadded:: 0.16.1

    Parameters
    ----------
    rho : numpy.ndarray
        Radial grid points :math:`\rho = \sqrt{s}`, must start at 0 and end
        at 1 and be strictly increasing.
    emission_density : numpy.ndarray
        Emission density :math:`S(\rho)` at each grid point (arbitrary units,
        must be >= 0). Only the shape matters; it is normalized internally.
    mode_m : numpy.ndarray of int
        Poloidal mode numbers :math:`m_k \ge 0`.
    mode_n : numpy.ndarray of int
        Toroidal mode numbers :math:`n_k` in units of the number of field
        periods (VMEC's ``xn``/``nfp``).
    rmnc : numpy.ndarray
        Cosine Fourier coefficients of :math:`R` in [cm] with shape
        ``(len(rho), len(mode_m))``.
    zmns : numpy.ndarray
        Sine Fourier coefficients of :math:`Z` in [cm] with shape
        ``(len(rho), len(mode_m))``.
    energy : openmc.stats.Univariate or Sequence[openmc.stats.Univariate]
        Energy distribution(s). Either a single distribution used at all
        radii, or one distribution per ``rho`` grid point (selected by
        stochastic interpolation between the two bracketing grid points).
    num_field_periods : int
        Number of field periods :math:`N_{fp}` (default: 1).
    rmns : numpy.ndarray, optional
        Sine Fourier coefficients of :math:`R` in [cm] for
        non-stellarator-symmetric equilibria.
    zmnc : numpy.ndarray, optional
        Cosine Fourier coefficients of :math:`Z` in [cm] for
        non-stellarator-symmetric equilibria.
    time : openmc.stats.Univariate, optional
        Time distribution of the source. If None, particles are born at
        :math:`t = 0`.
    strength : float
        Strength of the source (default: 1.0)
    constraints : dict
        Constraints on sampled source particles. See :class:`SourceBase` for
        valid keys and values.

    Attributes
    ----------
    rho : numpy.ndarray
        Radial grid points
    emission_density : numpy.ndarray
        Emission density at each grid point
    mode_m : numpy.ndarray
        Poloidal mode numbers
    mode_n : numpy.ndarray
        Toroidal mode numbers (per field period)
    rmnc : numpy.ndarray
        Cosine coefficients of R in [cm]
    zmns : numpy.ndarray
        Sine coefficients of Z in [cm]
    rmns : numpy.ndarray or None
        Sine coefficients of R in [cm]
    zmnc : numpy.ndarray or None
        Cosine coefficients of Z in [cm]
    num_field_periods : int
        Number of field periods
    energy : list of openmc.stats.Univariate
        Energy distribution(s)
    time : openmc.stats.Univariate or None
        Time distribution of the source
    strength : float
        Strength of the source
    type : str
        Indicator of source type: 'stellarator'
    constraints : dict
        Constraints on sampled source particles

    """

    def __init__(
        self,
        rho: Sequence[float],
        emission_density: Sequence[float],
        mode_m: Sequence[int],
        mode_n: Sequence[int],
        rmnc: Sequence[Sequence[float]],
        zmns: Sequence[Sequence[float]],
        energy: Univariate | Sequence[Univariate],
        num_field_periods: int = 1,
        rmns: Sequence[Sequence[float]] | None = None,
        zmnc: Sequence[Sequence[float]] | None = None,
        time: Univariate | None = None,
        strength: float = 1.0,
        constraints: dict[str, Any] | None = None
    ):
        super().__init__(strength=strength, constraints=constraints)
        self.rho = rho
        self.emission_density = emission_density
        self.mode_m = mode_m
        self.mode_n = mode_n
        self.rmnc = rmnc
        self.zmns = zmns
        self.rmns = rmns
        self.zmnc = zmnc
        self.num_field_periods = num_field_periods
        self.energy = energy
        self.time = time

        self._validate()

    def _validate(self):
        """Validate relationships between stellarator source parameters."""
        n_rho = len(self.rho)
        n_modes = len(self.mode_m)
        if len(self.emission_density) != n_rho:
            raise ValueError(
                f"emission_density (length {len(self.emission_density)}) must "
                f"have the same length as rho (length {n_rho})")
        if not np.any(self.emission_density > 0.0):
            raise ValueError("emission_density must contain a positive value")
        if len(self.mode_n) != n_modes:
            raise ValueError(
                f"mode_n (length {len(self.mode_n)}) must have the same "
                f"length as mode_m (length {n_modes})")
        for name in ('rmnc', 'zmns', 'rmns', 'zmnc'):
            coeff = getattr(self, name)
            if coeff is not None and coeff.shape != (n_rho, n_modes):
                raise ValueError(
                    f"{name} must have shape (len(rho), len(mode_m)) = "
                    f"({n_rho}, {n_modes}), got {coeff.shape}")
        if (self.rmns is None) != (self.zmnc is None):
            raise ValueError(
                "rmns and zmnc must both be given for non-stellarator-"
                "symmetric equilibria")
        if len(self.energy) not in (1, n_rho):
            raise ValueError(
                f"Number of energy distributions ({len(self.energy)}) must be "
                f"either 1 or equal to the number of rho grid points "
                f"({n_rho})")

    @property
    def type(self) -> str:
        return "stellarator"

    @property
    def rho(self) -> np.ndarray:
        return self._rho

    @rho.setter
    def rho(self, value: Sequence[float]):
        value = np.asarray(value, dtype=float)
        if value.ndim != 1 or len(value) < 2:
            raise ValueError("rho must be a 1-D array with at least 2 points")
        if value[0] != 0.0:
            raise ValueError("rho must start at 0")
        if value[-1] != 1.0:
            raise ValueError("rho must end at 1")
        if not np.all(np.diff(value) > 0):
            raise ValueError("rho must be strictly increasing")
        self._rho = value

    @property
    def emission_density(self) -> np.ndarray:
        return self._emission_density

    @emission_density.setter
    def emission_density(self, value: Sequence[float]):
        value = np.asarray(value, dtype=float)
        if value.ndim != 1:
            raise ValueError("emission_density must be a 1-D array")
        if np.any(value < 0):
            raise ValueError("emission_density values cannot be negative")
        self._emission_density = value

    @property
    def mode_m(self) -> np.ndarray:
        return self._mode_m

    @mode_m.setter
    def mode_m(self, value: Sequence[int]):
        value = np.asarray(value, dtype=int)
        if value.ndim != 1 or len(value) < 1:
            raise ValueError("mode_m must be a 1-D array with at least 1 mode")
        if np.any(value < 0):
            raise ValueError("mode_m values must be >= 0")
        self._mode_m = value

    @property
    def mode_n(self) -> np.ndarray:
        return self._mode_n

    @mode_n.setter
    def mode_n(self, value: Sequence[int]):
        value = np.asarray(value, dtype=int)
        if value.ndim != 1:
            raise ValueError("mode_n must be a 1-D array")
        self._mode_n = value

    @staticmethod
    def _check_coeff(name, value, none_ok=False):
        if value is None:
            if none_ok:
                return None
            raise ValueError(f"{name} must be given")
        value = np.asarray(value, dtype=float)
        if value.ndim != 2:
            raise ValueError(f"{name} must be a 2-D array with shape "
                             "(len(rho), len(mode_m))")
        return value

    @property
    def rmnc(self) -> np.ndarray:
        return self._rmnc

    @rmnc.setter
    def rmnc(self, value):
        self._rmnc = self._check_coeff('rmnc', value)

    @property
    def zmns(self) -> np.ndarray:
        return self._zmns

    @zmns.setter
    def zmns(self, value):
        self._zmns = self._check_coeff('zmns', value)

    @property
    def rmns(self) -> np.ndarray | None:
        return self._rmns

    @rmns.setter
    def rmns(self, value):
        self._rmns = self._check_coeff('rmns', value, none_ok=True)

    @property
    def zmnc(self) -> np.ndarray | None:
        return self._zmnc

    @zmnc.setter
    def zmnc(self, value):
        self._zmnc = self._check_coeff('zmnc', value, none_ok=True)

    @property
    def num_field_periods(self) -> int:
        return self._num_field_periods

    @num_field_periods.setter
    def num_field_periods(self, value: int):
        cv.check_type('num_field_periods', value, Integral)
        cv.check_greater_than('num_field_periods', value, 0)
        self._num_field_periods = int(value)

    @property
    def energy(self) -> list[Univariate]:
        return self._energy

    @energy.setter
    def energy(self, value: Univariate | Sequence[Univariate]):
        if isinstance(value, Univariate):
            self._energy = [value]
        else:
            cv.check_iterable_type('energy distributions', value, Univariate)
            self._energy = list(value)

    @property
    def time(self) -> Univariate | None:
        return self._time

    @time.setter
    def time(self, value: Univariate | None):
        if value is not None:
            cv.check_type('time distribution', value, Univariate)
        self._time = value

    @staticmethod
    def _evaluate_emission_density(emission_density, rho):
        """Evaluate an emission density callable or validate an array."""
        if callable(emission_density):
            return np.asarray(emission_density(rho), dtype=float)
        emission_density = np.asarray(emission_density, dtype=float)
        if emission_density.shape != rho.shape:
            raise ValueError(
                f"emission_density array (length {len(emission_density)}) "
                f"must have one value per radial surface ({len(rho)}); "
                "alternatively provide a callable S(rho)")
        return emission_density

    @classmethod
    def from_vmec(
        cls,
        wout: PathLike,
        emission_density,
        energy: Univariate | Sequence[Univariate],
        **kwargs
    ) -> Self:
        """Generate a stellarator source from a VMEC wout file.

        The VMEC radial grid is uniform in the normalized toroidal flux
        :math:`s`; it is relabeled here as :math:`\\rho_j = \\sqrt{s_j}`, which
        leaves the per-surface Fourier coefficients unchanged. VMEC's ``xn``
        convention (kernel :math:`\\cos(m u - x_n v)` with ``xn`` a multiple of
        the number of field periods) is converted to per-field-period mode
        numbers. Coefficients are converted from [m] to [cm].

        Classic-format (NetCDF3) wout files are read with
        :func:`scipy.io.netcdf_file` and NetCDF4-format files (which are
        HDF5-based) with :mod:`h5py`, so no additional packages are required.

        Parameters
        ----------
        wout : path-like
            Path to a VMEC ``wout_*.nc`` NetCDF output file.
        emission_density : callable or numpy.ndarray
            Either a callable ``S(rho)`` evaluated at the radial grid points,
            or an array with one value per VMEC radial surface (arbitrary
            units, must be >= 0).
        energy : openmc.stats.Univariate or Sequence[openmc.stats.Univariate]
            Energy distribution(s); see the class docstring.
        **kwargs
            Additional keyword arguments passed to the constructor
            (e.g. ``time``, ``strength``, ``constraints``).

        Returns
        -------
        openmc.StellaratorSource

        """
        ds = cls._read_wout(wout)
        rmnc = ds['rmnc']
        zmns = ds['zmns']
        xm = ds['xm'].astype(int)
        xn = ds['xn'].astype(int)
        nfp = int(ds['nfp'])
        lasym = bool(ds.get('lasym__logical__', 0))
        rmns = ds['rmns'] if lasym else None
        zmnc = ds['zmnc'] if lasym else None

        ns = rmnc.shape[0]
        s = np.linspace(0.0, 1.0, ns)
        rho = np.sqrt(s)

        # m to cm
        rmnc = rmnc * 100.0
        zmns = zmns * 100.0
        if lasym:
            rmns = rmns * 100.0
            zmnc = zmnc * 100.0

        return cls(
            rho=rho,
            emission_density=cls._evaluate_emission_density(
                emission_density, rho),
            mode_m=xm,
            mode_n=xn // nfp,
            rmnc=rmnc,
            zmns=zmns,
            rmns=rmns,
            zmnc=zmnc,
            num_field_periods=nfp,
            energy=energy,
            **kwargs
        )

    @staticmethod
    def _read_wout(path):
        """Read the needed variables from a VMEC wout file.

        Classic (NetCDF3) files are read with scipy; NetCDF4-format files are
        HDF5-based and read with h5py.
        """
        names = ('rmnc', 'zmns', 'rmns', 'zmnc', 'xm', 'xn', 'nfp',
                 'lasym__logical__')
        try:
            with netcdf_file(str(path), mmap=False) as ds:
                return {k: np.asarray(v[()], dtype=float)
                        for k, v in ds.variables.items() if k in names}
        except (OSError, TypeError, ValueError):
            with h5py.File(path, 'r') as ds:
                return {k: np.asarray(ds[k][()], dtype=float)
                        for k in names if k in ds}

    @classmethod
    def from_desc(
        cls,
        eq,
        emission_density,
        energy: Univariate | Sequence[Univariate],
        n_rho: int = 33,
        **kwargs
    ) -> Self:
        """Generate a stellarator source from a DESC equilibrium.

        DESC represents :math:`R` and :math:`Z` in a Fourier-Zernike basis,
        :math:`X(\\rho,\\theta,\\zeta) = \\sum_{lmn} X_{lmn}\\,
        \\mathcal{R}_l^{|m|}(\\rho)\\, \\mathcal{F}_m(\\theta)\\,
        \\mathcal{F}_n(N_{fp}\\zeta)`, where :math:`\\mathcal{R}_l^{|m|}` is the
        (unnormalized) Zernike radial polynomial and :math:`\\mathcal{F}_m(x)`
        is :math:`\\cos(|m|x)` for :math:`m \\ge 0` and :math:`\\sin(|m|x)`
        otherwise. DESC evaluates the radial basis on a uniform grid in
        :math:`\\rho` and converts the product-form coefficients to the
        combined VMEC-style form. Coefficients are converted from [m] to [cm].
        This method requires the optional ``desc-opt`` package, which is
        imported only when the method is called.

        Parameters
        ----------
        eq : desc.equilibrium.Equilibrium or desc.equilibrium.EquilibriaFamily or path-like
            A DESC equilibrium object, or a path to a DESC HDF5 output file
            (the last equilibrium is used for an equilibrium family).
        emission_density : callable or numpy.ndarray
            Either a callable ``S(rho)`` evaluated at the radial grid points,
            or an array of length ``n_rho`` (arbitrary units, must be >= 0).
        energy : openmc.stats.Univariate or Sequence[openmc.stats.Univariate]
            Energy distribution(s); see the class docstring.
        n_rho : int
            Number of radial surfaces to extract (default: 33).
        **kwargs
            Additional keyword arguments passed to the constructor
            (e.g. ``time``, ``strength``, ``constraints``).

        Returns
        -------
        openmc.StellaratorSource

        """
        from desc.basis import zernike_radial
        from desc.equilibrium import Equilibrium, EquilibriaFamily
        from desc.io import load
        from desc.vmec_utils import ptolemy_identity_rev

        cv.check_type('n_rho', n_rho, Integral)
        cv.check_greater_than('n_rho', n_rho, 1)
        rho = np.linspace(0.0, 1.0, n_rho)
        if isinstance(eq, (str, os.PathLike)):
            eq = load(eq, file_format='hdf5')
        if isinstance(eq, EquilibriaFamily):
            if not len(eq):
                raise ValueError("The DESC equilibrium family is empty")
            eq = eq[-1]
        cv.check_type('eq', eq, Equilibrium)

        def convert(name):
            basis = getattr(eq, f'{name}_basis')
            l, m, _ = basis.modes.T  # noqa: E741
            modes, inverse = np.unique(basis.modes[:, 1:], axis=0, return_inverse=True)
            radial = np.asarray(zernike_radial(rho[:, None], l, m))
            values = np.zeros((n_rho, len(modes)))
            # Sum radial modes sharing the same angular mode on each surface.
            np.add.at(values.T, inverse,
                      (radial * np.asarray(getattr(eq, f'{name}_lmn'))).T)
            m, n, sine, cosine = ptolemy_identity_rev(modes[:, 0], modes[:, 1], values)
            return list(zip(m, n)), sine, cosine

        r_modes, r_sine, r_cosine = convert('R')
        z_modes, z_sine, z_cosine = convert('Z')
        modes = sorted(set(r_modes) | set(z_modes))
        mode_m, mode_n = np.asarray(modes, dtype=int).T
        indices = {mode: k for k, mode in enumerate(modes)}
        rmnc, rmns, zmnc, zmns = (
            np.zeros((n_rho, len(modes))) for _ in range(4))
        # Align the R/Z mode lists and convert meters to centimeters.
        r_indices = [indices[mode] for mode in r_modes]
        z_indices = [indices[mode] for mode in z_modes]
        rmnc[:, r_indices], rmns[:, r_indices] = 100*r_cosine, 100*r_sine
        zmnc[:, z_indices], zmns[:, z_indices] = 100*z_cosine, 100*z_sine

        # Drop modes whose R and Z coefficients are exactly zero on every surface.
        keep = np.any((rmnc != 0) | (rmns != 0) | (zmnc != 0) | (zmns != 0), axis=0)
        mode_m, mode_n = mode_m[keep], mode_n[keep]
        rmnc, rmns, zmnc, zmns = (a[:, keep] for a in (rmnc, rmns, zmnc, zmns))

        # Drop the asymmetric tables for stellarator-symmetric equilibria
        sym = not (np.any(rmns) or np.any(zmnc))

        return cls(
            rho=rho,
            emission_density=cls._evaluate_emission_density(
                emission_density, rho),
            mode_m=mode_m,
            mode_n=mode_n,
            rmnc=rmnc,
            zmns=zmns,
            rmns=None if sym else rmns,
            zmnc=None if sym else zmnc,
            num_field_periods=eq.NFP,
            energy=energy,
            **kwargs
        )

    def populate_xml_element(self, element):
        """Add necessary stellarator source information to an XML element

        Returns
        -------
        element : lxml.etree._Element
            XML element containing source data

        """
        self._validate()

        ET.SubElement(element, "num_field_periods").text = str(self.num_field_periods)
        ET.SubElement(element, "rho").text = ' '.join(str(r) for r in self.rho)
        ET.SubElement(element, "emission_density").text = ' '.join(
            str(s) for s in self.emission_density)
        ET.SubElement(element, "mode_m").text = ' '.join(str(m) for m in self.mode_m)
        ET.SubElement(element, "mode_n").text = ' '.join(str(n) for n in self.mode_n)

        # Coefficient tables flattened row-major (surface index varies slowest)
        for name in ('rmnc', 'zmns', 'rmns', 'zmnc'):
            coeff = getattr(self, name)
            if coeff is not None:
                ET.SubElement(element, name).text = \
                    ' '.join(str(c) for c in coeff.ravel())

        # Energy distribution(s)
        for dist in self.energy:
            element.append(dist.to_xml_element('energy'))

        # Time distribution
        if self.time is not None:
            element.append(self.time.to_xml_element('time'))

    @classmethod
    def from_xml_element(cls, elem: ET.Element) -> Self:
        """Generate stellarator source from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.StellaratorSource
            Source generated from XML element

        """
        num_field_periods = int(get_text(elem, 'num_field_periods', 1))

        rho = np.array(get_text(elem, 'rho').split(), dtype=float)
        emission_density = np.array(
            get_text(elem, 'emission_density').split(), dtype=float)
        mode_m = np.array(get_text(elem, 'mode_m').split(), dtype=int)
        mode_n = np.array(get_text(elem, 'mode_n').split(), dtype=int)

        shape = (len(rho), len(mode_m))
        coeffs = {}
        for name in ('rmnc', 'zmns', 'rmns', 'zmnc'):
            text = get_text(elem, name)
            coeffs[name] = (np.array(text.split(), dtype=float).reshape(shape)
                            if text else None)

        # Read energy distributions
        energy = [Univariate.from_xml_element(e) for e in elem.findall('energy')]

        # Read time distribution
        time_elem = elem.find('time')
        time = Univariate.from_xml_element(time_elem) if time_elem is not None else None

        # Read constraints and strength
        constraints = cls._get_constraints(elem)
        strength = float(get_text(elem, 'strength', 1.0))

        return cls(
            rho=rho,
            emission_density=emission_density,
            mode_m=mode_m,
            mode_n=mode_n,
            rmnc=coeffs['rmnc'],
            zmns=coeffs['zmns'],
            rmns=coeffs['rmns'],
            zmnc=coeffs['zmnc'],
            num_field_periods=num_field_periods,
            energy=energy,
            time=time,
            strength=strength,
            constraints=constraints
        )


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
    particle : ParticleType or str or int
        Type of the particle (type, name, or PDG number)

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
        particle: ParticleType | str | int = ParticleType.NEUTRON
    ):

        self.r = tuple(r)
        self.u = tuple(u)
        self.E = float(E)
        self.time = float(time)
        self.wgt = float(wgt)
        self.delayed_group = delayed_group
        self.surf_id = surf_id
        self.particle = particle

    @property
    def particle(self) -> ParticleType:
        return self._particle

    @particle.setter
    def particle(self, particle):
        self._particle = ParticleType(particle)

    def __repr__(self):
        return f'<SourceParticle: {str(self.particle)} at E={self.E:.6e} eV>'

    def to_tuple(self) -> tuple:
        """Return source particle attributes as a tuple

        Returns
        -------
        tuple
            Source particle attributes

        """
        return (self.r, self.u, self.E, self.time, self.wgt,
                self.delayed_group, self.surf_id, self.particle.pdg_number)


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
    cv.check_iterable_type(
        "source particles", source_particles, SourceParticle)
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
                particle_type = ParticleType(particle.pdgcode)

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
                 str(sp.particle)) for sp in self]

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
            fh.attrs['version'] = np.array([_VERSION_STATEPOINT, 2])
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
    """Read a collision track file in HDF5 format.

    Parameters
    ----------
    filename : str or path-like
        Path to the HDF5 collision track file.

    Returns
    -------
    numpy.ndarray
        Structured array containing collision track data.

    See Also
    --------
    read_collision_track_mcpl
    read_collision_track_file
    """

    with h5py.File(filename, 'r') as file:
        data = file['collision_track_bank'][:]

    return data


def read_collision_track_mcpl(file_path):
    """Read a collision track file in MCPL format.

    Parameters
    ----------
    file_path : str or path-like
        Path to the MCPL collision track file.

    Returns
    -------
    numpy.ndarray
        Structured array of particle collision track information, including
        position, direction, energy, weight, reaction data, and identifiers.

    See Also
    --------
    read_collision_track_hdf5
    read_collision_track_file
    """
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
        if f'blob_{i}' in myfile.blobs:
            blob_data = myfile.blobs[f'blob_{i}']
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
            data['particle'].append(ParticleType(p.pdgcode))
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
    """Read a collision track file (HDF5 or MCPL) and return its data.

    Parameters
    ----------
    filename : str or path-like
        Path to the collision track file to read. Must end with
        ``.h5`` or ``.mcpl``.

    Returns
    -------
    numpy.ndarray
        Structured array containing collision track data.

    See Also
    --------
    read_collision_track_hdf5
    read_collision_track_mcpl
    """

    filename = Path(filename)
    if filename.suffix not in ('.h5', '.mcpl'):
        raise ValueError('Collision track file must have a .h5 or .mcpl extension.')

    if filename.suffix == '.h5':
        return read_collision_track_hdf5(filename)
    else:
        return read_collision_track_mcpl(filename)
