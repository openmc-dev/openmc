from abc import ABC
import textwrap

import lxml.etree as ET
import numpy as np
from numpy.typing import ArrayLike

import openmc
import openmc.checkvalue as cv
from .mixin import IDManagerMixin
from .mesh import RegularMesh, RectilinearMesh


# Field types
_FIELD_TYPES = {}


class FieldBase(IDManagerMixin, ABC):
    """Abstract base class for fields defined on a geometric mesh.

    While there is no abstract method currently, deriving from ABC enforces
    the contract of not instantiating FieldBase objects directly.

    Parameters
    ----------
    mesh : openmc.RegularMesh or openmc.RectilinearMesh
        Spatial mesh on which the field is defined.

    values : array-like
        Field values. Shape must be (n_elements,) for scalar fields
        or (n_elements, n_components) for vector fields.
    mapping : {'nodal', 'cell'}
        How field values map to the mesh. Default is 'cell'.

    Attributes
    ----------
    mesh : openmc.MeshBase
        Spatial mesh associated with the field.
    values : numpy.ndarray
        Array of field values.
    mapping : {'nodal', 'cell'}
        How field values map to the mesh.

    """
    # Shared accross all FieldBase subclasses so that IDs are globally unique
    next_id = 1
    used_ids = set()

    _field_type = None
    _n_components = None
    _supported_meshes = (RegularMesh, RectilinearMesh)

    def __init__(self, mesh, values, mapping="cell", field_id=None, name=''):
        self.mesh = mesh
        self.mapping = mapping
        self.values = values
        self.id = field_id
        self.name = name

    def __init_subclass__(cls, **kwargs):
        """Auto-register subclasses that define _field_type to _FIELD_TYPES."""
        super().__init_subclass__(**kwargs)
        if cls._field_type is not None:
            _FIELD_TYPES[cls._field_type] = cls

    @property
    def mesh(self) -> RegularMesh | RectilinearMesh:
        return self._mesh

    @mesh.setter
    def mesh(self, mesh: RegularMesh | RectilinearMesh):
        cv.check_type('mesh', mesh, self._supported_meshes)
        self._mesh = mesh

    @property
    def values(self) -> np.ndarray:
        return self._values

    @values.setter
    def values(self, values: ArrayLike):
        if not hasattr(self, "_mapping"):
            raise AttributeError("Mapping must be set before values.")

        # Check that _n_components is defined
        if self._n_components is None:
            raise NotImplementedError(
                f"{type(self).__name__} must define _n_components."
            )

        # Check that mesh is set
        if not hasattr(self, "_mesh"):
            raise AttributeError("Mesh must be set before values.")

        values = np.asarray(values, dtype=float)

        # Determine expected range shape
        if self.mapping == "nodal":
            n = self._mesh.n_vertices
        elif self.mapping == "cell":
            n = self._mesh.n_elements

        if self._n_components == 1:
            expected = (n,)
        else:
            expected = (n, self._n_components)

        if values.shape != expected:
            raise ValueError(
                f"{type(self).__name__} values has shape {values.shape}, "
                f"expected {expected} based on mesh with {n} element(s) "
                f"and {self._n_components} component(s)."
            )

        self._validate_data(values)
        self._values = values

    @property
    def mapping(self) -> str:
        return self._mapping

    @mapping.setter
    def mapping(self, mapping: str):
        cv.check_value("mapping", mapping, ("nodal", "cell"))
        self._mapping = mapping

    def _validate_data(self, values):
        """Subclass-specific validation.

        Parameters
        ----------
        values : numpy.ndarray
            Data array to validate.

        """
        pass

    def to_xml_element(self):
        """Create an XML element representing this field.

        Returns
        -------
        element : xml.etree._Element
            XML element for this field.

        """
        element = ET.Element("field")
        element.set("type", self._field_type)
        element.set("id", str(self.id))

        if self.name:
            element.set("name", self.name)

        mesh_elem = ET.SubElement(element, "mesh")
        mesh_elem.text = str(self.mesh.id)

        # Mapping
        mapping_elem = ET.SubElement(element, "mapping")
        mapping_elem.text = self.mapping

        values_elem = ET.SubElement(element, "values")
        flat = self.values.flatten(order="C")

        # Assumes two XML block levels for the indent
        indent = 8 * ' '
        # f"{v:.15e}" instead of str(v) for consistency?
        values_elem.text = f"\n{indent}".join(
            textwrap.wrap(" ".join([str(v) for v in flat]), 80)
        )
        return element

    @classmethod
    def from_xml_element(cls, elem, meshes):
        """Construct a field from an XML element.

        Parameters
        ----------
        elem : xml.etree._Element
            XML <field> element to parse.
        meshes : dict
            Dictionary mapping mesh IDs (int) to mesh instances

        Returns
        -------
        Field
            Reconstructed field instance

        """
        field_type = elem.get("type")
        field_id = int(elem.get("id"))
        name = elem.get("name", "")

        if cls is FieldBase:
            if field_type not in _FIELD_TYPES:
                raise TypeError(
                    f"Unknown field type '{field_type}'. "
                    f"Known types: {list(_FIELD_TYPES.keys())}"
                )
            subcls = _FIELD_TYPES[field_type]
        else:
            subcls = cls

        mesh_id = int(elem.find("mesh").text)
        if mesh_id not in meshes:
            raise ValueError(f"Mesh with id={mesh_id} not found.")
        mesh = meshes[mesh_id]

        # Mapping
        mapping_elem = elem.find("mapping")
        mapping = mapping_elem.text if mapping_elem is not None else "cell"
        values_text = elem.find("values").text
        flat = np.array([float(v) for v in values_text.split()])

        if mapping == "nodal":
            n = mesh.n_vertices
        elif mapping == "cell":
            n = mesh.n_elements
        else:
            raise ValueError(f"Unknown mapping '{mapping}'")

        expected_size = n * subcls._n_components
        if flat.size != expected_size:
            raise ValueError(f"Expected {expected_size} values, got {flat.size}.")

        if subcls._n_components > 1:
            flat = flat.reshape((n, subcls._n_components), order="C")

        kwargs = {
            "mesh": mesh,
            "mapping": mapping,
            "values": flat,
            "field_id": field_id,
            "name": name
        }
        return subcls(**kwargs)

    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        # Only checks that mesh IDs are consistent
        return (
            self.mesh.id == other.mesh.id
            and np.array_equal(self.values, other.values)
            and self.mapping == other.mapping
        )

    def __repr__(self):
        repr_str = f"{type(self).__name__}(id={self.id}"
        if self.name:
            repr_str += f", name='{self.name}'"
        repr_str += f", mapping='{self.mapping}'"
        repr_str += f", mesh_id={self.mesh.id}"
        repr_str += f", n_elements={self.mesh.n_elements})"
        return repr_str


class TemperatureField(FieldBase):
    """Temperature field defined on a mesh.

    Parameters
    ----------
    mesh : Mesh
        Spatial mesh associated with the field.
    values : iterable of float
        Temperature values (K) for each mesh element. Must be strictly positive.
    field_id : int, optional
        Unique identifier for this field.
    name : str, optional
        User-defined name for the field.
    mapping : {'nodal', 'cell'}
        How field values map to the mesh. Default is 'cell'.

    """
    _field_type = "temperature"
    _n_components = 1

    def _validate_data(self, values):
        if np.any(values <= 0.0):
            raise ValueError(
                "All temperature values must be strictly positive (Kelvin). "
                f"Found minimum value: {values.min():.4e}"
            )


class VelocityField(FieldBase):
    """Velocity field defined on a mesh.

    Parameters
    ----------
    mesh : Mesh
        Spatial mesh associated with the field.
    values : iterable of float
        Velocity vectors (vx, vy, vz) in cm/s for each mesh element.
    field_id : int, optional
        Unique identifier for this field.
    name : str, optional
        User-defined name for the field.
    mapping : {'nodal', 'cell'}
        How field values map to the mesh. Default is 'cell'.

    """
    _field_type = "velocity"
    _n_components = 3

    @property
    def magnitude(self):
        """Velocity magnitude at each mesh element."""
        return np.linalg.norm(self.values, axis=1)

    def max_speed(self):
        """Return the maximum speed in the field (cm/s)."""
        return float(self.magnitude.max())
