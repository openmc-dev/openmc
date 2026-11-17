from collections.abc import Iterable
from numbers import Real

import openmc
import openmc.checkvalue as cv


class ScalarField:
    """Scalar field defined on a geometric mesh.

    Attributes
    ----------
    mesh : Mesh
        Spatial mesh associated with the field
    values : iterable of float
        List of values associated with each mesh cell

    """
    def __init__(self, mesh, values):
        cv.check_type('values', values, Iterable)
        values = list(values)
        cv.check_type('values', values, Iterable, Real)

        # Check mesh compatibility
        if not isinstance(mesh, (openmc.RegularMesh, openmc.RectilinearMesh)):
            raise TypeError(
                f"{type(self)} only implemented for regular and "
                "rectilinear meshes.")

        # Check values/mesh size consistency
        if mesh.n_elements != len(values):
            raise ValueError(
                f"Inconsistent number of values '{len(values)}' compared to "
                f"the number of mesh cells '{mesh.n_elements}' declared in "
                f"an instance of {type(self)}.")

        # Assign
        self.mesh = mesh
        self.values = values

    def __eq__(self, other):
        if not isinstance(other, ScalarField):
            return NotImplemented
        return self.mesh.id == other.mesh.id and self.values == other.values


class TemperatureField(ScalarField):
    """Temperature field.

    Attributes
    ----------
    mesh : Mesh
        Spatial mesh associated with the field
    values : iterable of float
        List of values associated with each mesh cell

    """
    def __init__(self, mesh, values):
        super().__init__(mesh, values)

        # Check that values are non-negative
        for v in self.values:
            if v < 0:
                raise ValueError(
                    "Temperature values declared in the temperature field "
                    "must be non-negative.")
