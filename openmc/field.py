from abc import ABC, abstractmethod

import openmc


class ScalarField(ABC):
    """Scalar field defined on a geometric mesh.

    Attributes
    ----------
    mesh : Mesh
        Spatial mesh associated with the field
    values : iterable of float
        List of values associated with each mesh cell

    """
    def __init__(self, mesh, values):
        """Initialization.

        Parameters
        ----------
        mesh : Mesh
            Spatial mesh associated with the field
        values : iterable of float
            List of values associated with each mesh cell

        """
        self.mesh = mesh
        self.values = values



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
        """Initialization.

        Parameters
        ----------
        mesh : Mesh
            Spatial mesh associated with the field
        values : iterable of float
            List of values associated with each mesh cell

        """
        super().__init__(mesh, values)
