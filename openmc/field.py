import openmc


class ScalarField():
    """Scalar field defined on a geometric mesh.

    Attributes
    ----------
    mesh : Mesh
        Spatial mesh associated with the field
    values : iterable of float
        List of values associated with each mesh cell

    """
    def __init__(self, mesh, values):
        # Check mesh compatibility
        compatible_mesh_types = (openmc.RegularMesh, openmc.RectilinearMesh)
        compatible = False
        for t in compatible_mesh_types:
            if isinstance(mesh, t):
                compatible = True
        if not compatible:
            raise NotImplementedError(
                f"{type(self)} only implemented for regular and rectilinear meshes.")

        # Check values/mesh size consistency
        if mesh.n_elements != len(values):
            raise ValueError(
                f"Inconsistent number of values '{len(values)}' compared to the number "
                f"of mesh cells '{mesh.n_elements}' declared in an instance of "
                f"{type(self)}.")

        # Assign
        self.mesh = mesh
        self.values = values

    def __eq__(self, other):
        if not isinstance(other, ScalarField):
            return NotImplementedError()
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

        # Check that values are positive
        for v in values:
            if v < 0:
                raise ValueError(
                    "Temperature values declared in the temperature field must be "
                    "positive.")
