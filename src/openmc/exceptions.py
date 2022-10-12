class OpenMCError(Exception):
    """Root exception class for OpenMC."""


class GeometryError(OpenMCError):
    """Geometry-related error"""


class InvalidIDError(OpenMCError):
    """Use of an ID that is invalid."""


class AllocationError(OpenMCError):
    """Error related to memory allocation."""


class OutOfBoundsError(OpenMCError):
    """Index in array out of bounds."""


class DataError(OpenMCError):
    """Error relating to nuclear data."""


class PhysicsError(OpenMCError):
    """Error relating to performing physics."""


class InvalidArgumentError(OpenMCError):
    """Argument passed was invalid."""


class InvalidTypeError(OpenMCError):
    """Tried to perform an operation on the wrong type."""


class SetupError(OpenMCError):
    """Error while setting up a problem."""
