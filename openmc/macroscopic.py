import warnings

from openmc.checkvalue import check_type


class Macroscopic(str):
    """A Macroscopic object that can be used in a material.

    .. deprecated:: 0.16.1
            Use a str name of the MGXS object with :meth:`Material.add_macroscopic`.

    Parameters
    ----------
    name : str
        Name of the macroscopic data, e.g. UO2

    Attributes
    ----------
    name : str
        Name of the nuclide, e.g. UO2

    """

    def __new__(cls, name):
        warnings.warn(
            "The Macroscopic class is deprecated."
            "Use a str name of the MGXS object with :meth:`Material.add_macroscopic`.",
            category=DeprecationWarning, stacklevel=2
        )
        check_type('name', name, str)
        return super().__new__(cls, name)

    @property
    def name(self):
        return self
