from openmc.checkvalue import check_type


class Macroscopic(str):
    """A Macroscopic object that can be used in a material.

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
        check_type('name', name, str)
        return super().__new__(cls, name)

    @property
    def name(self):
        return self
