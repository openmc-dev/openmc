import warnings


class Nuclide(str):
    """A nuclide that can be used in a material.

    Parameters
    ----------
    name : str
        Name of the nuclide, e.g. 'U235'

    Attributes
    ----------
    name : str
        Name of the nuclide, e.g. 'U235'

    """

    def __new__(cls, name):
        # Initialize class attributes
        orig_name = name

        if '-' in name:
            name = name.replace('-', '')
            name = name.replace('Nat', '0')
            if name.endswith('m'):
                name = name[:-1] + '_m1'

            msg = ('OpenMC nuclides follow the GNDS naming convention. '
                   f'Nuclide "{orig_name}" is being renamed as "{name}".')
            warnings.warn(msg)

        return super().__new__(cls, name)

    @property
    def name(self):
        return self
