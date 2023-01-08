from enum import IntEnum

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
            value (str): The string representation of the particle type.

        Returns
        -------
            ParticleType: The corresponding ParticleType instance.
        """
        try:
            return cls[value.upper()]
        except KeyError:
            raise ValueError(f"Invalid string for creation of {cls.__name__}: {value}")

    def to_string(self):
        """
        Returns a string representation of the ParticleType instance.

        Returns:
            str: The lowercase name of the ParticleType instance.
        """
        return self.name.lower()
