from numbers import Integral

from openmc.data import gnds_name, zam, ATOMIC_SYMBOL


_PDG_NAME = {
    2112: 'neutron',
    22: 'photon',
    11: 'electron',
    -11: 'positron',
    2212: 'H1',
}

_ALIAS_PDG = {
    'neutron': 2112,
    'n': 2112,
    'photon': 22,
    'gamma': 22,
    'electron': 11,
    'positron': -11,
    'proton': 2212,
    'p': 2212,
    'h1': 2212,
    'deuteron': 1000010020,
    'd': 1000010020,
    'h2': 1000010020,
    'triton': 1000010030,
    't': 1000010030,
    'h3': 1000010030,
    'alpha': 1000020040,
    'he4': 1000020040,
}

_LEGACY_PARTICLE_INDEX = {
    0: 2112,
    1: 22,
    2: 11,
    3: -11,
}


class ParticleType:
    """Particle type defined by a PDG number.

    ParticleType uses the Particle Data Group (PDG) Monte Carlo numbering scheme
    to uniquely identify particle types. This includes elementary particles
    (neutrons, photons, etc.) and nuclear codes for isotopes.

    Parameters
    ----------
    value : str, int, or ParticleType
        The particle identifier. Can be:

        - A string name (e.g., 'neutron', 'photon', 'He4', 'U235')
        - An integer PDG number (e.g., 2112 for neutron)
        - A string with PDG prefix (e.g., 'pdg:2112')
        - An existing ParticleType instance

    Attributes
    ----------
    pdg_number : int
        The PDG number for this particle type
    zam : tuple of int or None
        For nuclear particles, the (Z, A, m) tuple where Z is atomic number,
        A is mass number, and m is metastable state. None for elementary particles.
    is_nucleus : bool
        Whether this particle is a nucleus (ion)

    Examples
    --------
    >>> neutron = ParticleType('neutron')
    >>> neutron.pdg_number
    2112
    >>> he4 = ParticleType('He4')
    >>> he4.zam
    (2, 4, 0)
    >>> ParticleType(2112) == ParticleType('neutron')
    True

    """

    __slots__ = ('_pdg_number',)

    def __init__(self, value: 'str | int | ParticleType'):
        if isinstance(value, ParticleType):
            pdg = value._pdg_number
        elif isinstance(value, str):
            pdg = self._pdg_number_from_string(value)
        elif isinstance(value, Integral):
            pdg = int(value)
            # Handle legacy particle indices (0, 1, 2, 3)
            if pdg in _LEGACY_PARTICLE_INDEX:
                pdg = _LEGACY_PARTICLE_INDEX[pdg]
        else:
            raise TypeError(f"Cannot create ParticleType from {type(value).__name__}")

        self._pdg_number = pdg

    def __eq__(self, other):
        if isinstance(other, ParticleType):
            return self._pdg_number == other._pdg_number
        if isinstance(other, Integral):
            return self._pdg_number == int(other)
        if isinstance(other, str):
            try:
                return self._pdg_number == ParticleType(other)._pdg_number
            except (ValueError, TypeError):
                return False
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self._pdg_number)

    def __int__(self) -> int:
        return self._pdg_number

    @property
    def pdg_number(self) -> int:
        return self._pdg_number

    @staticmethod
    def _pdg_number_from_string(value: str) -> int:
        """Parse a string to get a PDG number.

        Parameters
        ----------
        value : str
            Particle identifier string

        Returns
        -------
        int
            PDG number

        Raises
        ------
        ValueError
            If string cannot be parsed as a valid particle identifier

        """
        s = value.strip()
        if not s:
            raise ValueError('Particle identifier cannot be empty.')

        lower = s.lower()
        if lower.startswith('pdg:'):
            code_str = lower[4:]
            try:
                return int(code_str)
            except ValueError:
                raise ValueError(f'Invalid PDG number: {code_str}')

        if lower in _ALIAS_PDG:
            return _ALIAS_PDG[lower]

        # Assume it is a GNDS nuclide name
        Z, A, m = zam(s)
        if Z <= 0 or Z > 999 or A <= 0 or A > 999 or m < 0 or m > 9:
            raise ValueError('Invalid Z/A/m for nuclear PDG number.')
        return 1000000000 + Z * 10000 + A * 10 + m

    def __repr__(self) -> str:
        return f'<ParticleType: {str(self)} (PDG={self._pdg_number})>'

    def __str__(self) -> str:
        """Return a canonical string representation of the particle type.

        Returns
        -------
        str
            Canonical name (e.g., 'neutron', 'He4', 'pdg:12345')

        """
        if self._pdg_number in _PDG_NAME:
            return _PDG_NAME[self._pdg_number]

        if (zam_tuple := self.zam) is not None:
            Z, A, m = zam_tuple
            if Z <= 0 or Z > max(ATOMIC_SYMBOL) or A <= 0 or A > 999:
                raise ValueError(f"Invalid nuclear PDG number: {self._pdg_number}")
            return gnds_name(Z, A, m)

        return f'pdg:{self._pdg_number}'

    @property
    def zam(self) -> 'tuple[int, int, int] | None':
        """Return the (Z, A, m) tuple for nuclear particles.

        Returns
        -------
        tuple of int or None
            For nuclear particles, returns (Z, A, m) where Z is atomic number,
            A is mass number, and m is metastable state. Returns None for
            elementary particles.

        """
        if self._pdg_number < 1000000000:
            return None
        Z = (self._pdg_number // 10000) % 1000
        A = (self._pdg_number // 10) % 1000
        m = self._pdg_number % 10
        if Z <= 0 or A <= 0:
            return None
        else:
            return (Z, A, m)

    @property
    def is_nucleus(self) -> bool:
        """Return whether this particle is a nucleus.

        Returns
        -------
        bool
            True if the particle is a nucleus (ion), False otherwise

        """
        return self.zam is not None


# Define common particle constants
ParticleType.NEUTRON = ParticleType(2112)
ParticleType.PHOTON = ParticleType(22)
ParticleType.ELECTRON = ParticleType(11)
ParticleType.POSITRON = ParticleType(-11)
ParticleType.PROTON = ParticleType(2212)
ParticleType.DEUTERON = ParticleType(1000010020)
ParticleType.TRITON = ParticleType(1000010030)
ParticleType.ALPHA = ParticleType(1000020040)
