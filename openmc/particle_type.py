from enum import IntEnum

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


class ParticleType(IntEnum):
    """Class representing a particle type based on PDG numbers."""
    NEUTRON = 2112
    PHOTON = 22
    ELECTRON = 11
    POSITRON = -11
    PROTON = 2212
    DEUTERON = 1000010020
    TRITON = 1000010030
    ALPHA = 1000020040

    @classmethod
    def _missing_(cls, value):
        if isinstance(value, str):
            try:
                int_value = cls._pdg_number_from_string(value)
            except (TypeError, ValueError):
                return None
        else:
            try:
                int_value = int(value)
            except (TypeError, ValueError):
                return None

        if int_value in _LEGACY_PARTICLE_INDEX:
            int_value = _LEGACY_PARTICLE_INDEX[int_value]

        member = cls._value2member_map_.get(int_value)
        if member is not None:
            return member

        obj = int.__new__(cls, int_value)
        obj._value_ = int_value
        obj._name_ = None
        return obj

    @staticmethod
    def _pdg_number_from_string(value: str) -> int:
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
        return f'<ParticleType: {str(self)} (PDG={self.value})>'

    def __str__(self) -> str:
        """Returns a canonical string representation of the particle type."""
        if self.value in _PDG_NAME:
            return _PDG_NAME[self.value]

        if (zam := self.zam) is not None:
            Z, A, m = zam
            if Z <= 0 or Z > max(ATOMIC_SYMBOL) or A <= 0 or A > 999:
                raise ValueError(f"Invalid nuclear PDG number: {self.value}")
            return gnds_name(Z, A, m)

        return f'pdg:{self.value}'

    @property
    def zam(self) -> tuple[int, int, int] | None:
        """Returns the (Z, A, m) tuple for nuclear particles."""
        if self.value < 1000000000:
            return None
        Z = (self.value // 10000) % 1000
        A = (self.value // 10) % 1000
        m = self.value % 10
        if Z <= 0 or A <= 0:
            return None
        else:
            return (Z, A, m)

    @property
    def is_nucleus(self) -> bool:
        return self.zam is not None
