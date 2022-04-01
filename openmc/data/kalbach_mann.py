from collections.abc import Iterable
from numbers import Real, Integral
from warnings import warn

import numpy as np

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from openmc.stats import Tabular, Univariate, Discrete, Mixture
from .function import Tabulated1D, INTERPOLATION_SCHEME
from .angle_energy import AngleEnergy
from .data import EV_PER_MEV
from .endf import get_list_record, get_tab2_record


# Kalbach-Mann constants as defined in ENDF-6 manual BNL-203218-2018-INRE,
# Revision 215, File 6 description for LAW=1 and LANG=2.
_C1 = 0.04  # [1/MeV]
_C2 = 1.8E-6  # [1/MeV^3]
_C3 = 6.7E-7  # [1/MeV^4]
_ET1 = 130.  # [MeV]
_ET3 = 41.  # [MeV]
_M_NEUTRON = 1.
_M_PROTON = 1.
_M_DEUTERON = 1.
_M_TRITON = None
_M_3HE = None
_M_ALPHA = 0.
_SM_NEUTRON = 1/2.
_SM_PROTON = 1.
_SM_DEUTERON = 1.
_SM_TRITON = 1.
_SM_3HE = 1.
_SM_ALPHA = 2.

# Kalbach-Mann M coefficients
_TABULATED_PARTICLE_M = {
    1: _M_NEUTRON,
    1001: _M_PROTON,
    1002: _M_DEUTERON,
    1003: _M_TRITON,
    2003: _M_3HE,
    2004: _M_ALPHA
}

# Kalbach-Mann m coefficients
_TABULATED_PARTICLE_SM = {
    1: _SM_NEUTRON,
    1001: _SM_PROTON,
    1002: _SM_DEUTERON,
    1003: _SM_TRITON,
    2003: _SM_3HE,
    2004: _SM_ALPHA
}

# Breaking energy as defined in ENDF-6 manual BNL-203218-2018-INRE,
# Revision 215, Appendix H, Table 3.
_BREAKING_ENERGY_NEUTRON = 0.
_BREAKING_ENERGY_PROTON = 0.
_BREAKING_ENERGY_DEUTERON = 2.224566  # [MeV]
_BREAKING_ENERGY_TRITON = 8.481798  # [MeV]
_BREAKING_ENERGY_3HE = 7.718043  # [MeV]
_BREAKING_ENERGY_ALPHA = 28.29566  # [MeV]

_TABULATED_BREAKING_ENERGY = {
    1: _BREAKING_ENERGY_NEUTRON,
    1001: _BREAKING_ENERGY_PROTON,
    1002: _BREAKING_ENERGY_DEUTERON,
    1003: _BREAKING_ENERGY_TRITON,
    2003: _BREAKING_ENERGY_3HE,
    2004: _BREAKING_ENERGY_ALPHA
}

# Abundant IZA translation in merged library
_IZA_TRANSLATION = {
    6000: 6012,
}


class AtomicRepresentation(EqualityMixin):
    """Atomic representation of an isotope or a particle.

    Parameters
    ----------
    z: int
        Number of protons (atomic number)
    a: int
        Number of nucleons (mass number)

    Raises
    ------
    IOError:
        When the number of protons (z) declared is higher than the number
        of nucleons (a)

    Attributes
    ----------
    z: int
        Number of protons (atomic number)
    a: int
        Number of nucleons (mass number)
    n: int
        Number of neutrons
    breaking_energy: float
        Energy required to break the isotope or particle into their
        constituent nucleons from tabulated values
    M: float
        Kalbach-Mann M coefficient
    m: float
        Kalbach-Mann m coefficient
    iza: int
        ZA identifier defined as:
        iza = Z x 1000 + A,
        where Z is the number of protons and A the number of nucleons

    """
    def __init__(self, z, a):
        self._consistency_check(z, a)
        self._z = z
        self._a = a
        self.z = self._z
        self.a = self._a

    def __add__(self, other):
        """Adds two AtomicRepresentations.

        """
        z = self.z + other.z
        a = self.a + other.a
        return AtomicRepresentation(z=z, a=a)

    def __sub__(self, other):
        """Substracts two AtomicRepresentations.

        """
        z = self.z - other.z
        a = self.a - other.a
        return AtomicRepresentation(z=z, a=a)

    @property
    def a(self):
        self._consistency_check(self._z, self._a)
        return self._a

    @property
    def z(self):
        self._consistency_check(self._z, self._a)
        return self._z

    @property
    def n(self):
        return self.a - self.z

    @property
    def breaking_energy(self):
        breaking_energy = None
        if self.iza in _TABULATED_BREAKING_ENERGY:
            breaking_energy = _TABULATED_BREAKING_ENERGY[self.iza]
        return breaking_energy

    @property
    def M(self):
        M = None
        if self.iza in _TABULATED_PARTICLE_M:
            M = _TABULATED_PARTICLE_M[self.iza]
        return M

    @property
    def m(self):
        m = None
        if self.iza in _TABULATED_PARTICLE_SM:
            m = _TABULATED_PARTICLE_SM[self.iza]
        return m

    @property
    def iza(self):
        iza = self.z * 1000 + self.a
        return iza

    @a.setter
    def a(self, an):
        cv.check_type('a', an, Integral)
        cv.check_greater_than('a', an, 0, equality=True)
        self._consistency_check(self._z, an)
        self._a = an

    @z.setter
    def z(self, zn):
        cv.check_type('z', zn, Integral)
        cv.check_greater_than('z', zn, 0, equality=True)
        self._consistency_check(zn, self._a)
        self._z = zn

    @staticmethod
    def _consistency_check(z, a):
        """Simple consistency check.

        Parameters
        ----------
        z: int
            Number of protons (atomic number)
        a: int
            Number of nucleons (mass number)

        Raises
        ------
        IOError:
            When the number of protons (z) declared is higher than the number
            of nucleons (a)

        """
        if z > a:
            raise IOError(
                "Number of protons (%i) incompatible with number of "
                "nucleons (%i)" % (z, a)
            )

    @classmethod
    def from_iza(cls, iza):
        """Instantiates an AtomicRepresentation from a ZA identifier.

        The ZA identifier is defined as:
        iza = Z x 1000 + A,
        where Z is the number of protons and A the number of nucleons.

        Parameters
        ----------
        iza: int
            ZA identifier

        Returns
        -------
        AtomicRepresentation
            Atomic representation of the isotope/particle

        """
        if iza in _IZA_TRANSLATION.keys():
            iza = _IZA_TRANSLATION[iza]

        z = int(iza/1000)
        a = np.mod(iza, 1000)

        return cls(z, a)


def _calculate_separation_energy(compound, nucleus, particle):
    """Calculates the separation energy as defined in ENDF-6 manual
    BNL-203218-2018-INRE, Revision 215, File 6 description for LAW=1
    and LANG=2. This function can be used for the incident or emitted
    particle of the following reaction: A + a -> C -> B + b

    Parameters
    ----------
    compound: AtomicRepresentation
        Atomic representation of the compound (C)
    nucleus: AtomicRepresentation
        Atomic representation of the nucleus (A or B)
    particle: AtomicRepresentation
        Atomic representation of the particle (a or b)

    Returns
    -------
    separation_energy: float
        Separation energy in MeV

    """
    coef_1 = 15.68 * (compound.a - nucleus.a)
    coef_2 = 28.07 * ((compound.n - compound.z)**2 / float(compound.a) - \
        (nucleus.n - nucleus.z)**2 / float(nucleus.a))
    coef_3 = 18.56 * (compound.a**(2./3.) - nucleus.a**(2./3.))
    coef_4 = 33.22 * ((compound.n - compound.z)**2 / float(compound.a)**(4./3.) - \
        (nucleus.n - nucleus.z)**2 / float(nucleus.a)**(4./3.))
    coef_5 = 0.717 * (compound.z**2 / float(compound.a)**(1./3.) - \
        nucleus.z**2 / float(nucleus.a)**(1./3.))
    coef_6 = 1.211 * (compound.z**2 / float(compound.a) - \
        nucleus.z**2 / float(nucleus.a))

    separation_energy = coef_1 - coef_2 - coef_3 + coef_4 \
        - coef_5 + coef_6 - particle.breaking_energy

    return separation_energy


def _return_entrance_channel_energy(e_p, awr_t, awr_p):
    """Returns the entrance channel energy as defined in ENDF-6 manual
    BNL-203218-2018-INRE, Revision 215, File 6 description for LAW=1
    and LANG=2.

    Parameters
    ----------
    e_p: float
        Energy of the incident projectile in the laboratory system in eV
    awr_t: float
        Atomic weight ratio of the target
    awr_p: float
        Atomic weight ratio of the projectile

    Returns
    -------
    epsilon_p: float
        Entrance channel energy in eV

    """
    epsilon_p = e_p * awr_t / (awr_t + awr_p)
    return epsilon_p


def _return_emission_channel_energy(e_e, awr_r, awr_e):
    """Returns the emission channel energy as defined in ENDF-6 manual
    BNL-203218-2018-INRE, Revision 215, File 6 description for LAW=1
    and LANG=2.

    Parameters
    ----------
    e_e: float
        Energy of the emitted particle in the center of mass system in eV
    awr_r: float
        Atomic weight ratio of the residual nucleus
    awr_e: float
        Atomic weight ratio of the emitted particle

    Returns
    -------
    epsilon_e: float
        Emission channel energy in eV

    """
    epsilon_e = e_e * (awr_r + awr_e) / awr_r
    return epsilon_e


def _calculate_kalbach_slope(energy_projectile,
                             energy_emitted,
                             projectile,
                             target,
                             compound,
                             emitted,
                             residual):
    """Calculate the Kalbach slope for projectiles other than photons
    as defined in ENDF-6 manual BNL-203218-2018-INRE, Revision 215,
    File 6 description for LAW=1 and LANG=2.

    The entrance and emission channel energies are not calculated with
    the AWR number, but approximated with the number of mass instead.

    Parameters
    ----------
    energy_projectile: float
        Energy of the projectile in the laboratory system in eV
    energy_emitted: float
        Energy of the emitted particle in the center of mass system in eV
    projectile: AtomicRepresentation
        Atomic representation of the projectile
    target: AtomicRepresentation
        Atomic representation of the target
    compound: AtomicRepresentation
        Atomic representation of the compound
    emitted: AtomicRepresentation
        Atomic representation of the emitted particle
    residual: AtomicRepresentation
        Atomic representation of the residual nucleus

    Returns
    -------
    slope: float
        Kalbach-Mann slope

    """
    epsilon_a = _return_entrance_channel_energy(
        energy_projectile,
        target.a,
        projectile.a
    ) / EV_PER_MEV
    epsilon_b = _return_emission_channel_energy(
        energy_emitted,
        residual.a,
        emitted.a
    ) / EV_PER_MEV

    s_a = _calculate_separation_energy(compound, target, projectile)
    s_b = _calculate_separation_energy(compound, residual, emitted)

    e_a = epsilon_a + s_a
    e_b = epsilon_b + s_b

    r_1 = min(e_a, _ET1)
    r_3 = min(e_a, _ET3)

    x_1 = r_1 * e_b / e_a
    x_3 = r_3 * e_b / e_a

    slope = _C1 * x_1 \
        + _C2 * x_1**3 \
        + _C3 * projectile.M * emitted.m * x_3**4

    return slope


def return_kalbach_slope(energy_projectile,
                         energy_emitted,
                         iza_projectile,
                         iza_emitted,
                         iza_target):
    """Returns Kalbach-Mann slope from calculations.
    
    The associated reaction is defined as:
    A + a -> C -> B + b

    Where:

    - A is the targeted nucleus,
    - a is the projectile,
    - C is the compound,
    - B is the residual nucleus,
    - b is the emitted particle.

    This function uses the concept of ZA identifier defined as:
    iza = Z x 1000 + A,
    where Z is the number of protons and A the number of nucleons.

    The Kalbach-Mann slope calculation is done as defined in ENDF-6 manual
    BNL-203218-2018-INRE, Revision 215, File 6 description for LAW=1 and
    LANG=2. One exception to this, is that the entrance and emission channel
    energies are not calculated with the AWR number, but approximated with
    the number of mass instead.

    Parameters
    ----------
    energy_projectile: float
        Energy of the projectile in the laboratory system in eV
    energy_emitted: float
        Energy of the emitted particle in the center of mass system in eV
    iza_projectile: int
        ZA identifier of the projectile
    iza_emitted: int
        ZA identifier of the emitted particle
    iza_target: int
        ZA identifier of the targeted nucleus

    Raises
    ------
    NotImplementedError:
        When the ZA identifier of the projectile is not equal to 1
        (ie. other than a neutron).

    Returns
    -------
    slope: float
        Kalbach-Mann slope given with the same format as ACE file.

    """
    # TODO: develop for photons as projectile
    # TODO: test for other particles than neutron
    if iza_projectile != 1:
        raise NotImplementedError(
            "Developed and tested for neutron projectile only."
        )

    projectile = AtomicRepresentation.from_iza(iza_projectile)
    emitted = AtomicRepresentation.from_iza(iza_emitted)
    target = AtomicRepresentation.from_iza(iza_target)
    compound = projectile + target
    residual = compound - emitted

    slope = _calculate_kalbach_slope(
        energy_projectile,
        energy_emitted,
        projectile,
        target,
        compound,
        emitted,
        residual
    )

    return float("%7e" % slope)


class KalbachMann(AngleEnergy):
    """Kalbach-Mann distribution

    Parameters
    ----------
    breakpoints : Iterable of int
        Breakpoints defining interpolation regions
    interpolation : Iterable of int
        Interpolation codes
    energy : Iterable of float
        Incoming energies at which distributions exist
    energy_out : Iterable of openmc.stats.Univariate
        Distribution of outgoing energies corresponding to each incoming energy
    precompound : Iterable of openmc.data.Tabulated1D
        Precompound factor 'r' as a function of outgoing energy for each
        incoming energy
    slope : Iterable of openmc.data.Tabulated1D
        Kalbach-Chadwick angular distribution slope value 'a' as a function of
        outgoing energy for each incoming energy

    Attributes
    ----------
    breakpoints : Iterable of int
        Breakpoints defining interpolation regions
    interpolation : Iterable of int
        Interpolation codes
    energy : Iterable of float
        Incoming energies at which distributions exist
    energy_out : Iterable of openmc.stats.Univariate
        Distribution of outgoing energies corresponding to each incoming energy
    precompound : Iterable of openmc.data.Tabulated1D
        Precompound factor 'r' as a function of outgoing energy for each
        incoming energy
    slope : Iterable of openmc.data.Tabulated1D
        Kalbach-Chadwick angular distribution slope value 'a' as a function of
        outgoing energy for each incoming energy

    """

    def __init__(self, breakpoints, interpolation, energy, energy_out,
                 precompound, slope):
        super().__init__()
        self.breakpoints = breakpoints
        self.interpolation = interpolation
        self.energy = energy
        self.energy_out = energy_out
        self.precompound = precompound
        self.slope = slope

    @property
    def breakpoints(self):
        return self._breakpoints

    @property
    def interpolation(self):
        return self._interpolation

    @property
    def energy(self):
        return self._energy

    @property
    def energy_out(self):
        return self._energy_out

    @property
    def precompound(self):
        return self._precompound

    @property
    def slope(self):
        return self._slope

    @breakpoints.setter
    def breakpoints(self, breakpoints):
        cv.check_type('Kalbach-Mann breakpoints', breakpoints,
                      Iterable, Integral)
        self._breakpoints = breakpoints

    @interpolation.setter
    def interpolation(self, interpolation):
        cv.check_type('Kalbach-Mann interpolation', interpolation,
                      Iterable, Integral)
        self._interpolation = interpolation

    @energy.setter
    def energy(self, energy):
        cv.check_type('Kalbach-Mann incoming energy', energy,
                      Iterable, Real)
        self._energy = energy

    @energy_out.setter
    def energy_out(self, energy_out):
        cv.check_type('Kalbach-Mann distributions', energy_out,
                      Iterable, Univariate)
        self._energy_out = energy_out

    @precompound.setter
    def precompound(self, precompound):
        cv.check_type('Kalbach-Mann precompound factor', precompound,
                      Iterable, Tabulated1D)
        self._precompound = precompound

    @slope.setter
    def slope(self, slope):
        cv.check_type('Kalbach-Mann slope', slope, Iterable, Tabulated1D)
        self._slope = slope

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['type'] = np.string_('kalbach-mann')

        dset = group.create_dataset('energy', data=self.energy)
        dset.attrs['interpolation'] = np.vstack((self.breakpoints,
                                                 self.interpolation))

        # Determine total number of (E,p,r,a) tuples and create array
        n_tuple = sum(len(d) for d in self.energy_out)
        distribution = np.empty((5, n_tuple))

        # Create array for offsets
        offsets = np.empty(len(self.energy_out), dtype=int)
        interpolation = np.empty(len(self.energy_out), dtype=int)
        n_discrete_lines = np.empty(len(self.energy_out), dtype=int)
        j = 0

        # Populate offsets and distribution array
        for i, (eout, km_r, km_a) in enumerate(zip(
                self.energy_out, self.precompound, self.slope)):
            n = len(eout)
            offsets[i] = j

            if isinstance(eout, Mixture):
                discrete, continuous = eout.distribution
                n_discrete_lines[i] = m = len(discrete)
                interpolation[i] = 1 if continuous.interpolation == 'histogram' else 2
                distribution[0, j:j+m] = discrete.x
                distribution[1, j:j+m] = discrete.p
                distribution[2, j:j+m] = discrete.c
                distribution[0, j+m:j+n] = continuous.x
                distribution[1, j+m:j+n] = continuous.p
                distribution[2, j+m:j+n] = continuous.c
            else:
                if isinstance(eout, Tabular):
                    n_discrete_lines[i] = 0
                    interpolation[i] = 1 if eout.interpolation == 'histogram' else 2
                elif isinstance(eout, Discrete):
                    n_discrete_lines[i] = n
                    interpolation[i] = 1
                distribution[0, j:j+n] = eout.x
                distribution[1, j:j+n] = eout.p
                distribution[2, j:j+n] = eout.c

            distribution[3, j:j+n] = km_r.y
            distribution[4, j:j+n] = km_a.y
            j += n

        # Create dataset for distributions
        dset = group.create_dataset('distribution', data=distribution)

        # Write interpolation as attribute
        dset.attrs['offsets'] = offsets
        dset.attrs['interpolation'] = interpolation
        dset.attrs['n_discrete_lines'] = n_discrete_lines

    @classmethod
    def from_hdf5(cls, group):
        """Generate Kalbach-Mann distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.KalbachMann
            Kalbach-Mann energy distribution

        """
        interp_data = group['energy'].attrs['interpolation']
        energy_breakpoints = interp_data[0, :]
        energy_interpolation = interp_data[1, :]
        energy = group['energy'][()]

        data = group['distribution']
        offsets = data.attrs['offsets']
        interpolation = data.attrs['interpolation']
        n_discrete_lines = data.attrs['n_discrete_lines']

        energy_out = []
        precompound = []
        slope = []
        n_energy = len(energy)
        for i in range(n_energy):
            # Determine length of outgoing energy distribution and number of
            # discrete lines
            j = offsets[i]
            if i < n_energy - 1:
                n = offsets[i+1] - j
            else:
                n = data.shape[1] - j
            m = n_discrete_lines[i]

            # Create discrete distribution if lines are present
            if m > 0:
                eout_discrete = Discrete(data[0, j:j+m], data[1, j:j+m])
                eout_discrete.c = data[2, j:j+m]
                p_discrete = eout_discrete.c[-1]

            # Create continuous distribution
            if m < n:
                interp = INTERPOLATION_SCHEME[interpolation[i]]
                eout_continuous = Tabular(data[0, j+m:j+n], data[1, j+m:j+n], interp)
                eout_continuous.c = data[2, j+m:j+n]

            # If both continuous and discrete are present, create a mixture
            # distribution
            if m == 0:
                eout_i = eout_continuous
            elif m == n:
                eout_i = eout_discrete
            else:
                eout_i = Mixture([p_discrete, 1. - p_discrete],
                                 [eout_discrete, eout_continuous])

            # Precompound factor and slope are on rows 3 and 4, respectively
            km_r = Tabulated1D(data[0, j:j+n], data[3, j:j+n])
            km_a = Tabulated1D(data[0, j:j+n], data[4, j:j+n])

            energy_out.append(eout_i)
            precompound.append(km_r)
            slope.append(km_a)

        return cls(energy_breakpoints, energy_interpolation,
                   energy, energy_out, precompound, slope)

    @classmethod
    def from_ace(cls, ace, idx, ldis):
        """Generate Kalbach-Mann energy-angle distribution from ACE data

        Parameters
        ----------
        ace : openmc.data.ace.Table
            ACE table to read from
        idx : int
            Index in XSS array of the start of the energy distribution data
            (LDIS + LOCC - 1)
        ldis : int
            Index in XSS array of the start of the energy distribution block
            (e.g. JXS[11])

        Returns
        -------
        openmc.data.KalbachMann
            Kalbach-Mann energy-angle distribution

        """
        # Read number of interpolation regions and incoming energies
        n_regions = int(ace.xss[idx])
        n_energy_in = int(ace.xss[idx + 1 + 2*n_regions])

        # Get interpolation information
        idx += 1
        if n_regions > 0:
            breakpoints = ace.xss[idx:idx + n_regions].astype(int)
            interpolation = ace.xss[idx + n_regions:idx + 2*n_regions].astype(int)
        else:
            breakpoints = np.array([n_energy_in])
            interpolation = np.array([2])

        # Incoming energies at which distributions exist
        idx += 2*n_regions + 1
        energy = ace.xss[idx:idx + n_energy_in]*EV_PER_MEV

        # Location of distributions
        idx += n_energy_in
        loc_dist = ace.xss[idx:idx + n_energy_in].astype(int)

        # Initialize variables
        energy_out = []
        km_r = []
        km_a = []

        # Read each outgoing energy distribution
        for i in range(n_energy_in):
            idx = ldis + loc_dist[i] - 1

            # intt = interpolation scheme (1=hist, 2=lin-lin)
            INTTp = int(ace.xss[idx])
            intt = INTTp % 10
            n_discrete_lines = (INTTp - intt)//10
            if intt not in (1, 2):
                warn("Interpolation scheme for continuous tabular distribution "
                     "is not histogram or linear-linear.")
                intt = 2

            n_energy_out = int(ace.xss[idx + 1])
            data = ace.xss[idx + 2:idx + 2 + 5*n_energy_out].copy()
            data.shape = (5, n_energy_out)
            data[0, :] *= EV_PER_MEV

            # Create continuous distribution
            eout_continuous = Tabular(data[0][n_discrete_lines:],
                                      data[1][n_discrete_lines:]/EV_PER_MEV,
                                      INTERPOLATION_SCHEME[intt],
                                      ignore_negative=True)
            eout_continuous.c = data[2][n_discrete_lines:]
            if np.any(data[1][n_discrete_lines:] < 0.0):
                warn("Kalbach-Mann energy distribution has negative "
                     "probabilities.")

            # If discrete lines are present, create a mixture distribution
            if n_discrete_lines > 0:
                eout_discrete = Discrete(data[0][:n_discrete_lines],
                                         data[1][:n_discrete_lines])
                eout_discrete.c = data[2][:n_discrete_lines]
                if n_discrete_lines == n_energy_out:
                    eout_i = eout_discrete
                else:
                    p_discrete = min(sum(eout_discrete.p), 1.0)
                    eout_i = Mixture([p_discrete, 1. - p_discrete],
                                     [eout_discrete, eout_continuous])
            else:
                eout_i = eout_continuous

            energy_out.append(eout_i)
            km_r.append(Tabulated1D(data[0], data[3]))
            km_a.append(Tabulated1D(data[0], data[4]))

        return cls(breakpoints, interpolation, energy, energy_out, km_r, km_a)

    @classmethod
    def from_endf(cls, file_obj, iza_emitted, iza_target, projectile_mass):
        """Generate Kalbach-Mann distribution from an ENDF evaluation.

        If the projectile is a neutron, the slope is calculated when it is
        not given explicitly.

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of the Kalbach-Mann distribution
        iza_emitted : int
            ZA identifier of the emitted particle
        iza_target : int
            ZA identifier of the target
        projectile_mass: float
            Mass of the projectile

        Warns
        -----
        UserWarning
            If the mass of the projectile is not equal to 1 (other than
            a neutron), the slope is not calculated and set to 0 if missing.

        Returns
        -------
        openmc.data.KalbachMann
            Kalbach-Mann energy-angle distribution

        """
        params, tab2 = get_tab2_record(file_obj)
        lep = params[3]
        ne = params[5]
        energy = np.zeros(ne)
        n_discrete_energies = np.zeros(ne, dtype=int)
        energy_out = []
        precompound = []
        slope = []
        calculated_slope = []
        for i in range(ne):
            items, values = get_list_record(file_obj)
            energy[i] = items[1]
            n_discrete_energies[i] = items[2]
            # TODO: split out discrete energies
            n_angle = items[3]
            n_energy_out = items[5]
            values = np.asarray(values)
            values.shape = (n_energy_out, n_angle + 2)

            # Outgoing energy distribution at the i-th incoming energy
            eout_i = values[:, 0]
            eout_p_i = values[:, 1]
            energy_out_i = Tabular(eout_i, eout_p_i, INTERPOLATION_SCHEME[lep])
            energy_out.append(energy_out_i)

            # Precompound factors for Kalbach-Mann
            r_i = values[:, 2]

            # Slope factors for Kalbach-Mann
            if n_angle == 2:
                a_i = values[:, 3]
                calculated_slope.append(False)
            else:
                # Check if the projectile is not a neutron
                if not np.isclose(projectile_mass, 1.0, atol=1.0e-12, rtol=0.):
                    warn(
                        "Kalbach-Mann slope calculation is only available with "
                        "neutrons as projectile. Slope coefficients are set to 0."
                    )
                    a_i = np.zeros_like(r_i)
                    calculated_slope.append(False)

                else:
                    # TODO: retrieve IZA of the projectile
                    iza_projectile = 1
                    a_i = [return_kalbach_slope(energy_projectile=energy[i],
                                                energy_emitted=e,
                                                iza_projectile=iza_projectile,
                                                iza_emitted=iza_emitted,
                                                iza_target=iza_target)
                           for e in eout_i]
                    calculated_slope.append(True)

            precompound.append(Tabulated1D(eout_i, r_i))
            slope.append(Tabulated1D(eout_i, a_i))

        km_distribution = cls(tab2.breakpoints, tab2.interpolation, energy,
                              energy_out, precompound, slope)

        # List of bool to indicate slope calculation by OpenMC
        km_distribution._calculated_slope = calculated_slope

        return km_distribution
