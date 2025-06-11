from collections.abc import Iterable
from numbers import Real, Integral
from warnings import warn

import numpy as np

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from openmc.stats import Tabular, Univariate, Discrete, Mixture
from .function import Tabulated1D, INTERPOLATION_SCHEME
from .angle_energy import AngleEnergy
from .data import EV_PER_MEV, NEUTRON_MASS_EV
from .endf import get_list_record, get_tab2_record


class _AtomicRepresentation(EqualityMixin):
    """Atomic representation of an isotope or a particle.

    Parameters
    ----------
    z : int
        Number of protons (atomic number)
    a : int
        Number of nucleons (mass number)

    Raises
    ------
    ValueError
        When the number of protons (z) declared is higher than the number
        of nucleons (a)

    Attributes
    ----------
    z : int
        Number of protons (atomic number)
    a : int
        Number of nucleons (mass number)
    n : int
        Number of neutrons
    za : int
        ZA identifier, 1000*Z + A, where Z is the atomic number and A the mass
        number

    """
    def __init__(self, z, a):
        # Sanity checks on values
        cv.check_type('z', z, Integral)
        cv.check_greater_than('z', z, 0, equality=True)
        cv.check_type('a', a, Integral)
        cv.check_greater_than('a', a, 0, equality=True)
        if z > a:
            raise ValueError(f"Number of protons ({z}) must be less than or "
                             f"equal to number of nucleons ({a}).")

        self._z = z
        self._a = a

    def __add__(self, other):
        """Add two _AtomicRepresentations"""
        z = self.z + other.z
        a = self.a + other.a
        return _AtomicRepresentation(z=z, a=a)

    def __sub__(self, other):
        """Substract two _AtomicRepresentations"""
        z = self.z - other.z
        a = self.a - other.a
        return _AtomicRepresentation(z=z, a=a)

    @property
    def a(self):
        return self._a

    @property
    def z(self):
        return self._z

    @property
    def n(self):
        return self.a - self.z

    @property
    def za(self):
        return self.z * 1000 + self.a

    @classmethod
    def from_za(cls, za):
        """Instantiate an _AtomicRepresentation from a ZA identifier.

        Parameters
        ----------
        za : int
            ZA identifier, 1000*Z + A, where Z is the atomic number and A the
            mass number

        Returns
        -------
        _AtomicRepresentation
            Atomic representation of the isotope/particle

        """
        z, a = divmod(za, 1000)
        return cls(z, a)


def _separation_energy(compound, nucleus, particle):
    """Calculates the separation energy as defined in ENDF-6 manual
    BNL-203218-2018-INRE, Revision 215, File 6 description for LAW=1
    and LANG=2. This function can be used for the incident or emitted
    particle of the following reaction: A + a -> C -> B + b

    Parameters
    ----------
    compound : _AtomicRepresentation
        Atomic representation of the compound (C)
    nucleus : _AtomicRepresentation
        Atomic representation of the nucleus (A or B)
    particle : _AtomicRepresentation
        Atomic representation of the particle (a or b)

    Returns
    -------
    separation_energy : float
        Separation energy in MeV

    """
    # Determine A, Z, and N for compound and nucleus
    A_c = compound.a
    Z_c = compound.z
    N_c = compound.n
    A_a = nucleus.a
    Z_a = nucleus.z
    N_a = nucleus.n

    # Determine breakup energy of incident particle (ENDF-6 Formats Manual,
    # Appendix H, Table 3) in MeV
    za_to_breaking_energy = {
        1: 0.0,
        1001: 0.0,
        1002: 2.224566,
        1003: 8.481798,
        2003: 7.718043,
        2004: 28.29566
    }
    I_a = za_to_breaking_energy[particle.za]

    # Eq. 4 in in doi:10.1103/PhysRevC.37.2350 or ENDF-6 Formats Manual section
    # 6.2.3.2
    return (
        15.68 * (A_c - A_a) -
        28.07 * ((N_c - Z_c)**2 / A_c - (N_a - Z_a)**2 / A_a) -
        18.56 * (A_c**(2./3.) - A_a**(2./3.)) +
        33.22 * ((N_c - Z_c)**2 / A_c**(4./3.) - (N_a - Z_a)**2 / A_a**(4./3.)) -
        0.717 * (Z_c**2 / A_c**(1./3.) - Z_a**2 / A_a**(1./3.)) +
        1.211 * (Z_c**2 / A_c - Z_a**2 / A_a) -
        I_a
    )


def kalbach_slope(energy_projectile, energy_emitted, za_projectile,
                  za_emitted, za_target):
    """Returns Kalbach-Mann slope from calculations.

    The associated reaction is defined as:
    A + a -> C -> B + b

    Where:

    - A is the targeted nucleus,
    - a is the projectile,
    - C is the compound,
    - B is the residual nucleus,
    - b is the emitted particle.

    The Kalbach-Mann slope calculation is done as defined in ENDF-6 manual
    BNL-203218-2018-INRE, Revision 215, File 6 description for LAW=1 and
    LANG=2. One exception to this, is that the entrance and emission channel
    energies are not calculated with the AWR number, but approximated with
    the number of mass instead.

    .. versionadded:: 0.13.1

    Parameters
    ----------
    energy_projectile : float
        Energy of the projectile in the laboratory system in eV
    energy_emitted : float
        Energy of the emitted particle in the center of mass system in eV
    za_projectile : int
        ZA identifier of the projectile
    za_emitted : int
        ZA identifier of the emitted particle
    za_target : int
        ZA identifier of the targeted nucleus

    Raises
    ------
    NotImplementedError
        When the projectile is not a neutron

    Returns
    -------
    slope : float
        Kalbach-Mann slope given with the same format as ACE file.

    """

    if za_projectile == 0:
        # Calculate slope for photons using Eq. 3 in doi:10.1080/18811248.1995.9731830
        # or ENDF-6 Formats Manual section 6.2.3.2
        slope_n = kalbach_slope(energy_projectile, energy_emitted, 1,
                  za_emitted, za_target)
        return slope_n * np.sqrt(0.5*energy_projectile/NEUTRON_MASS_EV)*np.minimum(4,np.maximum(1,9.3/np.sqrt(energy_emitted/EV_PER_MEV))) 
        
    # Special handling of elemental carbon
    if za_emitted == 6000:
        za_emitted = 6012
    if za_target == 6000:
        za_target = 6012

    projectile = _AtomicRepresentation.from_za(za_projectile)
    emitted = _AtomicRepresentation.from_za(za_emitted)
    target = _AtomicRepresentation.from_za(za_target)
    compound = projectile + target
    residual = compound - emitted

    # Calculate entrance and emission channel energy in MeV, defined in section
    # 6.2.3.2 in the ENDF-6 Formats Manual
    epsilon_a = energy_projectile * target.a / (target.a + projectile.a) / EV_PER_MEV
    epsilon_b = energy_emitted * (residual.a + emitted.a) \
        / (residual.a * EV_PER_MEV)

    # Calculate separation energies using Eq. 4 in doi:10.1103/PhysRevC.37.2350
    # or ENDF-6 Formats Manual section 6.2.3.2
    s_a = _separation_energy(compound, target, projectile)
    s_b = _separation_energy(compound, residual, emitted)

    # See Eq. 10 in doi:10.1103/PhysRevC.37.2350 or section 6.2.3.2 in the
    # ENDF-6 Formats Manual
    za_to_M = {1: 1.0, 1001: 1.0, 1002: 1.0, 2004: 0.0}
    za_to_m = {1: 0.5, 1001: 1.0, 1002: 1.0, 1003: 1.0, 2003: 1.0, 2004: 2.0}
    M = za_to_M[projectile.za]
    m = za_to_m[emitted.za]
    e_a = epsilon_a + s_a
    e_b = epsilon_b + s_b
    r_1 = min(e_a, 130.)
    r_3 = min(e_a, 41.)
    x_1 = r_1 * e_b / e_a
    x_3 = r_3 * e_b / e_a
    return 0.04 * x_1 + 1.8e-6 * x_1**3 + 6.7e-7 * M * m * x_3**4


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
    is_photon : bool 
        Whether projectile is a photon, defaults to False

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
    is_photon : bool 
        Whether projectile particle is a photon        

    """

    def __init__(self, breakpoints, interpolation, energy, energy_out,
                 precompound, slope, is_photon=False):
        super().__init__()
        self.breakpoints = breakpoints
        self.interpolation = interpolation
        self.is_photon = is_photon
        self.energy = energy
        self.energy_out = energy_out
        self.precompound = precompound
        self.slope = slope

    @property
    def breakpoints(self):
        return self._breakpoints

    @breakpoints.setter
    def breakpoints(self, breakpoints):
        cv.check_type('Kalbach-Mann breakpoints', breakpoints,
                      Iterable, Integral)
        self._breakpoints = breakpoints

    @property
    def interpolation(self):
        return self._interpolation

    @interpolation.setter
    def interpolation(self, interpolation):
        cv.check_type('Kalbach-Mann interpolation', interpolation,
                      Iterable, Integral)
        self._interpolation = interpolation
        
    @property
    def is_photon(self):
        return self._is_photon

    @is_photon.setter
    def is_photon(self, is_photon):
        cv.check_type('Kalbach-Mann is_photon', is_photon, bool)
        self._is_photon = is_photon     

    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self, energy):
        cv.check_type('Kalbach-Mann incoming energy', energy,
                      Iterable, Real)
        self._energy = energy

    @property
    def energy_out(self):
        return self._energy_out

    @energy_out.setter
    def energy_out(self, energy_out):
        cv.check_type('Kalbach-Mann distributions', energy_out,
                      Iterable, Univariate)
        self._energy_out = energy_out

    @property
    def precompound(self):
        return self._precompound

    @precompound.setter
    def precompound(self, precompound):
        cv.check_type('Kalbach-Mann precompound factor', precompound,
                      Iterable, Tabulated1D)
        self._precompound = precompound

    @property
    def slope(self):
        return self._slope

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
        group.attrs['type'] = np.bytes_('kalbach-mann')
        group.attrs['is_photon'] = self.is_photon

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
        is_photon = bool(group.attrs.get("is_photon", False))
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
                   energy, energy_out, precompound, slope, is_photon=is_photon)

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
        is_photon : bool
            Whether projectile is photon

        Returns
        -------
        openmc.data.KalbachMann
            Kalbach-Mann energy-angle distribution

        """
        is_photon = bool(ace.data_type.value == 'u')
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

        return cls(breakpoints, interpolation, energy, energy_out, km_r, km_a, is_photon=is_photon)

    @classmethod
    def from_endf(cls, file_obj, za_emitted, za_target, za_projectile):
        """Generate Kalbach-Mann distribution from an ENDF evaluation.

        If the projectile is a neutron, the slope is calculated when it is
        not given explicitly.

        .. versionchanged:: 0.13.1
            Arguments changed to accommodate slope calculation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of the Kalbach-Mann distribution
        za_emitted : int
            ZA identifier of the emitted particle
        za_target : int
            ZA identifier of the target
        za_projectile : int
            ZA identifier of the projectile

        Returns
        -------
        openmc.data.KalbachMann
            Kalbach-Mann energy-angle distribution

        """
        is_photon = bool(za_projectile==0)
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
                a_i = [kalbach_slope(energy_projectile=energy[i],
                                     energy_emitted=e,
                                     za_projectile=za_projectile,
                                     za_emitted=za_emitted,
                                     za_target=za_target)
                       for e in eout_i]
                calculated_slope.append(True)

            precompound.append(Tabulated1D(eout_i, r_i))
            slope.append(Tabulated1D(eout_i, a_i))

        km_distribution = cls(tab2.breakpoints, tab2.interpolation, energy,
                              energy_out, precompound, slope, is_photon)

        # List of bool to indicate slope calculation by OpenMC
        km_distribution._calculated_slope = calculated_slope

        return km_distribution
