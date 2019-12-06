from collections.abc import Iterable
from numbers import Real, Integral
from warnings import warn

import numpy as np

import openmc.checkvalue as cv
from openmc.stats import Tabular, Univariate, Discrete, Mixture
from .function import Tabulated1D, INTERPOLATION_SCHEME
from .angle_energy import AngleEnergy
from .data import EV_PER_MEV
from .endf import get_list_record, get_tab2_record


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
            data[0,:] *= EV_PER_MEV

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
    def from_endf(cls, file_obj):
        """Generate Kalbach-Mann distribution from an ENDF evaluation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of the Kalbach-Mann distribution

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
            eout_i = values[:,0]
            eout_p_i = values[:,1]
            energy_out_i = Tabular(eout_i, eout_p_i, INTERPOLATION_SCHEME[lep])
            energy_out.append(energy_out_i)

            # Precompound and slope factors for Kalbach-Mann
            r_i = values[:,2]
            if n_angle == 2:
                a_i = values[:,3]
            else:
                a_i = np.zeros_like(r_i)
            precompound.append(Tabulated1D(eout_i, r_i))
            slope.append(Tabulated1D(eout_i, a_i))

        return cls(tab2.breakpoints, tab2.interpolation, energy,
                   energy_out, precompound, slope)
