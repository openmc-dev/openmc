from collections import Iterable
from numbers import Real, Integral

import numpy as np

import openmc.checkvalue as cv
from openmc.stats import Tabular, Univariate, Discrete, Mixture
from .container import interpolation_scheme
from .angle_energy import AngleEnergy


class CorrelatedAngleEnergy(AngleEnergy):
    """Correlated angle-energy distribution

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
    mu : Iterable of Iterable of openmc.stats.Univariate
        Distribution of scattering cosine for each incoming/outgoing energy

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
    mu : Iterable of Iterable of openmc.stats.Univariate
        Distribution of scattering cosine for each incoming/outgoing energy

    """

    def __init__(self, breakpoints, interpolation, energy, energy_out, mu):
        super(CorrelatedAngleEnergy, self).__init__()
        self.breakpoints = breakpoints
        self.interpolation = interpolation
        self.energy = energy
        self.energy_out = energy_out
        self.mu = mu

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
    def mu(self):
        return self._mu

    @breakpoints.setter
    def breakpoints(self, breakpoints):
        cv.check_type('correlated angle-energy breakpoints', breakpoints,
                      Iterable, Integral)
        self._breakpoints = breakpoints

    @interpolation.setter
    def interpolation(self, interpolation):
        cv.check_type('correlated angle-energy interpolation', interpolation,
                      Iterable, Integral)
        self._interpolation = interpolation

    @energy.setter
    def energy(self, energy):
        cv.check_type('correlated angle-energy incoming energy', energy,
                      Iterable, Real)
        self._energy = energy

    @energy_out.setter
    def energy_out(self, energy_out):
        cv.check_type('correlated angle-energy outgoing energy', energy_out,
                      Iterable, Univariate)
        self._energy_out = energy_out

    @mu.setter
    def mu(self, mu):
        cv.check_iterable_type('correlated angle-energy outgoing cosine',
                               mu, Univariate, 2, 2)
        self._mu = mu

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['type'] = np.string_('correlated')

        dset = group.create_dataset('energy', data=self.energy)
        dset.attrs['interpolation'] = np.vstack((self.breakpoints,
                                                 self.interpolation))

        # Determine total number of (E,p) pairs and create array
        n_tuple = sum(len(d.x) for d in self.energy_out)
        eout = np.empty((5, n_tuple))

        # Make sure all mu data is tabular
        mu_tabular = []
        for i, mu_i in enumerate(self.mu):
            mu_tabular.append([mu_ij if isinstance(mu_ij, (Tabular, Discrete)) else
                               mu_ij.to_tabular() for mu_ij in mu_i])

        # Determine total number of (mu,p) points and create array
        n_tuple = sum(sum(len(mu_ij.x) for mu_ij in mu_i)
                      for mu_i in mu_tabular)
        mu = np.empty((3, n_tuple))

        # Create array for offsets
        offsets = np.empty(len(self.energy_out), dtype=int)
        interpolation = np.empty(len(self.energy_out), dtype=int)
        n_discrete_lines = np.empty(len(self.energy_out), dtype=int)
        offset_e = 0
        offset_mu = 0

        # Populate offsets and eout array
        for i, d in enumerate(self.energy_out):
            n = len(d)
            offsets[i] = offset_e

            if isinstance(d, Mixture):
                discrete, continuous = d.distribution
                n_discrete_lines[i] = m = len(discrete)
                interpolation[i] = 1 if continuous.interpolation == 'histogram' else 2
                eout[0, offset_e:offset_e+m] = discrete.x
                eout[1, offset_e:offset_e+m] = discrete.p
                eout[2, offset_e:offset_e+m] = discrete.c
                eout[0, offset_e+m:offset_e+n] = continuous.x
                eout[1, offset_e+m:offset_e+n] = continuous.p
                eout[2, offset_e+m:offset_e+n] = continuous.c
            else:
                if isinstance(d, Tabular):
                    n_discrete_lines[i] = 0
                    interpolation[i] = 1 if d.interpolation == 'histogram' else 2
                elif isinstance(d, Discrete):
                    n_discrete_lines[i] = n
                    interpolation[i] = 1
                eout[0, offset_e:offset_e+n] = d.x
                eout[1, offset_e:offset_e+n] = d.p
                eout[2, offset_e:offset_e+n] = d.c

            for j, mu_ij in enumerate(mu_tabular[i]):
                if isinstance(mu_ij, Discrete):
                    eout[3, offset_e+j] = 0
                else:
                    eout[3, offset_e+j] = 1 if mu_ij.interpolation == 'histogram' else 2
                eout[4, offset_e+j] = offset_mu

                n_mu = len(mu_ij)
                mu[0, offset_mu:offset_mu+n_mu] = mu_ij.x
                mu[1, offset_mu:offset_mu+n_mu] = mu_ij.p
                mu[2, offset_mu:offset_mu+n_mu] = mu_ij.c

                offset_mu += n_mu

            offset_e += n

        # Create dataset for outgoing energy distributions
        dset = group.create_dataset('energy_out', data=eout)

        # Write interpolation on outgoing energy as attribute
        dset.attrs['offsets'] = offsets
        dset.attrs['interpolation'] = interpolation
        dset.attrs['n_discrete_lines'] = n_discrete_lines

        # Create dataset for outgoing angle distributions
        group.create_dataset('mu', data=mu)

    @classmethod
    def from_hdf5(cls, group):
        """Generate correlated angle-energy distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.CorrelatedAngleEnergy
            Correlated angle-energy distribution

        """
        interp_data = group['energy'].attrs['interpolation']
        energy_breakpoints = interp_data[0,:]
        energy_interpolation = interp_data[1,:]
        energy = group['energy'].value

        offsets = group['energy_out'].attrs['offsets']
        interpolation = group['energy_out'].attrs['interpolation']
        n_discrete_lines = group['energy_out'].attrs['n_discrete_lines']
        dset_eout = group['energy_out'].value
        energy_out = []

        dset_mu = group['mu'].value
        mu = []

        n_energy = len(energy)
        for i in range(n_energy):
            # Determine length of outgoing energy distribution and number of
            # discrete lines
            offset_e = offsets[i]
            if i < n_energy - 1:
                n = offsets[i+1] - offset_e
            else:
                n = dset_eout.shape[1] - offset_e
            m = n_discrete_lines[i]

            # Create discrete distribution if lines are present
            if m > 0:
                x = dset_eout[0, offset_e:offset_e+m]
                p = dset_eout[1, offset_e:offset_e+m]
                eout_discrete = Discrete(x, p)
                eout_discrete.c = dset_eout[2, offset_e:offset_e+m]
                p_discrete = eout_discrete.c[-1]

            # Create continuous distribution
            if m < n:
                interp = interpolation_scheme[interpolation[i]]

                x = dset_eout[0, offset_e+m:offset_e+n]
                p = dset_eout[1, offset_e+m:offset_e+n]
                eout_continuous = Tabular(x, p, interp, ignore_negative=True)
                eout_continuous.c = dset_eout[2, offset_e+m:offset_e+n]

            # If both continuous and discrete are present, create a mixture
            # distribution
            if m == 0:
                eout_i = eout_continuous
            elif m == n:
                eout_i = eout_discrete
            else:
                eout_i = Mixture([p_discrete, 1. - p_discrete],
                                 [eout_discrete, eout_continuous])

            # Read angular distributions
            mu_i = []
            for j in range(n):
                # Determine interpolation scheme
                interp_code = int(dset_eout[3, offsets[i] + j])

                # Determine offset and length
                offset_mu = int(dset_eout[4, offsets[i] + j])
                if offsets[i] + j < dset_eout.shape[1] - 1:
                    n_mu = int(dset_eout[4, offsets[i] + j + 1]) - offset_mu
                else:
                    n_mu = dset_mu.shape[1] - offset_mu

                # Get data
                x = dset_mu[0, offset_mu:offset_mu+n_mu]
                p = dset_mu[1, offset_mu:offset_mu+n_mu]
                c = dset_mu[2, offset_mu:offset_mu+n_mu]

                if interp_code == 0:
                    mu_ij = Discrete(x, p)
                else:
                    mu_ij = Tabular(x, p, interpolation_scheme[interp_code],
                                    ignore_negative=True)
                mu_ij.c = c
                mu_i.append(mu_ij)

                offset_mu += n_mu

            energy_out.append(eout_i)
            mu.append(mu_i)

            j += n

        return cls(energy_breakpoints, energy_interpolation,
                   energy, energy_out, mu)
