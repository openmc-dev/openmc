from collections import Iterable
from numbers import Real

import numpy as np

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from openmc.stats import Univariate, Tabular, Uniform
from .function import INTERPOLATION_SCHEME


class AngleDistribution(EqualityMixin):
    """Angle distribution as a function of incoming energy

    Parameters
    ----------
    energy : Iterable of float
        Incoming energies at which distributions exist
    mu : Iterable of openmc.stats.Univariate
        Distribution of scattering cosines corresponding to each incoming energy

    Attributes
    ----------
    energy : Iterable of float
        Incoming energies at which distributions exist
    mu : Iterable of openmc.stats.Univariate
        Distribution of scattering cosines corresponding to each incoming energy

    """

    def __init__(self, energy, mu):
        super(AngleDistribution, self).__init__()
        self.energy = energy
        self.mu = mu

    @property
    def energy(self):
        return self._energy

    @property
    def mu(self):
        return self._mu

    @energy.setter
    def energy(self, energy):
        cv.check_type('angle distribution incoming energy', energy,
                      Iterable, Real)
        self._energy = energy

    @mu.setter
    def mu(self, mu):
        cv.check_type('angle distribution scattering cosines', mu,
                      Iterable, Univariate)
        self._mu = mu

    def to_hdf5(self, group):
        """Write angle distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        dset = group.create_dataset('energy', data=self.energy)

        # Make sure all data is tabular
        mu_tabular = [mu_i if isinstance(mu_i, Tabular) else
                      mu_i.to_tabular() for mu_i in self.mu]

        # Determine total number of (mu,p) pairs and create array
        n_pairs = sum([len(mu_i.x) for mu_i in mu_tabular])
        pairs = np.empty((3, n_pairs))

        # Create array for offsets
        offsets = np.empty(len(mu_tabular), dtype=int)
        interpolation = np.empty(len(mu_tabular), dtype=int)
        j = 0

        # Populate offsets and pairs array
        for i, mu_i in enumerate(mu_tabular):
            n = len(mu_i.x)
            offsets[i] = j
            interpolation[i] = 1 if mu_i.interpolation == 'histogram' else 2
            pairs[0, j:j+n] = mu_i.x
            pairs[1, j:j+n] = mu_i.p
            pairs[2, j:j+n] = mu_i.c
            j += n

        # Create dataset for distributions
        dset = group.create_dataset('mu', data=pairs)

        # Write interpolation as attribute
        dset.attrs['offsets'] = offsets
        dset.attrs['interpolation'] = interpolation

    @classmethod
    def from_hdf5(cls, group):
        """Generate angular distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.AngleDistribution
            Angular distribution

        """
        energy = group['energy'].value
        data = group['mu']
        offsets = data.attrs['offsets']
        interpolation = data.attrs['interpolation']

        mu = []
        n_energy = len(energy)
        for i in range(n_energy):
            # Determine length of outgoing energy distribution and number of
            # discrete lines
            j = offsets[i]
            if i < n_energy - 1:
                n = offsets[i+1] - j
            else:
                n = data.shape[1] - j

            interp = INTERPOLATION_SCHEME[interpolation[i]]
            mu_i = Tabular(data[0, j:j+n], data[1, j:j+n], interp)
            mu_i.c = data[2, j:j+n]

            mu.append(mu_i)

        return cls(energy, mu)

    @classmethod
    def from_ace(cls, ace, location_dist, location_start):
        """Generate an angular distribution from ACE data

        Parameters
        ----------
        ace : openmc.data.ace.Table
            ACE table to read from
        location_dist : int
            Index in the XSS array corresponding to the start of a block,
            e.g. JXS(9).
        location_start : int
            Index in the XSS array corresponding to the start of an angle
            distribution array

        Returns
        -------
        openmc.data.AngleDistribution
            Angular distribution

        """
        # Set starting index for angle distribution
        idx = location_dist + location_start - 1

        # Number of energies at which angular distributions are tabulated
        n_energies = int(ace.xss[idx])
        idx += 1

        # Incoming energy grid
        energy = ace.xss[idx:idx + n_energies]
        idx += n_energies

        # Read locations for angular distributions
        lc = ace.xss[idx:idx + n_energies].astype(int)
        idx += n_energies

        mu = []
        for i in range(n_energies):
            if lc[i] > 0:
                # Equiprobable 32 bin distribution
                idx = location_dist + abs(lc[i]) - 1
                cos = ace.xss[idx:idx + 33]
                pdf = np.zeros(33)
                pdf[:32] = 1.0/(32.0*np.diff(cos))
                cdf = np.linspace(0.0, 1.0, 33)

                mu_i = Tabular(cos, pdf, 'histogram', ignore_negative=True)
                mu_i.c = cdf
            elif lc[i] < 0:
                # Tabular angular distribution
                idx = location_dist + abs(lc[i]) - 1
                intt = int(ace.xss[idx])
                n_points = int(ace.xss[idx + 1])
                data = ace.xss[idx + 2:idx + 2 + 3*n_points]
                data.shape = (3, n_points)

                mu_i = Tabular(data[0], data[1], INTERPOLATION_SCHEME[intt])
                mu_i.c = data[2]
            else:
                # Isotropic angular distribution
                mu_i = Uniform(-1., 1.)

            mu.append(mu_i)

        return cls(energy, mu)
