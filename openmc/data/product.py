from collections import Iterable
from numbers import Real
import sys

import numpy as np
from numpy.polynomial.polynomial import Polynomial

import openmc.checkvalue as cv
from .function import Tabulated1D
from .angle_energy import AngleEnergy

if sys.version_info[0] >= 3:
    basestring = str


class Product(object):
    """Secondary particle emitted in a nuclear reaction

    Parameters
    ----------
    particle : str, optional
        What particle the reaction product is. Defaults to 'neutron'.

    Attributes
    ----------
    applicability : Iterable of openmc.data.Tabulated1D
        Probability of sampling a given distribution for this product.
    decay_rate : float
        Decay rate in inverse seconds
    distribution : Iterable of openmc.data.AngleEnergy
        Distributions of energy and angle of product.
    emission_mode : {'prompt', 'delayed', 'total'}
        Indicate whether the particle is emitted immediately or whether it
        results from the decay of reaction product (e.g., neutron emitted from a
        delayed neutron precursor). A special value of 'total' is used when the
        yield represents particles from prompt and delayed sources.
    particle : str
        What particle the reaction product is.
    yield_ : float or openmc.data.Tabulated1D or numpy.polynomial.Polynomial
        Yield of secondary particle in the reaction.

    """

    def __init__(self, particle='neutron'):
        self.particle = particle
        self.decay_rate = 0.0
        self.emission_mode = 'prompt'
        self.distribution = []
        self.applicability = []
        self.yield_ = 1

    def __repr__(self):
        if isinstance(self.yield_, Real):
            return "<Product: {}, emission={}, yield={}>".format(
                self.particle, self.emission_mode, self.yield_)
        elif isinstance(self.yield_, Tabulated1D):
            if np.all(self.yield_.y == self.yield_.y[0]):
                return "<Product: {}, emission={}, yield={}>".format(
                    self.particle, self.emission_mode, self.yield_.y[0])
            else:
                return "<Product: {}, emission={}, yield=tabulated>".format(
                    self.particle, self.emission_mode)
        else:
            return "<Product: {}, emission={}, yield=polynomial>".format(
                self.particle, self.emission_mode)

    @property
    def applicability(self):
        return self._applicability

    @property
    def decay_rate(self):
        return self._decay_rate

    @property
    def distribution(self):
        return self._distribution

    @property
    def emission_mode(self):
        return self._emission_mode

    @property
    def particle(self):
        return self._particle

    @property
    def yield_(self):
        return self._yield

    @applicability.setter
    def applicability(self, applicability):
        cv.check_type('product distribution applicability', applicability,
                      Iterable, Tabulated1D)
        self._applicability = applicability

    @decay_rate.setter
    def decay_rate(self, decay_rate):
        cv.check_type('product decay rate', decay_rate, Real)
        cv.check_greater_than('product decay rate', decay_rate, 0.0, True)
        self._decay_rate = decay_rate

    @distribution.setter
    def distribution(self, distribution):
        cv.check_type('product angle-energy distribution', distribution,
                      Iterable, AngleEnergy)
        self._distribution = distribution

    @emission_mode.setter
    def emission_mode(self, emission_mode):
        cv.check_value('product emission mode', emission_mode,
                       ('prompt', 'delayed', 'total'))
        self._emission_mode = emission_mode

    @particle.setter
    def particle(self, particle):
        cv.check_type('product particle type', particle, basestring)
        self._particle = particle

    @yield_.setter
    def yield_(self, yield_):
        cv.check_type('product yield', yield_,
                      (Real, Tabulated1D, Polynomial))
        self._yield = yield_

    def to_hdf5(self, group):
        """Write product to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['particle'] = np.string_(self.particle)
        group.attrs['emission_mode'] = np.string_(self.emission_mode)
        if self.decay_rate > 0.0:
            group.attrs['decay_rate'] = self.decay_rate

        # Write yield
        if isinstance(self.yield_, Tabulated1D):
            self.yield_.to_hdf5(group, 'yield')
            dset = group['yield']
            dset.attrs['type'] = np.string_('tabulated')
        elif isinstance(self.yield_, Polynomial):
            dset = group.create_dataset('yield', data=self.yield_.coef)
            dset.attrs['type'] = np.string_('polynomial')
        else:
            dset = group.create_dataset('yield', data=float(self.yield_))
            dset.attrs['type'] = np.string_('constant')

        # Write applicability/distribution
        group.attrs['n_distribution'] = len(self.distribution)
        for i, d in enumerate(self.distribution):
            dgroup = group.create_group('distribution_{}'.format(i))
            if self.applicability:
                self.applicability[i].to_hdf5(dgroup, 'applicability')
            d.to_hdf5(dgroup)

    @classmethod
    def from_hdf5(cls, group):
        """Generate reaction product from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.Product
            Reaction product

        """
        particle = group.attrs['particle'].decode()
        p = cls(particle)

        p.emission_mode = group.attrs['emission_mode'].decode()
        if 'decay_rate' in group.attrs:
            p.decay_rate = group.attrs['decay_rate']

        # Read yield
        yield_type = group['yield'].attrs['type'].decode()
        if yield_type == 'constant':
            p.yield_ = group['yield'].value
        elif yield_type == 'polynomial':
            p.yield_ = Polynomial(group['yield'].value)
        elif yield_type == 'tabulated':
            p.yield_ = Tabulated1D.from_hdf5(group['yield'])

        # Read applicability/distribution
        n_distribution = group.attrs['n_distribution']
        distribution = []
        applicability = []
        for i in range(n_distribution):
            dgroup = group['distribution_{}'.format(i)]
            if 'applicability' in dgroup:
                applicability.append(Tabulated1D.from_hdf5(
                    dgroup['applicability']))
            distribution.append(AngleEnergy.from_hdf5(dgroup))

        p.distribution = distribution
        p.applicability = applicability

        return p
