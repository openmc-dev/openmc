from collections.abc import Iterable
from numbers import Real, Integral

import numpy as np

import openmc.checkvalue as cv
from openmc.stats import Tabular, Univariate, Discrete
from .angle_energy import AngleEnergy
from .endf import get_tab2_record, get_tab1_record


class LaboratoryAngleEnergy(AngleEnergy):
    """Laboratory angle-energy distribution

    Parameters
    ----------
    breakpoints : Iterable of int
        Breakpoints defining interpolation regions
    interpolation : Iterable of int
        Interpolation codes
    energy : Iterable of float
        Incoming energies at which distributions exist
    mu : Iterable of openmc.stats.Univariate
        Distribution of scattering cosines for each incoming energy
    energy_out : Iterable of Iterable of openmc.stats.Univariate
        Distribution of outgoing energies for each incoming energy/scattering
        cosine

    Attributes
    ----------
    breakpoints : Iterable of int
        Breakpoints defining interpolation regions
    interpolation : Iterable of int
        Interpolation codes
    energy : Iterable of float
        Incoming energies at which distributions exist
    mu : Iterable of openmc.stats.Univariate
        Distribution of scattering cosines for each incoming energy
    energy_out : Iterable of Iterable of openmc.stats.Univariate
        Distribution of outgoing energies for each incoming energy/scattering
        cosine

    """

    def __init__(self, breakpoints, interpolation, energy, mu, energy_out):
        super().__init__()
        self.breakpoints = breakpoints
        self.interpolation = interpolation
        self.energy = energy
        self.mu = mu
        self.energy_out = energy_out

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
    def mu(self):
        return self._mu

    @property
    def energy_out(self):
        return self._energy_out

    @breakpoints.setter
    def breakpoints(self, breakpoints):
        cv.check_type('laboratory angle-energy breakpoints', breakpoints,
                      Iterable, Integral)
        self._breakpoints = breakpoints

    @interpolation.setter
    def interpolation(self, interpolation):
        cv.check_type('laboratory angle-energy interpolation', interpolation,
                      Iterable, Integral)
        self._interpolation = interpolation

    @energy.setter
    def energy(self, energy):
        cv.check_type('laboratory angle-energy incoming energy', energy,
                      Iterable, Real)
        self._energy = energy

    @mu.setter
    def mu(self, mu):
        cv.check_type('laboratory angle-energy outgoing cosine', mu,
                      Iterable, Univariate)
        self._mu = mu

    @energy_out.setter
    def energy_out(self, energy_out):
        cv.check_iterable_type('laboratory angle-energy outgoing energy',
                               energy_out, Univariate, 2, 2)
        self._energy_out = energy_out

    @classmethod
    def from_endf(cls, file_obj):
        """Generate laboratory angle-energy distribution from an ENDF evaluation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for a correlated
            angle-energy distribution

        Returns
        -------
        openmc.data.LaboratoryAngleEnergy
            Laboratory angle-energy distribution

        """
        params, tab2 = get_tab2_record(file_obj)
        ne = params[5]
        energy = np.zeros(ne)
        mu = []
        energy_out = []
        for i in range(ne):
            params, _ = get_tab2_record(file_obj)
            energy[i] = params[1]
            n_mu = params[5]
            mu_i = np.zeros(n_mu)
            p_mu_i = np.zeros(n_mu)
            energy_out_i = []
            for j in range(n_mu):
                params, f = get_tab1_record(file_obj)
                mu_i[j] = params[1]
                p_mu_i[j] = sum(f.y)
                energy_out_i.append(Tabular(f.x, f.y))
            mu.append(Tabular(mu_i, p_mu_i))
            energy_out.append(energy_out_i)

        return cls(tab2.breakpoints, tab2.interpolation, energy, mu, energy_out)

    def to_hdf5(self, group):
        raise NotImplementedError
