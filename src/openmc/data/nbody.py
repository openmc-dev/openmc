from numbers import Real, Integral

import numpy as np

import openmc.checkvalue as cv
from .angle_energy import AngleEnergy
from .endf import get_cont_record

class NBodyPhaseSpace(AngleEnergy):
    """N-body phase space distribution

    Parameters
    ----------
    total_mass : float
        Total mass of product particles
    n_particles : int
        Number of product particles
    atomic_weight_ratio : float
        Atomic weight ratio of target nuclide
    q_value : float
        Q value for reaction in eV

    Attributes
    ----------
    total_mass : float
        Total mass of product particles
    n_particles : int
        Number of product particles
    atomic_weight_ratio : float
        Atomic weight ratio of target nuclide
    q_value : float
        Q value for reaction in eV

    """

    def __init__(self, total_mass, n_particles, atomic_weight_ratio, q_value):
        self.total_mass = total_mass
        self.n_particles = n_particles
        self.atomic_weight_ratio = atomic_weight_ratio
        self.q_value = q_value

    @property
    def total_mass(self):
        return self._total_mass

    @property
    def n_particles(self):
        return self._n_particles

    @property
    def atomic_weight_ratio(self):
        return self._atomic_weight_ratio

    @property
    def q_value(self):
        return self._q_value

    @total_mass.setter
    def total_mass(self, total_mass):
        name = 'N-body phase space total mass'
        cv.check_type(name, total_mass, Real)
        cv.check_greater_than(name, total_mass, 0.)
        self._total_mass = total_mass

    @n_particles.setter
    def n_particles(self, n_particles):
        name = 'N-body phase space number of particles'
        cv.check_type(name, n_particles, Integral)
        cv.check_greater_than(name, n_particles, 0)
        self._n_particles = n_particles

    @atomic_weight_ratio.setter
    def atomic_weight_ratio(self, atomic_weight_ratio):
        name = 'N-body phase space atomic weight ratio'
        cv.check_type(name, atomic_weight_ratio, Real)
        cv.check_greater_than(name, atomic_weight_ratio, 0.0)
        self._atomic_weight_ratio = atomic_weight_ratio

    @q_value.setter
    def q_value(self, q_value):
        name = 'N-body phase space Q value'
        cv.check_type(name, q_value, Real)
        self._q_value = q_value

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['type'] = np.string_('nbody')
        group.attrs['total_mass'] = self.total_mass
        group.attrs['n_particles'] = self.n_particles
        group.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio
        group.attrs['q_value'] = self.q_value

    @classmethod
    def from_hdf5(cls, group):
        """Generate N-body phase space distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.NBodyPhaseSpace
            N-body phase space distribution

        """
        total_mass = group.attrs['total_mass']
        n_particles = group.attrs['n_particles']
        awr = group.attrs['atomic_weight_ratio']
        q_value = group.attrs['q_value']
        return cls(total_mass, n_particles, awr, q_value)

    @classmethod
    def from_ace(cls, ace, idx, q_value):
        """Generate N-body phase space distribution from ACE data

        Parameters
        ----------
        ace : openmc.data.ace.Table
            ACE table to read from
        idx : int
            Index in XSS array of the start of the energy distribution data
            (LDIS + LOCC - 1)
        q_value : float
            Q-value for reaction in eV

        Returns
        -------
        openmc.data.NBodyPhaseSpace
            N-body phase space distribution

        """
        n_particles = int(ace.xss[idx])
        total_mass = ace.xss[idx + 1]
        return cls(total_mass, n_particles, ace.atomic_weight_ratio, q_value)

    @classmethod
    def from_endf(cls, file_obj):
        """Generate N-body phase space distribution from an ENDF evaluation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positions at the start of the N-body phase space
            distribution

        Returns
        -------
        openmc.data.NBodyPhaseSpace
            N-body phase space distribution

        """
        items = get_cont_record(file_obj)
        total_mass = items[0]
        n_particles = items[5]
        # TODO: get awr and Q value
        return cls(total_mass, n_particles, 1.0, 0.0)
