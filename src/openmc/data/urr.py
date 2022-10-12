from collections.abc import Iterable
from numbers import Integral, Real

import numpy as np

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from .data import EV_PER_MEV


class ProbabilityTables(EqualityMixin):
    r"""Unresolved resonance region probability tables.

    Parameters
    ----------
    energy : Iterable of float
        Energies in eV at which probability tables exist
    table : numpy.ndarray
        Probability tables for each energy. This array is of shape (N, 6, M)
        where N is the number of energies and M is the number of bands. The
        second dimension indicates whether the value is for the cumulative
        probability (0), total (1), elastic (2), fission (3), :math:`(n,\gamma)`
        (4), or heating number (5).
    interpolation : {2, 5}
        Interpolation scheme between tables
    inelastic_flag : int
        A value less than zero indicates that the inelastic cross section is
        zero within the unresolved energy range. A value greater than zero
        indicates the MT number for a reaction whose cross section is to be used
        in the unresolved range.
    absorption_flag : int
        A value less than zero indicates that the "other absorption" cross
        section is zero within the unresolved energy range. A value greater than
        zero indicates the MT number for a reaction whose cross section is to be
        used in the unresolved range.
    multiply_smooth : bool
        Indicate whether probability table values are cross sections (False) or
        whether they must be multiply by the corresponding "smooth" cross
        sections (True).

    Attributes
    ----------
    energy : Iterable of float
        Energies in eV at which probability tables exist
    table : numpy.ndarray
        Probability tables for each energy. This array is of shape (N, 6, M)
        where N is the number of energies and M is the number of bands. The
        second dimension indicates whether the value is for the cumulative
        probability (0), total (1), elastic (2), fission (3), :math:`(n,\gamma)`
        (4), or heating number (5).
    interpolation : {2, 5}
        Interpolation scheme between tables
    inelastic_flag : int
        A value less than zero indicates that the inelastic cross section is
        zero within the unresolved energy range. A value greater than zero
        indicates the MT number for a reaction whose cross section is to be used
        in the unresolved range.
    absorption_flag : int
        A value less than zero indicates that the "other absorption" cross
        section is zero within the unresolved energy range. A value greater than
        zero indicates the MT number for a reaction whose cross section is to be
        used in the unresolved range.
    multiply_smooth : bool
        Indicate whether probability table values are cross sections (False) or
        whether they must be multiply by the corresponding "smooth" cross
        sections (True).
    """

    def __init__(self, energy, table, interpolation, inelastic_flag=-1,
                 absorption_flag=-1, multiply_smooth=False):
        self.energy = energy
        self.table = table
        self.interpolation = interpolation
        self.inelastic_flag = inelastic_flag
        self.absorption_flag = absorption_flag
        self.multiply_smooth = multiply_smooth

    @property
    def absorption_flag(self):
        return self._absorption_flag

    @property
    def energy(self):
        return self._energy

    @property
    def inelastic_flag(self):
        return self._inelastic_flag

    @property
    def interpolation(self):
        return self._interpolation

    @property
    def multiply_smooth(self):
        return self._multiply_smooth

    @property
    def table(self):
        return self._table

    @absorption_flag.setter
    def absorption_flag(self, absorption_flag):
        cv.check_type('absorption flag', absorption_flag, Integral)
        self._absorption_flag = absorption_flag

    @energy.setter
    def energy(self, energy):
        cv.check_type('probability table energies', energy, Iterable, Real)
        self._energy = energy

    @inelastic_flag.setter
    def inelastic_flag(self, inelastic_flag):
        cv.check_type('inelastic flag', inelastic_flag, Integral)
        self._inelastic_flag = inelastic_flag

    @interpolation.setter
    def interpolation(self, interpolation):
        cv.check_value('interpolation', interpolation, [2, 5])
        self._interpolation = interpolation

    @multiply_smooth.setter
    def multiply_smooth(self, multiply_smooth):
        cv.check_type('multiply by smooth', multiply_smooth, bool)
        self._multiply_smooth = multiply_smooth

    @table.setter
    def table(self, table):
        cv.check_type('probability tables', table, np.ndarray)
        self._table = table

    def to_hdf5(self, group):
        """Write probability tables to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['interpolation'] = self.interpolation
        group.attrs['inelastic'] = self.inelastic_flag
        group.attrs['absorption'] = self.absorption_flag
        group.attrs['multiply_smooth'] = int(self.multiply_smooth)

        group.create_dataset('energy', data=self.energy)
        group.create_dataset('table', data=self.table)

    @classmethod
    def from_hdf5(cls, group):
        """Generate probability tables from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.ProbabilityTables
            Probability tables

        """
        interpolation = group.attrs['interpolation']
        inelastic_flag = group.attrs['inelastic']
        absorption_flag = group.attrs['absorption']
        multiply_smooth = bool(group.attrs['multiply_smooth'])

        energy = group['energy'][()]
        table = group['table'][()]

        return cls(energy, table, interpolation, inelastic_flag,
                   absorption_flag, multiply_smooth)

    @classmethod
    def from_ace(cls, ace):
        """Generate probability tables from an ACE table

        Parameters
        ----------
        ace : openmc.data.ace.Table
            ACE table to read from

        Returns
        -------
        openmc.data.ProbabilityTables
            Unresolved resonance region probability tables

        """
        # Check if URR probability tables are present
        idx = ace.jxs[23]
        if idx == 0:
            return None

        N = int(ace.xss[idx])      # Number of incident energies
        M = int(ace.xss[idx+1])    # Length of probability table
        interpolation = int(ace.xss[idx+2])
        inelastic_flag = int(ace.xss[idx+3])
        absorption_flag = int(ace.xss[idx+4])
        multiply_smooth = (int(ace.xss[idx+5]) == 1)
        idx += 6

        # Get energies at which tables exist
        energy = ace.xss[idx : idx+N]*EV_PER_MEV
        idx += N

        # Get probability tables
        table = ace.xss[idx : idx+N*6*M].copy()
        table.shape = (N, 6, M)

        # Convert units on heating numbers
        table[:,5,:] *= EV_PER_MEV

        return cls(energy, table, interpolation, inelastic_flag,
                   absorption_flag, multiply_smooth)
