from abc import ABCMeta, abstractmethod

import numpy as np

import openmc.data


class AngleEnergy(object):
    """Distribution in angle and energy of a secondary particle."""

    __metaclass = ABCMeta

    @abstractmethod
    def to_hdf5(self, group):
        pass

    @staticmethod
    def from_hdf5(group):
        """Generate angle-energy distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.AngleEnergy
            Angle-energy distribution

        """
        dist_type = group.attrs['type'].decode()
        if dist_type == 'uncorrelated':
            return openmc.data.UncorrelatedAngleEnergy.from_hdf5(group)
        elif dist_type == 'correlated':
            return openmc.data.CorrelatedAngleEnergy.from_hdf5(group)
        elif dist_type == 'kalbach-mann':
            return openmc.data.KalbachMann.from_hdf5(group)
        elif dist_type == 'nbody':
            return openmc.data.NBodyPhaseSpace.from_hdf5(group)
