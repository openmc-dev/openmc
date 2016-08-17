import numpy as np

import openmc.checkvalue as cv
from .angle_energy import AngleEnergy
from .energy_distribution import EnergyDistribution
from .angle_distribution import AngleDistribution


class UncorrelatedAngleEnergy(AngleEnergy):
    """Uncorrelated angle-energy distribution

    Parameters
    ----------
    angle : openmc.data.AngleDistribution
        Distribution of outgoing angles represented as scattering cosines
    energy : openmc.data.EnergyDistribution
        Distribution of outgoing energies

    Attributes
    ----------
    angle : openmc.data.AngleDistribution
        Distribution of outgoing angles represented as scattering cosines
    energy : openmc.data.EnergyDistribution
        Distribution of outgoing energies

    """

    def __init__(self, angle=None, energy=None):
        self._angle = None
        self._energy = None

        if angle is not None:
            self.angle = angle
        if energy is not None:
            self.energy = energy

    @property
    def angle(self):
        return self._angle

    @property
    def energy(self):
        return self._energy

    @angle.setter
    def angle(self, angle):
        cv.check_type('uncorrelated angle distribution', angle,
                      AngleDistribution)
        self._angle = angle

    @energy.setter
    def energy(self, energy):
        cv.check_type('uncorrelated energy distribution', energy,
                      EnergyDistribution)
        self._energy = energy

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['type'] = np.string_('uncorrelated')
        if self.angle is not None:
            angle_group = group.create_group('angle')
            self.angle.to_hdf5(angle_group)

        if self.energy is not None:
            energy_group = group.create_group('energy')
            self.energy.to_hdf5(energy_group)

    @classmethod
    def from_hdf5(cls, group):
        """Generate uncorrelated angle-energy distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.UncorrelatedAngleEnergy
            Uncorrelated angle-energy distribution

        """
        dist = cls()
        if 'angle' in group:
            dist.angle = AngleDistribution.from_hdf5(group['angle'])
        if 'energy' in group:
            dist.energy = EnergyDistribution.from_hdf5(group['energy'])
        return dist
