from abc import ABCMeta, abstractmethod
from io import StringIO

import openmc.data
from openmc.mixin import EqualityMixin


class AngleEnergy(EqualityMixin, metaclass=ABCMeta):
    """Distribution in angle and energy of a secondary particle."""
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
        elif dist_type == 'coherent_elastic':
            return openmc.data.CoherentElasticAE.from_hdf5(group)
        elif dist_type == 'incoherent_elastic':
            return openmc.data.IncoherentElasticAE.from_hdf5(group)
        elif dist_type == 'incoherent_elastic_discrete':
            return openmc.data.IncoherentElasticAEDiscrete.from_hdf5(group)
        elif dist_type == 'incoherent_inelastic_discrete':
            return openmc.data.IncoherentInelasticAEDiscrete.from_hdf5(group)
        elif dist_type == 'incoherent_inelastic':
            return openmc.data.IncoherentInelasticAE.from_hdf5(group)

    @staticmethod
    def from_ace(ace, location_dist, location_start, rx=None):
        """Generate an angle-energy distribution from ACE data

        Parameters
        ----------
        ace : openmc.data.ace.Table
            ACE table to read from
        location_dist : int
            Index in the XSS array corresponding to the start of a block,
            e.g. JXS(11) for the the DLW block.
        location_start : int
            Index in the XSS array corresponding to the start of an energy
            distribution array
        rx : Reaction
            Reaction this energy distribution will be associated with

        Returns
        -------
        distribution : openmc.data.AngleEnergy
            Secondary angle-energy distribution

        """
        # Set starting index for energy distribution
        idx = location_dist + location_start - 1

        law = int(ace.xss[idx + 1])
        location_data = int(ace.xss[idx + 2])

        # Position index for reading law data
        idx = location_dist + location_data - 1

        # Parse energy distribution data
        if law == 2:
            distribution = openmc.data.UncorrelatedAngleEnergy()
            distribution.energy = openmc.data.DiscretePhoton.from_ace(ace, idx)
        elif law in (3, 33):
            distribution = openmc.data.UncorrelatedAngleEnergy()
            distribution.energy = openmc.data.LevelInelastic.from_ace(ace, idx)
        elif law == 4:
            distribution = openmc.data.UncorrelatedAngleEnergy()
            distribution.energy = openmc.data.ContinuousTabular.from_ace(
                ace, idx, location_dist)
        elif law == 5:
            distribution = openmc.data.UncorrelatedAngleEnergy()
            distribution.energy = openmc.data.GeneralEvaporation.from_ace(ace, idx)
        elif law == 7:
            distribution = openmc.data.UncorrelatedAngleEnergy()
            distribution.energy = openmc.data.MaxwellEnergy.from_ace(ace, idx)
        elif law == 9:
            distribution = openmc.data.UncorrelatedAngleEnergy()
            distribution.energy = openmc.data.Evaporation.from_ace(ace, idx)
        elif law == 11:
            distribution = openmc.data.UncorrelatedAngleEnergy()
            distribution.energy = openmc.data.WattEnergy.from_ace(ace, idx)
        elif law == 44:
            distribution = openmc.data.KalbachMann.from_ace(
                ace, idx, location_dist)
        elif law == 61:
            distribution = openmc.data.CorrelatedAngleEnergy.from_ace(
                ace, idx, location_dist)
        elif law == 66:
            distribution = openmc.data.NBodyPhaseSpace.from_ace(
                ace, idx, rx.q_value)
        else:
            raise ValueError("Unsupported ACE secondary energy "
                             "distribution law {}".format(law))

        return distribution
