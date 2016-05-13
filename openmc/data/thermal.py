from collections import Iterable
from numbers import Real

import numpy as np

import openmc.checkvalue as cv


class CoherentElastic(object):
    """Coherent elastic scattering data from a crystalline material

    Parameters
    ----------
    bragg_edges : Iterable of float
        Bragg edge energies in MeV
    factors : Iterable of float
        Partial sum of structure factors, :math:`\sum\limits_{i=1}^{E_i<E} S_i`

    Attributes
    ----------
    bragg_edges : Iterable of float
        Bragg edge energies in MeV
    factors : Iterable of float
        Partial sum of structure factors, :math:`\sum\limits_{i=1}^{E_i<E} S_i`

    """

    def __init__(self, bragg_edges, factors):
        self.bragg_edges = bragg_edges
        self.factors = factors

    def __call__(self, E):
        if isinstance(E, Iterable):
            E = np.asarray(E)
        idx = np.searchsorted(self.bragg_edges, E)
        return self.factors[idx]/E

    def __len__(self):
        return len(self.bragg_edges)

    @property
    def bragg_edges(self):
        return self._bragg_edges

    @property
    def factors(self):
        return self._factors

    @bragg_edges.setter
    def bragg_edges(self, bragg_edges):
        cv.check_type('Bragg edges', bragg_edges, Iterable, Real)
        self._bragg_edges = np.asarray(bragg_edges)

    @factors.setter
    def factors(self, factors):
        cv.check_type('structure factor cumulative sums', factors,
                      Iterable, Real)
        self._factors = np.asarray(factors)

    def to_hdf5(self, group, name):
        """Write coherent elastic scattering to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to
        name : str
            Name of the dataset to create

        """
        dataset = group.create_dataset(name, data=np.vstack(
            [self.bragg_edges, self.factors]))
        dataset.attrs['type'] = np.string_('bragg')

    @classmethod
    def from_hdf5(cls, dataset):
        """Read coherent elastic scattering from an HDF5 dataset

        Parameters
        ----------
        group : h5py.Dataset
            HDF5 group to write to

        Returns
        -------
        openmc.data.CoherentElastic
            Coherent elastic scattering cross section

        """
        bragg_edges = dataset.value[0,:]
        factors = dataset.value[1,:]
        return cls(bragg_edges, factors)
