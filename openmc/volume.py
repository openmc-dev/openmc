from collections import Iterable, Mapping
from numbers import Real, Integral
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

from openmc import Cell, Union
import openmc.checkvalue as cv


class VolumeCalculation(object):
    """Stochastic volume calculation specifications and results.

    Parameters
    ----------
    cells : Iterable of Cell
        Cells to find volumes of
    samples : int
        Number of samples used to generate volume estimates
    lower_left : Iterable of float
        Lower-left coordinates of bounding box used to sample points. If this
        argument is not supplied, an attempt is made to automatically determine
        a bounding box.
    upper_right : Iterable of float
        Upper-right coordinates of bounding box used to sample points. If this
        argument is not supplied, an attempt is made to automatically determine
        a bounding box.

    Attributes
    ----------
    cell_ids : Iterable of int
        IDs of cells to find volumes of
    samples : int
        Number of samples used to generate volume estimates
    lower_left : Iterable of float
        Lower-left coordinates of bounding box used to sample points
    upper_right : Iterable of float
        Upper-right coordinates of bounding box used to sample points
    results : dict
        Dictionary whose keys are unique IDs of cells and values are
        dictionaries with calculated volumes and total number of atoms for each
        nuclide present in the cell.
    volumes : dict
        Dictionary whose keys are unique IDs of cells and values are the
        estimated volumes
    atoms_dataframe : pandas.DataFrame
        DataFrame showing the estimated number of atoms for each nuclide present
        in each cell specified.

    """
    def __init__(self, cells, samples, lower_left=None,
                 upper_right=None):
        self._results = None

        cv.check_type('cells', cells, Iterable, Cell)
        self.cell_ids = [c.id for c in cells]
        self.samples = samples

        if lower_left is not None:
            self.lower_left = lower_left
            if upper_right is None:
                raise ValueError('Both lower-left and upper-right coordinates '
                                 'should be specified')
            self.upper_right = upper_right
        else:
            ll, ur = Union(*[c.region for c in cells]).bounding_box
            if np.any(np.isinf(ll)) or np.any(np.isinf(ur)):
                raise ValueError('Could not automatically determine bounding box '
                                 'for stochastic volume calculation.')
            else:
                self.lower_left = ll
                self.upper_right = ur

    @property
    def cell_ids(self):
        return self._cell_ids

    @property
    def samples(self):
        return self._samples

    @property
    def lower_left(self):
        return self._lower_left

    @property
    def upper_right(self):
        return self._upper_right

    @property
    def results(self):
        return self._results

    @property
    def volumes(self):
        return {uid: results['volume'] for uid, results in self.results.items()}

    @property
    def atoms_dataframe(self):
        items = []
        columns = ['Cell', 'Nuclide', 'Atoms', 'Uncertainty']
        for cell_id, results in self.results.items():
            for name, atoms in results['atoms']:
                items.append((cell_id, name, atoms[0], atoms[1]))

        return pd.DataFrame.from_records(items, columns=columns)

    @cell_ids.setter
    def cell_ids(self, cell_ids):
        cv.check_type('cell IDs', cell_ids, Iterable, Real)
        self._cell_ids = cell_ids

    @samples.setter
    def samples(self, samples):
        cv.check_type('number of samples', samples, Integral)
        cv.check_greater_than('number of samples', samples, 0)
        self._samples = samples

    @lower_left.setter
    def lower_left(self, lower_left):
        name = 'lower-left bounding box coordinates',
        cv.check_type(name, lower_left, Iterable, Real)
        cv.check_length(name, lower_left, 3)
        self._lower_left = lower_left

    @upper_right.setter
    def upper_right(self, upper_right):
        name = 'upper-right bounding box coordinates'
        cv.check_type(name, upper_right, Iterable, Real)
        cv.check_length(name, upper_right, 3)
        self._upper_right = upper_right

    @results.setter
    def results(self, results):
        cv.check_type('results', results, Mapping)
        self._results = results

    @classmethod
    def from_hdf5(cls, filename):
        """Load stochastic volume calculation results from HDF5 file.

        Parameters
        ----------
        filename : str
            Path to volume.h5 file

        Returns
        -------
        openmc.VolumeCalculation
            Results of the stochastic volume calculation

        """
        import h5py

        with h5py.File(filename, 'r') as f:
            samples = f.attrs['samples']
            lower_left = f.attrs['lower_left']
            upper_right = f.attrs['upper_right']

            results = {}
            cell_ids = []
            for obj_name in f:
                if obj_name.startswith('cell_'):
                    cell_id = int(obj_name[5:])
                    cell_ids.append(cell_id)
                    group = f[obj_name]
                    volume = tuple(group['volume'].value)
                    nucnames = group['nuclides'].value
                    atoms = group['atoms'].value

                    atom_list = []
                    for name_i, atoms_i in zip(nucnames, atoms):
                        atom_list.append((name_i.decode(), tuple(atoms_i)))
                    results[cell_id] = {'volume': volume, 'atoms': atom_list}

        # Instantiate some throw-away cells that are used by the constructor to
        # assign IDs
        cells = [Cell(uid) for uid in cell_ids]

        # Instantiate the class and assign results
        vol = cls(cells, samples, lower_left, upper_right)
        vol.results = results
        return vol

    def to_xml(self):
        """Return XML representation of the volume calculation

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing volume calculation data

        """
        element = ET.Element("volume_calc")
        cell_elem = ET.SubElement(element, "cells")
        cell_elem.text = ' '.join(str(uid) for uid in self.cell_ids)
        samples_elem = ET.SubElement(element, "samples")
        samples_elem.text = str(self.samples)
        ll_elem = ET.SubElement(element, "lower_left")
        ll_elem.text = ' '.join(str(x) for x in self.lower_left)
        ur_elem = ET.SubElement(element, "upper_right")
        ur_elem.text = ' '.join(str(x) for x in self.upper_right)
        return element
