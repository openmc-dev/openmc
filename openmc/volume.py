from collections import Iterable, Mapping
from numbers import Real, Integral
from xml.etree import ElementTree as ET
from warnings import warn

import numpy as np
import pandas as pd

import openmc
import openmc.checkvalue as cv


class VolumeCalculation(object):
    """Stochastic volume calculation specifications and results.

    Parameters
    ----------
    domains : Iterable of openmc.Cell, openmc.Material, or openmc.Universe
        Domains to find volumes of
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
    ids : Iterable of int
        IDs of domains to find volumes of
    domain_type : {'cell', 'material', 'universe'}
        Type of each domain
    samples : int
        Number of samples used to generate volume estimates
    lower_left : Iterable of float
        Lower-left coordinates of bounding box used to sample points
    upper_right : Iterable of float
        Upper-right coordinates of bounding box used to sample points
    results : dict
        Dictionary whose keys are unique IDs of domains and values are
        dictionaries with calculated volumes and total number of atoms for each
        nuclide present in the domain.
    volumes : dict
        Dictionary whose keys are unique IDs of domains and values are the
        estimated volumes
    atoms_dataframe : pandas.DataFrame
        DataFrame showing the estimated number of atoms for each nuclide present
        in each domain specified.

    """
    def __init__(self, domains, samples, lower_left=None,
                 upper_right=None):
        self._results = None

        cv.check_type('domains', domains, Iterable,
                      (openmc.Cell, openmc.Material, openmc.Universe))
        if isinstance(domains[0], openmc.Cell):
            self._domain_type = 'cell'
        elif isinstance(domains[0], openmc.Material):
            self._domain_type = 'material'
        elif isinstance(domains[0], openmc.Universe):
            self._domain_type = 'universe'
        self.ids = [d.id for d in domains]

        self.samples = samples

        if lower_left is not None:
            if upper_right is None:
                raise ValueError('Both lower-left and upper-right coordinates '
                                 'should be specified')

            # For cell domains, try to compute bounding box and make sure
            # user-specified one is valid
            if self.domain_type == 'cell':
                for c in domains:
                    if c.region is None:
                        continue
                    ll, ur = c.region.bounding_box
                    if np.any(np.isinf(ll)) or np.any(np.isinf(ur)):
                        continue
                    if (np.any(np.asarray(lower_left) > ll) or
                        np.any(np.asarray(upper_right) < ur)):
                        warn("Specified bounding box is smaller than computed "
                             "bounding box for cell {}. Volume calculation may "
                             "be incorrect!".format(c.id))

            self.lower_left = lower_left
            self.upper_right = upper_right
        else:
            if self.domain_type == 'cell':
                ll, ur = openmc.Union(*[c.region for c in domains]).bounding_box
                if np.any(np.isinf(ll)) or np.any(np.isinf(ur)):
                    raise ValueError('Could not automatically determine bounding '
                                     'box for stochastic volume calculation.')
                else:
                    self.lower_left = ll
                    self.upper_right = ur
            else:
                raise ValueError('Could not automatically determine bounding box '
                                 'for stochastic volume calculation.')

    @property
    def ids(self):
        return self._ids

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
    def domain_type(self):
        return self._domain_type

    @property
    def volumes(self):
        return {uid: results['volume'] for uid, results in self.results.items()}

    @property
    def atoms_dataframe(self):
        items = []
        columns = [self.domain_type.capitalize(), 'Nuclide', 'Atoms',
                   'Uncertainty']
        for uid, results in self.results.items():
            for name, atoms in results['atoms']:
                items.append((uid, name, atoms[0], atoms[1]))

        return pd.DataFrame.from_records(items, columns=columns)

    @ids.setter
    def ids(self, ids):
        cv.check_type('domain IDs', ids, Iterable, Real)
        self._ids = ids

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
            domain_type = f.attrs['domain_type'].decode()
            samples = f.attrs['samples']
            lower_left = f.attrs['lower_left']
            upper_right = f.attrs['upper_right']

            results = {}
            ids = []
            for obj_name in f:
                if obj_name.startswith('domain_'):
                    domain_id = int(obj_name[7:])
                    ids.append(domain_id)
                    group = f[obj_name]
                    volume = tuple(group['volume'].value)
                    nucnames = group['nuclides'].value
                    atoms = group['atoms'].value

                    atom_list = []
                    for name_i, atoms_i in zip(nucnames, atoms):
                        atom_list.append((name_i.decode(), tuple(atoms_i)))
                    results[domain_id] = {'volume': volume, 'atoms': atom_list}

        # Instantiate some throw-away domains that are used by the constructor
        # to assign IDs
        if domain_type == 'cell':
            domains = [openmc.Cell(uid) for uid in ids]
        elif domain_type == 'material':
            domains = [openmc.Material(uid) for uid in ids]
        elif domain_type == 'universe':
            domains = [openmc.Universe(uid) for uid in ids]

        # Instantiate the class and assign results
        vol = cls(domains, samples, lower_left, upper_right)
        vol.results = results
        return vol

    def to_xml_element(self):
        """Return XML representation of the volume calculation

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing volume calculation data

        """
        element = ET.Element("volume_calc")
        dt_elem = ET.SubElement(element, "domain_type")
        dt_elem.text = self.domain_type
        id_elem = ET.SubElement(element, "domain_ids")
        id_elem.text = ' '.join(str(uid) for uid in self.ids)
        samples_elem = ET.SubElement(element, "samples")
        samples_elem.text = str(self.samples)
        ll_elem = ET.SubElement(element, "lower_left")
        ll_elem.text = ' '.join(str(x) for x in self.lower_left)
        ur_elem = ET.SubElement(element, "upper_right")
        ur_elem.text = ' '.join(str(x) for x in self.upper_right)
        return element
