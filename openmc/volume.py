from collections import OrderedDict
from collections.abc import Iterable, Mapping
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import warnings

import numpy as np
import pandas as pd
import h5py
from uncertainties import ufloat

import openmc
import openmc.checkvalue as cv

_VERSION_VOLUME = 1


class VolumeCalculation:
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
    threshold : float
        Threshold for the maximum standard deviation of volume in the calculation
    atoms : dict
        Dictionary mapping unique IDs of domains to a mapping of nuclides to
        total number of atoms for each nuclide present in the domain. For
        example, {10: {'U235': 1.0e22, 'U238': 5.0e22, ...}}.
    atoms_dataframe : pandas.DataFrame
        DataFrame showing the estimated number of atoms for each nuclide present
        in each domain specified.
    volumes : dict
        Dictionary mapping unique IDs of domains to estimated volumes in cm^3.
    threshold : float
        Threshold for the maxmimum standard deviation of volumes.
    trigger_type : {'variance', 'std_dev', 'rel_err'}
        Value type used to halt volume calculation
    iterations : int
        Number of iterations over samples (for calculations with a trigger).

    """
    def __init__(self, domains, samples, lower_left=None, upper_right=None):
        self._atoms = {}
        self._volumes = {}
        self._threshold = None
        self._trigger_type = None
        self._iterations = None

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
                    ll, ur = c.bounding_box
                    if np.any(np.isinf(ll)) or np.any(np.isinf(ur)):
                        continue
                    if (np.any(np.asarray(lower_left) > ll) or
                        np.any(np.asarray(upper_right) < ur)):
                        warnings.warn(
                            "Specified bounding box is smaller than computed "
                            "bounding box for cell {}. Volume calculation may "
                            "be incorrect!".format(c.id))

            self.lower_left = lower_left
            self.upper_right = upper_right
        else:
            if self.domain_type == 'cell':
                ll, ur = openmc.Union(c.region for c in domains).bounding_box
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
    def threshold(self):
        return self._threshold

    @property
    def trigger_type(self):
        return self._trigger_type

    @property
    def iterations(self):
        return self._iterations

    @property
    def domain_type(self):
        return self._domain_type

    @property
    def atoms(self):
        return self._atoms

    @property
    def volumes(self):
        return self._volumes

    @property
    def atoms_dataframe(self):
        items = []
        columns = [self.domain_type.capitalize(), 'Nuclide', 'Atoms']
        for uid, atoms_dict in self.atoms.items():
            for name, atoms in atoms_dict.items():
                items.append((uid, name, atoms))

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

    @threshold.setter
    def threshold(self, threshold):
        name = 'volume std. dev. threshold'
        cv.check_type(name, threshold, Real)
        cv.check_greater_than(name, threshold, 0.0)
        self._threshold = threshold

    @trigger_type.setter
    def trigger_type(self, trigger_type):
        cv.check_value('tally trigger type', trigger_type,
                       ('variance', 'std_dev', 'rel_err'))
        self._trigger_type = trigger_type

    @iterations.setter
    def iterations(self, iterations):
        name = 'volume calculation iterations'
        cv.check_type(name, iterations, Integral)
        cv.check_greater_than(name, iterations, 0)
        self._iterations = iterations

    @volumes.setter
    def volumes(self, volumes):
        cv.check_type('volumes', volumes, Mapping)
        self._volumes = volumes

    @atoms.setter
    def atoms(self, atoms):
        cv.check_type('atoms', atoms, Mapping)
        self._atoms = atoms

    def set_trigger(self, threshold, trigger_type):
        """Set a trigger on the voulme calculation

        Parameters
        ----------
        threshold : float
            Threshold for the maxmimum standard deviation of volumes
        trigger_type : {'variance', 'std_dev', 'rel_err'}
            Value type used to halt volume calculation
        """
        self.trigger_type = trigger_type
        self.threshold = threshold

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
        with h5py.File(filename, 'r') as f:
            cv.check_filetype_version(f, "volume", _VERSION_VOLUME)

            domain_type = f.attrs['domain_type'].decode()
            samples = f.attrs['samples']
            lower_left = f.attrs['lower_left']
            upper_right = f.attrs['upper_right']

            threshold = f.attrs.get('threshold')
            trigger_type = f.attrs.get('trigger_type')
            iterations = f.attrs.get('iterations', 1)

            volumes = {}
            atoms = {}
            ids = []
            for obj_name in f:
                if obj_name.startswith('domain_'):
                    domain_id = int(obj_name[7:])
                    ids.append(domain_id)
                    group = f[obj_name]
                    volume = ufloat(*group['volume'][()])
                    volumes[domain_id] = volume
                    nucnames = group['nuclides'][()]
                    atoms_ = group['atoms'][()]
                    atom_dict = OrderedDict()
                    for name_i, atoms_i in zip(nucnames, atoms_):
                        atom_dict[name_i.decode()] = ufloat(*atoms_i)
                    atoms[domain_id] = atom_dict

        # Instantiate some throw-away domains that are used by the constructor
        # to assign IDs
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', openmc.IDWarning)
            if domain_type == 'cell':
                domains = [openmc.Cell(uid) for uid in ids]
            elif domain_type == 'material':
                domains = [openmc.Material(uid) for uid in ids]
            elif domain_type == 'universe':
                domains = [openmc.Universe(uid) for uid in ids]

        # Instantiate the class and assign results
        vol = cls(domains, samples, lower_left, upper_right)

        if trigger_type is not None:
            vol.set_trigger(threshold, trigger_type.decode())

        vol.iterations = iterations
        vol.volumes = volumes
        vol.atoms = atoms
        return vol

    def load_results(self, filename):
        """Load stochastic volume calculation results from an HDF5 file.

        Parameters
        ----------
        filename : str
            Path to volume.h5 file

        """
        results = type(self).from_hdf5(filename)

        # Make sure properties match
        assert self.ids == results.ids
        assert np.all(self.lower_left == results.lower_left)
        assert np.all(self.upper_right == results.upper_right)

        # Copy results
        self.volumes = results.volumes
        self.atoms = results.atoms

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
        if self.threshold:
            trigger_elem = ET.SubElement(element, "threshold")
            trigger_elem.set("type", self.trigger_type)
            trigger_elem.set("threshold", str(self.threshold))
        return element
