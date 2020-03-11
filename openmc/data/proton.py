from collections import OrderedDict
from collections.abc import Mapping, Callable
from io import StringIO
from numbers import Integral, Real

import h5py
import numpy as np

from openmc.mixin import EqualityMixin
import openmc.checkvalue as cv
from . import HDF5_VERSION, HDF5_VERSION_MAJOR
from .ace import Table, get_metadata, get_table
from .data import ATOMIC_SYMBOL
from .endf import Evaluation, get_head_record, get_tab1_record, get_list_record
from .function import Tabulated1D

_REACTION_NAME = {2: '(p,elastic)', 4: '(p,level)', 5: '(p,misc)', 11: '(p,2nd)', 
                  16: '(p,2n)', 17: '(p,3n)', 18: '(p,fission)', 22: '(p,na)', 
                  24: '(p,2na)', 25: '(p,3na)', 28: '(p,np)', 29: '(p,n2a)',
                  30: '(p,2n2a)', 32: '(p,nd)', 33: '(p,nt)', 34: '(p,nHe-3)',
                  35: '(p,nd2a)', 36: '(p,nt2a)', 37: '(p,4n)', 41: '(p,2np)', 
                  42: '(p,3np)', 44: '(p,n2p)', 45: '(p,npa)', 91: '(p,nc)', 
                  102: '(p,gamma)', 103: '(p,p)', 104: '(p,d)', 105: '(p,t)', 
                  106: '(p,3He)', 107: '(p,a)', 108: '(p,2a)', 109: '(p,3a)', 
                  111: '(p,2p)', 112: '(p,pa)', 113: '(p,t2a)', 114: '(p,d2a)', 
                  115: '(p,pd)', 116: '(p,pt)', 117: '(p,da)',  203: '(p,Xp)',
                 204: '(p,Xd)', 205: '(p,Xt)', 206: '(p,X3He)', 207: '(p,Xa)',
                 301: 'heating', 444: 'damage-energy', 699: '(p,dc)', 749: '(p,tc)', 
                 799: '(p,3Hec)', 849: '(p,ac)', 891: '(p,2nc)', 
                 901: 'heating-local'}
_REACTION_NAME.update({i: '(p,n{})'.format(i - 50) for i in range(50, 91)})
_REACTION_NAME.update({i: '(p,d{})'.format(i - 650) for i in range(650, 699)})
_REACTION_NAME.update({i: '(p,t{})'.format(i - 700) for i in range(700, 749)})
_REACTION_NAME.update({i: '(p,3He{})'.format(i - 750) for i in range(750, 799)})
_REACTION_NAME.update({i: '(p,a{})'.format(i - 800) for i in range(800, 849)})
_REACTION_NAME.update({i: '(p,2n{})'.format(i - 875) for i in range(875, 891)})

class ProtonReaction(EqualityMixin):

    def __init__(self, mt):
        self.mt = mt
        self._xs = None

    def __repr__(self):
        if self.mt in _REACTION_NAME:
            return "<Photon Reaction: MT={} {}>".format(
                self.mt, _REACTION_NAME[self.mt][0])
        else:
            return "<Photon Reaction: MT={}>".format(self.mt)

    @property
    def xs(self):
        return self._xs

    @xs.setter
    def xs(self, xs):
        cv.check_type('reaction cross section', xs, Callable)
        self._xs = xs

    @classmethod
    def from_endf(cls, ev, mt):
        """Generate proton reaction from an ENDF evaluation

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF photo-atomic interaction data evaluation
        mt : int
            The MT value of the reaction to get data for

        Returns
        -------
        openmc.data.PhotonReaction
            Photon reaction data

        """

        rx = cls(mt)

        # Integrated cross sections
        if (3, mt) in ev.section:
            file_obj = StringIO(ev.section[3, mt])
            get_head_record(file_obj)
            rx.xs = get_tab1_record(file_obj)[1]

        
        return rx

    #  A little redundant
    def to_hdf5(self, group, energy):
        """Write proton reaction to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['mt'] = self.mt
        if self.mt in _REACTION_NAME:
            group.attrs['label'] = np.string_(_REACTION_NAME[self.mt])
        else:
            group.attrs['label'] = np.string_(self.mt)

        group.create_dataset('xs', data=self.xs(energy))

    @classmethod
    def from_hdf5(cls, group, energy):
        """Generate proton reaction from an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from
        energy : dict
            Array of energies at which cross sections are tabulated at.

        Returns
        -------
        openmc.data.PhotonReaction
            Reaction data

        """

        mt = group.attrs['mt']
        rx = cls(mt)

        # Read cross section
        xs = group['xs'][()]
        tabulated_xs = Tabulated1D(energy, xs)
        rx.xs = tabulated_xs

        return rx

class IncidentProton(EqualityMixin):
    """Photon interaction data.

    This class stores data derived from an ENDF-6 format proton interaction
    sublibrary. Instances of this class are not normally instantiated by the
    user but rather created using the factory methods
    :meth:`IncidentProton.from_hdf5`, :meth:`IncidentProton.from_ace`, and
    :meth:`IncidentProton.from_endf`.

    Parameters
    ----------
    atomic_number : int
        Number of protons in the target nucleus

    Attributes
    ----------
    atomic_number : int
        Number of protons in the target nucleus
    reactions : collections.OrderedDict
        Contains the cross sections for each photon reaction. The keys are MT
        values and the values are instances of :class:`ProtonReaction`.

    """

    def __init__(self, name, atomic_number, mass_number):
        self.name = name
        self.atomic_number = atomic_number
        self.mass_number = mass_number
        self.reactions = OrderedDict()
        self.energy = []

    def __contains__(self, mt):
        return mt in self.reactions

    def __getitem__(self, mt):
        if mt in self.reactions:
            return self.reactions[mt]
        else:
            raise KeyError('No reaction with MT={}.'.format(mt))

    def __repr__(self):
        return "<IncidentProton: {}>".format(self.name)

    def __iter__(self):
        return iter(self.reactions.values())

    @property
    def atomic_number(self):
        return self._atomic_number

    @property
    def name(self):
        return self._name
    
    @property
    def mass_number(self):
        return self._mass_number


    @name.setter
    def name(self, name):
        cv.check_type('name', name, str)
        self._name = name

    @atomic_number.setter
    def atomic_number(self, atomic_number):
        cv.check_type('atomic number', atomic_number, Integral)
        cv.check_greater_than('atomic number', atomic_number, 0, True)
        self._atomic_number = atomic_number

    @mass_number.setter
    def mass_number(self, mass_number):
        cv.check_type('mass number', mass_number, Integral)
        cv.check_greater_than('mass number', mass_number, 0, True)
        self._mass_number = mass_number


    @classmethod
    def from_endf(cls, ev_or_filename):
        """Generate incident proton data from an ENDF evaluation

        Parameters
        ----------
        ev_or_filename : openmc.data.endf.Evaluation or str
            ENDF evaluation to read from. If given as a string, it is assumed to
            be the filename for the ENDF file.
  

        Returns
        -------
        openmc.data.IncidentProton
            Proton interaction data

      """

        if isinstance(ev_or_filename, Evaluation):
            ev = ev_or_filename
        else:
            ev = Evaluation(ev_or_filename)


        atomic_number = ev.target['atomic_number']
        mass_number = ev.target['mass_number']

        element = ATOMIC_SYMBOL[atomic_number]
        name = '{}{}'.format(element, mass_number)

        data = cls(name, atomic_number, mass_number)

        # Read each reaction
        for mf, mt, nc, mod in ev.reaction_list:
            if mf == 3:
                data.reactions[mt] = ProtonReaction.from_endf(ev, mt)

        # Determine union energy grid
        union_grid = np.array([])
        for rx in data:
            union_grid = np.union1d(union_grid, rx.xs.x)

        data.energy = union_grid
        
        return data

    def export_to_hdf5(self, path, mode='a', libver='earliest'):
        """Export incident photon data to an HDF5 file.

        Parameters
        ----------
        path : str
            Path to write HDF5 file to
        mode : {'r', 'r+', 'w', 'x', 'a'}
            Mode that is used to open the HDF5 file. This is the second argument
            to the :class:`h5py.File` constructor.
        libver : {'earliest', 'latest'}
            Compatibility mode for the HDF5 file. 'latest' will produce files
            that are less backwards compatible but have performance benefits.

        """
        # Open file and write version
        f = h5py.File(str(path), mode, libver=libver)
        f.attrs['filetype'] = np.string_('data_proton')
        if 'version' not in f.attrs:
            f.attrs['version'] = np.array(HDF5_VERSION)

        group = f.create_group(self.name)
        group.attrs['Z'] = self.atomic_number
        group.attrs['A'] = self.mass_number

        # Determine union energy grid
        union_grid = np.array([])
        for rx in self:
            union_grid = np.union1d(union_grid, rx.xs.x)
        group.create_dataset('energy', data=union_grid)

        # Write cross sections
        rxs_group = group.create_group('reactions')
        for mt, rx in self.reactions.items():
            key = _REACTION_NAME[mt]
            # Only elastic cross section for now
            if mt == 2:
                rx_group = rxs_group.create_group('reaction_{:03}'.format(rx.mt))
            else:
                continue
             
            rx.to_hdf5(rx_group, union_grid)
        
        f.close()

    @classmethod
    def from_hdf5(cls, group_or_filename):
        """Generate proton interaction data from HDF5 group

        Parameters
        ----------
        group_or_filename : h5py.Group or str
            HDF5 group containing interaction data. If given as a string, it is
            assumed to be the filename for the HDF5 file, and the first group is
            used to read from.

        Returns
        -------
        openmc.data.IncidentProton
            Proton interaction data

        """
        if isinstance(group_or_filename, h5py.Group):
            group = group_or_filename
        else:
            h5file = h5py.File(str(group_or_filename), 'r')

            # Make sure version matches
            if 'version' in h5file.attrs:
                major, minor = h5file.attrs['version']
                # For now all versions of HDF5 data can be read
            else:
                raise IOError(
                    'HDF5 data does not indicate a version. Your installation of '
                    'the OpenMC Python API expects version {}.x data.'
                    .format(HDF5_VERSION_MAJOR))

            group = list(h5file.values())[0]

        name = group.name[1:]
        atomic_number = group.attrs['Z']
        mass_number = group.attrs['A']

        data = cls(name, atomic_number, mass_number)

        # Read energy grid
        data.energy = group['energy'][()]

        # Read reaction data
        rxs_group = group['reactions']
        for name, obj in sorted(rxs_group.items()):
            if name.startswith('reaction_'):
                rx = ProtonReaction.from_hdf5(obj, data.energy)
                data.reactions[rx.mt] = rx

        return data



