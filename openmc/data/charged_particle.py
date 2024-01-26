from collections import OrderedDict
from collections.abc import Callable
from io import StringIO
from numbers import Integral, Real
from warnings import warn
import os
import tempfile

import h5py
import numpy as np

from openmc.mixin import EqualityMixin
import openmc.checkvalue as cv
from . import HDF5_VERSION, HDF5_VERSION_MAJOR
from .ace import Table, get_metadata, get_table, Library
from .data import ATOMIC_SYMBOL, EV_PER_MEV
from .endf import Evaluation, get_head_record, get_tab1_record
from .function import Tabulated1D
from .njoy import make_ace

_REACTION_NAME = {
    1: "(p,total)",
    2: "(p,elastic)",
    3: "(p,non-elastic)",
    4: "(p,level)",
    5: "(p,misc)",
    11: "(p,2nd)",
    16: "(p,2n)",
    17: "(p,3n)",
    18: "(p,fission)",
    22: "(p,na)",
    24: "(p,2na)",
    25: "(p,3na)",
    28: "(p,np)",
    29: "(p,n2a)",
    30: "(p,2n2a)",
    32: "(p,nd)",
    33: "(p,nt)",
    34: "(p,nHe-3)",
    35: "(p,nd2a)",
    36: "(p,nt2a)",
    37: "(p,4n)",
    41: "(p,2np)",
    42: "(p,3np)",
    44: "(p,n2p)",
    45: "(p,npa)",
    91: "(p,nc)",
    101: "(p,disapr)",
    102: "(p,gamma)",
    103: "(p,p)",
    104: "(p,d)",
    105: "(p,t)",
    106: "(p,3He)",
    107: "(p,a)",
    108: "(p,2a)",
    109: "(p,3a)",
    111: "(p,2p)",
    112: "(p,pa)",
    113: "(p,t2a)",
    114: "(p,d2a)",
    115: "(p,pd)",
    116: "(p,pt)",
    117: "(p,da)",
    203: "(p,Xp)",
    204: "(p,Xd)",
    205: "(p,Xt)",
    206: "(p,X3He)",
    207: "(p,Xa)",
    301: "heating",
    444: "damage-energy",
    699: "(p,dc)",
    749: "(p,tc)",
    799: "(p,3Hec)",
    849: "(p,ac)",
    891: "(p,2nc)",
    901: "heating-local",
}
_REACTION_NAME.update({i: "(p,n{})".format(i - 50) for i in range(50, 91)})
_REACTION_NAME.update({i: "(p,d{})".format(i - 650) for i in range(650, 699)})
_REACTION_NAME.update({i: "(p,t{})".format(i - 700) for i in range(700, 749)})
_REACTION_NAME.update({i: "(p,3He{})".format(i - 750) for i in range(750, 799)})
_REACTION_NAME.update({i: "(p,a{})".format(i - 800) for i in range(800, 849)})
_REACTION_NAME.update({i: "(p,2n{})".format(i - 875) for i in range(875, 891)})


class ChargedParticleReaction(EqualityMixin):
    def __init__(self, mt):
        self.mt = mt
        self._xs = None
        self._redundant = False

    def __repr__(self):
        if self.mt in _REACTION_NAME:
            return "<Proton Reaction: MT={} {}>".format(
                self.mt, _REACTION_NAME[self.mt]
            )
        else:
            return "<Proton Reaction: MT={}>".format(self.mt)

    @property
    def xs(self):
        return self._xs

    @property
    def redundant(self):
        return self._redundant

    @xs.setter
    def xs(self, xs):
        cv.check_type("reaction cross section", xs, Callable)
        self._xs = xs

    @redundant.setter
    def redundant(self, redundant):
        cv.check_type("redundant", redundant, (bool, np.bool_))
        self._redundant = redundant

    @classmethod
    def from_endf(cls, ev, mt):
        """Generate charged particle reaction from an ENDF evaluation

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation
            ENDF photo-atomic interaction data evaluation
        mt : int
            The MT value of the reaction to get data for

        Returns
        -------
        openmc.data.ChargedParticleReaction
            ChargedParticle reaction data

        """

        rx = cls(mt)

        # Integrated cross sections
        if (3, mt) in ev.section:
            file_obj = StringIO(ev.section[3, mt])
            get_head_record(file_obj)
            rx.xs = get_tab1_record(file_obj)[1]

        return rx

    def to_hdf5(self, group, energy):
        """Write proton reaction to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs["mt"] = self.mt
        group.attrs["redundant"] = 1 if self.redundant else 0

        if self.mt in _REACTION_NAME:
            group.attrs["label"] = np.string_(_REACTION_NAME[self.mt])
        else:
            group.attrs["label"] = np.string_(self.mt)

        dset = group.create_dataset("xs", data=self.xs(energy))
        threshold_idx = getattr(self.xs, "_threshold_idx", 0)
        dset.attrs["threshold_idx"] = threshold_idx

    @classmethod
    def from_ace(cls, ace, i_reaction):
        # Get nuclide energy grid
        n_grid = ace.nxs[3]
        grid = ace.xss[ace.jxs[1] : ace.jxs[1] + n_grid] * EV_PER_MEV

        if i_reaction > 0:
            mt = int(ace.xss[ace.jxs[3] + i_reaction - 1])
            rx = cls(mt)

            # ==================================================================
            # CROSS SECTION

            # Get locator for cross-section data
            loc = int(ace.xss[ace.jxs[6] + i_reaction - 1])

            # Determine starting index on energy grid
            threshold_idx = int(ace.xss[ace.jxs[7] + loc - 1]) - 1

            # Determine number of energies in reaction
            n_energy = int(ace.xss[ace.jxs[7] + loc])
            energy = grid[threshold_idx : threshold_idx + n_energy]

            # Read reaction cross section
            xs = ace.xss[ace.jxs[7] + loc + 1 : ace.jxs[7] + loc + 1 + n_energy]

            # Fix negatives -- known issue for Y89 in JEFF 3.2
            if np.any(xs < 0.0):
                warn(
                    "Negative cross sections found for MT={} in {}. Setting "
                    "to zero.".format(rx.mt, ace.name)
                )
                xs[xs < 0.0] = 0.0

            tabulated_xs = Tabulated1D(energy, xs)
            tabulated_xs._threshold_idx = threshold_idx
            rx.xs = tabulated_xs

        else:
            # Elastic scattering
            mt = 2
            rx = cls(mt)

            # Get elastic cross section values
            elastic_xs = ace.xss[ace.jxs[1] + 3 * n_grid : ace.jxs[1] + 4 * n_grid]

            # Fix negatives -- known issue for Ti46,49,50 in JEFF 3.2
            if np.any(elastic_xs < 0.0):
                warn(
                    "Negative elastic scattering cross section found for {}. "
                    "Setting to zero.".format(ace.name)
                )
                elastic_xs[elastic_xs < 0.0] = 0.0

            tabulated_xs = Tabulated1D(grid, elastic_xs)
            tabulated_xs._threshold_idx = 0
            rx.xs = tabulated_xs

        return rx

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
        openmc.data.ChargedParticleReaction
            Reaction data

        """

        mt = group.attrs["mt"]
        rx = cls(mt)

        # Read cross section
        xs = group["xs"][()]
        tabulated_xs = Tabulated1D(energy, xs)
        rx.xs = tabulated_xs

        return rx


class IncidentChargedParticle(EqualityMixin):
    """Charged particle interaction data.

    This class stores data derived from an ENDF-6 format proton interaction
    sublibrary. Instances of this class are not normally instantiated by the
    user but rather created using the factory methods
    :meth:`IncidentChargedParticle.from_hdf5`, :meth:`IncidentChargedParticle.from_ace`, and
    :meth:`IncidentChargedParticle.from_endf`.

    Parameters
    ----------
    name : str
        Name of the nuclide using the GND naming convention
    atomic_number : int
        Number of protons in the target nucleus
    atomic_number : int
        Number of protons in the target nucleus
    mass_number : int
        Number of nucleons in the target nucleus
    metastable : int
        Metastable state of the target nucleus. A value of zero indicates ground
        state.
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide.

    Attributes
    ----------
    atomic_number : int
        Number of protons in the target nucleus
    atomic_symbol : str
        Atomic symbol of the nuclide, e.g., 'Zr'
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide.
    mass_number : int
        Number of nucleons in the target nucleus
    metastable : int
        Metastable state of the target nucleus. A value of zero indicates ground
        state.
    name : str
        Name of the nuclide using the GND naming convention
    reactions : collections.OrderedDict
        Contains the cross sections, secondary angle and energy distributions,
        and other associated data for each reaction. The keys are the MT values
        and the values are Reaction objects.

    """

    def __init__(
        self, name, atomic_number, mass_number, metastable, atomic_weight_ratio
    ):
        self.name = name
        self.atomic_number = atomic_number
        self.mass_number = mass_number
        self.reactions = OrderedDict()
        self.energy = []
        self.metastable = metastable
        self.atomic_weight_ratio = atomic_weight_ratio

    def __contains__(self, mt):
        return mt in self.reactions

    def __getitem__(self, mt):
        if mt in self.reactions:
            return self.reactions[mt]
        else:
            raise KeyError("No reaction with MT={}.".format(mt))

    def __repr__(self):
        return "<IncidentChargedParticle: {}>".format(self.name)

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

    @property
    def metastable(self):
        return self._metastable

    @property
    def atomic_weight_ratio(self):
        return self._atomic_weight_ratio

    @property
    def atomic_symbol(self):
        return ATOMIC_SYMBOL[self.atomic_number]

    @name.setter
    def name(self, name):
        cv.check_type("name", name, str)
        self._name = name

    @atomic_number.setter
    def atomic_number(self, atomic_number):
        cv.check_type("atomic number", atomic_number, Integral)
        cv.check_greater_than("atomic number", atomic_number, 0, True)
        self._atomic_number = atomic_number

    @mass_number.setter
    def mass_number(self, mass_number):
        cv.check_type("mass number", mass_number, Integral)
        cv.check_greater_than("mass number", mass_number, 0, True)
        self._mass_number = mass_number

    @metastable.setter
    def metastable(self, metastable):
        cv.check_type("metastable", metastable, Integral)
        cv.check_greater_than("metastable", metastable, 0, True)
        self._metastable = metastable

    @atomic_weight_ratio.setter
    def atomic_weight_ratio(self, atomic_weight_ratio):
        cv.check_type("atomic weight ratio", atomic_weight_ratio, Real)
        cv.check_greater_than("atomic weight ratio", atomic_weight_ratio, 0.0)
        self._atomic_weight_ratio = atomic_weight_ratio

    @classmethod
    def from_endf(cls, ev_or_filename):
        """Generate incident charged particle data from an ENDF evaluation

        Parameters
        ----------
        ev_or_filename : openmc.data.endf.Evaluation or str
            ENDF evaluation to read from. If given as a string, it is assumed to
            be the filename for the ENDF file.


        Returns
        -------
        openmc.data.IncidentChargedParticle
            Proton interaction data

        """

        if isinstance(ev_or_filename, Evaluation):
            ev = ev_or_filename
        else:
            ev = Evaluation(ev_or_filename)

        atomic_number = ev.target["atomic_number"]
        mass_number = ev.target["mass_number"]
        metastable = ev.target["isomeric_state"]
        atomic_weight_ratio = ev.target["mass"]

        element = ATOMIC_SYMBOL[atomic_number]
        name = "{}{}".format(element, mass_number)

        data = cls(name, atomic_number, mass_number, metastable, atomic_weight_ratio)

        # Read each reaction
        for mf, mt, nc, mod in ev.reaction_list:
            if mf == 3:
                data.reactions[mt] = ChargedParticleReaction.from_endf(ev, mt)

        data._evaluation = ev
        return data

    @classmethod
    def from_ace(cls, ace_or_filename, metastable_scheme="nndc"):
        """Generate incident proton continuous-energy data from an ACE table

        Parameters
        ----------
        ace_or_filename : openmc.data.ace.Table or str
            ACE table to read from. If the value is a string, it is assumed to
            be the filename for the ACE file.
        metastable_scheme : {'nndc', 'mcnp'}
            Determine how ZAID identifiers are to be interpreted in the case of
            a metastable nuclide. Because the normal ZAID (=1000*Z + A) does not
            encode metastable information, different conventions are used among
            different libraries. In MCNP libraries, the convention is to add 400
            for a metastable nuclide except for Am242m, for which 95242 is
            metastable and 95642 (or 1095242 in newer libraries) is the ground
            state. For NNDC libraries, ZAID is given as 1000*Z + A + 100*m.

        Returns
        -------
        openmc.data.IncidentChargedParticle
            Incident proton continuous-energy data

        """

        # First obtain the data for the first provided ACE table/file
        if isinstance(ace_or_filename, Table):
            ace = ace_or_filename
        else:
            ace = get_table(ace_or_filename)

        # If mass number hasn't been specified, make an educated guess
        zaid, xs = ace.name.split(".")
        if not xs.endswith("h"):
            raise TypeError(
                "{} is not a continuous-energy proton ACE table.".format(ace)
            )
        name, element, Z, mass_number, metastable = get_metadata(
            int(zaid), metastable_scheme
        )

        data = cls(name, Z, mass_number, metastable, ace.atomic_weight_ratio)

        # Read energy grid
        n_energy = ace.nxs[3]
        i = ace.jxs[1]
        energy = ace.xss[i : i + n_energy] * EV_PER_MEV
        data.energy = energy

        total_xs = ace.xss[i + n_energy : i + 2 * n_energy]
        absorption_xs = ace.xss[i + 2 * n_energy : i + 3 * n_energy]
        heating_number = ace.xss[i + 4 * n_energy : i + 5 * n_energy] * EV_PER_MEV

        # Create redundant reactions (total, absorption, and heating)
        total = ChargedParticleReaction(1)
        total.xs = Tabulated1D(energy, total_xs)
        total.redundant = True
        data.reactions[1] = total

        if np.count_nonzero(absorption_xs) > 0:
            absorption = ChargedParticleReaction(101)
            absorption.xs = Tabulated1D(energy, absorption_xs)
            absorption.redundant = True
            data.reactions[101] = absorption

        heating = ChargedParticleReaction(301)
        heating.xs = Tabulated1D(energy, heating_number * total_xs)
        heating.redundant = True
        data.reactions[301] = heating

        # Read each reaction
        n_reaction = ace.nxs[4] + 1
        for i in range(n_reaction):
            rx = ChargedParticleReaction.from_ace(ace, i)
            data.reactions[rx.mt] = rx

        return data

    def export_to_hdf5(self, path, mode="a", libver="earliest"):
        """Export incident proton data to an HDF5 file.

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

        # If data come from ENDF, don't allow exporting to HDF5
        if hasattr(self, "_evaluation"):
            raise NotImplementedError(
                "Cannot export incident neutron data that "
                "originated from an ENDF file."
            )

        # Open file and write version
        f = h5py.File(str(path), mode, libver=libver)
        f.attrs["filetype"] = np.string_("data_proton")
        if "version" not in f.attrs:
            f.attrs["version"] = np.array(HDF5_VERSION)

        group = f.create_group(self.name)
        group.attrs["Z"] = self.atomic_number
        group.attrs["A"] = self.mass_number
        group.attrs["metastable"] = self.metastable
        group.attrs["atomic_weight_ratio"] = self.atomic_weight_ratio

        # Determine union energy grid
        union_grid = np.array([])
        for rx in self:
            union_grid = np.union1d(union_grid, rx.xs.x)
        group.create_dataset("energy", data=union_grid)

        # Write cross sections
        rxs_group = group.create_group("reactions")
        for mt, rx in self.reactions.items():
            if not rx.redundant:
                rx_group = rxs_group.create_group("reaction_{:03}".format(rx.mt))
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
        openmc.data.IncidentChargedParticle
            Proton interaction data

        """
        if isinstance(group_or_filename, h5py.Group):
            group = group_or_filename
        else:
            h5file = h5py.File(str(group_or_filename), "r")

            # Make sure version matches
            if "version" in h5file.attrs:
                major, minor = h5file.attrs["version"]
                # For now all versions of HDF5 data can be read
            else:
                raise IOError(
                    "HDF5 data does not indicate a version. Your installation of "
                    "the OpenMC Python API expects version {}.x data.".format(
                        HDF5_VERSION_MAJOR
                    )
                )

            group = list(h5file.values())[0]

        name = group.name[1:]
        atomic_number = group.attrs["Z"]
        mass_number = group.attrs["A"]
        metastable = group.attrs["metastable"]
        atomic_weight_ratio = group.attrs["atomic_weight_ratio"]

        data = cls(name, atomic_number, mass_number, metastable, atomic_weight_ratio)

        # Read energy grid
        data.energy = group["energy"][()]

        # Read reaction data
        rxs_group = group["reactions"]
        for name, obj in sorted(rxs_group.items()):
            if name.startswith("reaction_"):
                rx = ChargedParticleReaction.from_hdf5(obj, data.energy)
                data.reactions[rx.mt] = rx

        return data

    @classmethod
    def from_njoy(cls, filename, evaluation=None, **kwargs):
        """Generate incident proton data by running NJOY.

        Parameters
        ----------
        filename : str
            Path to ENDF file
        evaluation : openmc.data.endf.Evaluation, optional
            If the ENDF file contains multiple material evaluations, this
            argument indicates which evaluation to use.
        **kwargs
            Keyword arguments passed to :func:`openmc.data.njoy.make_ace`

        Returns
        -------
        data : openmc.data.IncidentChargedParticle
            Incident proton continuous-energy data

        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run NJOY to create an ACE library
            kwargs.setdefault("output_dir", tmpdir)
            for key in ("acer", "pendf", "heatr", "broadr", "gaspr", "purr"):
                kwargs.setdefault(key, os.path.join(kwargs["output_dir"], key))
            kwargs["evaluation"] = evaluation
            kwargs["pendf"] = False
            kwargs["broadr"] = False
            kwargs["heatr"] = False
            kwargs["gaspr"] = False
            kwargs["purr"] = False
            make_ace(filename, [0.0], **kwargs)

            # Create instance from ACE tables within library
            lib = Library(kwargs["acer"])
            data = cls.from_ace(lib.tables[0])

        return data
