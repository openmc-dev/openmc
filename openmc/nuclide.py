from numbers import Integral, Real
import sys
import warnings
import os

from six import string_types

import openmc.checkvalue as cv
from openmc.plot_data import *
import openmc.data


class Nuclide(object):
    """A nuclide that can be used in a material.

    Parameters
    ----------
    name : str
        Name of the nuclide, e.g. U235

    Attributes
    ----------
    name : str
        Name of the nuclide, e.g. U235
    scattering : 'data' or 'iso-in-lab' or None
        The type of angular scattering distribution to use

    """

    def __init__(self, name=''):
        # Initialize class attributes
        self._name = ''
        self._scattering = None

        # Set the Material class attributes
        self.name = name

    def __eq__(self, other):
        if isinstance(other, Nuclide):
            if self.name != other.name:
                return False
            else:
                return True
        elif isinstance(other, string_types) and other == self.name:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self == other

    def __gt__(self, other):
        return repr(self) > repr(other)

    def __lt__(self, other):
        return not self > other

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        string = 'Nuclide    -    {0}\n'.format(self._name)
        if self.scattering is not None:
            string += '{0: <16}{1}{2}\n'.format('\tscattering', '=\t',
                                                self.scattering)
        return string

    @property
    def name(self):
        return self._name

    @property
    def scattering(self):
        return self._scattering

    @name.setter
    def name(self, name):
        cv.check_type('name', name, string_types)
        self._name = name

        if '-' in name:
            self._name = name.replace('-', '')
            self._name = self._name.replace('Nat', '0')
            if self._name.endswith('m'):
                self._name = self._name[:-1] + '_m1'

            msg = 'OpenMC nuclides follow the GND naming convention. Nuclide ' \
                  '"{}" is being renamed as "{}".'.format(name, self._name)
            warnings.warn(msg)

    @scattering.setter
    def scattering(self, scattering):
        if not scattering in ['data', 'iso-in-lab', None]:
            msg = 'Unable to set scattering for Nuclide to {0} which ' \
                  'is not "data", "iso-in-lab", or None'.format(scattering)
            raise ValueError(msg)

        self._scattering = scattering

    def plot_xs(self, types, divisor_types=None, temperature=294.,
                Erange=(1.E-5, 20.E6), sab_name=None, cross_sections=None,
                **kwargs):
        """Creates a figure of continuous-energy cross sections for this
        nuclide

        Parameters
        ----------
        types : Iterable of values of PLOT_TYPES
            The type of cross sections to include in the plot
        divisor_types : Iterable of values of PLOT_TYPES, optional
            Cross section types which will divide those produced by types
            before plotting. A type of 'unity' can be used to effectively not
            divide some types.
        temperature : float, optional
            Temperature in Kelvin to plot. If not specified, a default
            temperature of 294K will be plotted. Note that the nearest
            temperature in the library for each nuclide will be used as opposed
            to using any interpolation.
        Erange : tuple of floats
            Energy range (in eV) to plot the cross section within
        sab_name : str, optional
            Name of S(a,b) library to apply to MT=2 data when applicable.
        cross_sections : str, optional
            Location of cross_sections.xml file. Default is None.
        **kwargs
            All keyword arguments are passed to
            :func:`matplotlib.pyplot.figure`.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Matplotlib Figure of the generated macroscopic cross section

        """

        from matplotlib import pyplot as plt

        E, data = self.calculate_xs(types, temperature, sab_name,
                                    cross_sections)

        if divisor_types:
            cv.check_length('divisor types', divisor_types, len(types),
                            len(types))
            Ediv, data_div = self.calculate_xs(divisor_types, temperature,
                                               sab_name, cross_sections)

            # Create a new union grid, interpolate data and data_div on to that
            # grid, and then do the actual division
            Enum = E[:]
            E = np.union1d(Enum, Ediv)
            data_new = []

            for l in range(len(types)):
                data_new.append(openmc.data.Combination([data[l], data_div[l]],
                                                        [np.divide]))
                if divisor_types[l] != 'unity':
                    types[l] = types[l] + ' / ' + divisor_types[l]
            data = data_new

        # Generate the plot
        fig = plt.figure(**kwargs)
        ax = fig.add_subplot(111)
        for i in range(len(data)):
            to_plot = data[i](E)
            if np.sum(to_plot) > 0.:
                ax.loglog(E, to_plot, label=types[i])

        ax.set_xlabel('Energy [eV]')
        ax.set_ylabel('Microscopic Cross Section [b]')
        ax.legend(loc='best')
        ax.set_xlim(Erange)
        if self.name is not None:
            title = 'Microscopic Cross Section for ' + self.name
            ax.set_title(title)

        return fig

    def calculate_xs(self, types, temperature=294., sab_name=None,
                     cross_sections=None):
        """Calculates continuous-energy cross sections of a requested type

        Parameters
        ----------
        types : Iterable of values of PLOT_TYPES
            The type of cross sections to calculate
        temperature : float, optional
            Temperature in Kelvin to plot. If not specified, a default
            temperature of 294K will be plotted. Note that the nearest
            temperature in the library for each nuclide will be used as opposed
            to using any interpolation.
        sab_name : str, optional
            Name of S(a,b) library to apply to MT=2 data when applicable.
        cross_sections : str, optional
            Location of cross_sections.xml file. Default is None.

        Returns
        -------
        E : numpy.array
            Energies at which cross sections are calculated, in units of eV
        data : numpy.ndarray
            Cross sections calculated at the energy grid described by unionE

        """

        # Check types
        if cross_sections is not None:
            cv.check_type('cross_sections', cross_sections, str)
        cv.check_iterable_type('types', types, str)
        if sab_name:
            cv.check_type('sab_name', sab_name, str)

        # Parse the types
        mts = []
        ops = []
        yields = []
        for line in types:
            if line in openmc.plot_data.PLOT_TYPES:
                mts.append(PLOT_TYPES_MT[line])
                yields.append(PLOT_TYPES_YIELD[line])
                ops.append(PLOT_TYPES_OP[line])
            else:
                # Not a built-in type, we have to parse it ourselves
                raise NotImplementedError()

        # If cross_sections is None, get the cross sections from the
        # OPENMC_CROSS_SECTIONS environment variable
        if cross_sections is None:
            cross_sections = os.environ.get('OPENMC_CROSS_SECTIONS')

        # If a cross_sections library is present, check natural nuclides
        # against the nuclides in the library
        if cross_sections is not None:
            library = openmc.data.DataLibrary.from_xml(cross_sections)
        else:
            raise ValueError("cross_sections or OPENMC_CROSS_SECTIONS "
                             "environmental variable must be set")

        # Convert temperature to format needed for access in the library
        cv.check_type('temperature', temperature, Real)
        strT = "{}K".format(int(round(temperature)))
        T = temperature

        # Now we can create the data sets to be plotted
        E = []
        xs = []
        lib = library.get_by_material(self.name)
        if lib is not None:
            nuc = openmc.data.IncidentNeutron.from_hdf5(lib['path'])
            # Obtain the nearest temperature
            if strT in nuc.temperatures:
                nucT = strT
            else:
                data_Ts = nuc.temperatures
                for t in range(len(data_Ts)):
                    # Take off the "K" and convert to a float
                    data_Ts[t] = float(data_Ts[t][:-1])
                min_delta = np.finfo(np.float64).max
                closest_t = -1
                for t in data_Ts:
                    if abs(data_Ts[t] - T) < min_delta:
                        closest_t = t
                nucT = "{}K".format(int(round(data_Ts[closest_t])))

            # Prep S(a,b) data if needed
            if sab_name:
                sab = openmc.data.ThermalScattering.from_hdf5(sab_name)
                # Obtain the nearest temperature
                if strT in sab.temperatures:
                    sabT = strT
                else:
                    data_Ts = sab.temperatures
                    for t in range(len(data_Ts)):
                        # Take off the "K" and convert to a float
                        data_Ts[t] = float(data_Ts[t][:-1])
                    min_delta = np.finfo(np.float64).max
                    closest_t = -1
                    for t in data_Ts:
                        if abs(data_Ts[t] - T) < min_delta:
                            closest_t = t
                    sabT = "{}K".format(int(round(data_Ts[closest_t])))

                # Create an energy grid composed the S(a,b) and
                # the nuclide's grid
                grid = nuc.energy[nucT]
                sab_Emax = 0.
                sab_funcs = []
                if sab.elastic_xs:
                    elastic = sab.elastic_xs[sabT]
                    if isinstance(elastic, openmc.data.CoherentElastic):
                        grid = np.union1d(grid, elastic.bragg_edges)
                        if elastic.bragg_edges[-1] > sab_Emax:
                            sab_Emax = elastic.bragg_edges[-1]
                    elif isinstance(elastic, openmc.data.Tabulated1D):
                        grid = np.union1d(grid, elastic.x)
                        if elastic.x[-1] > sab_Emax:
                            sab_Emax = elastic.x[-1]
                    sab_funcs.append(elastic)
                if sab.inelastic_xs:
                    inelastic = sab.inelastic_xs[sabT]
                    grid = np.union1d(grid, inelastic.x)
                    if inelastic.x[-1] > sab_Emax:
                            sab_Emax = inelastic.x[-1]
                    sab_funcs.append(inelastic)
                E = grid
            else:
                E = nuc.energy[nucT]

            for i, mt_set in enumerate(mts):
                # Get the reaction xs data from the nuclide
                funcs = []
                op = ops[i]
                for mt, yield_check in zip(mt_set, yields[i]):
                    if mt == 2:
                        if sab_name:
                            # Then we need to do a piece-wise function of
                            # The S(a,b) and non-thermal data
                            sab_sum = openmc.data.Sum(sab_funcs)
                            pw_funcs = openmc.data.Regions1D(
                                [sab_sum, nuc[mt].xs[nucT]],
                                [sab_Emax])
                            funcs.append(pw_funcs)
                        else:
                            funcs.append(nuc[mt].xs[nucT])
                    elif mt in nuc:
                        if yield_check:
                            found_it = False
                            for prod in nuc[mt].products:
                                if prod.particle == 'neutron' and \
                                    prod.emission_mode == 'total':
                                    func = openmc.data.Combination(
                                        [nuc[mt].xs[nucT], prod.yield_],
                                        [np.multiply])
                                    funcs.append(func)
                                    found_it = True
                                    break
                            if not found_it:
                                for prod in nuc[mt].products:
                                    if prod.particle == 'neutron' and \
                                        prod.emission_mode == 'prompt':
                                        func = openmc.data.Combination(
                                            [nuc[mt].xs[nucT],
                                             prod.yield_], [np.multiply])
                                        funcs.append(func)
                                        found_it = True
                                        break
                            if not found_it:
                                # Assume the yield is 1
                                funcs.append(nuc[mt].xs[nucT])
                        else:
                            funcs.append(nuc[mt].xs[nucT])
                    elif mt == 0:
                        funcs.append(lambda x: 1.)
                    else:
                        funcs.append(lambda x: 0.)
                xs.append(openmc.data.Combination(funcs, op))
        else:
            raise ValueError(nuclide[0] + " not in library")

        return E, xs
