from numbers import Integral, Real

import numpy as np

import openmc.checkvalue as cv
import openmc.data

# Supported keywords for continuous-energy cross section plotting
PLOT_TYPES = ['total', 'scatter', 'elastic', 'inelastic', 'fission',
              'absorption', 'capture', 'nu-fission', 'nu-scatter', 'unity',
              'slowing-down power', 'damage']
# Supported keywoards for multi-group cross section plotting
PLOT_TYPES_MGXS = ['total', 'absorption', 'fission', 'kappa-fission',
                   'chi', 'chi-prompt', 'nu-fission', 'prompt-nu-fission',
                   'inverse-velocity', 'unity']

# Special MT values
UNITY_MT = -1
XI_MT = -2

# MTs to combine to generate associated plot_types
PLOT_TYPES_MT = {'total': (2, 3,),
                 'scatter': (2, 4, 11, 16, 17, 22, 23, 24, 25, 28, 29, 30,
                             32, 33, 34, 35, 36, 37, 41, 42, 44, 45, 152,
                             153, 154, 156, 157, 158, 159, 160, 161, 162,
                             163, 164, 165, 166, 167, 168, 169, 170, 171,
                             172, 173, 174, 175, 176, 177, 178, 179, 180,
                             181, 183, 184, 190, 194, 196, 198, 199, 200,
                             875, 891),
                 'elastic': (2,),
                 'inelastic': (4, 11, 16, 17, 22, 23, 24, 25, 28, 29, 30,
                               32, 33, 34, 35, 36, 37, 41, 42, 44, 45, 152,
                               153, 154, 156, 157, 158, 159, 160, 161, 162,
                               163, 164, 165, 166, 167, 168, 169, 170, 171,
                               172, 173, 174, 175, 176, 177, 178, 179, 180,
                               181, 183, 184, 190, 194, 196, 198, 199, 200,
                               875, 891),
                 'fission': (18,),
                 'absorption': (27,), 'capture': (101,),
                 'nu-fission': (18,),
                 'nu-scatter': (2, 4, 11, 16, 17, 22, 23, 24, 25, 28, 29, 30,
                                32, 33, 34, 35, 36, 37, 41, 42, 44, 45, 152,
                                153, 154, 156, 157, 158, 159, 160, 161, 162,
                                163, 164, 165, 166, 167, 168, 169, 170, 171,
                                172, 173, 174, 175, 176, 177, 178, 179, 180,
                                181, 183, 184, 190, 194, 196, 198, 199, 200,
                                875, 891),
                 'unity': (UNITY_MT,),
                 'slowing-down power': (2, 4, 11, 16, 17, 22, 23, 24, 25, 28,
                                        29, 30, 32, 33, 34, 35, 36, 37, 41, 42,
                                        44, 45, 152, 153, 154, 156, 157, 158,
                                        159, 160, 161, 162, 163, 164, 165, 166,
                                        167, 168, 169, 170, 171, 172, 173, 174,
                                        175, 176, 177, 178, 179, 180, 181, 183,
                                        184, 190, 194, 196, 198, 199, 200, 875,
                                        891, XI_MT),
                 'damage': (444,)}
# Operations to use when combining MTs the first np.add is used in reference
# to zero
PLOT_TYPES_OP = {'total': (np.add,),
                 'scatter': (np.add,) * (len(PLOT_TYPES_MT['scatter']) - 1),
                 'elastic': (),
                 'inelastic': (np.add,) * (len(PLOT_TYPES_MT['inelastic']) - 1),
                 'fission': (), 'absorption': (),
                 'capture': (), 'nu-fission': (),
                 'nu-scatter': (np.add,) * (len(PLOT_TYPES_MT['nu-scatter']) - 1),
                 'unity': (),
                 'slowing-down power':
                    (np.add,) * (len(PLOT_TYPES_MT['slowing-down power']) - 2) + (np.multiply,),
                 'damage': ()}

# Whether or not to multiply the reaction by the yield as well
PLOT_TYPES_YIELD = {'total': (False, False),
                    'scatter': (False,) * len(PLOT_TYPES_MT['scatter']),
                    'elastic': (False,),
                    'inelastic': (False,) * len(PLOT_TYPES_MT['inelastic']),
                    'fission': (False,), 'absorption': (False,),
                    'capture': (False,), 'nu-fission': (True,),
                    'nu-scatter': (True,) * len(PLOT_TYPES_MT['nu-scatter']),
                    'unity': (False,),
                    'slowing-down power':
                        (True,) * len(PLOT_TYPES_MT['slowing-down power']),
                    'damage': (False,)}

# Types of plots to plot linearly in y
PLOT_TYPES_LINEAR = {'nu-fission / fission', 'nu-scatter / scatter',
                     'nu-fission / absorption', 'fission / absorption'}


def plot_xs(this, types, divisor_types=None, temperature=294., axis=None,
            energy_range=(1.E-5, 20.E6), sab_name=None, cross_sections=None,
            enrichment=None, **kwargs):
    """Creates a figure of continuous-energy cross sections for this item

    Parameters
    ----------
    this : {openmc.Element, openmc.Nuclide, openmc.Material}
        Object to source data from
    types : Iterable of values of PLOT_TYPES
        The type of cross sections to include in the plot.
    divisor_types : Iterable of values of PLOT_TYPES, optional
        Cross section types which will divide those produced by types
        before plotting. A type of 'unity' can be used to effectively not
        divide some types.
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    axis : matplotlib.axes, optional
        A previously generated axis to use for plotting. If not specified,
        a new axis and figure will be generated.
    energy_range : tuple of floats
        Energy range (in eV) to plot the cross section within
    sab_name : str, optional
        Name of S(a,b) library to apply to MT=2 data when applicable; only used
        for items which are instances of openmc.Element or openmc.Nuclide
    cross_sections : str, optional
        Location of cross_sections.xml file. Default is None.
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None. This is only used for
        items which are instances of openmc.Element
    **kwargs
        All keyword arguments are passed to
        :func:`matplotlib.pyplot.figure`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        If axis is None, then a Matplotlib Figure of the generated
        cross section will be returned. Otherwise, a value of
        None will be returned as the figure and axes have already been
        generated.

    """

    from matplotlib import pyplot as plt

    if isinstance(this, openmc.Nuclide):
        data_type = 'nuclide'
    elif isinstance(this, openmc.Element):
        data_type = 'element'
    elif isinstance(this, openmc.Material):
        data_type = 'material'
    else:
        raise TypeError("Invalid type for plotting")

    E, data = calculate_xs(this, types, temperature, sab_name, cross_sections,
                           enrichment)

    if divisor_types:
        cv.check_length('divisor types', divisor_types, len(types),
                        len(types))
        Ediv, data_div = calculate_xs(this, divisor_types, temperature,
                                      sab_name, cross_sections, enrichment)

        # Create a new union grid, interpolate data and data_div on to that
        # grid, and then do the actual division
        Enum = E[:]
        E = np.union1d(Enum, Ediv)
        if data_type == 'nuclide':
            data_new = []
        else:
            data_new = np.zeros((len(types), len(E)))

        for line in range(len(types)):
            if data_type == 'nuclide':
                data_new.append(openmc.data.Combination([data[line],
                                                         data_div[line]],
                                                        [np.divide]))
            else:
                data_new[line, :] = \
                    np.divide(np.interp(E, Enum, data[line, :]),
                              np.interp(E, Ediv, data_div[line, :]))
            if divisor_types[line] != 'unity':
                types[line] = types[line] + ' / ' + divisor_types[line]
        data = data_new

    # Generate the plot
    if axis is None:
        fig = plt.figure(**kwargs)
        ax = fig.add_subplot(111)
    else:
        fig = None
        ax = axis
    # Set to loglog or semilogx depending on if we are plotting a data
    # type which we expect to vary linearly
    if set(types).issubset(PLOT_TYPES_LINEAR):
        plot_func = ax.semilogx
    else:
        plot_func = ax.loglog
    # Plot the data
    for i in range(len(data)):
        if data_type == 'nuclide':
            to_plot = data[i](E)
        else:
            to_plot = data[i, :]
        if np.sum(to_plot) > 0.:
            plot_func(E, to_plot, label=types[i])

    ax.set_xlabel('Energy [eV]')
    if divisor_types:
        if data_type == 'nuclide':
            ylabel = 'Nuclidic Microscopic Data'
        elif data_type == 'element':
            ylabel = 'Elemental Microscopic Data'
        elif data_type == 'material':
            ylabel = 'Macroscopic Data'
    else:
        if data_type == 'nuclide':
            ylabel = 'Microscopic Cross Section [b]'
        elif data_type == 'element':
            ylabel = 'Elemental Cross Section [b]'
        elif data_type == 'material':
            ylabel = 'Macroscopic Cross Section [1/cm]'
    ax.set_ylabel(ylabel)
    ax.legend(loc='best')
    ax.set_xlim(energy_range)
    if this.name is not None:
        title = 'Cross Section for ' + this.name
        ax.set_title(title)

    return fig


def plot_mgxs(this, types, divisor_types=None, temperature=294., axis=None,
              energy_range=(1.E-5, 20.E6), cross_sections=None,
              ce_cross_sections=None, enrichment=None, **kwargs):
    """Creates a figure of multi-group cross sections for this item

    Parameters
    ----------
    this : {openmc.Element, openmc.Nuclide, openmc.Material, openmc.Macroscopic}
        Object to source data from
    types : Iterable of values of PLOT_TYPES
        The type of cross sections to include in the plot.
    divisor_types : Iterable of values of PLOT_TYPES, optional
        Cross section types which will divide those produced by types
        before plotting. A type of 'unity' can be used to effectively not
        divide some types.
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    axis : matplotlib.axes, optional
        A previously generated axis to use for plotting. If not specified,
        a new axis and figure will be generated.
    energy_range : tuple of floats
        Energy range (in eV) to plot the cross section within
    cross_sections : str, optional
        Location of MGXS HDF5 Library file. Default is None.
    cross_sections : str, optional
        Location of continuous-energy cross_sections.xml file. Default is None.
        This is used only for expanding an openmc.Element object passed as this
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None. This is only used for
        items which are instances of openmc.Element
    **kwargs
        All keyword arguments are passed to
        :func:`matplotlib.pyplot.figure`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        If axis is None, then a Matplotlib Figure of the generated
        cross section will be returned. Otherwise, a value of
        None will be returned as the figure and axes have already been
        generated.

    """

    from matplotlib import pyplot as plt

    if isinstance(this, openmc.Nuclide):
        data_type = 'nuclide'
    elif isinstance(this, openmc.Element):
        data_type = 'element'
    elif isinstance(this, openmc.Material):
        data_type = 'material'
    elif isinstance(this, openmc.Macroscopic):
        data_type = 'macroscopic'
    else:
        raise TypeError("Invalid type for plotting")

    E, data = calculate_mgxs(this, types, temperature, cross_sections,
                             ce_cross_sections, enrichment)

    if divisor_types:
        cv.check_length('divisor types', divisor_types, len(types),
                        len(types))
        Ediv, data_div = calculate_mgxs(this, divisor_types, temperature,
                                        cross_sections, ce_cross_sections,
                                        enrichment)

        # Perform the division
        for line in range(len(types)):
            data[line, :] = np.divide(data[line, :], data_div[line, :])
            if divisor_types[line] != 'unity':
                types[line] = types[line] + ' / ' + divisor_types[line]

    # Modify the data such that the it plots the group values as distinct
    # horizontal lines in an efficient manner
    E_plot = []
    data_plot = []
    for line in range(len(types)):
        for g in range(len(E) - 1):
            if line == 0:
                E_plot.append(E[g])
                E_plot.append(E[g + 1])
                E_plot.append(None)
            data_plot.append(data[line, g])
            data_plot.append(data[line, g])
            data_plot.append(None)

    # Generate the plot
    if axis is None:
        fig = plt.figure(**kwargs)
        ax = fig.add_subplot(111)
    else:
        fig = None
        ax = axis
    # Set to loglog or semilogx depending on if we are plotting a data
    # type which we expect to vary linearly
    if set(types).issubset(PLOT_TYPES_LINEAR):
        plot_func = ax.semilogx
    else:
        plot_func = ax.loglog
    # Plot the data
    for i in range(len(data_plot)):
        plot_func(E_plot, data_plot[i, :], label=types[i])

    ax.set_xlabel('Energy [eV]')
    if divisor_types:
        if data_type == 'nuclide':
            ylabel = 'Nuclidic Microscopic Data'
        elif data_type == 'element':
            ylabel = 'Elemental Microscopic Data'
        elif data_type == 'material':
            ylabel = 'Macroscopic Data'
    else:
        if data_type == 'nuclide':
            ylabel = 'Microscopic Cross Section [b]'
        elif data_type == 'element':
            ylabel = 'Elemental Cross Section [b]'
        elif data_type == 'material':
            ylabel = 'Macroscopic Cross Section [1/cm]'
    ax.set_ylabel(ylabel)
    ax.legend(loc='best')
    ax.set_xlim(energy_range)
    if this.name is not None:
        title = 'Cross Section for ' + this.name
        ax.set_title(title)

    return fig


def calculate_xs(this, types, temperature=294., sab_name=None,
                 cross_sections=None, enrichment=None):
    """Calculates continuous-energy cross sections of a requested type

    Parameters
    ----------
    this : openmc.Element, openmc.Nuclide, or openmc.Material
        Object to source data from
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
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None
        (natural composition).

    Returns
    -------
    energy_grid : numpy.array
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Cross sections calculated at the energy grid described by energy_grid

    """

    # Check types
    cv.check_type('temperature', temperature, Real)
    if sab_name:
        cv.check_type('sab_name', sab_name, str)
    if enrichment:
        cv.check_type('enrichment', enrichment, Real)

    if isinstance(this, openmc.Nuclide):
        energy_grid, data = _calculate_xs_nuclide(this, types, temperature,
                                                  sab_name, cross_sections)
    elif isinstance(this, openmc.Element):
        energy_grid, data = _calculate_xs_element(this, types, temperature,
                                                  sab_name, cross_sections,
                                                  enrichment)
    elif isinstance(this, openmc.Material):
        energy_grid, data = _calculate_xs_material(this, types, temperature,
                                                   cross_sections)
    else:
        raise TypeError("Invalid type")

    return energy_grid, data


def _calculate_xs_element(this, types, temperature=294., sab_name=None,
                          cross_sections=None, enrichment=None):
    """Calculates continuous-energy cross sections of a requested type

    Parameters
    ----------
    this : openmc.Element
        Element object to source data from
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
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None
        (natural composition).

    Returns
    -------
    energy_grid : numpy.array
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Macroscopic cross sections calculated at the energy grid described
        by energy_grid

    """

    # Load the library
    library = openmc.data.DataLibrary.from_xml(cross_sections)

    # Expand elements in to nuclides with atomic densities
    nuclides = this.expand(100., 'ao', enrichment=enrichment,
                           cross_sections=cross_sections)

    # For ease of processing split out nuc and nuc_density
    nuc_fractions = [nuclide[1] for nuclide in nuclides]

    # Identify the nuclides which have S(a,b) data
    sabs = {}
    for nuclide in nuclides:
        sabs[nuclide[0].name] = None
    if sab_name:
        sab = openmc.data.ThermalScattering.from_hdf5(sab_name)
        for nuc in sab.nuclides:
            sabs[nuc] = library.get_by_material(sab_name)['path']

    # Now we can create the data sets to be plotted
    xs = []
    E = []
    for nuclide in nuclides:
        sab_tab = sabs[nuclide[0].name]
        temp_E, temp_xs = calculate_xs(nuclide[0], types, temperature, sab_tab,
                                       cross_sections)
        E.append(temp_E)
        xs.append(temp_xs)

    # Condense the data for every nuclide
    # First create a union energy grid
    energy_grid = E[0]
    for n in range(1, len(E)):
        energy_grid = np.union1d(energy_grid, E[n])

    # Now we can combine all the nuclidic data
    data = np.zeros((len(types), len(energy_grid)))
    for line in range(len(types)):
        if types[line] == 'unity':
            data[line, :] = 1.
        else:
            for n in range(len(nuclides)):
                data[line, :] += nuc_fractions[n] * xs[n][line](energy_grid)

    return energy_grid, data


def _calculate_xs_nuclide(this, types, temperature=294., sab_name=None,
                          cross_sections=None):
    """Calculates continuous-energy cross sections of a requested type

    Parameters
    ----------
    this : openmc.Nuclide
        Nuclide object to source data from
    types : Iterable of str or Integral
        The type of cross sections to calculate; values can either be those
        in openmc.PLOT_TYPES or integers which correspond to reaction
        channel (MT) numbers.
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
    energy_grid : numpy.array
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Cross sections calculated at the energy grid described by
        energy_grid

    """

    # Parse the types
    mts = []
    ops = []
    yields = []
    for line in types:
        if line in PLOT_TYPES:
            mts.append(PLOT_TYPES_MT[line])
            yields.append(PLOT_TYPES_YIELD[line])
            ops.append(PLOT_TYPES_OP[line])
        else:
            # Not a built-in type, we have to parse it ourselves
            cv.check_type('MT in types', line, Integral)
            cv.check_greater_than('MT in types', line, 0)
            mts.append((line,))
            yields.append((False,))
            ops.append(())

    # Load the library
    library = openmc.data.DataLibrary.from_xml(cross_sections)

    # Convert temperature to format needed for access in the library
    strT = "{}K".format(int(round(temperature)))
    T = temperature

    # Now we can create the data sets to be plotted
    energy_grid = []
    xs = []
    lib = library.get_by_material(this.name)
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
            energy_grid = grid
        else:
            energy_grid = nuc.energy[nucT]

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
                elif mt == UNITY_MT:
                    funcs.append(lambda x: 1.)
                elif mt == XI_MT:
                    awr = nuc.atomic_weight_ratio
                    alpha = ((awr - 1.) / (awr + 1.))**2
                    xi = 1. + alpha * np.log(alpha) / (1. - alpha)
                    funcs.append(lambda x: xi)
                else:
                    funcs.append(lambda x: 0.)
            xs.append(openmc.data.Combination(funcs, op))
    else:
        raise ValueError(this.name + " not in library")

    return energy_grid, xs


def _calculate_xs_material(this, types, temperature=294., cross_sections=None):
    """Calculates continuous-energy macroscopic cross sections of a
    requested type

    Parameters
    ----------
    this : openmc.Material
        Material object to source data from
    types : Iterable of values of PLOT_TYPES
        The type of cross sections to calculate
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    cross_sections : str, optional
        Location of cross_sections.xml file. Default is None.

    Returns
    -------
    energy_grid : numpy.array
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Macroscopic cross sections calculated at the energy grid described
        by energy_grid

    """

    if this.temperature is not None:
        T = this.temperature
    else:
        T = temperature

    # Load the library
    library = openmc.data.DataLibrary.from_xml(cross_sections)

    # Expand elements in to nuclides with atomic densities
    nuclides = this.get_nuclide_atom_densities(cross_sections)

    # For ease of processing split out nuc and nuc_density
    nuc_densities = [nuclide[1][1] for nuclide in nuclides.items()]

    # Identify the nuclides which have S(a,b) data
    sabs = {}
    for nuclide in nuclides.items():
        sabs[nuclide[0].name] = None
    for sab_name in this._sab:
        sab = openmc.data.ThermalScattering.from_hdf5(
            library.get_by_material(sab_name)['path'])
        for nuc in sab.nuclides:
            sabs[nuc] = library.get_by_material(sab_name)['path']

    # Now we can create the data sets to be plotted
    xs = []
    E = []
    for nuclide in nuclides.items():
        sab_tab = sabs[nuclide[0].name]
        temp_E, temp_xs = calculate_xs(nuclide[0], types, T, sab_tab,
                                       cross_sections)
        E.append(temp_E)
        xs.append(temp_xs)

    # Condense the data for every nuclide
    # First create a union energy grid
    energy_grid = E[0]
    for n in range(1, len(E)):
        energy_grid = np.union1d(energy_grid, E[n])

    # Now we can combine all the nuclidic data
    data = np.zeros((len(types), len(energy_grid)))
    for line in range(len(types)):
        if types[line] == 'unity':
            data[line, :] = 1.
        else:
            for n in range(len(nuclides)):
                data[line, :] += nuc_densities[n] * xs[n][line](energy_grid)

    return energy_grid, data


def calculate_mgxs(this, types, temperature=294., cross_sections=None,
                   ce_cross_sections=None, enrichment=None):
    """Calculates continuous-energy cross sections of a requested type

    Parameters
    ----------
    this : openmc.Element, openmc.Nuclide, or openmc.Material
        Object to source data from
    types : Iterable of values of PLOT_TYPES
        The type of cross sections to calculate
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    cross_sections : str, optional
        Location of MGXS HDF5 Library file. Default is None.
    ce_cross_sections : str, optional
        Location of continuous-energy cross_sections.xml file. Default is None.
        This is used only for expanding an openmc.Element object passed as this
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None
        (natural composition).

    Returns
    -------
    numpy.array
        Energy group structure in units of eV
    data : numpy.ndarray
        Multi-Group Cross sections

    """

    # Check types
    cv.check_type('temperature', temperature, Real)
    if enrichment:
        cv.check_type('enrichment', enrichment, Real)

    cv.check_type("cross_sections", cross_sections, str)
    library = openmc.MGXSLibrary.from_hdf5(cross_sections)

    if isinstance(this, (openmc.Nuclide, openmc.Macroscopic)):
        data = _calculate_mgxs_nuc_macro(this, types, library, temperature)
    elif isinstance(this, (openmc.Element, openmc.Material)):
        data = _calculate_mgxs_elem_mat(this, types, library, temperature,
                                        ce_cross_sections, enrichment)
    else:
        raise TypeError("Invalid type")

    return library.energy_groups.group_structure, data


def _calculate_mgxs_nuc_macro(this, types, library, temperature=294.):
    """Determines the multi-group cross sections of a nuclide or macroscopic
    object

    Parameters
    ----------
    this : {openmc.Nuclide, openmc.Macroscopic}
        Object to source data from
    types : Iterable of str
        The type of cross sections to calculate; values can either be those
        in openmc.PLOT_TYPES_MGXS
    library : openmc.MGXSLibrary
        MGXS Library containing the data of interest
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.

    Returns
    -------
    data : numpy.ndarray
        Cross sections values of the requested types

    """

    # Check the parameters
    for line in types:
        cv.check_type("line", line, str)
        cv.check_value("line", line, PLOT_TYPES_MGXS)

    xsdata = library[this.name]

    if xsdata is not None:
        # Obtain the nearest temperature
        data_Ts = library.temperatures
        for t in range(len(data_Ts)):
            # Take off the "K" and convert to a float
            data_Ts[t] = float(data_Ts[t][:-1])
        min_delta = np.finfo(np.float64).max
        closest_t = -1
        for t in data_Ts:
            if abs(data_Ts[t] - temperature) < min_delta:
                closest_t = t
        t = closest_t

        # Get the data
        data = np.zeros((len(types), library.energy_groups.num_groups))
        for line in types:
            if line == 'unity':
                data[line, :] = 1.
            else:
                attr = line.replace(' ', '_').replace('-', '_')
                data[line, :] = getattr(xsdata, attr)[t]
    else:
        raise ValueError("{} not present in provided MGXS "
                         "library".format(this.name))

    return data


def _calculate_mgxs_elem_mat(this, types, library, temperature=294.,
                             ce_cross_sections=None, enrichment=None):
    """Determines the multi-group cross sections of an element or material
    object

    Parameters
    ----------
    this : {openmc.Element, openmc.Material}
        Object to source data from
    types : Iterable of str
        The type of cross sections to calculate; values can either be those
        in openmc.PLOT_TYPES_MGXS
    library : openmc.MGXSLibrary
        MGXS Library containing the data of interest
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    ce_cross_sections : str, optional
        Location of continuous-energy cross_sections.xml file. Default is None.
        This is used only for expanding the elements
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None
        (natural composition).

    Returns
    -------
    data : numpy.ndarray
        Cross sections values of the requested types

    """

    # Check the parameters
    for line in types:
        cv.check_type("line", line, str)
        cv.check_value("line", line, PLOT_TYPES_MGXS)

    if isinstance(this, openmc.Material):
        if this.temperature is not None:
            T = this.temperature
        else:
            T = temperature
        # Expand elements in to nuclides with atomic densities
        nuclides = this.get_nuclide_atom_densities(ce_cross_sections)

        # For ease of processing split out nuc and nuc_density
        nuc_multiplier = [nuclide[1][1] for nuclide in nuclides.items()]
    else:
        T = temperature
        # Expand elements in to nuclides with atomic densities
        nuclides = this.expand(100., 'ao', enrichment=enrichment,
                               cross_sections=ce_cross_sections)

        # For ease of processing split out nuc and nuc_fractions
        nuc_multiplier = [nuclide[1] for nuclide in nuclides]

    nuc_data = []
    for nuclide in nuclides:
        nuc_data.append(_calculate_mgxs_nuc_macro(nuclide[0], types, library,
                                                  T))

    # Combine across the nuclides
    data = np.zeros((len(types), library.energy_groups.num_groups))
    for line in range(len(types)):
        if types[line] == 'unity':
            data[line, :] = 1.
        else:
            for n in range(len(nuclides)):
                data[line, :] += nuc_multiplier[n] * nuc_data[n][line, :]

    return data
