from __future__ import annotations
from itertools import chain
from numbers import Integral, Real
from typing import Dict, Iterable, List

import numpy as np

import openmc.checkvalue as cv
import openmc.data

# Supported keywords for continuous-energy cross section plotting
PLOT_TYPES = {'total', 'scatter', 'elastic', 'inelastic', 'fission',
              'absorption', 'capture', 'nu-fission', 'nu-scatter', 'unity',
              'slowing-down power', 'damage'}

# Supported keywords for multi-group cross section plotting
PLOT_TYPES_MGXS = {'total', 'absorption', 'scatter', 'fission',
                   'kappa-fission', 'nu-fission', 'prompt-nu-fission',
                   'deleyed-nu-fission', 'chi', 'chi-prompt', 'chi-delayed',
                   'inverse-velocity', 'beta', 'decay-rate', 'unity'}
# Create a dictionary which can be used to convert PLOT_TYPES_MGXS to the
# openmc.XSdata attribute name needed to access the data
_PLOT_MGXS_ATTR = {line: line.replace(' ', '_').replace('-', '_')
                   for line in PLOT_TYPES_MGXS}
_PLOT_MGXS_ATTR['scatter'] = 'scatter_matrix'

# Special MT values
UNITY_MT = -1
XI_MT = -2

# MTs to combine to generate associated plot_types
_INELASTIC = [mt for mt in openmc.data.SUM_RULES[3] if mt != 27]
PLOT_TYPES_MT = {
    'total': openmc.data.SUM_RULES[1],
    'scatter': [2] + _INELASTIC,
    'elastic': [2],
    'inelastic': _INELASTIC,
    'fission': [18],
    'absorption': [27],
    'capture': [101],
    'nu-fission': [18],
    'nu-scatter': [2] + _INELASTIC,
    'unity': [UNITY_MT],
    'slowing-down power': [2] + [XI_MT],
    'damage': [444]
}

# Types of plots to plot linearly in y
PLOT_TYPES_LINEAR = {'nu-fission / fission', 'nu-scatter / scatter',
                     'nu-fission / absorption', 'fission / absorption'}

# Minimum and maximum energies for plotting (units of eV)
_MIN_E = 1.e-5
_MAX_E = 20.e6


ELEMENT_NAMES = list(openmc.data.ELEMENT_SYMBOL.values())[1:]


def _get_legend_label(this, type):
    """Gets a label for the element or nuclide or material and reaction plotted"""
    if isinstance(this, str):
        if type in openmc.data.DADZ:
            if this in ELEMENT_NAMES:
                return f'{this} {type}'
            else:  # this is a nuclide so the legend can contain more information
                z, a, m = openmc.data.zam(this)
                da, dz = openmc.data.DADZ[type]
                gnds_name = openmc.data.gnds_name(z + dz, a + da, m)
                # makes a string with nuclide reaction and new nuclide
                # For example "Be9 (n,2n) Be8"
                return f'{this} {type} {gnds_name}'
        return f'{this} {type}'
    elif this.name == '':
        return f'Material {this.id} {type}'
    else:
        return f'{this.name} {type}'


def _get_yaxis_label(reactions, divisor_types):
    """Gets a y axis label for the type of data plotted"""

    heat_values = {"heating", "heating-local", "damage-energy"}

    # if all the types are heating a different stem and unit is needed
    if all(set(value).issubset(heat_values) for value in reactions.values()):
        stem = "Heating"
    elif all(isinstance(item, str) for item in reactions.keys()):
        for nuc_reactions in reactions.values():
            for reaction in nuc_reactions:
                if reaction in heat_values:
                    raise TypeError(
                        "Mixture of heating and Microscopic reactions. "
                        "Invalid type for plotting"
                    )
        stem = "Microscopic"
    elif all(isinstance(item, openmc.Material) for item in reactions.keys()):
        stem = 'Macroscopic'
    else:
        msg = "Mixture of openmc.Material and elements/nuclides. Invalid type for plotting"
        raise TypeError(msg)

    if divisor_types:
        mid, units = "Data", ""
    else:
        mid = "Cross Section"
        units = {
            "Macroscopic": "[1/cm]",
            "Microscopic": "[b]",
            "Heating": "[eV-barn]",
        }[stem]

    return f'{stem} {mid} {units}'

def _get_title(reactions):
    """Gets a title for the type of data plotted"""
    if len(reactions) == 1:
        this, = reactions
        name = this.name if isinstance(this, openmc.Material) else this
        return f'Cross Section Plot For {name}'
    else:
        return 'Cross Section Plot'


def plot_xs(
    reactions: Dict[str | openmc.Material, List[str]],
    divisor_types: Iterable[str] | None = None,
    temperature: float = 294.0,
    axis: "plt.Axes" | None = None,
    sab_name: str | None = None,
    ce_cross_sections: str | None = None,
    mg_cross_sections: str | None = None,
    enrichment: float | None = None,
    plot_CE: bool = True,
    orders: Iterable[int] | None = None,
    divisor_orders: Iterable[int] | None = None,
    energy_axis_units: str = "eV",
    **kwargs,
) -> "plt.Figure" | None:
    """Creates a figure of continuous-energy cross sections for this item.

    Parameters
    ----------
    reactions : dict
        keys can be either a nuclide or element in string form or an
        openmc.Material object. Values are a list of the types of
        cross sections to include in the plot.
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
    sab_name : str, optional
        Name of S(a,b) library to apply to MT=2 data when applicable.
    ce_cross_sections : str, optional
        Location of cross_sections.xml file. Default is None.
    mg_cross_sections : str, optional
        Location of MGXS HDF5 Library file. Default is None.
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None.
    plot_CE : bool, optional
        Denotes whether or not continuous-energy will be plotted. Defaults to
        plotting the continuous-energy data.
    orders : Iterable of Integral, optional
        The scattering order or delayed group index to use for the
        corresponding entry in types. Defaults to the 0th order for scattering
        and the total delayed neutron data. This only applies to plots of
        multi-group data.
    divisor_orders : Iterable of Integral, optional
        Same as orders, but for divisor_types
    **kwargs :
        All keyword arguments are passed to
        :func:`matplotlib.pyplot.figure`.
    energy_axis_units : {'eV', 'keV', 'MeV'}
        Units used on the plot energy axis

        .. versionadded:: 0.15.0

    Returns
    -------
    fig : matplotlib.figure.Figure
        If axis is None, then a Matplotlib Figure of the generated
        cross section will be returned. Otherwise, a value of
        None will be returned as the figure and axes have already been
        generated.

    """
    import matplotlib.pyplot as plt

    cv.check_type("plot_CE", plot_CE, bool)
    cv.check_value("energy_axis_units", energy_axis_units, {"eV", "keV", "MeV"})

    axis_scaling_factor = {"eV": 1.0, "keV": 1e-3, "MeV": 1e-6}

    # Generate the plot
    if axis is None:
        fig, ax = plt.subplots(**kwargs)
    else:
        fig = None
        ax = axis

    all_types = []

    for this, types in reactions.items():
        all_types = all_types + types

        if plot_CE:
            cv.check_type("this", this, (str, openmc.Material))
            # Calculate for the CE cross sections
            E, data = calculate_cexs(this, types, temperature, sab_name,
                                    ce_cross_sections, enrichment)
            if divisor_types:
                cv.check_length('divisor types', divisor_types, len(types))
                Ediv, data_div = calculate_cexs(this, divisor_types, temperature,
                                                sab_name, ce_cross_sections,
                                                enrichment)

                # Create a new union grid, interpolate data and data_div on to that
                # grid, and then do the actual division
                Enum = E[:]
                E = np.union1d(Enum, Ediv)
                data_new = np.zeros((len(types), len(E)))

                for line in range(len(types)):
                    data_new[line, :] = \
                        np.divide(np.interp(E, Enum, data[line, :]),
                                np.interp(E, Ediv, data_div[line, :]))
                    if divisor_types[line] != 'unity':
                        types[line] = types[line] + ' / ' + divisor_types[line]
                data = data_new
        else:
            # Calculate for MG cross sections
            E, data = calculate_mgxs(this, types, orders, temperature,
                                    mg_cross_sections, ce_cross_sections,
                                    enrichment)

            if divisor_types:
                cv.check_length('divisor types', divisor_types, len(types))
                Ediv, data_div = calculate_mgxs(this, divisor_types,
                                                divisor_orders, temperature,
                                                mg_cross_sections,
                                                ce_cross_sections, enrichment)

                # Perform the division
                for line in range(len(types)):
                    data[line, :] /= data_div[line, :]
                    if divisor_types[line] != 'unity':
                        types[line] += ' / ' + divisor_types[line]

        E *= axis_scaling_factor[energy_axis_units]

        # Plot the data
        for i in range(len(data)):
            data[i, :] = np.nan_to_num(data[i, :])
            if np.sum(data[i, :]) > 0.:
                ax.plot(E, data[i, :], label=_get_legend_label(this, types[i]))

    # Set to loglog or semilogx depending on if we are plotting a data
    # type which we expect to vary linearly
    if set(all_types).issubset(PLOT_TYPES_LINEAR):
        ax.set_xscale('log')
        ax.set_yscale('linear')
    else:
        ax.set_xscale('log')
        ax.set_yscale('log')

    ax.set_xlabel(f"Energy [{energy_axis_units}]")
    if plot_CE:
        ax.set_xlim(
            _MIN_E * axis_scaling_factor[energy_axis_units],
            _MAX_E * axis_scaling_factor[energy_axis_units],
        )
    else:
        ax.set_xlim(E[-1], E[0])

    ax.set_ylabel(_get_yaxis_label(reactions, divisor_types))
    ax.legend(loc='best')
    ax.set_title(_get_title(reactions))

    return fig


def calculate_cexs(this, types, temperature=294., sab_name=None,
                   cross_sections=None, enrichment=None, ncrystal_cfg=None):
    """Calculates continuous-energy cross sections of a requested type.

    Parameters
    ----------
    this : str or openmc.Material
        Object to source data from. Nuclides and elements should be input as a
        str
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
    ncrystal_cfg : str, optional
        Configuration string for NCrystal material.

    Returns
    -------
    energy_grid : numpy.ndarray
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Cross sections calculated at the energy grid described by energy_grid

    """

    # Check types
    cv.check_type('this', this, (str, openmc.Material))
    cv.check_type('temperature', temperature, Real)
    if sab_name:
        cv.check_type('sab_name', sab_name, str)
    if enrichment:
        cv.check_type('enrichment', enrichment, Real)

    if isinstance(this, str):
        if this in ELEMENT_NAMES:
            energy_grid, data = _calculate_cexs_elem_mat(
                this, types, temperature, cross_sections, sab_name, enrichment
            )
        else:
            energy_grid, xs = _calculate_cexs_nuclide(
                this, types, temperature, sab_name, cross_sections,
                ncrystal_cfg
            )

            # Convert xs (Iterable of Callable) to a grid of cross section values
            # calculated on the points in energy_grid for consistency with the
            # element and material functions.
            data = np.zeros((len(types), len(energy_grid)))
            for line in range(len(types)):
                data[line, :] = xs[line](energy_grid)
    else:
        energy_grid, data = _calculate_cexs_elem_mat(this, types, temperature,
                                                     cross_sections)

    return energy_grid, data


def _calculate_cexs_nuclide(this, types, temperature=294., sab_name=None,
                            cross_sections=None, ncrystal_cfg=None):
    """Calculates continuous-energy cross sections of a requested type.

    Parameters
    ----------
    this : str
        Nuclide object to source data from
    types : Iterable of str or Integral
        The type of cross sections to calculate; values can either be those
        in openmc.PLOT_TYPES or keys from openmc.data.REACTION_MT which
        correspond to a reaction description e.g '(n,2n)' or integers which
        correspond to reaction channel (MT) numbers.
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    sab_name : str, optional
        Name of S(a,b) library to apply to MT=2 data when applicable.
    cross_sections : str, optional
        Location of cross_sections.xml file. Default is None.
    ncrystal_cfg : str, optional
        Configuration string for NCrystal material.

    Returns
    -------
    energy_grid : numpy.ndarray
        Energies at which cross sections are calculated, in units of eV
    data : Iterable of Callable
        Requested cross section functions

    """

    # Load the library
    library = openmc.data.DataLibrary.from_xml(cross_sections)

    # Convert temperature to format needed for access in the library
    strT = f"{int(round(temperature))}K"
    T = temperature

    # Now we can create the data sets to be plotted
    energy_grid = []
    xs = []
    lib = library.get_by_material(this)
    if lib is not None:
        nuc = openmc.data.IncidentNeutron.from_hdf5(lib['path'])
        # Obtain the nearest temperature
        if strT in nuc.temperatures:
            nucT = strT
        else:
            delta_T = np.array(nuc.kTs) - T * openmc.data.K_BOLTZMANN
            closest_index = np.argmin(np.abs(delta_T))
            nucT = nuc.temperatures[closest_index]

        # Prep S(a,b) data if needed
        if sab_name:
            sab = openmc.data.ThermalScattering.from_hdf5(
                library.get_by_material(sab_name, data_type='thermal')['path'])
            # Obtain the nearest temperature
            if strT in sab.temperatures:
                sabT = strT
            else:
                delta_T = np.array(sab.kTs) - T * openmc.data.K_BOLTZMANN
                closest_index = np.argmin(np.abs(delta_T))
                sabT = sab.temperatures[closest_index]

            # Create an energy grid composed the S(a,b) and the nuclide's grid
            grid = nuc.energy[nucT]
            sab_Emax = 0.
            sab_funcs = []
            if sab.elastic is not None:
                elastic = sab.elastic.xs[sabT]
                if isinstance(elastic, openmc.data.CoherentElastic):
                    grid = np.union1d(grid, elastic.bragg_edges)
                    if elastic.bragg_edges[-1] > sab_Emax:
                        sab_Emax = elastic.bragg_edges[-1]
                elif isinstance(elastic, openmc.data.Tabulated1D):
                    grid = np.union1d(grid, elastic.x)
                    if elastic.x[-1] > sab_Emax:
                        sab_Emax = elastic.x[-1]
                sab_funcs.append(elastic)
            if sab.inelastic is not None:
                inelastic = sab.inelastic.xs[sabT]
                grid = np.union1d(grid, inelastic.x)
                if inelastic.x[-1] > sab_Emax:
                        sab_Emax = inelastic.x[-1]
                sab_funcs.append(inelastic)
            energy_grid = grid
        else:
            energy_grid = nuc.energy[nucT]

        # Parse the types
        mts = []
        ops = []
        yields = []
        for line in types:
            if line in PLOT_TYPES:
                tmp_mts = [mtj for mti in PLOT_TYPES_MT[line] for mtj in
                           nuc.get_reaction_components(mti)]
                mts.append(tmp_mts)
                if line.startswith('nu'):
                    yields.append(True)
                else:
                    yields.append(False)
                if XI_MT in tmp_mts:
                    ops.append((np.add,) * (len(tmp_mts) - 2) + (np.multiply,))
                else:
                    ops.append((np.add,) * (len(tmp_mts) - 1))
            elif line in openmc.data.REACTION_MT:
                mt_number = openmc.data.REACTION_MT[line]
                cv.check_type('MT in types', mt_number, Integral)
                cv.check_greater_than('MT in types', mt_number, 0)
                tmp_mts = nuc.get_reaction_components(mt_number)
                mts.append(tmp_mts)
                ops.append((np.add,) * (len(tmp_mts) - 1))
                yields.append(False)
            elif isinstance(line, int):
                # Not a built-in type, we have to parse it ourselves
                cv.check_type('MT in types', line, Integral)
                cv.check_greater_than('MT in types', line, 0)
                tmp_mts = nuc.get_reaction_components(line)
                mts.append(tmp_mts)
                ops.append((np.add,) * (len(tmp_mts) - 1))
                yields.append(False)
            else:
                raise TypeError("Invalid type", line)

        for i, mt_set in enumerate(mts):
            # Get the reaction xs data from the nuclide
            funcs = []
            op = ops[i]
            for mt in mt_set:
                if mt == 2:
                    if sab_name:
                        # Then we need to do a piece-wise function of
                        # The S(a,b) and non-thermal data
                        sab_sum = openmc.data.Sum(sab_funcs)
                        pw_funcs = openmc.data.Regions1D(
                            [sab_sum, nuc[mt].xs[nucT]],
                            [sab_Emax])
                        funcs.append(pw_funcs)
                    elif ncrystal_cfg:
                        import NCrystal
                        nc_scatter = NCrystal.createScatter(ncrystal_cfg)
                        nc_func = nc_scatter.crossSectionNonOriented
                        nc_emax = 5 # eV # this should be obtained from NCRYSTAL_MAX_ENERGY
                        energy_grid = np.union1d(np.geomspace(min(energy_grid),
                                                              1.1*nc_emax,
                                                              1000),energy_grid) # NCrystal does not have
                                                                                 # an intrinsic energy grid
                        pw_funcs = openmc.data.Regions1D(
                            [nc_func, nuc[mt].xs[nucT]],
                            [nc_emax])
                        funcs.append(pw_funcs)
                    else:
                        funcs.append(nuc[mt].xs[nucT])
                elif mt in nuc:
                    if yields[i]:
                        # Get the total yield first if available. This will be
                        # used primarily for fission.
                        for prod in chain(nuc[mt].products,
                                          nuc[mt].derived_products):
                            if prod.particle == 'neutron' and \
                                prod.emission_mode == 'total':
                                func = openmc.data.Combination(
                                    [nuc[mt].xs[nucT], prod.yield_],
                                    [np.multiply])
                                funcs.append(func)
                                break
                        else:
                            # Total doesn't exist so we have to create from
                            # prompt and delayed. This is used for scatter
                            # multiplication.
                            func = None
                            for prod in chain(nuc[mt].products,
                                              nuc[mt].derived_products):
                                if prod.particle == 'neutron' and \
                                    prod.emission_mode != 'total':
                                    if func:
                                        func = openmc.data.Combination(
                                            [prod.yield_, func], [np.add])
                                    else:
                                        func = prod.yield_
                            if func:
                                funcs.append(openmc.data.Combination(
                                    [func, nuc[mt].xs[nucT]], [np.multiply]))
                            else:
                                # If func is still None, then there were no
                                # products. In that case, assume the yield is
                                # one as its not provided for some summed
                                # reactions like MT=4
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
            funcs = funcs if funcs else [lambda x: 0.]
            xs.append(openmc.data.Combination(funcs, op))
    else:
        raise ValueError(this + " not in library")

    return energy_grid, xs


def _calculate_cexs_elem_mat(this, types, temperature=294.,
                             cross_sections=None, sab_name=None,
                             enrichment=None):
    """Calculates continuous-energy cross sections of a requested type.

    Parameters
    ----------
    this : openmc.Material or str
        Object to source data from. Element can be input as str
    types : Iterable of values of PLOT_TYPES
        The type of cross sections to calculate
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    cross_sections : str, optional
        Location of cross_sections.xml file. Default is None.
    sab_name : str, optional
        Name of S(a,b) library to apply to MT=2 data when applicable.
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None
        (natural composition).

    Returns
    -------
    energy_grid : numpy.ndarray
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Cross sections calculated at the energy grid described by energy_grid

    """

    if isinstance(this, openmc.Material):
        if this.temperature is not None:
            T = this.temperature
        else:
            T = temperature
    else:
        T = temperature

    # Load the library
    library = openmc.data.DataLibrary.from_xml(cross_sections)

    ncrystal_cfg = None
    if isinstance(this, openmc.Material):
        # Expand elements in to nuclides with atomic densities
        nuc_fractions = this.get_nuclide_atom_densities()
        # Create a dict of [nuclide name] = nuclide object to carry forward
        # with a common nuclides format between openmc.Material and Elements
        nuclides = {nuclide: nuclide for nuclide in nuc_fractions}
        # Add NCrystal cfg string if it exists
        ncrystal_cfg = this.ncrystal_cfg
    else:
        # Expand elements in to nuclides with atomic densities
        nuclides = openmc.Element(this).expand(1., 'ao', enrichment=enrichment,
                               cross_sections=cross_sections)
        # For ease of processing split out the nuclide and its fraction
        nuc_fractions = {nuclide[0]: nuclide[1] for nuclide in nuclides}
        # Create a dict of [nuclide name] = nuclide object to carry forward
        # with a common nuclides format between openmc.Material and Elements
        nuclides = {nuclide[0]: nuclide[0] for nuclide in nuclides}

    # Identify the nuclides which have S(a,b) data
    sabs = {}
    for nuclide in nuclides.items():
        sabs[nuclide[0]] = None
    if isinstance(this, openmc.Material):
        for sab_name, _ in this._sab:
            sab = openmc.data.ThermalScattering.from_hdf5(
                library.get_by_material(sab_name, data_type='thermal')['path'])
            for nuc in sab.nuclides:
                sabs[nuc] = sab_name
    else:
        if sab_name:
            sab = openmc.data.ThermalScattering.from_hdf5(
                library.get_by_material(sab_name, data_type='thermal')['path'])
            for nuc in sab.nuclides:
                sabs[nuc] = sab_name

    # Now we can create the data sets to be plotted
    xs = {}
    E = []
    for nuclide in nuclides.items():
        name = nuclide[0]
        nuc = nuclide[1]
        sab_name = sabs[name]
        temp_E, temp_xs = calculate_cexs(nuc, types, T, sab_name, cross_sections,
                                         ncrystal_cfg=ncrystal_cfg
                                         )
        E.append(temp_E)
        # Since the energy grids are different, store the cross sections as
        # a tabulated function so they can be calculated on any grid needed.
        xs[name] = [openmc.data.Tabulated1D(temp_E, temp_xs[line])
                    for line in range(len(types))]

    # Condense the data for every nuclide
    # First create a union energy grid
    energy_grid = E[0]
    for grid in E[1:]:
        energy_grid = np.union1d(energy_grid, grid)

    # Now we can combine all the nuclidic data
    data = np.zeros((len(types), len(energy_grid)))
    for line in range(len(types)):
        if types[line] == 'unity':
            data[line, :] = 1.
        else:
            for nuclide in nuclides.items():
                name = nuclide[0]
                data[line, :] += (nuc_fractions[name] *
                                  xs[name][line](energy_grid))

    return energy_grid, data


def calculate_mgxs(this, types, orders=None, temperature=294.,
                   cross_sections=None, ce_cross_sections=None,
                   enrichment=None):
    """Calculates multi-group cross sections of a requested type.

    If the data for the nuclide or macroscopic object in the library is
    represented as angle-dependent data then this method will return the
    geometric average cross section over all angles.

    Parameters
    ----------
    this : str or openmc.Material
        Object to source data from. Nuclides and elements can be input as a str
    types : Iterable of values of PLOT_TYPES_MGXS
        The type of cross sections to calculate
    orders : Iterable of Integral, optional
        The scattering order or delayed group index to use for the
        corresponding entry in types. Defaults to the 0th order for scattering
        and the total delayed neutron data.
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.
    cross_sections : str, optional
        Location of MGXS HDF5 Library file. Default is None.
    ce_cross_sections : str, optional
        Location of continuous-energy cross_sections.xml file. Default is None.
    enrichment : float, optional
        Enrichment for U235 in weight percent. For example, input 4.95 for
        4.95 weight percent enriched U. Default is None
        (natural composition).

    Returns
    -------
    energy_grid : numpy.ndarray
        Energies at which cross sections are calculated, in units of eV
    data : numpy.ndarray
        Cross sections calculated at the energy grid described by energy_grid

    """

    # Check types
    cv.check_type('temperature', temperature, Real)
    if enrichment:
        cv.check_type('enrichment', enrichment, Real)
    cv.check_iterable_type('types', types, str)

    cv.check_type("cross_sections", cross_sections, str)
    library = openmc.MGXSLibrary.from_hdf5(cross_sections)

    if this in ELEMENT_NAMES or isinstance(this, openmc.Material):
        mgxs = _calculate_mgxs_elem_mat(this, types, library, orders,
                                        temperature, ce_cross_sections,
                                        enrichment)
    elif isinstance(this, str):
        mgxs = _calculate_mgxs_nuc_macro(this, types, library, orders,
                                         temperature)
    else:
        raise TypeError("Invalid type")

    # Convert the data to the format needed
    data = np.zeros((len(types), 2 * library.energy_groups.num_groups))
    energy_grid = np.zeros(2 * library.energy_groups.num_groups)
    for g in range(library.energy_groups.num_groups):
        energy_grid[g * 2: g * 2 + 2] = \
            library.energy_groups.group_edges[g: g + 2]
    # Ensure the energy will show on a log-axis by replacing 0s with a
    # sufficiently small number
    energy_grid[0] = max(energy_grid[0], _MIN_E)

    for line in range(len(types)):
        for g in range(library.energy_groups.num_groups):
            data[line, g * 2: g * 2 + 2] = mgxs[line, g]

    return energy_grid[::-1], data


def _calculate_mgxs_nuc_macro(this, types, library, orders=None,
                              temperature=294.):
    """Determines the multi-group cross sections of a nuclide or macroscopic
    object.

    If the data for the nuclide or macroscopic object in the library is
    represented as angle-dependent data then this method will return the
    geometric average cross section over all angles.

    Parameters
    ----------
    this : str
        Object to source data from
    types : Iterable of str
        The type of cross sections to calculate; values can either be those
        in openmc.PLOT_TYPES_MGXS
    library : openmc.MGXSLibrary
        MGXS Library containing the data of interest
    orders : Iterable of Integral, optional
        The scattering order or delayed group index to use for the
        corresponding entry in types. Defaults to the 0th order for scattering
        and the total delayed neutron data.
    temperature : float, optional
        Temperature in Kelvin to plot. If not specified, a default
        temperature of 294K will be plotted. Note that the nearest
        temperature in the library for each nuclide will be used as opposed
        to using any interpolation.

    Returns
    -------
    data : numpy.ndarray
        Cross sections calculated at the energy grid described by energy_grid

    """

    # Check the parameters and grab order/delayed groups
    if orders:
        cv.check_iterable_type('orders', orders, Integral,
                               min_depth=len(types), max_depth=len(types))
    else:
        orders = [None] * len(types)
    for i, line in enumerate(types):
        cv.check_type("line", line, str)
        cv.check_value("line", line, PLOT_TYPES_MGXS)
        if orders[i]:
            cv.check_greater_than("order value", orders[i], 0, equality=True)

    xsdata = library.get_by_name(this)

    if xsdata is not None:
        # Obtain the nearest temperature
        t = np.abs(xsdata.temperatures - temperature).argmin()

        # Get the data
        data = np.zeros((len(types), library.energy_groups.num_groups))
        for i, line in enumerate(types):
            if 'fission' in line and not xsdata.fissionable:
                continue
            elif line == 'unity':
                data[i, :] = 1.
            else:
                # Now we have to get the cross section data and properly
                # treat it depending on the requested type.
                # First get the data in a generic fashion
                temp_data = getattr(xsdata, _PLOT_MGXS_ATTR[line])[t]
                shape = temp_data.shape[:]
                # If we have angular data, then want the geometric
                # average over all provided angles.  Since the angles are
                # equi-distant, un-weighted averaging will suffice
                if xsdata.representation == 'angle':
                    temp_data = np.mean(temp_data, axis=(0, 1))

                # Now we can look at the shape of the data to identify how
                # it should be modified to produce an array of values
                # with groups.
                if shape in (xsdata.xs_shapes["[G']"],
                             xsdata.xs_shapes["[G]"]):
                    # Then the data is already an array vs groups so copy
                    # and move along
                    data[i, :] = temp_data
                elif shape == xsdata.xs_shapes["[G][G']"]:
                    # Sum the data over outgoing groups to create our array vs
                    # groups
                    data[i, :] = np.sum(temp_data, axis=1)
                elif shape == xsdata.xs_shapes["[DG]"]:
                    # Then we have a constant vs groups with a value for each
                    # delayed group. The user-provided value of orders tells us
                    # which delayed group we want. If none are provided, then
                    # we sum all the delayed groups together.
                    if orders[i]:
                        if orders[i] < len(shape[0]):
                            data[i, :] = temp_data[orders[i]]
                    else:
                        data[i, :] = np.sum(temp_data[:])
                elif shape in (xsdata.xs_shapes["[DG][G']"],
                               xsdata.xs_shapes["[DG][G]"]):
                    # Then we have an array vs groups with values for each
                    # delayed group. The user-provided value of orders tells us
                    # which delayed group we want. If none are provided, then
                    # we sum all the delayed groups together.
                    if orders[i]:
                        if orders[i] < len(shape[0]):
                            data[i, :] = temp_data[orders[i], :]
                    else:
                        data[i, :] = np.sum(temp_data[:, :], axis=0)
                elif shape == xsdata.xs_shapes["[DG][G][G']"]:
                    # Then we have a delayed group matrix. We will first
                    # remove the outgoing group dependency
                    temp_data = np.sum(temp_data, axis=-1)
                    # And then proceed in exactly the same manner as the
                    # "[DG][G']" or "[DG][G]" shapes in the previous block.
                    if orders[i]:
                        if orders[i] < len(shape[0]):
                            data[i, :] = temp_data[orders[i], :]
                    else:
                        data[i, :] = np.sum(temp_data[:, :], axis=0)
                elif shape == xsdata.xs_shapes["[G][G'][Order]"]:
                    # This is a scattering matrix with angular data
                    # First remove the outgoing group dependence
                    temp_data = np.sum(temp_data, axis=1)
                    # The user either provided a specific order or we resort
                    # to the default 0th order
                    if orders[i]:
                        order = orders[i]
                    else:
                        order = 0
                    # If the order is available, store the data for that order
                    # if it is not available, then the expansion coefficient
                    # is zero and thus we already have the correct value.
                    if order < shape[1]:
                        data[i, :] = temp_data[:, order]
    else:
        raise ValueError(f"{this} not present in provided MGXS library")

    return data


def _calculate_mgxs_elem_mat(this, types, library, orders=None,
                             temperature=294., ce_cross_sections=None,
                             enrichment=None):
    """Determines the multi-group cross sections of an element or material
    object.

    If the data for the nuclide or macroscopic object in the library is
    represented as angle-dependent data then this method will return the
    geometric average cross section over all angles.

    Parameters
    ----------
    this : str or openmc.Material
        Object to source data from. Elements can be input as a str
    types : Iterable of str
        The type of cross sections to calculate; values can either be those
        in openmc.PLOT_TYPES_MGXS
    library : openmc.MGXSLibrary
        MGXS Library containing the data of interest
    orders : Iterable of Integral, optional
        The scattering order or delayed group index to use for the
        corresponding entry in types. Defaults to the 0th order for scattering
        and the total delayed neutron data.
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
        Cross sections calculated at the energy grid described by energy_grid

    """

    if isinstance(this, openmc.Material):
        if this.temperature is not None:
            T = this.temperature
        else:
            T = temperature

        # Check to see if we have nuclides/elements or a macroscopic object
        if this._macroscopic is not None:
            # We have macroscopics
            nuclides = {this._macroscopic: this.density}
        else:
            # Expand elements in to nuclides with atomic densities
            nuclides = this.get_nuclide_atom_densities()

        # For ease of processing split out nuc and nuc_density
        nuc_fraction = list(nuclides.values())
    else:
        T = temperature
        # Expand elements in to nuclides with atomic densities
        nuclides = openmc.Element(this).expand(100., 'ao', enrichment=enrichment,
                               cross_sections=ce_cross_sections)

        # For ease of processing split out nuc and nuc_fractions
        nuc_fraction = [nuclide[1] for nuclide in nuclides]

    nuc_data = []
    for nuclide in nuclides.items():
        nuc_data.append(_calculate_mgxs_nuc_macro(nuclide[0], types, library,
                                                  orders, T))

    # Combine across the nuclides
    data = np.zeros((len(types), library.energy_groups.num_groups))
    for line in range(len(types)):
        if types[line] == 'unity':
            data[line, :] = 1.
        else:
            for n in range(len(nuclides)):
                data[line, :] += nuc_fraction[n] * nuc_data[n][line, :]

    return data
