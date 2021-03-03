import openmc
import numpy as np

names = ['H', 'O', 'Zr', 'U235', 'U238']


def build_openmc_xs_lib(name, groups, temperatures, xsdict, micro=True):
    """Build an Openm XSdata based on dictionary values"""
    xsdata = openmc.XSdata(name, groups, temperatures=temperatures)
    xsdata.order = 0
    for tt in temperatures:
        xsdata.set_absorption(xsdict[tt]['absorption'][name], temperature=tt)
        xsdata.set_scatter_matrix(xsdict[tt]['scatter'][name], temperature=tt)
        xsdata.set_total(xsdict[tt]['total'][name], temperature=tt)
        if (name in xsdict[tt]['nu-fission'].keys()):
            xsdata.set_nu_fission(xsdict[tt]['nu-fission'][name],
                                  temperature=tt)
            xsdata.set_chi(np.array([1., 0.]), temperature=tt)
    return xsdata


def create_micro_xs_dict():
    """Returns micro xs library"""
    xs_micro = {}
    reactions = ['absorption', 'total', 'scatter', 'nu-fission']
    # chi is unnecessary  when energy bound is in thermal region
    # Temperature 300K
    # absorption
    xs_micro[300] = {r: {} for r in reactions}
    xs_micro[300]['absorption']['H'] = np.array([1.0285E-4, 0.0057])
    xs_micro[300]['absorption']['O'] = np.array([7.1654E-5, 3.0283E-6])
    xs_micro[300]['absorption']['Zr'] = np.array([4.5918E-5, 3.6303E-5])
    xs_micro[300]['absorption']['U235'] = np.array([0.0035, 0.1040])
    xs_micro[300]['absorption']['U238'] = np.array([0.0056, 0.0094])
    # nu-scatter matrix
    xs_micro[300]['scatter']['H'] = np.array([[[0.0910, 0.01469],
                                               [0.0, 0.3316]]])
    xs_micro[300]['scatter']['O'] = np.array([[[0.0814, 3.3235E-4],
                                               [0.0, 0.0960]]])
    xs_micro[300]['scatter']['Zr'] = np.array([[[0.0311, 2.6373E-5],
                                                [0.0, 0.0315]]])
    xs_micro[300]['scatter']['U235'] = np.array([[[0.0311, 2.6373E-5],
                                                  [0.0, 0.0315]]])
    xs_micro[300]['scatter']['U238'] = np.array([[[0.0551, 2.2341E-5],
                                                  [0.0, 0.0526]]])
    # nu-fission
    xs_micro[300]['nu-fission']['U235'] = np.array([0.0059, 0.2160])
    xs_micro[300]['nu-fission']['U238'] = np.array([0.0019, 1.4627E-7])
    # total
    xs_micro[300]['total']['H'] = xs_micro[300]['absorption']['H'] + \
                                  np.sum(xs_micro[300]['scatter']['H'][0], 1)
    xs_micro[300]['total']['O'] = xs_micro[300]['absorption']['O'] + \
                                  np.sum(xs_micro[300]['scatter']['O'][0], 1)

    xs_micro[300]['total']['Zr'] = xs_micro[300]['absorption']['Zr'] + \
                                   np.sum(xs_micro[300]['scatter']['Zr'][0], 1)

    xs_micro[300]['total']['U235'] = xs_micro[300]['absorption']['U235'] + \
                                     np.sum(xs_micro[300]['scatter']['U235'][0], 1)

    xs_micro[300]['total']['U238'] = xs_micro[300]['absorption']['U238'] + \
                                     np.sum(xs_micro[300]['scatter']['U238'][0], 1)

    # Temperature 600K
    xs_micro[600] = {r: {} for r in reactions}
    # absorption
    xs_micro[600]['absorption']['H'] = np.array([1.0356E-4, 0.0046])
    xs_micro[600]['absorption']['O'] = np.array([7.2678E-5, 2.4963E-6])
    xs_micro[600]['absorption']['Zr'] = np.array([4.7256E-5, 2.9757E-5])
    xs_micro[600]['absorption']['U235'] = np.array([0.0035, 0.0853])
    xs_micro[600]['absorption']['U238'] = np.array([0.0058, 0.0079])
    # nu-scatter matrix
    xs_micro[600]['scatter']['H'] = np.array([[[0.0910, 0.0138],
                                               [0.0, 0.3316]]])
    xs_micro[600]['scatter']['O'] = np.array([[[0.0814, 3.5367E-4],
                                               [0.0, 0.0959]]])
    xs_micro[600]['scatter']['Zr'] = np.array([[[0.0311, 3.2293E-5],
                                                [0.0, 0.0314]]])
    xs_micro[600]['scatter']['U235'] = np.array([[[0.0022, 1.9763E-6],
                                                  [9.1634E-8, 0.0039]]])
    xs_micro[600]['scatter']['U238'] = np.array([[[0.0556, 2.8803E-5],
                                                  [0.0, 0.0536]]])
    # nu-fission
    xs_micro[600]['nu-fission']['U235'] = np.array([0.0059, 0.1767])
    xs_micro[600]['nu-fission']['U238'] = np.array([0.0019, 1.2405E-7])
    # total
    xs_micro[600]['total']['H'] = xs_micro[600]['absorption']['H'] + \
                                  np.sum(xs_micro[600]['scatter']['H'][0], 1)
    xs_micro[600]['total']['O'] = xs_micro[600]['absorption']['O'] + \
                                  np.sum(xs_micro[600]['scatter']['O'][0], 1)

    xs_micro[600]['total']['Zr'] = xs_micro[600]['absorption']['Zr'] + \
                                   np.sum(xs_micro[600]['scatter']['Zr'][0], 1)

    xs_micro[600]['total']['U235'] = xs_micro[600]['absorption']['U235'] + \
                                     np.sum(xs_micro[600]['scatter']['U235'][0], 1)

    xs_micro[600]['total']['U238'] = xs_micro[600]['absorption']['U238'] + \
                                     np.sum(xs_micro[600]['scatter']['U238'][0], 1)

    # Temperature 900K
    xs_micro[900] = {r: {} for r in reactions}
    # absorption
    xs_micro[900]['absorption']['H'] = np.array([1.0529E-4, 0.0040])
    xs_micro[900]['absorption']['O'] = np.array([7.3055E-5, 2.1850E-6])
    xs_micro[900]['absorption']['Zr'] = np.array([4.7141E-5, 2.5941E-5])
    xs_micro[900]['absorption']['U235'] = np.array([0.0035, 0.0749])
    xs_micro[900]['absorption']['U238'] = np.array([0.0060, 0.0071])
    # total
    xs_micro[900]['total']['H'] = np.array([0.2982, 0.7332])
    xs_micro[900]['total']['O'] = np.array([0.0885, 0.1004])
    xs_micro[900]['total']['Zr'] = np.array([0.0370, 0.0317])
    xs_micro[900]['total']['U235'] = np.array([0.0061, 0.0789])
    xs_micro[900]['total']['U238'] = np.array([0.0707, 0.0613])
    # nu-scatter matrix
    xs_micro[900]['scatter']['H'] = np.array([[[0.0913, 0.0147],
                                               [0.0, 0.4020]]])
    xs_micro[900]['scatter']['O'] = np.array([[[0.0812, 4.0413E-4],
                                               [0.0, 0.0965]]])
    xs_micro[900]['scatter']['Zr'] = np.array([[[0.0311, 3.6735E-5],
                                                [0.0, 0.0314]]])
    xs_micro[900]['scatter']['U235'] = np.array([[[0.0022, 2.9034E-6],
                                                  [1.3117E-8, 0.0039]]])
    xs_micro[900]['scatter']['U238'] = np.array([[[0.0560, 3.7619E-5],
                                                  [0.0, 0.0538]]])
    # nu-fission
    xs_micro[900]['nu-fission']['U235'] = np.array([0.0059, 0.1545])
    xs_micro[900]['nu-fission']['U238'] = np.array([0.0019, 1.1017E-7])
    # total
    xs_micro[900]['total']['H'] = xs_micro[900]['absorption']['H'] + \
                                  np.sum(xs_micro[900]['scatter']['H'][0], 1)
    xs_micro[900]['total']['O'] = xs_micro[900]['absorption']['O'] + \
                                  np.sum(xs_micro[900]['scatter']['O'][0], 1)

    xs_micro[900]['total']['Zr'] = xs_micro[900]['absorption']['Zr'] + \
                                   np.sum(xs_micro[900]['scatter']['Zr'][0], 1)

    xs_micro[900]['total']['U235'] = xs_micro[900]['absorption']['U235'] + \
                                     np.sum(xs_micro[900]['scatter']['U235'][0], 1)

    xs_micro[900]['total']['U238'] = xs_micro[900]['absorption']['U238'] + \
                                     np.sum(xs_micro[900]['scatter']['U238'][0], 1)

    # roll axis for scatter matrix
    for t in xs_micro:
        for n in xs_micro[t]['scatter']:
            xs_micro[t]['scatter'][n] = np.rollaxis(xs_micro[t]['scatter'][n],
                                                    0, 3)
    return xs_micro


def create_macro_dict(xs_micro):
    """Create a dictionary with two group cross-section"""
    xs_macro = {}
    for t, d1 in xs_micro.items():
        xs_macro[t] = {}
        for r, d2 in d1.items():
            temp = []
            xs_macro[t][r] = {}
            for n, v in d2.items():
                temp.append(d2[n])
            # The name 'macro' is needed to store data at the same level
            # of a xs_macro dictionary as for xs_micro and use it in
            # function build_openmc_xs_lib
            xs_macro[t][r]['macro'] = sum(temp)
    return xs_macro


def create_openmc_2mg_libs(names):
    """Built a micro/macro two group openmc MGXS libraries"""
    # Initialized library params
    group_edges = [0.0, 0.625, 20.0e6]
    groups = openmc.mgxs.EnergyGroups(group_edges=group_edges)
    mg_cross_sections_file_micro = openmc.MGXSLibrary(groups)
    mg_cross_sections_file_macro = openmc.MGXSLibrary(groups)
    # Building a micro mg library
    micro_cs = create_micro_xs_dict()
    for name in names:
        mg_cross_sections_file_micro.add_xsdata(build_openmc_xs_lib(name,
                                                                    groups,
                                                                    [t for t in
                                                                     micro_cs],
                                                                    micro_cs))
    # Building a macro mg library
    macro_xs = create_macro_dict(micro_cs)
    mg_cross_sections_file_macro.add_xsdata(build_openmc_xs_lib('macro',
                                                                groups,
                                                                [t for t in
                                                                 macro_xs],
                                                                macro_xs))
    # Exporting library to hdf5 files
    mg_cross_sections_file_micro.export_to_hdf5('micro_2g.h5')
    mg_cross_sections_file_macro.export_to_hdf5('macro_2g.h5')
    # Returning the macro_xs dict is needed for analytical solution
    return macro_xs


def analytical_solution_2g_therm(xsmin, xsmax=None, wgt=1.0):
    """ Calculate eigenvalue based on analytical solution for eq Lf = (1/k)Qf
    in two group for infinity dilution media in assumption of group
    boundary in thermal spectra < 1.e+3 Ev
    Parameters:
    ----------
    xsmin : dict
        macro cross-sections dictionary with minimum range temperature
    xsmax : dict
        macro cross-sections dictionary with maximum range temperature
        by default: None not used for standalone temperature
    wgt : float
        weight for interpolation by default 1.0
    Returns:
    -------
    keff : np.float64
       analytical eigenvalue of critical eq matrix
    """
    if xsmax is None:
        sa = xsmin['absorption']['macro']
        ss12 = xsmin['scatter']['macro'][0][1][0]
        nsf = xsmin['nu-fission']['macro']
    else:
        sa = xsmin['absorption']['macro'] * wgt + \
             xsmax['absorption']['macro'] * (1 - wgt)
        ss12 = xsmin['scatter']['macro'][0][1][0] * wgt + \
               xsmax['scatter']['macro'][0][1][0] * (1 - wgt)
        nsf = xsmin['nu-fission']['macro'] * wgt + \
              xsmax['nu-fission']['macro'] * (1 - wgt)
    L = np.array([sa[0] + ss12, 0.0, -ss12, sa[1]]).reshape(2, 2)
    Q = np.array([nsf[0], nsf[1], 0.0, 0.0]).reshape(2, 2)
    arr = np.linalg.inv(L).dot(Q)
    return np.amax(np.linalg.eigvals(arr))


def build_inf_model(xsnames, xslibname, temperature, tempmethod='nearest'):
    """ Building an infinite medium for openmc multi-group testing
    Parameters:
    ----------
    xsnames : list of str()
        list with xs names
    xslibname:
        name of hdf5 file with cross-section library
    temperature : float
        value of a current temperature in K
    tempmethod : {'nearest', 'interpolation'}
        by default 'nearest'
    """
    inf_medium = openmc.Material(name='test material', material_id=1)
    inf_medium.set_density("sum")
    for xs in xsnames:
        inf_medium.add_nuclide(xs, 1)
    INF = 11.1
    # Instantiate a Materials collection and export to XML
    materials_file = openmc.Materials([inf_medium])
    materials_file.cross_sections = xslibname
    materials_file.export_to_xml()

    # Instantiate boundary Planes
    min_x = openmc.XPlane(boundary_type='reflective', x0=-INF)
    max_x = openmc.XPlane(boundary_type='reflective', x0=INF)
    min_y = openmc.YPlane(boundary_type='reflective', y0=-INF)
    max_y = openmc.YPlane(boundary_type='reflective', y0=INF)

    # Instantiate a Cell
    cell = openmc.Cell(cell_id=1, name='cell')
    cell.temperature = temperature
    # Register bounding Surfaces with the Cell
    cell.region = +min_x & -max_x & +min_y & -max_y

    # Fill the Cell with the Material
    cell.fill = inf_medium

    # Create root universe
    root_universe = openmc.Universe(name='root universe', cells=[cell])

    # Create Geometry and set root Universe
    openmc_geometry = openmc.Geometry(root_universe)

    # Export to "geometry.xml"
    openmc_geometry.export_to_xml()

    # OpenMC simulation parameters
    batches = 200
    inactive = 5
    particles = 5000

    # Instantiate a Settings object
    settings_file = openmc.Settings()
    settings_file.batches = batches
    settings_file.inactive = inactive
    settings_file.particles = particles
    settings_file.energy_mode = 'multi-group'
    settings_file.output = {'summary': False}
    # Create an initial uniform spatial source distribution over fissionable zones
    bounds = [-INF, -INF, -INF, INF, INF, INF]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    settings_file.temperature = {'method': tempmethod}
    settings_file.source = openmc.Source(space=uniform_dist)
    settings_file.export_to_xml()
