import numpy as np

import openmc.checkvalue as cv
import openmc.data
import h5py

import matplotlib
from mpl_toolkits.mplot3d import Axes3D
matplotlib.use("Agg")
import matplotlib.pyplot as plt
#import seaborn as sns

SPATIAL_PLOT_TYPES = ['pin_powers', 'assembly_powers', 'flux',
                      'adjoint_flux', 'precursors', 'kappa_fission',
                      'pin_cell_kappa_fission', 'assembly_kappa_fission',
                      'pin_cell_shape', 'assembly_shape']
SCALAR_PLOT_TYPES = ['core_power_density', 'reactivity', 'pnl', 'beta_eff']


def get_scalar_data(variable, log_file, variable_twin=None, directory='.'):
    """Creates a figure of a scalar variable from a transient solve.

    Parameters
    ----------
    variable : str
        Variable to plot.
    log_file : str
        Location of log file.
    variable_twin : str
        Variable to plot on twinx axis. Default is None.
    directory : str
        Directory to save pngs. Default is '.'.
    axis : matplotlib.axes, optional
        A previously generated axis to use for plotting. If not specified,
        a new axis and figure will be generated.
    **kwargs
        All keyword arguments are passed to
        :func:`matplotlib.pyplot.figure`.

    Returns
    -------
    fig : matplotlib.figure.Figure or None
        If axis is None, then a Matplotlib Figure of the generated spatial
        variable will be returned. Otherwise, a value of None will be returned
        as the figure and axes have already been generated.

    """

    # Open hdf5 log file
    f = h5py.File(log_file, 'r')

    # Retrieve time steps
    time_steps = f['time_steps'].keys()
    time_steps = np.sort(np.array([float(i) for i in time_steps]))

    # Retrieve the spatial data
    scalar_data = np.zeros((len(time_steps), 3))
    for i,step in enumerate(time_steps):
        data = f['time_steps'][str(step)].attrs[variable]
        scalar_data[i,0] = step
        scalar_data[i,1] = data

        if variable_twin is not None:
            data = f['time_steps'][str(step)].attrs[variable_twin]
            scalar_data[i,2] = data

    f.close()
    return scalar_data


def spatial_plot(variable, log_file, directory='.', plane='xy', plane_num=0,
                 group=0, axis=None, animation=False, bar=True, normalize=550., **kwargs):
    """Creates a figure of a spatial variable from a transient solve.

    Parameters
    ----------
    variable : str
        Variable to plot.
    log_file : str
        Location of log file.
    directory : str
        Directory to save pngs. Default is '.'.
    plane : str
        Cut plane to plot variable. Default is xy.
    plane_num : int
        Plane number to plot variable. Default is 0.
    group : int
        Energy group to plot variable. Default is 0.
    axis : matplotlib.axes, optional
        A previously generated axis to use for plotting. If not specified,
        a new axis and figure will be generated.
    animation : bool
        Whether to create an animation or create separate plots for each time
        step.
    **kwargs
        All keyword arguments are passed to
        :func:`matplotlib.pyplot.figure`.

    Returns
    -------
    fig : matplotlib.figure.Figure or None
        If animation is False, the plots for each time saved and None
        will be returned. If animation is True and axis is None, then a
        Matplotlib Figure of the generated spatial variable will be returned.
        Otherwise, a value of None will be returned as the figure and axes have
        already been generated.

    """

    # Open hdf5 log file
    f = h5py.File(log_file, 'r')

    # Retrieve time steps
    time_steps = f['time_steps'].keys()
    time_steps = np.sort(np.array([float(i) for i in time_steps]))
    num_outer_time_steps = f.attrs['num_outer_time_steps']

    # Retrieve mesh domain size and number of groups
    if variable in ['pin_powers', 'pin_cell_kappa_fission',
                    'pin_cell_shape']:
        domains = f['pin_cell_mesh'].attrs['dimension']
    elif variable in ['assembly_powers', 'assembly_kappa_fission',
                      'assembly_shape']:
        domains = f['assembly_mesh'].attrs['dimension']
    elif variable in ['flux', 'adjoint_flux', 'precursors',
                      'kappa_fission', 'current_balance', 'current_computed']:
        domains = f['mesh'].attrs['dimension']

    if 'precursors' in variable:
        num_groups = f.attrs['num_delayed_groups']
    elif 'powers' in variable:
        num_groups = 1
    else:
        num_groups = f['energy_groups'].attrs['num_groups']

    # Retrieve the spatial data
    #spatial_data = np.zeros((num_outer_time_steps, domains[2], domains[1],
    #                         domains[0], num_groups))
    outer_time_steps = []
    outer_i = 0
    spatial_data = []
    for i,step in enumerate(time_steps):
        if u'flux' in f['time_steps'][str(step)].keys():
            outer_time_steps.append(step)
            data = f['time_steps'][str(step)][variable][:]
            data.shape = (domains[2], domains[1], domains[0], num_groups)
            #spatial_data[outer_i] = data
            spatial_data.append(data)
            outer_i += 1

    # Get the limits used to scale all the figures
    data_lim = [-np.finfo(np.float).max, np.finfo(np.float).max]
    for i,step in enumerate(outer_time_steps):
        if plane == 'xy':
            data_lim[0] = np.max([data_lim[0], np.min(spatial_data[i][plane_num, :, :, group])])
            data_lim[1] = np.min([data_lim[1], np.max(spatial_data[i][plane_num, :, :, group])])
        elif plane == 'xz':
            data_lim[0] = np.max([data_lim[0], np.min(spatial_data[i, :, plane_num, :, group])])
            data_lim[1] = np.min([data_lim[1], np.max(spatial_data[i, :, plane_num, :, group])])
        elif plane == 'yz':
            data_lim[0] = np.max([data_lim[0], np.min(spatial_data[i, :, :, plane_num, group])])
            data_lim[1] = np.min([data_lim[1], np.max(spatial_data[i, :, :, plane_num, group])])

    # Plot the spatial data
    for i,step in enumerate(outer_time_steps):
        fig = plt.figure()
        if bar:
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = fig.add_subplot(111)
        if plane == 'xy':
            if bar:
                nx = domains[0]
                ny = domains[1]
                zeros = np.zeros(nx*ny)
                ones = np.ones(nx*ny)

                color_map = np.copy(spatial_data[i][plane_num, :, :, group].flatten())
                fracs = color_map.astype(float)/normalize
                norm = matplotlib.colors.Normalize(0.0, 1.0)
                colors = matplotlib.cm.jet(norm(fracs))

                cax = ax.bar3d(np.tile(range(nx), ny), np.repeat(range(ny), nx), zeros, ones, ones, spatial_data[i][plane_num, :, :, group].flatten(), color=colors)
                #cax = ax.bar3d(np.tile(range(nx), ny), np.repeat(range(ny), nx), zeros, ones, ones, spatial_data[i][plane_num, :, :, group].flatten())
                ax.set_zlim(0,normalize)
            else:
                cax = ax.imshow(spatial_data[i][plane_num, :, :, group],
                                vmin=data_lim[0], vmax=data_lim[1], interpolation='nearest')
        elif plane == 'xz':
            cax = ax.imshow(spatial_data[i, :, plane_num, :, group],
                            vmin=data_lim[0], vmax=data_lim[1], interpolation='nearest')
        elif plane == 'yz':
            cax = ax.imshow(spatial_data[i, :, :, plane_num, group],
                            vmin=data_lim[0], vmax=data_lim[1], interpolation='nearest')
        #fig.colorbar(cax)
        ax.set_title('{} time step {} @ t = {:.5f} s'.format(variable.replace('_', ' '), i, step))
        plt.savefig('{}/{}_{}_plane_{}_{:.5f}_s.png'.format(directory, variable, plane, plane_num, step))
        plt.close()

    f.close()
    return None

def scalar_plot(variable, log_file, variable_twin=None, directory='.',
                axis=None, **kwargs):
    """Creates a figure of a scalar variable from a transient solve.

    Parameters
    ----------
    variable : str
        Variable to plot.
    log_file : str
        Location of log file.
    variable_twin : str
        Variable to plot on twinx axis. Default is None.
    directory : str
        Directory to save pngs. Default is '.'.
    axis : matplotlib.axes, optional
        A previously generated axis to use for plotting. If not specified,
        a new axis and figure will be generated.
    **kwargs
        All keyword arguments are passed to
        :func:`matplotlib.pyplot.figure`.

    Returns
    -------
    fig : matplotlib.figure.Figure or None
        If axis is None, then a Matplotlib Figure of the generated spatial
        variable will be returned. Otherwise, a value of None will be returned
        as the figure and axes have already been generated.

    """

    scalar_data = get_scalar_data(variable, log_file, variable_twin, directory)

    # Plot the spatial data
    if axis == None:
        fig, ax = plt.subplots()
    else:
        fig = None
        ax = axis

    ax.plot(scalar_data[:,0], scalar_data[:,1])
    ax.set_xlabel('time (s)', color='k')
    ax.set_ylabel(variable.replace('_', ' '), color='k')
    ax.tick_params('y', colors='k')

    if variable_twin is not None:
        twin_ax = ax.twinx()
        twin_ax.plot(scalar_data[:,0], scalar_data[:,2], color='r')
        twin_ax.set_ylabel(variable_twin.replace('_', ' '), color='r')
        twin_ax.tick_params('y', colors='r')
        plt.gcf().tight_layout()

    if variable_twin is None:
        plt.savefig('{}/{}.png'.format(directory, variable))
    else:
        plt.savefig('{}/{}-{}.png'.format(directory, variable, variable_twin))

    return fig, ax
