from collections.abc import Iterable
import os
from pathlib import Path
from numbers import Integral
import time
import warnings
import subprocess

import h5py
import numpy as np

import openmc
from openmc.checkvalue import check_type, check_value, check_iterable_type, \
    check_length
import openmc.deplete as dep
from openmc.data.library import DataLibrary
from openmc.exceptions import DataError, InvalidIDError, SetupError


class Model:
    """Model container.

    This class can be used to store instances of :class:`openmc.Geometry`,
    :class:`openmc.Materials`, :class:`openmc.Settings`,
    :class:`openmc.Tallies`, and :class:`openmc.Plots`, thus making a complete
    model. The :meth:`Model.export_to_xml` method will export XML files for all
    attributes that have been set. If the :attr:`Model.materials` attribute is
    not set, it will attempt to create a ``materials.xml`` file based on all
    materials appearing in the geometry.

    Parameters
    ----------
    geometry : openmc.Geometry, optional
        Geometry information
    materials : openmc.Materials, optional
        Materials information
    settings : openmc.Settings, optional
        Settings information
    tallies : openmc.Tallies, optional
        Tallies information
    plots : openmc.Plots, optional
        Plot information
    chain_file : str or Path, optional
        Path to the depletion chain XML file.  Defaults to the chain
        found under the ``depletion_chain`` in the
        :envvar:`OPENMC_CROSS_SECTIONS` environment variable if it exists. If a
        str is provided it will be converted to a Path object.
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV].
        If not given, values will be pulled from the ``chain_file``.

    Attributes
    ----------
    geometry : openmc.Geometry
        Geometry information
    materials : openmc.Materials
        Materials information
    settings : openmc.Settings
        Settings information
    tallies : openmc.Tallies
        Tallies information
    plots : openmc.Plots
        Plot information
    chain_file : str or Path
        Path to the depletion chain XML file.  Defaults to the chain
        found under the ``depletion_chain`` in the
        :envvar:`OPENMC_CROSS_SECTIONS` environment variable if it exists. If a
        str is provided it will be converted to a Path object.
    fission_q : dict
        Dictionary of nuclides and their fission Q values [eV].
        If not given, values will be pulled from the ``chain_file``.

    """

    def __init__(self, geometry=None, materials=None, settings=None,
                 tallies=None, plots=None, chain_file=None, fission_q=None):
        self.geometry = openmc.Geometry()
        self.materials = openmc.Materials()
        self.settings = openmc.Settings()
        self.tallies = openmc.Tallies()
        self.plots = openmc.Plots()

        if geometry is not None:
            self.geometry = geometry
        if materials is not None:
            self.materials = materials
        if settings is not None:
            self.settings = settings
        if tallies is not None:
            self.tallies = tallies
        if plots is not None:
            self.plots = plots

        self.chain_file = chain_file
        self.fission_q = fission_q

        if materials is None:
            mats = self.geometry.get_all_materials().values()
        else:
            mats = self.materials
        self._materials_by_id = {mat.id: mat for mat in mats}
        self._materials_by_name = {mat.name: mat for mat in mats}
        cells = self.geometry.get_all_cells()
        self._cells_by_id = {cell.id: cell for cell in cells.values()}
        self._cells_by_name = {cell.name: cell for cell in cells.values()}

    @property
    def geometry(self):
        return self._geometry

    @property
    def materials(self):
        return self._materials

    @property
    def settings(self):
        return self._settings

    @property
    def tallies(self):
        return self._tallies

    @property
    def plots(self):
        return self._plots

    @property
    def chain_file(self):
        return self._chain_file

    @property
    def fission_q(self):
        return self._fission_q

    @property
    def C_init(self):
        return openmc.lib.LIB_INIT

    @geometry.setter
    def geometry(self, geometry):
        check_type('geometry', geometry, openmc.Geometry)
        self._geometry = geometry

    @materials.setter
    def materials(self, materials):
        check_type('materials', materials, Iterable, openmc.Material)
        if isinstance(materials, openmc.Materials):
            self._materials = materials
        else:
            del self._materials[:]
            for mat in materials:
                self._materials.append(mat)

    @settings.setter
    def settings(self, settings):
        check_type('settings', settings, openmc.Settings)
        self._settings = settings

    @tallies.setter
    def tallies(self, tallies):
        check_type('tallies', tallies, Iterable, openmc.Tally)
        if isinstance(tallies, openmc.Tallies):
            self._tallies = tallies
        else:
            del self._tallies[:]
            for tally in tallies:
                self._tallies.append(tally)

    @plots.setter
    def plots(self, plots):
        check_type('plots', plots, Iterable, openmc.Plot)
        if isinstance(plots, openmc.Plots):
            self._plots = plots
        else:
            del self._plots[:]
            for plot in plots:
                self._plots.append(plot)

    @chain_file.setter
    def chain_file(self, chain_file):
        check_type('chain_file', chain_file, (type(None), str, Path))
        if isinstance(chain_file, str):
            self._chain_file = Path(chain_file).resolve()
        elif isinstance(chain_file, Path):
            self._chain_file = chain_file.resolve()
        else:
            self._chain_file = None

    @fission_q.setter
    def fission_q(self, fission_q):
        check_type('fission_q', fission_q, (type(None), dict))
        self._fission_q = fission_q

    @classmethod
    def from_xml(cls, geometry='geometry.xml', materials='materials.xml',
                 settings='settings.xml'):
        """Create model from existing XML files
        When initializing this way, the user must manually load plots, tallies,
        the chain_file and fission_q attributes.

        Parameters
        ----------
        geometry : str
            Path to geometry.xml file
        materials : str
            Path to materials.xml file
        settings : str
            Path to settings.xml file

        Returns
        -------
        openmc.model.Model
            Model created from XML files

        """
        materials = openmc.Materials.from_xml(materials)
        geometry = openmc.Geometry.from_xml(geometry, materials)
        settings = openmc.Settings.from_xml(settings)
        return cls(geometry, materials, settings)

    def init_C_api(self, threads=None, geometry_debug=False, restart_file=None,
                   tracks=False, output=True, event_based=False):
        """Initializes the model in memory via the C-API

        Parameters
        ----------
        threads : int, optional
            Number of OpenMP threads. If OpenMC is compiled with OpenMP threading
            enabled, the default is implementation-dependent but is usually equal to
            the number of hardware threads available (or a value set by the
            :envvar:`OMP_NUM_THREADS` environment variable).
        geometry_debug : bool, optional
            Turn on geometry debugging during simulation. Defaults to False.
        restart_file : str, optional
            Path to restart file to use
        tracks : bool, optional
            Write tracks for all particles. Defaults to False.
        output : bool
            Capture OpenMC output from standard out
        event_based : bool, optional
            Turns on event-based parallelism, instead of default history-based
        """

        # TODO: right now the only way to set most of the above parameters via
        # the C-API are at initialization time despite use-cases existing to
        # set them for individual runs. For now this functionality is exposed
        # where it exists (here in init), but in the future the functionality
        # should be exposed so that it can be accessed via model.run(...)

        # TODO: the output flag is not yet implemented for the openmc.lib.init
        # command. This will be the subject of future work.

        args = []

        if isinstance(threads, Integral) and threads > 0:
            args += ['-s', str(threads)]

        if geometry_debug:
            args.append('-g')

        if event_based:
            args.append('-e')

        if isinstance(restart_file, str):
            args += ['-r', restart_file]

        if tracks:
            args.append('-t')

        self.clear_C_api()

        if dep.comm.rank == 0:
            self.export_to_xml()
        dep.comm.barrier()
        # TODO: Implement the output flag somewhere on the C++ side
        openmc.lib.init(args=args, intracomm=dep.comm)

    def clear_C_api(self):
        """Finalize simulation and free memory allocated for the C-API"""
        openmc.lib.finalize()

    def deplete(self, timesteps, chain_file=None, method='cecm',
                fission_q=None, final_step=True, directory='.',
                **kwargs):
        """Deplete model using specified timesteps/power

        Parameters
        ----------
        timesteps : iterable of float
            Array of timesteps in units of [s]. Note that values are not
            cumulative.
        chain_file : str, optional
            Path to the depletion chain XML file.  Defaults to the chain
            found under the ``depletion_chain`` in the
            :envvar:`OPENMC_CROSS_SECTIONS` environment variable if it exists.
        method : str, optional
             Integration method used for depletion (e.g., 'cecm', 'predictor').
             Defaults to 'cecm'.
        fission_q : dict, optional
            Dictionary of nuclides and their fission Q values [eV].
            If not given, values will be pulled from the ``chain_file``.
            Defaults to pulling from the ``chain_file``.
        final_step : bool, optional
            Indicate whether or not a transport solve should be run at the end
            of the last timestep. Defaults to running this transport solve.
        directory : str, optional
            Directory to write XML files to. If it doesn't exist already, it
            will be created. Defaults to the current working directory
        **kwargs
            Keyword arguments passed to integration function (e.g.,
            :func:`openmc.deplete.integrator.cecm`)

        """

        # To keep Model.deplete(...) API compatibility, we will allow the
        # chain_file and fission_q params to be set if provided while we set
        # the depletion operator
        if chain_file is not None:
            this_chain_file = Path(chain_file).resolve()
            warnings.warn("The chain_file argument of Model.deplete(...) "
                          "has been deprecated and may be removed in a "
                          "future version. The Model.chain_file should be"
                          "used instead.", DeprecationWarning)
        else:
            this_chain_file = self.chain_file
        if fission_q is not None:
            warnings.warn("The fission_q argument of Model.deplete(...) "
                          "has been deprecated and may be removed in a "
                          "future version. The Model.fission_q should be"
                          "used instead.", DeprecationWarning)
            this_fission_q = fission_q
        else:
            this_fission_q = fission_q

        # Create directory if required
        d = Path(directory)
        if not d.is_dir():
            d.mkdir(parents=True)
        start_dir = Path.cwd()
        os.chdir(d)

        depletion_operator = \
            dep.Operator(self.geometry, self.settings,
                         str(this_chain_file.absolute()),
                         fission_q=this_fission_q)
        # Tell depletion_operator.finalize NOT to clear C-API memory when it is
        # done
        depletion_operator.cleanup_when_done = False

        # Set up the integrator
        integrator_class = dep.integrators.integrator_factory(method)
        integrator = integrator_class(depletion_operator,
                                      timesteps, **kwargs)

        # Now perform the depletion
        integrator.integrate(final_step)

        # If we did not perform a transport calculation on the final step, then
        # make the code update the C-API material inventory
        if not final_step:
            depletion_operator._update_materials()

        # Now make the python Materials match the C-API material data
        for mat_id, mat in self._materials_by_id.items():
            if mat.depletable:
                # Get the C data
                c_mat = openmc.lib.materials[mat_id]
                nuclides, densities = c_mat._get_densities()
                # And now we can remove isotopes and add these ones in
                atom_density = 0.
                for nuc, density in zip(nuclides, densities):
                    mat.remove_nuclide(nuc)  # Replace if it's there
                    mat.add_nuclide(nuc, density)
                    atom_density += density
                mat.set_density('atom/b-cm', atom_density)

        os.chdir(start_dir)

    def export_to_xml(self, directory='.'):
        """Export model to XML files.

        Parameters
        ----------
        directory : str
            Directory to write XML files to. If it doesn't exist already, it
            will be created.

        """
        # Create directory if required
        d = Path(directory)
        if not d.is_dir():
            d.mkdir(parents=True)

        self.settings.export_to_xml(d)
        self.geometry.export_to_xml(d)

        # If a materials collection was specified, export it. Otherwise, look
        # for all materials in the geometry and use that to automatically build
        # a collection.
        if self.materials:
            self.materials.export_to_xml(d)
        else:
            materials = openmc.Materials(self.geometry.get_all_materials()
                                         .values())
            materials.export_to_xml(d)

        if self.tallies:
            self.tallies.export_to_xml(d)
        if self.plots:
            self.plots.export_to_xml(d)

    def import_properties(self, filename):
        """Import physical properties

        .. versionchanged:: 0.12.3
            This method now updates values as loaded in memory with the C-API

        Parameters
        ----------
        filename : str
            Path to properties HDF5 file

        See Also
        --------
        openmc.lib.export_properties

        """
        cells = self.geometry.get_all_cells()
        materials = self.geometry.get_all_materials()

        with h5py.File(filename, 'r') as fh:
            cells_group = fh['geometry/cells']

            # Make sure number of cells matches
            n_cells = fh['geometry'].attrs['n_cells']
            if n_cells != len(cells):
                raise ValueError("Number of cells in properties file doesn't "
                                 "match current model.")

            # Update temperatures for cells filled with materials
            for name, group in cells_group.items():
                cell_id = int(name.split()[1])
                cell = cells[cell_id]
                if cell.fill_type in ('material', 'distribmat'):
                    cell.temperature = group['temperature'][()]
                    if self.C_init:
                        C_cell = openmc.lib.cells[cell_id]
                        C_cell.set_temperature(group['temperature'][()])

            # Make sure number of materials matches
            mats_group = fh['materials']
            n_cells = mats_group.attrs['n_materials']
            if n_cells != len(materials):
                raise ValueError("Number of materials in properties file "
                                 "doesn't match current model.")

            # Update material densities
            for name, group in mats_group.items():
                mat_id = int(name.split()[1])
                atom_density = group.attrs['atom_density']
                materials[mat_id].set_density('atom/b-cm', atom_density)
                if self.C_init:
                    C_mat = openmc.lib.materials[mat_id]
                    C_mat.set_density(atom_density, 'atom/b-cm')

    def run(self, particles=None, threads=None, geometry_debug=False,
            restart_file=None, tracks=False, output=True, cwd='.',
            openmc_exec='openmc', mpi_args=None, event_based=False):
        """Runs OpenMC. If the C-API has been initialized, then the C-API is
        used, otherwise, this method creates the XML files and runs OpenMC via
        a system call. In both cases this method returns the path to the last
        statepoint file generated.

        .. versionchanged:: 0.12
            Instead of returning the final k-effective value, this function now
            returns the path to the final statepoint written.

        .. versionchanged:: 0.12.3
            This method can utilize the C-API for execution

        Parameters
        ----------
        particles : int, optional
            Number of particles to simulate per generation.
        threads : int, optional
            Number of OpenMP threads. If OpenMC is compiled with OpenMP
            threading enabled, the default is implementation-dependent but is
            usually equal to the number of hardware threads available (or a
            value set by the :envvar:`OMP_NUM_THREADS` environment variable).
        geometry_debug : bool, optional
            Turn on geometry debugging during simulation. Defaults to False.
        restart_file : str, optional
            Path to restart file to use
        tracks : bool, optional
            Write tracks for all particles. Defaults to False.
        output : bool, optional
            Capture OpenMC output from standard out
        cwd : str, optional
            Path to working directory to run in. Defaults to the current
            working directory.
        openmc_exec : str, optional
            Path to OpenMC executable. Defaults to 'openmc'.
        mpi_args : list of str, optional
            MPI execute command and any additional MPI arguments to pass,
            e.g. ['mpiexec', '-n', '8'].
        event_based : bool, optional
            Turns on event-based parallelism, instead of default history-based

        Returns
        -------
        Path
            Path to the last statepoint written by this run
            (None if no statepoint was written)

        """

        # Setting tstart here ensures we don't pick up any pre-existing
        # statepoint files in the output directory
        tstart = time.time()
        last_statepoint = None

        # Operate in the provided working directory
        with openmc.change_directory(Path(cwd)):
            if self.C_init:
                # Handle the run options as applicable
                # First dont allow ones that must be set via init
                for arg_name, arg, default in zip(
                    ['threads', 'geometry_debug', 'restart_file', 'tracks'],
                    [threads, geometry_debug, restart_file, tracks],
                    [None, False, None, False]):
                    if arg != default:
                        msg = "{} must be set via Model.c_init(...)".format(
                            arg_name)
                        raise ValueError(msg)

                if particles is not None:
                    init_particles = openmc.lib.settings.particles
                    if isinstance(particles, Integral) and particles > 0:
                        openmc.lib.settings.particles = particles

                # If we dont want output, make the verbosity quiet
                if not output:
                    init_verbosity = openmc.lib.settings.verbosity
                    if not output:
                        openmc.lib.settings.verbosity = 1

                # Event-based can be set at init-time or on a case-basis.
                # Handle the argument here.
                # TODO This will be dealt with in a future change to the C-API

                # Then run using the C-API
                openmc.lib.run()

                # Reset changes for the openmc.run kwargs handling
                if particles is not None:
                    openmc.lib.settings.particles = init_particles
                if not output:
                    # Then re-set the initial verbosity
                    openmc.lib.settings.verbosity = init_verbosity

            else:
                # Then run via the command line
                self.export_to_xml()
                openmc.run(particles, threads, geometry_debug, restart_file,
                           tracks, output, Path('.'), openmc_exec, mpi_args,
                           event_based)

            # Get output directory and return the last statepoint written
            if self.settings.output and 'path' in self.settings.output:
                output_dir = Path(self.settings.output['path'])
            else:
                output_dir = Path.cwd()
            for sp in output_dir.glob('statepoint.*.h5'):
                mtime = sp.stat().st_mtime
                if mtime >= tstart:  # >= allows for poor clock resolution
                    tstart = mtime
                    last_statepoint = sp
        return last_statepoint

    def calculate_volume(self, threads=None, output=True, cwd='.',
                         openmc_exec='openmc', mpi_args=None,
                         apply_volumes=True):
        """Runs an OpenMC stochastic volume calculation and, if requested,
        applies volumes to the model

        Parameters
        ----------
        threads : int, optional
            Number of OpenMP threads. If OpenMC is compiled with OpenMP
            threading enabled, the default is implementation-dependent but is
            usually equal to the number of hardware threads available (or a
            value set by the :envvar:`OMP_NUM_THREADS` environment variable).
            This currenty only applies to the case when not using the C-API.
        output : bool, optional
            Capture OpenMC output from standard out
        openmc_exec : str, optional
            Path to OpenMC executable. Defaults to 'openmc'.
            This only applies to the case when not using the C-API.
        mpi_args : list of str, optional
            MPI execute command and any additional MPI arguments to pass,
            e.g. ['mpiexec', '-n', '8'].
            This only applies to the case when not using the C-API.
        cwd : str, optional
            Path to working directory to run in. Defaults to the current
            working directory.
        apply_volumes : bool, optional
            Whether apply the volume calculation results from this calculation
            to the model. Defaults to applying the volumes.
        """

        if len(self.settings.volume_calculation) == 0:
            # Then there is no volume calculation specified
            raise ValueError("The Settings.volume_calculation attribute must"
                             " be specified before executing this method!")

        with openmc.change_directory(Path(cwd)):
            if self.C_init:
                if threads is not None:
                    msg = "Threads must be set via Model.c_init(...)"
                    raise ValueError(msg)
                if mpi_args is not None:
                    msg = "The MPI environment must be set otherwise such as" \
                        "with the call to mpi4py"
                    raise ValueError(msg)

                # Apply the output settings
                if not output:
                    init_verbosity = openmc.lib.settings.verbosity
                    if not output:
                        openmc.lib.settings.verbosity = 1

                # Compute the volumes
                openmc.lib.calculate_volumes()

                # Reset the output verbosity
                if not output:
                    openmc.lib.settings.verbosity = init_verbosity
            else:
                openmc.calculate_volumes(threads=threads, output=output,
                                         openmc_exec=openmc_exec,
                                         mpi_args=mpi_args)

            # Now we apply the volumes
            if apply_volumes:
                # Load the results
                f_names = \
                    ["volume_{}.h5".format(i + 1)
                     for i in range(len(self.settings.volume_calculations))]
                vol_calcs = [openmc.VolumeCalculation.load_results(f_name)
                             for f_name in f_names]
                # And now we can add them to the model
                for vol_calc in vol_calcs:
                    # First add them to the Python side
                    self.geometry.add_volume_information(vol_calc)

                    # And now repeat for the C-API
                    if vol_calc.domain == 'material':
                        # Then we can do this in the C-API
                        for domain_id in vol_calc.domains:
                            self.update_material_volumes(
                                [domain_id], vol_calc.volumes[domain_id])

    def plot_geometry(self, output=True, cwd='.', openmc_exec='openmc',
                      convert=True, convert_exec='convert'):
        """Creates plot images as specified by the Model.plots attribute

        If convert is True, this function requires that a program is installed
        to convert PPM files to PNG files. Typically, that would be
        `ImageMagick <https://www.imagemagick.org>`_ which includes a
        `convert` command.

        Parameters
        ----------
        output : bool, optional
            Capture OpenMC output from standard out
        cwd : str, optional
            Path to working directory to run in. Defaults to the current
            working directory.
        openmc_exec : str, optional
            Path to OpenMC executable. Defaults to 'openmc'.
            This only applies to the case when not using the C-API.
        convert : bool, optional
            Whether or not to attempt to convert from PPM to PNG
        convert_exec : str, optional
            Command that can convert PPM files into PNG files
        """

        if len(self.plots) == 0:
            # Then there is no volume calculation specified
            raise ValueError("The Model.plots attribute must be specified "
                             "before executing this method!")

        with openmc.change_directory(Path(cwd)):
            if self.C_init:
                # Apply the output settings
                if not output:
                    init_verbosity = openmc.lib.settings.verbosity
                    if not output:
                        openmc.lib.settings.verbosity = 1

                # Compute the volumes
                openmc.lib.plot_geometry()

                # Reset the output verbosity
                if not output:
                    openmc.lib.settings.verbosity = init_verbosity
            else:
                openmc.plot_geometry(output=output, openmc_exec=openmc_exec)

            if convert:
                for p in self.plots:
                    if p.filename is not None:
                        ppm_file = f'{p.filename}.ppm'
                    else:
                        ppm_file = f'plot_{p.id}.ppm'
                    png_file = ppm_file.replace('.ppm', '.png')
                    subprocess.check_call([convert_exec, ppm_file, png_file])

    def _change_py_C_attribs(self, names_or_ids, value, obj_type, attrib_name,
                             density_units='atom/b-cm'):
        # Method to do the same work whether it is a cell or material and
        # a temperature or volume
        check_type('names_or_ids', names_or_ids, Iterable, (Integral, str))
        check_type('obj_type', obj_type, str)
        obj_type = obj_type.lower()
        check_value('obj_type', obj_type, ('material', 'cell'))
        check_value('attrib_name', attrib_name,
                    ('temperature', 'volume', 'density', 'rotation',
                     'translation'))
        # The C-API only allows setting density units of atom/b-cm and g/cm3
        check_value('density_units', density_units, ('atom/b-cm', 'g/cm3'))
        # The C-API has no way to set cell volume so lets raise an exception
        if obj_type == 'cell' and attrib_name == 'volume':
            raise NotImplementedError(
                'Setting a Cell volume is not supported!')
        # Same with setting temperatures, TODO: update C-API for this
        if obj_type == 'material' and attrib_name == 'temperature':
            raise NotImplementedError(
                'Setting a Material temperature is not yet supported!')
        # And some items just dont make sense
        if obj_type == 'cell' and attrib_name == 'density':
            raise ValueError('Cannot set a Cell density!')
        if obj_type == 'material' and attrib_name in ('rotation',
                                                      'translation'):
            raise ValueError('Cannot set a material rotation/translation!')

        # Set the
        if obj_type == 'cell':
            by_name = self._cells_by_name
            by_id = self._cells_by_id
            C_by_id = openmc.lib.cells
        else:
            by_name = self._materials_by_name
            by_id = self._materials_by_id
            C_by_id = openmc.lib.materials

        # Get the list of ids to use if converting from names and accepting
        # only values that have actual ids
        ids = [None] * len(names_or_ids)
        for i, name_or_id in enumerate(names_or_ids):
            if isinstance(name_or_id, Integral):
                if name_or_id in by_id:
                    ids[i] = int(name_or_id)
                else:
                    msg = '{} ID {} is not present in the model!'.format(
                        obj_type.capitalize(), name_or_id)
                    raise InvalidIDError(msg)
            elif isinstance(name_or_id, str):
                if name_or_id in by_name:
                    ids[i] = by_name[name_or_id].id
                else:
                    msg = '{} {} is not present in the model!'.format(
                        obj_type.capitalize(), name_or_id)
                    raise InvalidIDError(msg)

        # Now perform the change to both python and C-API
        for id_ in ids:
            obj = by_id[id_]
            if attrib_name == 'rotation':
                obj.rotation = value
            elif attrib_name == 'translation':
                obj.translation = value
            elif attrib_name == 'volume':
                obj.volume = value
            elif attrib_name == 'temperature':
                obj.temperature = value
            elif attrib_name == 'density':
                obj.set_density(density_units, value)
            # Next lets keep what is in C-API memory up to date as well
            if self.C_init:
                C_obj = C_by_id[id_]
                if attrib_name == 'rotation':
                    C_obj.rotation = value
                elif attrib_name == 'translation':
                    C_obj.translation = value
                elif attrib_name == 'volume':
                    C_obj.volume = value
                elif attrib_name == 'temperature':
                    C_obj.set_temperature(value)
                elif attrib_name == 'density':
                    C_obj.set_density(value, density_units)

    def rotate_cells(self, names_or_ids, vector):
        """Rotate the identified cell(s) by the specified rotation vector.
        The rotation is only applied to cells filled with a universe.

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The cell names (if str) or id (if int) that are to be translated
            or rotated. This parameter can include a mix of names and ids.
        vector : Iterable of float
            The rotation vector of length 3 to apply. This array specifies the
            angles in degrees about the x, y, and z axes, respectively.

        """

        self._change_py_C_attribs(names_or_ids, vector, 'cell', 'rotation')

    def translate_cells(self, names_or_ids, vector):
        """Translate the identified cell(s) by the specified translation vector.
        The translation is only applied to cells filled with a universe.

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The cell names (if str) or id (if int) that are to be translated
            or rotated. This parameter can include a mix of names and ids.
        vector : Iterable of float
            The translation vector of length 3 to apply. This array specifies
            the x, y, and z dimensions of the translation.

        """

        self._change_py_C_attribs(names_or_ids, vector, 'cell', 'translation')

    def update_densities(self, names_or_ids, density, density_units='atom/b-cm'):
        """Update the density of a given set of materials to a new value

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The material names (if str) or id (if int) that are to be updated.
            This parameter can include a mix of names and ids.
        density : float
            The density to apply in the units specified by `density_units`
        density_units : {'atom/b-cm', 'g/cm3'}, optional
            Units for `density`. Defaults to 'atom/b-cm'

        """

        self._change_py_C_attribs(names_or_ids, density, 'material', 'density',
                                  density_units)

    def update_cell_temperatures(self, names_or_ids, temperature):
        """Update the temperature of a set of cells to the given value

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The cell names (if str) or id (if int) that are to be updated.
            This parameter can include a mix of names and ids.
        temperature : float
            The temperature to apply in units of Kelvin

        """

        self._change_py_C_attribs(names_or_ids, temperature, 'cell',
                                  'temperature')

    def update_material_temperatures(self, names_or_ids, temperature):
        """Update the temperature of a set of materials to the given value

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The material names (if str) or id (if int) that are to be updated.
            This parameter can include a mix of names and ids.
        temperature : float
            The temperature to apply in units of Kelvin

        """

        self._change_py_C_attribs(names_or_ids, temperature, 'material',
                                  'temperature')

    def update_material_volumes(self, names_or_ids, volume):
        """Update the volume of a set of materials to the given value

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The material names (if str) or id (if int) that are to be updated.
            This parameter can include a mix of names and ids.
        volume : float
            The volume to apply in units of cm^3

        """

        self._change_py_C_attribs(names_or_ids, volume, 'material', 'volume')
