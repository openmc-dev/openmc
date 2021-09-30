from collections.abc import Iterable
import operator
import os
from pathlib import Path
from numbers import Integral
import time
import warnings
import subprocess
from contextlib import contextmanager

import h5py

import openmc
from openmc.checkvalue import check_type, check_value
from openmc.exceptions import InvalidIDError


@contextmanager
def _change_directory(working_dir):
    """A context manager for executing in a provided working directory"""
    start_dir = Path().cwd()
    Path.mkdir(working_dir, exist_ok=True)
    os.chdir(working_dir)
    try:
        yield
    finally:
        os.chdir(start_dir)


class Model:
    """Model container.

    This class can be used to store instances of :class:`openmc.Geometry`,
    :class:`openmc.Materials`, :class:`openmc.Settings`,
    :class:`openmc.Tallies`, and :class:`openmc.Plots`, thus making a complete
    model. The :meth:`Model.export_to_xml` method will export XML files for all
    attributes that have been set. If the :attr:`Model.materials` attribute is
    not set, it will attempt to create a ``materials.xml`` file based on all
    materials appearing in the geometry.

    .. versionchanged:: 0.13.0

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
    intracomm : mpi4py.MPI.Intracomm or openmc.DummyCommunicator
        MPI intracommunicator; this defaults to the mpi4py world communicator
        from if present, or a DummyCommunicator otherwise. If an alternative
        communicator is desired, this parameter should be modified accordingly.

    """

    def __init__(self, geometry=None, materials=None, settings=None,
                 tallies=None, plots=None, intracomm=None):
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

        self.intracomm = openmc.comm

        # Store dictionaries to the materials and cells by ID and names
        if materials is None:
            mats = self.geometry.get_all_materials().values()
        else:
            mats = self.materials
        cells = self.geometry.get_all_cells()
        # Get the ID maps
        self._materials_by_id = {mat.id: mat for mat in mats}
        self._cells_by_id = {cell.id: cell for cell in cells.values()}

        # Get the names maps, but since names are not unique, store a list for
        # each name key. In this way when the user requests a change by a name,
        # the change will be applied to all of the same name.
        self._cells_by_name = {}
        for cell in cells.values():
            if cell.name not in self._cells_by_name:
                self._cells_by_name[cell.name] = [cell]
            else:
                self._cells_by_name[cell.name].append(cell)
        self._materials_by_name = {}
        for mat in mats:
            if mat.name not in self._materials_by_name:
                self._materials_by_name[mat.name] = [mat]
            else:
                self._materials_by_name[mat.name].append(mat)

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
    def intracomm(self):
        return self._intracomm

    @property
    def is_initialized(self):
        # TODO: Replace openmc.lib.is_initialized with a direct ctypes access
        # to simulation::initialized
        return openmc.lib.is_initialized

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

    @intracomm.setter
    def intracomm(self, intracomm):
        check_type('intracomm', intracomm, type(openmc.comm))
        self._intracomm = intracomm

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

    def init_lib(self, threads=None, geometry_debug=False, restart_file=None,
                 tracks=False, output=True, event_based=False):
        """Initializes the model in memory via the C API

        .. versionadded:: 0.13.0

        Parameters
        ----------
        threads : int, optional
            Number of OpenMP threads. If OpenMC is compiled with OpenMP
            threading enabled, the default is implementation-dependent but is
            usually equal to the number of hardware threads available
            (or a value set by the :envvar:`OMP_NUM_THREADS` environment
            variable).
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
        # the C API are at initialization time despite use-cases existing to
        # set them for individual runs. For now this functionality is exposed
        # where it exists (here in init), but in the future the functionality
        # should be exposed so that it can be accessed via model.run(...)

        args = openmc.process_CLI_arguments(
            volume=False, geometry_debug=geometry_debug,
            restart_file=restart_file, threads=threads, tracks=tracks,
            event_based=event_based)
        # Args adds the openmc_exec command in the first entry; remove it
        args = args[1:]

        self.finalize_lib()

        if self.intracomm.rank == 0:
            self.export_to_xml()
        self.intracomm.barrier()

        if isinstance(self.intracomm, openmc.DummyCommunicator):
            # openmc.lib.init does not accept DummyCommunicator, and importing
            # the DummyCommunicator class there is overkill. Filter it here
            intracomm = None
        else:
            intracomm = self.intracomm
        openmc.lib.init(args=args, intracomm=intracomm, output=output)

    def finalize_lib(self):
        """Finalize simulation and free memory allocated for the C API

        .. versionadded:: 0.13.0

        """

        openmc.lib.finalize()

    def deplete(self, timesteps, method='cecm', final_step=True,
                operator_kwargs=None, directory='.', output=True,
                **integrator_kwargs):
        """Deplete model using specified timesteps/power

        .. versionchanged:: 0.13.0

        Parameters
        ----------
        timesteps : iterable of float
            Array of timesteps in units of [s]. Note that values are not
            cumulative.
        method : str, optional
             Integration method used for depletion (e.g., 'cecm', 'predictor').
             Defaults to 'cecm'.
        final_step : bool, optional
            Indicate whether or not a transport solve should be run at the end
            of the last timestep. Defaults to running this transport solve.
        operator_kwargs : dict
            Keyword arguments passed to the depletion Operator initializer
            (e.g., :func:`openmc.deplete.Operator`)
        directory : str, optional
            Directory to write XML files to. If it doesn't exist already, it
            will be created. Defaults to the current working directory
        output : bool
            Capture OpenMC output from standard out
        integrator_kwargs : dict
            Remaining keyword arguments passed to the depletion Integrator
            initializer (e.g., :func:`openmc.deplete.integrator.cecm`).

        """

        if operator_kwargs is None:
            op_kwargs = {}
        elif isinstance(operator_kwargs, dict):
            op_kwargs = operator_kwargs
        else:
            msg = "operator_kwargs must be a dict or None"
            raise ValueError(msg)

        # Import openmc.deplete here so the Model can be used even if the
        # shared library is unavailable.
        import openmc.deplete as dep

        # Store whether or not the library was initialized when we started
        started_initialized = self.is_initialized

        with _change_directory(Path(directory)):
            with openmc.lib.quiet_dll(output):
                depletion_operator = \
                    dep.Operator(self.geometry, self.settings, **op_kwargs)

            # Tell depletion_operator.finalize NOT to clear C API memory when
            # it is done
            depletion_operator.cleanup_when_done = False

            # Set up the integrator
            check_value('method', method,
                        dep.integrators.integrator_by_name.keys())
            integrator_class = dep.integrators.integrator_by_name[method]
            integrator = integrator_class(depletion_operator, timesteps,
                                          **integrator_kwargs)

            # Now perform the depletion
            with openmc.lib.quiet_dll(output):
                integrator.integrate(final_step)

            # Now make the python Materials match the C API material data
            for mat_id, mat in self._materials_by_id.items():
                if mat.depletable:
                    # Get the C data
                    c_mat = openmc.lib.materials[mat_id]
                    nuclides, densities = c_mat._get_densities()
                    # And now we can remove isotopes and add these ones in
                    mat.nuclides.clear()
                    for nuc, density in zip(nuclides, densities):
                        mat.add_nuclide(nuc, density)
                    mat.set_density('atom/b-cm', sum(densities))

            # If we didnt start intialized, we should cleanup after ourselves
            if not started_initialized:
                depletion_operator.cleanup_when_done = True
                depletion_operator.finalize()

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

        .. versionchanged:: 0.13.0
            This method now updates values as loaded in memory with the C API

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
                    if self.is_initialized:
                        lib_cell = openmc.lib.cells[cell_id]
                        lib_cell.set_temperature(group['temperature'][()])

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
                if self.is_initialized:
                    C_mat = openmc.lib.materials[mat_id]
                    C_mat.set_density(atom_density, 'atom/b-cm')

    def run(self, particles=None, threads=None, geometry_debug=False,
            restart_file=None, tracks=False, output=True, cwd='.',
            openmc_exec='openmc', mpi_args=None, event_based=False):
        """Runs OpenMC. If the C API has been initialized, then the C API is
        used, otherwise, this method creates the XML files and runs OpenMC via
        a system call. In both cases this method returns the path to the last
        statepoint file generated.

        .. versionchanged:: 0.12
            Instead of returning the final k-effective value, this function now
            returns the path to the final statepoint written.

        .. versionchanged:: 0.13.0
            This method can utilize the C API for execution

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
        with _change_directory(Path(cwd)):
            if self.is_initialized:
                # Handle the run options as applicable
                # First dont allow ones that must be set via init
                for arg_name, arg, default in zip(
                    ['threads', 'geometry_debug', 'restart_file', 'tracks'],
                    [threads, geometry_debug, restart_file, tracks],
                    [None, False, None, False]):
                    if arg != default:
                        msg = f"{arg_name} must be set via Model.is_initialized(...)"
                        raise ValueError(msg)

                if particles is not None:
                    init_particles = openmc.lib.settings.particles
                    if isinstance(particles, Integral) and particles > 0:
                        openmc.lib.settings.particles = particles

                # Event-based can be set at init-time or on a case-basis.
                # Handle the argument here.
                # TODO This will be dealt with in a future change to the C API

                # Then run using the C API
                openmc.lib.run(output)

                # Reset changes for the openmc.run kwargs handling
                if particles is not None:
                    openmc.lib.settings.particles = init_particles

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

    def calculate_volumes(self, threads=None, output=True, cwd='.',
                          openmc_exec='openmc', mpi_args=None,
                          apply_volumes=True):
        """Runs an OpenMC stochastic volume calculation and, if requested,
        applies volumes to the model

        .. versionadded:: 0.13.0

        Parameters
        ----------
        threads : int, optional
            Number of OpenMP threads. If OpenMC is compiled with OpenMP
            threading enabled, the default is implementation-dependent but is
            usually equal to the number of hardware threads available (or a
            value set by the :envvar:`OMP_NUM_THREADS` environment variable).
            This currenty only applies to the case when not using the C API.
        output : bool, optional
            Capture OpenMC output from standard out
        openmc_exec : str, optional
            Path to OpenMC executable. Defaults to 'openmc'.
            This only applies to the case when not using the C API.
        mpi_args : list of str, optional
            MPI execute command and any additional MPI arguments to pass,
            e.g. ['mpiexec', '-n', '8'].
            This only applies to the case when not using the C API.
        cwd : str, optional
            Path to working directory to run in. Defaults to the current
            working directory.
        apply_volumes : bool, optional
            Whether apply the volume calculation results from this calculation
            to the model. Defaults to applying the volumes.
        """

        if len(self.settings.volume_calculations) == 0:
            # Then there is no volume calculation specified
            raise ValueError("The Settings.volume_calculation attribute must"
                             " be specified before executing this method!")

        with _change_directory(Path(cwd)):
            if self.is_initialized:
                if threads is not None:
                    msg = "Threads must be set via Model.is_initialized(...)"
                    raise ValueError(msg)
                if mpi_args is not None:
                    msg = "The MPI environment must be set otherwise such as" \
                        "with the call to mpi4py"
                    raise ValueError(msg)

                # Compute the volumes
                openmc.lib.calculate_volumes(output)

            else:
                self.export_to_xml()
                openmc.calculate_volumes(threads=threads, output=output,
                                         openmc_exec=openmc_exec,
                                         mpi_args=mpi_args)

            # Now we apply the volumes
            if apply_volumes:
                # Load the results and add them to the model
                for i, vol_calc in enumerate(self.settings.volume_calculations):
                    f_name = f"volume_{i + 1}.h5"
                    vol_calc.load_results(f_name)
                    # First add them to the Python side
                    self.geometry.add_volume_information(vol_calc)

                    # And now repeat for the C API
                    if self.is_initialized and vol_calc.domain_type == 'material':
                        # Then we can do this in the C API
                        for domain_id in vol_calc.ids:
                            openmc.lib.materials[domain_id].volume = \
                                vol_calc.volumes[domain_id].n

    def plot_geometry(self, output=True, cwd='.', openmc_exec='openmc',
                      convert=True, convert_exec='convert'):
        """Creates plot images as specified by the Model.plots attribute

        If convert is True, this function requires that a program is installed
        to convert PPM files to PNG files. Typically, that would be
        `ImageMagick <https://www.imagemagick.org>`_ which includes a
        `convert` command.

        .. versionadded:: 0.13.0

        Parameters
        ----------
        output : bool, optional
            Capture OpenMC output from standard out
        cwd : str, optional
            Path to working directory to run in. Defaults to the current
            working directory.
        openmc_exec : str, optional
            Path to OpenMC executable. Defaults to 'openmc'.
            This only applies to the case when not using the C API.
        convert : bool, optional
            Whether or not to attempt to convert from PPM to PNG
        convert_exec : str, optional
            Command that can convert PPM files into PNG files
        """

        if len(self.plots) == 0:
            # Then there is no volume calculation specified
            raise ValueError("The Model.plots attribute must be specified "
                             "before executing this method!")

        with _change_directory(Path(cwd)):
            # TODO: openmc.is_initialized doesnt read plots.xml unless it is
            # in plot mode so the following will not work. Commented out for
            # now and replacing with non-C API code
            self.export_to_xml()
            openmc.plot_geometry(output=output, openmc_exec=openmc_exec)

            # if self.is_initialized:
            #     # Compute the volumes
            #     openmc.lib.plot_geometry(output)
            # else:
            #     self.export_to_xml()
            #     openmc.plot_geometry(output=output, openmc_exec=openmc_exec)

            if convert:
                for p in self.plots:
                    if p.filename is not None:
                        ppm_file = f'{p.filename}.ppm'
                    else:
                        ppm_file = f'plot_{p.id}.ppm'
                    png_file = ppm_file.replace('.ppm', '.png')
                    subprocess.check_call([convert_exec, ppm_file, png_file])

    def _change_py_lib_attribs(self, names_or_ids, value, obj_type,
                               attrib_name, density_units='atom/b-cm'):
        # Method to do the same work whether it is a cell or material and
        # a temperature or volume
        check_type('names_or_ids', names_or_ids, Iterable, (Integral, str))
        check_type('obj_type', obj_type, str)
        obj_type = obj_type.lower()
        check_value('obj_type', obj_type, ('material', 'cell'))
        check_value('attrib_name', attrib_name,
                    ('temperature', 'volume', 'density', 'rotation',
                     'translation'))
        # The C API only allows setting density units of atom/b-cm and g/cm3
        check_value('density_units', density_units, ('atom/b-cm', 'g/cm3'))
        # The C API has no way to set cell volume so lets raise an exception
        if obj_type == 'cell' and attrib_name == 'volume':
            raise NotImplementedError(
                'Setting a Cell volume is not supported!')
        # Same with setting temperatures, TODO: update C API for this
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
            obj_by_id = openmc.lib.cells
        else:
            by_name = self._materials_by_name
            by_id = self._materials_by_id
            obj_by_id = openmc.lib.materials

        # Get the list of ids to use if converting from names and accepting
        # only values that have actual ids
        ids = []
        for name_or_id in names_or_ids:
            if isinstance(name_or_id, Integral):
                if name_or_id in by_id:
                    ids.append(int(name_or_id))
                else:
                    cap_obj = obj_type.capitalize()
                    msg = f'{cap_obj} ID {name_or_id} " \
                        "is not present in the model!'
                    raise InvalidIDError(msg)
            elif isinstance(name_or_id, str):
                if name_or_id in by_name:
                    # Then by_name[name_or_id] is a list so we need to add all
                    # entries
                    ids.extend([obj.id for obj in by_name[name_or_id]])
                else:
                    cap_obj = obj_type.capitalize()
                    msg = f'{cap_obj} {name_or_id} " \
                        "is not present in the model!'
                    raise InvalidIDError(msg)

        # Now perform the change to both python and C API
        for id_ in ids:
            obj = by_id[id_]
            if attrib_name == 'density':
                obj.set_density(density_units, value)
            else:
                setattr(obj, attrib_name, value)
            # Next lets keep what is in C API memory up to date as well
            if self.is_initialized:
                lib_obj = obj_by_id[id_]
                if attrib_name == 'temperature':
                    lib_obj.set_temperature(value)
                elif attrib_name == 'density':
                    lib_obj.set_density(value, density_units)
                else:
                    setattr(lib_obj, attrib_name, value)

    def rotate_cells(self, names_or_ids, vector):
        """Rotate the identified cell(s) by the specified rotation vector.
        The rotation is only applied to cells filled with a universe.

        .. note:: If applying this change to a name that is not unique, then
        the change will be applied to all objects of that name.

        .. versionadded:: 0.13.0

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The cell names (if str) or id (if int) that are to be translated
            or rotated. This parameter can include a mix of names and ids.
        vector : Iterable of float
            The rotation vector of length 3 to apply. This array specifies the
            angles in degrees about the x, y, and z axes, respectively.

        """

        self._change_py_lib_attribs(names_or_ids, vector, 'cell', 'rotation')

    def translate_cells(self, names_or_ids, vector):
        """Translate the identified cell(s) by the specified translation vector.
        The translation is only applied to cells filled with a universe.

        .. note:: If applying this change to a name that is not unique, then
        the change will be applied to all objects of that name.

        .. versionadded:: 0.13.0

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The cell names (if str) or id (if int) that are to be translated
            or rotated. This parameter can include a mix of names and ids.
        vector : Iterable of float
            The translation vector of length 3 to apply. This array specifies
            the x, y, and z dimensions of the translation.

        """

        self._change_py_lib_attribs(names_or_ids, vector, 'cell', 'translation')

    def update_densities(self, names_or_ids, density, density_units='atom/b-cm'):
        """Update the density of a given set of materials to a new value

        .. note:: If applying this change to a name that is not unique, then
        the change will be applied to all objects of that name.

        .. versionadded:: 0.13.0

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

        self._change_py_lib_attribs(names_or_ids, density, 'material', 'density',
                                  density_units)

    def update_cell_temperatures(self, names_or_ids, temperature):
        """Update the temperature of a set of cells to the given value

        .. note:: If applying this change to a name that is not unique, then
        the change will be applied to all objects of that name.

        .. versionadded:: 0.13.0

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The cell names (if str) or id (if int) that are to be updated.
            This parameter can include a mix of names and ids.
        temperature : float
            The temperature to apply in units of Kelvin

        """

        self._change_py_lib_attribs(names_or_ids, temperature, 'cell',
                                  'temperature')

    def update_material_temperatures(self, names_or_ids, temperature):
        """Update the temperature of a set of materials to the given value

        .. note:: If applying this change to a name that is not unique, then
        the change will be applied to all objects of that name.

        .. versionadded:: 0.13.0

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The material names (if str) or id (if int) that are to be updated.
            This parameter can include a mix of names and ids.
        temperature : float
            The temperature to apply in units of Kelvin

        """

        self._change_py_lib_attribs(names_or_ids, temperature, 'material',
                                  'temperature')

    def update_material_volumes(self, names_or_ids, volume):
        """Update the volume of a set of materials to the given value

        .. note:: If applying this change to a name that is not unique, then
        the change will be applied to all objects of that name.

        .. versionadded:: 0.13.0

        Parameters
        ----------
        names_or_ids : Iterable of str or int
            The material names (if str) or id (if int) that are to be updated.
            This parameter can include a mix of names and ids.
        volume : float
            The volume to apply in units of cm^3

        """

        self._change_py_lib_attribs(names_or_ids, volume, 'material', 'volume')
