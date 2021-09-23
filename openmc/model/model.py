from collections.abc import Iterable
from pathlib import Path
import time
import warnings

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
        self.depletion_operator = None

        if self.materials is None:
            mats = self.geometry.get_all_materials().values()
        else:
            mats = self.materials
        self._materials_by_id = {mat.id: mat for mat in mats}
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
    def depletion_operator(self):
        return self._depletion_operator

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
            self._chain_file = Path(chain_file)
        else:
            self._chain_file = chain_file

    @fission_q.setter
    def fission_q(self, fission_q):
        check_type('fission_q', fission_q, (type(None), dict))
        self._fission_q = fission_q

    @depletion_operator.setter
    def depletion_operator(self, depletion_operator):
        check_type('depletion_operator', depletion_operator,
                   (type(None), dep.Operator))
        self._depletion_operator = depletion_operator

    @classmethod
    def from_xml(cls, geometry='geometry.xml', materials='materials.xml',
                 settings='settings.xml'):
        """Create model from existing XML files

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

    def init_C_api(self, use_depletion_operator=False):
        """Initializes the model in memory via the C-API

        Parameters
        ----------
        directory : str
            Directory to write XML files to. If it doesn't exist already, it
            will be created.
        use_depletion_operator : bool, optional
            If True, the model will be loaded using the depletion operator
            including all isotopes necessary from fission. This parameter will
            use the :attr:`Model.chain_file` and :attr:`Model.fission_q`
            attributes.
        """

        if use_depletion_operator:
            # Create OpenMC transport operator
            self.depletion_operator = \
                dep.Operator(self.geometry, self.settings,
                             str(self.chain_file), fission_q=self.fission_q)
        else:
            openmc.lib.hard_reset()
        if dep.comm.rank == 0:
            self.export_to_xml()
        dep.comm.barrier()
        openmc.lib.init(intracomm=dep.comm)

    def clear_C_api(self):
        """Finalize simulation and free memory allocated for the C-API"""
        openmc.lib.finalize()

    def deplete(self, timesteps, chain_file=None, method='cecm',
                fission_q=None, final_step=True, **kwargs):
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
        method : str
             Integration method used for depletion (e.g., 'cecm', 'predictor')
        fission_q : dict, optional
            Dictionary of nuclides and their fission Q values [eV].
            If not given, values will be pulled from the ``chain_file``.
        final_step : bool, optional
            Indicate whether or not a transport solve should be run at the end
            of the last timestep.

            .. versionadded:: 0.12.3
        **kwargs
            Keyword arguments passed to integration function (e.g.,
            :func:`openmc.deplete.integrator.cecm`)

        """

        if self.C_init and self.depletion_operator is not None:
            # Then the user has properly initialized the information and we can
            # just carry forward
            pass
        elif self.C_init and self.depletion_operator is None:
            # Then the user has initialzed the C-API but without the depletion
            # isotopes loaded. We would have to reset and reload data, but
            # doing so could lose user information. Therefore let us just
            # provide an error and quit.
            msg = "Model.deplete(...) cannot be called after " \
                "Model.init_C_api(...) if the use_depletion_operator " \
                "argument to Model.init_C_api(...) is False."
            raise SetupError(msg)
        else:
            # To get here, the C-API is not initialized. So we can do that now
            # To keep Model.deplete(...) API compatibility, we will allow the
            # chain_file and fission_q params to be set since we havent loaded
            # the API anyways
            if chain_file is not None:
                self.chain_file = chain_file
                warnings.warn("The chain_file argument of Model.deplete(...) "
                              "has been deprecated and may be removed in a "
                              "future version. The Model.chain_file should be"
                              "used instead.", DeprecationWarning)
            if fission_q is not None:
                warnings.warn("The fission_q argument of Model.deplete(...) "
                              "has been deprecated and may be removed in a "
                              "future version. The Model.fission_q should be"
                              "used instead.", DeprecationWarning)
                self.fission_q = fission_q
            self.init_C_api(use_depletion_operator=True)

        # Set up the integrator
        integrator_class = dep.integrators.integrator_factory(method)
        integrator = integrator_class(self.depletion_operator,
                                      timesteps, **kwargs)

        # Now perform the depletion
        integrator.integrate(final_step)

        # If we did not perform a transport calculation on the final step, then
        # make the code update the C-API material inventory
        if not final_step:
            self.depletion_operator._update_materials()

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

            # Make sure number of materials matches
            mats_group = fh['materials']
            n_cells = mats_group.attrs['n_materials']
            if n_cells != len(materials):
                raise ValueError("Number of materials in properties file doesn't "
                                 "match current model.")

            # Update material densities
            for name, group in mats_group.items():
                mat_id = int(name.split()[1])
                atom_density = group.attrs['atom_density']
                materials[mat_id].set_density('atom/b-cm', atom_density)

    def run(self, **kwargs):
        """Runs OpenMC. If the C-API has been initialized, then the C-API is
        used, otherwise, this method creates the XML files and runs OpenMC via
        a system cal. In both cases this method returns the path to the last
        statepoint file generated.

        .. versionchanged:: 0.12
            Instead of returning the final k-effective value, this function now
            returns the path to the final statepoint written.

        .. versionchanged:: 0.12.3
            This method can utilize the C-API for execution

        Parameters
        ----------
        **kwargs
            Keyword arguments passed to :func:`openmc.run`

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

        if self.C_init:
            # Then run using the C-API
            openmc.lib.run()
        else:
            # Then run via the command line
            self.export_to_xml()
            openmc.run(**kwargs)

        # Get output directory and return the last statepoint written this run
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

    def _move_cell(self, cell_names_or_ids, vector, attrib_name):
        # Method to do the same work whether it is a rotation or translation
        check_type('cell_names_or_ids', cell_names_or_ids, Iterable,
                   (np.int, int, str))
        check_type('vector', vector, Iterable, (np.float, float))
        check_length('vector', vector, 3)
        check_value('attrib_name', attrib_name, ('rotation', 'translation'))

        # Get the list of cell ids to use y converting from names and accepting
        # only values that have actual ids
        cell_ids = [None] * len(cell_names_or_ids)
        for c, cell_name_or_id in enumerate(cell_names_or_ids):
            if isinstance(cell_name_or_id, (int, np.int)):
                if cell_name_or_id in self._cells_by_id:
                    cell_ids[c] = int(cell_name_or_id)
                msg = 'Cell ID {} is not present in the model!'.format(
                    cell_name_or_id)
                raise InvalidIDError(msg)
            elif isinstance(cell_name_or_id, str):
                if cell_name_or_id in self._cells_by_name:
                    cell_ids[c] = self._cells_by_name[cell_name_or_id]
                else:
                    msg = 'Cell {} is not present in the model!'.format(
                        cell_name_or_id)
                    raise InvalidIDError(msg)

        # Now perform the motion
        for cell_id in cell_ids:
            cell = self._cells_by_id[cell_id]
            if attrib_name == 'rotation':
                cell.rotation = vector
            elif attrib_name == 'translation':
                cell.translation = vector
            # Next lets keep what is in C-API memory up to date as well
            if self.C_init:
                C_cell = openmc.lib.cells[cell_id]
                if attrib_name == 'rotation':
                    C_cell.rotation = vector
                elif attrib_name == 'translation':
                    C_cell.translation = vector

    def rotate_cells(self, cell_names_or_ids, vector):
        """Rotate the identified cell(s) by the specified rotation vector.
        The rotation is only applied to cells filled with a universe.

        Parameters
        ----------
        cell_names_or_ids : Iterable of str or int
            The cell names (if str) or id (if int) that are to be translated
            or rotated. This parameter can include a mix of names and ids.
        vector : Iterable of float
            The rotation vector of length 3 to apply. This array specifies the
            angles in degrees about the x, y, and z axes, respectively.

        """

        self._move_cell(cell_names_or_ids, vector, 'rotation')

    def translate_cells(self, cell_names_or_ids, vector):
        """Translate the identified cell(s) by the specified translation vector.
        The translation is only applied to cells filled with a universe.

        Parameters
        ----------
        cell_names_or_ids : Iterable of str or int
            The cell names (if str) or id (if int) that are to be translated
            or rotated. This parameter can include a mix of names and ids.
        vector : Iterable of float
            The translation vector of length 3 to apply. This array specifies
            the x, y, and z dimensions of the translation.

        """

        self._move_cell(cell_names_or_ids, vector, 'translation')
