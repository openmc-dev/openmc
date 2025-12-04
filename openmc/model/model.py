from __future__ import annotations
from collections.abc import Callable, Iterable, Sequence
import copy
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path
import math
from numbers import Integral, Real
import random
import re
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Any, Protocol
import warnings

import h5py
import lxml.etree as ET
import numpy as np
from scipy.optimize import curve_fit

import openmc
import openmc._xml as xml
from openmc.dummy_comm import DummyCommunicator
from openmc.executor import _process_CLI_arguments
from openmc.checkvalue import check_type, check_value, PathLike
from openmc.exceptions import InvalidIDError
from openmc.plots import add_plot_params, _BASIS_INDICES
from openmc.utility_funcs import change_directory


# Protocol for a function that is passed to search_keff
class ModelModifier(Protocol):
    def __call__(self, val: float, **kwargs: Any) -> None:
        ...


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
        The model information can now be loaded in to OpenMC directly via
        openmc.lib

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

    """

    def __init__(
        self,
        geometry: openmc.Geometry | None = None,
        materials: openmc.Materials | None = None,
        settings: openmc.Settings | None = None,
        tallies: openmc.Tallies | None = None,
        plots: openmc.Plots | None = None,
    ):
        self.geometry = openmc.Geometry() if geometry is None else geometry
        self.materials = openmc.Materials() if materials is None else materials
        self.settings = openmc.Settings() if settings is None else settings
        self.tallies = openmc.Tallies() if tallies is None else tallies
        self.plots = openmc.Plots() if plots is None else plots

    @property
    def geometry(self) -> openmc.Geometry:
        return self._geometry

    @geometry.setter
    def geometry(self, geometry):
        check_type('geometry', geometry, openmc.Geometry)
        self._geometry = geometry

    @property
    def materials(self) -> openmc.Materials:
        return self._materials

    @materials.setter
    def materials(self, materials):
        check_type('materials', materials, Iterable, openmc.Material)
        if isinstance(materials, openmc.Materials):
            self._materials = materials
        else:
            if not hasattr(self, '_materials'):
                self._materials = openmc.Materials()
            del self._materials[:]
            for mat in materials:
                self._materials.append(mat)

    @property
    def settings(self) -> openmc.Settings:
        return self._settings

    @settings.setter
    def settings(self, settings):
        check_type('settings', settings, openmc.Settings)
        self._settings = settings

    @property
    def tallies(self) -> openmc.Tallies:
        return self._tallies

    @tallies.setter
    def tallies(self, tallies):
        check_type('tallies', tallies, Iterable, openmc.Tally)
        if isinstance(tallies, openmc.Tallies):
            self._tallies = tallies
        else:
            if not hasattr(self, '_tallies'):
                self._tallies = openmc.Tallies()
            del self._tallies[:]
            for tally in tallies:
                self._tallies.append(tally)

    @property
    def plots(self) -> openmc.Plots:
        return self._plots

    @plots.setter
    def plots(self, plots):
        check_type('plots', plots, Iterable, openmc.PlotBase)
        if isinstance(plots, openmc.Plots):
            self._plots = plots
        else:
            if not hasattr(self, '_plots'):
                self._plots = openmc.Plots()
            del self._plots[:]
            for plot in plots:
                self._plots.append(plot)

    @property
    def bounding_box(self) -> openmc.BoundingBox:
        return self.geometry.bounding_box

    @property
    def is_initialized(self) -> bool:
        try:
            import openmc.lib
            return openmc.lib.is_initialized
        except ImportError:
            return False

    @property
    @cache
    def _materials_by_id(self) -> dict:
        """Dictionary mapping material ID --> material"""
        if self.materials:
            mats = self.materials
        else:
            mats = self.geometry.get_all_materials().values()
        return {mat.id: mat for mat in mats}

    @property
    @cache
    def _cells_by_id(self) -> dict:
        """Dictionary mapping cell ID --> cell"""
        cells = self.geometry.get_all_cells()
        return {cell.id: cell for cell in cells.values()}

    @property
    @cache
    def _cells_by_name(self) -> dict[int, openmc.Cell]:
        # Get the names maps, but since names are not unique, store a set for
        # each name key. In this way when the user requests a change by a name,
        # the change will be applied to all of the same name.
        result = {}
        for cell in self.geometry.get_all_cells().values():
            if cell.name not in result:
                result[cell.name] = set()
            result[cell.name].add(cell)
        return result

    @property
    @cache
    def _materials_by_name(self) -> dict[int, openmc.Material]:
        if self.materials is None:
            mats = self.geometry.get_all_materials().values()
        else:
            mats = self.materials
        result = {}
        for mat in mats:
            if mat.name not in result:
                result[mat.name] = set()
            result[mat.name].add(mat)
        return result

    # TODO: This should really get incorporated in lower-level calls to
    # get_all_materials, but right now it requires information from the Model object
    def _get_all_materials(self) -> dict[int, openmc.Material]:
        """Get all materials including those in DAGMC universes

        Returns
        -------
        dict
            Dictionary mapping material ID to material instances
        """
        # Get all materials from the Geometry object
        materials = self.geometry.get_all_materials()

        # Account for materials in DAGMC universes
        for cell in self.geometry.get_all_cells().values():
            if isinstance(cell.fill, openmc.DAGMCUniverse):
                names = cell.fill.material_names
                materials.update({
                    mat.id: mat for mat in self.materials if mat.name in names
                })

        return materials

    def add_kinetics_parameters_tallies(self, num_groups: int | None = None, nuclides: Sequence[str] | None = None):
        """Add tallies for calculating kinetics parameters using the IFP method.

        This method adds tallies to the model for calculating two kinetics
        parameters, the generation time and the effective delayed neutron
        fraction (beta effective). After a model is run, these parameters can be
        determined through the :meth:`openmc.StatePoint.get_kinetics_parameters` method.

        Parameters
        ----------
        num_groups : int, optional
            Number of precursor groups to filter the delayed neutron fraction.
            If None, only the total effective delayed neutron fraction is
            tallied.
        nuclides : int, optional 
            Nuclides to calculate separate kinetic parameters for.
            If None, do not separate kinetic parameters per nuclide.
        """
        if not any('ifp-time-numerator' in t.scores for t in self.tallies):
            gen_time_tally = openmc.Tally(name='IFP time numerator')
            gen_time_tally.scores = ['ifp-time-numerator']
            self.tallies.append(gen_time_tally)
        if not any('ifp-beta-numerator' in t.scores for t in self.tallies):
            beta_tally = openmc.Tally(name='IFP beta numerator')
            beta_tally.scores = ['ifp-beta-numerator']
            if num_groups is not None:
                beta_tally.filters = [openmc.DelayedGroupFilter(list(range(1, num_groups + 1)))]
            if nuclides is not None:
                beta_tally.nuclides = list(nuclides)
            self.tallies.append(beta_tally)
        if not any('ifp-denominator' in t.scores for t in self.tallies):
            denom_tally = openmc.Tally(name='IFP denominator')
            denom_tally.scores = ['ifp-denominator']
            self.tallies.append(denom_tally)

    @classmethod
    def from_xml(
        cls,
        geometry: PathLike = "geometry.xml",
        materials: PathLike = "materials.xml",
        settings: PathLike = "settings.xml",
        tallies: PathLike = "tallies.xml",
        plots: PathLike = "plots.xml",
    ) -> Model:
        """Create model from existing XML files

        Parameters
        ----------
        geometry : PathLike
            Path to geometry.xml file
        materials : PathLike
            Path to materials.xml file
        settings : PathLike
            Path to settings.xml file
        tallies : PathLike
            Path to tallies.xml file

            .. versionadded:: 0.13.0
        plots : PathLike
            Path to plots.xml file

            .. versionadded:: 0.13.0

        Returns
        -------
        openmc.model.Model
            Model created from XML files

        """
        materials = openmc.Materials.from_xml(materials)
        geometry = openmc.Geometry.from_xml(geometry, materials)
        settings = openmc.Settings.from_xml(settings)
        tallies = openmc.Tallies.from_xml(
            tallies) if Path(tallies).exists() else None
        plots = openmc.Plots.from_xml(plots) if Path(plots).exists() else None
        return cls(geometry, materials, settings, tallies, plots)

    @classmethod
    def from_model_xml(cls, path: PathLike = "model.xml") -> Model:
        """Create model from single XML file

        .. versionadded:: 0.13.3

        Parameters
        ----------
        path : PathLike
            Path to model.xml file
        """
        parser = ET.XMLParser(huge_tree=True)
        tree = ET.parse(path, parser=parser)
        root = tree.getroot()

        model = cls()

        meshes = {}
        model.settings = openmc.Settings.from_xml_element(
            root.find('settings'), meshes)
        model.materials = openmc.Materials.from_xml_element(
            root.find('materials'))
        model.geometry = openmc.Geometry.from_xml_element(
            root.find('geometry'), model.materials)

        if root.find('tallies') is not None:
            model.tallies = openmc.Tallies.from_xml_element(
                root.find('tallies'), meshes)

        if root.find('plots') is not None:
            model.plots = openmc.Plots.from_xml_element(root.find('plots'))

        return model

    def init_lib(
        self,
        threads: int | None = None,
        geometry_debug: bool = False,
        restart_file: PathLike | None = None,
        tracks: bool = False,
        output: bool = True,
        event_based: bool | None = None,
        intracomm=None,
        directory: PathLike | None = None,
    ):
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
        restart_file : PathLike, optional
            Path to restart file to use
        tracks : bool, optional
            Enables the writing of particles tracks. The number of particle
            tracks written to tracks.h5 is limited to 1000 unless
            Settings.max_tracks is set. Defaults to False.
        output : bool
            Capture OpenMC output from standard out
        event_based : None or bool, optional
            Turns on event-based parallelism if True. If None, the value in
            the Settings will be used.
        intracomm : mpi4py.MPI.Intracomm or None, optional
            MPI intracommunicator
        directory : PathLike or None, optional
            Directory to write XML files to. Defaults to None.
        """

        import openmc.lib

        # TODO: right now the only way to set most of the above parameters via
        # the C API are at initialization time despite use-cases existing to
        # set them for individual runs. For now this functionality is exposed
        # where it exists (here in init), but in the future the functionality
        # should be exposed so that it can be accessed via model.run(...)

        args = _process_CLI_arguments(
            volume=False, geometry_debug=geometry_debug,
            restart_file=restart_file, threads=threads, tracks=tracks,
            event_based=event_based, path_input=directory)

        # Args adds the openmc_exec command in the first entry; remove it
        args = args[1:]

        self.finalize_lib()

        # The Model object needs to be aware of the communicator so it can
        # use it in certain cases, therefore lets store the communicator
        if intracomm is not None:
            self._intracomm = intracomm
        else:
            self._intracomm = DummyCommunicator()

        if self._intracomm.rank == 0:
            if directory is not None:
                self.export_to_xml(directory=directory)
            else:
                self.export_to_xml()
        self._intracomm.barrier()

        # We cannot pass DummyCommunicator to openmc.lib.init so pass instead
        # the user-provided intracomm which will either be None or an mpi4py
        # communicator
        openmc.lib.init(args=args, intracomm=intracomm, output=output)

    def sync_dagmc_universes(self):
        """Synchronize all DAGMC universes in the current geometry.

        This method iterates over all DAGMC universes in the geometry and
        synchronizes their cells with the current material assignments. Requires
        that the model has been initialized via :meth:`Model.init_lib`.

        .. versionadded:: 0.15.1

        """
        if self.is_initialized:
            if self.materials:
                materials = self.materials
            else:
                materials = list(self.geometry.get_all_materials().values())
            for univ in self.geometry.get_all_universes().values():
                if isinstance(univ, openmc.DAGMCUniverse):
                    univ.sync_dagmc_cells(materials)
        else:
            raise ValueError("The model must be initialized before calling "
                             "this method")

    def finalize_lib(self):
        """Finalize simulation and free memory allocated for the C API

        .. versionadded:: 0.13.0

        """

        import openmc.lib

        openmc.lib.finalize()

    def deplete(
        self,
        method: str = "cecm",
        final_step: bool = True,
        operator_kwargs: dict | None = None,
        directory: PathLike = ".",
        output: bool = True,
        **integrator_kwargs,
    ):
        """Deplete model using specified timesteps/power

        .. versionchanged:: 0.13.0
            The *final_step*, *operator_kwargs*, *directory*, and *output*
            arguments were added.

        Parameters
        ----------
        timesteps : iterable of float or iterable of tuple
            Array of timesteps. Note that values are not cumulative. The units are
            specified by the `timestep_units` argument when `timesteps` is an
            iterable of float. Alternatively, units can be specified for each step
            by passing an iterable of (value, unit) tuples.
        method : str
             Integration method used for depletion (e.g., 'cecm', 'predictor').
             Defaults to 'cecm'.
        final_step : bool, optional
            Indicate whether or not a transport solve should be run at the end
            of the last timestep. Defaults to running this transport solve.
        operator_kwargs : dict
            Keyword arguments passed to the depletion operator initializer
            (e.g., :func:`openmc.deplete.Operator`)
        directory : PathLike, optional
            Directory to write XML files to. If it doesn't exist already, it
            will be created. Defaults to the current working directory
        output : bool
            Capture OpenMC output from standard out
        integrator_kwargs : dict
            Remaining keyword arguments passed to the depletion integrator
            (e.g., :class:`openmc.deplete.CECMIntegrator`).

        """

        if operator_kwargs is None:
            op_kwargs = {}
        elif isinstance(operator_kwargs, dict):
            op_kwargs = operator_kwargs
        else:
            raise ValueError("operator_kwargs must be a dict or None")

        # Import openmc.deplete here so the Model can be used even if the
        # shared library is unavailable.
        import openmc.deplete as dep

        # Store whether or not the library was initialized when we started
        started_initialized = self.is_initialized

        with change_directory(directory):
            with openmc.lib.quiet_dll(output):
                # TODO: Support use of IndependentOperator too
                depletion_operator = dep.CoupledOperator(self, **op_kwargs)

            # Tell depletion_operator.finalize NOT to clear C API memory when
            # it is done
            depletion_operator.cleanup_when_done = False

            # Set up the integrator
            check_value('method', method,
                        dep.integrators.integrator_by_name.keys())
            integrator_class = dep.integrators.integrator_by_name[method]
            integrator = integrator_class(depletion_operator, **integrator_kwargs)

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

    def export_to_xml(self, directory: PathLike = '.', remove_surfs: bool = False,
                      nuclides_to_ignore: Iterable[str] | None = None):
        """Export model to separate XML files.

        Parameters
        ----------
        directory : PathLike
            Directory to write XML files to. If it doesn't exist already, it
            will be created.
        remove_surfs : bool
            Whether or not to remove redundant surfaces from the geometry when
            exporting.

            .. versionadded:: 0.13.1
        nuclides_to_ignore : list of str
            Nuclides to ignore when exporting to XML.

        """
        # Create directory if required
        d = Path(directory)
        if not d.is_dir():
            d.mkdir(parents=True, exist_ok=True)

        self.settings.export_to_xml(d)
        self.geometry.export_to_xml(d, remove_surfs=remove_surfs)

        # If a materials collection was specified, export it. Otherwise, look
        # for all materials in the geometry and use that to automatically build
        # a collection.
        if self.materials:
            self.materials.export_to_xml(d, nuclides_to_ignore=nuclides_to_ignore)
        else:
            materials = openmc.Materials(self.geometry.get_all_materials()
                                         .values())
            materials.export_to_xml(d, nuclides_to_ignore=nuclides_to_ignore)

        if self.tallies:
            self.tallies.export_to_xml(d)
        if self.plots:
            self.plots.export_to_xml(d)

    def export_to_model_xml(self, path: PathLike = 'model.xml', remove_surfs: bool = False,
                            nuclides_to_ignore: Iterable[str] | None = None):
        """Export model to a single XML file.

        .. versionadded:: 0.13.3

        Parameters
        ----------
        path : str or PathLike
            Location of the XML file to write (default is 'model.xml'). Can be a
            directory or file path.
        remove_surfs : bool
            Whether or not to remove redundant surfaces from the geometry when
            exporting.
        nuclides_to_ignore : list of str
            Nuclides to ignore when exporting to XML.

        """
        xml_path = Path(path)
        # if the provided path doesn't end with the XML extension, assume the
        # input path is meant to be a directory. If the directory does not
        # exist, create it and place a 'model.xml' file there.
        if not str(xml_path).endswith('.xml'):
            if not xml_path.exists():
                xml_path.mkdir(parents=True, exist_ok=True)
            elif not xml_path.is_dir():
                raise FileExistsError(f"File exists and is not a directory: '{xml_path}'")
            xml_path /= 'model.xml'
        # if this is an XML file location and the file's parent directory does
        # not exist, create it before continuing
        elif not xml_path.parent.exists():
            xml_path.parent.mkdir(parents=True, exist_ok=True)

        if remove_surfs:
            warnings.warn("remove_surfs kwarg will be deprecated soon, please "
                          "set the Geometry.merge_surfaces attribute instead.")
            self.geometry.merge_surfaces = True

        # provide a memo to track which meshes have been written
        mesh_memo = set()
        settings_element = self.settings.to_xml_element(mesh_memo)
        geometry_element = self.geometry.to_xml_element()

        xml.clean_indentation(geometry_element, level=1)
        xml.clean_indentation(settings_element, level=1)

        # If a materials collection was specified, export it. Otherwise, look
        # for all materials in the geometry and use that to automatically build
        # a collection.
        if self.materials:
            materials = self.materials
        else:
            materials = openmc.Materials(self.geometry.get_all_materials()
                                         .values())

        with open(xml_path, 'w', encoding='utf-8', errors='xmlcharrefreplace') as fh:
            # write the XML header
            fh.write("<?xml version='1.0' encoding='utf-8'?>\n")
            fh.write("<model>\n")
            # Write the materials collection to the open XML file first.
            # This will write the XML header also
            materials._write_xml(fh, False, level=1,
                                 nuclides_to_ignore=nuclides_to_ignore)
            # Write remaining elements as a tree
            fh.write(ET.tostring(geometry_element, encoding="unicode"))
            fh.write(ET.tostring(settings_element, encoding="unicode"))

            if self.tallies:
                tallies_element = self.tallies.to_xml_element(mesh_memo)
                xml.clean_indentation(
                    tallies_element, level=1, trailing_indent=self.plots)
                fh.write(ET.tostring(tallies_element, encoding="unicode"))
            if self.plots:
                plots_element = self.plots.to_xml_element()
                xml.clean_indentation(
                    plots_element, level=1, trailing_indent=False)
                fh.write(ET.tostring(plots_element, encoding="unicode"))
            fh.write("</model>\n")

    def import_properties(self, filename: PathLike):
        """Import physical properties

        .. versionchanged:: 0.13.0
            This method now updates values as loaded in memory with the C API

        Parameters
        ----------
        filename : PathLike
            Path to properties HDF5 file

        See Also
        --------
        openmc.lib.export_properties

        """
        import openmc.lib

        cells = self.geometry.get_all_cells()
        materials = self.geometry.get_all_materials()

        with h5py.File(filename, 'r') as fh:
            cells_group = fh['geometry/cells']

            # Make sure number of cells matches
            n_cells = fh['geometry'].attrs['n_cells']
            if n_cells != len(cells):
                raise ValueError("Number of cells in properties file doesn't "
                                 "match current model.")

            # Update temperatures and densities for cells filled with materials
            for name, group in cells_group.items():
                cell_id = int(name.split()[1])
                cell = cells[cell_id]
                if cell.fill_type in ('material', 'distribmat'):
                    temperature = group['temperature'][()]
                    cell.temperature = temperature
                    if self.is_initialized:
                        lib_cell = openmc.lib.cells[cell_id]
                        if temperature.size > 1:
                            for i, T in enumerate(temperature):
                                lib_cell.set_temperature(T, i)
                        else:
                            lib_cell.set_temperature(temperature[0])

                    if group['density']:
                      density = group['density'][()]
                      if density.size > 1:
                          cell.density = [rho for rho in density]
                      else:
                          cell.density = density
                      if self.is_initialized:
                          lib_cell = openmc.lib.cells[cell_id]
                          if density.size > 1:
                              for i, rho in enumerate(density):
                                  lib_cell.set_density(rho, i)
                          else:
                              lib_cell.set_density(density[0])

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

    def run(
        self,
        particles: int | None = None,
        threads: int | None = None,
        geometry_debug: bool = False,
        restart_file: PathLike | None = None,
        tracks: bool = False,
        output: bool = True,
        cwd: PathLike = ".",
        openmc_exec: PathLike = "openmc",
        mpi_args: Iterable[str] = None,
        event_based: bool | None = None,
        export_model_xml: bool = True,
        apply_tally_results: bool = False,
        **export_kwargs,
    ) -> Path:
        """Run OpenMC

        If the C API has been initialized, then the C API is used, otherwise,
        this method creates the XML files and runs OpenMC via a system call. In
        both cases this method returns the path to the last statepoint file
        generated.

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
        restart_file : str or PathLike
            Path to restart file to use
        tracks : bool, optional
            Enables the writing of particles tracks. The number of particle
            tracks written to tracks.h5 is limited to 1000 unless
            Settings.max_tracks is set. Defaults to False.
        output : bool, optional
            Capture OpenMC output from standard out
        cwd : PathLike, optional
            Path to working directory to run in. Defaults to the current working
            directory.
        openmc_exec : str, optional
            Path to OpenMC executable. Defaults to 'openmc'.
        mpi_args : list of str, optional
            MPI execute command and any additional MPI arguments to pass, e.g.
            ['mpiexec', '-n', '8'].
        event_based : None or bool, optional
            Turns on event-based parallelism if True. If None, the value in the
            Settings will be used.
        export_model_xml : bool, optional
            Exports a single model.xml file rather than separate files. Defaults
            to True.

            .. versionadded:: 0.13.3
        apply_tally_results : bool
            Whether to apply results of the final statepoint file to the
            model's tally objects.

            .. versionadded:: 0.15.1
        **export_kwargs
            Keyword arguments passed to either :meth:`Model.export_to_model_xml`
            or :meth:`Model.export_to_xml`.

        Returns
        -------
        Path
            Path to the last statepoint written by this run (None if no
            statepoint was written)

        """

        # Setting tstart here ensures we don't pick up any pre-existing
        # statepoint files in the output directory -- just in case there are
        # differences between the system clock and the filesystem, we get the
        # time of a just-created temporary file
        with NamedTemporaryFile() as fp:
            tstart = Path(fp.name).stat().st_mtime
        last_statepoint = None

        # Operate in the provided working directory
        with change_directory(cwd):
            if self.is_initialized:
                # Handle the run options as applicable
                # First dont allow ones that must be set via init
                for arg_name, arg, default in zip(
                    ['threads', 'geometry_debug', 'restart_file', 'tracks'],
                    [threads, geometry_debug, restart_file, tracks],
                    [None, False, None, False]
                ):
                    if arg != default:
                        msg = f"{arg_name} must be set via Model.is_initialized(...)"
                        raise ValueError(msg)

                init_particles = openmc.lib.settings.particles
                if particles is not None:
                    if isinstance(particles, Integral) and particles > 0:
                        openmc.lib.settings.particles = particles

                init_event_based = openmc.lib.settings.event_based
                if event_based is not None:
                    openmc.lib.settings.event_based = event_based

                # Then run using the C API
                openmc.lib.run(output)

                # Reset changes for the openmc.run kwargs handling
                openmc.lib.settings.particles = init_particles
                openmc.lib.settings.event_based = init_event_based

            else:
                # Then run via the command line
                if export_model_xml:
                    self.export_to_model_xml(**export_kwargs)
                else:
                    self.export_to_xml(**export_kwargs)
                path_input = export_kwargs.get("path", None)
                openmc.run(particles, threads, geometry_debug, restart_file,
                           tracks, output, Path('.'), openmc_exec, mpi_args,
                           event_based, path_input)

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

        if apply_tally_results:
            self.apply_tally_results(last_statepoint)

        return last_statepoint

    def calculate_volumes(
        self,
        threads: int | None = None,
        output: bool = True,
        cwd: PathLike = ".",
        openmc_exec: PathLike = "openmc",
        mpi_args: list[str] | None = None,
        apply_volumes: bool = True,
        export_model_xml: bool = True,
        **export_kwargs,
    ):
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
        export_model_xml : bool, optional
            Exports a single model.xml file rather than separate files. Defaults
            to True.
        **export_kwargs
            Keyword arguments passed to either :meth:`Model.export_to_model_xml`
            or :meth:`Model.export_to_xml`.

        """

        if len(self.settings.volume_calculations) == 0:
            # Then there is no volume calculation specified
            raise ValueError("The Settings.volume_calculations attribute must"
                             " be specified before executing this method!")

        with change_directory(cwd):
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
                if export_model_xml:
                    self.export_to_model_xml(**export_kwargs)
                else:
                    self.export_to_xml(**export_kwargs)
                path_input = export_kwargs.get("path", None)
                openmc.calculate_volumes(
                    threads=threads, output=output, openmc_exec=openmc_exec,
                    mpi_args=mpi_args, path_input=path_input
                )

            # Now we apply the volumes
            if apply_volumes:
                # Load the results and add them to the model
                for i, vol_calc in enumerate(self.settings.volume_calculations):
                    vol_calc.load_results(f"volume_{i + 1}.h5")
                    # First add them to the Python side
                    if vol_calc.domain_type == "material" and self.materials:
                        for material in self.materials:
                            if material.id in vol_calc.volumes:
                                material.add_volume_information(vol_calc)
                    else:
                        self.geometry.add_volume_information(vol_calc)

                    # And now repeat for the C API
                    if self.is_initialized and vol_calc.domain_type == 'material':
                        # Then we can do this in the C API
                        for domain_id in vol_calc.ids:
                            openmc.lib.materials[domain_id].volume = \
                                vol_calc.volumes[domain_id].n


    def _set_plot_defaults(
        self,
        origin: Sequence[float] | None,
        width: Sequence[float] | None,
        pixels: int | Sequence[int],
        basis: str
    ):
        x, y, _ = _BASIS_INDICES[basis]

        bb = self.bounding_box
        # checks to see if bounding box contains -inf or inf values
        if np.isinf(bb.extent[basis]).any():
            if origin is None:
                origin = (0, 0, 0)
            if width is None:
                width = (10, 10)
        else:
            if origin is None:
                # if nan values in the bb.center they get replaced with 0.0
                # this happens when the bounding_box contains inf values
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    origin = np.nan_to_num(bb.center)
            if width is None:
                bb_width = bb.width
                width = (bb_width[x], bb_width[y])

        if isinstance(pixels, int):
            aspect_ratio = width[0] / width[1]
            pixels_y = math.sqrt(pixels / aspect_ratio)
            pixels = (int(pixels / pixels_y), int(pixels_y))

        return origin, width, pixels

    def id_map(
        self,
        origin: Sequence[float] | None = None,
        width: Sequence[float] | None = None,
        pixels: int | Sequence[int] = 40000,
        basis: str = 'xy',
        **init_kwargs
    ) -> np.ndarray:
        """Generate an ID map for domains based on the plot parameters

        If the model is not yet initialized, it will be initialized with
        openmc.lib. If the model is initialized, the model will remain
        initialized after this method call exits.

        .. versionadded:: 0.15.3

        Parameters
        ----------
        origin : Sequence[float], optional
            Origin of the plot. If unspecified, this argument defaults to the
            center of the bounding box if the bounding box does not contain inf
            values for the provided basis, otherwise (0.0, 0.0, 0.0).
        width : Sequence[float], optional
            Width of the plot. If unspecified, this argument defaults to the
            width of the bounding box if the bounding box does not contain inf
            values for the provided basis, otherwise (10.0, 10.0).
        pixels : int | Sequence[int], optional
            If an iterable of ints is provided then this directly sets the
            number of pixels to use in each basis direction. If a single int is
            provided then this sets the total number of pixels in the plot and
            the number of pixels in each basis direction is calculated from this
            total and the image aspect ratio based on the width argument.
        basis : {'xy', 'yz', 'xz'}, optional
            Basis of the plot.
        **init_kwargs
            Keyword arguments passed to :meth:`Model.init_lib`.

        Returns
        -------
        id_map : numpy.ndarray
            A NumPy array with shape (vertical pixels, horizontal pixels, 3) of
            OpenMC property IDs with dtype int32. The last dimension of the
            array contains cell IDs, cell instances, and material IDs (in that
            order).
        """
        import openmc.lib

        origin, width, pixels = self._set_plot_defaults(
            origin, width, pixels, basis)

        # initialize the openmc.lib.plot._PlotBase object
        plot_obj = openmc.lib.plot._PlotBase()
        plot_obj.origin = origin
        plot_obj.width = width[0]
        plot_obj.height = width[1]
        plot_obj.h_res = pixels[0]
        plot_obj.v_res = pixels[1]
        plot_obj.basis = basis

        # Silence output by default. Also set arguments to start in volume
        # calculation mode to avoid loading cross sections
        init_kwargs.setdefault('output', False)
        init_kwargs.setdefault('args', ['-c'])

        with openmc.lib.TemporarySession(self, **init_kwargs):
            return openmc.lib.id_map(plot_obj)

    @add_plot_params
    def plot(
        self,
        origin: Sequence[float] | None = None,
        width: Sequence[float] | None = None,
        pixels: int | Sequence[int] = 40000,
        basis: str = 'xy',
        color_by: str = 'cell',
        colors: dict | None = None,
        seed: int | None = None,
        openmc_exec: PathLike = 'openmc',
        axes=None,
        legend: bool = False,
        axis_units: str = 'cm',
        outline: bool | str = False,
        show_overlaps: bool = False,
        overlap_color: Sequence[int] | str | None = None,
        n_samples: int | None = None,
        plane_tolerance: float = 1.,
        legend_kwargs: dict | None = None,
        source_kwargs: dict | None = None,
        contour_kwargs: dict | None = None,
        **kwargs,
    ):
        """Display a slice plot of the model.

        .. versionadded:: 0.15.1
        """
        import matplotlib.image as mpimg
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt

        check_type('n_samples', n_samples, int | None)
        check_type('plane_tolerance', plane_tolerance, Real)
        if legend_kwargs is None:
            legend_kwargs = {}
        legend_kwargs.setdefault('bbox_to_anchor', (1.05, 1))
        legend_kwargs.setdefault('loc', 2)
        legend_kwargs.setdefault('borderaxespad', 0.0)
        if source_kwargs is None:
            source_kwargs = {}
        source_kwargs.setdefault('marker', 'x')

        # Set indices using basis and create axis labels
        x, y, z = _BASIS_INDICES[basis]
        xlabel, ylabel = f'{basis[0]} [{axis_units}]', f'{basis[1]} [{axis_units}]'

        # Determine extents of plot
        origin, width, pixels = self._set_plot_defaults(
            origin, width, pixels, basis)

        axis_scaling_factor = {'km': 0.00001, 'm': 0.01, 'cm': 1, 'mm': 10}

        x_min = (origin[x] - 0.5*width[0]) * axis_scaling_factor[axis_units]
        x_max = (origin[x] + 0.5*width[0]) * axis_scaling_factor[axis_units]
        y_min = (origin[y] - 0.5*width[1]) * axis_scaling_factor[axis_units]
        y_max = (origin[y] + 0.5*width[1]) * axis_scaling_factor[axis_units]

        # Determine whether any materials contains macroscopic data and if so,
        # set energy mode accordingly
        _energy_mode = self.settings._energy_mode
        for mat in self.geometry.get_all_materials().values():
            if mat._macroscopic is not None:
                self.settings.energy_mode = 'multi-group'
                break

        with TemporaryDirectory() as tmpdir:
            _plot_seed = self.settings.plot_seed
            if seed is not None:
                self.settings.plot_seed = seed

            # Create plot object matching passed arguments
            plot = openmc.Plot()
            plot.origin = origin
            plot.width = width
            plot.pixels = pixels
            plot.basis = basis
            plot.color_by = color_by
            plot.show_overlaps = show_overlaps
            if overlap_color is not None:
                plot.overlap_color = overlap_color
            if colors is not None:
                plot.colors = colors
            self.plots.append(plot)

            # Run OpenMC in geometry plotting mode
            self.plot_geometry(False, cwd=tmpdir, openmc_exec=openmc_exec)

            # Undo changes to model
            self.plots.pop()
            self.settings._plot_seed = _plot_seed
            self.settings._energy_mode = _energy_mode

            # Read image from file
            img_path = Path(tmpdir) / f'plot_{plot.id}.png'
            if not img_path.is_file():
                img_path = img_path.with_suffix('.ppm')
            img = mpimg.imread(str(img_path))

            # Create a figure sized such that the size of the axes within
            # exactly matches the number of pixels specified
            if axes is None:
                px = 1/plt.rcParams['figure.dpi']
                fig, axes = plt.subplots()
                axes.set_xlabel(xlabel)
                axes.set_ylabel(ylabel)
                params = fig.subplotpars
                width = pixels[0]*px/(params.right - params.left)
                height = pixels[1]*px/(params.top - params.bottom)
                fig.set_size_inches(width, height)

            if outline:
                # Combine R, G, B values into a single int
                rgb = (img * 256).astype(int)
                image_value = (rgb[..., 0] << 16) + \
                    (rgb[..., 1] << 8) + (rgb[..., 2])

                # Set default arguments for contour()
                if contour_kwargs is None:
                    contour_kwargs = {}
                contour_kwargs.setdefault('colors', 'k')
                contour_kwargs.setdefault('linestyles', 'solid')
                contour_kwargs.setdefault('algorithm', 'serial')

                axes.contour(
                    image_value,
                    origin="upper",
                    levels=np.unique(image_value),
                    extent=(x_min, x_max, y_min, y_max),
                    **contour_kwargs
                )

            # add legend showing which colors represent which material
            # or cell if that was requested
            if legend:
                if plot.colors == {}:
                    raise ValueError("Must pass 'colors' dictionary if you "
                                     "are adding a legend via legend=True.")

                if color_by == "cell":
                    expected_key_type = openmc.Cell
                else:
                    expected_key_type = openmc.Material

                patches = []
                for key, color in plot.colors.items():

                    if isinstance(key, int):
                        raise TypeError(
                            "Cannot use IDs in colors dict for auto legend.")
                    elif not isinstance(key, expected_key_type):
                        raise TypeError(
                            "Color dict key type does not match color_by")

                    # this works whether we're doing cells or materials
                    label = key.name if key.name != '' else key.id

                    # matplotlib takes RGB on 0-1 scale rather than 0-255. at
                    # this point PlotBase has already checked that 3-tuple
                    # based colors are already valid, so if the length is three
                    # then we know it just needs to be converted to the 0-1
                    # format.
                    if len(color) == 3 and not isinstance(color, str):
                        scaled_color = (
                            color[0]/255, color[1]/255, color[2]/255)
                    else:
                        scaled_color = color

                    key_patch = mpatches.Patch(color=scaled_color, label=label)
                    patches.append(key_patch)

                axes.legend(handles=patches, **legend_kwargs)

            # Plot image and return the axes
            if outline != 'only':
                axes.imshow(img, extent=(x_min, x_max, y_min, y_max), **kwargs)


        if n_samples:
            # Sample external source particles
            particles = self.sample_external_source(n_samples)

            # Get points within tolerance of the slice plane
            slice_value = origin[z]
            xs = []
            ys = []
            tol = plane_tolerance
            for particle in particles:
                if (slice_value - tol < particle.r[z] < slice_value + tol):
                    xs.append(particle.r[x])
                    ys.append(particle.r[y])
            axes.scatter(xs, ys, **source_kwargs)

        return axes

    def sample_external_source(
        self,
        n_samples: int = 1000,
        prn_seed: int | None = None,
        **init_kwargs
    ) -> openmc.ParticleList:
        """Sample external source and return source particles.

        .. versionadded:: 0.15.1

        Parameters
        ----------
        n_samples : int
            Number of samples
        prn_seed : int
            Pseudorandom number generator (PRNG) seed; if None, one will be
            generated randomly.
        **init_kwargs
            Keyword arguments passed to :func:`openmc.lib.init`

        Returns
        -------
        openmc.ParticleList
            List of samples source particles
        """
        import openmc.lib

        # Silence output by default. Also set arguments to start in volume
        # calculation mode to avoid loading cross sections
        init_kwargs.setdefault('output', False)
        init_kwargs.setdefault('args', ['-c'])

        with openmc.lib.TemporarySession(self, **init_kwargs):
            return openmc.lib.sample_external_source(
                n_samples=n_samples, prn_seed=prn_seed
            )

    def apply_tally_results(self, statepoint: PathLike | openmc.StatePoint):
        """Apply results from a statepoint to tally objects on the Model

        Parameters
        ----------
        statepoint : PathLike or openmc.StatePoint
            Statepoint file used to update tally results
        """
        self.tallies.add_results(statepoint)

    def plot_geometry(
        self,
        output: bool = True,
        cwd: PathLike = ".",
        openmc_exec: PathLike = "openmc",
        export_model_xml: bool = True,
        **export_kwargs,
    ):
        """Creates plot images as specified by the Model.plots attribute

        .. versionadded:: 0.13.0

        Parameters
        ----------
        output : bool, optional
            Capture OpenMC output from standard out
        cwd : PathLike, optional
            Path to working directory to run in. Defaults to the current
            working directory.
        openmc_exec : PathLike, optional
            Path to OpenMC executable. Defaults to 'openmc'.
            This only applies to the case when not using the C API.
        export_model_xml : bool, optional
            Exports a single model.xml file rather than separate files. Defaults
            to True.
        **export_kwargs
            Keyword arguments passed to either :meth:`Model.export_to_model_xml`
            or :meth:`Model.export_to_xml`.

        """

        if len(self.plots) == 0:
            # Then there is no volume calculation specified
            raise ValueError("The Model.plots attribute must be specified "
                             "before executing this method!")

        with change_directory(cwd):
            if self.is_initialized:
                # Compute the volumes
                openmc.lib.plot_geometry(output)
            else:
                if export_model_xml:
                    self.export_to_model_xml(**export_kwargs)
                else:
                    self.export_to_xml(**export_kwargs)
                path_input = export_kwargs.get("path", None)
                openmc.plot_geometry(output=output, openmc_exec=openmc_exec,
                                     path_input=path_input)

    def _change_py_lib_attribs(
        self,
        names_or_ids: Iterable[str] | Iterable[int],
        value: float | Iterable[float],
        obj_type: str,
        attrib_name: str,
        density_units: str = "atom/b-cm",
    ):
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
        # The C API has no way to set cell volume or material temperature
        # so lets raise exceptions as needed
        if obj_type == 'cell' and attrib_name == 'volume':
            raise NotImplementedError(
                'Setting a Cell volume is not supported!')
        if obj_type == 'material' and attrib_name == 'temperature':
            raise NotImplementedError(
                'Setting a material temperature is not supported!')

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
            if self.is_initialized:
                obj_by_id = openmc.lib.cells
        else:
            by_name = self._materials_by_name
            by_id = self._materials_by_id
            if self.is_initialized:
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
                if attrib_name == 'density':
                    lib_obj.set_density(value, density_units)
                elif attrib_name == 'temperature':
                    lib_obj.set_temperature(value)
                else:
                    setattr(lib_obj, attrib_name, value)

    def rotate_cells(
        self, names_or_ids: Iterable[str] | Iterable[int], vector: Iterable[float]
    ):
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

    def translate_cells(
        self, names_or_ids: Iterable[str] | Iterable[int], vector: Iterable[float]
    ):
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

        self._change_py_lib_attribs(names_or_ids, vector, 'cell',
                                    'translation')

    def update_densities(
        self,
        names_or_ids: Iterable[str] | Iterable[int],
        density: float,
        density_units: str = "atom/b-cm",
    ):
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

        self._change_py_lib_attribs(names_or_ids, density, 'material',
                                    'density', density_units)

    def update_cell_temperatures(
        self, names_or_ids: Iterable[str] | Iterable[int], temperature: float
    ):
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

    def update_material_volumes(
        self, names_or_ids: Iterable[str] | Iterable[int], volume: float
    ):
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

    def differentiate_depletable_mats(self, diff_volume_method: str = None):
        """Assign distribmats for each depletable material

        .. versionadded:: 0.14.0

        .. versionchanged:: 0.15.1
            diff_volume_method default is None, do not set volumes on the new
            material ovjects. Is now a convenience method for
            differentiate_mats(diff_volume_method, depletable_only=True)

        Parameters
        ----------
        diff_volume_method : str
            Specifies how the volumes of the new materials should be found.
            - None: Do not assign volumes to the new materials (Default)
            - 'divide equally': Divide the original material volume equally between the new materials
            - 'match cell': Set the volume of the material to the volume of the cell they fill
        """
        self.differentiate_mats(diff_volume_method, depletable_only=True)

    def differentiate_mats(self, diff_volume_method: str = None, depletable_only: bool = True):
        """Assign distribmats for each material

        .. versionadded:: 0.15.1

        Parameters
        ----------
        diff_volume_method : str
            Specifies how the volumes of the new materials should be found.
            - None: Do not assign volumes to the new materials (Default)
            - 'divide equally': Divide the original material volume equally between the new materials
            - 'match cell': Set the volume of the material to the volume of the cell they fill
        depletable_only : bool
            Default is True, only depletable materials will be differentiated. If False, all materials will be
            differentiated.
        """
        check_value('volume differentiation method', diff_volume_method, ("divide equally", "match cell", None))

        # Count the number of instances for each cell and material
        self.geometry.determine_paths(instances_only=True)

        # Get list of materials
        if self.materials:
            materials = self.materials
        else:
            materials = list(self.geometry.get_all_materials().values())

        # Find all or depletable_only materials which have multiple instance
        distribmats = set()
        for mat in materials:
            # Differentiate all materials with multiple instances
            diff_mat = mat.num_instances > 1
            # If depletable_only is True, differentiate only depletable materials
            if depletable_only:
                diff_mat = diff_mat and mat.depletable
            if diff_mat:
                # Assign volumes to the materials according to requirements
                if diff_volume_method == "divide equally":
                    if mat.volume is None:
                        raise RuntimeError(
                            "Volume not specified for "
                            f"material with ID={mat.id}.")
                    else:
                        mat.volume /= mat.num_instances
                elif diff_volume_method == "match cell":
                    for cell in self.geometry.get_all_material_cells().values():
                        if cell.fill == mat:
                            if not cell.volume:
                                raise ValueError(
                                    f"Volume of cell ID={cell.id} not specified. "
                                    "Set volumes of cells prior to using "
                                    "diff_volume_method='match cell'.")
                distribmats.add(mat)

        if not distribmats:
            return

        # Assign distribmats to cells
        for cell in self.geometry.get_all_material_cells().values():
            if cell.fill in distribmats:
                mat = cell.fill

                # Clone materials
                if cell.num_instances > 1:
                    cell.fill = [mat.clone() for _ in range(cell.num_instances)]
                else:
                    cell.fill = mat.clone()

                # For 'match cell', assign volumes based on the cells
                if diff_volume_method == 'match cell':
                    if cell.fill_type == 'distribmat':
                        for clone_mat in cell.fill:
                            clone_mat.volume = cell.volume
                    else:
                        cell.fill.volume = cell.volume

        if self.materials is not None:
            self.materials = openmc.Materials(
                self.geometry.get_all_materials().values()
            )

    def _generate_infinite_medium_mgxs(
        self,
        groups: openmc.mgxs.EnergyGroups,
        nparticles: int,
        mgxs_path: PathLike,
        correction: str | None,
        directory: PathLike,
    ):
        """Generate a MGXS library by running multiple OpenMC simulations, each
        representing an infinite medium simulation of a single isolated
        material. A discrete source is used to sample particles, with an equal
        strength spread across each of the energy groups. This is a highly naive
        method that ignores all spatial self shielding effects and all resonance
        shielding effects between materials.

        Parameters
        ----------
        groups : openmc.mgxs.EnergyGroups
            Energy group structure for the MGXS.
        nparticles : int
            Number of particles to simulate per batch when generating MGXS.
        mgxs_path : str
            Filename for the MGXS HDF5 file.
        correction : str
            Transport correction to apply to the MGXS. Options are None and
            "P0".
        directory : str
            Directory to run the simulation in, so as to contain XML files.
        """
        warnings.warn("The infinite medium method of generating MGXS may hang "
                      "if a material has a k-infinity > 1.0.")
        mgxs_sets = []
        for material in self.materials:
            model = openmc.Model()

            # Set materials on the model
            model.materials = [material]

            # Settings
            model.settings.batches = 100
            model.settings.particles = nparticles
            model.settings.run_mode = 'fixed source'

            # Make a discrete source that is uniform over the bins of the group structure
            n_groups = groups.num_groups
            midpoints = []
            strengths = []
            for i in range(n_groups):
                bounds = groups.get_group_bounds(i+1)
                midpoints.append((bounds[0] + bounds[1]) / 2.0)
                strengths.append(1.0)

            energy_distribution = openmc.stats.Discrete(x=midpoints, p=strengths)
            model.settings.source = openmc.IndependentSource(
                space=openmc.stats.Point(), energy=energy_distribution)
            model.settings.output = {'summary': True, 'tallies': False}

            # Geometry
            box = openmc.model.RectangularPrism(
                100000.0, 100000.0, boundary_type='reflective')
            name = material.name
            infinite_cell = openmc.Cell(name=name, fill=material, region=-box)
            infinite_universe = openmc.Universe(name=name, cells=[infinite_cell])
            model.geometry.root_universe = infinite_universe

            # Add MGXS Tallies

            # Initialize MGXS library with a finished OpenMC geometry object
            mgxs_lib = openmc.mgxs.Library(model.geometry)

            # Pick energy group structure
            mgxs_lib.energy_groups = groups

            # Disable transport correction
            mgxs_lib.correction = correction

            # Specify needed cross sections for random ray
            if correction == 'P0':
                mgxs_lib.mgxs_types = [
                    'nu-transport', 'absorption', 'nu-fission', 'fission',
                    'consistent nu-scatter matrix', 'multiplicity matrix', 'chi'
                ]
            elif correction is None:
                mgxs_lib.mgxs_types = [
                    'total', 'absorption', 'nu-fission', 'fission',
                    'consistent nu-scatter matrix', 'multiplicity matrix', 'chi'
                ]

            # Specify a "cell" domain type for the cross section tally filters
            mgxs_lib.domain_type = "material"

            # Specify the cell domains over which to compute multi-group cross sections
            mgxs_lib.domains = model.geometry.get_all_materials().values()

            # Do not compute cross sections on a nuclide-by-nuclide basis
            mgxs_lib.by_nuclide = False

            # Check the library - if no errors are raised, then the library is satisfactory.
            mgxs_lib.check_library_for_openmc_mgxs()

            # Construct all tallies needed for the multi-group cross section library
            mgxs_lib.build_library()

            # Create a "tallies.xml" file for the MGXS Library
            mgxs_lib.add_to_tallies_file(model.tallies, merge=True)

            # Run
            statepoint_filename = model.run(cwd=directory)

            # Load MGXS
            with openmc.StatePoint(statepoint_filename) as sp:
                mgxs_lib.load_from_statepoint(sp)

            # Create a MGXS File which can then be written to disk
            mgxs_set = mgxs_lib.get_xsdata(domain=material, xsdata_name=name)
            mgxs_sets.append(mgxs_set)

        # Write the file to disk
        mgxs_file = openmc.MGXSLibrary(energy_groups=groups)
        for mgxs_set in mgxs_sets:
            mgxs_file.add_xsdata(mgxs_set)
        mgxs_file.export_to_hdf5(mgxs_path)

    @staticmethod
    def _create_stochastic_slab_geometry(
        materials: Sequence[openmc.Material],
        cell_thickness: float = 1.0,
        num_repeats: int = 100,
    ) -> tuple[openmc.Geometry, openmc.stats.Box]:
        """Create a geometry representing a stochastic "sandwich" of materials in a
        layered slab geometry. To reduce the impact of the order of materials in
        the slab, the materials are applied to 'num_repeats' different randomly
        positioned layers of 'cell_thickness' each.

        Parameters
        ----------
        materials : list of openmc.Material
            List of materials to assign. Each material will appear exactly num_repeats times,
            then the ordering is randomly shuffled.
        cell_thickness : float, optional
            Thickness of each lattice cell in x (default 1.0 cm).
        num_repeats : int, optional
            Number of repeats for each material (default 100).

        Returns
        -------
        geometry : openmc.Geometry
            The constructed geometry.
        box : openmc.stats.Box
            A spatial sampling distribution covering the full slab domain.
        """
        if not materials:
            raise ValueError("At least one material must be provided.")

        num_materials = len(materials)
        total_cells = num_materials * num_repeats
        total_width = total_cells * cell_thickness

        # Generate an infinite cell/universe for each material
        universes = []
        for i in range(num_materials):
            cell = openmc.Cell(fill=materials[i])
            universes.append(openmc.Universe(cells=[cell]))

        # Make a list of randomized material idx assignments for the stochastic slab
        assignments = list(range(num_materials)) * num_repeats
        random.seed(42)
        random.shuffle(assignments)

        # Create a list of the (randomized) universe assignments to be used
        # when defining the problem lattice.
        lattice_entries = [universes[m] for m in assignments]

        # Create the RectLattice for the 1D material variation in x.
        lattice = openmc.RectLattice()
        lattice.pitch = (cell_thickness, total_width, total_width)
        lattice.lower_left = (0.0, 0.0, 0.0)
        lattice.universes = [[lattice_entries]]
        lattice.outer = universes[0]

        # Define the six outer surfaces with reflective boundary conditions
        rpp = openmc.model.RectangularParallelepiped(
            0.0, total_width, 0.0, total_width, 0.0, total_width,
            boundary_type='reflective'
        )

        # Create an outer cell that fills with the lattice.
        outer_cell = openmc.Cell(fill=lattice, region=-rpp)

        # Build the geometry
        geometry = openmc.Geometry([outer_cell])

        # Define the spatial distribution that covers the full cubic domain
        box = openmc.stats.Box(*outer_cell.bounding_box)

        return geometry, box

    def _generate_stochastic_slab_mgxs(
        self,
        groups: openmc.mgxs.EnergyGroups,
        nparticles: int,
        mgxs_path: PathLike,
        correction: str | None,
        directory: PathLike,
    ) -> None:
        """Generate MGXS assuming a stochastic "sandwich" of materials in a layered
        slab geometry. While geometry-specific spatial shielding effects are not
        captured, this method can be useful when the geometry has materials only
        found far from the source region that the "material_wise" method would
        not be capable of generating cross sections for. Conversely, this method
        will generate cross sections for all materials in the problem regardless
        of type. If this is a fixed source problem, a discrete source is used to
        sample particles, with an equal strength spread across each of the
        energy groups.

        Parameters
        ----------
        groups : openmc.mgxs.EnergyGroups
            Energy group structure for the MGXS.
        nparticles : int
            Number of particles to simulate per batch when generating MGXS.
        mgxs_path : str
            Filename for the MGXS HDF5 file.
        correction : str
            Transport correction to apply to the MGXS. Options are None and
            "P0".
        directory : str
            Directory to run the simulation in, so as to contain XML files.
        """
        model = openmc.Model()
        model.materials = self.materials

        # Settings
        model.settings.batches = 200
        model.settings.inactive = 100
        model.settings.particles = nparticles
        model.settings.output = {'summary': True, 'tallies': False}
        model.settings.run_mode = self.settings.run_mode

        # Stochastic slab geometry
        model.geometry, spatial_distribution = Model._create_stochastic_slab_geometry(
            model.materials)

        # Make a discrete source that is uniform over the bins of the group structure
        n_groups = groups.num_groups
        midpoints = []
        strengths = []
        for i in range(n_groups):
            bounds = groups.get_group_bounds(i+1)
            midpoints.append((bounds[0] + bounds[1]) / 2.0)
            strengths.append(1.0)

        energy_distribution = openmc.stats.Discrete(x=midpoints, p=strengths)
        model.settings.source = [openmc.IndependentSource(
            space=spatial_distribution, energy=energy_distribution, strength=1.0)]

        model.settings.output = {'summary': True, 'tallies': False}

        # Add MGXS Tallies

        # Initialize MGXS library with a finished OpenMC geometry object
        mgxs_lib = openmc.mgxs.Library(model.geometry)

        # Pick energy group structure
        mgxs_lib.energy_groups = groups

        # Disable transport correction
        mgxs_lib.correction = correction

       # Specify needed cross sections for random ray
        if correction == 'P0':
            mgxs_lib.mgxs_types = ['nu-transport', 'absorption', 'nu-fission', 'fission',
                                   'consistent nu-scatter matrix', 'multiplicity matrix', 'chi']
        elif correction is None:
            mgxs_lib.mgxs_types = ['total', 'absorption', 'nu-fission', 'fission',
                                   'consistent nu-scatter matrix', 'multiplicity matrix', 'chi']

        # Specify a "cell" domain type for the cross section tally filters
        mgxs_lib.domain_type = "material"

        # Specify the cell domains over which to compute multi-group cross sections
        mgxs_lib.domains = model.geometry.get_all_materials().values()

        # Do not compute cross sections on a nuclide-by-nuclide basis
        mgxs_lib.by_nuclide = False

        # Check the library - if no errors are raised, then the library is satisfactory.
        mgxs_lib.check_library_for_openmc_mgxs()

        # Construct all tallies needed for the multi-group cross section library
        mgxs_lib.build_library()

        # Create a "tallies.xml" file for the MGXS Library
        mgxs_lib.add_to_tallies_file(model.tallies, merge=True)

        # Run
        statepoint_filename = model.run(cwd=directory)

        # Load MGXS
        with openmc.StatePoint(statepoint_filename) as sp:
            mgxs_lib.load_from_statepoint(sp)

        names = [mat.name for mat in mgxs_lib.domains]

        # Create a MGXS File which can then be written to disk
        mgxs_file = mgxs_lib.create_mg_library(xs_type='macro', xsdata_names=names)
        mgxs_file.export_to_hdf5(mgxs_path)

    def _generate_material_wise_mgxs(
        self,
        groups: openmc.mgxs.EnergyGroups,
        nparticles: int,
        mgxs_path: PathLike,
        correction: str | None,
        directory: PathLike,
    ) -> None:
        """Generate a material-wise MGXS library for the model by running the
        original continuous energy OpenMC simulation of the full material
        geometry and source, and tally MGXS data for each material. This method
        accurately conserves reaction rates totaled over the entire simulation
        domain. However, when the geometry has materials only found far from the
        source region, it is possible the Monte Carlo solver may not be able to
        score any tallies to these material types, thus resulting in zero cross
        section values for these materials. For such cases, the "stochastic
        slab" method may be more appropriate.

        Parameters
        ----------
        groups : openmc.mgxs.EnergyGroups
            Energy group structure for the MGXS.
        nparticles : int
            Number of particles to simulate per batch when generating MGXS.
        mgxs_path : PathLike
            Filename for the MGXS HDF5 file.
        correction : str
            Transport correction to apply to the MGXS. Options are None and
            "P0".
        directory : PathLike
            Directory to run the simulation in, so as to contain XML files.
        """
        model = copy.deepcopy(self)
        model.tallies = openmc.Tallies()

        # Settings
        model.settings.batches = 200
        model.settings.inactive = 100
        model.settings.particles = nparticles
        model.settings.output = {'summary': True, 'tallies': False}

        # Add MGXS Tallies

        # Initialize MGXS library with a finished OpenMC geometry object
        mgxs_lib = openmc.mgxs.Library(model.geometry)

        # Pick energy group structure
        mgxs_lib.energy_groups = groups

        # Disable transport correction
        mgxs_lib.correction = correction

        # Specify needed cross sections for random ray
        if correction == 'P0':
            mgxs_lib.mgxs_types = [
                'nu-transport', 'absorption', 'nu-fission', 'fission',
                'consistent nu-scatter matrix', 'multiplicity matrix', 'chi'
            ]
        elif correction is None:
            mgxs_lib.mgxs_types = [
                'total', 'absorption', 'nu-fission', 'fission',
                'consistent nu-scatter matrix', 'multiplicity matrix', 'chi'
            ]

        # Specify a "cell" domain type for the cross section tally filters
        mgxs_lib.domain_type = "material"

        # Specify the cell domains over which to compute multi-group cross sections
        mgxs_lib.domains = model.geometry.get_all_materials().values()

        # Do not compute cross sections on a nuclide-by-nuclide basis
        mgxs_lib.by_nuclide = False

        # Check the library - if no errors are raised, then the library is satisfactory.
        mgxs_lib.check_library_for_openmc_mgxs()

        # Construct all tallies needed for the multi-group cross section library
        mgxs_lib.build_library()

        # Create a "tallies.xml" file for the MGXS Library
        mgxs_lib.add_to_tallies_file(model.tallies, merge=True)

        # Run
        statepoint_filename = model.run(cwd=directory)

        # Load MGXS
        with openmc.StatePoint(statepoint_filename) as sp:
            mgxs_lib.load_from_statepoint(sp)

        names = [mat.name for mat in mgxs_lib.domains]

        # Create a MGXS File which can then be written to disk
        mgxs_file = mgxs_lib.create_mg_library(
            xs_type='macro', xsdata_names=names)
        mgxs_file.export_to_hdf5(mgxs_path)

    def convert_to_multigroup(
        self,
        method: str = "material_wise",
        groups: str = "CASMO-2",
        nparticles: int = 2000,
        overwrite_mgxs_library: bool = False,
        mgxs_path: PathLike = "mgxs.h5",
        correction: str | None = None,
    ):
        """Convert all materials from continuous energy to multigroup.

        If no MGXS data library file is found, generate one using one or more
        continuous energy Monte Carlo simulations.

        Parameters
        ----------
        method : {"material_wise", "stochastic_slab", "infinite_medium"}, optional
            Method to generate the MGXS.
        groups : openmc.mgxs.EnergyGroups or str, optional
            Energy group structure for the MGXS or the name of the group
            structure (based on keys from openmc.mgxs.GROUP_STRUCTURES).
        mgxs_path : str, optional
            Filename of the mgxs.h5 library file.
        correction : str, optional
            Transport correction to apply to the MGXS. Options are None and
            "P0".
        """
        if isinstance(groups, str):
            groups = openmc.mgxs.EnergyGroups(groups)

        # Do all work (including MGXS generation) in a temporary directory
        # to avoid polluting the working directory with residual XML files
        with TemporaryDirectory() as tmpdir:

            # Determine if there are DAGMC universes in the model. If so, we need to synchronize
            # the dagmc materials with cells.
            # TODO: Can this be done without having to init/finalize?
            for univ in self.geometry.get_all_universes().values():
                if isinstance(univ, openmc.DAGMCUniverse):
                    self.init_lib(directory=tmpdir)
                    self.sync_dagmc_universes()
                    self.finalize_lib()
                    break

            # Make sure all materials have a name, and that the name is a valid HDF5
            # dataset name
            for material in self.materials:
                if not material.name or not material.name.strip():
                    material.name = f"material {material.id}"
                material.name = re.sub(r'[^a-zA-Z0-9]', '_', material.name)

            # If needed, generate the needed MGXS data library file
            if not Path(mgxs_path).is_file() or overwrite_mgxs_library:
                if method == "infinite_medium":
                    self._generate_infinite_medium_mgxs(
                        groups, nparticles, mgxs_path, correction, tmpdir)
                elif method == "material_wise":
                    self._generate_material_wise_mgxs(
                        groups, nparticles, mgxs_path, correction, tmpdir)
                elif method == "stochastic_slab":
                    self._generate_stochastic_slab_mgxs(
                        groups, nparticles, mgxs_path, correction, tmpdir)
                else:
                    raise ValueError(
                        f'MGXS generation method "{method}" not recognized')
            else:
                print(f'Existing MGXS library file "{mgxs_path}" will be used')

            # Convert all continuous energy materials to multigroup
            self.materials.cross_sections = mgxs_path
            for material in self.materials:
                material.set_density('macro', 1.0)
                material._nuclides = []
                material._sab = []
                material.add_macroscopic(material.name)

            self.settings.energy_mode = 'multi-group'

    def convert_to_random_ray(self):
        """Convert a multigroup model to use random ray.

        This method determines values for the needed settings and adds them to
        the settings.random_ray dictionary so as to enable random ray mode. The
        settings that are populated are:

        - 'ray_source' (openmc.IndependentSource): Where random ray starting
          points are sampled from.
        - 'distance_inactive' (float): The "dead zone" distance at the beginning
          of the ray.
        - 'distance_active' (float): The "active" distance of the ray
        - 'particles' (int): Number of rays to simulate

        The method will determine reasonable defaults for each of the above
        variables based on analysis of the model's geometry. The function will
        have no effect if the random ray dictionary is already defined in the
        model settings.
        """
        # If the random ray dictionary is already set, don't overwrite it
        if self.settings.random_ray:
            warnings.warn("Random ray conversion skipped as "
                          "settings.random_ray dictionary is already set.")
            return

        if self.settings.energy_mode != 'multi-group':
            raise ValueError(
                "Random ray conversion failed: energy mode must be "
                "'multi-group'. Use convert_to_multigroup() first."
            )

        # Helper function for detecting infinity
        def _replace_infinity(value):
            if np.isinf(value):
                return 1.0 if value > 0 else -1.0
            return value

        # Get a bounding box for sampling rays. We can utilize the geometry's bounding box
        # though for 2D problems we need to detect the infinities and replace them with an
        # arbitrary finite value.
        bounding_box = self.geometry.bounding_box
        lower_left = [_replace_infinity(v) for v in bounding_box.lower_left]
        upper_right = [_replace_infinity(v) for v in bounding_box.upper_right]
        uniform_dist_ray = openmc.stats.Box(lower_left, upper_right)
        rr_source = openmc.IndependentSource(space=uniform_dist_ray)
        self.settings.random_ray['ray_source'] = rr_source

        # For the dead zone and active length, a reasonable guess is the larger of either:
        # 1) The maximum chord length through the geometry (as defined by its bounding box)
        # 2) 30 cm
        # Then, set the active length to be 5x longer than the dead zone length, for the sake of efficiency.
        chord_length = np.array(upper_right) - np.array(lower_left)
        max_length = max(np.linalg.norm(chord_length), 30.0)

        self.settings.random_ray['distance_inactive'] = max_length
        self.settings.random_ray['distance_active'] = 5 * max_length

        # Take a wild guess as to how many rays are needed
        self.settings.particles = 2 * int(max_length)

    def keff_search(
        self,
        func: ModelModifier,
        x0: float,
        x1: float,
        target: float = 1.0,
        k_tol: float = 1e-4,
        sigma_final: float = 3e-4,
        p: float = 0.5,
        q: float = 0.95,
        memory: int = 4,
        x_min: float | None = None,
        x_max: float | None = None,
        b0: int | None = None,
        b_min: int = 20,
        b_max: int | None = None,
        maxiter: int = 50,
        output: bool = False,
        func_kwargs: dict[str, Any] | None = None,
        run_kwargs: dict[str, Any] | None = None,
    ) -> SearchResult:
        r"""Perform a keff search on a model parametrized by a single variable.

        This method uses the GRsecant method described in a paper by `Price and
        Roskoff <https://doi.org/10.1016/j.pnucene.2023.104731>`_. The GRsecant
        method is a modification of the secant method that accounts for
        uncertainties in the function evaluations. The method uses a weighted
        linear fit of the most recent function evaluations to predict the next
        point to evaluate. It also adaptively changes the number of batches to
        meet the target uncertainty value at each iteration.

        The target uncertainty for iteration :math:`n+1` is determined by the
        following equation (following Eq. (8) in the paper):

        .. math::
            \sigma_{i+1} = q \sigma_\text{final} \left ( \frac{ \min \left \{
            \left\lvert k_i - k_\text{target} \right\rvert : k=0,1,\dots,n
            \right \} }{k_\text{tol}} \right )^p

        where :math:`q` is a multiplicative factor less than 1, given as the
        ``sigma_factor`` parameter below.

        Parameters
        ----------
        func : ModelModifier
            Function that takes the parameter to be searched and makes a
            modification to the model.
        x0 : float
            First guess for the parameter passed to `func`
        x1 : float
            Second guess for the parameter passed to `func`
        target : float, optional
            keff value to search for
        k_tol : float, optional
            Stopping criterion on the function value; the absolute value must be
            within ``k_tol`` of zero to be accepted.
        sigma_final : float, optional
            Maximum accepted k-effective uncertainty for the stopping criterion.
        p : float, optional
            Exponent used in the stopping criterion.
        q : float, optional
            Multiplicative factor used in the stopping criterion.
        memory : int, optional
            Number of most-recent points used in the weighted linear fit of
            ``f(x) = a + b x`` to predict the next point.
        x_min : float, optional
            Minimum allowed value for the parameter ``x``.
        x_max : float, optional
            Maximum allowed value for the parameter ``x``.
        b0 : int, optional
            Number of active batches to use for the initial function
            evaluations. If None, uses the model's current setting.
        b_min : int, optional
            Minimum number of active batches to use in a function evaluation.
        b_max : int, optional
            Maximum number of active batches to use in a function evaluation.
        maxiter : int, optional
            Maximum number of iterations to perform.
        output : bool, optional
            Whether or not to display output showing iteration progress.
        func_kwargs : dict, optional
            Keyword-based arguments to pass to the `func` function.
        run_kwargs : dict, optional
            Keyword arguments to pass to :meth:`openmc.Model.run` or
            :meth:`openmc.lib.run`.

        Returns
        -------
        SearchResult
            Result object containing the estimated root (parameter value) and
            evaluation history (parameters, means, standard deviations, and
            batches), plus convergence status and termination reason.

        """
        import openmc.lib

        check_type('model modifier', func, Callable)
        check_type('target', target, Real)
        if memory < 2:
            raise ValueError("memory must be ≥ 2")
        func_kwargs = {} if func_kwargs is None else dict(func_kwargs)
        run_kwargs = {} if run_kwargs is None else dict(run_kwargs)
        run_kwargs.setdefault('output', False)

        # Create lists to store the history of evaluations
        xs: list[float] = []
        fs: list[float] = []
        ss: list[float] = []
        gs: list[int] = []
        count = 0

        # Helper function to evaluate f and store results
        def eval_at(x: float, batches: int) -> tuple[float, float]:
            # Modify the model with the current guess
            func(x, **func_kwargs)

            # Change the number of batches and run the model
            batches += self.settings.inactive
            if openmc.lib.is_initialized:
                openmc.lib.settings.set_batches(batches)
                openmc.lib.reset()
                openmc.lib.run(**run_kwargs)
                sp_filepath = f'statepoint.{batches}.h5'
            else:
                self.settings.batches = batches
                sp_filepath = self.run(**run_kwargs)

            # Extract keff and its uncertainty
            with openmc.StatePoint(sp_filepath) as sp:
                keff = sp.keff

            if output:
                nonlocal count
                count += 1
                print(f'Iteration {count}: {batches=}, {x=:.6g}, {keff=:.5f}')

            xs.append(float(x))
            fs.append(float(keff.n - target))
            ss.append(float(keff.s))
            gs.append(int(batches))
            return fs[-1], ss[-1]

        # Default b0 to current model settings if not explicitly provided
        if b0 is None:
            b0 = self.settings.batches - self.settings.inactive

        # Perform the search (inlined GRsecant) in a temporary directory
        with TemporaryDirectory() as tmpdir:
            if not openmc.lib.is_initialized:
                run_kwargs.setdefault('cwd', tmpdir)

            # ---- Seed with two evaluations
            f0, s0 = eval_at(x0, b0)
            if abs(f0) <= k_tol and s0 <= sigma_final:
                return SearchResult(x0, xs, fs, ss, gs, True, "converged")
            f1, s1 = eval_at(x1, b0)
            if abs(f1) <= k_tol and s1 <= sigma_final:
                return SearchResult(x1, xs, fs, ss, gs, True, "converged")

            for _ in range(maxiter - 2):
                # ------ Step 1: propose next x via GRsecant
                m = min(memory, len(xs))

                # Perform a curve fit on f(x) = a + bx accounting for
                # uncertainties. This is equivalent to minimizing the function
                # in Equation (A.14)
                (a, b), _ = curve_fit(
                    lambda x, a, b: a + b*x,
                    xs[-m:], fs[-m:], sigma=ss[-m:], absolute_sigma=True
                )
                x_new = float(-a / b)

                # Clamp x_new to the bounds if provided
                if x_min is not None:
                    x_new = max(x_new, x_min)
                if x_max is not None:
                    x_new = min(x_new, x_max)

                # ------ Step 2: choose target σ for next run (Eq. 8 + clamp)

                min_abs_f = float(np.min(np.abs(fs)))
                base = q * sigma_final
                ratio = min_abs_f / k_tol if k_tol > 0 else 1.0
                sig = base * (ratio ** p)
                sig_target = max(sig, base)

                # ------ Step 3: choose generations to hit σ_target (Appendix C)

                # Use at least two past points for regression
                if len(gs) >= 2 and np.var(np.log(gs)) > 0.0:
                    # Perform a curve fit based on Eq. (C.3) to solve for ln(k).
                    # Note that unlike in the paper, we do not leave r as an
                    # undetermined parameter and choose r=0.5.
                    (ln_k,), _ = curve_fit(
                        lambda ln_b, ln_k: ln_k - 0.5*ln_b,
                        np.log(gs[-4:]), np.log(ss[-4:]),
                    )
                    k = float(np.exp(ln_k))
                else:
                    k = float(ss[-1] * math.sqrt(gs[-1]))

                b_new = (k / sig_target) ** 2

                # Clamp and round up to integer
                b_new = max(b_min, math.ceil(b_new))
                if b_max is not None:
                    b_new = min(b_new, b_max)

                # Evaluate at proposed x with batches determined above
                f_new, s_new = eval_at(x_new, b_new)

                # Termination based on both criteria (|f| and σ)
                if abs(f_new) <= k_tol and s_new <= sigma_final:
                    return SearchResult(x_new, xs, fs, ss, gs, True, "converged")

            return SearchResult(xs[-1], xs, fs, ss, gs, False, "maxiter")


@dataclass
class SearchResult:
    """Result of a GRsecant keff search.

    Attributes
    ----------
    root : float
        Estimated parameter value where f(x) = 0 at termination.
    parameters : list[float]
        Parameter values (x) evaluated during the search, in order.
    keffs : list[float]
        Estimated keff values for each evaluation.
    stdevs : list[float]
        One-sigma uncertainties of keff for each evaluation.
    batches : list[int]
        Number of active batches used for each evaluation.
    converged : bool
        Whether both |f| <= k_tol and sigma <= sigma_final were met.
    flag : str
        Reason for termination (e.g., "converged", "maxiter").
    """
    root: float
    parameters: list[float] = field(repr=False)
    means: list[float] = field(repr=False)
    stdevs: list[float] = field(repr=False)
    batches: list[int] = field(repr=False)
    converged: bool
    flag: str

    @property
    def function_calls(self) -> int:
        """Number of function evaluations performed."""
        return len(self.parameters)

    @property
    def total_batches(self) -> int:
        """Total number of active batches used across all evaluations."""
        return sum(self.batches)


