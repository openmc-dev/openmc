from __future__ import annotations
from collections.abc import Sequence
from contextlib import nullcontext
import copy
from datetime import datetime
import json
from numbers import Integral
from pathlib import Path

import numpy as np
import openmc
from . import IndependentOperator, PredictorIntegrator
from .chain import Chain
from .microxs import get_microxs_and_flux, write_microxs_hdf5, read_microxs_hdf5
from .results import Results
from ..checkvalue import PathLike
from ..mpi import comm
from openmc.lib import TemporarySession


def get_activation_materials(
    model: openmc.Model,
    mmv_list: list[openmc.MeshMaterialVolumes]
) -> openmc.Materials:
    """Get a list of activation materials for each mesh element/material.

    When performing a mesh-based R2S calculation, a unique material is needed
    for each activation region, which is a combination of a mesh element and a
    material within that mesh element. This function generates a list of such
    materials, each with a unique name and volume corresponding to the mesh
    element and material.

    Parameters
    ----------
    model : openmc.Model
        The full model containing the geometry and materials.
    mmv_list : list of openmc.MeshMaterialVolumes
        List of mesh material volumes objects, one per mesh, containing the
        materials and their volumes for each mesh element.

    Returns
    -------
    openmc.Materials
        A list of materials, each corresponding to a unique mesh element and
        material combination across all meshes.

    """
    # Get all materials in the model
    material_dict = model._get_all_materials()

    # Create a new activation material for each element-material combination
    # across all meshes
    materials = openmc.Materials()
    for mesh_idx, mmv in enumerate(mmv_list):
        mat_ids = mmv._materials[mmv._materials > -1]
        volumes = mmv._volumes[mmv._materials > -1]
        elems, _ = np.where(mmv._materials > -1)

        for elem, mat_id, vol in zip(elems, mat_ids, volumes):
            mat = material_dict[mat_id]
            new_mat = mat.clone()
            new_mat.depletable = True
            new_mat.name = f'Mesh {mesh_idx}, Element {elem}, Material {mat_id}'
            new_mat.volume = vol
            materials.append(new_mat)

    return materials


class R2SManager:
    """Manager for Rigorous 2-Step (R2S) method calculations.

    This class is responsible for managing the materials and sources needed for
    mesh-based or cell-based R2S calculations. It provides methods to get
    activation materials and decay photon sources based on the mesh/cells and
    materials in the OpenMC model. Multiple meshes can be specified as domains,
    in which case each element--material combination of each mesh is treated as
    an activation region (meshes are assumed to be non-overlapping).

    This class supports the use of a different models for the neutron and photon
    transport calculation. However, for cell-based calculations, it assumes that
    the only changes in the model are material assignments. For mesh-based
    calculations, it checks material assignments in the photon model and any
    element--material combinations that don't appear in the photon model are
    skipped.

    Parameters
    ----------
    neutron_model : openmc.Model
        The OpenMC model to use for neutron transport.
    domains : openmc.MeshBase or Sequence[openmc.MeshBase] or Sequence[openmc.Cell]
        The mesh(es) or a sequence of cells that represent the spatial units
        over which the R2S calculation will be performed. When a single
        :class:`~openmc.MeshBase` or a sequence of meshes is given, each
        element--material combination across all meshes is treated as an
        activation region.
    photon_model : openmc.Model, optional
        The OpenMC model to use for photon transport calculations. If None, a
        shallow copy of the neutron_model will be created and used.

    Attributes
    ----------
    domains : list of openmc.MeshBase or Sequence[openmc.Cell]
        The meshes or a sequence of cells that represent the spatial units over
        which the R2S calculation will be performed.
    neutron_model : openmc.Model
        The OpenMC model used for neutron transport.
    photon_model : openmc.Model
        The OpenMC model used for photon transport calculations.
    method : {'mesh-based', 'cell-based'}
        Indicates whether the R2S calculation uses mesh elements ('mesh-based')
        as the spatial discretization or a list of cells ('cell-based').
    results : dict
        A dictionary that stores results from the R2S calculation.

    """
    def __init__(
        self,
        neutron_model: openmc.Model,
        domains: openmc.MeshBase | Sequence[openmc.MeshBase] | Sequence[openmc.Cell],
        photon_model: openmc.Model | None = None,
    ):
        self.neutron_model = neutron_model
        if photon_model is None:
            # Create a shallow copy of the neutron model for photon transport
            self.photon_model = openmc.Model(
                geometry=copy.copy(neutron_model.geometry),
                materials=copy.copy(neutron_model.materials),
                settings=copy.copy(neutron_model.settings),
                tallies=copy.copy(neutron_model.tallies),
                plots=copy.copy(neutron_model.plots),
            )
        else:
            self.photon_model = photon_model
        if isinstance(domains, openmc.MeshBase):
            self.method = 'mesh-based'
            self.domains = [domains]
        elif isinstance(domains, Sequence) and len(domains) > 0 and \
                isinstance(domains[0], openmc.MeshBase):
            self.method = 'mesh-based'
            self.domains = list(domains)
        else:
            self.method = 'cell-based'
            self.domains = list(domains)
        self.results = {}

    def run(
        self,
        timesteps: Sequence[float] | Sequence[tuple[float, str]],
        source_rates: float | Sequence[float],
        timestep_units: str = 's',
        photon_time_indices: Sequence[int] | None = None,
        output_dir: PathLike | None = None,
        bounding_boxes: dict[int, openmc.BoundingBox] | None = None,
        chain_file: PathLike | Chain | None = None,
        micro_kwargs: dict | None = None,
        mat_vol_kwargs: dict | None = None,
        run_kwargs: dict | None = None,
        operator_kwargs: dict | None = None,
        by_parent_nuclide: bool = False,
    ):
        """Run the R2S calculation.

        Parameters
        ----------
        timesteps : Sequence[float] or Sequence[tuple[float, str]]
            Sequence of timesteps. Note that values are not cumulative. The
            units are specified by the `timestep_units` argument when
            `timesteps` is an iterable of float. Alternatively, units can be
            specified for each step by passing an iterable of (value, unit)
            tuples.
        source_rates : float or Sequence[float]
            Source rate in [neutron/sec] for each interval in `timesteps`.
        timestep_units : {'s', 'min', 'h', 'd', 'a'}, optional
            Units for values specified in the `timesteps` argument when passing
            float values. 's' means seconds, 'min' means minutes, 'h' means
            hours, 'd' means days, and 'a' means years (Julian).
        photon_time_indices : Sequence[int], optional
            Sequence of time indices at which photon transport should be run;
            represented as indices into the array of times formed by the
            timesteps. For example, if two timesteps are specified, the array of
            times would contain three entries, and [2] would indicate computing
            photon results at the last time. A value of None indicates to run
            photon transport at each time that has a decay photon source.
        output_dir : PathLike, optional
            Path to directory where R2S calculation outputs will be saved. If
            not provided, a timestamped directory 'r2s_YYYY-MM-DDTHH-MM-SS' is
            created. Subdirectories will be created for the neutron transport,
            activation, and photon transport steps.
        bounding_boxes : dict[int, openmc.BoundingBox], optional
            Dictionary mapping cell IDs to bounding boxes used for spatial
            source sampling in cell-based R2S calculations. Required if method
            is 'cell-based'.
        chain_file : PathLike or openmc.deplete.Chain, optional
            Path to the depletion chain XML file or depletion chain object to
            use during activation. If not provided, the default configured
            chain file will be used.
        micro_kwargs : dict, optional
            Additional keyword arguments passed to
            :func:`openmc.deplete.get_microxs_and_flux` during the neutron
            transport step.
        mat_vol_kwargs : dict, optional
            Additional keyword arguments passed to
            :meth:`openmc.MeshBase.material_volumes`.
        run_kwargs : dict, optional
            Additional keyword arguments passed to :meth:`openmc.Model.run`
            during the neutron and photon transport step. By default, output is
            disabled.
        operator_kwargs : dict, optional
            Additional keyword arguments passed to
            :class:`openmc.deplete.IndependentOperator`.
        by_parent_nuclide : bool, optional
            Whether to score photon tallies separately for each parent
            radionuclide. A :class:`~openmc.ParentNuclideFilter` is added to
            tallies that do not already contain one, with bins determined from
            the prepared decay photon sources.

        Returns
        -------
        Path
            Path to the output directory containing all calculation results
        """

        if output_dir is None:
            # Create timestamped output directory and broadcast to all ranks for
            # consistency (different ranks may have slightly different times)
            stamp = datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
            output_dir = Path(comm.bcast(f'r2s_{stamp}'))
        else:
            output_dir = Path(output_dir)

        # Set run_kwargs for the neutron transport step
        if micro_kwargs is None:
            micro_kwargs = {}
        if run_kwargs is None:
            run_kwargs = {}
        if operator_kwargs is None:
            operator_kwargs = {}
        run_kwargs.setdefault('output', False)
        micro_kwargs.setdefault('run_kwargs', run_kwargs)

        # DecaySpectrum distributions are resolved in the C++ solver using
        # OPENMC_CHAIN_FILE. If a Chain object was passed, write an XML
        # representation alongside the R2S outputs.
        if isinstance(chain_file, Chain):
            output_dir.mkdir(parents=True, exist_ok=True)
            chain_path = output_dir / 'chain.xml'
            if comm.rank == 0:
                chain_file.export_to_xml(chain_path)
            comm.barrier()
        else:
            chain_path = chain_file

        chain_context = (
            openmc.config.patch('chain_file', chain_path)
            if chain_path is not None else nullcontext()
        )
        with chain_context:
            self.step1_neutron_transport(
                output_dir / 'neutron_transport', mat_vol_kwargs, micro_kwargs
            )
            self.step2_activation(
                timesteps, source_rates, timestep_units,
                output_dir / 'activation', operator_kwargs=operator_kwargs
            )
            self.step3_photon_source(
                photon_time_indices, bounding_boxes, output_dir / 'photon_transport',
                mat_vol_kwargs=mat_vol_kwargs,
            )
            self.step4_photon_transport(
                output_dir / 'photon_transport', run_kwargs=run_kwargs,
                by_parent_nuclide=by_parent_nuclide,
            )

        return output_dir

    def step1_neutron_transport(
        self,
        output_dir: PathLike = "neutron_transport",
        mat_vol_kwargs: dict | None = None,
        micro_kwargs: dict | None = None
    ):
        """Run the neutron transport step.

        This step computes the material volume fractions on each mesh, creates
        mesh-material filters, and retrieves the fluxes and microscopic cross
        sections for each mesh/material combination via a single transport
        solve. This step will populate the 'fluxes' and 'micros' keys in the
        results dictionary. Any tallies assigned to ``neutron_model`` are run
        in the same transport solve and their results are stored under the
        'neutron_tallies' key. For a mesh-based calculation, this step will
        also populate the 'mesh_material_volumes' key (a list of
        :class:`~openmc.MeshMaterialVolumes`, one per mesh).

        Parameters
        ----------
        output_dir : PathLike, optional
            The directory where the results will be saved.
        mat_vol_kwargs : dict, optional
            Additional keyword arguments based to
            :meth:`openmc.MeshBase.material_volumes`.
        micro_kwargs : dict, optional
            Additional keyword arguments passed to
            :func:`openmc.deplete.get_microxs_and_flux`.

        """

        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        if self.method == 'mesh-based':
            # Compute material volume fractions on each mesh
            if mat_vol_kwargs is None:
                mat_vol_kwargs = {}
            mat_vol_kwargs.setdefault('bounding_boxes', True)

            mmv_list = []
            domain_filters = []
            for i, mesh in enumerate(self.domains):
                mmv = comm.bcast(
                    mesh.material_volumes(self.neutron_model, **mat_vol_kwargs))
                mmv_list.append(mmv)

                # Save results to file
                if comm.rank == 0:
                    mmv.save(output_dir / f'mesh_material_volumes_{i}.npz')

                # Create mesh-material filter for this mesh
                domain_filters.append(
                    openmc.MeshMaterialFilter.from_volumes(mesh, mmv))

            self.results['mesh_material_volumes'] = mmv_list
            domains = domain_filters
        else:
            domains: Sequence[openmc.Cell] = self.domains

            # Check to make sure that each cell is filled with a material and
            # that the volume has been set

            # TODO: If volumes are not set, run volume calculation for cells
            for cell in domains:
                if cell.fill is None:
                    raise ValueError(
                        f"Cell {cell.id} is not filled with a materials. "
                        "Please set the fill material for each cell before "
                        "running the R2S calculation."
                    )
                if cell.volume is None:
                    raise ValueError(
                        f"Cell {cell.id} does not have a volume set. "
                        "Please set the volume for each cell before running "
                        "the R2S calculation."
                    )

        # Set default keyword arguments for microxs and flux calculation
        if micro_kwargs is None:
            micro_kwargs = {}
        micro_kwargs.setdefault('include_model_tallies', True)
        micro_kwargs.setdefault('path_statepoint', output_dir / 'statepoint.h5')
        micro_kwargs.setdefault('path_input', output_dir / 'model.xml')
        statepoint_path = Path(micro_kwargs['path_statepoint']).resolve()
        micro_kwargs['path_statepoint'] = statepoint_path

        # Determine neutron tally IDs based on whether model tallies are included
        neutron_tally_ids = []
        if micro_kwargs['include_model_tallies']:
            neutron_tally_ids = [tally.id for tally in self.neutron_model.tallies]

        # Run neutron transport and get fluxes and micros. Run via openmc.lib to
        # maintain a consistent parallelism strategy with the activation step.
        with TemporarySession(output=False):
            self.results['fluxes'], self.results['micros'] = get_microxs_and_flux(
                self.neutron_model, domains, **micro_kwargs)

        with openmc.StatePoint(statepoint_path) as sp:
            self.results['neutron_tallies'] = [
                sp.tallies[tally_id] for tally_id in neutron_tally_ids
            ]

        # Save the IDs separately from the statepoint so that neutron tally
        # results can be restored by load_results
        if comm.rank == 0:
            with open(output_dir / 'tally_ids.json', 'w') as f:
                json.dump(neutron_tally_ids, f)

        # Save flux and micros to file
        if comm.rank == 0:
            np.save(output_dir / 'fluxes.npy', self.results['fluxes'])
            write_microxs_hdf5(self.results['micros'], output_dir / 'micros.h5')

    def step2_activation(
        self,
        timesteps: Sequence[float] | Sequence[tuple[float, str]],
        source_rates: float | Sequence[float],
        timestep_units: str = 's',
        output_dir: PathLike = 'activation',
        operator_kwargs: dict | None = None,
    ):
        """Run the activation step.

        This step creates a unique copy of each activation material based on the
        mesh elements or cells, then solves the depletion equations for each
        material using the fluxes and microscopic cross sections obtained in the
        neutron transport step. This step will populate the 'depletion_results'
        and 'activation_materials' keys in the results dictionary.

        Parameters
        ----------
        timesteps : Sequence[float] or Sequence[tuple[float, str]]
            Sequence of timesteps. Note that values are not cumulative. The
            units are specified by the `timestep_units` argument when
            `timesteps` is an iterable of float. Alternatively, units can be
            specified for each step by passing an iterable of (value, unit)
            tuples.
        source_rates : float | Sequence[float]
            Source rate in [neutron/sec] for each interval in `timesteps`.
        timestep_units : {'s', 'min', 'h', 'd', 'a'}, optional
            Units for values specified in the `timesteps` argument when passing
            float values. 's' means seconds, 'min' means minutes, 'h' means
            hours, 'd' means days, and 'a' means years (Julian).
        output_dir : PathLike, optional
            Path to directory where activation calculation outputs will be
            saved.
        operator_kwargs : dict, optional
            Additional keyword arguments passed to
            :class:`openmc.deplete.IndependentOperator`.
        """

        if self.method == 'mesh-based':
            # Get unique material for each (mesh, material) combination
            mmv_list = self.results['mesh_material_volumes']
            self.results['activation_materials'] = get_activation_materials(
                self.neutron_model, mmv_list)
        else:
            # Create unique material for each cell
            activation_mats = openmc.Materials()
            for cell in self.domains:
                mat = cell.fill.clone()
                mat.name = f'Cell {cell.id}'
                mat.depletable = True
                mat.volume = cell.volume
                activation_mats.append(mat)
            self.results['activation_materials'] = activation_mats

        # Save activation materials to file
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        self.results['activation_materials'].export_to_xml(
            output_dir / 'materials.xml')

        # Create depletion operator for the activation materials
        if operator_kwargs is None:
            operator_kwargs = {}
        operator_kwargs.setdefault('normalization_mode', 'source-rate')
        op = IndependentOperator(
            self.results['activation_materials'],
            self.results['fluxes'],
            self.results['micros'],
            **operator_kwargs
        )

        # Create time integrator and solve depletion equations
        integrator = PredictorIntegrator(
            op, timesteps, source_rates=source_rates, timestep_units=timestep_units
        )
        output_path = output_dir / 'depletion_results.h5'
        integrator.integrate(final_step=False, path=output_path)
        comm.barrier()

        # Get depletion results
        self.results['depletion_results'] = Results(output_path)

    def step3_photon_source(
        self,
        time_indices: Sequence[int] | None = None,
        bounding_boxes: dict[int, openmc.BoundingBox] | None = None,
        output_dir: PathLike = 'photon_transport',
        mat_vol_kwargs: dict | None = None,
    ):
        """Create decay photon sources.

        This step creates decay photon sources from the activated materials for
        each specified time. In mesh-based mode, the sources are created using
        mesh material volumes, while in cell-based mode, they are created using
        bounding boxes for each cell. This step will populate the
        'photon_sources' key in the results dictionary.

        Parameters
        ----------
        time_indices : Sequence[int], optional
            Sequence of time indices at which photon transport should be run;
            represented as indices into the array of times formed by the
            timesteps. For example, if two timesteps are specified, the array of
            times would contain three entries, and [2] would indicate computing
            photon results at the last time. A value of None indicates to run
            photon transport at each time that has a decay photon source.
        bounding_boxes : dict[int, openmc.BoundingBox], optional
            Dictionary mapping cell IDs to bounding boxes used for spatial
            source sampling in cell-based R2S calculations. Required if method
            is 'cell-based'.
        output_dir : PathLike, optional
            Path to directory where photon source outputs will be saved.
        mat_vol_kwargs : dict, optional
            Additional keyword arguments passed to
            :meth:`openmc.MeshBase.material_volumes`.
        """

        # Do not retain sources from an earlier successful call if this source
        # preparation attempt fails.
        self.results.pop('photon_sources', None)

        # TODO: Automatically determine bounding box for each cell
        if bounding_boxes is None and self.method == 'cell-based':
            raise ValueError("bounding_boxes must be provided for cell-based "
                             "R2S calculations.")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Determine and validate time indices before preparing source data.
        n_steps = len(self.results['depletion_results'])
        implicit_time_indices = time_indices is None
        if time_indices is None:
            time_indices = list(range(n_steps))
        else:
            time_indices = list(time_indices)
            if not time_indices:
                raise ValueError('time_indices must contain at least one index')

            normalized_indices = []
            for index in time_indices:
                if isinstance(index, bool) or not isinstance(index, Integral):
                    raise TypeError('time_indices must contain only integers')
                index = int(index)
                if index < -n_steps or index >= n_steps:
                    raise IndexError(
                        f'Photon time index {index} is out of range for '
                        f'{n_steps} depletion results')
                normalized_index = index % n_steps
                normalized_indices.append(normalized_index)

            # Remove duplicates while preserving order
            time_indices = list(dict.fromkeys(normalized_indices))

        # Check whether the photon model is different
        neutron_univ = self.neutron_model.geometry.root_universe
        photon_univ = self.photon_model.geometry.root_universe
        different_photon_model = (neutron_univ != photon_univ)

        # For mesh-based calculations, compute material volume fractions for the
        # photon model if it is different from the neutron model to account for
        # potential material changes
        if self.method == 'mesh-based' and different_photon_model:
            if mat_vol_kwargs is None:
                mat_vol_kwargs = {}
            photon_mmv_list = []
            for i, mesh in enumerate(self.domains):
                photon_mmv = comm.bcast(
                    mesh.material_volumes(self.photon_model, **mat_vol_kwargs))
                photon_mmv_list.append(photon_mmv)

                # Save photon MMV results to file
                if comm.rank == 0:
                    photon_mmv.save(
                        output_dir / f'mesh_material_volumes_{i}.npz')

            self.results['mesh_material_volumes_photon'] = photon_mmv_list

        # Get dictionary of cells in the photon model
        if different_photon_model:
            photon_cells = self.photon_model.geometry.get_all_cells()

        # Determine eligible work items upfront (independent of time index).
        if self.method == 'mesh-based':
            work_items = self._get_mesh_work_items()
        else:
            work_items = []
            for cell, original_mat in zip(
                    self.domains, self.results['activation_materials']):
                if different_photon_model:
                    if cell.id not in photon_cells or \
                            cell.fill.id != photon_cells[cell.id].fill.id:
                        continue
                work_items.append((cell, original_mat, bounding_boxes[cell.id]))

        # Create decay photon sources for each time index
        photon_sources = {
            time_index: self._create_photon_sources(time_index, work_items)
            for time_index in time_indices
        }

        # Determine if any times have no decay photon sources. If the user
        # didn't specify any specific time indices, remove those times from the
        # photon_sources dictionary. If the user did specify time indices, raise
        # an error if any of those times have no decay photon sources.
        empty_indices = [
            time_index for time_index, sources in photon_sources.items()
            if not sources
        ]
        if implicit_time_indices:
            for time_index in empty_indices:
                del photon_sources[time_index]
            if not photon_sources:
                raise RuntimeError(
                    'No decay photon sources were found at any depletion time')
        elif empty_indices:
            indices = ', '.join(str(index) for index in empty_indices)
            raise RuntimeError(
                f'No decay photon source was found for requested time '
                f'indices: {indices}')

        self.results['photon_sources'] = photon_sources

    def step4_photon_transport(
        self,
        output_dir: PathLike = 'photon_transport',
        run_kwargs: dict | None = None,
        by_parent_nuclide: bool = False,
    ):
        """Run photon transport using prepared decay photon sources.

        This step runs a photon transport calculation for each source list
        created by :meth:`step3_photon_source`. It will populate the
        'photon_tallies' key in the results dictionary.

        Parameters
        ----------
        output_dir : PathLike, optional
            Path to directory where photon transport outputs will be saved.
        run_kwargs : dict, optional
            Additional keyword arguments passed to :meth:`openmc.Model.run`.
            By default, output is disabled.
        by_parent_nuclide : bool, optional
            Whether to score photon tallies separately for each parent
            radionuclide. A :class:`~openmc.ParentNuclideFilter` is added to
            tallies that do not already contain one, with bins determined from
            the prepared decay photon sources.
        """
        if 'photon_sources' not in self.results:
            raise RuntimeError(
                'Photon sources must be created with step3_photon_source '
                'before running photon transport.')
        photon_sources = self.results['photon_sources']
        if not photon_sources:
            raise RuntimeError(
                'No decay photon sources are available for transport')

        if by_parent_nuclide:
            radionuclides = sorted({
                nuclide
                for sources in photon_sources.values()
                for source in sources
                for nuclide in source.energy.nuclides
            })

            if radionuclides:
                parent_filter = openmc.ParentNuclideFilter(radionuclides)
                for tally in self.photon_model.tallies:
                    if not tally.contains_filter(openmc.ParentNuclideFilter):
                        tally.filters.append(parent_filter)

        if run_kwargs is None:
            run_kwargs = {}
        run_kwargs.setdefault('output', False)

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if comm.rank == 0:
            tally_ids = [tally.id for tally in self.photon_model.tallies]
            with open(output_dir / 'tally_ids.json', 'w') as f:
                json.dump(tally_ids, f)

        self.results['photon_tallies'] = {}

        # Ensure photon transport is enabled in settings.
        self.photon_model.settings.photon_transport = True

        for time_index, sources in photon_sources.items():
            self.photon_model.settings.source = sources

            # Run photon transport calculation
            photon_dir = Path(output_dir) / f'time_{time_index}'
            with TemporarySession(self.photon_model, cwd=photon_dir):
                statepoint_path = self.photon_model.run(**run_kwargs)

            # Store tally results
            with openmc.StatePoint(statepoint_path) as sp:
                self.results['photon_tallies'][time_index] = [
                    sp.tallies[tally.id] for tally in self.photon_model.tallies
                ]

    def _get_mesh_work_items(self):
        """Enumerate mesh-based work items across all meshes.

        Returns a list of (index_mat, mat_id, bbox) tuples for each eligible
        mesh element--material combination, where index_mat is the index into
        the activation materials list, mat_id is the material ID, and bbox is
        the bounding box for that mesh element--material combination.

        Returns
        -------
        list of tuple
            Each tuple is (index_mat, mat_id, bbox).
        """
        mmv_list = self.results['mesh_material_volumes']
        photon_mmv_list = self.results.get('mesh_material_volumes_photon')

        work_items = []
        index_mat = 0
        for mesh_idx, mat_vols in enumerate(mmv_list):
            photon_mat_vols = photon_mmv_list[mesh_idx] \
                if photon_mmv_list is not None else None

            n_elements = mat_vols.num_elements
            for index_elem in range(n_elements):
                if photon_mat_vols is not None:
                    photon_materials = {
                        mat_id
                        for mat_id, _ in photon_mat_vols.by_element(index_elem)
                        if mat_id is not None
                    }

                for mat_id, _, bbox in mat_vols.by_element(
                        index_elem, include_bboxes=True):
                    if mat_id is None:
                        continue
                    if photon_mat_vols is not None \
                            and mat_id not in photon_materials:
                        index_mat += 1
                        continue
                    work_items.append((index_mat, mat_id, bbox))
                    index_mat += 1

        return work_items

    def _create_photon_sources(self, time_index, work_items):
        """Create decay photon sources for a set of regions.

        Builds :class:`openmc.IndependentSource` objects with
        :class:`openmc.stats.DecaySpectrum` energy distributions that will be
        serialized to XML and resolved against the depletion chain by the C++
        solver.

        Parameters
        ----------
        time_index : int
            Index into depletion results.
        work_items : list of tuple
            For mesh-based: list of (index_mat, mat_id, bbox).
            For cell-based: list of (cell, original_mat, bbox).

        Returns
        -------
        list of openmc.IndependentSource
            Photon sources for each activated region.
        """
        step_result = self.results['depletion_results'][time_index]
        materials = self.results['activation_materials']
        mesh_based = self.method == 'mesh-based'
        if mesh_based:
            mat_dict = self.neutron_model._get_all_materials()

        sources = []
        for item in work_items:
            if mesh_based:
                index_mat, domain_id, bbox = item
                original_mat = materials[index_mat]
                domain = mat_dict[domain_id]
            else:
                cell, original_mat, bbox = item
                domain = cell

            activated_mat = step_result.get_material(str(original_mat.id))
            nuclides = activated_mat.get_nuclide_atom_densities()
            if not nuclides:
                continue

            # Eliminate nuclides with zero density
            nuclides = {nuclide: density for nuclide, density in nuclides.items()
                        if density > 0}

            energy = openmc.stats.DecaySpectrum(nuclides, activated_mat.volume)
            energy.clip(inplace=True)
            if not energy.nuclides:
                continue

            sources.append(openmc.IndependentSource(
                space=openmc.stats.Box(bbox.lower_left, bbox.upper_right),
                energy=energy,
                particle='photon',
                constraints={'domains': [domain]},
            ))

        return sources

    def load_results(self, path: PathLike):
        """Load results from a previous R2S calculation.

        Parameters
        ----------
        path : PathLike
            Path to the directory containing the R2S calculation results.

        """
        path = Path(path)

        # Load neutron transport results
        neutron_dir = path / 'neutron_transport'
        if self.method == 'mesh-based':
            mmv_files = sorted(neutron_dir.glob('mesh_material_volumes*.npz'),
                               key=lambda p: int(p.stem.split('_')[-1])
                               if p.stem[-1].isdigit() else 0)
            if mmv_files:
                self.results['mesh_material_volumes'] = [
                    openmc.MeshMaterialVolumes.from_npz(f) for f in mmv_files
                ]
        fluxes_file = neutron_dir / 'fluxes.npy'
        if fluxes_file.exists():
            self.results['fluxes'] = list(np.load(fluxes_file, allow_pickle=True))
            micros_dict = read_microxs_hdf5(neutron_dir / 'micros.h5')
            self.results['micros'] = [
                micros_dict[f'domain_{i}'] for i in range(len(micros_dict))
            ]

        # Load neutron tally results from the neutron transport statepoint
        tally_ids_path = neutron_dir / 'tally_ids.json'
        statepoint_path = neutron_dir / 'statepoint.h5'
        if tally_ids_path.exists() and statepoint_path.exists():
            with tally_ids_path.open('r') as f:
                tally_ids = json.load(f)
            with openmc.StatePoint(statepoint_path) as sp:
                self.results['neutron_tallies'] = [
                    sp.tallies[tally_id] for tally_id in tally_ids
                ]

        # Load activation results
        activation_dir = path / 'activation'
        activation_results = activation_dir / 'depletion_results.h5'
        if activation_results.exists():
            self.results['depletion_results'] = Results(activation_results)
        activation_mats_file = activation_dir / 'materials.xml'
        if activation_mats_file.exists():
            self.results['activation_materials'] = \
                openmc.Materials.from_xml(activation_mats_file)

        # Load photon transport results
        photon_dir = path / 'photon_transport'

        # Load photon mesh material volumes if they exist (for mesh-based calculations)
        if self.method == 'mesh-based':
            photon_mmv_files = sorted(
                photon_dir.glob('mesh_material_volumes*.npz'),
                key=lambda p: int(p.stem.split('_')[-1])
                if p.stem[-1].isdigit() else 0)
            if photon_mmv_files:
                self.results['mesh_material_volumes_photon'] = [
                    openmc.MeshMaterialVolumes.from_npz(f)
                    for f in photon_mmv_files
                ]

        # Load tally IDs from JSON file
        tally_ids_path = photon_dir / 'tally_ids.json'
        if tally_ids_path.exists():
            with tally_ids_path.open('r') as f:
                tally_ids = json.load(f)
            self.results['photon_tallies'] = {}

            # For each photon transport calc, load the statepoint and get the
            # tally results based on tally_ids
            for time_dir in photon_dir.glob('time_*'):
                time_index = int(time_dir.name.split('_')[1])
                for sp_path in time_dir.glob('statepoint.*.h5'):
                    with openmc.StatePoint(sp_path) as sp:
                        self.results['photon_tallies'][time_index] = [
                            sp.tallies[tally_id] for tally_id in tally_ids
                        ]
