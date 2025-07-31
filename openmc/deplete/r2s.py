from __future__ import annotations
from collections.abc import Sequence
from datetime import datetime
import json
from pathlib import Path

import numpy as np
import openmc
from . import IndependentOperator, PredictorIntegrator
from .microxs import get_microxs_and_flux, write_microxs_hdf5, read_microxs_hdf5
from .results import Results
from ..checkvalue import PathLike
from ..mesh import _get_all_materials


def get_activation_materials(
    model: openmc.Model, mmv: openmc.MeshMaterialVolumes
) -> list[openmc.Material]:
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
    mmv : openmc.MeshMaterialVolumes
        The mesh material volumes object containing the materials and their
        volumes for each mesh element.

    Returns
    -------
    list of openmc.Material
        A list of materials, each corresponding to a unique mesh element and
        material combination.

    """
    # Get the material ID, volume, and element index for each element-material
    # combination
    mat_ids = mmv._materials[mmv._materials > -1]
    volumes = mmv._volumes[mmv._materials > -1]
    elems, _ = np.where(mmv._materials > -1)

    # Get all materials in the model
    material_dict = _get_all_materials(model)

    # Create a new activation material for each element-material combination
    materials = []
    for elem, mat_id, vol in zip(elems, mat_ids, volumes):
        mat = material_dict[mat_id]
        new_mat = mat.clone()
        new_mat.depletable = True
        new_mat.name = f'Element {elem}, Material {mat_id}'
        new_mat.volume = vol
        materials.append(new_mat)

    return materials


def get_decay_photon_source_mesh(
    model: openmc.Model,
    mesh: openmc.MeshBase,
    materials: list[openmc.Material],
    results: Results,
    mat_vols: openmc.MeshMaterialVolumes,
    time_index: int = -1,
) -> list[openmc.MeshSource]:
    """Create decay photon source for a mesh-based R2S calculation.

    This function creates N :class:`MeshSource` objects where N is the maximum
    number of unique materials that appears in a single mesh element. For each
    mesh element-material combination, and IndependentSource instance is created
    with a spatial constraint limited the sampled decay photons to the correct
    region.

    Parameters
    ----------
    model : openmc.Model
        The full model containing the geometry and materials.
    mesh : openmc.MeshBase
        The mesh object defining the spatial regions over which activation is
        performed.
    materials : list of openmc.Material
        List of materials that are activated in the model.
    results : openmc.deplete.Results
        The results object containing the depletion results for the materials.
    mat_vols : openmc.MeshMaterialVolumes
        The mesh material volumes object containing the materials and their
        volumes for each mesh element.

    Returns
    -------
    list of openmc.MeshSource
        A list of MeshSource objects, each containing IndependentSource
        instances for the decay photons in the corresponding mesh element.

    """
    mat_dict = _get_all_materials(model)

    # Some MeshSource objects will have empty positions; create a "null source"
    # that is used for this case
    null_source = openmc.IndependentSource(particle='photon', strength=0.0)

    # List to hold sources for each MeshSource (length = N)
    source_lists = []

    # Index in the overall list of activated materials
    index_mat = 0

    # Total number of mesh elements
    n_elements = mat_vols.num_elements

    for index_elem in range(n_elements):
        for j, (mat_id, _) in enumerate(mat_vols.by_element(index_elem)):
            # Skip void volume
            if mat_id is None:
                continue

            # Check whether a new MeshSource object is needed
            if j >= len(source_lists):
                source_lists.append([null_source]*n_elements)

            # Get activated material composition
            original_mat = materials[index_mat]
            activated_mat = results[time_index].get_material(str(original_mat.id))

            # Create decay photon source source
            energy = activated_mat.get_decay_photon_energy()
            source_lists[j][index_elem] = openmc.IndependentSource(
                energy=energy,
                particle='photon',
                strength=energy.integral(),
                constraints={'domains': [mat_dict[mat_id]]}
            )

            # Increment index of activated material
            index_mat += 1

    # Return list of mesh sources
    return [openmc.MeshSource(mesh, sources) for sources in source_lists]


class R2SManager:
    """Manager for Rigorous 2-Step (R2S) method calculations.

    This class is responsible for managing the materials and sources needed for
    mesh-based or cell-based R2S calculations. It provides methods to get
    activation materials and decay photon sources based on the mesh/cells and
    materials in the OpenMC model.

    Parameters
    ----------
    model : openmc.Model
        The OpenMC model containing the geometry and materials.
    domains : openmc.MeshBase or Sequence[openmc.Cell]
        The mesh or a sequence of cells that represent the spatial units over
        which the R2S calculation will be performed.

    Attributes
    ----------
    domains : openmc.MeshBase or Sequence[openmc.Cell]
        The mesh or a sequence of cells that represent the spatial units over
        which the R2S calculation will be performed.
    model : openmc.Model
        The OpenMC model containing the geometry and materials.
    method : {'mesh-based', 'cell-based'}
        Indicates whether the R2S calculation uses mesh elements ('mesh-based')
        as the spatial discetization or a list of a cells ('cell-based').
    results : dict
        A dictionary that stores results from the R2S calculation.

    """
    def __init__(
        self,
        model: openmc.Model,
        domains: openmc.MeshBase | Sequence[openmc.Cell],
    ):
        self.model = model
        if isinstance(domains, openmc.MeshBase):
            self.method = 'mesh-based'
        else:
            self.method = 'cell-based'
        self.domains = domains
        self.results = {}

    def run(
        self,
        timesteps: Sequence[float] | Sequence[tuple[float, str]],
        source_rates: float | Sequence[float],
        timestep_units: str = 's',
        photon_time_indices: Sequence[int] | None = None,
        photon_tallies: Sequence[openmc.Tally] | None = None,
        photon_settings: openmc.Settings | None = None,
        output_dir: PathLike | None = None,
        bounding_boxes: dict[int, openmc.BoundingBox] | None = None,
        micro_kwargs: dict | None = None,
        mat_vol_kwargs: dict | None = None,
        run_kwargs: dict | None = None,
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
            photon transport for each time.
        photon_tallies : Sequence[openmc.Tally], optional
            A sequence of tallies to be used in the photon transport step. If
            None, no tallies are present.
        photon_settings : openmc.Settings, optional
            Custom settings to use for the photon transport calculation. If
            provided, the original model settings will be temporarily replaced
            during the photon transport calculations and restored afterward.
        output_dir : PathLike, optional
            Path to directory where R2S calculation outputs will be saved. If
            not provided, a timestamped directory 'r2s_YYYY-MM-DDTHH-MM-SS' is
            created. Subdirectories will be created for the neutron transport,
            activation, and photon transport steps.
        bounding_boxes : dict[int, openmc.BoundingBox], optional
            Dictionary mapping cell IDs to bounding boxes used for spatial
            source sampling in cell-based R2S calculations. Required if method
            is 'cell-based'.
        micro_kwargs : dict, optional
            Additional keyword arguments passed to
            :func:`openmc.deplete.get_microxs_and_flux` during the neutron
            transport step.
        mat_vol_kwargs : dict, optional
            Additional keyword arguments passed to
            :meth:`openmc.MeshBase.material_volumes` during the neutron
            transport step.
        run_kwargs : dict, optional
            Additional keyword arguments passed to :meth:`openmc.Model.run`
            during the neutron and photon transport step. By default, output is
            disabled.

        Returns
        -------
        Path
            Path to the output directory containing all calculation results
        """

        if output_dir is None:
            stamp = datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
            output_dir = Path(f'r2s_{stamp}')

        # Set run_kwargs for the neutron transport step
        if micro_kwargs is None:
            micro_kwargs = {}
        if run_kwargs is None:
            run_kwargs = {}
        run_kwargs.setdefault('output', False)
        micro_kwargs.setdefault('run_kwargs', run_kwargs)

        self.step1_neutron_transport(
            output_dir / 'neutron_transport', mat_vol_kwargs, micro_kwargs
        )
        self.step2_activation(
            timesteps, source_rates, timestep_units, output_dir / 'activation'
        )
        self.step3_photon_transport(
            photon_time_indices, photon_tallies, bounding_boxes,
            output_dir / 'photon_transport', run_kwargs=run_kwargs,
            settings=photon_settings
        )

        return output_dir

    def step1_neutron_transport(
        self,
        output_dir: PathLike = "neutron_transport",
        mat_vol_kwargs: dict | None = None,
        micro_kwargs: dict | None = None
    ):
        """Run the neutron transport step.

        This step computes the material volume fractions on the mesh, creates a
        mesh-material filter, and retrieves the fluxes and microscopic cross
        sections for each mesh/material combination. This step will populate the
        'fluxes' and 'micros' keys in the results dictionary. For a mesh-based
        calculation, it will also populate the 'mesh_material_volumes' key.

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

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if self.method == 'mesh-based':
            # Compute material volume fractions on the mesh
            self.results['mesh_material_volumes'] = mmv = \
                self.domains.material_volumes(self.model, **mat_vol_kwargs)

            # Save results to file
            mmv.save(output_dir / 'mesh_material_volumes.npz')

            # Create mesh-material filter based on what combos were found
            domains = openmc.MeshMaterialFilter.from_volumes(self.domains, mmv)
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
        micro_kwargs.setdefault('path_statepoint', output_dir / 'statepoint.h5')
        micro_kwargs.setdefault('path_input', output_dir / 'model.xml')

        # Run neutron transport and get fluxes and micros
        self.results['fluxes'], self.results['micros'] = get_microxs_and_flux(
            self.model, domains, **micro_kwargs)

        # Save flux and micros to file
        np.save(output_dir / 'fluxes.npy', self.results['fluxes'])
        write_microxs_hdf5(self.results['micros'], output_dir / 'micros.h5')

    def step2_activation(
        self,
        timesteps: Sequence[float] | Sequence[tuple[float, str]],
        source_rates: float | Sequence[float],
        timestep_units: str = 's',
        output_dir: PathLike = 'activation',
    ):
        """Run the activation step.

        This step creates a unique copy of each activation material based on the
        mesh elements or cells, then solves the depletion equations for each
        material using the fluxes and microscopic cross sections obtained in the
        neutron transport step.  This step will populate the 'depletion_results'
        key in the results dictionary.

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
        """

        if self.method == 'mesh-based':
            # Get unique material for each (mesh, material) combination
            mmv = self.results['mesh_material_volumes']
            self.results['activation_materials'] = get_activation_materials(self.model, mmv)
        else:
            # Create unique material for each cell
            activation_mats = []
            for cell in self.domains:
                mat = cell.fill.clone()
                mat.name = f'Cell {cell.id}'
                mat.depletable = True
                mat.volume = cell.volume
                activation_mats.append(mat)
            self.results['activation_materials'] = activation_mats

        # Create depletion operator for the activation materials
        op = IndependentOperator(
            self.results['activation_materials'],
            self.results['fluxes'],
            self.results['micros'],
            normalization_mode='source-rate'
        )

        # Create time integrator and solve depletion equations
        integrator = PredictorIntegrator(
            op, timesteps, source_rates=source_rates,
            timestep_units=timestep_units
        )
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / 'depletion_results.h5'
        integrator.integrate(final_step=False, path=output_path)

        # Get depletion results
        self.results['depletion_results'] = Results(output_path)

    def step3_photon_transport(
        self,
        time_indices: Sequence[int] | None = None,
        tallies: Sequence[openmc.Tally] | None = None,
        bounding_boxes: dict[int, openmc.BoundingBox] | None = None,
        output_dir: PathLike = 'photon_transport',
        run_kwargs: dict | None = None,
        settings: openmc.Settings | None = None,
    ):
        """Run the photon transport step.

        This step performs photon transport calculations using decay photon
        sources created from the activated materials. For each specified time,
        it creates appropriate photon sources and runs a transport calculation.
        In mesh-based mode, the sources are created using the mesh material
        volumes, while in cell-based mode, they are created using bounding boxes
        for each cell.  This step will populate the 'photon_tallies' key in the
        results dictionary.

        Parameters
        ----------
        time_indices : Sequence[int], optional
            Sequence of time indices at which photon transport should be run;
            represented as indices into the array of times formed by the
            timesteps. For example, if two timesteps are specified, the array of
            times would contain three entries, and [2] would indicate computing
            photon results at the last time. A value of None indicates to run
            photon transport for each time.
        tallies : Sequence[openmc.Tally], optional
            A sequence of tallies to be used in the photon transport step. If
            None, no tallies are present.
        bounding_boxes : dict[int, openmc.BoundingBox], optional
            Dictionary mapping cell IDs to bounding boxes used for spatial
            source sampling in cell-based R2S calculations. Required if method
            is 'cell-based'.
        output_dir : PathLike, optional
            Path to directory where photon transport outputs will be saved.
        run_kwargs : dict, optional
            Additional keyword arguments passed to :meth:`openmc.Model.run`
            during the photon transport step. By default, output is disabled.
        settings : openmc.Settings, optional
            Custom settings to use for the photon transport calculation. If
            provided, the original model settings will be temporarily replaced
            during the photon transport calculations and restored afterward.
        """

        # TODO: Automatically determine bounding box for each cell
        # TODO: Voiding/changing materials between neutron/photon steps

        # Set default run arguments if not provided
        if run_kwargs is None:
            run_kwargs = {}
        run_kwargs.setdefault('output', False)

        # Add photon tallies to model if provided
        if tallies is None:
            tallies = []
        self.model.tallies = tallies

        # Save original settings to restore later
        original_settings = self.model.settings
        if settings is not None:
            self.model.settings = settings

        # Write out JSON file with tally IDs that can be used for loading
        # results
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        tally_ids = [tally.id for tally in tallies]
        with open(output_dir / 'tally_ids.json', 'w') as f:
            json.dump(tally_ids, f)

        self.results['photon_tallies'] = {}

        try:
            for time_index in time_indices:
                # Create decay photon source
                if self.method == 'mesh-based':
                    self.model.settings.source = get_decay_photon_source_mesh(
                        self.model,
                        self.domains,
                        self.results['activation_materials'],
                        self.results['depletion_results'],
                        self.results['mesh_material_volumes'],
                        time_index=time_index,
                    )
                else:
                    sources = []
                    results = self.results['depletion_results']
                    for cell, original_mat in zip(self.domains, self.results['activation_materials']):
                        bounding_box = bounding_boxes[cell.id]

                        # Get activated material composition
                        activated_mat = results[time_index].get_material(str(original_mat.id))

                        # Create decay photon source source
                        space = openmc.stats.Box(*bounding_box)
                        energy = activated_mat.get_decay_photon_energy()
                        source = openmc.IndependentSource(
                            space=space,
                            energy=energy,
                            particle='photon',
                            strength=energy.integral(),
                            domains=[cell]
                        )
                        sources.append(source)
                    self.model.settings.source = sources

                # Convert time_index (which may be negative) to a normal index
                if time_index < 0:
                    time_index = len(self.results['depletion_results']) + time_index

                # Run photon transport calculation
                run_kwargs['cwd'] = Path(output_dir) / f'time_{time_index}'
                statepoint_path = self.model.run(**run_kwargs)

                # Store tally results
                with openmc.StatePoint(statepoint_path) as sp:
                    self.results['photon_tallies'][time_index] = [
                        sp.tallies[tally.id] for tally in tallies
                    ]
        finally:
            # Restore original settings
            self.model.settings = original_settings

    def load_results(self, path: PathLike):
        """Load results from a previous R2S calculation.

        Parameters
        ----------
        path : PathLike
            Path to the directory containing the R2S calculation results.

        """
        path = Path(path)

        # Load neutron transport results
        if self.method == 'mesh-based':
            neutron_dir = path / 'neutron_transport'
            mmv_file = neutron_dir / 'mesh_material_volumes.npz'
            if mmv_file.exists():
                self.results['mesh_material_volumes'] = \
                    openmc.MeshMaterialVolumes.from_npz(mmv_file)
        fluxes_file = neutron_dir / 'fluxes.npy'
        if fluxes_file.exists():
            self.results['fluxes'] = list(np.load(fluxes_file, allow_pickle=True))
            micros_dict = read_microxs_hdf5(neutron_dir / 'micros.h5')
            self.results['micros'] = [
                micros_dict[f'domain_{i}'] for i in range(len(micros_dict))
            ]

        # Load activation results
        activation_dir = path / 'activation'
        activation_results = activation_dir / 'depletion_results.h5'
        if activation_results.exists():
            self.results['depletion_results'] = Results()

        # Load photon transport results
        photon_dir = path / 'photon_transport'

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
