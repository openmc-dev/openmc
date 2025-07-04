from __future__ import annotations
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path

import numpy as np
import openmc
from . import get_microxs_and_flux, IndependentOperator, PredictorIntegrator
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

    """
    def __init__(
        self,
        model: openmc.Model,
        domains: openmc.MeshBase,
    ):
        self.model = model
        if isinstance(domains, openmc.MeshBase):
            self.method = 'mesh-based'
        else:
            self.method = 'cell-based'
        self.domains = domains
        # TODO: Photon settings
        # TODO: Think about MPI
        # TODO: Option to use CoupledOperator? Needed for true burnup
        # TODO: Voiding/changing materials between neutron/photon steps
        self.results = {}

    def run(
        self,
        timesteps: Sequence[float] | Sequence[tuple[float, str]],
        source_rates: float | Sequence[float],
        cooling_times: Sequence[int],
        dose_tallies: Sequence[openmc.Tally] | None = None,
        timestep_units: str = 's',
        output_dir: PathLike | None = None,
        bounding_boxes: dict[int, openmc.BoundingBox] | None = None,
        micro_kwargs: dict | None = None,
        mat_vol_kwargs: dict | None = None,
        run_kwargs: dict | None = None,
    ):
        """Run the R2S calculation."""

        if output_dir is None:
            stamp = datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
            output_dir = Path(f'r2s_{stamp}')

        self.step1_neutron_transport(
            output_dir / 'neutron_transport', mat_vol_kwargs, micro_kwargs
        )
        self.step2_activation(
            timesteps, source_rates, timestep_units, output_dir / 'activation'
        )
        self.step3_photon_transport(
            cooling_times, dose_tallies, bounding_boxes,
            output_dir / 'photon_transport', run_kwargs=run_kwargs
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

        """

        # TODO: Energy discretization
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if self.method == 'mesh-based':
            # Compute material volume fractions on the mesh
            self.results['mesh_material_volumes'] = mmv = \
                self.domains.material_volumes(self.model, **mat_vol_kwargs)

            # Save results to file
            mmv.save(output_dir / 'mesh_material_volumes.h5')

            # Create mesh-material filter based on what combos were found
            domains = openmc.MeshMaterialFilter.from_volumes(self.domains, mmv)
        else:
            # TODO: Run volume calculation for cells
            domains = self.domains

        # Get fluxes and microscopic cross sections on each mesh/material
        if micro_kwargs is None:
            micro_kwargs = {}
        micro_kwargs.setdefault('path_statepoint', output_dir / 'statepoint.h5')

        # Run neutron transport and get fluxes and micros
        self.results['fluxes'], self.results['micros'] = get_microxs_and_flux(
            self.model, domains, **micro_kwargs)

        # Export model to output directory
        self.model.export_to_model_xml(output_dir / 'model.xml')

        # Save flux and micros to file
        np.save(output_dir / 'fluxes.npy', self.results['fluxes'])
        for i, micros in enumerate(self.results['micros']):
            micros.to_csv(output_dir / f'micros_{i}.csv')

    def step2_activation(
        self,
        timesteps: Sequence[float] | Sequence[tuple[float, str]],
        source_rates: float | Sequence[float],
        timestep_units: str = 's',
        output_dir: PathLike = 'activation',
    ):
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
        cooling_times: Sequence[int],
        dose_tallies: Sequence[openmc.Tally] | None = None,
        bounding_boxes: dict[int, openmc.BoundingBox] | None = None,
        output_dir: PathLike = 'photon_transport',
        run_kwargs: dict | None = None,
    ):
        # TODO: Automatically determine bounding box for each cell

        # Set default run arguments if not provided
        if run_kwargs is None:
            run_kwargs = {}
        run_kwargs.setdefault('output', False)

        # Add dose tallies to model if provided
        if dose_tallies is not None:
            self.model.tallies.extend(dose_tallies)

        for time_index in cooling_times:
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
            self.model.run(**run_kwargs)
