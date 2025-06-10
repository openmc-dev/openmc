from __future__ import annotations

import numpy as np
import openmc
from .results import Results
from ..mesh import _get_all_materials


def get_activation_materials(
        model: openmc.Model,
        mmv: openmc.MeshMaterialVolumes
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
        The OpenMC model containing the materials.
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


# TODO: put as part of a R2S class
def get_decay_photon_source(
    model: openmc.Model,
    mesh: openmc.MeshBase,
    materials: list[openmc.Material],
    results: Results,
    mat_vols: openmc.MeshMaterialVolumes,
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
        The OpenMC model containing the geometry and materials.
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
            activated_mat = results[-1].get_material(str(original_mat.id))

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
