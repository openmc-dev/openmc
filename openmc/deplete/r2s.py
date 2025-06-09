from __future__ import annotations

import numpy as np
import openmc
from .results import Results


def get_material_copies(model: openmc.Model, mmv: openmc.MeshMaterialVolumes):
    mat_ids = mmv._materials[mmv._materials > -1]
    volumes = mmv._volumes[mmv._materials > -1]
    elems, _ = np.where(mmv._materials > -1)
    # TODO: Handle DAGMC case
    material_dict = model.geometry.get_all_materials()
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
def get_mesh_source(
    model: openmc.Model,
    mesh: openmc.MeshBase,
    materials: list[openmc.Material],
    results: Results,
    mat_vols: openmc.MeshMaterialVolumes,
) -> list[openmc.MeshSource]:
    """Create decay photon source for a mesh-based R2S calculation.

    This function creates N :class:`MeshSource` objects where N is the maximum
    number of unique materials that appears in a single mesh element."""
    mat_dict = model.geometry.get_all_materials()

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
