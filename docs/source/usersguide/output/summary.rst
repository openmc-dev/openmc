.. _usersguide_summary:

===================
Summary File Format
===================

The current revision of the summary file format is 1.

**/filetype** (*char[]*)

    String indicating the type of file.

**/revision** (*int*)

    Revision of the summary file format. Any time a change is made in the
    format, this integer is incremented.

**/version_major** (*int*)

    Major version number for OpenMC

**/version_minor** (*int*)

    Minor version number for OpenMC

**/version_release** (*int*)

    Release version number for OpenMC

**/date_and_time** (*char[]*)

    Date and time the state point was written.

**/n_procs** (*int*)

    Number of MPI processes used.

**/n_particles** (*int8_t*)

    Number of particles used per generation.

**/n_batches** (*int*)

    Number of batches to simulate.

if (run_mode == MODE_EIGENVALUE)

    **/n_inactive** (*int*)

        Number of inactive batches.

    **/n_active** (*int*)

        Number of active batches.

    **/gen_per_batch** (*int*)

        Number of generations per batch.

end if

**/geometry/n_cells** (*int*)

**/geometry/n_surfaces** (*int*)

**/geometry/n_universes** (*int*)

**/geometry/n_lattices** (*int*)

do i = 1, n_cells

    **/geometry/cells/cell <uid>/index** (*int*)

    **/geometry/cells/cell <uid>/name** (*char[]*)

    **/geometry/cells/cell <uid>/universe** (*int*)

    **/geometry/cells/cell <uid>/fill_type** (*char[]*)

    **/geometry/cells/cell <uid>/material** (*int*)

    **/geometry/cells/cell <uid>/maps** (*int*)

    **/geometry/cells/cell <uid>/offset** (*int[]*)

    **/geometry/cells/cell <uid>/translated** (*int*)

    **/geometry/cells/cell <uid>/translation** (*double[]*)

    **/geometry/cells/cell <uid>/rotated** (*int*)

    **/geometry/cells/cell <uid>/rotation** (*double[]*)

    **/geometry/cells/cell <uid>/lattice** (*int*)

    **/geometry/cells/cell <uid>/surfaces** (*int[]*)

end do

do i = 1, n_surfaces

    **/geometry/surfaces/surface <uid>/index** (*int*)

    **/geometry/surfaces/surface <uid>/name** (*char[]*)

    **/geometry/surfaces/surface <uid>/type** (*char[]*)

    **/geometry/surfaces/surface <uid>/coefficients** (*double[]*)

    **/geometry/surfaces/surface <uid>/boundary_condition** (*char[]*)

end do

do i = 1, n_universes

    **/geometry/universes/universe <uid>/index** (*int*)

    **/geometry/universes/universe <uid>/cells** (*int[]*)

end do

do i = 1, n_lattices

    **/geometry/lattices/lattice <uid>/index** (*int*)

    **/geometry/lattices/lattice <uid>/name** (*char[]*)

    **/geometry/lattices/lattice <uid>/type** (*char[]*)

    **/geometry/lattices/lattice <uid>/pitch** (*double[]*)

    **/geometry/lattices/lattice <uid>/outer** (*int*)

    **/geometry/lattices/lattice <uid>/offset_size** (*int[]*)

    **/geometry/lattices/lattice <uid>/maps** (*int*)

    **/geometry/lattices/lattice <uid>/offsets** (*int[]*)

    **/geometry/lattices/lattice <uid>/universes** (*int[]*)

    if (rectangular lattice)

        **/geometry/lattices/lattice <uid>/dimension** (*int[]*)

        **/geometry/lattices/lattice <uid>/lower_left** (*double[]*)

    elseif (hexagonal lattice)

        **/geometry/lattices/lattice <uid>/n_rings** (*int*)

        **/geometry/lattices/lattice <uid>/n_axial** (*int*)

        **/geometry/lattices/lattice <uid>/center** (*double[]*)

    end if

end do

**/n_materials** (*int*)

do i = 1, n_materials

    **/materials/material <uid>/index** (*int*)

    **/materials/material <uid>/name** (*char[]*)

    **/materials/material <uid>/atom_density** (*double[]*)

    **/materials/material <uid>/nuclides** (*int[]*)

    **/materials/material <uid>/nuclide_densities** (*double[]*)

    **/materials/material <uid>/sab_names** (*char[][]*)

end do

**/tallies/n_tallies** (*int*)

**/tallies/n_meshes** (*int*)

do i = 1, n_meshes

   **/tallies/mesh <uid>/index** (*int*)

   **/tallies/mesh <uid>/type** (*char[]*)

   **/tallies/mesh <uid>/dimension** (*int[]*)

   **/tallies/mesh <uid>/lower_left** (*double[]*)

   **/tallies/mesh <uid>/upper_right** (*double[]*)

   **/tallies/mesh <uid>/width** (*double[]*)

end do

do i = 1, n_tallies

    **/tallies/tally <uid>/index** (*int*)

    **/tallies/tally <uid>/name** (*char[]*)

    **/tallies/tally <uid>/total_score_bins** (*int*)

    **/tallies/tally <uid>/total_filter_bins** (*int*)

    **/tallies/tally <uid>/n_filters** (*int*)

    do j = 1, n_filters

        **/tallies/tally <uid>/filter j/type** (*char[]*)

        **/tallies/tally <uid>/filter j/n_bins** (*int*)

        **/tallies/tally <uid>/filter j/bins** (*int[]* or *double[]*)

        **/tallies/tally <uid>/filter j/type_name** (*char[]*)

    end do

    **/tallies/tally <uid>/nuclides** (*char[][]*)

    **/tallies/tally <uid>/n_score_bins** (*int*)

    **/tallies/tally <uid>/score_bins** (*char[][]*)

end do
