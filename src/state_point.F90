module state_point

!===============================================================================
! STATE_POINT -- This module handles writing and reading state point
! files. State points are contain complete tally results, source sites, and
! various other data. They can be used to restart a run or to reconstruct
! confidence intervals for tallies (this requires post-processing via Python
! scripts).
!
! State points can be written at any batch during a simulation, or at specified
! intervals, using the <state_point ... /> tag.
!===============================================================================


  use constants
  use error,              only: fatal_error, warning
  use geometry_header,    only: Cell, Universe, Lattice, Surface
  use global
  use output,             only: write_message, time_stamp
  use source,             only: write_source_bank, read_source_bank
  use string,             only: to_str, zero_padded, count_digits
  use material_header,    only: Material
  use output_interface
  use tally_header,       only: TallyObject
  use tally,              only: write_tally_result, read_tally_result
  use dict_header,        only: ElemKeyValueII

#ifdef HDF5
  use hdf5
  use hdf5_interface,   only: dims1, hdf5_rank, dset, dspace, hdf5_err
#endif

#ifdef MPI
  use mpi
#endif

  implicit none

  type(BinaryOutput)        :: sp      ! Statepoint/source output file
  type(BinaryOutput)        :: fh      ! general file handle

contains

!===============================================================================
! WRITE_STATE_POINT
!===============================================================================

  subroutine write_state_point()

    character(MAX_FILE_LEN)    :: filename
    integer                    :: i, j, k, m, n
    integer                    :: n_x, n_y, n_z
    integer, allocatable       :: temp_array(:)
    integer, allocatable       :: temp_array2(:)
    integer, allocatable       :: lattice_universes(:,:,:)
    type(Cell),     pointer    :: c => null()
    type(Universe), pointer    :: u => null()
    type(Lattice),  pointer    :: lat => null()
    type(Surface),  pointer    :: s => null()
    type(Material), pointer    :: mat => null()
    type(TallyObject), pointer :: t => null()
    type(ElemKeyValueII), pointer :: current => null()
    type(ElemKeyValueII), pointer :: next => null()

    ! Set filename for state point
    filename = trim(path_output) // 'statepoint.' // &
        & zero_padded(current_batch, count_digits(n_batches))

    if (dd_run) then
      filename = trim(filename) // '.domain_' // &
          & zero_padded(domain_decomp % meshbin, &
                        count_digits(domain_decomp % n_domains))
    end if

    ! Append appropriate extension
#ifdef HDF5
    filename = trim(filename) // '.h5'
#else
    filename = trim(filename) // '.binary'
#endif

    ! Write message
    call write_message("Creating state point " // trim(filename) // "...", 1)

    if (master .or. (dd_run .and. domain_decomp % local_master)) then
      ! Create statepoint file
      call sp % file_create(filename)

      ! Write file type
      call sp % write_data(FILETYPE_STATEPOINT, "filetype")

      ! Write revision number for state point file
      call sp % write_data(REVISION_STATEPOINT, "revision")

      ! Write OpenMC version
      call sp % write_data(VERSION_MAJOR, "version_major")
      call sp % write_data(VERSION_MINOR, "version_minor")
      call sp % write_data(VERSION_RELEASE, "version_release")

      ! Write current date and time
      call sp % write_data(time_stamp(), "date_and_time")

      ! Write path to input
      call sp % write_data(path_input, "path")

      ! Write out random number seed
      call sp % write_data(seed, "seed")

      ! Write run information
      call sp % write_data(run_mode, "run_mode")
      call sp % write_data(n_particles, "n_particles")

      ! Write out current batch number
      call sp % write_data(current_batch, "current_batch")

      ! Write out information for eigenvalue run
      if (run_mode == MODE_EIGENVALUE) then
        call sp % write_data(n_inactive, "n_inactive")
        call sp % write_data(gen_per_batch, "gen_per_batch")
        call sp % write_data(k_generation, "k_generation", &
             length=current_batch*gen_per_batch)
        call sp % write_data(entropy, "entropy", length=current_batch*gen_per_batch)
        call sp % write_data(k_col_abs, "k_col_abs")
        call sp % write_data(k_col_tra, "k_col_tra")
        call sp % write_data(k_abs_tra, "k_abs_tra")
        call sp % write_data(k_combined, "k_combined", length=2)

        ! Write out CMFD info
        if (cmfd_on) then
          call sp % write_data(1, "cmfd_on")
          call sp % write_data(cmfd % indices, "indices", length=4, group="cmfd")
          call sp % write_data(cmfd % k_cmfd, "k_cmfd", length=current_batch, &
               group="cmfd")
          call sp % write_data(cmfd % cmfd_src, "cmfd_src", &
               length=(/cmfd % indices(4), cmfd % indices(1), &
               cmfd % indices(2), cmfd % indices(3)/), &
               group="cmfd")
          call sp % write_data(cmfd % entropy, "cmfd_entropy", &
                          length=current_batch, group="cmfd")
          call sp % write_data(cmfd % balance, "cmfd_balance", &
               length=current_batch, group="cmfd")
          call sp % write_data(cmfd % dom, "cmfd_dominance", &
               length = current_batch, group="cmfd")
          call sp % write_data(cmfd % src_cmp, "cmfd_srccmp", &
               length = current_batch, group="cmfd")
        else
          call sp % write_data(0, "cmfd_on")
        end if
      end if

      ! Begin writing geometry/material information
      call sp % write_data(n_cells, "n_cells", group="geometry")
      call sp % write_data(n_universes, "n_universes", group="geometry")
      call sp % write_data(n_lattices, "n_lattices", group="geometry")
      call sp % write_data(n_surfaces, "n_surfaces", group="geometry")
      call sp % write_data(n_materials, "n_materials", group="geometry")

      if (n_lattices > 0) then
        ! Print list of lattice IDs
        allocate(temp_array(n_lattices))
        do i = 1, n_lattices        
          lat => lattices(i)
          temp_array(i) = lat % id
        end do
        call sp % write_data(temp_array, "lattice_ids", &
             group="geometry", length=n_lattices)
        deallocate(temp_array)
      end if

      ! Print list of universe IDs
      allocate(temp_array(n_universes))
      do i = 1, n_universes        
        u => universes(i)
        temp_array(i) = u % id
      end do
      call sp % write_data(temp_array, "universe_ids", &
           group="geometry", length=n_universes)
      deallocate(temp_array)

      ! Print list of surface IDs
      current => surface_dict % keys()
      i = 1
      allocate(temp_array(n_surfaces))
      allocate(temp_array2(n_surfaces))
      do while (associated(current))
        temp_array(i) = current % key
        temp_array2(i) = current % value
        ! Move to next surface
        next => current % next
        deallocate(current)
        current => next
        i = i + 1
      end do
      call sp % write_data(temp_array, "surface_keys", &
           group="geometry", length=n_surfaces)
      call sp % write_data(temp_array2, "surface_ids", &
           group="geometry", length=n_surfaces)
      deallocate(temp_array)
      deallocate(temp_array2)

      ! Print list of material IDs
      allocate(temp_array(n_materials))
      do i = 1, n_materials        
        mat => materials(i)
        temp_array(i) = mat % id
      end do
      call sp % write_data(temp_array, "material_ids", &
           group="geometry", length=n_materials)
      deallocate(temp_array)

      ! Print list of cell keys-> IDs
      current => cell_dict % keys()
      i = 1
      allocate(temp_array(n_cells))
      allocate(temp_array2(n_cells))
      do while (associated(current))
        temp_array(i) = current % key
        temp_array2(i) = current % value
        ! Move to next universe
        next => current % next
        deallocate(current)
        current => next
        i = i + 1
      end do
      call sp % write_data(temp_array, "cell_keys", &
           group="geometry", length=n_cells)
      call sp % write_data(temp_array2, "cell_ids", &
           group="geometry", length=n_cells)
      deallocate(temp_array)
      deallocate(temp_array2)

      ! =======================================================================
      ! WRITE INFORMATION ON CELLS

      ! Create a cell group (nothing directly written in this group) then close
#ifdef HDF5
      call sp % open_group("geometry/cells")
      call sp % close_group()
#endif

      ! Write information on each cell
      CELL_LOOP: do i = 1, n_cells
        c => cells(i)    
        ! Write information on what fills this cell
        call sp % write_data(c % type, "fill_type", &
               group="geometry/cells/cell " // trim(to_str(c % id)))

        call sp % write_data(c % n_surfaces, "n_surfaces", &
               group="geometry/cells/cell " // trim(to_str(c % id)))
        call sp % write_data(c % surfaces, "surfaces", & 
               length= c % n_surfaces, &
               group="geometry/cells/cell " // trim(to_str(c % id)))

        select case (c % type)
        case (CELL_NORMAL)          
          if (c % material == MATERIAL_VOID) then
            call sp % write_data(-1, "material", &
                 group="geometry/cells/cell " // trim(to_str(c % id)))
          else
            call sp % write_data(materials(c % material) % id, "material", &
                 group="geometry/cells/cell " // trim(to_str(c % id)))
          end if
        case (CELL_FILL)
          call sp % write_data(universes(c % fill) % id, "fill", &
               group="geometry/cells/cell " // trim(to_str(c % id))) 
          call sp % write_data(size(c % offset), "maps", &
               group="geometry/cells/cell " // trim(to_str(c % id)))
          if (size(c % offset) > 0) then
            call sp % write_data(c % offset, "offset", length=size(c%offset), &
                 group="geometry/cells/cell " // trim(to_str(c % id)))
          end if
        case (CELL_LATTICE)
          call sp % write_data(lattices(c % fill) % id, "lattice", &
               group="geometry/cells/cell " // trim(to_str(c % id))) 
        end select

      end do CELL_LOOP

      ! =======================================================================
      ! WRITE INFORMATION ON UNIVERSES

      ! Create universes group (nothing directly written here) then close
#ifdef HDF5
      call sp % open_group("geometry/universes")
      call sp % close_group()
#endif

      ! Write information on each universe
      UNIVERSE_LOOP: do i = 1, n_universes
        u => universes(i)
        call sp % write_data(u % n_cells, "n_cells", &
             group="geometry/universes/universe " // trim(to_str(u % id)))

        ! Write list of cells in this universe
        if (u % n_cells > 0) then
          call sp % write_data(u % cells, "cells", length=u % n_cells, &
               group="geometry/universes/universe " // trim(to_str(u % id)))
        end if

      end do UNIVERSE_LOOP

      ! =======================================================================
      ! WRITE INFORMATION ON LATTICES

      ! Create lattices group (nothing directly written here) then close
#ifdef HDF5
      call sp % open_group("geometry/lattices")
      call sp % close_group()
#endif

      ! Write information on each lattice
      LATTICE_LOOP: do i = 1, n_lattices
        lat => lattices(i)

        ! Write lattice type
        select case(lat % type)
        case (LATTICE_RECT)
          call sp % write_data(1, "type", &
               group="geometry/lattices/lattice " // trim(to_str(lat % id)))
        case (LATTICE_HEX)
          call sp % write_data(2, "type", &
               group="geometry/lattices/lattice " // trim(to_str(lat % id)))
        end select

        ! Write lattice dimensions, number of offset maps, and offsets
        if (lat % n_dimension == 2) then
          call sp % write_data((/lat % dimension(1),lat % dimension(2),1/), & 
               "dimension", length=3, &
               group="geometry/lattices/lattice " // trim(to_str(lat % id)))
        else          
         call sp % write_data(lat % dimension, "dimension", &
              length=3, &
              group="geometry/lattices/lattice " // trim(to_str(lat % id)))
        end if

        call sp % write_data(size(lat % offset,1), "maps", &
             group="geometry/lattices/lattice " // trim(to_str(lat % id)))
        call sp % write_data(size(lat % offset), "offset_size", &
             group="geometry/lattices/lattice " // trim(to_str(lat % id)))
        if (size(lat % offset) > 0) then
            call sp % write_data(lat % offset, "offset", &
                 length=shape(lat % offset), &
                 group="geometry/lattices/lattice " // trim(to_str(lat % id)))
        end if

        ! Determine dimensions of lattice
        n_x = lat % dimension(1)
        n_y = lat % dimension(2)
        if (lat % n_dimension == 3) then
          n_z = lat % dimension(3)
        else
          n_z = 1
        end if

        ! Write lattice universes
        allocate(lattice_universes(n_x, n_y, n_z))
        do j = 1, n_x
          do k = 1, n_y
            do m = 1, n_z
              lattice_universes(j,k,m) = universes(lat % universes(j,k,m)) % id
            end do
          end do
        end do
        call sp % write_data(lattice_universes, "universes", &
             length=(/n_x, n_y, n_z/), &
             group="geometry/lattices/lattice " // trim(to_str(lat % id)))
        deallocate(lattice_universes)

      end do LATTICE_LOOP

      ! =======================================================================
      ! WRITE INFORMATION ON SURFACES

      ! Create surfaces group (nothing directly written here) then close
#ifdef HDF5
      call sp % open_group("geometry/surfaces")
      call sp % close_group()
#endif

      ! Write information on each surface
      SURFACE_LOOP: do i = 1, n_surfaces
        s => surfaces(i)
        call sp % write_data(s % type, "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
        call sp % write_data(s % bc, "bc", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))

        call sp % write_data(size(s % coeffs), "n_coeffs", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
        if (size(s % coeffs) > 0) then
          call sp % write_data(s % coeffs, "coeffs", length=size(s % coeffs), &
               group="geometry/surfaces/surface " // trim(to_str(s % id)))
        end if

        if (allocated(s % neighbor_pos)) then
          call sp % write_data(size(s % neighbor_pos), "n_neighbor_pos", &
               group="geometry/surfaces/surface " // trim(to_str(s % id)))
          call sp % write_data(s % neighbor_pos, "neighbor_pos",  &
               length=size(s % neighbor_pos), &
               group="geometry/surfaces/surface " // trim(to_str(s % id)))
        else
          j = 0
          call sp % write_data(j, "n_neighbor_pos", &
               group="geometry/surfaces/surface " // trim(to_str(s % id)))
        endif

        if (allocated(s % neighbor_neg)) then
          call sp % write_data(size(s % neighbor_neg), "n_neighbor_neg", &
               group="geometry/surfaces/surface " // trim(to_str(s % id)))
          call sp % write_data(s % neighbor_neg, "neighbor_neg", &
               length=size(s % neighbor_neg), &
               group="geometry/surfaces/surface " // trim(to_str(s % id)))
        else
          j = 0
          call sp % write_data(j, "n_neighbor_neg", &
               group="geometry/surfaces/surface " // trim(to_str(s % id)))
        endif
      end do SURFACE_LOOP

      ! =======================================================================
      ! WRITE INFORMATION ON MATERIALS

      ! Create materials group (nothing directly written here) then close
#ifdef HDF5
      call sp % open_group("geometry/materials")
      call sp % close_group()
#endif

      ! Write information on each material
      MATERIAL_LOOP: do i = 1, n_materials
        mat => materials(i)
        call sp % write_data(mat % n_nuclides, "n_nuclides", &
             group="geometry/materials/material " // trim(to_str(mat % id)))
        call sp % write_data(mat % nuclide, "nuclide", & 
             length=mat % n_nuclides, &
             group="geometry/materials/material " // trim(to_str(mat % id)))

        j = 0
        if (mat % distrib_dens) then
          j = 1
        end if 
        call sp % write_data(j, "distrib_dens", &
             group="geometry/materials/material " // trim(to_str(mat % id)))

        j = 0
        if (mat % distrib_comp) then
          j = 1
        end if
        call sp % write_data(j, "distrib_comp", &
             group="geometry/materials/material " // trim(to_str(mat % id)))

        call sp % write_data(mat % density % num, "n_density", &
             group="geometry/materials/material " // trim(to_str(mat % id)))
        call sp % write_data(mat % density % density, "density", &
             length=mat % density % num, &
             group="geometry/materials/material " // trim(to_str(mat % id)))

        call sp % write_data(mat % n_comp, "n_comp", &
             group="geometry/materials/material " // trim(to_str(mat % id)))

        j = 0
        if (mat % otf_compositions) then
          j = 1
        end if
        call sp % write_data(j, "otf_compositions", &
             group="geometry/materials/material " // trim(to_str(mat % id)))

        if (.not. mat % otf_compositions) then

#ifdef HDF5
          call sp % open_group("geometry/materials/material "&
                 // trim(to_str(mat % id)) // "/compositions/")
          call sp % close_group()
#endif
          COMPOSITION_LOOP: do j = 1, mat % n_comp
            call sp % write_data(mat % comp(j) % atom_density, "atom_density ", &
                 length=mat % n_nuclides, & 
                 group="geometry/materials/material " &
                 // trim(to_str(mat % id)) // "/compositions/" &
                 // trim(to_str(j)))
          end do COMPOSITION_LOOP

        end if

        call sp % write_data(mat % n_sab, "n_sab", &
             group="geometry/materials/material " // trim(to_str(mat % id)))
        if (mat % n_sab > 0) then
          call sp % write_data(mat % i_sab_nuclides, "i_sab_nuclides", length=mat % n_sab, &
               group="geometry/materials/material " // trim(to_str(mat % id)))
          call sp % write_data(mat % i_sab_tables, "i_sab_tables", length=mat % n_sab, &
               group="geometry/materials/material " // trim(to_str(mat % id)))
        endif

      end do MATERIAL_LOOP

      ! Write number of meshes
      call sp % write_data(n_meshes, "n_meshes", group="tallies")

      ! Write information for meshes
      MESH_LOOP: do i = 1, n_meshes
        call sp % write_data(meshes(i) % id, "id", &
             group="tallies/mesh" // to_str(i))
        call sp % write_data(meshes(i) % type, "type", &
             group="tallies/mesh" // to_str(i))
        call sp % write_data(meshes(i) % n_dimension, "n_dimension", &
             group="tallies/mesh" // to_str(i))
        call sp % write_data(meshes(i) % dimension, "dimension", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
        call sp % write_data(meshes(i) % lower_left, "lower_left", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
        call sp % write_data(meshes(i) % upper_right, "upper_right", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
        call sp % write_data(meshes(i) % width, "width", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
      end do MESH_LOOP

      ! Write number of tallies
      call sp % write_data(n_tallies, "n_tallies", group="tallies")

      ! Write all tally information except results
      TALLY_METADATA: do i = 1, n_tallies
        !Get pointer to tally
        t => tallies(i)

        ! Write id
        call sp % write_data(t % id, "id", group="tallies/tally" // to_str(i))

        ! Write label                       
        call sp % write_data(len(t % label), "labellen", group="tallies/tally" // to_str(i))
        if (len(t % label) > 0) then
          call sp % write_data(t % label, "label", group="tallies/tally" // to_str(i))
        endif

        ! Write estimator type                                                                      
        call sp % write_data(t % estimator, "estimator", group="tallies/tally" // to_str(i))

        ! Write number of realizations
        call sp % write_data(t % n_realizations, "n_realizations", &
             group="tallies/tally" // to_str(i))

        ! Write size of each tally
        call sp % write_data(t % total_score_bins, "total_score_bins", &
             group="tallies/tally" // to_str(i))
        call sp % write_data(t % total_filter_bins, "total_filter_bins", &
             group="tallies/tally" // to_str(i))

        ! Write on-the-fly allocation tally info
        if (t % on_the_fly_allocation) then
          n = t % next_filter_idx - 1
          call sp % write_data(n, "otf_size_results_filters", &
               group="tallies/tally" // to_str(i))
          ! Write otf filter bin mapping
          allocate(temp_array(n))
          do j = 1, n
            temp_array(j) = t % reverse_filter_index_map % get_key(j)
          end do
          call sp % write_data(temp_array, "otf_filter_bin_map", &
               group="tallies/tally" // to_str(i), &
               length=n)
          deallocate(temp_array)
        else
          call sp % write_data(NONE, "otf_size_results_filters", &
               group="tallies/tally" // to_str(i))
        end if

        ! Write number of filters
        call sp % write_data(t % n_filters, "n_filters", &
             group="tallies/tally" // to_str(i))

        ! Write filter information
        FILTER_LOOP: do j = 1, t % n_filters

          ! Write type of filter
          call sp % write_data(t % filters(j) % type, "type", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j))
          ! Write offset for this filter
          call sp % write_data(t % filters(j) % offset, "offset", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j))

          ! Write number of bins for this filter
          call sp % write_data(t % filters(j) % n_bins, "n_bins", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j))

          ! Write bins
          if (t % filters(j) % type == FILTER_ENERGYIN .or. &
              t % filters(j) % type == FILTER_ENERGYOUT) then
            call sp % write_data(t % filters(j) % real_bins, "bins", &
                 group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
                 length=size(t % filters(j) % real_bins))
          else
            call sp % write_data(t % filters(j) % int_bins, "bins", &
                 group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
                 length=size(t % filters(j) % int_bins))
          end if

        end do FILTER_LOOP

        ! Write number of nuclide bins
        call sp % write_data(t % n_nuclide_bins, "n_nuclide_bins", &
             group="tallies/tally" // to_str(i))

        ! Set up nuclide bin array and then write
        allocate(temp_array(t % n_nuclide_bins))
        NUCLIDE_LOOP: do j = 1, t % n_nuclide_bins
          if (t % nuclide_bins(j) > 0) then
            temp_array(j) = nuclides(t % nuclide_bins(j)) % zaid
          else
            temp_array(j) = t % nuclide_bins(j)
          end if
        end do NUCLIDE_LOOP
        call sp % write_data(temp_array, "nuclide_bins", &
             group="tallies/tally" // to_str(i), length=t % n_nuclide_bins)
        deallocate(temp_array)

        ! Write number of score bins, score bins, and moment order
        call sp % write_data(t % n_score_bins, "n_score_bins", &
             group="tallies/tally" // to_str(i))
        call sp % write_data(t % score_bins, "score_bins", &
             group="tallies/tally" // to_str(i), length=t % n_score_bins)
        call sp % write_data(t % moment_order, "moment_order", &
             group="tallies/tally" // to_str(i), length=t % n_score_bins)

        ! Write number of user score bins
        call sp % write_data(t % n_user_score_bins, "n_user_score_bins", &
             group="tallies/tally" // to_str(i))

      end do TALLY_METADATA

      ! Indicate where source bank is stored in statepoint
      if (source_separate) then
        call sp % write_data(0, "source_present")
      else
        call sp % write_data(1, "source_present")
      end if

    end if

    ! Check for the no-tally-reduction method
    if (.not. reduce_tallies) then
      ! If using the no-tally-reduction method, we need to collect tally
      ! results before writing them to the state point file.

      if (dd_run) then
        call fatal_error('no_reduce not implemented with domain decomposition')
      end if

      call write_tally_results_nr()

    elseif (master .or. (dd_run .and. domain_decomp % local_master)) then

      ! Write number of global realizations
      call sp % write_data(n_realizations, "n_realizations")

      ! Write global tallies
      call sp % write_data(N_GLOBAL_TALLIES, "n_global_tallies")
      call write_tally_result(sp, global_tallies, "global_tallies", &
           n1=N_GLOBAL_TALLIES, n2=1)

      ! Write tallies
      if (tallies_on) then

        ! Indicate that tallies are on
        call sp % write_data(1, "tallies_present", group="tallies")

        ! Write all tally results
        TALLY_RESULTS: do i = 1, n_tallies

          ! Set point to current tally
          t => tallies(i)

          ! Write sum and sum_sq for each bin
          if (t % on_the_fly_allocation) then

            n = t % next_filter_idx - 1
            call write_tally_result(sp, t % results(:,1:n), "results", &
                 group="tallies/tally" // to_str(i), &
                 n1=size(t % results, 1), n2=n)
          else

            call write_tally_result(sp, t % results, "results", &
                 group="tallies/tally" // to_str(i), &
                 n1=size(t % results, 1), n2=size(t % results, 2))

          endif
        end do TALLY_RESULTS

      else

        ! Indicate tallies are off
        call sp % write_data(0, "tallies_present", group="tallies")

      end if

      ! Close the file for serial writing
      call sp % file_close()

    end if

  end subroutine write_state_point

!===============================================================================
! WRITE_SOURCE_POINT
!===============================================================================

  subroutine write_source_point()

    type(BinaryOutput) :: sp
    character(MAX_FILE_LEN) :: filename

    ! Check to write out source for a specified batch
    if (sourcepoint_batch % contains(current_batch)) then

      ! Create or open up file
      if (source_separate) then

        ! Set filename
        filename = trim(path_output) // 'source.' // &
            & zero_padded(current_batch, count_digits(n_batches))

        if (dd_run) then
          filename = trim(filename) // '.domain_' // &
              & zero_padded(domain_decomp % meshbin, &
                            count_digits(domain_decomp % n_domains))
        end if

#ifdef HDF5
        filename = trim(filename) // '.h5'
#else
        filename = trim(filename) // '.binary'
#endif

        ! Write message for new file creation
        call write_message("Creating source file " // trim(filename) // "...", &
             &1)

        ! Create separate source file
        call sp % file_create(filename, serial = .false.)

        ! Write file type
        call sp % write_data(FILETYPE_SOURCE, "filetype")

      else

        ! Set filename for state point
        filename = trim(path_output) // 'statepoint.' // &
            & zero_padded(current_batch, count_digits(n_batches))

        if (dd_run) then
          filename = trim(filename) // '.domain_' // &
              & zero_padded(domain_decomp % meshbin, &
                            count_digits(domain_decomp % n_domains))
        end if

#ifdef HDF5
        filename = trim(filename) // '.h5'
#else
        filename = trim(filename) // '.binary'
#endif

        ! Reopen statepoint file in parallel
        call sp % file_open(filename, 'w', serial = .false.)

      end if

      ! Write out source
      call write_source_bank(sp)

      ! Close file
      call sp % file_close()

    end if

    ! Also check to write source separately in overwritten file
    if (source_latest) then

      ! Set filename
      filename = trim(path_output) // 'source'
#ifdef HDF5
      filename = trim(filename) // '.h5'
#else
      filename = trim(filename) // '.binary'
#endif

      ! Write message for new file creation
      call write_message("Creating source file " // trim(filename) // "...", 1)

      ! Always create this file because it will be overwritten
      call sp % file_create(filename, serial = .false.)

      ! Write file type
      call sp % write_data(FILETYPE_SOURCE, "filetype")

      ! Write out source
      call write_source_bank(sp)

      ! Close file
      call sp % file_close()

    end if

  end subroutine write_source_point

!===============================================================================
! WRITE_TALLY_RESULTS_NR
!===============================================================================

  subroutine write_tally_results_nr()

    integer :: i      ! loop index
    integer :: n      ! number of filter bins
    integer :: m      ! number of score bins
    integer :: n_bins ! total number of bins
    real(8), allocatable :: tally_temp(:,:,:) ! contiguous array of results
    real(8), target :: global_temp(2,N_GLOBAL_TALLIES)
#ifdef MPI
    real(8) :: dummy  ! temporary receive buffer for non-root reduces
#endif
    type(TallyObject), pointer :: t => null()
    type(TallyResult), allocatable :: tallyresult_temp(:,:)

    ! ==========================================================================
    ! COLLECT AND WRITE GLOBAL TALLIES

    if (master) then
      ! Write number of realizations
      call sp % write_data(n_realizations, "n_realizations")

      ! Write number of global tallies
      call sp % write_data(N_GLOBAL_TALLIES, "n_global_tallies")
    end if

    ! Copy global tallies into temporary array for reducing
    n_bins = 2 * N_GLOBAL_TALLIES
    global_temp(1,:) = global_tallies(:) % sum
    global_temp(2,:) = global_tallies(:) % sum_sq

    if (master) then
      ! The MPI_IN_PLACE specifier allows the master to copy values into a
      ! receive buffer without having a temporary variable
#ifdef MPI
      call MPI_REDUCE(MPI_IN_PLACE, global_temp, n_bins, MPI_REAL8, MPI_SUM, &
           0, MPI_COMM_WORLD, mpi_err)
#endif

      ! Transfer values to value on master
      if (current_batch == n_batches) then
        global_tallies(:) % sum    = global_temp(1,:)
        global_tallies(:) % sum_sq = global_temp(2,:)
      end if

      ! Put reduced value in temporary tally result
      allocate(tallyresult_temp(N_GLOBAL_TALLIES, 1))
      tallyresult_temp(:,1) % sum    = global_temp(1,:)
      tallyresult_temp(:,1) % sum_sq = global_temp(2,:)


      ! Write out global tallies sum and sum_sq
      call write_tally_result(sp, tallyresult_temp, "global_tallies", &
           n1=N_GLOBAL_TALLIES, n2=1)

      ! Deallocate temporary tally result
      deallocate(tallyresult_temp)
    else
      ! Receive buffer not significant at other processors
#ifdef MPI
      call MPI_REDUCE(global_temp, dummy, n_bins, MPI_REAL8, MPI_SUM, &
           0, MPI_COMM_WORLD, mpi_err)
#endif
    end if

    if (tallies_on) then
      ! Indicate that tallies are on
      if (master) then
        call sp % write_data(1, "tallies_present", group="tallies")
      end if

      ! Write all tally results
      TALLY_RESULTS: do i = 1, n_tallies
        t => tallies(i)

        ! Determine size of tally results array
        m = size(t % results, 1)
        n = size(t % results, 2)
        n_bins = m*n*2

        ! Allocate array for storing sums and sums of squares, but
        ! contiguously in memory for each
        allocate(tally_temp(2,m,n))
        tally_temp(1,:,:) = t % results(:,:) % sum
        tally_temp(2,:,:) = t % results(:,:) % sum_sq

        if (master) then
          ! The MPI_IN_PLACE specifier allows the master to copy values into
          ! a receive buffer without having a temporary variable
#ifdef MPI
          call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
#endif
          ! At the end of the simulation, store the results back in the
          ! regular TallyResults array
          if (current_batch == n_batches) then
            t % results(:,:) % sum = tally_temp(1,:,:)
            t % results(:,:) % sum_sq = tally_temp(2,:,:)
          end if

         ! Put in temporary tally result
         allocate(tallyresult_temp(m,n))
         tallyresult_temp(:,:) % sum    = tally_temp(1,:,:)
         tallyresult_temp(:,:) % sum_sq = tally_temp(2,:,:)

         ! Write reduced tally results to file
          call write_tally_result(sp, t % results, "results", &
               group="tallies/tally" // to_str(i), n1=m, n2=n)

          ! Deallocate temporary tally result
          deallocate(tallyresult_temp)
        else
          ! Receive buffer not significant at other processors
#ifdef MPI
          call MPI_REDUCE(tally_temp, dummy, n_bins, MPI_REAL8, MPI_SUM, &
               0, MPI_COMM_WORLD, mpi_err)
#endif
        end if

        ! Deallocate temporary copy of tally results
        deallocate(tally_temp)
      end do TALLY_RESULTS
    else
      if (master) then
        ! Indicate that tallies are off
        call sp % write_data(0, "tallies_present", group="tallies")
      end if
    end if

  end subroutine write_tally_results_nr

!===============================================================================
! LOAD_STATE_POINT
!===============================================================================

  subroutine load_state_point()

    character(MAX_FILE_LEN) :: path_temp
    character(19)           :: current_time
    character(52)           :: label
    integer                 :: i
    integer                 :: j
    integer                 :: k
    integer                 :: n_univ
    integer                 :: n_cell
    integer                 :: n_lat
    integer                 :: n_surf
    integer                 :: n_mat
    integer                 :: dist_dens
    integer                 :: dist_comp
    integer                 :: length(4)
    integer                 :: int_array(3)
    integer, allocatable    :: temp_array(:)
    integer, allocatable    :: temp_array3D(:,:,:)
    integer, allocatable    :: temp_array4D(:,:,:,:)
    logical                 :: source_present
    real(8)                 :: real_array(3) 
    real(8), allocatable    :: temp_real_array(:)
    type(TallyObject), pointer :: t => null()

    ! Write message
    call write_message("Loading state point " // trim(path_state_point) &
         &// "...", 1)

    ! Open file for reading
    call sp % file_open(path_state_point, 'r', serial = .false.)

    ! Read filetype
    call sp % read_data(int_array(1), "filetype")

    ! Read revision number for state point file and make sure it matches with
    ! current version
    call sp % read_data(int_array(1), "revision")
    if (int_array(1) /= REVISION_STATEPOINT) then
      call fatal_error("State point version does not match current version &
           &in OpenMC.")
    end if

    ! Read OpenMC version
    call sp % read_data(int_array(1), "version_major")
    call sp % read_data(int_array(2), "version_minor")
    call sp % read_data(int_array(3), "version_release")
    if (int_array(1) /= VERSION_MAJOR .or. int_array(2) /= VERSION_MINOR &
        .or. int_array(3) /= VERSION_RELEASE) then
      if (master) call warning("State point file was created with a different &
           &version of OpenMC.")
    end if

    ! Read date and time
    call sp % read_data(current_time, "date_and_time")

    ! Read path to input
    call sp % read_data(path_temp, "path")

    ! Read and overwrite random number seed
    call sp % read_data(seed, "seed")

    ! Read and overwrite run information except number of batches
    call sp % read_data(run_mode, "run_mode")
    call sp % read_data(n_particles, "n_particles")

    ! Read batch number to restart at
    call sp % read_data(restart_batch, "current_batch")

    if (restart_batch > n_batches) then
      call fatal_error("The number batches specified in settings.xml is fewer &
           & than the number of batches in the given statepoint file.")
    end if

    ! Read information specific to eigenvalue run
    if (run_mode == MODE_EIGENVALUE) then
      call sp % read_data(int_array(1), "n_inactive")
      call sp % read_data(gen_per_batch, "gen_per_batch")
      call sp % read_data(k_generation, "k_generation", &
           length=restart_batch*gen_per_batch)
      call sp % read_data(entropy, "entropy", length=restart_batch*gen_per_batch)
      call sp % read_data(k_col_abs, "k_col_abs")
      call sp % read_data(k_col_tra, "k_col_tra")
      call sp % read_data(k_abs_tra, "k_abs_tra")
      call sp % read_data(real_array(1:2), "k_combined", length=2)

      ! Take maximum of statepoint n_inactive and input n_inactive
      n_inactive = max(n_inactive, int_array(1))

      ! Read in to see if CMFD was on
      call sp % read_data(int_array(1), "cmfd_on")

      ! Write out CMFD info
      if (int_array(1) == 1) then
        call sp % read_data(cmfd % indices, "indices", length=4, group="cmfd")
        call sp % read_data(cmfd % k_cmfd, "k_cmfd", length=restart_batch, &
             group="cmfd")
        length = cmfd % indices([4,1,2,3])
        call sp % read_data(cmfd % cmfd_src, "cmfd_src", &
             length=length, group="cmfd")
        call sp % read_data(cmfd % entropy, "cmfd_entropy", &
                       length=restart_batch, group="cmfd")
        call sp % read_data(cmfd % balance, "cmfd_balance", &
             length=restart_batch, group="cmfd")
        call sp % read_data(cmfd % dom, "cmfd_dominance", &
             length = restart_batch, group="cmfd")
        call sp % read_data(cmfd % src_cmp, "cmfd_srccmp", &
             length = restart_batch, group="cmfd")
      end if
    end if

#ifdef HDF5
        ! Do not bother reading geometry
#else
    ! Begin reading geometry information
    ! Set all values to temp variables. We already reread the geometry
    call sp % read_data(n_cell, "n_cells", group="geometry")
    call sp % read_data(n_univ, "n_universes", group="geometry")
    call sp % read_data(n_lat, "n_lattices", group="geometry")
    call sp % read_data(n_surf, "n_surfaces", group="geometry")
    call sp % read_data(n_mat, "n_materials", group="geometry")

    if (n_lat > 0) then
      ! Read list of lattice IDs
      allocate(temp_array(n_lat))
      call sp % read_data(temp_array, "lattice_ids", &
           group="geometry", length=n_lat)
      deallocate(temp_array)
    end if

    ! Read list of universe IDs
    allocate(temp_array(n_univ))
    call sp % read_data(temp_array, "universe_ids", &
         group="geometry", length=n_univ)
    deallocate(temp_array)

    ! Read list of surface IDs
    allocate(temp_array(n_surf))
    call sp % read_data(temp_array, "surface_keys", &
         group="geometry", length=n_surf)
    call sp % read_data(temp_array, "surface_ids", &
         group="geometry", length=n_surf)
    deallocate(temp_array)

    ! Read list of material IDs
    allocate(temp_array(n_mat))
    call sp % read_data(temp_array, "material_ids", &
         group="geometry", length=n_mat)
    deallocate(temp_array)

    ! Read list of cell keys-> IDs
    allocate(temp_array(n_cell))
    call sp % read_data(temp_array, "cell_keys", &
         group="geometry", length=n_cell)
    call sp % read_data(temp_array, "cell_ids", &
         group="geometry", length=n_cell)
    deallocate(temp_array)

    ! =========================================================================
    ! READ INFORMATION ON CELLS

    ! Read information on each cell
    CELL_LOOP: do i = 1, n_cell 
      ! Write information on what fills this cell
      call sp % read_data(j, "fill_type", &
             group="geometry/cells/cell ")

      call sp % read_data(k, "n_surfaces", &
             group="geometry/cells/cell ")
      if (k > 0) then
        allocate(temp_array(k))
        call sp % read_data(temp_array, "surfaces", &
               group="geometry/cells/cell ", length=k)
        deallocate(temp_array)
      endif

      select case (j)
      case (CELL_NORMAL)          
          call sp % read_data(k, "material", &
               group="geometry/cells/cell ")
      case (CELL_FILL)
        call sp % read_data(k, "fill", &
             group="geometry/cells/cell ") 
        call sp % read_data(k, "maps", &
             group="geometry/cells/cell ")
        if (k > 0) then
          allocate(temp_array(k))
          call sp % read_data(temp_array, "offset", length=k, &
               group="geometry/cells/cell ")
          deallocate(temp_array)
        end if
      case (CELL_LATTICE)
        call sp % read_data(j, "lattice", &
             group="geometry/cells/cell ") 
      end select

    end do CELL_LOOP

    ! =========================================================================
    ! Read INFORMATION ON UNIVERSES

    ! Read information on each universe
    UNIVERSE_LOOP: do i = 1, n_univ
      call sp % read_data(j, "n_cells", &
           group="geometry/universes/universe ")

      ! Read list of cells in this universe
      if (j > 0) then
        allocate(temp_array(j))
        call sp % read_data(temp_array, "cells", length=j, &
             group="geometry/universes/universe ")
        deallocate(temp_array)
      end if
  

    end do UNIVERSE_LOOP

    ! =========================================================================
    ! READ INFORMATION ON LATTICES

    ! Read information on each lattice
    LATTICE_LOOP: do i = 1, n_lat

      ! Read lattice type
      call sp % read_data(j, "type", &
           group="geometry/lattices/lattice ")
      
      ! Read lattice dimensions, number of offset maps, and offsets

      call sp % read_data(int_array, "dimension", &
           length=3, &
           group="geometry/lattices/lattice ")

      call sp % read_data(j, "maps", &
           group="geometry/lattices/lattice ")

      call sp % read_data(j, "offset_size", &
           group="geometry/lattices/lattice ")
      if (j > 0) then
          allocate(temp_array4D(j,int_array(1), int_array(2), int_array(3)))
          call sp % read_data(temp_array4D, "offset", &
               length=(/j,int_array(1),int_array(2),int_array(3)/), &
               group="geometry/lattices/lattice ")
          deallocate(temp_array4D)
      end if

      ! Determine dimensions of lattice

      ! Read lattice universes
      allocate(temp_array3D(int_array(1), int_array(2), int_array(3)))
      call sp % read_data(temp_array3D, "universes", &
           length=(/int_array(1),int_array(2),int_array(3)/), &
           group="geometry/lattices/lattice ")
      deallocate(temp_array3D)

    end do LATTICE_LOOP   
    
    ! =========================================================================
    ! READ INFORMATION ON SURFACES
    
    ! Read information on each surface
    SURFACE_LOOP: do i = 1, n_surf
      call sp % read_data(j, "type", &
           group="geometry/surfaces/surface ")
      call sp % read_data(j, "bc", &
           group="geometry/surfaces/surface ")

      call sp % read_data(j, "n_coeffs", &
           group="geometry/surfaces/surface ")
      if (j > 0) then
        allocate(temp_real_array(j))
        call sp % read_data(temp_real_array, "coeffs", length=j, &
             group="geometry/surfaces/surface ")
        deallocate(temp_real_array)
      end if

      call sp % read_data(j, "n_neighbor_pos", &
           group="geometry/surfaces/surface ")
      if (j > 0) then
        allocate(temp_array(j))
        call sp % read_data(temp_array, "neighbor_pos", length=j, &
             group="geometry/surfaces/surface ")
        deallocate(temp_array)
      end if

      call sp % read_data(j, "n_neighbor_neg", &
           group="geometry/surfaces/surface ")
      if (j > 0) then
        allocate(temp_array(j))
        call sp % read_data(temp_array, "neighbor_neg", length=j, &
             group="geometry/surfaces/surface ")
        deallocate(temp_array)
      end if

    end do SURFACE_LOOP

    ! =========================================================================
    ! READ INFORMATION ON MATERIALS

    ! Read information on each material
    MATERIAL_LOOP: do i = 1, n_mat
      call sp % read_data(k, "n_nuclides", &
           group="geometry/materials/material ")
      allocate(temp_array(k))
      call sp % read_data(temp_array, "nuclide", length=k, &
           group="geometry/materials/material ")
      deallocate(temp_array)

      call sp % read_data(dist_dens, "distrib_dens", &
           group="geometry/materials/material ")
      call sp % read_data(dist_comp, "distrib_comp", &
           group="geometry/materials/material ")

      call sp % read_data(j, "n_density", &
           group="geometry/materials/material ")
      allocate(temp_real_array(j))
      call sp % read_data(temp_real_array, "density", length=j, &
           group="geometry/materials/material ")
      deallocate(temp_real_array)

      call sp % read_data(k, "n_comp", &
             group="geometry/materials/material ")

      call sp % read_data(j, "otf_compositions", &
           group="geometry/materials/material " // &
                 trim(to_str(materials(i) % id)))
      
      if (j == 0) then
        allocate(temp_real_array(k))
        COMPOSITION_LOOP: do j = 1, k
          call sp % read_data(temp_real_array, "atom_density ", &
               length=size(temp_real_array), & 
               group="geometry/materials/material ")
        end do COMPOSITION_LOOP
        deallocate(temp_real_array)
      end if

      call sp % read_data(j, "n_sab", &
           group="geometry/materials/material ")
      if (j > 0) then
        allocate(temp_array(j))
        call sp % read_data(temp_array, "i_sab_nuclides", length=j, &
             group="geometry/materials/material ")
        call sp % read_data(temp_array, "i_sab_tables", length=j, &
             group="geometry/materials/material ")
        deallocate(temp_array)
      end if

    end do MATERIAL_LOOP
#endif

    ! Read number of meshes
    call sp % read_data(n_meshes, "n_meshes", group="tallies")

    ! Read and overwrite mesh information
    MESH_LOOP: do i = 1, n_meshes
      call sp % read_data(meshes(i) % id, "id", &
           group="tallies/mesh" // to_str(i))
      call sp % read_data(meshes(i) % type, "type", &
           group="tallies/mesh" // to_str(i))
      call sp % read_data(meshes(i) % n_dimension, "n_dimension", &
           group="tallies/mesh" // to_str(i))
      call sp % read_data(meshes(i) % dimension, "dimension", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension)
      call sp % read_data(meshes(i) % lower_left, "lower_left", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension)
      call sp % read_data(meshes(i) % upper_right, "upper_right", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension)
      call sp % read_data(meshes(i) % width, "width", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension)
    end do MESH_LOOP

    ! Read and overwrite number of tallies
    call sp % read_data(n_tallies, "n_tallies", group="tallies")

    ! Read in tally metadata
    TALLY_METADATA: do i = 1, n_tallies

      ! Get pointer to tally
      t => tallies(i)
       
      ! Read tally id
      call sp % read_data(t % id, "id", group="tallies/tally" // to_str(i))

      ! Read tally label length
      call sp % read_data(j, "labellen", group="tallies/tally" // to_str(i))

      ! Read tally label
      if (j > 0) then
        call sp % read_data(label, "label", group="tallies/tally" // to_str(i))
      end if
        
      ! Read estimator type                                                                      
      call sp % read_data(t % estimator, "estimator", & 
           group="tallies/tally" // to_str(i))
                                  
      ! Read number of realizations
      call sp % read_data(t % n_realizations, "n_realizations", &
           group="tallies/tally" // to_str(i))

      ! Read size of tally results
      call sp % read_data(int_array(1), "total_score_bins", &
           group="tallies/tally" // to_str(i))
       
      call sp % read_data(int_array(2), "total_filter_bins", &
           group="tallies/tally" // to_str(i))

      ! Check size of tally results array
      if (int_array(1) /= t % total_score_bins .and. &
          int_array(2) /= t % total_filter_bins) then
        call fatal_error("Input file tally structure is different from &
             &restart.")
      end if

      ! Read number of filters
      call sp % read_data(t % n_filters, "n_filters", &
           group="tallies/tally" // to_str(i))

      ! Read filter information
      FILTER_LOOP: do j = 1, t % n_filters

        ! Read type of filter
        call sp % read_data(t % filters(j) % type, "type", &
             group="tallies/tally" // trim(to_str(i)) & 
             // "/filter" // to_str(j))

          ! Write offset for this filter
          call sp % read_data(t % filters(j) % offset, "offset", &
               group="tallies/tally" // trim(to_str(i)) &
               // "/filter" // to_str(j))

        ! Read number of bins for this filter
        call sp % read_data(t % filters(j) % n_bins, "n_bins", &
             group="tallies/tally" // trim(to_str(i)) &
             // "/filter" // to_str(j))

        ! Read bins
        if (t % filters(j) % type == FILTER_ENERGYIN .or. &
            t % filters(j) % type == FILTER_ENERGYOUT) then
          call sp % read_data(t % filters(j) % real_bins, "bins", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" & 
               // to_str(j), length=size(t % filters(j) % real_bins))
        else
          call sp % read_data(t % filters(j) % int_bins, "bins", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" & 
               // to_str(j), length=size(t % filters(j) % int_bins))
        end if

      end do FILTER_LOOP

      ! Read number of nuclide bins
      call sp % read_data(t % n_nuclide_bins, "n_nuclide_bins", &
           group="tallies/tally" // to_str(i))

      ! Set up nuclide bin array and then write
      allocate(temp_array(t % n_nuclide_bins))
      call sp % read_data(temp_array, "nuclide_bins", &
           group="tallies/tally" // to_str(i), length=t % n_nuclide_bins)
      NUCLIDE_LOOP: do j = 1, t % n_nuclide_bins
        if (temp_array(j) > 0) then
          nuclides(t % nuclide_bins(j)) % zaid = temp_array(j)
        else
          t % nuclide_bins(j) = temp_array(j)
        end if
      end do NUCLIDE_LOOP
      deallocate(temp_array)

      ! Write number of score bins, score bins, and scatt order
      call sp % read_data(t % n_score_bins, "n_score_bins", &
           group="tallies/tally" // to_str(i))
      call sp % read_data(t % score_bins, "score_bins", &
           group="tallies/tally" // to_str(i), length=t % n_score_bins)
      call sp % read_data(t % moment_order, "moment_order", &
           group="tallies/tally" // to_str(i), length=t % n_score_bins)

      ! Write number of user score bins
      call sp % read_data(t % n_user_score_bins, "n_user_score_bins", &
           group="tallies/tally" // to_str(i))

    end do TALLY_METADATA

    ! Check for source in statepoint if needed
    call sp % read_data(int_array(1), "source_present")
    if (int_array(1) == 1) then
      source_present = .true.
    else
      source_present = .false.
    end if

    ! Check to make sure source bank is present
    if (path_source_point == path_state_point .and. .not. source_present) then
      call fatal_error("Source bank must be contained in statepoint restart &
           &file")
    end if

    ! Read tallies to master
    if (master) then

      ! Read number of realizations for global tallies
      call sp % read_data(n_realizations, "n_realizations", collect=.false.)

      ! Read number of global tallies
      call sp % read_data(int_array(1), "n_global_tallies", collect=.false.)
      if (int_array(1) /= N_GLOBAL_TALLIES) then
        call fatal_error("Number of global tallies does not match in state &
             &point.")
      end if

      ! Read global tally data
      call read_tally_result(sp, global_tallies, "global_tallies", &
           n1=N_GLOBAL_TALLIES, n2=1)

      ! Check if tally results are present
      call sp % read_data(int_array(1), "tallies_present", group="tallies", collect=.false.)

      ! Read in sum and sum squared
      if (int_array(1) == 1) then
        TALLY_RESULTS: do i = 1, n_tallies

          ! Set pointer to tally
          t => tallies(i)

          ! Read sum and sum_sq for each bin
          call read_tally_result(sp, t % results, "results", &
               group="tallies/tally" // to_str(i), &
               n1=size(t % results, 1), n2=size(t % results, 2))

        end do TALLY_RESULTS
      end if
    end if

    ! Read source if in eigenvalue mode
    if (run_mode == MODE_EIGENVALUE) then

      ! Check if source was written out separately
      if (.not. source_present) then

        ! Close statepoint file
        call sp % file_close()

        ! Write message
        call write_message("Loading source file " // trim(path_source_point) &
             &// "...", 1)

        ! Open source file
        call sp % file_open(path_source_point, 'r', serial = .false.)

        ! Read file type
        call sp % read_data(int_array(1), "filetype")

      end if

      ! Write out source
      call read_source_bank(sp)

    end if

    ! Close file
    call sp % file_close()

  end subroutine load_state_point

  subroutine read_source
! TODO write this routine
! TODO what if n_particles does not match source bank
  end subroutine read_source


!===============================================================================
! WRITE_DISTRIBMAT_COMPS
!===============================================================================

  subroutine write_distribmat_comps()

    character(MAX_FILE_LEN)    :: filename
    integer :: i, j
    integer :: idx
    type(BinaryOutput) :: fh
    type(Material), pointer :: mat => null()

#ifdef HDF5
    integer(HID_T) :: file_id
#endif

    ! Create files and write headers (master only)
    if (master) then
      do i = 1, n_materials
        mat => materials(i)
        if (mat % n_comp > 1) then

          filename = trim(path_output) // 'material'
          filename = trim(filename) // '.m' // &
              & zero_padded(mat % id, count_digits(n_materials))
#ifdef HDF5
          filename = trim(filename) // '.h5'
#else
          filename = trim(filename) // '.binary'
#endif
          ! Create file and write header
          call fh % file_create(filename, record_len = 8)
          call fh % write_data(mat % n_nuclides, 'n_nuclides', record=1)
          call fh % write_data(mat % n_comp, 'n_instances', record=2)
          call fh % file_close()
#ifdef HDF5
          ! Create the full dataset initially so all other procs can write to it
          hdf5_rank = 1
          dims1(1) = mat % n_comp * mat % n_nuclides
          call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, hdf5_err)
          call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)
          call h5dcreate_f(file_id, 'comps', H5T_NATIVE_DOUBLE, dspace, dset, hdf5_err)
          call h5dclose_f(dset, hdf5_err)
          call h5sclose_f(dspace, hdf5_err)
          call h5fclose_f(file_id, hdf5_err)
#endif
        end if
      end do
    end if

#ifdef MPI
    ! For parallel io we need to wait for master to create the files
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

    ! Write compositions
        ! Only master processes write the compositions
    if (master .or. (dd_run .and. domain_decomp % local_master)) then
      do i = 1, n_materials
        mat => materials(i)
        if (mat % n_comp > 1) then
          
          filename = trim(path_output) // 'material'
          filename = trim(filename) // '.m' // &
              & zero_padded(mat % id, count_digits(n_materials))
#ifdef HDF5
          filename = trim(filename) // '.h5'
#else
          filename = trim(filename) // '.binary'
#endif
          if (master) then
            call write_message("Writing distributed material " // trim(filename) &
                               // "...", 1)
          end if

          call fh % file_open(filename, 'w', serial = .false., &
                              direct_access = .true., &
                              record_len = 8 * mat % n_nuclides)
    
          do j = 1, mat % n_comp

            if (mat % otf_compositions) then
              if (mat % reverse_comp_index_map % has_key(j)) then
                idx = mat % reverse_comp_index_map % get_key(j)
              else
                cycle
              end if
            else
              idx = j
            end if

            call fh % write_data(mat % comp(j) % atom_density, "comps", &
                 length=mat % n_nuclides, record=idx, offset=16, collect=.false.)
          end do 

          call fh % file_close()

        end if
      end do
    end if

  end subroutine write_distribmat_comps

end module state_point
