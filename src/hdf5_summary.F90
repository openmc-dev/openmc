module hdf5_summary

#ifdef HDF5

  use ace_header,      only: Reaction, UrrData, Nuclide
  use constants
  use endf,            only: reaction_name
  use geometry_header, only: Cell, Surface, Universe, Lattice
  use global
  use hdf5_interface
  use material_header, only: Material
  use output_interface
  use output,          only: time_stamp
  use string,          only: to_str
  use tally_header,    only: TallyObject

contains

!===============================================================================
! HDF5_WRITE_SUMMARY
!===============================================================================

  subroutine hdf5_write_summary()

    character(MAX_FILE_LEN) :: filename = "summary"

    ! Create a new file using default properties.
    call file_create(filename, "serial") 

    ! Write header information
    call hdf5_write_header()

    ! Write eigenvalue information
    if (run_mode == MODE_EIGENVALUE) then

      ! Write number of particles
      call write_data(n_particles, "n_particles")

      ! Use H5LT interface to write n_batches, n_inactive, and n_active
      call write_data(n_batches, "n_batches")
      call write_data(n_inactive, "n_inactive")
      call write_data(n_active, "n_active")
      call write_data(gen_per_batch, "gen_per_batch")

      ! Add description of each variable
      call write_attribute_string("n_particles", &
           "description", "Number of particles per generation")
      call write_attribute_string("n_batches", &
           "description", "Total number of batches")
      call write_attribute_string("n_inactive", &
           "description", "Number of inactive batches")
      call write_attribute_string("n_active", &
           "description", "Number of active batches")
      call write_attribute_string("gen_per_batch", &
           "description", "Number of generations per batch")
    end if

    call hdf5_write_geometry()
    call hdf5_write_materials()
    call hdf5_write_nuclides()
    if (n_tallies > 0) then
      call hdf5_write_tallies()
    end if

    ! Terminate access to the file.
    call h5fclose_f(hdf5_output_file, hdf5_err)

  end subroutine hdf5_write_summary

!===============================================================================
! HDF5_WRITE_HEADER
!===============================================================================

  subroutine hdf5_write_header()

    ! Write version information
    call write_data(VERSION_MAJOR, "version_major")
    call write_data(VERSION_MINOR, "version_minor")
    call write_data(VERSION_RELEASE, "version_release") 

    ! Write current date and time
    call write_data(time_stamp(), "date_and_time")

    ! Write MPI information
    call write_data(n_procs, "n_procs")
    call write_attribute_string("n_procs", "description", &
         "Number of MPI processes") 

  end subroutine hdf5_write_header

!===============================================================================
! HDF5_WRITE_GEOMETRY
!===============================================================================

  subroutine hdf5_write_geometry()

    integer          :: i, j, k, m
    integer          :: n_x, n_y, n_z
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: dims3(3)
    integer(HID_T)   :: geometry_group
    integer(HID_T)   :: cell_group
    integer(HID_T)   :: surface_group
    integer(HID_T)   :: universe_group
    integer(HID_T)   :: lattice_group
    integer(HID_T)   :: temp_group
    integer, allocatable :: lattice_universes(:,:,:)
    type(Cell),     pointer :: c => null()
    type(Surface),  pointer :: s => null()
    type(Universe), pointer :: u => null()
    type(Lattice),  pointer :: lat => null()

    ! Use H5LT interface to write number of geometry objects
    call write_data(n_cells, "n_cells", group="geometry")
    call write_data(n_surfaces, "n_surfaces", group="geometry")
    call write_data(n_universes, "n_universes", group="geometry")
    call write_data(n_lattices, "n_lattices", group="geometry")

    ! ==========================================================================
    ! WRITE INFORMATION ON CELLS

    ! Create a cell group (nothing directly written in this group) then close
    call hdf5_open_group("geometry/cells")
    call hdf5_close_group()

    ! Write information on each cell
    CELL_LOOP: do i = 1, n_cells
      c => cells(i)

      ! Write universe for this cell
      call write_data(universes(c % universe) % id, "universe", &
           group="geometry/cells/cell " // trim(to_str(c % id)))

      ! Write information on what fills this cell
      select case (c % type)
      case (CELL_NORMAL)
        call write_data("normal", "fill_type", &
             group="geometry/cells/cell " // trim(to_str(c % id)))
        if (c % material == MATERIAL_VOID) then
          call write_data(-1, "material", &
               group="geometry/cells/cell " // trim(to_str(c % id)))
        else
          call write_data(materials(c % material) % id, "material", &
               group="geometry/cells/cell " // trim(to_str(c % id)))
        end if
      case (CELL_FILL)
        call write_data("universe", "fill_type", &
             group="geometry/cells/cell " // trim(to_str(c % id)))
        call write_data(universes(c % fill) % id, "material", &
             group="geometry/cells/cell " // trim(to_str(c % id))) 
      case (CELL_LATTICE)
        call write_data("lattice", "fill_type", &
             group="geometry/cells/cell " // trim(to_str(c % id)))
        call write_data(lattices(c % fill) % id, "lattice", &
             group="geometry/cells/cell " // trim(to_str(c % id))) 
      end select

      ! Write list of bounding surfaces
      if (c % n_surfaces > 0) then
        call write_data(c % surfaces, "surfaces", length= c % n_surfaces, &
             group="geometry/cells/cell " // trim(to_str(c % id)))
      end if

    end do CELL_LOOP

    ! ==========================================================================
    ! WRITE INFORMATION ON SURFACES

    ! Create surfaces group (nothing directly written here) then close
    call hdf5_open_group("geometry/surfaces")
    call hdf5_close_group()

    ! Write information on each surface
    SURFACE_LOOP: do i = 1, n_surfaces
      s => surfaces(i)

      ! Write surface type
      select case (s % type)
      case (SURF_PX)
        call write_data("X Plane", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_PY)
        call write_data("Y Plane", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_PZ)
        call write_data("Z Plane", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_PLANE)
        call write_data("Plane", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CYL_X)
        call write_data("X Cylinder", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CYL_Y)
        call write_data("Y Cylinder", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CYL_Z)
        call write_data("Z Cylinder", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_SPHERE)
        call write_data("Sphere", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CONE_X)
        call write_data("X Cone", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CONE_Y)
        call write_data("Y Cone", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CONE_Z)
        call write_data("Z Cone", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      end select

      ! Write coefficients for surface
      call write_data(s % coeffs, "coefficients", length=size(s % coeffs), &
           group="geometry/surfaces/surface " // trim(to_str(s % id)))

      ! Write positive neighbors
      if (allocated(s % neighbor_pos)) then
        call write_data(s % neighbor_pos, "neighbors_positive", &
             length=size(s % neighbor_pos), &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      end if

      ! Write negative neighbors
      if (allocated(s % neighbor_neg)) then
        call write_data(s % neighbor_neg, "neighbors_negative", &
             length=size(s % neighbor_neg), &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      end if

      ! Write boundary condition
      select case (s % bc)
      case (BC_TRANSMIT)
        call write_data("transmission", "boundary_condition", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (BC_VACUUM)
        call write_data("vacuum", "boundary_condition", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (BC_REFLECT)
        call write_data("reflective", "boundary_condition", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (BC_PERIODIC)
        call write_data("periodic", "boundary_condition", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      end select

    end do SURFACE_LOOP

    ! ==========================================================================
    ! WRITE INFORMATION ON UNIVERSES

    ! Create universes group (nothing directly written here) then close
    call hdf5_open_group("geometry/universes")
    call hdf5_close_group()

    ! Write information on each universe
    UNIVERSE_LOOP: do i = 1, n_universes
      u => universes(i)

      ! Write list of cells in this universe
      if (u % n_cells > 0) then
        call write_data(u % cells, "cells", length=u % n_cells, &
             group="geometry/universes/universe " // trim(to_str(u % id)))
      end if

    end do UNIVERSE_LOOP

    ! ==========================================================================
    ! WRITE INFORMATION ON LATTICES

    ! Create lattices group (nothing directly written here) then close
    call hdf5_open_group("geometry/lattices")
    call hdf5_close_group()

    ! Write information on each lattice
    LATTICE_LOOP: do i = 1, n_lattices
      lat => lattices(i)

      ! Write lattice type
      select case(lat % type)
      case (LATTICE_RECT)
        call write_data("rectangular", "type", &
             group="geometry/lattices/lattice " // trim(to_str(lat % id)))
      case (LATTICE_HEX)
        call write_data("hexagonal", "type", &
             group="geometry/lattices/lattice " // trim(to_str(lat % id)))
      end select

      ! Write lattice dimensions, lower left corner, and width of element
      call write_data(lat % dimension, "dimension", &
           length=lat % n_dimension, &
           group="geometry/lattices/lattice " // trim(to_str(lat % id)))
      call write_data(lat % lower_left, "lower_left", &
           length=lat % n_dimension, &
           group="geometry/lattices/lattice " // trim(to_str(lat % id)))
      call write_data(lat % width, "width", &
           length=lat % n_dimension, &
           group="geometry/lattices/lattice " // trim(to_str(lat % id)))

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
      call write_data(lattice_universes, "universes", &
           length=(/n_x, n_y, n_z/), &
           group="geometry/lattices/lattice " // trim(to_str(lat % id)))
      deallocate(lattice_universes)

    end do LATTICE_LOOP
stop
  end subroutine hdf5_write_geometry

!===============================================================================
! HDF5_WRITE_MATERIALS
!===============================================================================

  subroutine hdf5_write_materials()

    integer          :: i
    integer          :: j
    integer(HSIZE_T) :: dims(1)
    integer(HID_T)   :: materials_group
    integer(HID_T)   :: temp_group
    integer, allocatable :: zaids(:)
    type(Material), pointer :: m => null()

    ! Create group for materials
    call h5gcreate_f(hdf5_output_file, "/materials", materials_group, hdf5_err)

    ! Use H5LT interface to write number of materials
    call hdf5_write_integer(materials_group, "n_materials", n_materials)

    ! Write information on each material
    do i = 1, n_materials
      m => materials(i)

      ! Create group for i-th universe
      call h5gcreate_f(materials_group, "material " // trim(to_str(m % id)), &
           temp_group, hdf5_err)

      ! Write atom density with units
      call hdf5_write_double(temp_group, "atom_density", m % density)
      call h5ltset_attribute_string_f(temp_group, "atom_density", &
           "units", "atom/barn-cm", hdf5_err)

      ! Copy ZAID for each nuclide to temporary array
      allocate(zaids(m % n_nuclides))
      do j = 1, m % n_nuclides
        zaids(j) = nuclides(m % nuclide(j)) % zaid
      end do

      ! Write temporary array to 'nuclides'
      dims(1) = m % n_nuclides
      call h5ltmake_dataset_int_f(temp_group, "nuclides", 1, &
           dims, zaids, hdf5_err)

      ! Deallocate temporary array
      deallocate(zaids)

      ! Write atom densities
      call h5ltmake_dataset_double_f(temp_group, "nuclide_densities", 1, &
           dims, m % atom_density, hdf5_err)

      ! Write S(a,b) information if present
      if (m % n_sab > 0) then
        dims(1) = m % n_sab
        call h5ltmake_dataset_int_f(temp_group, "i_sab_nuclides", 1, &
             dims, m % i_sab_nuclides, hdf5_err)
        call h5ltmake_dataset_int_f(temp_group, "i_sab_tables", 1, &
             dims, m % i_sab_tables, hdf5_err)
      end if

      ! Close group for i-th material
      call h5gclose_f(temp_group, hdf5_err)
    end do

    ! Close materials group
    call h5gclose_f(materials_group, hdf5_err)

  end subroutine hdf5_write_materials

!===============================================================================
! HDF5_WRITE_TALLIES
!===============================================================================

  subroutine hdf5_write_tallies()

    integer           :: i, j
    integer(HSIZE_T)  :: dims(1)
    integer(HID_T)    :: tallies_group
    integer(HID_T)    :: temp_group
    integer(HID_T)    :: filter_group     ! group for i-th filter
    integer, allocatable :: temp_array(:) ! nuclide bin array
    type(TallyObject), pointer :: t => null()

    ! Create group for tallies
    call h5gcreate_f(hdf5_output_file, "tallies", tallies_group, hdf5_err)

    ! Write total number of meshes
    call hdf5_write_integer(tallies_group, "n_meshes", n_meshes)

    ! Write information for meshes
    MESH_LOOP: do i = 1, n_meshes
      ! Create temporary group for each mesh
      call h5gcreate_f(tallies_group, "mesh" // to_str(i), &
           temp_group, hdf5_err)

      ! Write type and number of dimensions
      call hdf5_write_integer(temp_group, "type", meshes(i) % type)
      call hdf5_write_integer(temp_group, "n_dimension", &
           meshes(i) % n_dimension)

      ! Write mesh information
      dims(1) = meshes(i) % n_dimension
      call h5ltmake_dataset_int_f(temp_group, "dimension", 1, &
           dims, meshes(i) % dimension, hdf5_err)
      call h5ltmake_dataset_double_f(temp_group, "lower_left", 1, &
           dims, meshes(i) % lower_left, hdf5_err)
      call h5ltmake_dataset_double_f(temp_group, "upper_right", 1, &
           dims, meshes(i) % upper_right, hdf5_err)
      call h5ltmake_dataset_double_f(temp_group, "width", 1, &
           dims, meshes(i) % width, hdf5_err)

      ! Close temporary group for mesh
      call h5gclose_f(temp_group, hdf5_err)
    end do MESH_LOOP

    ! Write number of tallies
    call hdf5_write_integer(tallies_group, "n_tallies", n_tallies)

    TALLY_METADATA: do i = 1, n_tallies
      ! Get pointer to tally
      t => tallies(i)

      ! Create group for this tally
      call h5gcreate_f(tallies_group, "tally" // to_str(i), &
           temp_group, hdf5_err)

      ! Write size of each tally
      call hdf5_write_integer(temp_group, "total_score_bins", &
           t % total_score_bins)
      call hdf5_write_integer(temp_group, "total_filter_bins", &
           t % total_filter_bins)

      ! Write number of filters
      call hdf5_write_integer(temp_group, "n_filters", t % n_filters)

      FILTER_LOOP: do j = 1, t % n_filters
        ! Create filter group
        call h5gcreate_f(temp_group, "filter" // to_str(j), filter_group, &
             hdf5_err)

        ! Write type of filter
        call hdf5_write_integer(filter_group, "type", t % filters(j) % type)

        ! Write number of bins for this filter
        call hdf5_write_integer(filter_group, "n_bins", t % filters(j) % n_bins)

        ! Write filter bins
        if (t % filters(j) % type == FILTER_ENERGYIN .or. &
             t % filters(j) % type == FILTER_ENERGYOUT) then
          dims(1) = size(t % filters(j) % real_bins)
          call h5ltmake_dataset_double_f(filter_group, "bins", 1, &
               dims, t % filters(j) % real_bins, hdf5_err)
        else
          dims(1) = size(t % filters(j) % int_bins)
          call h5ltmake_dataset_int_f(filter_group, "bins", 1, &
               dims, t % filters(j) % int_bins, hdf5_err)
        end if

        ! Write name of type
        select case (t % filters(j) % type)
        case(FILTER_UNIVERSE)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "universe", hdf5_err)
        case(FILTER_MATERIAL)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "material", hdf5_err)
        case(FILTER_CELL)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "cell", hdf5_err)
        case(FILTER_CELLBORN)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "cellborn", hdf5_err)
        case(FILTER_SURFACE)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "surface", hdf5_err)
        case(FILTER_MESH)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "mesh", hdf5_err)
        case(FILTER_ENERGYIN)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "energy", hdf5_err)
        case(FILTER_ENERGYOUT)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "energyout", hdf5_err)
        end select

        ! Close group for this filter
        call h5gclose_f(filter_group, hdf5_err)
      end do FILTER_LOOP

      ! Write number of nuclide bins
      call hdf5_write_integer(temp_group, "n_nuclide_bins", &
           t % n_nuclide_bins)


      ! Create temporary array for nuclide bins
      allocate(temp_array(t % n_nuclide_bins))
      NUCLIDE_LOOP: do j = 1, t % n_nuclide_bins
        if (t % nuclide_bins(j) > 0) then
          temp_array(j) = nuclides(t % nuclide_bins(j)) % zaid
        else
          temp_array(j) = t % nuclide_bins(j)
        end if
      end do NUCLIDE_LOOP

      ! Write and deallocate nuclide bins
      dims(1) = t % n_nuclide_bins
      call h5ltmake_dataset_int_f(temp_group, "nuclide_bins", 1, &
           dims, temp_array, hdf5_err)
      deallocate(temp_array)

      ! Write number of score bins
      call hdf5_write_integer(temp_group, "n_score_bins", &
           t % n_score_bins)
      dims(1) = t % n_score_bins
      call h5ltmake_dataset_int_f(temp_group, "score_bins", 1, &
           dims, t % score_bins, hdf5_err)

      ! Close tally group
      call h5gclose_f(temp_group, hdf5_err)
    end do TALLY_METADATA

    ! Close tallies group
    call h5gclose_f(tallies_group, hdf5_err)

  end subroutine hdf5_write_tallies

!===============================================================================
! HDF5_WRITE_NUCLIDES
!===============================================================================

  subroutine hdf5_write_nuclides()

    integer        :: i, j
    integer        :: size_total
    integer        :: size_xs
    integer        :: size_angle
    integer        :: size_energy
    integer(HID_T) :: group
    integer(HID_T) :: nuclide_group
    integer(HID_T) :: reactions_group
    integer(HID_T) :: rxn_group
    type(Nuclide),  pointer :: nuc => null()
    type(Reaction), pointer :: rxn => null()
    type(UrrData),  pointer :: urr => null()

    ! Create group for nuclides
    call h5gcreate_f(hdf5_output_file, "/nuclides", group, hdf5_err)

    ! Use H5LT interface to write number of nuclides
    call hdf5_write_integer(group, "n_nuclides", n_nuclides_total)

    ! Write information on each nuclide
    do i = 1, n_nuclides_total
      nuc => nuclides(i)

      ! Determine size of cross-sections
      size_xs = (5 + nuc % n_reaction) * nuc % n_grid * 8
      size_total = size_xs

      ! Create group for i-th nuclide
      call h5gcreate_f(group, trim(nuc % name), nuclide_group, hdf5_err)

      ! Write some basic attributes
      call hdf5_write_integer(nuclide_group, "zaid", nuc % zaid)
      call hdf5_write_double(nuclide_group, "awr", nuc % awr)
      call hdf5_write_double(nuclide_group, "kT", nuc % kT)
      call hdf5_write_integer(nuclide_group, "n_grid", nuc % n_grid)
      call hdf5_write_integer(nuclide_group, "n_reactions", nuc % n_reaction)
      call hdf5_write_integer(nuclide_group, "n_fission", nuc % n_fission)
      call hdf5_write_integer(nuclide_group, "size_xs", size_xs)

      ! =======================================================================
      ! WRITE INFORMATION ON EACH REACTION

      ! Create overall group 
      call h5gcreate_f(nuclide_group, "reactions", reactions_group, hdf5_err)

      do j = 1, nuc % n_reaction
        ! Information on each reaction
        rxn => nuc % reactions(j)

        ! Determine size of angle distribution
        if (rxn % has_angle_dist) then
          size_angle = rxn % adist % n_energy * 16 + size(rxn % adist % data) * 8
        else
          size_angle = 0
        end if

        ! Determine size of energy distribution
        if (rxn % has_energy_dist) then
          size_energy = size(rxn % edist % data) * 8
        else
          size_energy = 0
        end if

        ! Create reaction group 
        call h5gcreate_f(reactions_group, reaction_name(rxn % MT), &
             rxn_group, hdf5_err)

        ! Write information on reaction
        call hdf5_write_double(rxn_group, "Q_value", rxn % Q_value)
        call hdf5_write_integer(rxn_group, "multiplicity", rxn % multiplicity)
        call hdf5_write_integer(rxn_group, "threshold", rxn % threshold)
        call hdf5_write_integer(rxn_group, "size_angle", size_angle)
        call hdf5_write_integer(rxn_group, "size_energy", size_energy)

        call h5gclose_f(rxn_group, hdf5_err)

        ! Accumulate data size
        size_total = size_total + size_angle + size_energy
      end do

      ! Close overall group for reactions
      call h5gclose_f(reactions_group, hdf5_err)

      ! =======================================================================
      ! WRITE INFORMATION ON URR PROBABILITY TABLES

      if (nuc % urr_present) then
        urr => nuc % urr_data
        call hdf5_write_integer(nuclide_group, "urr_n_energy", urr % n_energy)
        call hdf5_write_integer(nuclide_group, "urr_n_prob", urr % n_prob)
        call hdf5_write_integer(nuclide_group, "urr_interp", urr % interp)
        call hdf5_write_integer(nuclide_group, "urr_inelastic", urr % inelastic_flag)
        call hdf5_write_integer(nuclide_group, "urr_absorption", urr % absorption_flag)
        call hdf5_write_double(nuclide_group, "urr_min_E", urr % energy(1))
        call hdf5_write_double(nuclide_group, "urr_max_E", urr % energy(urr % n_energy))
      end if

      ! Write total memory used
      call hdf5_write_integer(nuclide_group, "size_total", size_total)

      ! Close group for i-th nuclide
      call h5gclose_f(nuclide_group, hdf5_err)
    end do

    ! Close group for nuclides
    call h5gclose_f(group, hdf5_err)

  end subroutine hdf5_write_nuclides

!===============================================================================
! HDF5_WRITE_TIMING
!===============================================================================

  subroutine hdf5_write_timing()

    integer(HID_T)   :: timing_group
    integer(8)       :: total_particles
    real(8)          :: speed

    ! Create group for timing
    call h5gcreate_f(hdf5_output_file, "/timing", timing_group, hdf5_err)

    ! Write timing data
    call hdf5_write_double(timing_group, "time_initialize", time_initialize % elapsed)
    call hdf5_write_double(timing_group, "time_read_xs", time_read_xs % elapsed)
    call hdf5_write_double(timing_group, "time_unionize", time_unionize % elapsed)
    call hdf5_write_double(timing_group, "time_transport", time_transport % elapsed)
    call hdf5_write_double(timing_group, "time_bank", time_bank % elapsed)
    call hdf5_write_double(timing_group, "time_bank_sample", time_bank_sample % elapsed)
    call hdf5_write_double(timing_group, "time_bank_sendrecv", time_bank_sendrecv % elapsed)
    call hdf5_write_double(timing_group, "time_tallies", time_tallies % elapsed)
    call hdf5_write_double(timing_group, "time_inactive", time_inactive % elapsed)
    call hdf5_write_double(timing_group, "time_active", time_active % elapsed)
    call hdf5_write_double(timing_group, "time_finalize", time_finalize % elapsed)
    call hdf5_write_double(timing_group, "time_total", time_total % elapsed)

    ! Add descriptions to timing data
    call h5ltset_attribute_string_f(timing_group, "time_initialize", &
         "description", "Total time elapsed for initialization (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_read_xs", &
         "description", "Time reading cross-section libraries (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_unionize", &
         "description", "Time unionizing energy grid (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_transport", &
         "description", "Time in transport only (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_bank", &
         "description", "Total time synchronizing fission bank (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_bank_sample", &
         "description", "Time between generations sampling source sites (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_bank_sendrecv", &
         "description", "Time between generations SEND/RECVing source sites (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_tallies", &
         "description", "Time between batches accumulating tallies (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_inactive", &
         "description", "Total time in inactive batches (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_active", &
         "description", "Total time in active batches (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_finalize", &
         "description", "Total time for finalization (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_total", &
         "description", "Total time elapsed (s)", hdf5_err)

    ! Write calculation rate
    total_particles = n_particles * n_batches * gen_per_batch
    speed = real(total_particles) / (time_inactive % elapsed + &
         time_active % elapsed)
    call hdf5_write_double(timing_group, "neutrons_per_second", speed)

    ! Close timing group
    call h5gclose_f(timing_group, hdf5_err)

  end subroutine hdf5_write_timing

#endif

end module hdf5_summary
