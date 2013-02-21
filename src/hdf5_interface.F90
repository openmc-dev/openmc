module hdf5_interface

  use ace_header,      only: Reaction, UrrData
  use bank_header,     only: Bank
  use constants
  use endf,            only: reaction_name
  use error,           only: fatal_error
  use geometry_header, only: Cell, Surface, Universe, Lattice
  use global
  use material_header, only: Material
  use output,          only: write_message, time_stamp
  use string,          only: to_str
  use tally_header,    only: TallyObject

#ifdef HDF5
  use hdf5
  use h5lt
  use, intrinsic :: ISO_C_BINDING
#endif

  implicit none

#ifdef HDF5

contains

!===============================================================================
! HDF5_INITIALIZE
!===============================================================================

  subroutine hdf5_initialize()

    type(TallyResult), target :: tmp(2)          ! temporary TallyResult
    type(Bank),        target :: tmpb(2)         ! temporary Bank
    integer(HID_T)            :: coordinates_t   ! HDF5 type for 3 reals
    integer(HSIZE_T)          :: dims(1) = (/3/) ! size of coordinates

    ! Initialize FORTRAN interface.
    call h5open_f(hdf5_err)

    ! Create the compound datatype for TallyResult
    call h5tcreate_f(H5T_COMPOUND_F, h5offsetof(c_loc(tmp(1)), &
         c_loc(tmp(2))), hdf5_tallyresult_t, hdf5_err)
    call h5tinsert_f(hdf5_tallyresult_t, "sum", h5offsetof(c_loc(tmp(1)), &
         c_loc(tmp(1)%sum)), H5T_NATIVE_DOUBLE, hdf5_err)
    call h5tinsert_f(hdf5_tallyresult_t, "sum_sq", h5offsetof(c_loc(tmp(1)), &
         c_loc(tmp(1)%sum_sq)), H5T_NATIVE_DOUBLE, hdf5_err)

    ! Create compound type for xyz and uvw
    call h5tarray_create_f(H5T_NATIVE_DOUBLE, 1, dims, coordinates_t, hdf5_err)

    ! Create the compound datatype for Bank
    call h5tcreate_f(H5T_COMPOUND_F, h5offsetof(c_loc(tmpb(1)), &
         c_loc(tmpb(2))), hdf5_bank_t, hdf5_err)
    call h5tinsert_f(hdf5_bank_t, "wgt", h5offsetof(c_loc(tmpb(1)), &
         c_loc(tmpb(1)%wgt)), H5T_NATIVE_DOUBLE, hdf5_err)
    call h5tinsert_f(hdf5_bank_t, "xyz", h5offsetof(c_loc(tmpb(1)), &
         c_loc(tmpb(1)%xyz)), coordinates_t, hdf5_err)
    call h5tinsert_f(hdf5_bank_t, "uvw", h5offsetof(c_loc(tmpb(1)), &
         c_loc(tmpb(1)%uvw)), coordinates_t, hdf5_err)
    call h5tinsert_f(hdf5_bank_t, "E", h5offsetof(c_loc(tmpb(1)), &
         c_loc(tmpb(1)%E)), H5T_NATIVE_DOUBLE, hdf5_err)

    ! Determine type for integer(8)
    hdf5_integer8_t = h5kind_to_type(8, H5_INTEGER_KIND)

  end subroutine hdf5_initialize

!===============================================================================
! HDF5_FINALIZE
!===============================================================================

  subroutine hdf5_finalize()

    ! Release compound datatypes
    call h5tclose_f(hdf5_tallyresult_t, hdf5_err)

    ! Close FORTRAN interface.
    call h5close_f(hdf5_err)

  end subroutine hdf5_finalize

!===============================================================================
! HDF5_WRITE_SUMMARY
!===============================================================================

  subroutine hdf5_write_summary()

    character(MAX_FILE_LEN) :: filename = "summary.h5"

    ! Create a new file using default properties.
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5_output_file, hdf5_err)

    ! Write header information
    call hdf5_write_header()

    ! Write eigenvalue information
    if (run_mode == MODE_EIGENVALUE) then
      ! Need to write integer(8)'s using double instead since there is no H5LT
      ! call for making a dataset of type long
      call hdf5_write_double(hdf5_output_file, "n_particles", real(n_particles,8))

      ! Use H5LT interface to write n_batches, n_inactive, and n_active
      call hdf5_write_integer(hdf5_output_file, "n_batches", n_batches)
      call hdf5_write_integer(hdf5_output_file, "n_inactive", n_inactive)
      call hdf5_write_integer(hdf5_output_file, "n_active", n_active)
      call hdf5_write_integer(hdf5_output_file, "gen_per_batch", gen_per_batch)

      ! Add description of each variable
      call h5ltset_attribute_string_f(hdf5_output_file, "n_particles", &
           "description", "Number of particles per generation", hdf5_err)
      call h5ltset_attribute_string_f(hdf5_output_file, "n_batches", &
           "description", "Total number of batches", hdf5_err)
      call h5ltset_attribute_string_f(hdf5_output_file, "n_inactive", &
           "description", "Number of inactive batches", hdf5_err)
      call h5ltset_attribute_string_f(hdf5_output_file, "n_active", &
           "description", "Number of active batches", hdf5_err)
      call h5ltset_attribute_string_f(hdf5_output_file, "gen_per_batch", &
           "description", "Number of generations per batch", hdf5_err)
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
    call hdf5_write_integer(hdf5_output_file, "version_major", VERSION_MAJOR)
    call hdf5_write_integer(hdf5_output_file, "version_minor", VERSION_MINOR)
    call hdf5_write_integer(hdf5_output_file, "version_release", VERSION_RELEASE)

    ! Write current date and time
    call h5ltmake_dataset_string_f(hdf5_output_file, "date_and_time", &
         time_stamp(), hdf5_err)

    ! Write MPI information
    call hdf5_write_integer(hdf5_output_file, "n_procs", n_procs)
    call h5ltset_attribute_string_f(hdf5_output_file, "n_procs", &
         "description", "Number of MPI processes", hdf5_err)

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

    ! Create group for geometry
    call h5gcreate_f(hdf5_output_file, "/geometry", geometry_group, hdf5_err)

    ! Use H5LT interface to write number of geometry objects
    call hdf5_write_integer(geometry_group, "n_cells", n_cells)
    call hdf5_write_integer(geometry_group, "n_surfaces", n_surfaces)
    call hdf5_write_integer(geometry_group, "n_universes", n_universes)
    call hdf5_write_integer(geometry_group, "n_lattices", n_lattices)

    ! ==========================================================================
    ! WRITE INFORMATION ON CELLS

    call h5gcreate_f(geometry_group, "cells", cell_group, hdf5_err)

    ! Write information on each cell
    do i = 1, n_cells
      c => cells(i)

      ! Create group for i-th cell
      call h5gcreate_f(cell_group, "cell " // trim(to_str(c % id)), &
           temp_group, hdf5_err)

      ! Write universe for this cell
      call hdf5_write_integer(temp_group, "universe", &
           universes(c % universe) % id)

      ! Write information on what fills this cell
      select case (c % type)
      case (CELL_NORMAL)
        call h5ltmake_dataset_string_f(temp_group, "fill_type", "normal", &
             hdf5_err)
        if (c % material == MATERIAL_VOID) then
          call hdf5_write_integer(temp_group, "material", -1)
        else
          call hdf5_write_integer(temp_group, "material", &
               materials(c % material) % id)
        end if
      case (CELL_FILL)
        call h5ltmake_dataset_string_f(temp_group, "fill_type", "universe", &
             hdf5_err) 
        call hdf5_write_integer(temp_group, "material", &
             universes(c % fill) % id)
      case (CELL_LATTICE)
        call h5ltmake_dataset_string_f(temp_group, "fill_type", "lattice", &
             hdf5_err) 
        call hdf5_write_integer(temp_group, "lattice", &
             lattices(c % fill) % id)
      end select

      ! Write list of bounding surfaces
      if (c % n_surfaces > 0) then
        dims(1) = c % n_surfaces
        call h5ltmake_dataset_int_f(temp_group, "surfaces", 1, &
             dims, c % surfaces, hdf5_err)
      end if

      ! Close group for i-th cell
      call h5gclose_f(temp_group, hdf5_err)
    end do

    call h5gclose_f(cell_group, hdf5_err)

    ! ==========================================================================
    ! WRITE INFORMATION ON SURFACES

    call h5gcreate_f(geometry_group, "surfaces", surface_group, hdf5_err)

    ! Write information on each surface
    do i = 1, n_surfaces
      s => surfaces(i)

      ! Create group for i-th surface
      call h5gcreate_f(surface_group, "surface " // trim(to_str(s % id)), &
           temp_group, hdf5_err)

      ! Write surface type
      select case (s % type)
      case (SURF_PX)
        call h5ltmake_dataset_string_f(temp_group, "type", "X Plane", hdf5_err)
      case (SURF_PY)
        call h5ltmake_dataset_string_f(temp_group, "type", "Y Plane", hdf5_err)
      case (SURF_PZ)
        call h5ltmake_dataset_string_f(temp_group, "type", "Z Plane", hdf5_err)
      case (SURF_PLANE)
        call h5ltmake_dataset_string_f(temp_group, "type", "Plane", hdf5_err)
      case (SURF_CYL_X)
        call h5ltmake_dataset_string_f(temp_group, "type", "X Cylinder", hdf5_err)
      case (SURF_CYL_Y)
        call h5ltmake_dataset_string_f(temp_group, "type", "Y Cylinder", hdf5_err)
      case (SURF_CYL_Z)
        call h5ltmake_dataset_string_f(temp_group, "type", "Z Cylinder", hdf5_err)
      case (SURF_SPHERE)
        call h5ltmake_dataset_string_f(temp_group, "type", "Sphere", hdf5_err)
      case (SURF_CONE_X)
        call h5ltmake_dataset_string_f(temp_group, "type", "X Cone", hdf5_err)
      case (SURF_CONE_Y)
        call h5ltmake_dataset_string_f(temp_group, "type", "Y Cone", hdf5_err)
      case (SURF_CONE_Z)
        call h5ltmake_dataset_string_f(temp_group, "type", "Z Cone", hdf5_err)
      end select

      ! Write coefficients for surface
      dims(1) = size(s % coeffs)
      call h5ltmake_dataset_double_f(temp_group, "coefficients", 1, dims, &
           s % coeffs, hdf5_err)

      ! Write positive neighbors
      if (allocated(s % neighbor_pos)) then
        dims(1) = size(s % neighbor_pos)
        call h5ltmake_dataset_int_f(temp_group, "neighbors_positive", 1, dims, &
             s % neighbor_pos, hdf5_err)
      end if

      ! Write negative neighbors
      if (allocated(s % neighbor_neg)) then
        dims(1) = size(s % neighbor_neg)
        call h5ltmake_dataset_int_f(temp_group, "neighbors_negative", 1, dims, &
             s % neighbor_neg, hdf5_err)
      end if

      ! Write boundary condition
      select case (s % bc)
      case (BC_TRANSMIT)
        call h5ltmake_dataset_string_f(temp_group, "boundary_condition", &
             "transmission", hdf5_err)
      case (BC_VACUUM)
        call h5ltmake_dataset_string_f(temp_group, "boundary_condition", &
             "vacuum", hdf5_err)
      case (BC_REFLECT)
        call h5ltmake_dataset_string_f(temp_group, "boundary_condition", &
             "reflective", hdf5_err)
      case (BC_PERIODIC)
        call h5ltmake_dataset_string_f(temp_group, "boundary_condition", &
             "periodic", hdf5_err)
      end select

      ! Close group for i-th surface
      call h5gclose_f(temp_group, hdf5_err)
    end do

    call h5gclose_f(surface_group, hdf5_err)

    ! ==========================================================================
    ! WRITE INFORMATION ON UNIVERSES

    call h5gcreate_f(geometry_group, "universes", universe_group, hdf5_err)

    ! Write information on each universe
    do i = 1, n_universes
      u => universes(i)

      ! Create group for i-th universe
      call h5gcreate_f(universe_group, "universe " // trim(to_str(u % id)), &
           temp_group, hdf5_err)

      ! Write list of cells in this universe
      if (u % n_cells > 0) then
        dims(1) = u % n_cells
        call h5ltmake_dataset_int_f(temp_group, "cells", 1, dims, &
             u % cells, hdf5_err)
      end if

      ! Close group for i-th universe
      call h5gclose_f(temp_group, hdf5_err)
    end do

    call h5gclose_f(universe_group, hdf5_err)

    ! ==========================================================================
    ! WRITE INFORMATION ON LATTICES

    call h5gcreate_f(geometry_group, "lattices", lattice_group, hdf5_err)

    ! Write information on each lattice
    do i = 1, n_lattices
      lat => lattices(i)

      ! Create group for i-th lattice
      call h5gcreate_f(lattice_group, "lattice " // trim(to_str(lat % id)), &
           temp_group, hdf5_err)

      ! Write lattice type
      select case(lat % type)
      case (LATTICE_RECT)
        call h5ltmake_dataset_string_f(temp_group, "type", "rectangular", hdf5_err)
      case (LATTICE_HEX)
        call h5ltmake_dataset_string_f(temp_group, "type", "hexagonal", hdf5_err)
      end select

      ! Write lattice dimensions, lower left corner, and width of element
      dims(1) = lat % n_dimension
      call h5ltmake_dataset_int_f(temp_group, "dimension", 1, dims, &
           lat % dimension, hdf5_err)
      call h5ltmake_dataset_double_f(temp_group, "lower_left", 1, dims, &
           lat % lower_left, hdf5_err)
      call h5ltmake_dataset_double_f(temp_group, "width", 1, dims, &
           lat % width, hdf5_err)

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
      dims3 = (/ n_x, n_y, n_z /)
      call h5ltmake_dataset_int_f(temp_group, "universes", 3, dims3, &
           lattice_universes, hdf5_err)
      deallocate(lattice_universes)

      ! Close group for i-th lattice
      call h5gclose_f(temp_group, hdf5_err)
    end do

    call h5gclose_f(lattice_group, hdf5_err)

    ! Close geometry group
    call h5gclose_f(geometry_group, hdf5_err)

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
      call h5gcreate_f(tallies_group, "mesh " // to_str(meshes(i) % id), &
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
      call h5gcreate_f(tallies_group, "tally " // to_str(t % id), &
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
        call h5gcreate_f(temp_group, "filter " // to_str(j), filter_group, &
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

!===============================================================================
! HDF5_WRITE_STATE_POINT
!===============================================================================

  subroutine hdf5_write_state_point()

    integer              :: i, j             ! loop indices
    integer(HSIZE_T)     :: dims(1)          ! dimensions of 1D arrays
    integer(HSIZE_T)     :: dims2(2)         ! dimensions of 2D arrays
    integer(HSIZE_T)     :: dims4(4)         ! dimensions of 4D arrays
    integer(HID_T)       :: hdf5_state_point ! identifier for state point file
    integer(HID_T)       :: tallies_group    ! "tallies" group
    integer(HID_T)       :: temp_group       ! group for i-th tally or mesh
    integer(HID_T)       :: filter_group     ! group for i-th filter
    integer(HID_T)       :: cmfd_group       ! group for cmfd output
    integer(HID_T)       :: dspace           ! identifier for dataspace
    integer(HID_T)       :: dset             ! identifier for dataset
    integer, allocatable :: temp_array(:)    ! nuclide bin array
    type(c_ptr)          :: f_ptr            ! Pointer for h5dwrite
    type(TallyObject), pointer :: t => null()

    path_state_point = 'statepoint.' // trim(to_str(current_batch)) // '.h5'

    ! Write message
    message = "Creating HDF5 state point " // trim(path_state_point) // "..."
    call write_message(1)

    ! Only master process should continue from here
    if (.not. master) return

    ! Create a new state point file using default properties.
    call h5fcreate_f(path_state_point, H5F_ACC_TRUNC_F, hdf5_state_point, &
         hdf5_err)

    ! Write revision number for state point file
    call hdf5_write_integer(hdf5_state_point, "revision_statepoint", &
         REVISION_STATEPOINT)

    ! Write OpenMC version
    call hdf5_write_integer(hdf5_state_point, "version_major", VERSION_MAJOR)
    call hdf5_write_integer(hdf5_state_point, "version_minor", VERSION_MINOR)
    call hdf5_write_integer(hdf5_state_point, "version_release", VERSION_RELEASE)

    ! Write current date and time
    call h5ltmake_dataset_string_f(hdf5_state_point, "date_and_time", &
         time_stamp(), hdf5_err)

    ! Write path to input
    call h5ltmake_dataset_string_f(hdf5_state_point, "path", &
         path_input, hdf5_err)

    ! Write out random number seed
    call hdf5_write_long(hdf5_state_point, "seed", seed)

    ! Write run information
    call hdf5_write_integer(hdf5_state_point, "run_mode", run_mode)
    call hdf5_write_long(hdf5_state_point, "n_particles", n_particles)
    call hdf5_write_integer(hdf5_state_point, "n_batches", n_batches)

    ! Write out current batch number
    call hdf5_write_integer(hdf5_state_point, "current_batch", current_batch)

    ! Write out information for eigenvalue run
    if (run_mode == MODE_EIGENVALUE) then
      call hdf5_write_integer(hdf5_state_point, "n_inactive", n_inactive)
      call hdf5_write_integer(hdf5_state_point, "gen_per_batch", gen_per_batch)

      ! Write out keff and entropy
      dims(1) = current_batch
      call h5ltmake_dataset_double_f(hdf5_state_point, "k_batch", 1, &
           dims, k_batch, hdf5_err)
      dims(1) = current_batch*gen_per_batch
      call h5ltmake_dataset_double_f(hdf5_state_point, "entropy", 1, &
           dims, entropy, hdf5_err)
      call hdf5_write_double(hdf5_state_point, "k_col_abs", k_col_abs)
      call hdf5_write_double(hdf5_state_point, "k_col_tra", k_col_tra)
      call hdf5_write_double(hdf5_state_point, "k_abs_tra", k_abs_tra)
      dims(1) = 2
      call h5ltmake_dataset_double_f(hdf5_state_point, "k_combined", 1, &
           dims, k_combined, hdf5_err)
    end if

    ! Create group for tallies
    call h5gcreate_f(hdf5_state_point, "tallies", tallies_group, hdf5_err)

    ! Write total number of meshes
    call hdf5_write_integer(tallies_group, "n_meshes", n_meshes)

    ! Write information for meshes
    MESH_LOOP: do i = 1, n_meshes
      ! Create temporary group for each mesh
      call h5gcreate_f(tallies_group, "mesh " // to_str(meshes(i) % id), &
           temp_group, hdf5_err)

      ! Write id, type, and number of dimensions
      call hdf5_write_integer(temp_group, "id", meshes(i) % id)
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
      call h5gcreate_f(tallies_group, "tally " // to_str(t % id), &
           temp_group, hdf5_err)

      ! Write id
      call hdf5_write_integer(temp_group, "id", t % id)

      ! Write number of realizations
      call hdf5_write_integer(temp_group, "n_realizations", &
           t % n_realizations)

      ! Write size of each tally
      call hdf5_write_integer(temp_group, "total_score_bins", &
           t % total_score_bins)
      call hdf5_write_integer(temp_group, "total_filter_bins", &
           t % total_filter_bins)

      ! Write number of filters
      call hdf5_write_integer(temp_group, "n_filters", t % n_filters)

      FILTER_LOOP: do j = 1, t % n_filters
        ! Create filter group
        call h5gcreate_f(temp_group, "filter " // to_str(j), filter_group, &
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

      ! Write score bins and scattering order
      dims(1) = t % n_score_bins
      call h5ltmake_dataset_int_f(temp_group, "score_bins", 1, &
           dims, t % score_bins, hdf5_err)
      call h5ltmake_dataset_int_f(temp_group, "scatt_order", 1, &
           dims, t % scatt_order, hdf5_err)

      ! Write number of user score bins
      call hdf5_write_integer(temp_group, "n_user_score_bins", &
           t % n_user_score_bins)

      ! Close tally group
      call h5gclose_f(temp_group, hdf5_err)
    end do TALLY_METADATA

    ! Write number of realizations for global tallies
    call hdf5_write_integer(hdf5_state_point, "n_realizations", n_realizations)

    ! Write out global tallies sum and sum_sq
    call hdf5_write_integer(hdf5_state_point, "n_global_tallies", &
         N_GLOBAL_TALLIES)

    ! Write global tallies
    dims(1) = N_GLOBAL_TALLIES
    call h5screate_simple_f(1, dims, dspace, hdf5_err)
    call h5dcreate_f(hdf5_state_point, "global_tallies", hdf5_tallyresult_t, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(global_tallies(1))
    CALL h5dwrite_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

    if (tallies_on) then
      ! Indicate that tallies are on
      call hdf5_write_integer(tallies_group, "tallies_present", 1)

      ! Write tally sum and sum_sq
      TALLY_RESULTS: do i = 1, n_tallies
        ! Get pointer to tally
        t => tallies(i)

        ! Open group for the i-th tally
        call h5gopen_f(tallies_group, "tally " // to_str(t % id), &
             temp_group, hdf5_err)

        ! Write sum and sum_sq for each bin
        dims2 = shape(t % results)
        call h5screate_simple_f(2, dims2, dspace, hdf5_err)
        call h5dcreate_f(temp_group, "results", hdf5_tallyresult_t, &
             dspace, dset, hdf5_err)
        f_ptr = c_loc(t % results(1, 1))
        CALL h5dwrite_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)
        call h5dclose_f(dset, hdf5_err)
        call h5sclose_f(dspace, hdf5_err)

        ! Close group for the i-th tally
        call h5gclose_f(temp_group, hdf5_err)
      end do TALLY_RESULTS
    else
      ! Indicate that tallies are off
      call hdf5_write_integer(tallies_group, "tallies_present", 0)
    end if

    ! Close tallies group
    call h5gclose_f(tallies_group, hdf5_err)

    ! TODO: Use parallel HDF5 to write source bank

    ! Write source bank
    dims(1) = work
    call h5screate_simple_f(1, dims, dspace, hdf5_err)
    call h5dcreate_f(hdf5_state_point, "source_bank", hdf5_bank_t, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(source_bank(1))
    CALL h5dwrite_f(dset, hdf5_bank_t, f_ptr, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

    ! Write out CMFD info if active
    if (cmfd_on) then

      ! Create CMFD group
      call h5gcreate_f(hdf5_state_point, "cmfd", cmfd_group, hdf5_err)

      ! write out openmc source
      dims4 = shape(cmfd % openmc_src) 
      call h5ltmake_dataset_double_f(cmfd_group, "openmc_src", 4, &
           dims4, cmfd % openmc_src, hdf5_err)

      ! write out cmfd source
      dims4 = shape(cmfd % cmfd_src) 
      call h5ltmake_dataset_double_f(cmfd_group, "cmfd_src", 4, &
           dims4, cmfd % cmfd_src, hdf5_err)

      ! write out keff
      call hdf5_write_double(cmfd_group, "cmfd_keff", cmfd % keff) 

      ! Close CMFD group
      call h5gclose_f(tallies_group, hdf5_err)

    end if

    ! Close HDF5 state point file
    call h5fclose_f(hdf5_state_point, hdf5_err)

  end subroutine hdf5_write_state_point

!===============================================================================
! HDF5_LOAD_STATE_POINT
!===============================================================================

  subroutine hdf5_load_state_point()

    integer          :: i                ! loop index
    integer          :: mode             ! specified run mode
    integer          :: temp             ! statepoint revision for comparison
    integer(HID_T)   :: hdf5_state_point ! identifier for state point file
    integer(HID_T)   :: tally_group      ! identifier for tally group
    integer(HID_T)   :: dset             ! identifier for dataset
    integer(HSIZE_T) :: dims(1)          ! dimensions of 1D arrays
    type(c_ptr)      :: f_ptr            ! c pointer for h5dread

    ! Write message
    message = "Loading HDF5 state point " // trim(path_state_point) // "..."
    call write_message(1)

    ! Load state point file
    call h5fopen_f(path_state_point, H5F_ACC_RDONLY_F, &
         hdf5_state_point, hdf5_err)

    ! Raad revision number for state point file and make sure it matches with
    ! current version
    call hdf5_read_integer(hdf5_state_point, "revision_statepoint", temp)
    if (temp /= REVISION_STATEPOINT) then
      message = "State point binary version does not match current version " &
           // "in OpenMC."
      call fatal_error()
    end if

    ! Read and overwrite random number seed
    call hdf5_read_long(hdf5_state_point, "seed", seed)

    ! Read and overwrite run information
    call hdf5_read_integer(hdf5_state_point, "run_mode", mode)
    call hdf5_read_long(hdf5_state_point, "n_particles", n_particles)
    call hdf5_read_integer(hdf5_state_point, "n_batches", n_batches)

    ! Read batch number to restart at
    call hdf5_read_integer(hdf5_state_point, "current_batch", restart_batch)

    ! Read information specific to eigenvalue run
    if (mode == MODE_EIGENVALUE) then
      call hdf5_read_integer(hdf5_state_point, "n_inactive", n_inactive)
      call hdf5_read_integer(hdf5_state_point, "gen_per_batch", gen_per_batch)

      dims(1) = restart_batch
      call h5ltread_dataset_double_f(hdf5_state_point, "k_batch", &
           k_batch(1:restart_batch), dims, hdf5_err)
      call h5ltread_dataset_double_f(hdf5_state_point, "entropy", &
           entropy(1:restart_batch*gen_per_batch), dims, hdf5_err)
      call hdf5_read_double(hdf5_state_point, "k_col_abs", k_col_abs)
      call hdf5_read_double(hdf5_state_point, "k_col_tra", k_col_tra)
      call hdf5_read_double(hdf5_state_point, "k_abs_tra", k_abs_tra)
    end if

    ! Read number of realizations for global tallies
    call hdf5_read_integer(hdf5_state_point, "n_realizations", n_realizations)

    if (master) then
      ! Open global tallies dataset
      call h5dopen_f(hdf5_state_point, "global_tallies", dset, hdf5_err)

      ! Read global tallies
      f_ptr = c_loc(global_tallies(1))
      call h5dread_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)

      ! Close global tallies dataset
      call h5dclose_f(dset, hdf5_err)

      TALLIES_LOOP: do i = 1, n_tallies
        ! Open tally group
        call h5gopen_f(hdf5_state_point, "tallies/tally " // &
             to_str(tallies(i) % id), tally_group, hdf5_err)

        ! Read number of realizations
        call hdf5_read_integer(tally_group, "n_realizations", &
             tallies(i) % n_realizations)

        ! Open dataset for tally results
        call h5dopen_f(tally_group, "results", dset, hdf5_err)

        ! Read sum and sum_sq for each tally bin
        f_ptr = c_loc(tallies(i) % results(1,1))
        call h5dread_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)

        ! Close dataset for tally results
        call h5dclose_f(dset, hdf5_err)

        ! Close tally group
        call h5gclose_f(tally_group, hdf5_err)
      end do TALLIES_LOOP
    end if

    ! TODO: Use parallel HDF5 to read source bank in parallel

    ! Open dataset for source bank
    call h5dopen_f(hdf5_state_point, "source_bank", dset, hdf5_err)

    ! Read source bank
    f_ptr = c_loc(source_bank(1))
    call h5dread_f(dset, hdf5_bank_t, f_ptr, hdf5_err)

    ! Close dataset for source bank
    call h5dclose_f(dset, hdf5_err)

    ! Close HDF5 state point file
    call h5fclose_f(hdf5_state_point, hdf5_err)

  end subroutine hdf5_load_state_point

!===============================================================================
! HDF5_WRITE_INTEGER
!===============================================================================

  subroutine hdf5_write_integer(group, name, buffer)

    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    integer,        intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltmake_dataset_int_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_write_integer

!===============================================================================
! HDF5_WRITE_LONG
!===============================================================================

  subroutine hdf5_write_long(group, name, buffer)

    integer(HID_T),     intent(in) :: group
    character(*),       intent(in) :: name
    integer(8), target, intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
    integer(HID_T)   :: dspace
    integer(HID_T)   :: dset
    type(c_ptr)      :: f_ptr

    ! Create dataspace and dataset
    call h5screate_simple_f(rank, dims, dspace, hdf5_err)
    call h5dcreate_f(group, name, hdf5_integer8_t, dspace, dset, hdf5_err)

    ! Write eight-byte integer
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, hdf5_integer8_t, f_ptr, hdf5_err)

    ! Close dataspace and dataset for long integer
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

  end subroutine hdf5_write_long

!===============================================================================
! HDF5_WRITE_DOUBLE
!===============================================================================

  subroutine hdf5_write_double(group, name, buffer)

    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    real(8),        intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltmake_dataset_double_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_write_double

!===============================================================================
! HDF5_READ_INTEGER
!===============================================================================

  subroutine hdf5_read_integer(group, name, buffer)

    integer(HID_T), intent(in)    :: group
    character(*),   intent(in)    :: name
    integer,        intent(inout) :: buffer

    integer          :: buffer_copy(1)
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltread_dataset_int_f(group, name, buffer_copy, dims, hdf5_err)
    buffer = buffer_copy(1)

  end subroutine hdf5_read_integer

!===============================================================================
! HDF5_READ_LONG
!===============================================================================

  subroutine hdf5_read_long(group, name, buffer)

    integer(HID_T),     intent(in)  :: group
    character(*),       intent(in)  :: name
    integer(8), target, intent(out) :: buffer

    integer(HID_T) :: dset
    type(c_ptr)    :: f_ptr

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Get pointer to buffer
    f_ptr = c_loc(buffer)

    ! Read data from dataset
    call h5dread_f(dset, hdf5_integer8_t, f_ptr, hdf5_err)

    ! Close dataset
    call h5dclose_f(dset, hdf5_err)

  end subroutine hdf5_read_long

!===============================================================================
! HDF5_READ_DOUBLE
!===============================================================================

  subroutine hdf5_read_double(group, name, buffer)

    integer(HID_T), intent(in)  :: group
    character(*),   intent(in)  :: name
    real(8),        intent(out) :: buffer

    real(8)          :: buffer_copy(1)
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltread_dataset_double_f(group, name, buffer_copy, dims, hdf5_err)
    buffer = buffer_copy(1)

  end subroutine hdf5_read_double

#endif

end module hdf5_interface
