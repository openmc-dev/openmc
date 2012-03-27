module hdf5_interface

  use ace_header,      only: Reaction, UrrData
  use constants
  use endf,            only: reaction_name
  use geometry_header, only: Cell, Surface, Universe, Lattice
  use global
  use material_header, only: Material
  use string,          only: to_str

#ifdef HDF5
  use hdf5
  use h5lt
#endif

  implicit none

#ifdef HDF5

contains

!===============================================================================
! HDF5_CREATE_OUTPUT
!===============================================================================

  subroutine hdf5_create_output()

    character(MAX_FILE_LEN) :: filename ! File name
    integer :: error  ! Error flag

    ! set output file name
    filename = trim(path_input) // "output.h5"

    ! Initialize FORTRAN interface.
    call h5open_f (error)

    ! Create a new file using default properties.
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, hdf5_output_file, error)

  end subroutine hdf5_create_output

!===============================================================================
! HDF5_OPEN_OUTPUT
!===============================================================================

  subroutine hdf5_open_output()

    character(MAX_FILE_LEN) :: filename ! File name
    integer :: error ! Error flag

    ! set output file name
    filename = trim(path_input) // "output.h5"

    ! Initialize FORTRAN interface.
    call h5open_f(hdf5_err)

    ! Create a new file using default properties.
    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, hdf5_output_file, error)

  end subroutine hdf5_open_output

!===============================================================================
! HDF5_WRITE_HEADER
!===============================================================================

  subroutine hdf5_write_header()

    character(8)  :: date_
    character(10) :: time_
    character(19) :: current_time
     
    ! Write version information
    call hdf5_make_integer(hdf5_output_file, "version_major", VERSION_MAJOR)
    call hdf5_make_integer(hdf5_output_file, "version_minor", VERSION_MINOR)
    call hdf5_make_integer(hdf5_output_file, "version_release", VERSION_RELEASE)

    ! Write current date and time
    call date_and_time(DATE=date_, TIME=time_)
    current_time = date_(1:4) // "-" // date_(5:6) // "-" // date_(7:8) // &
         " " // time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)
    call h5ltmake_dataset_string_f(hdf5_output_file, "date_and_time", &
         current_time, hdf5_err)

    ! Write MPI information
    call hdf5_make_integer(hdf5_output_file, "n_procs", n_procs)
    call h5ltset_attribute_string_f(hdf5_output_file, "n_procs", &
         "description", "Number of MPI processes", hdf5_err)

  end subroutine hdf5_write_header

!===============================================================================
! HDF5_WRITE_SUMMARY
!===============================================================================

  subroutine hdf5_write_summary()

    ! Write criticality information
    if (problem_type == PROB_CRITICALITY) then
       ! Need to write integer(8)'s using double instead since there is no H5LT
       ! call for making a dataset of type long
       call hdf5_make_double(hdf5_output_file, "n_particles", real(n_particles,8))

       ! Use H5LT interface to write n_cycles, n_inactive, and n_active
       call hdf5_make_integer(hdf5_output_file, "n_batches", n_batches)
       call hdf5_make_integer(hdf5_output_file, "n_inactive", n_inactive)
       call hdf5_make_integer(hdf5_output_file, "n_active", n_active)
       call hdf5_make_integer(hdf5_output_file, "gen_per_batch", gen_per_batch)

       ! Add description of each variable
       call h5ltset_attribute_string_f(hdf5_output_file, "n_particles", &
            "description", "Number of particles per cycle", hdf5_err)
       call h5ltset_attribute_string_f(hdf5_output_file, "n_batches", &
            "description", "Total number of batches", hdf5_err)
       call h5ltset_attribute_string_f(hdf5_output_file, "n_inactive", &
            "description", "Number of inactive cycles", hdf5_err)
       call h5ltset_attribute_string_f(hdf5_output_file, "n_active", &
            "description", "Number of active cycles", hdf5_err)
       call h5ltset_attribute_string_f(hdf5_output_file, "gen_per_batch", &
            "description", "Number of generations per batch", hdf5_err)
    end if

    call hdf5_write_geometry()
    call hdf5_write_materials()
    call hdf5_write_nuclides()
    if (n_tallies > 0) then
       call hdf5_write_tallies()
    end if

  end subroutine hdf5_write_summary

!===============================================================================
! HDF5_WRITE_RESULTS
!===============================================================================

  subroutine hdf5_write_results()

    call hdf5_write_timing()
    call hdf5_write_global_tallies()

  end subroutine hdf5_write_results

!===============================================================================
! HDF5_WRITE_GEOMETRY
!===============================================================================

  subroutine hdf5_write_geometry()

    integer          :: i, j, k
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: dims2(2)
    integer(HID_T)   :: geometry_group
    integer(HID_T)   :: cell_group
    integer(HID_T)   :: surface_group
    integer(HID_T)   :: universe_group
    integer(HID_T)   :: lattice_group
    integer(HID_T)   :: temp_group
    integer, allocatable :: lattice_universes(:,:)
    type(Cell),     pointer :: c => null()
    type(Surface),  pointer :: s => null()
    type(Universe), pointer :: u => null()
    type(Lattice),  pointer :: l => null()
     
    ! Create group for geometry
    call h5gcreate_f(hdf5_output_file, "/geometry", geometry_group, hdf5_err)

    ! Use H5LT interface to write number of geometry objects
    call hdf5_make_integer(geometry_group, "n_cells", n_cells)
    call hdf5_make_integer(geometry_group, "n_surfaces", n_surfaces)
    call hdf5_make_integer(geometry_group, "n_universes", n_universes)
    call hdf5_make_integer(geometry_group, "n_lattices", n_lattices)

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
       call hdf5_make_integer(temp_group, "universe", &
            universes(c % universe) % id)
       
       ! Write information on what fills this cell
       select case (c % type)
       case (CELL_NORMAL)
          call h5ltmake_dataset_string_f(temp_group, "fill_type", "normal", &
               hdf5_err) 
          call hdf5_make_integer(temp_group, "material", &
               materials(c % material) % id)
       case (CELL_FILL)
          call h5ltmake_dataset_string_f(temp_group, "fill_type", "universe", &
               hdf5_err) 
          call hdf5_make_integer(temp_group, "material", &
               universes(c % fill) % id)
       case (CELL_LATTICE)
          call h5ltmake_dataset_string_f(temp_group, "fill_type", "lattice", &
               hdf5_err) 
          call hdf5_make_integer(temp_group, "lattice", &
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
       case (SURF_BOX_X)
       case (SURF_BOX_Y)
       case (SURF_BOX_Z)
       case (SURF_BOX)
       case (SURF_GQ)
          call h5ltmake_dataset_string_f(temp_group, "type", "General Quadratic", hdf5_err)
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
       l => lattices(i)

       ! Create group for i-th lattice
       call h5gcreate_f(lattice_group, "lattice " // trim(to_str(l % id)), &
            temp_group, hdf5_err)

       ! Write lattice type
       select case(l % type)
       case (LATTICE_RECT)
          call h5ltmake_dataset_string_f(temp_group, "type", "rectangular", hdf5_err)
       case (LATTICE_HEX)
          call h5ltmake_dataset_string_f(temp_group, "type", "hexagonal", hdf5_err)
       end select

       ! Write lattice dimensions, lower left corner, and width of element
       dims(1) = 2
       call h5ltmake_dataset_int_f(temp_group, "n_elements", 1, dims, &
            (/ l % n_x, l % n_y /), hdf5_err)
       call h5ltmake_dataset_double_f(temp_group, "lower_left", 1, dims, &
            (/ l % x0, l % y0 /), hdf5_err)
       call h5ltmake_dataset_double_f(temp_group, "element_width", 1, dims, &
            (/ l % width_x, l % width_y /), hdf5_err)

       ! Write lattice elements
       allocate(lattice_universes(l % n_x, l % n_y))
       do j = 1, l % n_x
          do k = 1, l % n_y
             lattice_universes(j,k) = universes(l % element(j,k)) % id
          end do
       end do
       dims2 = (/ l % n_x, l % n_y /)
       call h5ltmake_dataset_int_f(temp_group, "elements", 2, dims2, &
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
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: size_string = 12
    integer(HID_T)   :: materials_group
    integer(HID_T)   :: temp_group
    integer(HID_T)   :: string_type
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: dset_id
    type(Material), pointer :: m => null()
     
    ! Create group for materials
    call h5gcreate_f(hdf5_output_file, "/materials", materials_group, hdf5_err)

    ! Use H5LT interface to write number of materials
    call hdf5_make_integer(materials_group, "n_materials", n_materials)

    ! Write information on each material
    do i = 1, n_materials
       m => materials(i)

       ! Create group for i-th universe
       call h5gcreate_f(materials_group, "material " // trim(to_str(m % id)), &
            temp_group, hdf5_err)

       ! Write atom density with units
       call hdf5_make_double(temp_group, "atom_density", m % density)
       call h5ltset_attribute_string_f(temp_group, "atom_density", &
            "units", "atom/barn-cm", hdf5_err)

       ! Create string type of length 12
       call h5tcopy_f(H5T_C_S1, string_type, hdf5_err)
       call h5tset_size_f(string_type, size_string, hdf5_err)

       ! Create dataspace and dataset for writing nuclides
       dims(1) = m % n_nuclides
       call h5screate_simple_f(1, dims, dspace_id, hdf5_err)
       call h5dcreate_f(temp_group, "nuclides", string_type, dspace_id, &
            dset_id, hdf5_err)

       ! Write list of nuclides
       call h5dwrite_f(dset_id, string_type, m % names, dims, hdf5_err)

       ! Close dataspace and dataset for nuclides
       call h5dclose_f(dset_id, hdf5_err)
       call h5sclose_f(dspace_id, hdf5_err)

       ! Write atom densities
       call h5ltmake_dataset_double_f(temp_group, "nuclide_densities", 1, &
            dims, m % atom_density, hdf5_err)

       ! Write S(a,b) information if present
       if (m % has_sab_table) then
          call h5ltmake_dataset_string_f(temp_group, "sab_table", &
               m % sab_name, hdf5_err)
       end if

       ! Close group for i-th material
       call h5gclose_f(temp_group, hdf5_err)
    end do

    ! Close materials group
    call h5gclose_f(materials_group, hdf5_err)

  end subroutine hdf5_write_materials

!===============================================================================
! HDF5_WRITE_GLOBAL_TALLIES
!===============================================================================

  subroutine hdf5_write_global_tallies()

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/ 2 /)

    call h5ltmake_dataset_double_f(hdf5_output_file, "k_analog", &
         rank, dims, (/ global_tallies(K_ANALOG) % sum, &
         global_tallies(K_ANALOG) % sum_sq /), hdf5_err)
    call h5ltmake_dataset_double_f(hdf5_output_file, "k_collision", &
         rank, dims, (/ global_tallies(K_COLLISION) % sum, &
         global_tallies(K_COLLISION) % sum_sq /), hdf5_err)
    call h5ltmake_dataset_double_f(hdf5_output_file, "k_tracklength", &
         rank, dims, (/ global_tallies(K_TRACKLENGTH) % sum, &
         global_tallies(K_TRACKLENGTH) % sum_sq /), hdf5_err)
    call h5ltmake_dataset_double_f(hdf5_output_file, "leakage", &
         rank, dims, (/ global_tallies(LEAKAGE) % sum, &
         global_tallies(LEAKAGE) % sum_sq /), hdf5_err)

  end subroutine hdf5_write_global_tallies

!===============================================================================
! HDF5_WRITE_TALLIES
!===============================================================================

  subroutine hdf5_write_tallies()

    integer           :: i, j
    integer(SIZE_T)   :: num_elements
    integer(HSIZE_T)  :: dims(1)
    integer(HSIZE_T)  :: coord(1,1)
    integer(HSIZE_T)  :: size_string = 12
    integer(HID_T)    :: string_type
    integer(HID_T)    :: tallies_group
    integer(HID_T)    :: temp_group
    integer(HID_T)    :: dspace
    integer(HID_T)    :: subspace
    integer(HID_T)    :: dset
    character(12)     :: string
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()
     
    ! Create group for tallies
    call h5gcreate_f(hdf5_output_file, "/tallies", tallies_group, hdf5_err)

    ! Use H5LT interface to write number of tallies
    call hdf5_make_integer(tallies_group, "n_tallies", n_tallies)

    ! Write information on each material
    do i = 1, n_tallies
       t => tallies(i)

       ! Create group for i-th universe
       call h5gcreate_f(tallies_group, "tally " // trim(to_str(t % id)), &
            temp_group, hdf5_err)

       ! =======================================================================
       ! WRITE INFORMATION ON TALLY FILTERS

       ! Write filters
       call hdf5_make_integer(temp_group, "n_filters", t % n_filters)
       if (t % n_filters > 0) then
          dims(1) = t % n_filters
          call h5ltmake_dataset_int_f(temp_group, "filters", 1, &
               dims, t % filters, hdf5_err)
       end if

       ! Write universe_bins if present
       if (t % n_filter_bins(FILTER_UNIVERSE) > 0) then
          dims(1) = t % n_filter_bins(FILTER_UNIVERSE)
          call h5ltmake_dataset_int_f(temp_group, "universe_bins", 1, &
               dims, t % universe_bins(:) % scalar, hdf5_err)
       end if

       ! Write material_bins if present
       if (t % n_filter_bins(FILTER_MATERIAL) > 0) then
          dims(1) = t % n_filter_bins(FILTER_MATERIAL)
          call h5ltmake_dataset_int_f(temp_group, "material_bins", 1, &
               dims, t % material_bins(:) % scalar, hdf5_err)
       end if

       ! Write cell_bins if present
       if (t % n_filter_bins(FILTER_CELL) > 0) then
          dims(1) = t % n_filter_bins(FILTER_CELL)
          call h5ltmake_dataset_int_f(temp_group, "cell_bins", 1, &
               dims, t % cell_bins(:) % scalar, hdf5_err)
       end if

       ! Write cellborn_bins if present
       if (t % n_filter_bins(FILTER_CELLBORN) > 0) then
          dims(1) = t % n_filter_bins(FILTER_CELLBORN)
          call h5ltmake_dataset_int_f(temp_group, "cellborn_bins", 1, &
               dims, t % cellborn_bins(:) % scalar, hdf5_err)
       end if

       ! Write surface_bins if present
       if (t % n_filter_bins(FILTER_SURFACE) > 0) then
          dims(1) = t % n_filter_bins(FILTER_SURFACE)
          call h5ltmake_dataset_int_f(temp_group, "surface_bins", 1, &
               dims, t % surface_bins(:) % scalar, hdf5_err)
       end if

       ! Write incoming energy filter
       if (t % n_filter_bins(FILTER_ENERGYIN) > 0) then
          dims(1) = t % n_filter_bins(FILTER_ENERGYIN) + 1
          call h5ltmake_dataset_double_f(temp_group, "energy_in", 1, &
               dims, t % energy_in, hdf5_err)
       end if

       ! Write outgoing energy filter
       if (t % n_filter_bins(FILTER_ENERGYOUT) > 0) then
          dims(1) = t % n_filter_bins(FILTER_ENERGYOUT) + 1
          call h5ltmake_dataset_double_f(temp_group, "energy_out", 1, &
               dims, t % energy_out, hdf5_err)
       end if

       ! Write mesh information
       if (t % n_filter_bins(FILTER_MESH) > 0) then
          m => meshes(t % mesh)
          dims(1) = m % n_dimension

          ! Write mesh dimensions
          call h5ltmake_dataset_int_f(temp_group, "mesh_dimensions", 1, &
               dims, m % dimension, hdf5_err)

          ! Write mesh lower-left corner
          call h5ltmake_dataset_double_f(temp_group, "mesh_lower_left", 1, &
               dims, m % lower_left, hdf5_err)
          
          ! Write mesh element width
          call h5ltmake_dataset_double_f(temp_group, "mesh_element_width", 1, &
               dims, m % width, hdf5_err)
       end if
          
       ! =======================================================================
       ! WRITE INFORMATION ON TALLY SCORES

       ! Create string type of length 12
       call h5tcopy_f(H5T_C_S1, string_type, hdf5_err)
       call h5tset_size_f(string_type, size_string, hdf5_err)

       ! Create dataspace and dataset for writing nuclides
       dims(1) = t % n_score_bins
       call h5screate_simple_f(1, dims, dspace, hdf5_err)
       call h5dcreate_f(temp_group, "scores", string_type, dspace, &
            dset, hdf5_err)

       ! Create subspace for writing one element
       dims(1) = 1
       call h5screate_simple_f(1, dims, subspace, hdf5_err)

       ! Write list of nuclides
       num_elements = 1
       do j = 1, t % n_score_bins
          ! Determine string for this scoring bin
          select case (t % score_bins(j) % scalar)
          case (SCORE_FLUX)
             string = 'flux'
          case (SCORE_TOTAL)
             string = 'total'
          case (SCORE_SCATTER)
             string = 'scatter'
          case (SCORE_NU_SCATTER)
             string = 'nu-scatter'
          case (SCORE_SCATTER_1)
             string = '1st moment'
          case (SCORE_SCATTER_2)
             string = '2nd moment'
          case (SCORE_SCATTER_3)
             string = '3rd moment'
          case (SCORE_N_1N)
             string = '(n,1n)'
          case (SCORE_N_2N)
             string = '(n,2n)'
          case (SCORE_N_3N)
             string = '(n,3n)'
          case (SCORE_N_4N)
             string = '(n,4n)'
          case (SCORE_ABSORPTION)
             string = 'absorption'
          case (SCORE_FISSION)
             string = 'fission'
          case (SCORE_NU_FISSION)
             string = 'nu-fission'
          case (SCORE_CURRENT)
             string = 'current'
          end select

          ! Select block within array
          coord(1,1) = j
          call h5sselect_elements_f(dspace, H5S_SELECT_SET_F, 1, &
               num_elements, coord, hdf5_err)

          ! Write data to partial array
          call h5dwrite_f(dset, string_type, string, dims, hdf5_err, &
               subspace, dspace)
       end do

       ! Close dataspace and dataset for nuclides
       call h5dclose_f(dset, hdf5_err)
       call h5sclose_f(subspace, hdf5_err)
       call h5sclose_f(dspace, hdf5_err)

       ! Close group for i-th material
       call h5gclose_f(temp_group, hdf5_err)
    end do

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
    call hdf5_make_integer(group, "n_nuclides", n_nuclides_total)

    ! Write information on each nuclide
    do i = 1, n_nuclides_total
       nuc => nuclides(i)

       ! Determine size of cross-sections
       size_xs = (5 + nuc % n_reaction) * nuc % n_grid * 8
       size_total = size_xs

       ! Create group for i-th nuclide
       call h5gcreate_f(group, trim(nuc % name), nuclide_group, hdf5_err)

       ! Write some basic attributes
       call hdf5_make_integer(nuclide_group, "zaid", nuc % zaid)
       call hdf5_make_double(nuclide_group, "awr", nuc % awr)
       call hdf5_make_double(nuclide_group, "kT", nuc % kT)
       call hdf5_make_integer(nuclide_group, "n_grid", nuc % n_grid)
       call hdf5_make_integer(nuclide_group, "n_reactions", nuc % n_reaction)
       call hdf5_make_integer(nuclide_group, "n_fission", nuc % n_fission)
       call hdf5_make_integer(nuclide_group, "size_xs", size_xs)

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
          call hdf5_make_double(rxn_group, "Q_value", rxn % Q_value)
          call hdf5_make_integer(rxn_group, "n_neutrons", rxn % TY)
          call hdf5_make_integer(rxn_group, "grid_index", rxn % IE)
          call hdf5_make_integer(rxn_group, "size_angle", size_angle)
          call hdf5_make_integer(rxn_group, "size_energy", size_energy)

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
          call hdf5_make_integer(nuclide_group, "urr_n_energy", urr % n_energy)
          call hdf5_make_integer(nuclide_group, "urr_n_prob", urr % n_prob)
          call hdf5_make_integer(nuclide_group, "urr_interp", urr % interp)
          call hdf5_make_integer(nuclide_group, "urr_inelastic", urr % inelastic_flag)
          call hdf5_make_integer(nuclide_group, "urr_absorption", urr % absorption_flag)
          call hdf5_make_double(nuclide_group, "urr_min_E", urr % energy(1))
          call hdf5_make_double(nuclide_group, "urr_max_E", urr % energy(urr % n_energy))
       end if

       ! Write total memory used
       call hdf5_make_integer(nuclide_group, "size_total", size_total)

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
    call hdf5_make_double(timing_group, "time_initialize", time_initialize % elapsed)
    call hdf5_make_double(timing_group, "time_read_xs", time_read_xs % elapsed)
    call hdf5_make_double(timing_group, "time_unionize", time_unionize % elapsed)
    call hdf5_make_double(timing_group, "time_transport", time_transport % elapsed)
    call hdf5_make_double(timing_group, "time_intercycle", time_intercycle % elapsed)
    call hdf5_make_double(timing_group, "time_tallies", time_ic_tallies % elapsed)
    call hdf5_make_double(timing_group, "time_sample", time_ic_sample % elapsed)
    call hdf5_make_double(timing_group, "time_sendrecv", time_ic_sendrecv % elapsed)
    call hdf5_make_double(timing_group, "time_inactive", time_inactive % elapsed)
    call hdf5_make_double(timing_group, "time_active", time_active % elapsed)
    call hdf5_make_double(timing_group, "time_finalize", time_finalize % elapsed)
    call hdf5_make_double(timing_group, "time_total", time_total % elapsed)

    ! Add descriptions to timing data
    call h5ltset_attribute_string_f(timing_group, "time_initialize", &
         "description", "Total time elapsed for initialization (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_read_xs", &
         "description", "Time reading cross-section libraries (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_unionize", &
         "description", "Time unionizing energy grid (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_transport", &
         "description", "Time in transport only (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_intercycle", &
         "description", "Total time between generations (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_tallies", &
         "description", "Time between cycles accumulating tallies (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_sample", &
         "description", "Time between cycles sampling source sites (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_sendrecv", &
         "description", "Time between cycles SEND/RECVing source sites (s)", hdf5_err)
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
    call hdf5_make_double(timing_group, "neutrons_per_second", speed)

    ! Close timing group
    call h5gclose_f(timing_group, hdf5_err)

  end subroutine hdf5_write_timing

!===============================================================================
! HDF5_CLOSE_OUTPUT
!===============================================================================

  subroutine hdf5_close_output()

    ! Terminate access to the file.
    call h5fclose_f(hdf5_output_file, hdf5_err)

    ! Close FORTRAN interface.
    call h5close_f(hdf5_err)

  end subroutine hdf5_close_output

!===============================================================================
! HDF5_MAKE_INTEGER
!===============================================================================

  subroutine hdf5_make_integer(group, name, buffer)

    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    integer,        intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltmake_dataset_int_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_make_integer

!===============================================================================
! HDF5_MAKE_DOUBLE
!===============================================================================

  subroutine hdf5_make_double(group, name, buffer)

    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    real(8),        intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltmake_dataset_double_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_make_double

#endif

end module hdf5_interface
