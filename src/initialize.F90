module initialize

  use ace,              only: read_xs
  use bank_header,      only: Bank
  use constants
  use dict_header,      only: DictIntInt, ElemKeyValueII
  use energy_grid,      only: unionized_grid
  use error,            only: fatal_error, warning
  use geometry,         only: neighbor_lists
  use geometry_header,  only: Cell, Universe, Lattice, BASE_UNIVERSE
  use global
  use input_xml,        only: read_input_xml, read_cross_sections_xml,         &
                              cells_in_univ_dict, read_plots_xml
  use output,           only: title, header, write_summary, print_version,     &
                              print_usage, write_xs_summary, print_plot,       &
                              write_message
  use output_interface
  use random_lcg,       only: initialize_prng
  use source,           only: initialize_source
  use state_point,      only: load_state_point
  use string,           only: to_str, str_to_int, starts_with, ends_with
  use tally_header,     only: TallyObject, TallyResult
  use tally_initialize, only: configure_tallies

#ifdef MPI
  use mpi
#endif

#ifdef _OPENMP
  use omp_lib
#endif

#ifdef HDF5
  use hdf5_interface
  use hdf5_summary,     only: hdf5_write_summary
#endif

  implicit none

contains

!===============================================================================
! INITIALIZE_RUN takes care of all initialization tasks, i.e. reading
! from command line, reading xml input files, initializing random
! number seeds, reading cross sections, initializing starting source,
! setting up timers, etc.
!===============================================================================

  subroutine initialize_run()

    ! Start total and initialization timer
    call time_total % start()
    call time_initialize % start()

#ifdef MPI
    ! Setup MPI
    call initialize_mpi()
#endif

#ifdef HDF5
    ! Initialize HDF5 interface
    call hdf5_initialize()
#endif

    ! Read command line arguments
    call read_command_line()

    if (master) then
      ! Display title and initialization header
      call title()
      call header("INITIALIZATION", level=1)
    end if

    ! Read XML input files
    call read_input_xml()

    ! Initialize random number generator -- this has to be done after the input
    ! files have been read in case the user specified a seed for the random
    ! number generator

    call initialize_prng()

    ! Read plots.xml if it exists -- this has to be done separate from the other
    ! XML files because we need the PRNG to be initialized first
    if (run_mode == MODE_PLOTTING) call read_plots_xml()

    ! Set up universe structures
    call prepare_universes()

    ! Use dictionaries to redefine index pointers
    call adjust_indices()

    ! After reading input and basic geometry setup is complete, build lists of
    ! neighboring cells for efficient tracking
    call neighbor_lists()

    if (run_mode /= MODE_PLOTTING) then
      ! With the AWRs from the xs_listings, change all material specifications
      ! so that they contain atom percents summing to 1
      call normalize_ao()

      ! Read ACE-format cross sections
      call time_read_xs % start()
      call read_xs()
      call time_read_xs % stop()

      ! Construct unionized energy grid from cross-sections
      if (grid_method == GRID_UNION) then
        call time_unionize % start()
        call unionized_grid()
        call time_unionize % stop()
      end if

      ! Allocate and setup tally stride, matching_bins, and tally maps
      call configure_tallies()

      ! Determine how much work each processor should do
      call calculate_work()

      ! Allocate banks and create source particles -- for a fixed source
      ! calculation, the external source distribution is sampled during the
      ! run, not at initialization
      if (run_mode == MODE_EIGENVALUE) then
        call allocate_banks()
        if (.not. restart_run) call initialize_source()
      end if

      ! If this is a restart run, load the state point data and binary source
      ! file
      if (restart_run) call load_state_point()
    end if

    if (master) then
      if (run_mode == MODE_PLOTTING) then
        ! Display plotting information
        call print_plot()
      else
        ! Write summary information
#ifdef HDF5
        if (output_summary) call hdf5_write_summary()
#else
        if (output_summary) call write_summary()
#endif

        ! Write cross section information
        if (output_xs) call write_xs_summary()
      end if
    end if

    ! Check for particle restart run
    if (particle_restart_run) run_mode = MODE_PARTICLE

    ! Warn if overlap checking is on
    if (master .and. check_overlaps) then
      message = ""
      call write_message()
      message = "Cell overlap checking is ON"
      call warning()
    end if

    ! Stop initialization timer
    call time_initialize % stop()

  end subroutine initialize_run

#ifdef MPI
!===============================================================================
! INITIALIZE_MPI starts up the Message Passing Interface (MPI) and determines
! the number of processors the problem is being run with as well as the rank of
! each processor.
!===============================================================================

  subroutine initialize_mpi()

    integer                   :: bank_blocks(4)  ! Count for each datatype
    integer                   :: bank_types(4)   ! Datatypes
    integer(MPI_ADDRESS_KIND) :: bank_disp(4)    ! Displacements
    integer                   :: temp_type       ! temporary derived type
    integer                   :: result_blocks(1) ! Count for each datatype
    integer                   :: result_types(1)  ! Datatypes
    integer(MPI_ADDRESS_KIND) :: result_disp(1)   ! Displacements
    integer(MPI_ADDRESS_KIND) :: result_base_disp ! Base displacement
    integer(MPI_ADDRESS_KIND) :: lower_bound     ! Lower bound for TallyResult
    integer(MPI_ADDRESS_KIND) :: extent          ! Extent for TallyResult
    type(Bank)       :: b
    type(TallyResult) :: tr

    ! Indicate that MPI is turned on
    mpi_enabled = .true.

    ! Initialize MPI
    call MPI_INIT(mpi_err)

    ! Determine number of processors and rank of each processor
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

    ! Determine master
    if (rank == 0) then
      master = .true.
    else
      master = .false.
    end if

    ! ==========================================================================
    ! CREATE MPI_BANK TYPE

    ! Determine displacements for MPI_BANK type
    call MPI_GET_ADDRESS(b % wgt, bank_disp(1), mpi_err)
    call MPI_GET_ADDRESS(b % xyz, bank_disp(2), mpi_err)
    call MPI_GET_ADDRESS(b % uvw, bank_disp(3), mpi_err)
    call MPI_GET_ADDRESS(b % E,   bank_disp(4), mpi_err)

    ! Adjust displacements 
    bank_disp = bank_disp - bank_disp(1)

    ! Define MPI_BANK for fission sites
    bank_blocks = (/ 1, 3, 3, 1 /)
    bank_types = (/ MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8 /)
    call MPI_TYPE_CREATE_STRUCT(4, bank_blocks, bank_disp, & 
         bank_types, MPI_BANK, mpi_err)
    call MPI_TYPE_COMMIT(MPI_BANK, mpi_err)

    ! ==========================================================================
    ! CREATE MPI_TALLYRESULT TYPE

    ! Determine displacements for MPI_BANK type
    call MPI_GET_ADDRESS(tr % value, result_base_disp, mpi_err)
    call MPI_GET_ADDRESS(tr % sum, result_disp(1), mpi_err)

    ! Adjust displacements
    result_disp = result_disp - result_base_disp

    ! Define temporary type for TallyResult
    result_blocks = (/ 2 /)
    result_types = (/ MPI_REAL8 /)
    call MPI_TYPE_CREATE_STRUCT(1, result_blocks, result_disp, result_types, &
         temp_type, mpi_err)

    ! Adjust lower-bound and extent of type for tally score
    lower_bound = 0
    extent      = result_disp(1) + 16
    call MPI_TYPE_CREATE_RESIZED(temp_type, lower_bound, extent, &
         MPI_TALLYRESULT, mpi_err)

    ! Commit derived type for tally scores
    call MPI_TYPE_COMMIT(MPI_TALLYRESULT, mpi_err)

    ! Free temporary MPI type
    call MPI_TYPE_FREE(temp_type, mpi_err)

  end subroutine initialize_mpi
#endif

#ifdef HDF5

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

#endif

!===============================================================================
! READ_COMMAND_LINE reads all parameters from the command line
!===============================================================================

  subroutine read_command_line()

    integer :: i         ! loop index
    integer :: argc      ! number of command line arguments
    integer :: last_flag ! index of last flag
    integer :: filetype
    character(MAX_FILE_LEN) :: pwd      ! present working directory
    character(MAX_WORD_LEN), allocatable :: argv(:) ! command line arguments
    type(BinaryOutput) :: sp

    ! Get working directory
    call GET_ENVIRONMENT_VARIABLE("PWD", pwd)

    ! Check number of command line arguments and allocate argv
    argc = COMMAND_ARGUMENT_COUNT()

    ! Allocate and retrieve command arguments
    allocate(argv(argc))
    do i = 1, argc
      call GET_COMMAND_ARGUMENT(i, argv(i))
    end do

    ! Process command arguments
    last_flag = 0
    i = 1
    do while (i <= argc)
      ! Check for flags
      if (starts_with(argv(i), "-")) then
        select case (argv(i))
        case ('-p', '-plot', '--plot')
          run_mode = MODE_PLOTTING
          check_overlaps = .true.

        case ('-n', '-particles', '--particles')
          ! Read number of particles per cycle
          i = i + 1
          n_particles = str_to_int(argv(i))

          ! Check that number specified was valid
          if (n_particles == ERROR_INT) then
            message = "Must specify integer after " // trim(argv(i-1)) // &
                 " command-line flag."
            call fatal_error()
          end if
        case ('-r', '-restart', '--restart')
          ! Read path for state point/particle restart
          i = i + 1

          ! Check what type of file this is
          call sp % file_open(argv(i), 'r', serial = .false.)
          call sp % read_data(filetype, 'filetype')
          call sp % file_close()

          ! Set path and flag for type of run
          select case (filetype)
          case (FILETYPE_STATEPOINT)
            path_state_point = argv(i)
            restart_run = .true.
          case (FILETYPE_PARTICLE_RESTART)
            path_particle_restart = argv(i)
            particle_restart_run = .true.
          case default
            message = "Unrecognized file after restart flag."
            call fatal_error()
          end select

          ! If its a restart run check for additional source file
          if (restart_run .and. i + 1 <= argc) then

            ! Increment arg
            i = i + 1

            ! Check if it has extension we can read
            if ((ends_with(argv(i), '.binary') .or. &
                 ends_with(argv(i), '.h5'))) then

              ! Check file type is a source file
              call sp % file_open(argv(i), 'r', serial = .false.)
              call sp % read_data(filetype, 'filetype')
              call sp % file_close()
              if (filetype /= FILETYPE_SOURCE) then
                message = "Second file after restart flag must be a source file"
                call fatal_error()
              end if

              ! It is a source file
              path_source_point = argv(i)

            else ! Different option is specified not a source file

              ! Source is in statepoint file
              path_source_point = path_state_point

              ! Set argument back
              i = i - 1

            end if

          else ! No command line arg after statepoint

            ! Source is assumed to be in statepoint file
            path_source_point = path_state_point

          end if

        case ('-g', '-geometry-debug', '--geometry-debug')
          check_overlaps = .true.

        case ('-s', '--threads')
          ! Read number of threads
          i = i + 1

#ifdef _OPENMP          
          ! Read and set number of OpenMP threads
          n_threads = str_to_int(argv(i))
          if (n_threads < 1) then
            message = "Invalid number of threads specified on command line."
            call fatal_error()
          end if
          call omp_set_num_threads(n_threads)
#else
          message = "Ignoring number of threads specified on command line."
          call warning()
#endif

        case ('-?', '-help', '--help')
          call print_usage()
          stop
        case ('-v', '-version', '--version')
          call print_version()
          stop
        case ('-eps_tol', '-ksp_gmres_restart')
          ! Handle options that would be based to PETSC
          i = i + 1
        case ('-t', '-track', '--track')
          write_all_tracks = .true.
          i = i + 1
        case default
          message = "Unknown command line option: " // argv(i)
          call fatal_error()
        end select

        last_flag = i
      end if

      ! Increment counter
      i = i + 1
    end do

    ! Determine directory where XML input files are
    if (argc > 0 .and. last_flag < argc) then
      path_input = argv(last_flag + 1)
      ! Need to add working directory if the given path is a relative path
      if (.not. starts_with(path_input, "/")) then
        path_input = trim(pwd) // "/" // trim(path_input)
      end if
    else
      path_input = pwd
    end if

    ! Add slash at end of directory if it isn't there
    if (.not. ends_with(path_input, "/")) then
      path_input = trim(path_input) // "/"
    end if

    ! Free memory from argv
    deallocate(argv)

    ! TODO: Check that directory exists

  end subroutine read_command_line

!===============================================================================
! PREPARE_UNIVERSES allocates the universes array and determines the cells array
! for each universe.
!===============================================================================

  subroutine prepare_universes()

    integer              :: i                     ! index in cells array
    integer              :: i_univ                ! index in universes array
    integer              :: n_cells_in_univ       ! number of cells in a universe
    integer, allocatable :: index_cell_in_univ(:) ! the index in the univ%cells
                                                  ! array for each universe
    type(ElemKeyValueII), pointer :: pair_list => null()
    type(ElemKeyValueII), pointer :: current => null()
    type(ElemKeyValueII), pointer :: next => null()
    type(Universe),       pointer :: univ => null()
    type(Cell),           pointer :: c => null()

    allocate(universes(n_universes))

    ! We also need to allocate the cell count lists for each universe. The logic
    ! for this is a little more convoluted. In universe_dict, the (key,value)
    ! pairs are the id of the universe and the index in the array. In
    ! cells_in_univ_dict, it's the id of the universe and the number of cells.

    pair_list => universe_dict % keys()
    current => pair_list
    do while (associated(current))
      ! Find index of universe in universes array
      i_univ = current % value
      univ => universes(i_univ)
      univ % id = current % key

      ! Check for lowest level universe
      if (univ % id == 0) BASE_UNIVERSE = i_univ

      ! Find cell count for this universe
      n_cells_in_univ = cells_in_univ_dict % get_key(univ % id)

      ! Allocate cell list for universe
      allocate(univ % cells(n_cells_in_univ))
      univ % n_cells = n_cells_in_univ

      ! Move to next universe
      next => current % next
      deallocate(current)
      current => next
    end do

    ! Also allocate a list for keeping track of where cells have been assigned
    ! in each universe

    allocate(index_cell_in_univ(n_universes))
    index_cell_in_univ = 0

    do i = 1, n_cells
      c => cells(i)

      ! Get pointer to corresponding universe
      i_univ = universe_dict % get_key(c % universe)
      univ => universes(i_univ)

      ! Increment the index for the cells array within the Universe object and
      ! then store the index of the Cell object in that array
      index_cell_in_univ(i_univ) = index_cell_in_univ(i_univ) + 1
      univ % cells(index_cell_in_univ(i_univ)) = i
    end do

    ! Clear dictionary
    call cells_in_univ_dict % clear()

  end subroutine prepare_universes

!===============================================================================
! ADJUST_INDICES changes the values for 'surfaces' for each cell and the
! material index assigned to each to the indices in the surfaces and material
! array rather than the unique IDs assigned to each surface and material. Also
! assigns boundary conditions to surfaces based on those read into the bc_dict
! dictionary
!===============================================================================

  subroutine adjust_indices()

    integer :: i             ! index for various purposes
    integer :: j             ! index for various purposes
    integer :: k             ! loop index for lattices
    integer :: m             ! loop index for lattices
    integer :: mid, lid      ! material and lattice IDs
    integer :: n_x, n_y, n_z ! size of lattice
    integer :: i_array       ! index in surfaces/materials array 
    integer :: id            ! user-specified id
    type(Cell),        pointer :: c => null()
    type(Lattice),     pointer :: lat => null()
    type(TallyObject), pointer :: t => null()

    do i = 1, n_cells
      ! =======================================================================
      ! ADJUST SURFACE LIST FOR EACH CELL

      c => cells(i)
      do j = 1, c % n_surfaces
        id = c % surfaces(j)
        if (id < OP_DIFFERENCE) then
          if (surface_dict % has_key(abs(id))) then
            i_array = surface_dict % get_key(abs(id))
            c % surfaces(j) = sign(i_array, id)
          else
            message = "Could not find surface " // trim(to_str(abs(id))) // &
                 " specified on cell " // trim(to_str(c % id))
            call fatal_error()
          end if
        end if
      end do

      ! =======================================================================
      ! ADJUST UNIVERSE INDEX FOR EACH CELL

      id = c % universe
      if (universe_dict % has_key(id)) then
        c % universe = universe_dict % get_key(id)
      else
        message = "Could not find universe " // trim(to_str(id)) // &
             " specified on cell " // trim(to_str(c % id))
        call fatal_error()
      end if

      ! =======================================================================
      ! ADJUST MATERIAL/FILL POINTERS FOR EACH CELL

      id = c % material
      if (id == MATERIAL_VOID) then
        c % type = CELL_NORMAL
      elseif (id /= 0) then
        if (material_dict % has_key(id)) then
          c % type = CELL_NORMAL
          c % material = material_dict % get_key(id)
        else
          message = "Could not find material " // trim(to_str(id)) // &
               " specified on cell " // trim(to_str(c % id))
          call fatal_error()
        end if
      else
        id = c % fill
        if (universe_dict % has_key(id)) then
          c % type = CELL_FILL
          c % fill = universe_dict % get_key(id)
        elseif (lattice_dict % has_key(id)) then
          lid = lattice_dict % get_key(id)
          mid = lattices(lid) % outside
          c % type = CELL_LATTICE
          c % fill = lid
          if (mid == MATERIAL_VOID) then
            c % material = mid
          else if (material_dict % has_key(mid)) then
            c % material = material_dict % get_key(mid)
          else
            message = "Could not find material " // trim(to_str(mid)) // &
               " specified on lattice " // trim(to_str(lid))
            call fatal_error()
          end if
          
        else
          message = "Specified fill " // trim(to_str(id)) // " on cell " // &
               trim(to_str(c % id)) // " is neither a universe nor a lattice."
          call fatal_error()
        end if
      end if
    end do

    ! ==========================================================================
    ! ADJUST UNIVERSE INDICES FOR EACH LATTICE

    do i = 1, n_lattices
      lat => lattices(i)
      n_x = lat % dimension(1)
      n_y = lat % dimension(2)
      if (lat % n_dimension == 3) then
        n_z = lat % dimension(3)
      else
        n_z = 1
      end if

      do m = 1, n_z
        do k = 1, n_y
          do j = 1, n_x
            id = lat % universes(j,k,m)
            if (universe_dict % has_key(id)) then
              lat % universes(j,k,m) = universe_dict % get_key(id)
            else
              message = "Invalid universe number " // trim(to_str(id)) &
                   // " specified on lattice " // trim(to_str(lat % id))
              call fatal_error()
            end if
          end do
        end do
      end do

    end do

    TALLY_LOOP: do i = 1, n_tallies
      t => tallies(i)

      ! =======================================================================
      ! ADJUST INDICES FOR EACH TALLY FILTER

      FILTER_LOOP: do j = 1, t % n_filters

        select case (t % filters(j) % type)
        case (FILTER_CELL, FILTER_CELLBORN)

          do k = 1, t % filters(j) % n_bins
            id = t % filters(j) % int_bins(k)
            if (cell_dict % has_key(id)) then
              t % filters(j) % int_bins(k) = cell_dict % get_key(id)
            else
              message = "Could not find cell " // trim(to_str(id)) // &
                   " specified on tally " // trim(to_str(t % id))
              call fatal_error()
            end if
          end do

        case (FILTER_SURFACE)

          ! Check if this is a surface filter only for surface currents
          if (any(t % score_bins == SCORE_CURRENT)) cycle FILTER_LOOP

          do k = 1, t % filters(j) % n_bins
            id = t % filters(j) % int_bins(k)
            if (surface_dict % has_key(id)) then
              t % filters(j) % int_bins(k) = surface_dict % get_key(id)
            else
              message = "Could not find surface " // trim(to_str(id)) // &
                   " specified on tally " // trim(to_str(t % id))
              call fatal_error()
            end if
          end do

        case (FILTER_UNIVERSE)

          do k = 1, t % filters(j) % n_bins
            id = t % filters(j) % int_bins(k)
            if (universe_dict % has_key(id)) then
              t % filters(j) % int_bins(k) = universe_dict % get_key(id)
            else
              message = "Could not find universe " // trim(to_str(id)) // &
                   " specified on tally " // trim(to_str(t % id))
              call fatal_error()
            end if
          end do

        case (FILTER_MATERIAL)

          do k = 1, t % filters(j) % n_bins
            id = t % filters(j) % int_bins(k)
            if (material_dict % has_key(id)) then
              t % filters(j) % int_bins(k) = material_dict % get_key(id)
            else
              message = "Could not find material " // trim(to_str(id)) // &
                   " specified on tally " // trim(to_str(t % id))
              call fatal_error()
            end if
          end do

        case (FILTER_MESH)

          ! The mesh filter already has been set to the index in meshes rather
          ! than the user-specified id, so it doesn't need to be changed.

        end select

      end do FILTER_LOOP

    end do TALLY_LOOP

  end subroutine adjust_indices

!===============================================================================
! NORMALIZE_AO normalizes the atom or weight percentages for each material
!===============================================================================

  subroutine normalize_ao()

    integer        :: index_list      ! index in xs_listings array
    integer        :: i               ! index in materials array
    integer        :: j               ! index over nuclides in material
    real(8)        :: sum_percent     ! summation
    real(8)        :: awr             ! atomic weight ratio
    real(8)        :: x               ! atom percent
    logical        :: percent_in_atom ! nuclides specified in atom percent?
    logical        :: density_in_atom ! density specified in atom/b-cm?
    type(Material), pointer :: mat => null()

    ! first find the index in the xs_listings array for each nuclide in each
    ! material
    do i = 1, n_materials
      mat => materials(i)

      percent_in_atom = (mat % atom_density(1) > ZERO)
      density_in_atom = (mat % density > ZERO)

      sum_percent = ZERO
      do j = 1, mat % n_nuclides
        ! determine atomic weight ratio
        index_list = xs_listing_dict % get_key(mat % names(j))
        awr = xs_listings(index_list) % awr

        ! if given weight percent, convert all values so that they are divided
        ! by awr. thus, when a sum is done over the values, it's actually
        ! sum(w/awr)
        if (.not. percent_in_atom) then
          mat % atom_density(j) = -mat % atom_density(j) / awr
        end if
      end do

      ! determine normalized atom percents. if given atom percents, this is
      ! straightforward. if given weight percents, the value is w/awr and is
      ! divided by sum(w/awr)
      sum_percent = sum(mat % atom_density)
      mat % atom_density = mat % atom_density / sum_percent

      ! Change density in g/cm^3 to atom/b-cm. Since all values are now in atom
      ! percent, the sum needs to be re-evaluated as 1/sum(x*awr)
      if (.not. density_in_atom) then
        sum_percent = ZERO
        do j = 1, mat % n_nuclides
          index_list = xs_listing_dict % get_key(mat % names(j))
          awr = xs_listings(index_list) % awr
          x = mat % atom_density(j)
          sum_percent = sum_percent + x*awr
        end do
        sum_percent = ONE / sum_percent
        mat % density = -mat % density * N_AVOGADRO & 
             / MASS_NEUTRON * sum_percent
      end if

      ! Calculate nuclide atom densities
      mat % atom_density = mat % density * mat % atom_density
    end do

  end subroutine normalize_ao

!===============================================================================
! CALCULATE_WORK determines how many particles each processor should simulate
!===============================================================================

  subroutine calculate_work()

    integer    :: i         ! loop index
    integer    :: remainder ! Number of processors with one extra particle
    integer(8) :: i_bank    ! Running count of number of particles
    integer(8) :: min_work  ! Minimum number of particles on each proc
    integer(8) :: work_i    ! Number of particles on rank i

    allocate(work_index(0:n_procs))

    ! Determine minimum amount of particles to simulate on each processor
    min_work = n_particles/n_procs

    ! Determine number of processors that have one extra particle
    remainder = int(mod(n_particles, int(n_procs,8)), 4)

    i_bank = 0
    work_index(0) = 0
    do i = 0, n_procs - 1
      ! Number of particles for rank i
      if (i < remainder) then
        work_i = min_work + 1
      else
        work_i = min_work
      end if

      ! Set number of particles
      if (rank == i) work = work_i

      ! Set index into source bank for rank i
      i_bank = i_bank + work_i
      work_index(i+1) = i_bank
    end do

  end subroutine calculate_work

!===============================================================================
! ALLOCATE_BANKS allocates memory for the fission and source banks
!===============================================================================

  subroutine allocate_banks()

    integer :: alloc_err  ! allocation error code

    ! Allocate source bank
    allocate(source_bank(work), STAT=alloc_err)

    ! Check for allocation errors 
    if (alloc_err /= 0) then
      message = "Failed to allocate source bank."
      call fatal_error()
    end if

#ifdef _OPENMP
    ! If OpenMP is being used, each thread needs its own private fission
    ! bank. Since the private fission banks need to be combined at the end of a
    ! generation, there is also a 'master_fission_bank' that is used to collect
    ! the sites from each thread.

    n_threads = omp_get_max_threads()

!$omp parallel
    thread_id = omp_get_thread_num()

    if (thread_id == 0) then
       allocate(fission_bank(3*work))
    else
       allocate(fission_bank(3*work/n_threads))
    end if
!$omp end parallel
    allocate(master_fission_bank(3*work), STAT=alloc_err)
#else
    allocate(fission_bank(3*work), STAT=alloc_err)
#endif

    ! Check for allocation errors 
    if (alloc_err /= 0) then
      message = "Failed to allocate fission bank."
      call fatal_error()
    end if

  end subroutine allocate_banks

end module initialize
