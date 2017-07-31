module initialize

  use, intrinsic :: ISO_C_BINDING, only: c_loc

  use hdf5
#ifdef _OPENMP
  use omp_lib
#endif

  use bank_header,     only: Bank
  use constants
  use dict_header,     only: DictIntInt, ElemKeyValueII
  use set_header,      only: SetInt
  use error,           only: fatal_error, warning
  use geometry,        only: neighbor_lists, count_instance, calc_offsets, &
                             maximum_levels
  use geometry_header, only: Cell, Universe, Lattice, RectLattice, HexLattice,&
                             root_universe
  use global
  use hdf5_interface,  only: file_open, read_attribute, file_close, &
                             hdf5_bank_t, hdf5_integer8_t
  use input_xml,       only: read_input_xml, read_plots_xml
  use material_header, only: Material
  use message_passing
  use mgxs_data,       only: read_mgxs, create_macro_xs
  use output,          only: print_version, write_message, print_usage, &
                             print_plot
  use random_lcg,      only: initialize_prng
  use state_point,     only: load_state_point
  use string,          only: to_str, starts_with, ends_with, str_to_int
  use summary,         only: write_summary
  use tally_header,    only: TallyObject
  use tally_initialize,only: configure_tallies
  use tally_filter
  use tally,           only: init_tally_routines

  implicit none

contains

!===============================================================================
! OPENMC_INIT takes care of all initialization tasks, i.e. reading
! from command line, reading xml input files, initializing random
! number seeds, reading cross sections, initializing starting source,
! setting up timers, etc.
!===============================================================================

  subroutine openmc_init(intracomm) bind(C)
    integer, intent(in), optional :: intracomm  ! MPI intracommunicator

    ! Copy the communicator to a new variable. This is done to avoid changing
    ! the signature of this subroutine. If MPI is being used but no communicator
    ! was passed, assume MPI_COMM_WORLD.
#ifdef MPI
#ifdef MPIF08
    type(MPI_Comm), intent(in) :: comm     ! MPI intracommunicator
    if (present(intracomm)) then
      comm % MPI_VAL = intracomm
    else
      comm = MPI_COMM_WORLD
    end if
#else
    integer :: comm
    if (present(intracomm)) then
      comm = intracomm
    else
      comm = MPI_COMM_WORLD
    end if
#endif
#endif

    ! Start total and initialization timer
    call time_total%start()
    call time_initialize%start()

#ifdef MPI
    ! Setup MPI
    call initialize_mpi(comm)
#endif

    ! Initialize HDF5 interface
    call hdf5_initialize()

    ! Read command line arguments
    call read_command_line()

    ! Read XML input files
    call read_input_xml()

    ! Initialize random number generator -- this has to be done after the input
    ! files have been read in case the user specified a seed for the random
    ! number generator

    call initialize_prng()

    ! Read plots.xml if it exists -- this has to be done separate from the other
    ! XML files because we need the PRNG to be initialized first
    if (run_mode == MODE_PLOTTING) call read_plots_xml()

    ! Use dictionaries to redefine index pointers
    call adjust_indices()

    ! Initialize distribcell_filters
    call prepare_distribcell()

    ! After reading input and basic geometry setup is complete, build lists of
    ! neighboring cells for efficient tracking
    call neighbor_lists()

    ! Check to make sure there are not too many nested coordinate levels in the
    ! geometry since the coordinate list is statically allocated for performance
    ! reasons
    if (maximum_levels(universes(root_universe)) > MAX_COORD) then
      call fatal_error("Too many nested coordinate levels in the geometry. &
           &Try increasing the maximum number of coordinate levels by &
           &providing the CMake -Dmaxcoord= option.")
    end if

    if (run_mode /= MODE_PLOTTING) then
      ! Allocate and setup tally stride, filter_matches, and tally maps
      call configure_tallies()

      ! Set up tally procedure pointers
      call init_tally_routines()

      ! Determine how much work each processor should do
      call calculate_work()

      ! Allocate source bank, and for eigenvalue simulations also allocate the
      ! fission bank
      call allocate_banks()

      ! If this is a restart run, load the state point data and binary source
      ! file
      if (restart_run) call load_state_point()
    end if

    if (master) then
      if (run_mode == MODE_PLOTTING) then
        ! Display plotting information
        if (verbosity >= 5) call print_plot()
      else
        ! Write summary information
        if (output_summary) call write_summary()
      end if
    end if

    ! Check for particle restart run
    if (particle_restart_run) run_mode = MODE_PARTICLE

    ! Warn if overlap checking is on
    if (master .and. check_overlaps .and. run_mode /= MODE_PLOTTING) then
      call warning("Cell overlap checking is ON.")
    end if

    ! Stop initialization timer
    call time_initialize%stop()

  end subroutine openmc_init

#ifdef MPI
!===============================================================================
! INITIALIZE_MPI starts up the Message Passing Interface (MPI) and determines
! the number of processors the problem is being run with as well as the rank of
! each processor.
!===============================================================================

  subroutine initialize_mpi(intracomm)
#ifdef MPIF08
    type(MPI_Comm), intent(in) :: intracomm  ! MPI intracommunicator
#else
    integer, intent(in) :: intracomm         ! MPI intracommunicator
#endif

    integer                   :: mpi_err          ! MPI error code
    integer                   :: bank_blocks(5)   ! Count for each datatype
#ifdef MPIF08
    type(MPI_Datatype)        :: bank_types(5)
#else
    integer                   :: bank_types(5)    ! Datatypes
#endif
    integer(MPI_ADDRESS_KIND) :: bank_disp(5)     ! Displacements
    logical    :: init_called
    type(Bank) :: b

    ! Indicate that MPI is turned on
    mpi_enabled = .true.

    ! Initialize MPI
    call MPI_INITIALIZED(init_called, mpi_err)
    if (.not. init_called) call MPI_INIT(mpi_err)

    ! Determine number of processors and rank of each processor
    mpi_intracomm = intracomm
    call MPI_COMM_SIZE(mpi_intracomm, n_procs, mpi_err)
    call MPI_COMM_RANK(mpi_intracomm, rank, mpi_err)

    ! Determine master
    if (rank == 0) then
      master = .true.
    else
      master = .false.
    end if

    ! ==========================================================================
    ! CREATE MPI_BANK TYPE

    ! Determine displacements for MPI_BANK type
    call MPI_GET_ADDRESS(b % wgt,           bank_disp(1), mpi_err)
    call MPI_GET_ADDRESS(b % xyz,           bank_disp(2), mpi_err)
    call MPI_GET_ADDRESS(b % uvw,           bank_disp(3), mpi_err)
    call MPI_GET_ADDRESS(b % E,             bank_disp(4), mpi_err)
    call MPI_GET_ADDRESS(b % delayed_group, bank_disp(5), mpi_err)

    ! Adjust displacements
    bank_disp = bank_disp - bank_disp(1)

    ! Define MPI_BANK for fission sites
    bank_blocks = (/ 1, 3, 3, 1, 1 /)
    bank_types = (/ MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_INTEGER /)
    call MPI_TYPE_CREATE_STRUCT(5, bank_blocks, bank_disp, &
         bank_types, MPI_BANK, mpi_err)
    call MPI_TYPE_COMMIT(MPI_BANK, mpi_err)

  end subroutine initialize_mpi
#endif

!===============================================================================
! HDF5_INITIALIZE
!===============================================================================

  subroutine hdf5_initialize()

    type(Bank),        target :: tmpb(2)         ! temporary Bank
    integer                   :: hdf5_err
    integer(HID_T)            :: coordinates_t   ! HDF5 type for 3 reals
    integer(HSIZE_T)          :: dims(1) = (/3/) ! size of coordinates

    ! Initialize FORTRAN interface.
    call h5open_f(hdf5_err)

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
    call h5tinsert_f(hdf5_bank_t, "delayed_group", h5offsetof(c_loc(tmpb(1)), &
         c_loc(tmpb(1)%delayed_group)), H5T_NATIVE_INTEGER, hdf5_err)

    ! Determine type for integer(8)
    hdf5_integer8_t = h5kind_to_type(8, H5_INTEGER_KIND)

  end subroutine hdf5_initialize

!===============================================================================
! READ_COMMAND_LINE reads all parameters from the command line
!===============================================================================

  subroutine read_command_line()

    integer :: i         ! loop index
    integer :: argc      ! number of command line arguments
    integer :: last_flag ! index of last flag
    character(MAX_WORD_LEN) :: filetype
    integer(HID_T) :: file_id
    character(MAX_WORD_LEN), allocatable :: argv(:) ! command line arguments

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
            call fatal_error("Must specify integer after " // trim(argv(i-1)) &
                 &// " command-line flag.")
          end if
        case ('-r', '-restart', '--restart')
          ! Read path for state point/particle restart
          i = i + 1

          ! Check what type of file this is
          file_id = file_open(argv(i), 'r', parallel=.true.)
          call read_attribute(filetype, file_id, 'filetype')
          call file_close(file_id)

          ! Set path and flag for type of run
          select case (trim(filetype))
          case ('statepoint')
            path_state_point = argv(i)
            restart_run = .true.
          case ('particle restart')
            path_particle_restart = argv(i)
            particle_restart_run = .true.
          case default
            call fatal_error("Unrecognized file after restart flag: " // filetype // ".")
          end select

          ! If its a restart run check for additional source file
          if (restart_run .and. i + 1 <= argc) then

            ! Increment arg
            i = i + 1

            ! Check if it has extension we can read
            if (ends_with(argv(i), '.h5')) then

              ! Check file type is a source file
              file_id = file_open(argv(i), 'r', parallel=.true.)
              call read_attribute(filetype, file_id, 'filetype')
              call file_close(file_id)
              if (filetype /= 'source') then
                call fatal_error("Second file after restart flag must be a &
                     &source file")
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

        case ('-c', '--volume')
          run_mode = MODE_VOLUME

        case ('-s', '--threads')
          ! Read number of threads
          i = i + 1

#ifdef _OPENMP
          ! Read and set number of OpenMP threads
          n_threads = int(str_to_int(argv(i)), 4)
          if (n_threads < 1) then
            call fatal_error("Invalid number of threads specified on command &
                 &line.")
          end if
          call omp_set_num_threads(n_threads)
#else
          if (master) call warning("Ignoring number of threads specified on &
               &command line.")
#endif

        case ('-?', '-h', '-help', '--help')
          call print_usage()
          stop
        case ('-v', '-version', '--version')
          call print_version()
          stop
        case ('-t', '-track', '--track')
          write_all_tracks = .true.
        case default
          call fatal_error("Unknown command line option: " // argv(i))
        end select

        last_flag = i
      end if

      ! Increment counter
      i = i + 1
    end do

    ! Determine directory where XML input files are
    if (argc > 0 .and. last_flag < argc) then
      path_input = argv(last_flag + 1)
    else
      path_input = ''
    end if

    ! Add slash at end of directory if it isn't there
    if (.not. ends_with(path_input, "/") .and. len_trim(path_input) > 0) then
      path_input = trim(path_input) // "/"
    end if

    ! Free memory from argv
    deallocate(argv)

    ! TODO: Check that directory exists

  end subroutine read_command_line

!===============================================================================
! ADJUST_INDICES changes the values for 'surfaces' for each cell and the
! material index assigned to each to the indices in the surfaces and material
! array rather than the unique IDs assigned to each surface and material. Also
! assigns boundary conditions to surfaces based on those read into the bc_dict
! dictionary
!===============================================================================

  subroutine adjust_indices()

    integer :: i                      ! index for various purposes
    integer :: j                      ! index for various purposes
    integer :: k                      ! loop index for lattices
    integer :: m                      ! loop index for lattices
    integer :: lid                    ! lattice IDs
    integer :: i_array                ! index in surfaces/materials array
    integer :: id                     ! user-specified id
    type(Cell),        pointer :: c => null()
    class(Lattice),    pointer :: lat => null()

    do i = 1, n_cells
      ! =======================================================================
      ! ADJUST REGION SPECIFICATION FOR EACH CELL

      c => cells(i)
      do j = 1, size(c%region)
        id = c%region(j)
        ! Make sure that only regions are checked. Since OP_UNION is the
        ! operator with the lowest integer value, anything below it must denote
        ! a half-space
        if (id < OP_UNION) then
          if (surface_dict%has_key(abs(id))) then
            i_array = surface_dict%get_key(abs(id))
            c%region(j) = sign(i_array, id)
          else
            call fatal_error("Could not find surface " // trim(to_str(abs(id)))&
                 &// " specified on cell " // trim(to_str(c%id)))
          end if
        end if
      end do

      ! Also adjust the indices in the reverse Polish notation
      do j = 1, size(c%rpn)
        id = c%rpn(j)
        ! Again, make sure that only regions are checked
        if (id < OP_UNION) then
          i_array = surface_dict%get_key(abs(id))
          c%rpn(j) = sign(i_array, id)
        end if
      end do

      ! =======================================================================
      ! ADJUST UNIVERSE INDEX FOR EACH CELL

      id = c%universe
      if (universe_dict%has_key(id)) then
        c%universe = universe_dict%get_key(id)
      else
        call fatal_error("Could not find universe " // trim(to_str(id)) &
             &// " specified on cell " // trim(to_str(c%id)))
      end if

      ! =======================================================================
      ! ADJUST MATERIAL/FILL POINTERS FOR EACH CELL

      if (c % material(1) == NONE) then
        id = c % fill
        if (universe_dict % has_key(id)) then
          c % type = FILL_UNIVERSE
          c % fill = universe_dict % get_key(id)
        elseif (lattice_dict % has_key(id)) then
          lid = lattice_dict % get_key(id)
          c % type = FILL_LATTICE
          c % fill = lid
        else
          call fatal_error("Specified fill " // trim(to_str(id)) // " on cell "&
               // trim(to_str(c % id)) // " is neither a universe nor a &
               &lattice.")
        end if
      else
        do j = 1, size(c % material)
          id = c % material(j)
          if (id == MATERIAL_VOID) then
            c % type = FILL_MATERIAL
          else if (material_dict % has_key(id)) then
            c % type = FILL_MATERIAL
            c % material(j) = material_dict % get_key(id)
          else
            call fatal_error("Could not find material " // trim(to_str(id)) &
                 // " specified on cell " // trim(to_str(c % id)))
          end if
        end do
      end if
    end do

    ! ==========================================================================
    ! ADJUST UNIVERSE INDICES FOR EACH LATTICE

    do i = 1, n_lattices
      lat => lattices(i)%obj
      select type (lat)

      type is (RectLattice)
        do m = 1, lat%n_cells(3)
          do k = 1, lat%n_cells(2)
            do j = 1, lat%n_cells(1)
              id = lat%universes(j,k,m)
              if (universe_dict%has_key(id)) then
                lat%universes(j,k,m) = universe_dict%get_key(id)
              else
                call fatal_error("Invalid universe number " &
                     &// trim(to_str(id)) // " specified on lattice " &
                     &// trim(to_str(lat%id)))
              end if
            end do
          end do
        end do

      type is (HexLattice)
        do m = 1, lat%n_axial
          do k = 1, 2*lat%n_rings - 1
            do j = 1, 2*lat%n_rings - 1
              if (j + k < lat%n_rings + 1) then
                cycle
              else if (j + k > 3*lat%n_rings - 1) then
                cycle
              end if
              id = lat%universes(j, k, m)
              if (universe_dict%has_key(id)) then
                lat%universes(j, k, m) = universe_dict%get_key(id)
              else
                call fatal_error("Invalid universe number " &
                     &// trim(to_str(id)) // " specified on lattice " &
                     &// trim(to_str(lat%id)))
              end if
            end do
          end do
        end do

      end select

      if (lat%outer /= NO_OUTER_UNIVERSE) then
        if (universe_dict%has_key(lat%outer)) then
          lat%outer = universe_dict%get_key(lat%outer)
        else
          call fatal_error("Invalid universe number " &
               &// trim(to_str(lat%outer)) &
               &// " specified on lattice " // trim(to_str(lat%id)))
        end if
      end if

    end do

    ! =======================================================================
    ! ADJUST INDICES FOR EACH TALLY FILTER

    FILTER_LOOP: do i = 1, n_filters

      select type(filt => filters(i) % obj)
      type is (SurfaceFilter)
        ! Check if this is a surface filter only for surface currents
        if (.not. filt % current) call filt % initialize()
      class default
        call filt % initialize()
      end select

    end do FILTER_LOOP

  end subroutine adjust_indices

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
      call fatal_error("Failed to allocate source bank.")
    end if

    if (run_mode == MODE_EIGENVALUE) then
#ifdef _OPENMP
      ! If OpenMP is being used, each thread needs its own private fission
      ! bank. Since the private fission banks need to be combined at the end of
      ! a generation, there is also a 'master_fission_bank' that is used to
      ! collect the sites from each thread.

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
        call fatal_error("Failed to allocate fission bank.")
      end if
    end if

  end subroutine allocate_banks

!===============================================================================
! PREPARE_DISTRIBCELL initializes any distribcell filters present and sets the
! offsets for distribcells
!===============================================================================

  subroutine prepare_distribcell()

    integer :: i, j                ! Tally, filter loop counters
    logical :: distribcell_active  ! Does simulation use distribcell?
    integer, allocatable :: univ_list(:)              ! Target offsets
    integer, allocatable :: counts(:,:)               ! Target count
    logical, allocatable :: found(:,:)                ! Target found

    ! Assume distribcell is not needed until proven otherwise.
    distribcell_active = .false.

    ! We need distribcell if any tallies have distribcell filters.
    do i = 1, n_tallies
      do j = 1, size(tallies(i) % filter)
        select type(filt => filters(tallies(i) % filter(j)) % obj)
        type is (DistribcellFilter)
          distribcell_active = .true.
        end select
      end do
    end do

    ! We also need distribcell if any distributed materials or distributed
    ! temperatues are present.
    if (.not. distribcell_active) then
      do i = 1, n_cells
        if (size(cells(i) % material) > 1 .or. size(cells(i) % sqrtkT) > 1) then
          distribcell_active = .true.
          exit
        end if
      end do
    end if

    ! If distribcell isn't used in this simulation then no more work left to do.
    if (.not. distribcell_active) return

    ! Count the number of instances of each cell.
    call count_instance(universes(root_universe))

    ! Set the number of bins in all distribcell filters.
    do i = 1, n_tallies
      do j = 1, size(tallies(i) % filter)
        select type(filt => filters(tallies(i) % filter(j)) % obj)
        type is (DistribcellFilter)
          ! Set the number of bins to the number of instances of the cell.
          filt % n_bins = cells(filt % cell) % instances
        end select
      end do
    end do

    ! Make sure the number of materials and temperatures matches the number of
    ! cell instances.
    do i = 1, n_cells
      associate (c => cells(i))
        if (size(c % material) > 1) then
          if (size(c % material) /= c % instances) then
            call fatal_error("Cell " // trim(to_str(c % id)) // " was &
                 &specified with " // trim(to_str(size(c % material))) &
                 // " materials but has " // trim(to_str(c % instances)) &
                 // " distributed instances. The number of materials must &
                 &equal one or the number of instances.")
          end if
        end if
        if (size(c % sqrtkT) > 1) then
          if (size(c % sqrtkT) /= c % instances) then
            call fatal_error("Cell " // trim(to_str(c % id)) // " was &
                 &specified with " // trim(to_str(size(c % sqrtkT))) &
                 // " temperatures but has " // trim(to_str(c % instances)) &
                 // " distributed instances. The number of temperatures must &
                 &equal one or the number of instances.")
          end if
        end if
      end associate
    end do

    ! Allocate offset maps at each level in the geometry
    call allocate_offsets(univ_list, counts, found)

    ! Calculate offsets for each target distribcell
    do i = 1, n_maps
      do j = 1, n_universes
        call calc_offsets(univ_list(i), i, universes(j), counts, found)
      end do
    end do

  end subroutine prepare_distribcell

!===============================================================================
! ALLOCATE_OFFSETS determines the number of maps needed and allocates required
! memory for distribcell offset tables
!===============================================================================

  recursive subroutine allocate_offsets(univ_list, counts, found)

    integer, intent(out), allocatable     :: univ_list(:) ! Target offsets
    integer, intent(out), allocatable     :: counts(:,:)  ! Target count
    logical, intent(out), allocatable     :: found(:,:)   ! Target found

    integer      :: i, j, k   ! Loop counters
    type(SetInt) :: cell_list ! distribells to track

    ! Begin gathering list of cells in distribcell tallies
    n_maps = 0

    ! List all cells referenced in distribcell filters.
    do i = 1, n_tallies
      do j = 1, size(tallies(i) % filter)
        select type(filt => filters(tallies(i) % filter(j)) % obj)
        type is (DistribcellFilter)
          call cell_list % add(filt % cell)
        end select
      end do
    end do

    ! List all cells with multiple (distributed) materials or temperatures.
    do i = 1, n_cells
      if (size(cells(i) % material) > 1 .or. size(cells(i) % sqrtkT) > 1) then
        call cell_list % add(i)
      end if
    end do

    ! Compute the number of unique universes containing these distribcells
    ! to determine the number of offset tables to allocate
    do i = 1, n_universes
      do j = 1, size(universes(i) % cells)
        if (cell_list % contains(universes(i) % cells(j))) then
          n_maps = n_maps + 1
        end if
      end do
    end do

    ! Allocate the list of offset tables for each unique universe
    allocate(univ_list(n_maps))

    ! Allocate list to accumulate target distribcell counts in each universe
    allocate(counts(n_universes, n_maps))
    counts(:,:) = 0

    ! Allocate list to track if target distribcells are found in each universe
    allocate(found(n_universes, n_maps))
    found(:,:) = .false.


    ! Search through universes for distributed cells and assign each one a
    ! unique distribcell array index.
    k = 1
    do i = 1, n_universes
      do j = 1, size(universes(i) % cells)
        if (cell_list % contains(universes(i) % cells(j))) then
          cells(universes(i) % cells(j)) % distribcell_index = k
          univ_list(k) = universes(i) % id
          k = k + 1
        end if
      end do
    end do

    ! Allocate the offset tables for lattices
    do i = 1, n_lattices
      associate(lat => lattices(i) % obj)
        select type(lat)

        type is (RectLattice)
          allocate(lat % offset(n_maps, lat % n_cells(1), lat % n_cells(2), &
                   lat % n_cells(3)))
        type is (HexLattice)
          allocate(lat % offset(n_maps, 2 * lat % n_rings - 1, &
               2 * lat % n_rings - 1, lat % n_axial))
        end select

        lat % offset(:, :, :, :) = 0
      end associate
    end do

    ! Allocate offset table for fill cells
    do i = 1, n_cells
      if (cells(i) % type /= FILL_MATERIAL) then
        allocate(cells(i) % offset(n_maps))
      end if
    end do

    ! Free up memory
    call cell_list % clear()

  end subroutine allocate_offsets

end module initialize
