module initialize

  use, intrinsic :: ISO_C_BINDING, only: c_loc

  use hdf5
#ifdef _OPENMP
  use omp_lib
#endif

  use bank_header,     only: Bank
  use constants
  use set_header,      only: SetInt
  use error,           only: fatal_error, warning, write_message
  use geometry_header, only: Cell, Universe, Lattice, RectLattice, HexLattice,&
                             root_universe
  use hdf5_interface,  only: file_open, read_attribute, file_close, &
                             hdf5_bank_t, hdf5_integer8_t
  use input_xml,       only: read_input_xml
  use material_header, only: Material
  use message_passing
  use mgxs_data,       only: read_mgxs, create_macro_xs
  use output,          only: print_version, print_usage
  use random_lcg,      only: openmc_set_seed
  use settings
#ifdef _OPENMP
  use simulation_header, only: n_threads
#endif
  use string,          only: to_str, starts_with, ends_with, str_to_int
  use tally_header,    only: TallyObject
  use tally_filter
  use timer_header

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

#ifdef _OPENMP
    character(MAX_WORD_LEN) :: envvar
#endif

    ! Copy the communicator to a new variable. This is done to avoid changing
    ! the signature of this subroutine. If MPI is being used but no communicator
    ! was passed, assume MPI_COMM_WORLD.
#ifdef OPENMC_MPI
#ifdef OPENMC_MPIF08
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

#ifdef OPENMC_MPI
    ! Setup MPI
    call initialize_mpi(comm)
#endif

#ifdef _OPENMP
    ! Change schedule of main parallel-do loop if OMP_SCHEDULE is set
    call get_environment_variable("OMP_SCHEDULE", envvar)
    if (len_trim(envvar) == 0) then
      call omp_set_schedule(omp_sched_static, 0)
    end if
#endif

    ! Initialize HDF5 interface
    call hdf5_initialize()

    ! Read command line arguments
    call read_command_line()

    ! Initialize random number generator -- if the user specifies a seed, it
    ! will be re-initialized later
    call openmc_set_seed(DEFAULT_SEED)

    ! Read XML input files
    call read_input_xml()

    ! Check for particle restart run
    if (particle_restart_run) run_mode = MODE_PARTICLE

    ! Stop initialization timer
    call time_initialize%stop()

  end subroutine openmc_init

#ifdef OPENMC_MPI
!===============================================================================
! INITIALIZE_MPI starts up the Message Passing Interface (MPI) and determines
! the number of processors the problem is being run with as well as the rank of
! each processor.
!===============================================================================

  subroutine initialize_mpi(intracomm)
#ifdef OPENMC_MPIF08
    type(MPI_Comm), intent(in) :: intracomm  ! MPI intracommunicator
#else
    integer, intent(in) :: intracomm         ! MPI intracommunicator
#endif

    integer                   :: mpi_err          ! MPI error code
    integer                   :: bank_blocks(5)   ! Count for each datatype
#ifdef OPENMC_MPIF08
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

end module initialize
