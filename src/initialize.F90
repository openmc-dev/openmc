module initialize

  use ace,              only: read_xs
  use bank_header,      only: Bank
  use constants
  use datatypes,        only: dict_create, dict_get_key, dict_has_key,         &
                              dict_keys
  use datatypes_header, only: ListKeyValueII, DictionaryII
  use energy_grid,      only: unionized_grid
  use error,            only: fatal_error
  use geometry,         only: neighbor_lists
  use geometry_header,  only: Cell, Universe, Lattice, BASE_UNIVERSE
  use global
  use input_xml,        only: read_input_xml, read_cross_sections_xml,         &
                              cells_in_univ_dict, read_plots_xml
  use output,           only: title, header, print_summary, print_geometry,    &
                              print_plot, create_summary_file, print_usage,    &
                              create_xs_summary_file, print_version
  use random_lcg,       only: initialize_prng
  use source,           only: allocate_banks, initialize_source
  use string,           only: to_str, str_to_int, starts_with, ends_with,      &
                              lower_case
  use tally,            only: create_tally_map
  use tally_header,     only: TallyObject
  use timing,           only: timer_start, timer_stop

#ifdef MPI
  use mpi
#endif

#ifdef HDF5
  use hdf5_interface,   only: hdf5_create_output, hdf5_write_header,            &
                              hdf5_write_summary
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
    call timer_start(time_total)
    call timer_start(time_initialize)

    ! Setup MPI
    call setup_mpi()

    ! Read command line arguments
    call read_command_line()

    if (master) then
       ! Create output files
       call create_summary_file()
       call create_xs_summary_file()

#ifdef HDF5
       ! Open HDF5 output file for writing and write header information
       call hdf5_create_output()
       call hdf5_write_header()
#endif

       ! Display title and initialization header
       call title()
       call header("INITIALIZATION", level=1)
    end if

    ! set up dictionaries
    call create_dictionaries()

    ! Read XML input files
    call read_input_xml()

    ! Initialize random number generator
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
       ! Read cross section summary file to determine what files contain
       ! cross-sections
       call read_cross_sections_xml()

       ! With the AWRs from the xs_listings, change all material specifications
       ! so that they contain atom percents summing to 1
       call normalize_ao()

       ! Read ACE-format cross sections
       call timer_start(time_read_xs)
       call read_xs()
       call timer_stop(time_read_xs)

       ! Construct unionized energy grid from cross-sections
       if (grid_method == GRID_UNION) then
          call timer_start(time_unionize)
          call unionized_grid()
          call timer_stop(time_unionize)
       end if

       ! Create tally map
       call create_tally_map()

       ! allocate banks and create source particles
       call allocate_banks()
       call initialize_source()
    end if

    ! stop timer for initialization
    if (master) then
       if (run_mode == MODE_PLOTTING) then
          call print_geometry()
          call print_plot()
       else
          call print_summary()
#ifdef HDF5
          call hdf5_write_summary()
#endif
       end if
    end if

    ! Stop initialization timer
    call timer_stop(time_initialize)

  end subroutine initialize_run

!===============================================================================
! SETUP_MPI initilizes the Message Passing Interface (MPI) and determines the
! number of processors the problem is being run with as well as the rank of each
! processor.
!===============================================================================

  subroutine setup_mpi()

#ifdef MPI
    integer        :: bank_blocks(5) ! Count for each datatype
    integer        :: bank_types(5)  ! Datatypes
    integer(MPI_ADDRESS_KIND) :: bank_disp(5)   ! Displacements
    type(Bank)     :: b

    mpi_enabled = .true.

    ! Initialize MPI
    call MPI_INIT(mpi_err)
    if (mpi_err /= MPI_SUCCESS) then
       message = "Failed to initialize MPI."
       call fatal_error()
    end if

    ! Determine number of processors
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, mpi_err)
    if (mpi_err /= MPI_SUCCESS) then
       message = "Could not determine number of processors."
       call fatal_error()
    end if

    ! Determine rank of each processor
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    if (mpi_err /= MPI_SUCCESS) then
       message = "Could not determine MPI rank."
       call fatal_error()
    end if

    ! Determine master
    if (rank == 0) then
       master = .true.
    else
       master = .false.
    end if

    ! Determine displacements for MPI_BANK type
    call MPI_GET_ADDRESS(b % id,  bank_disp(1), mpi_err)
    call MPI_GET_ADDRESS(b % wgt, bank_disp(2), mpi_err)
    call MPI_GET_ADDRESS(b % xyz, bank_disp(3), mpi_err)
    call MPI_GET_ADDRESS(b % uvw, bank_disp(4), mpi_err)
    call MPI_GET_ADDRESS(b % E,   bank_disp(5), mpi_err)

    ! Adjust displacements 
    bank_disp = bank_disp - bank_disp(1)
    
    ! Define MPI_BANK for fission sites
    bank_blocks = (/ 1, 1, 3, 3, 1 /)
    bank_types = (/ MPI_INTEGER8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8 /)
    call MPI_TYPE_CREATE_STRUCT(5, bank_blocks, bank_disp, & 
         bank_types, MPI_BANK, mpi_err)
    call MPI_TYPE_COMMIT(MPI_BANK, mpi_err)

#else
    ! if no MPI, set processor to master
    mpi_enabled = .false.
    rank = 0
    n_procs = 1
    master = .true.
#endif

  end subroutine setup_mpi

!===============================================================================
! READ_COMMAND_LINE reads all parameters from the command line
!===============================================================================

  subroutine read_command_line()

    integer :: i         ! loop index
    integer :: argc      ! number of command line arguments
    integer :: last_flag ! index of last flag
    character(MAX_FILE_LEN) :: pwd      ! present working directory
    character(MAX_WORD_LEN), allocatable :: argv(:) ! command line arguments
    
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
          case ('-n', '-n_particles', '--n_particles')
             i = i + 1
             ! Read number of particles per cycle
             n_particles = str_to_int(argv(i))

             ! Check that number specified was valid
             if (n_particles == ERROR_INT) then
                message = "Must specify integer after " // trim(argv(i-1)) // &
                     " command-line flag."
                call fatal_error()
             end if
          case ('-?', '-help', '--help')
             call print_usage()
             stop
          case ('-v', '-version', '--version')
             call print_version()
             stop
          case ('-eps_tol', '-ksp_gmres_restart')
             ! Handle options that would be based to PETSC
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
! CREATE_DICTIONARIES initializes the various dictionary variables. It would be
! nice to avoid this and just have a check at the beginning of
! dictionary_add_key.
!===============================================================================

  subroutine create_dictionaries()

    ! Create all global dictionaries
    call dict_create(cell_dict)
    call dict_create(universe_dict)
    call dict_create(lattice_dict)
    call dict_create(surface_dict)
    call dict_create(material_dict)
    call dict_create(mesh_dict)
    call dict_create(tally_dict)
    call dict_create(xs_listing_dict)
    call dict_create(nuclide_dict)
    call dict_create(sab_dict)

    ! Create special dictionary used in input_xml
    call dict_create(cells_in_univ_dict)
    
  end subroutine create_dictionaries

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
    type(ListKeyValueII), pointer :: key_list => null()
    type(Universe),       pointer :: univ => null()
    type(Cell),           pointer :: c => null()

    allocate(universes(n_universes))

    ! We also need to allocate the cell count lists for each universe. The logic
    ! for this is a little more convoluted. In universe_dict, the (key,value)
    ! pairs are the id of the universe and the index in the array. In
    ! cells_in_univ_dict, it's the id of the universe and the number of cells.

    key_list => dict_keys(universe_dict)
    do while (associated(key_list))
       ! find index of universe in universes array
       i_univ = key_list % data % value
       univ => universes(i_univ)
       univ % id = key_list % data % key

       ! check for lowest level universe
       if (univ % id == 0) BASE_UNIVERSE = i_univ
       
       ! find cell count for this universe
       n_cells_in_univ = dict_get_key(cells_in_univ_dict, univ % id)

       ! allocate cell list for universe
       allocate(univ % cells(n_cells_in_univ))
       univ % n_cells = n_cells_in_univ
       
       ! move to next universe
       key_list => key_list % next
    end do

    ! Also allocate a list for keeping track of where cells have been assigned
    ! in each universe

    allocate(index_cell_in_univ(n_universes))
    index_cell_in_univ = 0

    do i = 1, n_cells
       c => cells(i)

       ! get pointer to corresponding universe
       i_univ = dict_get_key(universe_dict, c % universe)
       univ => universes(i_univ)

       ! increment the index for the cells array within the Universe object and
       ! then store the index of the Cell object in that array
       index_cell_in_univ(i_univ) = index_cell_in_univ(i_univ) + 1
       univ % cells(index_cell_in_univ(i_univ)) = i
    end do

  end subroutine prepare_universes

!===============================================================================
! ADJUST_INDICES changes the values for 'surfaces' for each cell and the
! material index assigned to each to the indices in the surfaces and material
! array rather than the unique IDs assigned to each surface and material. Also
! assigns boundary conditions to surfaces based on those read into the bc_dict
! dictionary
!===============================================================================

  subroutine adjust_indices()

    integer :: i       ! index for various purposes
    integer :: j       ! index for various purposes
    integer :: k       ! loop index for lattices
    integer :: i_array ! index in surfaces/materials array 
    integer :: id      ! user-specified id
    type(Cell),        pointer :: c => null()
    type(Lattice),     pointer :: l => null()
    type(TallyObject), pointer :: t => null()
    
    do i = 1, n_cells
       ! =======================================================================
       ! ADJUST SURFACE LIST FOR EACH CELL

       c => cells(i)
       do j = 1, c % n_surfaces
          id = c % surfaces(j)
          if (id < OP_DIFFERENCE) then
             if (dict_has_key(surface_dict, abs(id))) then
                i_array = dict_get_key(surface_dict, abs(id))
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
       if (dict_has_key(universe_dict, id)) then
          c % universe = dict_get_key(universe_dict, id)
       else
          message = "Could not find universe " // trim(to_str(id)) // &
               " specified on cell " // trim(to_str(c % id))
          call fatal_error()
       end if

       ! =======================================================================
       ! ADJUST MATERIAL/FILL POINTERS FOR EACH CELL

       id = c % material
       if (id /= 0) then
          if (dict_has_key(material_dict, id)) then
             c % type = CELL_NORMAL
             c % material = dict_get_key(material_dict, id)
          else
             message = "Could not find material " // trim(to_str(id)) // &
                   " specified on cell " // trim(to_str(c % id))
             call fatal_error()
          end if
       else
          id = c % fill
          if (dict_has_key(universe_dict, id)) then
             c % type = CELL_FILL
             c % fill = dict_get_key(universe_dict, id)
          elseif (dict_has_key(lattice_dict, id)) then
             c % type = CELL_LATTICE
             c % fill = dict_get_key(lattice_dict, id)
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
       l => lattices(i)
       do j = 1, l % n_x
          do k = 1, l % n_y
             id = l % element(j,k)
             if (dict_has_key(universe_dict, id)) then
                l % element(j,k) = dict_get_key(universe_dict, id)
             else
                message = "Invalid universe number " // trim(to_str(id)) &
                     // " specified on lattice " // trim(to_str(l % id))
                call fatal_error()
             end if
          end do
       end do
    end do

    do i = 1, n_tallies
       t => tallies(i)

       ! =======================================================================
       ! ADJUST CELL INDICES FOR EACH TALLY

       if (t % n_filter_bins(FILTER_CELL) > 0) then
          do j = 1, t % n_filter_bins(FILTER_CELL)
             id = t % cell_bins(j) % scalar
             if (dict_has_key(cell_dict, id)) then
                t % cell_bins(j) % scalar = dict_get_key(cell_dict, id)
             else
                message = "Could not find cell " // trim(to_str(id)) // &
                     " specified on tally " // trim(to_str(t % id))
                call fatal_error()
             end if
          end do
       end if

       ! =======================================================================
       ! ADJUST SURFACE INDICES FOR EACH TALLY

       if (t % n_filter_bins(FILTER_SURFACE) > 0) then
          do j = 1, t % n_filter_bins(FILTER_SURFACE)
             id = t % surface_bins(j) % scalar
             if (dict_has_key(surface_dict, id)) then
                t % surface_bins(j) % scalar = dict_get_key(surface_dict, id)
             else
                message = "Could not find surface " // trim(to_str(id)) // &
                     " specified on tally " // trim(to_str(t % id))
                call fatal_error()
             end if
          end do
       end if

       ! =======================================================================
       ! ADJUST UNIVERSE INDICES FOR EACH TALLY

       if (t % n_filter_bins(FILTER_UNIVERSE) > 0) then
          do j = 1, t % n_filter_bins(FILTER_UNIVERSE)
             id = t % universe_bins(j) % scalar
             if (dict_has_key(universe_dict, id)) then
                t % universe_bins(j) % scalar = dict_get_key(universe_dict, id)
             else
                message = "Could not find universe " // trim(to_str(id)) // &
                     " specified on tally " // trim(to_str(t % id))
                call fatal_error()
             end if
          end do
       end if

       ! =======================================================================
       ! ADJUST MATERIAL INDICES FOR EACH TALLY

       if (t % n_filter_bins(FILTER_MATERIAL) > 0) then
          do j = 1, t % n_filter_bins(FILTER_MATERIAL)
             id = t % material_bins(j) % scalar
             if (dict_has_key(material_dict, id)) then
                t % material_bins(j) % scalar = dict_get_key(material_dict, id)
             else
                message = "Could not find material " // trim(to_str(id)) // &
                     " specified on tally " // trim(to_str(t % id))
                call fatal_error()
             end if
          end do
       end if

       ! =======================================================================
       ! ADJUST CELLBORN INDICES FOR EACH TALLY

       if (t % n_filter_bins(FILTER_CELLBORN) > 0) then
          do j = 1, t % n_filter_bins(FILTER_CELLBORN)
             id = t % cellborn_bins(j) % scalar
             if (dict_has_key(cell_dict, id)) then
                t % cellborn_bins(j) % scalar = dict_get_key(cell_dict, id)
             else
                message = "Could not find material " // trim(to_str(id)) // &
                     " specified on tally " // trim(to_str(t % id))
                call fatal_error()
             end if
          end do
       end if

       ! =======================================================================
       ! ADJUST MESH INDICES FOR EACH TALLY

       if (t % n_filter_bins(FILTER_MESH) > 0) then
          id = t % mesh
          if (dict_has_key(mesh_dict, id)) then
             t % mesh = dict_get_key(mesh_dict, id)
          else
             message = "Could not find mesh " // trim(to_str(id)) // &
                  " specified on tally " // trim(to_str(t % id))
             call fatal_error()
          end if
       end if
    end do

  end subroutine adjust_indices

!===============================================================================
! NORMALIZE_AO normalizes the atom or weight percentages for each material
!===============================================================================

  subroutine normalize_ao()

    integer        :: index_list      ! index in xs_listings array
    integer        :: i               ! index in materials array
    integer        :: j               ! index over nuclides in material
    integer        :: n               ! length of string
    real(8)        :: sum_percent     ! 
    real(8)        :: awr             ! atomic weight ratio
    real(8)        :: x               ! atom percent
    logical        :: percent_in_atom ! nuclides specified in atom percent?
    logical        :: density_in_atom ! density specified in atom/b-cm?
    character(12)  :: key             ! name of nuclide, e.g. 92235.03c
    type(Material), pointer :: mat => null()
    
    ! first find the index in the xs_listings array for each nuclide in each
    ! material
    do i = 1, n_materials
       mat => materials(i)

       ! Check to make sure either all atom percents or all weight percents are
       ! given
       if (.not. (all(mat%atom_percent > ZERO) .or. & 
            all(mat%atom_percent < ZERO))) then
          message = "Cannot mix atom and weight percents in material " // &
               to_str(mat % id)
          call fatal_error()
       end if

       percent_in_atom = (mat%atom_percent(1) > ZERO)
       density_in_atom = (mat%density > ZERO)

       sum_percent = ZERO
       do j = 1, mat % n_nuclides
          ! Set indices for nuclides
          key = mat % names(j)

          ! Check to make sure cross-section is continuous energy neutron table
          n = len_trim(key)
          if (key(n:n) /= 'c') then
             message = "Cross-section table " // trim(key) // & 
                  " is not a continuous-energy neutron table."
             call fatal_error()
          end if

          if (dict_has_key(xs_listing_dict, key)) then
             index_list = dict_get_key(xs_listing_dict, key)
             mat % xs_listing(j) = index_list
          else
             message = "Cannot find cross-section " // trim(key) // &
                  " in specified cross_sections.xml file."
             call fatal_error()
          end if

          ! determine atomic weight ratio
          awr = xs_listings(index_list) % awr

          ! if given weight percent, convert all values so that they are divided
          ! by awr. thus, when a sum is done over the values, it's actually
          ! sum(w/awr)
          if (.not. percent_in_atom) then
             mat % atom_percent(j) = -mat % atom_percent(j) / awr
          end if
       end do

       ! determine normalized atom percents. if given atom percents, this is
       ! straightforward. if given weight percents, the value is w/awr and is
       ! divided by sum(w/awr)
       sum_percent = sum(mat%atom_percent)
       mat % atom_percent = mat % atom_percent / sum_percent

       ! Change density in g/cm^3 to atom/b-cm. Since all values are now in atom
       ! percent, the sum needs to be re-evaluated as 1/sum(x*awr)
       if (.not. density_in_atom) then
          sum_percent = ZERO
          do j = 1, mat % n_nuclides
             index_list = mat % xs_listing(j)
             awr = xs_listings(index_list) % awr
             x = mat % atom_percent(j)
             sum_percent = sum_percent + x*awr
          end do
          sum_percent = ONE / sum_percent
          mat % density = -mat % density * N_AVOGADRO & 
               / MASS_NEUTRON * sum_percent
       end if

       ! Calculate nuclide atom densities and deallocate atom_percent array
       ! since it is no longer needed past this point
       mat % atom_density = mat % density * mat % atom_percent
       deallocate(mat % atom_percent)
    end do

  end subroutine normalize_ao

end module initialize
