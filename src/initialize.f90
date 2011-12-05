module initialize

  use constants
  use cross_section,    only: read_xs
  use datatypes,        only: dict_create, dict_add_key, dict_get_key,         &
                              dict_has_key, dict_keys
  use datatypes_header, only: ListKeyValueII, DictionaryII
  use energy_grid,      only: unionized_grid, original_indices
  use error,            only: fatal_error
  use geometry,         only: neighbor_lists
  use geometry_header,  only: Cell, Surface, Universe, Lattice, BASE_UNIVERSE
  use global
  use input_xml,        only: read_input_xml, read_cross_sections_xml,         &
                              cells_in_univ_dict
  use logging,          only: create_log
  use mpi_routines,     only: setup_mpi
  use output,           only: title, header, print_summary, print_geometry,    &
                              print_plot
  use random_lcg,       only: initialize_prng
  use source,           only: initialize_source
  use string,           only: int_to_str, starts_with, ends_with, lower_case
  use tally,            only: create_tally_map, TallyObject
  use timing,           only: timer_start, timer_stop

  implicit none

  type(DictionaryII), pointer :: build_dict => null()

contains

!===============================================================================
! INITIALIZE_RUN takes care of all initialization tasks, i.e. reading
! from command line, reading xml input files, initializing random
! number seeds, reading cross sections, initializing starting source,
! setting up timers, etc.
!===============================================================================

  subroutine initialize_run()

    type(Universe), pointer :: univ

    ! Start initialization timer
    call timer_start(time_initialize)

    ! Setup MPI
    call setup_mpi()

    ! Read command line arguments
    call read_command_line()
    if (master) call create_log()

    ! Print the OpenMC title and version/date/time information
    if (master) call title()

    ! Print initialization header block
    if (master) call header("INITIALIZATION", 1)

    ! Initialize random number generator
    call initialize_prng()

    ! set up dictionaries
    call create_dictionaries()

    ! Read XML input files
    call read_input_xml()

    ! Set up universe structures
    call prepare_universes()

    ! Use dictionaries to redefine index pointers
    call adjust_indices()

    ! determine at which level universes are and link cells to parenting cells
    univ => universes(BASE_UNIVERSE)
    call build_universe(univ, 0, 0)

    ! After reading input and basic geometry setup is complete, build lists of
    ! neighboring cells for efficient tracking
    call neighbor_lists()

    if (.not. plotting) then
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
       call timer_start(time_unionize)
       call unionized_grid()
       call original_indices()
       call timer_stop(time_unionize)

       ! Create tally map
       call create_tally_map()

       ! create source particles
       call initialize_source()
    end if

    ! stop timer for initialization
    if (master) then
       if (plotting) then
          call print_geometry()
          call print_plot()
       else
          call print_summary()
       end if
    end if

    ! Stop initialization timer
    call timer_stop(time_initialize)

  end subroutine initialize_run

!===============================================================================
! READ_COMMAND_LINE reads all parameters from the command line
!===============================================================================

  subroutine read_command_line()

    integer :: i         ! loop index
    integer :: argc      ! number of command line arguments
    integer :: last_flag ! index of last flag
    character(MAX_FILE_LEN) :: pwd      ! present working directory
    character(MAX_WORD_LEN) :: argv(10) ! command line arguments
    
    ! Get working directory
    call GET_ENVIRONMENT_VARIABLE("PWD", pwd)

    ! Check number of command line arguments
    argc = COMMAND_ARGUMENT_COUNT()

    ! Get all command line arguments
    last_flag = 0
    do i = 1, argc
       call GET_COMMAND_ARGUMENT(i, argv(i))

       ! Check for flags
       if (starts_with(argv(i), "-")) then
          if (argv(i)(2:5) == 'plot') then
             plotting = .true.
          end if

          last_flag = i
       end if
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

    ! Create special dictionary for building universes
    call dict_create(build_dict)
    
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
    ! ucount_dict, it's the id of the universe and the number of cells.

    key_list => dict_keys(universe_dict)
    do while (associated(key_list))
       ! find index of universe in universes array
       i_univ = key_list%data%value
       univ => universes(i_univ)
       univ % id = key_list%data%key

       ! check for lowest level universe
       if (univ % id == 0) BASE_UNIVERSE = i_univ
       
       ! find cell count for this universe
       n_cells_in_univ = dict_get_key(cells_in_univ_dict, key_list%data%key)

       ! allocate cell list for universe
       allocate(univ % cells(n_cells_in_univ))
       univ % n_cells = n_cells_in_univ
       
       ! move to next universe
       key_list => key_list%next
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

    integer                 :: i            ! index in cells array
    integer                 :: j            ! index over surface list
    integer                 :: k
    integer                 :: i_array      ! index in surfaces/materials array 
    integer                 :: id           ! user-specified id
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
                message = "Could not find surface " // trim(int_to_str(abs(id))) // &
                     & " specified on cell " // trim(int_to_str(c % id))
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
          message = "Could not find universe " // trim(int_to_str(id)) // &
               " specified on cell " // trim(int_to_str(c % id))
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
             message = "Could not find material " // trim(int_to_str(id)) // &
                   " specified on cell " // trim(int_to_str(c % id))
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
             message = "Specified fill " // trim(int_to_str(id)) // " on cell " // &
                  trim(int_to_str(c % id)) // " is neither a universe nor a lattice."
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
                message = "Invalid universe number " // trim(int_to_str(id)) &
                     // " specified on lattice " // trim(int_to_str(l % id))
                call fatal_error()
             end if
          end do
       end do
    end do

    do i = 1, n_tallies
       t => tallies(i)

       ! =======================================================================
       ! ADJUST CELL INDICES FOR EACH TALLY

       if (t % n_bins(T_CELL) > 0) then
          do j = 1, t % n_bins(T_CELL)
             id = t % cell_bins(j) % scalar
             if (dict_has_key(cell_dict, id)) then
                t % cell_bins(j) % scalar = dict_get_key(cell_dict, id)
             else
                message = "Could not find cell " // trim(int_to_str(id)) // &
                     & " specified on tally " // trim(int_to_str(t % id))
                call fatal_error()
             end if
          end do
       end if

       ! =======================================================================
       ! ADJUST SURFACE INDICES FOR EACH TALLY

       if (t % n_bins(T_SURFACE) > 0) then
          do j = 1, t % n_bins(T_SURFACE)
             id = t % surface_bins(j) % scalar
             if (dict_has_key(surface_dict, id)) then
                t % surface_bins(j) % scalar = dict_get_key(surface_dict, id)
             else
                message = "Could not find surface " // trim(int_to_str(id)) // &
                     & " specified on tally " // trim(int_to_str(t % id))
                call fatal_error()
             end if
          end do
       end if

       ! =======================================================================
       ! ADJUST UNIVERSE INDICES FOR EACH TALLY

       if (t % n_bins(T_UNIVERSE) > 0) then
          do j = 1, t % n_bins(T_UNIVERSE)
             id = t % universe_bins(j) % scalar
             if (dict_has_key(universe_dict, id)) then
                t % universe_bins(j) % scalar = dict_get_key(universe_dict, id)
             else
                message = "Could not find universe " // trim(int_to_str(id)) // &
                     & " specified on tally " // trim(int_to_str(t % id))
                call fatal_error()
             end if
          end do
       end if

       ! =======================================================================
       ! ADJUST MATERIAL INDICES FOR EACH TALLY

       if (t % n_bins(T_MATERIAL) > 0) then
          do j = 1, t % n_bins(T_MATERIAL)
             id = t % material_bins(j) % scalar
             if (dict_has_key(material_dict, id)) then
                t % material_bins(j) % scalar = dict_get_key(material_dict, id)
             else
                message = "Could not find material " // trim(int_to_str(id)) // &
                     & " specified on tally " // trim(int_to_str(t % id))
                call fatal_error()
             end if
          end do
       end if

       ! =======================================================================
       ! ADJUST CELLBORN INDICES FOR EACH TALLY

       if (t % n_bins(T_CELLBORN) > 0) then
          do j = 1, t % n_bins(T_CELLBORN)
             id = t % cellborn_bins(j) % scalar
             if (dict_has_key(cell_dict, id)) then
                t % cellborn_bins(j) % scalar = dict_get_key(cell_dict, id)
             else
                message = "Could not find material " // trim(int_to_str(id)) // &
                     & " specified on tally " // trim(int_to_str(t % id))
                call fatal_error()
             end if
          end do
       end if

       ! =======================================================================
       ! ADJUST MESH INDICES FOR EACH TALLY

       if (t % n_bins(T_MESH) > 0) then
          id = t % mesh
          if (dict_has_key(mesh_dict, id)) then
             t % mesh = dict_get_key(mesh_dict, id)
          else
             message = "Could not find mesh " // trim(int_to_str(id)) // &
                  & " specified on tally " // trim(int_to_str(t % id))
             call fatal_error()
          end if
       end if
    end do

  end subroutine adjust_indices

!===============================================================================
! BUILD_UNIVERSE determines what level each universe is at and determines what
! the parent cell of each cell in a subuniverse is.
!===============================================================================

  recursive subroutine build_universe(univ, parent, level)

    type(Universe), pointer :: univ   ! univese pointer
    integer,     intent(in) :: parent ! cell containing universe
    integer,     intent(in) :: level  ! level of universe

    integer :: i      ! index for cells in universe
    integer :: x,y    ! indices for lattice positions
    integer :: i_cell ! index in cells array
    integer :: universe_num
    type(Cell),     pointer :: c => null()
    type(Universe), pointer :: subuniverse => null()
    type(Lattice),  pointer :: lat => null()

    ! set level of the universe
    univ % level = level

    ! loop over all cells in the universe
    do i = 1, univ % n_cells
       i_cell = univ % cells(i)
       c => cells(i_cell)
       c % parent = parent

       ! if this cell is filled with another universe, recursively
       ! call this subroutine
       if (c % type == CELL_FILL) then
          subuniverse => universes(c % fill)
          call build_universe(subuniverse, i_cell, level + 1)
       end if

       ! if this cell is filled by a lattice, need to build the
       ! universe for each unique lattice element
       if (c % type == CELL_LATTICE) then
          lat => lattices(c % fill)
          do x = 1, lat % n_x
             do y = 1, lat % n_y
                lat => lattices(cells(i_cell) % fill)
                universe_num = lat % element(x,y)
                if (.not. dict_has_key(build_dict, universe_num)) then
                   call dict_add_key(build_dict, universe_num, 0)
                   subuniverse => universes(universe_num)
                   call build_universe(subuniverse, i_cell, level + 1)
                end if
             end do
          end do
       end if

    end do

  end subroutine build_universe

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
            & all(mat%atom_percent < ZERO))) then
          message = "Cannot mix atom and weight percents in material " // &
               & int_to_str(mat % id)
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
               & / MASS_NEUTRON * sum_percent
       end if

       ! Calculate nuclide atom densities and deallocate atom_percent array
       ! since it is no longer needed past this point
       mat % atom_density = mat % density * mat % atom_percent
       deallocate(mat % atom_percent)
    end do

  end subroutine normalize_ao

end module initialize
