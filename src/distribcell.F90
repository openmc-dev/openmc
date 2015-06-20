module distribcell

  use constants
  use error,                  only: fatal_error, warning
  use geometry_header,        only: Cell, Surface, Universe, Lattice, &
                                    RectLattice, HexLattice, BASE_UNIVERSE
  use global
  use output,                 only: write_message, find_offset
  use particle_header,        only: LocalCoord, Particle
  use string,                 only: to_str
  use tally_header,           only: TallyResult, TallyMapItem, &
                                    TallyMapElement, TallyFilter

  implicit none

contains

!===============================================================================
! GET_DISTRIBCELL_INSTANCE sums up the distribcell offsets for all coordinate
! levels to determin the instance of the current material
!===============================================================================

  function get_distribcell_instance(p, map_num, target_cell) result(instance)

    type(Particle), intent(in) :: p           ! Particle containing coordinates
    integer,        intent(in) :: map_num     ! Distribcell map number
    integer,        intent(in) :: target_cell ! Target cell to search for
    integer                    :: instance

    integer                   :: offset
    type(LocalCoord), pointer :: coord

    instance = NO_BIN_FOUND

    ! Loop through each layer and sum the offsets for the given cell map
    coord => p % coord0
    offset = 0
    do while(associated(coord))
      if (cells(coord % cell) % type == CELL_FILL) then
        offset = offset + cells(coord % cell) % offset(map_num)
      elseif (cells(coord % cell) % type == CELL_LATTICE) then
        offset = offset + lattices(coord % next % lattice) % obj % &
             offset(map_num, coord % next % lattice_x, &
                    coord % next % lattice_y, coord % next % lattice_z)
      elseif (coord % cell == target_cell) then
        instance = offset + 1
        exit
      end if
      coord => coord % next
    end do
    nullify(coord)

  end function get_distribcell_instance

!===============================================================================
! PREPARE_DISTRIBCELL initializes any distribcell filters present and sets the
! offsets for distribcells
!===============================================================================

  subroutine prepare_distribcell()

    integer :: i, j       ! Tally, filter loop counters
    integer :: n_filt     ! Number of filters originally in tally
    logical :: count_all  ! Count all cells
    type(TallyObject),    pointer :: tally            ! Current tally
    type(Universe),       pointer :: univ             ! Pointer to universe
    type(Cell),           pointer :: c                ! Pointer to cell
    type(Material),       pointer :: mat              ! Pointer to material
    integer, allocatable :: univ_list(:)              ! Target offsets
    integer, allocatable :: counts(:,:)               ! Target count
    logical, allocatable :: found(:,:)                ! Target found

    count_all = .false.

    ! Loop over tallies    
    do i = 1, n_tallies

      ! Get pointer to tally
      tally => tallies(i)

      n_filt = tally % n_filters

      ! Loop over the filters to determine how many additional filters
      ! need to be added to this tally
      do j = 1, tally % n_filters

        ! Determine type of filter
        if (tally % filters(j) % type == FILTER_DISTRIBCELL) then
          count_all = .true.
          if (size(tally % filters(j) % int_bins) > 1) then
            call fatal_error("A distribcell filter was specified with &
                             &multiple bins. This feature is not supported.")
          end if      
        end if

      end do

    end do
    
    if (count_all) then
      
      univ => universes(BASE_UNIVERSE)

      ! sum the number of occurrences of all cells
      call count_instance(univ)

      ! Loop over tallies    
      do i = 1, n_tallies    

        ! Get pointer to tally
        tally => tallies(i)      

        ! Initialize the filters
        do j = 1, tally % n_filters

          ! Set the number of bins to the number of instances of the cell
          if (tally % filters(j) % type == FILTER_DISTRIBCELL) then
            c => cells(tally % filters(j) % int_bins(1))
            tally % filters(j) % n_bins = c % instances
          end if

        end do
      end do

    end if

    ! Allocate offset maps at each level in the geometry
    call allocate_offsets(univ_list, counts, found)

    ! Verify correct xml input of distributed materials
    call verify_distribmats()

    ! Calculate offsets for each target distribcell
    do i = 1, n_maps
      do j = 1, n_universes  
        univ => universes(j)
        call calc_offsets(univ_list(i), i, univ, counts, found)
      end do
    end do

    ! Assign all non-distributed materials to the same map
    do i = 1, n_materials
      mat => materials(i)
      if (.not. mat % distrib_comp) then
        mat % distribmap = n_maps
      end if
    end do

    ! Deallocate temporary target variable arrays
    deallocate(counts)
    deallocate(found)
    deallocate(univ_list)
  
  end subroutine prepare_distribcell

!===============================================================================
! ALLOCATE_OFFSETS determines the number of maps needed and allocates required
! memory for distribcell offset tables
!===============================================================================

  recursive subroutine allocate_offsets(univ_list, counts, found)

    integer, intent(out), allocatable     :: univ_list(:) ! Target offsets
    integer, intent(out), allocatable     :: counts(:,:)  ! Target count
    logical, intent(out), allocatable     :: found(:,:)   ! Target found

    integer :: mapnum                           ! Map index
    integer :: i, j, l, m                       ! Loop counters
    type(SetInt)               :: cell_list     ! distribells to track    
    type(Universe),    pointer :: univ          ! pointer to universe
    class(Lattice),    pointer :: lat           ! pointer to lattice
    type(TallyObject), pointer :: tally         ! pointer to tally
    type(TallyFilter), pointer :: filter        ! pointer to filter
    type(Material),    pointer :: mat           ! pointer to material
    type(Cell),        pointer :: c             ! pointer to cell
    
    ! Begin gathering list of cells in distribcell tallies
    n_maps = 0
    
    ! Populate list of distribcells to track from tallies
    do i = 1, n_tallies
      tally => tallies(i)

      do j = 1, tally % n_filters
        filter => tally % filters(j)        

        if (filter % type == FILTER_DISTRIBCELL) then
          if (.not. cell_list % contains(filter % int_bins(1))) then
            call cell_list % add(filter % int_bins(1))
          end if
        end if 

      end do
    end do

    ! Include list of distribcells to track from materials
    do i = 1, n_materials
        mat => materials(i)

        ! Skip non-distributed materials
        if (.not. (mat % distrib_dens .or. mat % distrib_comp)) cycle

        ! Find the cell this distributed material is assigned to
        do j = 1, n_cells
          c => cells(j)

          ! Skip non-normal cells
          if (.not. (c % type == CELL_NORMAL)) cycle

          ! Check if this cell was assigned this material
          if (c % material == material_dict % get_key(mat % id)) then

            ! Save the cell id on the material, and enforce that no two cells 
            ! share a distributed material
            if (mat % distribcell == NONE) then
              mat % distribcell = c % id
              call cell_list % add(c % id)
            else
              call fatal_error("Two cells, " // trim(to_str(mat % distribcell)) // &
                               " and " // trim(to_str(c % id)) // &
                               ", cannot share the same distributed material.")
            end if

          end if

        end do

        if (mat % distribcell < 1) then        
          if (master) call warning("No cell found for distributed material: " // &
              to_str(mat % id))
        end if        

    end do

    ! Compute the number of unique universes containing these distribcells
    ! to determine the number of offset tables to allocate
    do i = 1, n_universes
      univ => universes(i)
      do j = 1, univ % n_cells
        if (cell_list % contains(univ % cells(j))) then
          n_maps = n_maps + 1
        end if
      end do
    end do

    ! Create an extra map for non-distributed materials
    n_maps = n_maps + 1
    
    ! Allocate the list of offset tables for each unique universe
    allocate(univ_list(n_maps))

    ! Allocate list to accumulate target distribcell counts in each universe
    allocate(counts(n_universes, n_maps))

    ! Allocate list to track if target distribcells are found in each universe
    allocate(found(n_universes, n_maps))

    counts(:,:) = 0
    found(:,:) = .false.
    mapnum = 1

    do i = 1, n_universes
      univ => universes(i)

      do j = 1, univ % n_cells
      
        if (cell_list % contains(univ % cells(j))) then
          
          ! Loop over all tallies    
          do l = 1, n_tallies
            tally => tallies(l)
            
            do m = 1, tally % n_filters
              filter => tally % filters(m)
              
              ! Loop over only distribcell filters
              ! If filter points to cell we just found, set offset index
              if (filter % type == FILTER_DISTRIBCELL) then                  
                if (filter % int_bins(1) == univ % cells(j)) then
                  filter % offset = mapnum
                end if
              end if

            end do
          end do
     
          ! Loop over all materials
          do l = 1, n_materials
            mat => materials(l)

            ! Skip non-distributed mats
            if (mat % distribcell == NONE) cycle

            ! If this material is in the current cell, store its corresponding
            ! map index
            if (cell_dict % get_key(mat % distribcell) == univ % cells(j)) then
              mat % distribmap = mapnum
            end if

          end do
   
          univ_list(mapnum) = univ % id
          mapnum = mapnum + 1
        end if
      end do
    end do
    
    ! Allocate the offset tables for lattices    
    do i = 1, n_lattices
      lat => lattices(i) % obj

      select type(lat)

      type is (RectLattice)
        allocate(lat % offset(n_maps, lat % n_cells(1), lat % n_cells(2), &
                 lat % n_cells(3)))
      type is (HexLattice)
        allocate(lat % offset(n_maps, 2 * lat % n_rings - 1, &
             2 * lat % n_rings - 1, lat % n_axial))
      end select

      lat % offset(:, :, :, :) = 0

    end do

    ! Allocate offset table for fill cells
    do i = 1, n_cells
      if (cells(i) % material == NONE) then
        allocate(cells(i) % offset(n_maps))
      end if
    end do

  end subroutine allocate_offsets

!===============================================================================
! VERIFY_DISTRIBMATS verifies that all inputs are correct and then initializes
! the mappings for distributed materials
!===============================================================================

  subroutine verify_distribmats()

    integer :: i                                    ! Primary loop index
    integer :: j                                    ! Additional loop index
    integer :: num                                  ! size storage
    real(8) :: density                              ! density to use
    real(8),allocatable  :: atom_density(:)         ! composition to use
    type(Cell),     pointer, save :: c => null()    ! pointer to cell
    type(Material), pointer, save :: mat => null()  ! pointer to material

    call write_message("Verifying distributed materials...", 7)

    ! Verify that all distributed materials have a composition / density length
    ! equal to either 1 or the number of instance
    do i = 1, n_materials

      mat => materials(i)
      call write_message("Verifying mat " // trim(to_str(mat % id)) // "...", 9)

      if (mat % distrib_dens) then

        ! Skip unused mats
        if (mat % distribcell < 1) cycle

        c => cells(cell_dict % get_key(mat % distribcell))

        num = mat % density % num
        ! Ensure that there are a sensible number of densities specified
        if (.not.(num == 1 .or. num == c % instances)) then  

          call fatal_error("Invalid number of densities specified for " // &
                           "material " // to_str(mat % id))

        end if

        ! If num == 1, set all densities equal to the one given
        if (num == 1) then

          density = mat % density % density(1)
          deallocate(mat % density % density)
          allocate(mat % density % density(c % instances))
          do j = 1, c % instances
            mat % density % density(j) = density          
          end do          
          mat % density % num = c % instances

        end if

        ! Distribute the density by creating a fake composition distribution
        ! of the composition provided. Later the normalize_ao function will
        ! distribute the densities on the composition
        mat % distrib_comp = .true.
        mat % n_comp = c % instances
        allocate(atom_density(mat % n_nuclides))
        atom_density = mat % comp(1) % atom_density
        deallocate(mat % comp(1) % atom_density)
        deallocate(mat % comp)
        allocate(mat % comp(c % instances))
        do j = 1, c % instances
          allocate(mat % comp(j) % atom_density(mat % n_nuclides))
          mat % comp(j) % atom_density = atom_density          
        end do
        deallocate(atom_density)

      else if (mat % distrib_comp) then

        ! Skip unused mats
        if (mat % distribcell < 1) cycle

        c => cells(cell_dict % get_key(mat % distribcell))

        num = mat % n_comp

        ! Ensure that there are a sensible number of compositions specified
        if (.not.(num == 1 .or. num == c % instances)) then 

          call fatal_error("Invalid number of compositions specified for " // &
                           "material " // to_str(mat % id))

        end if

        ! If num == 1, set all compositions equal to the one given
        if (num == 1 .and. .not. mat % otf_compositions) then

          mat % n_comp = c % instances
          allocate(atom_density(mat % n_nuclides))
          atom_density = mat % comp(1) % atom_density
          deallocate(mat % comp(1) % atom_density)
          deallocate(mat % comp)
          allocate(mat % comp(c % instances))
          do j = 1, c % instances          
            allocate(mat % comp(j) % atom_density(mat % n_nuclides))
            mat % comp(j) % atom_density = atom_density          
          end do
          deallocate(atom_density)

        end if

      else 
        mat % distribmap = n_maps  
      end if

    end do

  end subroutine verify_distribmats

!===============================================================================
! CALC_OFFSETS calculates and stores the offsets in all fill cells. This
! routine is called once upon initialization.
!===============================================================================

  subroutine calc_offsets(goal, map, univ, counts, found)

    integer, intent(in)        :: goal         ! target universe ID
    integer, intent(in)        :: map          ! map index in vector of maps
    type(Universe), intent(in) :: univ         ! universe searching in
    integer, intent(inout)     :: counts(:,:)  ! target count
    logical, intent(inout)     :: found(:,:)   ! target found

    integer :: i                          ! index over cells
    integer :: j, k, m                    ! indices in lattice
    integer :: n                          ! number of cells to search
    integer :: offset                     ! total offset for a given cell
    integer :: cell_index                 ! index in cells array
    type(Cell),     pointer :: c          ! pointer to current cell
    type(Universe), pointer :: next_univ  ! next universe to cycle through
    class(Lattice), pointer :: lat        ! pointer to current lattice

    n = univ % n_cells
    offset = 0

    do i = 1, n

      cell_index = univ % cells(i)

      ! get pointer to cell
      c => cells(cell_index)

      ! ====================================================================
      ! AT LOWEST UNIVERSE, TERMINATE SEARCH
      if (c % type == CELL_NORMAL) then

      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_FILL) then
        ! Set offset for the cell on this level
        c % offset(map) = offset

        ! Count contents of this cell
        next_univ => universes(c % fill)
        offset = offset + count_target(next_univ, counts, found, goal, map)

        ! Move into the next universe
        next_univ => universes(c % fill)
        c => cells(cell_index)

      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_LATTICE) then

        ! Set current lattice
        lat => lattices(c % fill) % obj

        select type (lat)

        type is (RectLattice)

          ! Loop over lattice coordinates
          do j = 1, lat % n_cells(1)
            do k = 1, lat % n_cells(2)
              do m = 1, lat % n_cells(3)
                lat % offset(map, j, k, m) = offset
                next_univ => universes(lat % universes(j, k, m))
                offset = offset + &
                     count_target(next_univ, counts, found, goal, map)
              end do
            end do
          end do

        type is (HexLattice)

          ! Loop over lattice coordinates
          do m = 1, lat % n_axial
            do k = 1, 2*lat % n_rings - 1
              do j = 1, 2*lat % n_rings - 1
                ! This array location is never used
                if (j + k < lat % n_rings + 1) then
                  cycle
                ! This array location is never used
                else if (j + k > 3*lat % n_rings - 1) then
                  cycle
                else
                  lat % offset(map, j, k, m) = offset
                  next_univ => universes(lat % universes(j, k, m))
                  offset = offset + &
                       count_target(next_univ, counts, found, goal, map)
                end if
              end do
            end do
          end do
        end select

      end if
    end do

  end subroutine calc_offsets

!===============================================================================
! COUNT_TARGET recursively totals the numbers of occurances of a given
! universe ID beginning with the universe given.
!===============================================================================

  recursive function count_target(univ, counts, found, goal, map) result(count)

    type(Universe), intent(in) :: univ         ! universe to search through
    integer, intent(inout)     :: counts(:,:)  ! target count
    logical, intent(inout)     :: found(:,:)   ! target found
    integer, intent(in)        :: goal         ! target universe ID
    integer, intent(in)        :: map          ! current map

    integer :: i                           ! index over cells
    integer :: j, k, m                     ! indices in lattice
    integer :: n                           ! number of cells to search
    integer :: cell_index                  ! index in cells array
    integer :: count                       ! number of times target located
    type(Cell),     pointer :: c           ! pointer to current cell
    type(Universe), pointer :: next_univ   ! next univ to loop through
    class(Lattice), pointer :: lat         ! pointer to current lattice

    ! Don't research places already checked
    if (found(universe_dict % get_key(univ % id), map)) then
      count = counts(universe_dict % get_key(univ % id), map)
      return
    end if

    ! If this is the target, it can't contain itself.
    ! Count = 1, then quit
    if (univ % id == goal) then
      count = 1
      counts(universe_dict % get_key(univ % id), map) = 1
      found(universe_dict % get_key(univ % id), map) = .true.
      return
    end if

    count = 0
    n = univ % n_cells

    do i = 1, n

      cell_index = univ % cells(i)

      ! get pointer to cell
      c => cells(cell_index)

      ! ====================================================================
      ! AT LOWEST UNIVERSE, TERMINATE SEARCH
      if (c % type == CELL_NORMAL) then

      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_FILL) then

        next_univ => universes(c % fill)

        ! Found target - stop since target cannot contain itself
        if (next_univ % id == goal) then
          count = count + 1
          return
        end if

        count = count + count_target(next_univ, counts, found, goal, map)
        c => cells(cell_index)

      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_LATTICE) then

        ! Set current lattice
        lat => lattices(c % fill) % obj

        select type (lat)

        type is (RectLattice)

          ! Loop over lattice coordinates
          do j = 1, lat % n_cells(1)
            do k = 1, lat % n_cells(2)
              do m = 1, lat % n_cells(3)
                next_univ => universes(lat % universes(j, k, m))

                ! Found target - stop since target cannot contain itself
                if (next_univ % id == goal) then
                  count = count + 1
                  cycle
                end if

                count = count + &
                     count_target(next_univ, counts, found, goal, map)

              end do
            end do
          end do

          type is (HexLattice)

            ! Loop over lattice coordinates
            do m = 1, lat % n_axial
              do k = 1, 2*lat % n_rings - 1
                do j = 1, 2*lat % n_rings - 1
                  ! This array location is never used
                  if (j + k < lat % n_rings + 1) then
                    cycle
                  ! This array location is never used
                  else if (j + k > 3*lat % n_rings - 1) then
                    cycle
                  else
                    next_univ => universes(lat % universes(j, k, m))

                    ! Found target - stop since target cannot contain itself
                    if (next_univ % id == goal) then
                      count = count + 1
                      cycle
                    end if

                    count = count + &
                         count_target(next_univ, counts, found, goal, map)
                  end if
                end do
              end do
            end do

          end select

      end if
    end do

    counts(universe_dict % get_key(univ % id), map) = count
    found(universe_dict % get_key(univ % id), map) = .true.

  end function count_target

!===============================================================================
! COUNT_INSTANCE recursively totals the number of occurrences of all cells
! beginning with the universe given.
!===============================================================================

  recursive subroutine count_instance(univ)

    type(Universe), intent(in) :: univ  ! universe to search through

    integer :: i                          ! index over cells
    integer :: j, k, m                    ! indices in lattice
    integer :: n                          ! number of cells to search
    integer :: cell_index                 ! index in cells array
    type(Cell),     pointer :: c          ! pointer to current cell
    type(Universe), pointer :: next_univ  ! next universe to loop through
    class(Lattice), pointer :: lat        ! pointer to current lattice

    n = univ % n_cells

    do i = 1, n

      cell_index = univ % cells(i)

      ! get pointer to cell
      c => cells(cell_index)
      c % instances = c % instances + 1

      ! ====================================================================
      ! AT LOWEST UNIVERSE, TERMINATE SEARCH
      if (c % type == CELL_NORMAL) then

      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_FILL) then

        next_univ => universes(c % fill)

        call count_instance(next_univ)
        c => cells(cell_index)

      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_LATTICE) then

        ! Set current lattice
        lat => lattices(c % fill) % obj

        select type (lat)

        type is (RectLattice)

          ! Loop over lattice coordinates
          do j = 1, lat % n_cells(1)
            do k = 1, lat % n_cells(2)
              do m = 1, lat % n_cells(3)
                next_univ => universes(lat % universes(j, k, m))
                call count_instance(next_univ)
              end do
            end do
          end do

        type is (HexLattice)

          ! Loop over lattice coordinates
          do m = 1, lat % n_axial
            do k = 1, 2*lat % n_rings - 1
              do j = 1, 2*lat % n_rings - 1
                ! This array location is never used
                if (j + k < lat % n_rings + 1) then
                  cycle
                ! This array location is never used
                else if (j + k > 3*lat % n_rings - 1) then
                  cycle
                else
                  next_univ => universes(lat % universes(j, k, m))
                  call count_instance(next_univ)
                end if
              end do
            end do
          end do

        end select

      end if
    end do

  end subroutine count_instance

!===============================================================================
! DISTRIBUTION_HELP prints the set of human readable paths for all distributed
! materials.  This assists in building/debugging a model with distribcells.
!===============================================================================

  subroutine distribution_help()

    integer :: i      ! materials loop
    integer :: j      ! instances loop
    integer :: offset ! offset parameter for path generation
    type(Material), pointer :: mat   ! pointer to material
    type(Cell),     pointer :: c     ! pointer to cell
    type(Universe), pointer :: univ  ! pointer to universe
    character(MAX_FILE_LEN) :: path  ! path of summary file
    character(100)          :: label ! user-specified identifier

    ! Create filename for log file
    path = trim(path_output) // "distribution.out"

    ! Open log file for writing
    open(UNIT=UNIT_HELP, FILE=path, STATUS='replace', ACTION='write')

    do i = 1, n_materials

      mat => materials(i)

      if (mat % distrib_dens .or. mat % distrib_comp) then

        ! Skip unused mats
        if (mat % distribcell < 1) cycle

        c => cells(cell_dict % get_key(mat % distribcell))

        write(UNIT_HELP,*) 'Distributed Material:', mat % id
        write(UNIT_HELP,*) 'Number of Instances:', c % instances

        do j = 1, c % instances

          offset = 0
          label = ''
          univ => universes(BASE_UNIVERSE)
          call find_offset(mat % distribmap, material_dict % get_key(mat % id), univ, j - 1, offset, label, .true.)
          write(UNIT_HELP,*) label

        end do

      end if

    end do

  end subroutine distribution_help
 
end module distribcell
