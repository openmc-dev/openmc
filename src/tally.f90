module tally

  use constants
  use cross_section, only: get_macro_xs
  use error,         only: fatal_error
  use global
  use output,        only: message
  use search,        only: binary_search
  use string,        only: int_to_str
  use tally_header,  only: TallyScore, TallyMapItem, TallyMapElement

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! CALCULATE_KEFF
!===============================================================================

  subroutine calculate_keff(i_cycle)

    integer, intent(in) :: i_cycle ! index of current cycle

    integer(8)              :: total_bank ! total number of source sites
    integer                 :: n          ! active cycle number
    real(8)                 :: kcoll      ! keff collision estimator         
    real(8), save           :: k1 = 0.    ! accumulated keff
    real(8), save           :: k2 = 0.    ! accumulated keff**2
    real(8)                 :: std        ! stdev of keff over active cycles
    character(MAX_LINE_LEN) :: msg        ! output/error message
#ifdef MPI
    integer :: ierr
#endif

    msg = "Calculate cycle keff..."
    call message(msg, 8)

    ! set k1 and k2 at beginning of run
    if (i_cycle == 1) then
       k1 = ZERO
       k2 = ZERO
    end if

#ifdef MPI
    ! Collect number bank sites onto master process
    call MPI_REDUCE(n_bank, total_bank, 1, MPI_INTEGER8, MPI_SUM, 0, &
         & MPI_COMM_WORLD, ierr)
#else
    total_bank = n_bank
#endif

    ! Collect statistics and print output
    if (master) then
       kcoll = real(total_bank)/real(n_particles)*keff
       if (i_cycle > n_inactive) then
          n = i_cycle - n_inactive
          k1 = k1 + kcoll
          k2 = k2 + kcoll**2
          keff = k1/n
          std  = sqrt((k2/n-keff**2)/n)
          if (i_cycle > n_inactive+1) then
             write(6,101) i_cycle, kcoll, keff, std
          else
             write(6,100) i_cycle, kcoll
          end if
       else
          write(6,100) i_cycle, kcoll
          keff = kcoll
       end if
    end if

#ifdef MPI
    call MPI_BCAST(keff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
#endif

100 format (2X,I4,2X,F8.5)
101 format (2X,I4,2X,F8.5,9X,F8.5,1X,F8.5)

  end subroutine calculate_keff

!===============================================================================
! CREATE_TALLY_MAP creates a map that allows a quick determination of which
! tallies and bins need to be scored to when a particle makes a collision.
!===============================================================================

  subroutine create_tally_map()

    integer :: i
    integer :: j
    integer :: index
    integer :: n
    integer :: score_bins
    type(TallyObject), pointer :: t => null()

    ! allocate tally map array
    allocate(tally_maps(TALLY_MAP_TYPES))

    ! allocate list of items for each different filter type
    allocate(tally_maps(MAP_CELL)     % items(n_cells))
    allocate(tally_maps(MAP_SURFACE)  % items(n_surfaces))
    allocate(tally_maps(MAP_UNIVERSE) % items(n_universes))
    allocate(tally_maps(MAP_MATERIAL) % items(n_materials))
    allocate(tally_maps(MAP_MESH)     % items(100)) ! TODO: Change this
    allocate(tally_maps(MAP_BORNIN)   % items(n_cells))

    do i = 1, n_tallies
       t => tallies(i)

       ! initialize number of scoring bins
       score_bins = 1

       ! Add map elements for cell bins
       if (associated(t % cell_bins)) then
          n = size(t % cell_bins)
          do j = 1, n
             index = t % cell_bins(j) % scalar
             call add_map_element(tally_maps(MAP_CELL) % items(index), i, j)
          end do
          score_bins = score_bins * n
       end if

       ! Add map elements for surface bins
       if (associated(t % surface_bins)) then
          n = size(t % surface_bins)
          do j = 1, n
             index = t % surface_bins(j) % scalar
             call add_map_element(tally_maps(MAP_SURFACE) % items(index), i, j)
          end do
          score_bins = score_bins * n
       end if

       ! Add map elements for universe bins
       if (associated(t % universe_bins)) then
          n = size(t % universe_bins)
          do j = 1, n
             index = t % universe_bins(j) % scalar
             call add_map_element(tally_maps(MAP_UNIVERSE) % items(index), i, j)
          end do
          score_bins = score_bins * n
       end if

       ! Add map elements for material bins
       if (associated(t % material_bins)) then
          n = size(t % material_bins)
          do j = 1, n
             index = t % material_bins(j) % scalar
             call add_map_element(tally_maps(MAP_MATERIAL) % items(index), i, j)
          end do
          score_bins = score_bins * n
       end if

       ! Add map elements for bornin bins
       if (associated(t % bornin_bins)) then
          n = size(t % bornin_bins)
          do j = 1, size(t % bornin_bins)
             index = t % bornin_bins(j) % scalar
             call add_map_element(tally_maps(MAP_BORNIN) % items(index), i, j)
          end do
          score_bins = score_bins * n
       end if

       ! TODO: Determine size of mesh to increase number of scoring bins

       ! determine if there are subdivisions for incoming or outgoing energy to
       ! adjust the number of scoring bins appropriately
       if (allocated(t % energy_in)) then
          score_bins = score_bins * (size(t % energy_in) - 1)
       end if

       if (allocated(t % energy_out)) then
          score_bins = score_bins * (size(t % energy_out) - 1)
       end if

       ! Finally add scoring bins for the macro tallies
       if (associated(t % macro_bins)) then
          score_bins = score_bins * size(t % macro_bins)
       end if

       ! Allocate scores for tally
       allocate(t % score(score_bins))

    end do

  end subroutine create_tally_map

!===============================================================================
! ADD_MAP_ELEMENT adds a pair of tally and bin indices to the list for a given
! cell/surface/etc.
!===============================================================================

  subroutine add_map_element(item, index_tally, index_bin)

    type(TallyMapItem), intent(inout) :: item
    integer, intent(in) :: index_tally ! index in tallies array
    integer, intent(in) :: index_bin   ! index in bins array

    integer :: n
    type(TallyMapElement), allocatable :: temp(:)

    if (.not. allocated(item % elements)) then
       allocate(item % elements(1))
       item % elements(1) % index_tally = index_tally
       item % elements(1) % index_bin   = index_bin
    else
       ! determine size of elements array
       n = size(item % elements)

       ! allocate temporary storage and copy elements
       allocate(temp(n+1))
       temp(1:n) = item % elements

       ! move allocation back to main array
       call move_alloc(FROM=temp, TO=item%elements)

       ! set new element
       item % elements(n+1) % index_tally = index_tally
       item % elements(n+1) % index_bin   = index_bin
    end if

  end subroutine add_map_element

!===============================================================================
! SCORE_TALLY
!===============================================================================

  subroutine score_tally(p)

    type(Particle), pointer :: p     ! particle

    integer :: i
    integer :: n
    integer :: cell_bin
    integer :: surface_bin
    integer :: universe_bin
    integer :: material_bin
    integer :: bornin_bin
    integer :: energyin_bin
    integer :: energyout_bin
    type(TallyObject), pointer :: t

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    do i = 1, n_tallies
       t => tallies(i)
       
       ! determine next cell bin
       if (associated(t % cell_bins)) then
          cell_bin = get_next_bin(MAP_CELL, p % cell, i)
          if (cell_bin == NO_BIN_FOUND) cycle
       else
          cell_bin = 1
       end if

       ! determine next surface bin
       if (associated(t % surface_bins)) then
          surface_bin = get_next_bin(MAP_SURFACE, p % surface, i)
          if (surface_bin == NO_BIN_FOUND) cycle
       else
          surface_bin = 1
       end if

       ! determine next universe bin
       if (associated(t % universe_bins)) then
          universe_bin = get_next_bin(MAP_UNIVERSE, p % universe, i)
          if (universe_bin == NO_BIN_FOUND) cycle
       else
          universe_bin = 1
       end if

       ! determine next material bin
       if (associated(t % material_bins)) then
          material_bin = get_next_bin(MAP_MATERIAL, p % material, i)
          if (material_bin == NO_BIN_FOUND) cycle
       else
          material_bin = 1
       end if

       ! determine next bornin bin
       if (associated(t % bornin_bins)) then
          ! bornin_bin = get_next_bin(MAP_BORNIN, p % bornin, i)
          if (bornin_bin == NO_BIN_FOUND) cycle
       else
          bornin_bin = 1
       end if

       ! determine incoming energy bin
       if (allocated(t % energy_in)) then
          ! check if energy of the particle is within energy bins
          n = size(t % energy_in)
          if (p % E < t % energy_in(1) .or. p % E > t % energy_in(n)) cycle

          ! search to find incoming energy bin
          energyin_bin = binary_search(t % energy_in, n, p % E)
       else
          energyin_bin = 1
       end if

       ! determine outgoing energy bin
       if (allocated(t % energy_out)) then
          ! check if energy of the particle is within energy bins
          n = size(t % energy_out)
          if (p % E < t % energy_out(1) .or. p % E > t % energy_out(n)) cycle

          ! search to find incoming energy bin
          energyout_bin = binary_search(t % energy_out, n, p % E)
       else
          energyout_bin = 1
       end if

    end do

  end subroutine score_tally

!===============================================================================
! GET_NEXT_BIN
!===============================================================================

  function get_next_bin(i_map, i_cell, i_tally) result(bin)

    integer, intent(in) :: i_map
    integer, intent(in) :: i_cell
    integer, intent(in) :: i_tally
    integer             :: bin

    integer :: index_tally
    integer :: index_bin
    integer, allocatable, save :: position(:)
    logical, save :: first_entry = .true.
    integer       :: n

    ! initialize position to zero
    if (first_entry) then
       allocate(position(TALLY_MAP_TYPES))
       position = 0
       first_entry = .false.
    end if

    n = size(tally_maps(i_map) % items(i_cell) % elements)

    do
       ! Increment position in elements
       position(i_map) = position(i_map) + 1

       ! If we've reached the end of the array, there is no more bin to score to
       if (position(i_map) > n) then
          position(i_map) = 0
          bin = NO_BIN_FOUND
          return
       end if

       index_tally = tally_maps(i_map) % items(i_cell) % &
            elements(position(i_map)) % index_tally
       index_bin = tally_maps(i_map) % items(i_cell) % &
            elements(position(i_map)) % index_bin

       if (index_tally > i_tally) then
          ! Since the index being checked against is greater than the index we
          ! need (and the tally indices were added to elements sequentially), we
          ! know that no more bins will be scoring bins for this tally
          position(i_map) = 0
          bin = NO_BIN_FOUND
          return
       elseif (index_tally == i_tally) then
          ! Found a match
          bin = index_bin
          return
       end if

    end do

  end function get_next_bin

!===============================================================================
! ADD_TO_SCORE
!===============================================================================

  subroutine add_to_score(score, val)

    type(TallyScore), intent(inout) :: score
    real(8),          intent(in)    :: val
    
!!$    score % n_events = score % n_events + 1
!!$    score % val      = score % val + val
!!$    score % val_sq   = score % val_sq + val*val
    
  end subroutine add_to_score

end module tally
