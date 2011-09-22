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
    integer :: filter_bins
    integer :: score_bins
    character(MAX_LINE_LEN) :: msg        ! output/error message
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
       filter_bins = 1

       ! Add map elements for cell bins
       n = t % n_cell_bins
       if (n > 0) then
          do j = 1, n
             index = t % cell_bins(j) % scalar
             call add_map_element(tally_maps(MAP_CELL) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for surface bins
       n = t % n_surface_bins
       if (n > 0) then
          do j = 1, n
             index = t % surface_bins(j) % scalar
             call add_map_element(tally_maps(MAP_SURFACE) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for universe bins
       n = t % n_universe_bins
       if (n > 0) then
          do j = 1, n
             index = t % universe_bins(j) % scalar
             call add_map_element(tally_maps(MAP_UNIVERSE) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for material bins
       n = t % n_material_bins
       if (n > 0) then
          do j = 1, n
             index = t % material_bins(j) % scalar
             call add_map_element(tally_maps(MAP_MATERIAL) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for bornin bins
       n = t % n_bornin_bins
       if (n > 0) then
          do j = 1, n
             index = t % bornin_bins(j) % scalar
             call add_map_element(tally_maps(MAP_BORNIN) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! TODO: Determine size of mesh to increase number of scoring bins

       ! determine if there are subdivisions for incoming or outgoing energy to
       ! adjust the number of scoring bins appropriately
       n = t % n_energy_in
       if (n > 0) then
          filter_bins = filter_bins * n
       end if

       n = t % n_energy_out
       if (n > 0) then
          filter_bins = filter_bins * n
       end if

       ! Finally add scoring bins for the macro tallies
       n = t % n_macro_bins
       if (n > 0) then
          score_bins = n
       else
          msg = "Must have macro tally bins!"
          call fatal_error(msg)
       end if

       ! Allocate scores for tally
       allocate(t % scores(filter_bins, score_bins))

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
! SCORE_TALLY contains the main logic for scoring user-specified tallies
!===============================================================================

  subroutine score_tally(p)

    type(Particle), pointer :: p     ! particle

    integer :: i
    integer :: j
    integer :: n
    integer :: cell_bin
    integer :: surface_bin
    integer :: universe_bin
    integer :: material_bin
    integer :: bornin_bin
    integer :: energyin_bin
    integer :: energyout_bin

    integer :: score_index
    real(8) :: score
    type(TallyObject), pointer :: t

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    do i = 1, n_tallies
       t => tallies(i)
       
       ! =======================================================================
       ! DETERMINE SCORING BIN COMBINATION

       ! determine next cell bin
       if (t % n_cell_bins > 0) then
          cell_bin = get_next_bin(MAP_CELL, p % cell, i)
          if (cell_bin == NO_BIN_FOUND) cycle
       else
          cell_bin = 1
       end if

       ! determine next surface bin
       if (t % n_surface_bins > 0) then
          surface_bin = get_next_bin(MAP_SURFACE, p % surface, i)
          if (surface_bin == NO_BIN_FOUND) cycle
       else
          surface_bin = 1
       end if

       ! determine next universe bin
       if (t % n_universe_bins > 0) then
          universe_bin = get_next_bin(MAP_UNIVERSE, p % universe, i)
          if (universe_bin == NO_BIN_FOUND) cycle
       else
          universe_bin = 1
       end if

       ! determine next material bin
       if (t % n_material_bins > 0) then
          material_bin = get_next_bin(MAP_MATERIAL, p % material, i)
          if (material_bin == NO_BIN_FOUND) cycle
       else
          material_bin = 1
       end if

       ! determine next bornin bin
       if (t % n_bornin_bins > 0) then
          bornin_bin = get_next_bin(MAP_BORNIN, p % cell_born, i)
          if (bornin_bin == NO_BIN_FOUND) cycle
       else
          bornin_bin = 1
       end if

       ! determine incoming energy bin
       n = t % n_energy_in
       if (n > 0) then
          ! check if energy of the particle is within energy bins
          if (p % E < t % energy_in(1) .or. p % E > t % energy_in(n)) cycle

          ! search to find incoming energy bin
          energyin_bin = binary_search(t % energy_in, n, p % E)
       else
          energyin_bin = 1
       end if

       ! determine outgoing energy bin
       n = t % n_energy_out
       if (n > 0) then
          ! check if energy of the particle is within energy bins
          if (p % E < t % energy_out(1) .or. p % E > t % energy_out(n)) cycle

          ! search to find incoming energy bin
          energyout_bin = binary_search(t % energy_out, n, p % E)
       else
          energyout_bin = 1
       end if

       ! =======================================================================
       ! DETERMINE INDEX IN SCORES ARRAY

       ! If we have made it here, we have a scoring combination of bins for this
       ! tally -- now we need to determine where in the scores array we should
       ! be accumulating the tally values

       score_index = bins_to_index(t, cell_bin, surface_bin, universe_bin, &
            material_bin, bornin_bin, energyin_bin, energyout_bin)

       ! =======================================================================
       ! CALCULATE SCORES AND ACCUMULATE TALLY

       do j = 1, t % n_macro_bins
          ! Determine score
          select case(t % macro_bins(j) % scalar)
          case (MACRO_FLUX)
             score = p % wgt / material_xs % total
          case (MACRO_TOTAL)
             score = p % wgt
          case (MACRO_SCATTER)
             score = p % wgt * material_xs % scatter / material_xs % total
          case (MACRO_ABSORPTION)
             score = p % wgt * material_xs % absorption / material_xs % total
          case (MACRO_FISSION)
             score = p % wgt * material_xs % fission / material_xs % total
          case (MACRO_NU_FISSION)
             score = p % wgt * material_xs % nu_fission / material_xs % total
          end select

          ! Add score to tally
          call add_to_score(t % scores(score_index, j), score)

       end do

    end do

  end subroutine score_tally

!===============================================================================
! GET_NEXT_BIN determines the next scoring bin for a particular filter variable
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
! BINS_TO_INDEX takes a combination of cell/surface/etc bins and returns the
! index in the scores array for the given tally
!===============================================================================

  function bins_to_index(t, c_bin, s_bin, u_bin, m_bin, b_bin, &
       ei_bin, eo_bin) result (index)

    type(TallyObject), pointer :: t
    integer, intent(in) :: c_bin
    integer, intent(in) :: s_bin
    integer, intent(in) :: u_bin
    integer, intent(in) :: m_bin
    integer, intent(in) :: b_bin
    integer, intent(in) :: ei_bin
    integer, intent(in) :: eo_bin
    integer             :: index

    integer :: product

    index = 0
    product = 1
    
    index = index + (c_bin - 1)
    if (t % n_cell_bins > 0) then
       product = product * t % n_cell_bins
    end if

    index = index + (s_bin - 1) * product
    if (t % n_surface_bins > 0) then
       product = product * t % n_surface_bins
    end if

    index = index + (u_bin - 1) * product
    if (t % n_universe_bins > 0) then
       product = product * t % n_universe_bins
    end if

    index = index + (m_bin - 1) * product
    if (t % n_material_bins > 0) then
       product = product * t % n_material_bins
    end if

    index = index + (b_bin - 1) * product
    if (t % n_bornin_bins > 0) then
       product = product * t % n_bornin_bins
    end if

    index = index + (ei_bin - 1) * product
    if (t % n_energy_in > 9) then
       product = product * t % n_energy_in
    end if

    ! We need to shift the index by one since the array in fortran starts at
    ! unity, not zero
    index = index + (eo_bin - 1) * product + 1
    
  end function bins_to_index

!===============================================================================
! ADD_TO_SCORE accumulates a scoring contribution to a specific tally bin and
! specific response function. Note that we don't need to add the square of the
! contribution since that is done at the cycle level, not the history level
!===============================================================================

  subroutine add_to_score(score, val)

    type(TallyScore), intent(inout) :: score
    real(8),          intent(in)    :: val
    
    score % n_events    = score % n_events    + 1
    score % val_history = score % val_history + val
    
  end subroutine add_to_score

end module tally
