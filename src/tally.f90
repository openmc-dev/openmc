module tally

  use constants
  use cross_section, only: get_macro_xs
  use error,         only: fatal_error
  use global
  use mesh,          only: get_mesh_bin, bin_to_mesh_indices
  use output,        only: message, header
  use search,        only: binary_search
  use string,        only: int_to_str, real_to_str
  use tally_header,  only: TallyScore, TallyMapItem, TallyMapElement

#ifdef MPI
  use mpi
  use mpi_routines,  only: reduce_tallies
#endif

  implicit none

  integer, allocatable :: position(:)

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
          keff_std = std
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
! tallies and bins need to be scored to when a particle makes a collision. This
! subroutine also sets the stride attribute for each tally as well as allocating
! storage for the scores array.
!===============================================================================

  subroutine create_tally_map()

    integer :: i                   ! loop index for tallies
    integer :: j                   ! loop index for filter arrays
    integer :: index               ! filter bin entries
    integer :: n                   ! number of bins
    integer :: filter_bins         ! running total of number of filter bins
    integer :: score_bins          ! number of scoring bins
    character(MAX_LINE_LEN) :: msg ! output/error message
    type(TallyObject), pointer :: t => null()

    ! allocate tally map array -- note that we don't need a tally map for the
    ! energy_in and energy_out filters
    allocate(tally_maps(TALLY_TYPES - 3))

    ! allocate list of items for each different filter type
    allocate(tally_maps(T_UNIVERSE) % items(n_universes))
    allocate(tally_maps(T_MATERIAL) % items(n_materials))
    allocate(tally_maps(T_CELL)     % items(n_cells))
    allocate(tally_maps(T_CELLBORN) % items(n_cells))
    allocate(tally_maps(T_SURFACE)  % items(n_surfaces))

    ! Allocate and initialize tally map positioning for finding bins
    allocate(position(TALLY_TYPES - 3))
    position = 0

    do i = 1, n_tallies
       t => tallies(i)

       ! initialize number of scoring bins
       filter_bins = 1

       ! determine if there are subdivisions for incoming or outgoing energy to
       ! adjust the number of filter bins appropriately
       n = t % n_bins(T_ENERGYOUT)
       t % stride(T_ENERGYOUT) = filter_bins
       if (n > 0) then
          filter_bins = filter_bins * n
       end if

       n = t % n_bins(T_ENERGYIN)
       t % stride(T_ENERGYIN) = filter_bins
       if (n > 0) then
          filter_bins = filter_bins * n
       end if

       ! Add map elements for mesh bins
       n = t % n_bins(T_MESH)
       t % stride(T_MESH) = filter_bins
       if (n > 0) then
          filter_bins = filter_bins * n
       end if

       ! Add map elements for surface bins
       n = t % n_bins(T_SURFACE)
       t % stride(T_SURFACE) = filter_bins
       if (n > 0) then
          do j = 1, n
             index = t % surface_bins(j) % scalar
             call add_map_element(tally_maps(T_SURFACE) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for cellborn bins
       n = t % n_bins(T_CELLBORN)
       t % stride(T_CELLBORN) = filter_bins
       if (n > 0) then
          do j = 1, n
             index = t % cellborn_bins(j) % scalar
             call add_map_element(tally_maps(T_CELLBORN) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for cell bins
       n = t % n_bins(T_CELL)
       t % stride(T_CELL) = filter_bins
       if (n > 0) then
          do j = 1, n
             index = t % cell_bins(j) % scalar
             call add_map_element(tally_maps(T_CELL) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for material bins
       n = t % n_bins(T_MATERIAL)
       t % stride(T_MATERIAL) = filter_bins
       if (n > 0) then
          do j = 1, n
             index = t % material_bins(j) % scalar
             call add_map_element(tally_maps(T_MATERIAL) % items(index), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for universe bins
       n = t % n_bins(T_UNIVERSE)
       t % stride(T_UNIVERSE) = filter_bins
       if (n > 0) then
          do j = 1, n
             index = t % universe_bins(j) % scalar
             call add_map_element(tally_maps(T_UNIVERSE) % items(index), i, j)
          end do
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
       t % n_total_bins = filter_bins
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

    integer :: n                       ! size of elements array
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

  subroutine score_tally(p, scatter)

    type(Particle), pointer :: p
    logical, intent(in)     :: scatter

    integer :: i
    integer :: j
    integer :: n
    integer :: bins(TALLY_TYPES)
    integer :: score_index
    real(8) :: score
    real(8) :: last_wgt
    real(8) :: wgt
    logical :: in_mesh
    logical :: has_outgoing_energy
    type(TallyObject),    pointer :: t
    type(StructuredMesh), pointer :: m

    ! Copy particle's pre- and post-collision weight
    last_wgt = p % last_wgt
    wgt = p % wgt

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    do i = 1, n_tallies
       t => tallies(i)
       
       ! =======================================================================
       ! DETERMINE SCORING BIN COMBINATION

       ! determine next universe bin
       if (t % n_bins(T_UNIVERSE) > 0) then
          bins(T_UNIVERSE) = get_next_bin(T_UNIVERSE, p % universe, i)
          if (bins(T_UNIVERSE) == NO_BIN_FOUND) cycle
       else
          bins(T_UNIVERSE) = 1
       end if

       ! determine next material bin
       if (t % n_bins(T_MATERIAL) > 0) then
          bins(T_MATERIAL) = get_next_bin(T_MATERIAL, p % material, i)
          if (bins(T_MATERIAL) == NO_BIN_FOUND) cycle
       else
          bins(T_MATERIAL) = 1
       end if

       ! determine next cell bin
       if (t % n_bins(T_CELL) > 0) then
          bins(T_CELL) = get_next_bin(T_CELL, p % cell, i)
          if (bins(T_CELL) == NO_BIN_FOUND) cycle
       else
          bins(T_CELL) = 1
       end if

       ! determine next cellborn bin
       if (t % n_bins(T_CELLBORN) > 0) then
          bins(T_CELLBORN) = get_next_bin(T_CELLBORN, p % cell_born, i)
          if (bins(T_CELLBORN) == NO_BIN_FOUND) cycle
       else
          bins(T_CELLBORN) = 1
       end if

       ! determine next surface bin
       if (t % n_bins(T_SURFACE) > 0) then
          bins(T_SURFACE) = get_next_bin(T_SURFACE, p % surface, i)
          if (bins(T_SURFACE) == NO_BIN_FOUND) cycle
       else
          bins(T_SURFACE) = 1
       end if

       ! determine mesh bin
       if (t % n_bins(T_MESH) > 0) then
          m => meshes(t % mesh)
          call get_mesh_bin(m, p % xyz, bins(T_MESH), in_mesh)
          if (.not. in_mesh) cycle
       else
          bins(T_MESH) = 1
       end if

       ! determine incoming energy bin
       n = t % n_bins(T_ENERGYIN)
       if (n > 0) then
          ! check if energy of the particle is within energy bins
          if (p % last_E < t % energy_in(1) .or. &
               p % last_E > t % energy_in(n + 1)) cycle

          ! search to find incoming energy bin
          bins(T_ENERGYIN) = binary_search(t % energy_in, n + 1, p % last_E)
       else
          bins(T_ENERGYIN) = 1
       end if

       ! determine outgoing energy bin
       n = t % n_bins(T_ENERGYOUT)
       if (n > 0) then
          ! check if energy of the particle is within energy bins
          if (p % E < t % energy_out(1) .or. p % E > t % energy_out(n + 1)) cycle

          ! search to find incoming energy bin
          bins(T_ENERGYOUT) = binary_search(t % energy_out, n + 1, p % E)
          has_outgoing_energy = .true.
       else
          bins(T_ENERGYOUT) = 1
          has_outgoing_energy = .false.
       end if

       ! =======================================================================
       ! CALCULATE SCORES AND ACCUMULATE TALLY

       ! If we have made it here, we have a scoring combination of bins for this
       ! tally -- now we need to determine where in the scores array we should
       ! be accumulating the tally values

       ! Determine scoring index for this filter combination
       score_index = sum((bins - 1) * t % stride) + 1

       ! Determine score for each bin
       do j = 1, t % n_macro_bins
          if (has_outgoing_energy) then
             ! If this tally has an outgoing energy filter, the only supported
             ! reaction is scattering. For all other reactions, about

             if (.not. scatter) return
             
             ! Make sure bin is scattering -- since scattering has already
             ! occured, we do not need to multiply by the scattering cross
             ! section

             if (t % macro_bins(j) % scalar == MACRO_SCATTER) then
                score = last_wgt
             else
                ! call fatal_error
             end if
          else
             ! For tallies with no outgoing energy filter, the score is
             ! calculated normally depending on the quantity specified

             select case(t % macro_bins(j) % scalar)
             case (MACRO_FLUX)
                score = last_wgt / material_xs % total
             case (MACRO_TOTAL)
                score = last_wgt
             case (MACRO_SCATTER)
                score = last_wgt * (material_xs % total - material_xs % absorption) &
                     / material_xs % total
             case (MACRO_ABSORPTION)
                score = last_wgt * material_xs % absorption / material_xs % total
             case (MACRO_FISSION)
                score = last_wgt * material_xs % fission / material_xs % total
             case (MACRO_NU_FISSION)
                score = last_wgt * material_xs % nu_fission / material_xs % total
             end select
          end if
             
          ! Add score to tally
          call add_to_score(t % scores(score_index, j), score)

       end do

       ! Reset tally map positioning
       position = 0

    end do

  end subroutine score_tally

!===============================================================================
! GET_NEXT_BIN determines the next scoring bin for a particular filter variable
!===============================================================================

  function get_next_bin(i_map, i_item, i_tally) result(bin)

    integer, intent(in) :: i_map
    integer, intent(in) :: i_item
    integer, intent(in) :: i_tally
    integer             :: bin

    integer :: index_tally
    integer :: index_bin
    integer :: n

    ! If there are no scoring bins for this item, then return immediately
    if (.not. allocated(tally_maps(i_map) % items(i_item) % elements)) then
       bin = NO_BIN_FOUND
       return
    end if

    ! Check how many elements there are for this item
    n = size(tally_maps(i_map) % items(i_item) % elements)

    do
       ! Increment position in elements
       position(i_map) = position(i_map) + 1

       ! If we've reached the end of the array, there is no more bin to score to
       if (position(i_map) > n) then
          position(i_map) = 0
          bin = NO_BIN_FOUND
          return
       end if

       index_tally = tally_maps(i_map) % items(i_item) % &
            elements(position(i_map)) % index_tally
       index_bin = tally_maps(i_map) % items(i_item) % &
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

!===============================================================================
! SYNCHRONIZE_TALLIES accumulates the sum of the contributions from each history
! within the cycle to a new random variable
!===============================================================================

  subroutine synchronize_tallies()

    integer :: i   ! index in tallies array
    integer :: j   ! index over filter bins
    integer :: k   ! index over scoring bins
    real(8) :: val ! value of accumulated tally
    type(TallyObject), pointer :: t

#ifdef MPI
    call reduce_tallies()
    if (.not. master) return
#endif

    do i = 1, n_tallies
       t => tallies(i)

       ! Loop over all filter and scoring bins
       do j = 1, t % n_total_bins
          do k = 1, t % n_macro_bins
             ! Add the sum and square of the sum of contributions from each
             ! history within a cycle to the variables val and val_sq. This will
             ! later allow us to calculate a variance on the tallies

             val = t % scores(j,k) % val_history / n_particles
             t % scores(j,k) % val    = t % scores(j,k) % val    + val
             t % scores(j,k) % val_sq = t % scores(j,k) % val_sq + val*val

             ! Reset the within-cycle accumulation variable

             t % scores(j,k) % val_history = ZERO
          end do
       end do

    end do

  end subroutine synchronize_tallies

!===============================================================================
! WRITE_TALLIES creates an output file and writes out the mean values of all
! tallies and their standard deviations
!===============================================================================

  subroutine write_tallies()

    integer :: i                       ! index in tallies array
    integer :: j                       ! level in tally hierarchy
    integer :: k                       ! loop index for scoring bins
    integer :: bins(TALLY_TYPES) = 0   ! bins corresponding to each filter
    integer :: indent                  ! number of spaces to preceed output
    integer :: io_error                ! error in opening/writing file
    integer :: last_filter             ! lowest level filter type
    integer :: score_index             ! index in scores array for filters
    logical :: file_exists             ! does tallies.out file already exists? 
    logical :: has_filter(TALLY_TYPES) ! does tally have this filter?
    character(MAX_LINE_LEN) :: filename                 ! name of output file
    character(15)           :: filter_name(TALLY_TYPES) ! names of tally filters
    character(20)           :: macro_name(6)            ! names of macro scores
    character(80)           :: space = " "              ! spaces
    type(TallyObject), pointer :: t

    ! Initialize names for tally filter types
    filter_name(T_UNIVERSE)  = "Universe"
    filter_name(T_MATERIAL)  = "Material"
    filter_name(T_CELL)      = "Cell"
    filter_name(T_CELLBORN)  = "Birth Cell"
    filter_name(T_SURFACE)   = "Surface"
    filter_name(T_MESH)      = "Mesh"
    filter_name(T_ENERGYIN)  = "Incoming Energy"
    filter_name(T_ENERGYOUT) = "Outgoing Energy"

    ! Initialize names for macro scores
    macro_name(abs(MACRO_FLUX))       = "Flux"
    macro_name(abs(MACRO_TOTAL))      = "Total Reaction Rate"
    macro_name(abs(MACRO_SCATTER))    = "Scattering Rate"
    macro_name(abs(MACRO_ABSORPTION)) = "Absorption Rate"
    macro_name(abs(MACRO_FISSION))    = "Fission Rate"
    macro_name(abs(MACRO_NU_FISSION)) = "Nu-Fission Rate"

    ! Create filename for tally output
    filename = trim(path_input) // "tallies.out"

    ! Check if tally file already exists
    inquire(FILE=filename, EXIST=file_exists)
    if (file_exists) then
       ! Possibly make backup of old tally file
    end if

    ! Open tally file for writing
    open(FILE=filename, UNIT=UNIT_TALLY, STATUS='replace', &
         ACTION='write', IOSTAT=io_error)

    do i = 1, n_tallies
       t => tallies(i)

       ! Write header block
       call header("TALLY " // trim(int_to_str(t % uid)), 3, UNIT_TALLY)

       ! First determine which filters this tally has
       do j = 1, TALLY_TYPES
          if (t % n_bins(j) > 0) then
             has_filter(j) = .true.
             last_filter = j
          else
             has_filter(j) = .false.
          end if
       end do

       ! WARNING: Admittedly, the logic for moving for printing scores is
       ! extremely confusing and took quite a bit of time to get correct. The
       ! logic is structured this way since it is not practical to have a do
       ! loop for each filter variable (given that only a few filters are likely
       ! to be used for a given tally.

       ! Initialize bins, filter level, and indentation
       bins = 0
       j = 1
       indent = 0

       print_bin: do
          find_bin: do
             ! Increment bin combination
             bins(j) = bins(j) + 1

             ! =================================================================
             ! REACHED END OF BINS FOR THIS FILTER, MOVE TO NEXT FILTER

             if ((has_filter(j) .and. bins(j) > t % n_bins(j)) .or. &
                  ((.not. has_filter(j)) .and. bins(j) > 1)) then
                if (j == 1) then
                   ! This means we are done with all bin combinations
                   exit print_bin
                else
                   bins(j) = 0
                   j = j - 1
                   if (has_filter(j)) indent = indent - 2
                end if

             ! =================================================================
             ! VALID BIN -- WRITE FILTER INFORMATION OR EXIT TO WRITE SCORES

             else
                if (j == last_filter) then
                   exit find_bin
                else
                   if (has_filter(j)) then
                      ! Print current filter information
                      write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') space(1:indent), &
                           trim(filter_name(j)), trim(get_label(t, j, bins(j)))
                      indent = indent + 2
                   end if
                   j = j + 1
                end if
             end if

          end do find_bin

          ! Print filter information
          write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') space(1:indent), &
               trim(filter_name(j)), trim(get_label(t, j, bins(j)))

          ! Determine scoring index for this bin combination -- note that unlike
          ! in the score_tally subroutine, we have to use max(bins,1) since all
          ! bins below the lowest filter level will be zeros

          score_index = sum((max(bins,1) - 1) * t % stride) + 1

          ! Write scores for this filter bin combination
          indent = indent + 2
          do k = 1, t % n_macro_bins
             write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') & 
                  space(1:indent), macro_name(abs(t % macro_bins(k) % scalar)), &
                  real_to_str(t % scores(score_index,k) % val), &
                  trim(real_to_str(t % scores(score_index,k) % val_sq))
          end do
          indent = indent - 2

       end do print_bin

    end do

    close(UNIT=UNIT_TALLY)

  end subroutine write_tallies

!===============================================================================
! GET_LABEL returns a label for a cell/surface/etc given a tally, filter type,
! and corresponding bin
!===============================================================================

  function get_label(t, filter_type, bin) result(label)

    type(TallyObject), pointer :: t           ! tally object
    integer, intent(in)        :: filter_type ! type of filter
    integer, intent(in)        :: bin         ! bin in filter array
    character(30)              :: label         ! user-specified identifier

    integer              :: index  ! index in cells/surfaces/etc array
    integer, allocatable :: ijk(:) ! indices in mesh
    real(8)              :: E0     ! lower bound for energy bin
    real(8)              :: E1     ! upper bound for energy bin
    type(StructuredMesh), pointer :: m

    select case(filter_type)
    case (T_UNIVERSE)
       index = t % universe_bins(bin) % scalar
       label = int_to_str(universes(index) % uid)
    case (T_MATERIAL)
       index = t % material_bins(bin) % scalar
       label = int_to_str(materials(index) % uid)
    case (T_CELL)
       index = t % cell_bins(bin) % scalar
       label = int_to_str(cells(index) % uid)
    case (T_CELLBORN)
       index = t % cellborn_bins(bin) % scalar
       label = int_to_str(cells(index) % uid)
    case (T_SURFACE)
       index = t % surface_bins(bin) % scalar
       label = int_to_str(surfaces(index) % uid)
    case (T_MESH)
       m => meshes(t % mesh)
       allocate(ijk(m % n_dimension))
       call bin_to_mesh_indices(m, bin, ijk)
       if (m % n_dimension == 2) then
          label = "Index (" // trim(int_to_str(ijk(1))) // ", " // &
               trim(int_to_str(ijk(2))) // ")"
       elseif (m % n_dimension == 3) then
          label = "Index (" // trim(int_to_str(ijk(1))) // ", " // &
               trim(int_to_str(ijk(2))) // ", " // trim(int_to_str(ijk(3))) // ")"
       end if
    case (T_ENERGYIN)
       E0 = t % energy_in(bin)
       E1 = t % energy_in(bin + 1)
       label = "[" // trim(real_to_str(E0)) // ", " // trim(real_to_str(E1)) // ")"
    case (T_ENERGYOUT)
       E0 = t % energy_out(bin)
       E1 = t % energy_out(bin + 1)
       label = "[" // trim(real_to_str(E0)) // ", " // trim(real_to_str(E1)) // ")"
    end select

  end function get_label

!===============================================================================
! TALLY_STATISTICS computes the mean and standard deviation of the mean of each
! tally and stores them in the val and val_sq attributes of the TallyScores
! respectively
!===============================================================================

  subroutine tally_statistics()

    integer :: i    ! index in tallies array
    integer :: j    ! loop index for filter bins
    integer :: k    ! loop index for scoring bins
    integer :: n    ! number of active cycles
    real(8) :: val  ! sum(x)
    real(8) :: val2 ! sum(x*x)
    real(8) :: mean ! mean value
    real(8) :: std  ! standard deviation of the mean
    type(TallyObject), pointer :: t

    ! Number of active cycles
    n = n_cycles - n_inactive

    do i = 1, n_tallies
       t => tallies(i)

       do j = 1, t % n_total_bins
          do k = 1, t % n_macro_bins
             ! Copy values from tallies
             val  = t % scores(j,k) % val
             val2 = t % scores(j,k) % val_sq

             ! Calculate mean and standard deviation
             mean = val/n
             std = sqrt((val2/n - mean*mean)/n)

             ! Copy back into TallyScore
             t % scores(j,k) % val    = mean
             t % scores(j,k) % val_sq = std
          end do
       end do

    end do

  end subroutine tally_statistics

end module tally
