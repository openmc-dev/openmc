module tally

  use ISO_FORTRAN_ENV

  use constants
  use error,         only: fatal_error
  use global
  use mesh,          only: get_mesh_bin, bin_to_mesh_indices, get_mesh_indices
  use mesh_header,   only: StructuredMesh
  use output,        only: write_message, header
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
! CALCULATE_KEFF calculates the single cycle estimate of keff as well as the
! mean and standard deviation of the mean for active cycles and displays them
!===============================================================================

  subroutine calculate_keff(i_cycle)

    integer, intent(in) :: i_cycle ! index of current cycle

    integer(8)    :: total_bank ! total number of source sites
    integer       :: n          ! active cycle number
    real(8)       :: k_cycle    ! single cycle estimate of keff
    real(8), save :: k_sum      ! accumulated keff
    real(8), save :: k_sum_sq   ! accumulated keff**2

    message = "Calculate cycle keff..."
    call write_message(8)

    ! initialize sum and square of sum at beginning of run
    if (i_cycle == 1) then
       k_sum = ZERO
       k_sum_sq = ZERO
    end if

#ifdef MPI
    ! Collect number bank sites onto master process
    call MPI_REDUCE(n_bank, total_bank, 1, MPI_INTEGER8, MPI_SUM, 0, &
         & MPI_COMM_WORLD, mpi_err)
#else
    total_bank = n_bank
#endif

    ! Collect statistics and print output
    if (master) then
       ! Since the creation of bank sites was originally weighted by the last
       ! cycle keff, we need to multiply by that keff to get the current cycle's
       ! value

       k_cycle = real(total_bank)/real(n_particles)*keff

       if (i_cycle > n_inactive) then
          ! Active cycle number
          n = i_cycle - n_inactive

          ! Accumulate cycle estimate of k
          k_sum =    k_sum    + k_cycle
          k_sum_sq = k_sum_sq + k_cycle*k_cycle

          ! Determine mean and standard deviation of mean
          keff = k_sum/n
          keff_std = sqrt((k_sum_sq/n - keff*keff)/n)

          ! Display output for this cycle
          if (i_cycle > n_inactive+1) then
             write(UNIT=OUTPUT_UNIT, FMT=101) i_cycle, k_cycle, keff, keff_std
          else
             write(UNIT=OUTPUT_UNIT, FMT=100) i_cycle, k_cycle
          end if
       else
          ! Display output for inactive cycle
          write(UNIT=OUTPUT_UNIT, FMT=100) i_cycle, k_cycle
          keff = k_cycle
       end if
    end if

#ifdef MPI
    ! Broadcast new keff value to all processors
    call MPI_BCAST(keff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
#endif

100 format (2X,I5,2X,F8.5)
101 format (2X,I5,2X,F8.5,9X,F8.5,1X,F8.5)

  end subroutine calculate_keff

!===============================================================================
! CREATE_TALLY_MAP creates a map that allows a quick determination of which
! tallies and bins need to be scored to when a particle makes a collision. This
! subroutine also sets the stride attribute for each tally as well as allocating
! storage for the scores array.
!===============================================================================

  subroutine create_tally_map()

    integer :: i           ! loop index for tallies
    integer :: j           ! loop index for filter arrays
    integer :: i_item      ! filter bin entries
    integer :: n           ! number of bins
    integer :: filter_bins ! running total of number of filter bins
    integer :: score_bins  ! number of scoring bins
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()

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

       ! initialize number of filter bins
       filter_bins = 1

       ! handle surface current tallies separately
       if (t % surface_current) then
          m => meshes(t % mesh)

          t % stride(TS_SURFACE) = filter_bins
          ! Set stride for surface/direction
          if (m % n_dimension == 2) then
             filter_bins = filter_bins * 4
          elseif (m % n_dimension == 3) then
             filter_bins = filter_bins * 6
          end if

          ! Add filter for incoming energy
          n = t % n_bins(T_ENERGYIN)
          t % stride(TS_ENERGYIN) = filter_bins
          if (n > 0) then
             filter_bins = filter_bins * n
          end if
          
          ! account for z direction
          t % stride(TS_MESH_Z) = filter_bins
          filter_bins = filter_bins * (m % dimension(3) + 1)

          ! account for y direction
          t % stride(TS_MESH_Y) = filter_bins
          filter_bins = filter_bins * (m % dimension(2) + 1)

          ! account for z direction
          t % stride(TS_MESH_X) = filter_bins
          filter_bins = filter_bins * (m % dimension(1) + 1)

          ! Finally add scoring bins for the macro tallies and allocate scores
          ! for tally
          score_bins = t % n_macro_bins
          t % n_total_bins = filter_bins
          allocate(t % scores(filter_bins, score_bins))

          ! All the remaining logic is for non-surface-current tallies so we can
          ! safely skip it
          cycle
       end if

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
             i_item = t % surface_bins(j) % scalar
             call add_map_element(tally_maps(T_SURFACE) % items(i_item), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for cellborn bins
       n = t % n_bins(T_CELLBORN)
       t % stride(T_CELLBORN) = filter_bins
       if (n > 0) then
          do j = 1, n
             i_item = t % cellborn_bins(j) % scalar
             call add_map_element(tally_maps(T_CELLBORN) % items(i_item), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for cell bins
       n = t % n_bins(T_CELL)
       t % stride(T_CELL) = filter_bins
       if (n > 0) then
          do j = 1, n
             i_item = t % cell_bins(j) % scalar
             call add_map_element(tally_maps(T_CELL) % items(i_item), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for material bins
       n = t % n_bins(T_MATERIAL)
       t % stride(T_MATERIAL) = filter_bins
       if (n > 0) then
          do j = 1, n
             i_item = t % material_bins(j) % scalar
             call add_map_element(tally_maps(T_MATERIAL) % items(i_item), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for universe bins
       n = t % n_bins(T_UNIVERSE)
       t % stride(T_UNIVERSE) = filter_bins
       if (n > 0) then
          do j = 1, n
             i_item = t % universe_bins(j) % scalar
             call add_map_element(tally_maps(T_UNIVERSE) % items(i_item), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Finally add scoring bins for the macro tallies
       n = t % n_macro_bins
       if (n > 0) then
          score_bins = n
       else
          message = "Must have macro tally bins!"
          call fatal_error()
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

  subroutine score_tally(p, scattered, fissioned)

    type(Particle), pointer :: p
    logical, intent(in)     :: scattered
    logical, intent(in)     :: fissioned

    integer :: i
    integer :: j
    integer :: k
    integer :: n
    integer :: bins(TALLY_TYPES)
    integer :: bin_energyout
    integer :: score_index
    integer :: score_index0
    integer :: macro_bin
    integer :: mesh_bin
    real(8) :: score
    real(8) :: last_wgt
    real(8) :: wgt
    real(8) :: mu
    real(8) :: E_out
    logical :: has_energyout_bin
    logical :: analog
    type(TallyObject),    pointer :: t
    type(StructuredMesh), pointer :: m

    ! Copy particle's pre- and post-collision weight and angle
    last_wgt = p % last_wgt
    wgt = p % wgt
    mu = p % mu

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    do i = 1, n_tallies
       t => tallies(i)

       ! Surface current tallies are treated separately
       if (t % surface_current) cycle

       ! =======================================================================
       ! DETERMINE SCORING BIN COMBINATION

       ! determine mesh bin
       if (t % n_bins(T_MESH) > 0) then
          m => meshes(t % mesh)

          ! Determine if we're in the mesh first
          call get_mesh_bin(m, p % coord0 % xyz, mesh_bin)
          if (mesh_bin == NO_BIN_FOUND) cycle
          bins(T_MESH) = mesh_bin
       else
          bins(T_MESH) = 1
       end if

       ! determine next universe bin
       if (t % n_bins(T_UNIVERSE) > 0) then
          bins(T_UNIVERSE) = get_next_bin(T_UNIVERSE, p % coord % universe, i)
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
          bins(T_CELL) = get_next_bin(T_CELL, p % coord % cell, i)
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
          has_energyout_bin = .true.
       else
          bins(T_ENERGYOUT) = 1
          has_energyout_bin = .false.
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
          ! determine what type of macro bin
          macro_bin = t % macro_bins(j) % scalar

          ! determine if we need outgoing angle
          analog = (macro_bin == MACRO_NU_SCATTER .or. &
               macro_bin == MACRO_SCATTER_1 .or. macro_bin == MACRO_SCATTER_2 .or. &
               macro_bin == MACRO_SCATTER_3 .or. macro_bin == MACRO_N_1N .or. &
               macro_bin == MACRO_N_2N .or. macro_bin == MACRO_N_3N .or. &
               macro_bin == MACRO_N_4N)

          if (has_energyout_bin .or. analog) then
             ! If this tally has an outgoing energy filter, the only supported
             ! reaction is scattering or nu-fission. Note that some macro
             ! quantities can only be scored in analog such as scattering
             ! moments and (n,xn)

             if ((.not. scattered) .and. (.not. fissioned)) cycle

             if (scattered) then
                ! since scattering has already occured, we do not need to
                ! multiply by the scattering cross section

                select case (macro_bin)
                case (MACRO_SCATTER)
                   score = last_wgt
                case (MACRO_NU_SCATTER)
                   score = wgt
                case (MACRO_SCATTER_1)
                   score = last_wgt * mu
                case (MACRO_SCATTER_2)
                   score = last_wgt * 0.5*(3.0*mu*mu - ONE)
                case (MACRO_SCATTER_3)
                   score = last_wgt * 0.5*(5.0*mu*mu*mu - 3.0*mu)
                case (MACRO_N_1N)
                   if (wgt == last_wgt) then
                      score = last_wgt
                   else
                      cycle
                   end if
                case (MACRO_N_2N)
                   if (int(wgt/last_wgt) == 2) then
                      score = last_wgt
                   else
                      cycle
                   end if
                case (MACRO_N_3N)
                   if (int(wgt/last_wgt) == 3) then
                      score = last_wgt
                   else
                      cycle
                   end if
                case (MACRO_N_4N)
                   if (int(wgt/last_wgt) == 4) then
                      score = last_wgt
                   else
                      cycle
                   end if
                case (MACRO_NU_FISSION)
                   cycle
                case default
                   message = "Invalid macro reaction on analog tally."
                   call fatal_error()
                end select
             end if

             if (fissioned) then

                if (macro_bin /= MACRO_NU_FISSION) cycle

                ! Normally, we only need to make contributions to one scoring
                ! bin. However, in the case of fission, since multiple fission
                ! neutrons were emitted with different energies, multiple
                ! outgoing energy bins may have been scored to. The following
                ! logic treats this special case and scores to multiple bins

                ! save original outgoing energy bin and score index
                bin_energyout = bins(T_ENERGYOUT)
                score_index0 = score_index

                ! Since the creation of fission sites is weighted such that it
                ! is expected to create n_particles sites, we need to multiply
                ! the score by keff to get the true nu-fission rate. Otherwise,
                ! the sum of all nu-fission rates would be ~1.0.

                score = keff

                ! loop over number of particles banked
                do k = 1, p % n_bank
                   ! determine outgoing energy from fission bank
                   E_out = fission_bank(n_bank - p % n_bank + k) % E

                   ! change outgoing energy bin
                   bins(T_ENERGYOUT) = binary_search(t % energy_out, n + 1, E_out)

                   ! determine scoring index
                   score_index = sum((bins - 1) * t % stride) + 1

                   ! Add score to tally
                   call add_to_score(t % scores(score_index, j), score)
                end do

                ! reset outgoing energy bin and score index
                bins(T_ENERGYOUT) = bin_energyout 
                score_index = score_index0
                cycle

             end if
          else
             ! For tallies with no outgoing energy filter, the score is
             ! calculated normally depending on the quantity specified

             select case (macro_bin)
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
! SCORE_SURFACE_CURRENT tallies surface crossings in a mesh tally by manually
! determining which mesh surfaces were crossed
!===============================================================================

  subroutine score_surface_current(p)

    type(Particle),    pointer :: p

    integer :: i                 ! loop indices
    integer :: j                 ! loop indices
    integer :: k                 ! loop indices
    integer :: ijk0(3)           ! indices of starting coordinates
    integer :: ijk1(3)           ! indices of ending coordinates
    integer :: n_cross           ! number of surface crossings
    integer :: n                 ! number of incoming energy bins
    integer :: bins(TALLY_TYPES) ! scoring bin combination
    integer :: score_index       ! index of scoring bin
    real(8) :: uvw(3)            ! cosine of angle of particle
    real(8) :: xyz0(3)           ! starting/intermediate coordinates
    real(8) :: xyz1(3)           ! ending coordinates of particle
    real(8) :: xyz_cross(3)      ! coordinates of bounding surfaces
    real(8) :: d(3)              ! distance to each bounding surface
    real(8) :: distance          ! actual distance traveled
    logical :: start_in_mesh     ! particle's starting xyz in mesh?
    logical :: end_in_mesh       ! particle's ending xyz in mesh?
    logical :: x_same            ! same starting/ending x index (i)
    logical :: y_same            ! same starting/ending y index (j)
    logical :: z_same            ! same starting/ending z index (k)
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()

    bins = 1

    do i = 1, n_tallies
       ! Copy starting and ending location of particle
       xyz0 = p % last_xyz
       xyz1 = p % coord0 % xyz

       ! Get pointer to tally
       t => tallies(i)

       ! Skip non-surface-current tallies
       if (.not. t % surface_current) cycle

       ! Determine indices for starting and ending location
       m => meshes(t % mesh)
       call get_mesh_indices(m, xyz0, ijk0, start_in_mesh)
       call get_mesh_indices(m, xyz1, ijk1, end_in_mesh)

       ! Check to make sure start or end is in mesh
       if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) cycle

       ! Calculate number of surface crossings
       n_cross = sum(abs(ijk1 - ijk0))
       if (n_cross == 0) cycle

       ! Copy particle's direction
       uvw = p % coord0 % uvw

       ! determine incoming energy bin
       n = t % n_bins(T_ENERGYIN)
       if (n > 0) then
          ! check if energy of the particle is within energy bins
          if (p % last_E < t % energy_in(1) .or. &
               p % last_E > t % energy_in(n + 1)) cycle

          ! search to find incoming energy bin
          bins(TS_ENERGYIN) = binary_search(t % energy_in, n + 1, p % last_E)
       else
          bins(TS_ENERGYIN) = 1
       end if

       ! =======================================================================
       ! SPECIAL CASES WHERE TWO INDICES ARE THE SAME

       x_same = (ijk0(1) == ijk1(1))
       y_same = (ijk0(2) == ijk1(2))
       z_same = (ijk0(3) == ijk1(3))

       if (x_same .and. y_same) then
          ! Only z crossings
          if (uvw(3) > 0) then
             do j = ijk0(3), ijk1(3) - 1
                ijk0(3) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = OUT_TOP
                   bins(1:3) = ijk0 + 1
                   score_index = sum((bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(score_index, 1), p % wgt)
                end if
             end do
          else
             do j = ijk0(3) - 1, ijk1(3), -1
                ijk0(3) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = IN_TOP
                   bins(1:3) = ijk0 + 1
                   score_index = sum((bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(score_index, 1), p % wgt)
                end if
             end do
          end if
          cycle
       elseif (x_same .and. z_same) then
          ! Only y crossings
          if (uvw(2) > 0) then
             do j = ijk0(2), ijk1(2) - 1
                ijk0(2) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = OUT_FRONT
                   bins(1:3) = ijk0 + 1
                   score_index = sum((bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(score_index, 1), p % wgt)
                end if
             end do
          else
             do j = ijk0(2) - 1, ijk1(2), -1
                ijk0(2) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = IN_FRONT
                   bins(1:3) = ijk0 + 1
                   score_index = sum((bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(score_index, 1), p % wgt)
                end if
             end do
          end if
          cycle
       elseif (y_same .and. z_same) then
          ! Only x crossings
          if (uvw(1) > 0) then
             do j = ijk0(1), ijk1(1) - 1
                ijk0(1) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = OUT_RIGHT
                   bins(1:3) = ijk0 + 1
                   score_index = sum((bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(score_index, 1), p % wgt)
                end if
             end do
          else
             do j = ijk0(1) - 1, ijk1(1), -1
                ijk0(1) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = IN_RIGHT
                   bins(1:3) = ijk0 + 1
                   score_index = sum((bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(score_index, 1), p % wgt)
                end if
             end do
          end if
          cycle
       end if

       ! =======================================================================
       ! GENERIC CASE

       ! Bounding coordinates
       do j = 1, 3
          if (uvw(j) > 0) then
             xyz_cross(j) = m % origin(j) + ijk0(j) * m % width(j)
          else
             xyz_cross(j) = m % origin(j) + (ijk0(j) - 1) * m % width(j)
          end if
       end do

       do k = 1, n_cross
          ! Reset scoring bin index
          bins(TS_SURFACE) = 0

          ! Calculate distance to each bounding surface. We need to treat
          ! special case where the cosine of the angle is zero since this would
          ! result in a divide-by-zero.
          
          do j = 1, 3
             if (uvw(j) == 0) then
                d(j) = INFINITY
             else
                d(j) = (xyz_cross(j) - xyz0(j))/uvw(j)
             end if
          end do

          ! Determine the closest bounding surface of the mesh cell by
          ! calculating the minimum distance
          
          distance = minval(d)

          ! Now use the minimum distance and diretion of the particle to
          ! determine which surface was crossed
          
          if (distance == d(1)) then
             if (uvw(1) > 0) then
                ! Crossing into right mesh cell -- this is treated as outgoing
                ! current from (i,j,k)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = OUT_RIGHT
                   bins(1:3) = ijk0 + 1
                end if
                ijk0(1) = ijk0(1) + 1
                xyz_cross(1) = xyz_cross(1) + m % width(1)
             else
                ! Crossing into left mesh cell -- this is treated as incoming
                ! current in (i-1,j,k)
                ijk0(1) = ijk0(1) - 1
                xyz_cross(1) = xyz_cross(1) - m % width(1)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = IN_RIGHT
                   bins(1:3) = ijk0 + 1
                end if
             end if
          elseif (distance == d(2)) then
             if (uvw(2) > 0) then
                ! Crossing into front mesh cell -- this is treated as outgoing
                ! current in (i,j,k)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = OUT_FRONT
                   bins(1:3) = ijk0 + 1
                end if
                ijk0(2) = ijk0(2) + 1
                xyz_cross(2) = xyz_cross(2) + m % width(2)
             else
                ! Crossing into back mesh cell -- this is treated as incoming
                ! current in (i,j-1,k)
                ijk0(2) = ijk0(2) - 1
                xyz_cross(2) = xyz_cross(2) - m % width(2)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = IN_FRONT
                   bins(1:3) = ijk0 + 1
                end if
             end if
          else if (distance == d(3)) then
             if (uvw(3) > 0) then
                ! Crossing into top mesh cell -- this is treated as outgoing
                ! current in (i,j,k)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = OUT_TOP
                   bins(1:3) = ijk0 + 1
                end if
                ijk0(3) = ijk0(3) + 1
                xyz_cross(3) = xyz_cross(3) + m % width(3)
             else
                ! Crossing into bottom mesh cell -- this is treated as incoming
                ! current in (i,j,k-1)
                ijk0(3) = ijk0(3) - 1
                xyz_cross(3) = xyz_cross(3) - m % width(3)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(TS_SURFACE) = IN_TOP
                   bins(1:3) = ijk0 + 1
                end if
             end if
          end if

          ! Determine scoring index
          if (bins(TS_SURFACE) > 0) then
             score_index = sum((bins - 1) * t % stride) + 1

             ! Check for errors
             if (score_index <= 0 .or. score_index > t % n_total_bins) then
                message = "Score index outside range."
                call fatal_error()
             end if

             ! Add to surface current tally
             call add_to_score(t % scores(score_index, 1), p % wgt)
          end if

          ! Calculate new coordinates
          xyz0 = xyz0 + distance * uvw
       end do

    end do

  end subroutine score_surface_current

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
    character(MAX_FILE_LEN) :: filename                  ! name of output file
    character(15)           :: filter_name(TALLY_TYPES)  ! names of tally filters
    character(27)           :: macro_name(N_MACRO_TYPES) ! names of macro scores
    type(TallyObject), pointer :: t

    ! Skip if there are no tallies
    if (n_tallies == 0) return

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
    macro_name(abs(MACRO_NU_SCATTER)) = "Scattering Production Rate"
    macro_name(abs(MACRO_SCATTER_1))  = "First Scattering Moment"
    macro_name(abs(MACRO_SCATTER_2))  = "Second Scattering Moment"
    macro_name(abs(MACRO_SCATTER_3))  = "Third Scattering Moment"
    macro_name(abs(MACRO_N_1N))       = "(n,1n) Rate"
    macro_name(abs(MACRO_N_2N))       = "(n,2n) Rate"
    macro_name(abs(MACRO_N_3N))       = "(n,3n) Rate"
    macro_name(abs(MACRO_N_4N))       = "(n,4n) Rate"
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
       call header("TALLY " // trim(int_to_str(t % id)), unit=UNIT_TALLY, level=3)

       ! Handle surface current tallies separately
       if (t % surface_current) then
          call write_surface_current(t)
          cycle
       end if

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
                      write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                           trim(filter_name(j)), trim(get_label(t, j, bins(j)))
                      indent = indent + 2
                   end if
                   j = j + 1
                end if
             end if

          end do find_bin

          ! Print filter information
          write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
               trim(filter_name(j)), trim(get_label(t, j, bins(j)))

          ! Determine scoring index for this bin combination -- note that unlike
          ! in the score_tally subroutine, we have to use max(bins,1) since all
          ! bins below the lowest filter level will be zeros

          score_index = sum((max(bins,1) - 1) * t % stride) + 1

          ! Write scores for this filter bin combination
          indent = indent + 2
          do k = 1, t % n_macro_bins
             write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') & 
                  repeat(" ", indent), macro_name(abs(t % macro_bins(k) % scalar)), &
                  real_to_str(t % scores(score_index,k) % val), &
                  trim(real_to_str(t % scores(score_index,k) % val_sq))
          end do
          indent = indent - 2

       end do print_bin

    end do

    close(UNIT=UNIT_TALLY)

  end subroutine write_tallies

!===============================================================================
! WRITE_SURFACE_CURRENT writes out surface current tallies over a mesh to the
! tallies.out file.
!===============================================================================

  subroutine write_surface_current(t)

    type(TallyObject), pointer :: t

    integer :: i                 ! mesh index for x
    integer :: j                 ! mesh index for y
    integer :: k                 ! mesh index for z
    integer :: l                 ! mesh index for energy
    integer :: bins(TALLY_TYPES) ! bin combination
    integer :: n                 ! number of incoming energy bins
    integer :: len1              ! length of string 
    integer :: len2              ! length of string 
    integer :: score_index       ! index in scores array for filters
    logical :: print_ebin        ! should incoming energy bin be displayed?
    character(MAX_LINE_LEN) :: string
    type(StructuredMesh), pointer :: m => null()

    ! Get pointer to mesh
    m => meshes(t % mesh)

    ! initialize bins array
    bins = 1

    ! determine how many energy in bins there are
    n = t % n_bins(T_ENERGYIN)
    if (n > 0) then
       print_ebin = .true.
    else
       print_ebin = .false.
       n = 1
    end if

    do i = 1, m % dimension(1)
       string = "Mesh Index (" // trim(int_to_str(i)) // ", "
       len1 = len_trim(string)
       do j = 1, m % dimension(2)
          string = string(1:len1+1) // trim(int_to_str(j)) // ", "
          len2 = len_trim(string)
          do k = 1, m % dimension(3)
             ! Write mesh cell index
             string = string(1:len2+1) // trim(int_to_str(k)) // ")"
             write(UNIT=UNIT_TALLY, FMT='(1X,A)') trim(string)

             do l = 1, n
                ! Write incoming energy bin
                if (print_ebin) then
                   write(UNIT=UNIT_TALLY, FMT='(3X,A,1X,A)') &
                        "Incoming Energy", trim(get_label(t, T_ENERGYIN, l))
                end if

                ! Set incoming energy bin
                bins(TS_ENERGYIN) = l

                ! Left Surface
                bins(1:3) = (/ i-1, j, k /) + 1
                bins(TS_SURFACE) = IN_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Left", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                bins(TS_SURFACE) = OUT_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Left", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                ! Right Surface
                bins(1:3) = (/ i, j, k /) + 1
                bins(TS_SURFACE) = IN_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Right", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                bins(TS_SURFACE) = OUT_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Right", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                ! Back Surface
                bins(1:3) = (/ i, j-1, k /) + 1
                bins(TS_SURFACE) = IN_FRONT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Back", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                bins(TS_SURFACE) = OUT_FRONT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Back", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                ! Front Surface
                bins(1:3) = (/ i, j, k /) + 1
                bins(TS_SURFACE) = IN_FRONT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Front", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                bins(TS_SURFACE) = OUT_FRONT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Front", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                ! Bottom Surface
                bins(1:3) = (/ i, j, k-1 /) + 1
                bins(TS_SURFACE) = IN_TOP
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Bottom", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                bins(TS_SURFACE) = OUT_TOP
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Bottom", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                ! Top Surface
                bins(1:3) = (/ i, j, k /) + 1
                bins(TS_SURFACE) = IN_TOP
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Top", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))

                bins(TS_SURFACE) = OUT_TOP
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Top", &
                     real_to_str(t % scores(score_index,1) % val), &
                     trim(real_to_str(t % scores(score_index,1) % val_sq))
             end do

          end do
       end do
    end do

  end subroutine write_surface_current

!===============================================================================
! GET_LABEL returns a label for a cell/surface/etc given a tally, filter type,
! and corresponding bin
!===============================================================================

  function get_label(t, filter_type, bin) result(label)

    type(TallyObject), pointer :: t           ! tally object
    integer, intent(in)        :: filter_type ! type of filter
    integer, intent(in)        :: bin         ! bin in filter array
    character(30)              :: label         ! user-specified identifier

    integer              :: i      ! index in cells/surfaces/etc array
    integer, allocatable :: ijk(:) ! indices in mesh
    real(8)              :: E0     ! lower bound for energy bin
    real(8)              :: E1     ! upper bound for energy bin
    type(StructuredMesh), pointer :: m

    select case(filter_type)
    case (T_UNIVERSE)
       i = t % universe_bins(bin) % scalar
       label = int_to_str(universes(i) % id)
    case (T_MATERIAL)
       i = t % material_bins(bin) % scalar
       label = int_to_str(materials(i) % id)
    case (T_CELL)
       i = t % cell_bins(bin) % scalar
       label = int_to_str(cells(i) % id)
    case (T_CELLBORN)
       i = t % cellborn_bins(bin) % scalar
       label = int_to_str(cells(i) % id)
    case (T_SURFACE)
       i = t % surface_bins(bin) % scalar
       label = int_to_str(surfaces(i) % id)
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
