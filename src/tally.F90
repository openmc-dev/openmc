module tally

  use constants
  use error,         only: fatal_error
  use global
  use mesh,          only: get_mesh_bin, bin_to_mesh_indices, get_mesh_indices
  use mesh_header,   only: StructuredMesh
  use output,        only: header
  use search,        only: binary_search
  use string,        only: to_str
  use tally_header,  only: TallyScore, TallyMapItem, TallyMapElement

#ifdef MPI
  use mpi
#endif

  implicit none

  integer, allocatable :: position(:)

contains

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
    allocate(tally_maps(N_FILTER_TYPES - 3))

    ! allocate list of items for each different filter type
    allocate(tally_maps(FILTER_UNIVERSE) % items(n_universes))
    allocate(tally_maps(FILTER_MATERIAL) % items(n_materials))
    allocate(tally_maps(FILTER_CELL)     % items(n_cells))
    allocate(tally_maps(FILTER_CELLBORN) % items(n_cells))
    allocate(tally_maps(FILTER_SURFACE)  % items(n_surfaces))

    ! Allocate and initialize tally map positioning for finding bins
    allocate(position(N_FILTER_TYPES - 3))
    position = 0

    do i = 1, n_tallies
       t => tallies(i)

       ! initialize number of filter bins
       filter_bins = 1

       ! handle surface current tallies separately
       if (t % type == TALLY_SURFACE_CURRENT) then
          m => meshes(t % mesh)

          t % stride(SURF_FILTER_SURFACE) = filter_bins
          ! Set stride for surface/direction
          if (m % n_dimension == 2) then
             filter_bins = filter_bins * 4
          elseif (m % n_dimension == 3) then
             filter_bins = filter_bins * 6
          end if

          ! Add filter for incoming energy
          n = t % n_filter_bins(FILTER_ENERGYIN)
          t % stride(SURF_FILTER_ENERGYIN) = filter_bins
          if (n > 0) then
             filter_bins = filter_bins * n
          end if
          
          ! account for z direction
          t % stride(SURF_FILTER_MESH_Z) = filter_bins
          filter_bins = filter_bins * (m % dimension(3) + 1)

          ! account for y direction
          t % stride(SURF_FILTER_MESH_Y) = filter_bins
          filter_bins = filter_bins * (m % dimension(2) + 1)

          ! account for z direction
          t % stride(SURF_FILTER_MESH_X) = filter_bins
          filter_bins = filter_bins * (m % dimension(1) + 1)

          ! Finally add scoring bins for the tallies and allocate the main
          ! scores array for tally
          score_bins = t % n_score_bins
          t % n_total_bins = filter_bins
          allocate(t % scores(filter_bins, score_bins))

          ! All the remaining logic is for non-surface-current tallies so we can
          ! safely skip it
          cycle
       end if

       ! determine if there are subdivisions for incoming or outgoing energy to
       ! adjust the number of filter bins appropriately
       n = t % n_filter_bins(FILTER_ENERGYOUT)
       t % stride(FILTER_ENERGYOUT) = filter_bins
       if (n > 0) then
          filter_bins = filter_bins * n
       end if

       n = t % n_filter_bins(FILTER_ENERGYIN)
       t % stride(FILTER_ENERGYIN) = filter_bins
       if (n > 0) then
          filter_bins = filter_bins * n
       end if

       ! Add map elements for mesh bins
       n = t % n_filter_bins(FILTER_MESH)
       t % stride(FILTER_MESH) = filter_bins
       if (n > 0) then
          filter_bins = filter_bins * n
       end if

       ! Add map elements for surface bins
       n = t % n_filter_bins(FILTER_SURFACE)
       t % stride(FILTER_SURFACE) = filter_bins
       if (n > 0) then
          do j = 1, n
             i_item = t % surface_bins(j) % scalar
             call add_map_element(tally_maps(FILTER_SURFACE) % items(i_item), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for cellborn bins
       n = t % n_filter_bins(FILTER_CELLBORN)
       t % stride(FILTER_CELLBORN) = filter_bins
       if (n > 0) then
          do j = 1, n
             i_item = t % cellborn_bins(j) % scalar
             call add_map_element(tally_maps(FILTER_CELLBORN) % items(i_item), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for cell bins
       n = t % n_filter_bins(FILTER_CELL)
       t % stride(FILTER_CELL) = filter_bins
       if (n > 0) then
          do j = 1, n
             i_item = t % cell_bins(j) % scalar
             call add_map_element(tally_maps(FILTER_CELL) % items(i_item), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for material bins
       n = t % n_filter_bins(FILTER_MATERIAL)
       t % stride(FILTER_MATERIAL) = filter_bins
       if (n > 0) then
          do j = 1, n
             i_item = t % material_bins(j) % scalar
             call add_map_element(tally_maps(FILTER_MATERIAL) % items(i_item), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Add map elements for universe bins
       n = t % n_filter_bins(FILTER_UNIVERSE)
       t % stride(FILTER_UNIVERSE) = filter_bins
       if (n > 0) then
          do j = 1, n
             i_item = t % universe_bins(j) % scalar
             call add_map_element(tally_maps(FILTER_UNIVERSE) % items(i_item), i, j)
          end do
          filter_bins = filter_bins * n
       end if

       ! Finally add scoring bins for the tally
       n = t % n_score_bins
       if (n > 0) then
          score_bins = n
       else
          message = "Must have scoring bins!"
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
! SCORE_ANALOG_TALLY keeps track of how many events occur in a specified cell,
! energy range, etc. Note that since these are "analog" tallies, they are only
! triggered at every collision, not every event
!===============================================================================

  subroutine score_analog_tally()

    integer :: i                    ! loop index for analog tallies
    integer :: j                    ! loop index for scoring bins
    integer :: bins(N_FILTER_TYPES) ! scoring bin combination
    integer :: score_index          ! single index for single bin
    integer :: score_bin            ! scoring bin, e.g. SCORE_FLUX
    real(8) :: score                ! analog tally score
    real(8) :: last_wgt             ! pre-collision particle weight
    real(8) :: wgt                  ! post-collision particle weight
    real(8) :: mu                   ! cosine of angle of collision
    logical :: found_bin            ! scoring bin found?
    type(TallyObject), pointer :: t => null()

    ! Copy particle's pre- and post-collision weight and angle
    last_wgt = p % last_wgt
    wgt = p % wgt
    mu = p % mu

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    do i = 1, n_analog_tallies
       t => tallies(analog_tallies(i))

       ! =======================================================================
       ! DETERMINE SCORING BIN COMBINATION

       call get_scoring_bins(analog_tallies(i), bins, found_bin)
       if (.not. found_bin) cycle

       ! =======================================================================
       ! CALCULATE SCORES AND ACCUMULATE TALLY

       ! If we have made it here, we have a scoring combination of bins for this
       ! tally -- now we need to determine where in the scores array we should
       ! be accumulating the tally values

       ! Determine scoring index for this filter combination
       score_index = sum((bins - 1) * t % stride) + 1

       ! Determine score for each bin
       do j = 1, t % n_score_bins
          ! determine what type of score bin
          score_bin = t % score_bins(j) % scalar

          select case (score_bin)
          case (SCORE_FLUX)
             ! All events score to a flux bin. We actually use a collision
             ! estimator since there is no way to count 'events' exactly for the
             ! flux

             score = last_wgt / material_xs % total

          case (SCORE_TOTAL)
             ! All events will score to the total reaction rate. We can just use
             ! the weight of the particle entering the collision as the score

             score = last_wgt

          case (SCORE_SCATTER)
             ! Skip any event where the particle didn't scatter
             if (p % event /= EVENT_SCATTER) cycle

             ! Since only scattering events make it here, again we can use the
             ! weight entering the collision as the estimator for the reaction
             ! rate

             score = last_wgt

          case (SCORE_NU_SCATTER) 
             ! Skip any event where the particle didn't scatter
             if (p % event /= EVENT_SCATTER) cycle

             ! For scattering production, we need to use the post-collision
             ! weight as the estimate for the number of neutrons exiting a
             ! reaction with neutrons in the exit channel

             score = wgt

          case (SCORE_SCATTER_1)
             ! Skip any event where the particle didn't scatter
             if (p % event /= EVENT_SCATTER) cycle

             ! The first scattering moment can be determined by using the rate
             ! of scattering reactions multiplied by the cosine of the change in
             ! neutron's angle due to the collision

             score = last_wgt * mu

          case (SCORE_SCATTER_2)
             ! Skip any event where the particle didn't scatter
             if (p % event /= EVENT_SCATTER) cycle

             ! The second scattering moment can be determined in a similar
             ! manner to the first scattering moment

             score = last_wgt * 0.5*(3.0*mu*mu - ONE)

          case (SCORE_SCATTER_3)
             ! Skip any event where the particle didn't scatter
             if (p % event /= EVENT_SCATTER) cycle

             ! The first scattering moment can be determined by using the rate
             ! of scattering reactions multiplied by the cosine of the change in
             ! neutron's angle due to the collision

             score = last_wgt * 0.5*(5.0*mu*mu*mu - 3.0*mu)

          case (SCORE_DIFFUSION)
             ! Skip any event where the particle didn't scatter
             if (p % event /= EVENT_SCATTER) cycle
             
             ! Temporarily store the scattering cross section
             score = material_xs % total - material_xs % absorption

             ! Since this only gets tallied at every scattering event, the flux
             ! estimator is 1/Sigma_s. Therefore, the diffusion coefficient
             ! times flux is 1/(3*Sigma_s*(Sigma_t - mu*Sigma_s)).

             score = last_wgt / (3.0_8 * score * (material_xs % total - &
                  mu * score))

          case (SCORE_N_1N)
             ! Skip any event where the particle didn't scatter
             if (p % event /= EVENT_SCATTER) cycle

             ! Skip any events where weight of particle changed
             if (wgt /= last_wgt) cycle

             ! All events that reach this point are (n,1n) reactions
             score = last_wgt

          case (SCORE_N_2N)
             ! Skip any event where the particle didn't scatter
             if (p % event /= EVENT_SCATTER) cycle

             ! Skip any scattering events where the weight didn't double
             if (nint(wgt/last_wgt) /= 2) cycle

             ! All events that reach this point are (n,2n) reactions
             score = last_wgt

          case (SCORE_N_3N)
             ! Skip any event where the particle didn't scatter
             if (p % event /= EVENT_SCATTER) cycle

             ! Skip any scattering events where the weight didn't double
             if (nint(wgt/last_wgt) /= 3) cycle

             ! All events that reach this point are (n,3n) reactions
             score = last_wgt

          case (SCORE_N_4N)
             ! Skip any event where the particle didn't scatter
             if (p % event /= EVENT_SCATTER) cycle

             ! Skip any scattering events where the weight didn't double
             if (nint(wgt/last_wgt) /= 4) cycle

             ! All events that reach this point are (n,3n) reactions
             score = last_wgt

          case (SCORE_ABSORPTION)
             ! Skip any event where the particle wasn't absorbed
             if (p % event == EVENT_SCATTER) cycle

             ! All fission and absorption events will contribute here, so we can
             ! just use the particle's weight entering the collision

             score = last_wgt

          case (SCORE_FISSION)
             ! Skip any non-fission events
             if (p % event /= EVENT_FISSION) cycle

             ! All fission events will contribute, so again we can use
             ! particle's weight entering the collision as the estimate for the
             ! fission reaction rate

             score = last_wgt

          case (SCORE_NU_FISSION)
             ! Skip any non-fission events
             if (p % event /= EVENT_FISSION) cycle

             if (t % n_filter_bins(FILTER_ENERGYOUT) > 0) then
                ! Normally, we only need to make contributions to one scoring
                ! bin. However, in the case of fission, since multiple fission
                ! neutrons were emitted with different energies, multiple outgoing
                ! energy bins may have been scored to. The following logic treats
                ! this special case and scores to multiple bins

                call score_fission_eout(t, bins, j)
                cycle

             else
                ! If there is no outgoing energy filter, than we only need to
                ! score to one bin. For the score to be 'analog', we need to
                ! score the number of particles that were banked in the fission
                ! bank. Since this was weighted by 1/keff, we multiply by keff
                ! to get the proper score.

                score = keff * p % n_bank

             end if

          case default
             message = "Invalid score type on tally " // to_str(t % id) // "."
             call fatal_error()
          end select
             
          ! Add score to tally
          call add_to_score(t % scores(score_index, j), score)

       end do

       ! If the user has specified that we can assume all tallies are spatially
       ! separate, this implies that once a tally has been scored to, we needn't
       ! check the others. This cuts down on overhead when there are many
       ! tallies specified

       if (assume_separate) return

       ! Reset tally map positioning
       position = 0

    end do

  end subroutine score_analog_tally

!===============================================================================
! SCORE_FISSION_EOUT handles a special case where we need to store neutron
! production rate with an outgoing energy filter (think of a fission matrix). In
! this case, we may need to score to multiple bins if there were multiple
! neutrons produced with different energies.
!===============================================================================

  subroutine score_fission_eout(t, bins, j)

    type(TallyObject), pointer :: t
    integer, intent(inout)     :: bins(N_FILTER_TYPES)
    integer, intent(in)        :: j ! index for score

    integer :: k             ! loop index for bank sites
    integer :: bin_energyout ! original outgoing energy bin
    integer :: score_index   ! index for scoring bin combination
    real(8) :: score         ! actualy score
    real(8) :: E_out         ! energy of fission bank site

    ! save original outgoing energy bin and score index
    bin_energyout = bins(FILTER_ENERGYOUT)

    ! Since the creation of fission sites is weighted such that it is
    ! expected to create n_particles sites, we need to multiply the
    ! score by keff to get the true nu-fission rate. Otherwise, the sum
    ! of all nu-fission rates would be ~1.0.

    score = keff

    ! loop over number of particles banked
    do k = 1, p % n_bank
       ! determine outgoing energy from fission bank
       E_out = fission_bank(n_bank - p % n_bank + k) % E

       ! change outgoing energy bin
       bins(FILTER_ENERGYOUT) = binary_search(t % energy_out, &
            size(t % energy_out), E_out)

       ! determine scoring index
       score_index = sum((bins - 1) * t % stride) + 1

       ! Add score to tally
       call add_to_score(t % scores(score_index, j), score)
    end do

    ! reset outgoing energy bin and score index
    bins(FILTER_ENERGYOUT) = bin_energyout

  end subroutine score_fission_eout

!===============================================================================
! SCORE_TRACKLENGTH_TALLY calculates fluxes and reaction rates based on the
! track-length estimate of the flux. This is triggered at every event (surface
! crossing, lattice crossing, or collision) and thus cannot be done for tallies
! that require post-collision information.
!===============================================================================

  subroutine score_tracklength_tally(distance)

    real(8), intent(in) :: distance

    integer :: i                    ! loop index for tracklength tallies
    integer :: j                    ! loop index for scoring bins
    integer :: bins(N_FILTER_TYPES) ! scoring bin combination
    integer :: score_index          ! single index for single bin
    integer :: score_bin            ! scoring bin, e.g. SCORE_FLUX
    real(8) :: flux                 ! tracklength estimate of flux
    real(8) :: score                ! actual score (e.g., flux*xs)
    logical :: found_bin            ! scoring bin found?
    type(TallyObject), pointer :: t => null()

    ! Determine track-length estimate of flux
    flux = p % wgt * distance

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    do i = 1, n_tracklength_tallies
       t => tallies(tracklength_tallies(i))

       ! =======================================================================
       ! DETERMINE SCORING BIN COMBINATION

       call get_scoring_bins(tracklength_tallies(i), bins, found_bin)
       if (.not. found_bin) cycle

       ! =======================================================================
       ! CALCULATE SCORES AND ACCUMULATE TALLY

       ! If we have made it here, we have a scoring combination of bins for this
       ! tally -- now we need to determine where in the scores array we should
       ! be accumulating the tally values

       ! Determine scoring index for this filter combination
       score_index = sum((bins - 1) * t % stride) + 1

       ! Determine score for each bin
       do j = 1, t % n_score_bins
          ! determine what type of score bin
          score_bin = t % score_bins(j) % scalar

          ! Determine cross section 
          select case(score_bin)
          case (SCORE_FLUX)
             score = flux
          case (SCORE_TOTAL)
             score = material_xs % total * flux
          case (SCORE_SCATTER)
             score = (material_xs % total - material_xs % absorption) * flux
          case (SCORE_ABSORPTION)
             score = material_xs % absorption * flux
          case (SCORE_FISSION)
             score = material_xs % fission * flux
          case (SCORE_NU_FISSION)
             score = material_xs % nu_fission * flux
          case default
             message = "Invalid score type on tally " // to_str(t % id) // "."
             call fatal_error()
          end select

          ! Add score to tally
          call add_to_score(t % scores(score_index, j), score)

       end do

       ! If the user has specified that we can assume all tallies are spatially
       ! separate, this implies that once a tally has been scored to, we needn't
       ! check the others. This cuts down on overhead when there are many
       ! tallies specified

       if (assume_separate) return

       ! Reset tally map positioning
       position = 0

    end do

  end subroutine score_tracklength_tally

!===============================================================================
! GET_SCORING_BINS determines a combination of filter bins that should be scored
! for a tally based on the particle's current attributes.
!===============================================================================

  subroutine get_scoring_bins(index_tally, bins, found_bin)

    integer, intent(in)     :: index_tally
    integer, intent(out)    :: bins(N_FILTER_TYPES)
    logical, intent(out)    :: found_bin

    integer :: i        ! loop index for filters
    integer :: n        ! number of bins for single filter
    integer :: mesh_bin ! index for mesh bin
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()

    found_bin = .true.
    t => tallies(index_tally)
    bins = 1

    FILTER_LOOP: do i = 1, t % n_filters

       select case (t % filters(i))
       case (FILTER_MESH)
          ! determine mesh bin
          m => meshes(t % mesh)

          ! Determine if we're in the mesh first
          call get_mesh_bin(m, p % coord0 % xyz, mesh_bin)
          if (mesh_bin == NO_BIN_FOUND) then
             found_bin = .false.
             return
          end if
          bins(FILTER_MESH) = mesh_bin

       case (FILTER_UNIVERSE)
          ! determine next universe bin
          ! TODO: Account for multiple universes when performing this filter
          bins(FILTER_UNIVERSE) = get_next_bin(FILTER_UNIVERSE, &
               p % coord % universe, index_tally)
          if (bins(FILTER_UNIVERSE) == NO_BIN_FOUND) then
             found_bin = .false.
             return
          end if

       case (FILTER_MATERIAL)
          bins(FILTER_MATERIAL) = get_next_bin(FILTER_MATERIAL, &
               p % material, index_tally)
          if (bins(FILTER_MATERIAL) == NO_BIN_FOUND) then
             found_bin = .false.
             return
          end if

       case (FILTER_CELL)
          ! determine next cell bin
          ! TODO: Account for cells in multiple levels when performing this filter
          bins(FILTER_CELL) = get_next_bin(FILTER_CELL, &
               p % coord % cell, index_tally)
          if (bins(FILTER_CELL) == NO_BIN_FOUND) then
             found_bin = .false.
             return
          end if

       case (FILTER_CELLBORN)
          ! determine next cellborn bin
          bins(FILTER_CELLBORN) = get_next_bin(FILTER_CELLBORN, &
               p % cell_born, index_tally)
          if (bins(FILTER_CELLBORN) == NO_BIN_FOUND) then
             found_bin = .false.
             return
          end if

       case (FILTER_SURFACE)
          ! determine next surface bin
          bins(FILTER_SURFACE) = get_next_bin(FILTER_SURFACE, &
               p % surface, index_tally)
          if (bins(FILTER_SURFACE) == NO_BIN_FOUND) then
             found_bin = .false.
             return
          end if

       case (FILTER_ENERGYIN)
          ! determine incoming energy bin
          n = t % n_filter_bins(FILTER_ENERGYIN)
          ! check if energy of the particle is within energy bins
          if (p % last_E < t % energy_in(1) .or. &
               p % last_E > t % energy_in(n + 1)) then
             found_bin = .false.
             return
          end if

          ! search to find incoming energy bin
          bins(FILTER_ENERGYIN) = binary_search(t % energy_in, n + 1, p % last_E)

       case (FILTER_ENERGYOUT)
          ! determine outgoing energy bin
          n = t % n_filter_bins(FILTER_ENERGYOUT)
          ! check if energy of the particle is within energy bins
          if (p % E < t % energy_out(1) .or. p % E > t % energy_out(n + 1)) then
             found_bin = .false.
             return
          end if

          ! search to find incoming energy bin
          bins(FILTER_ENERGYOUT) = binary_search(t % energy_out, n + 1, p % E)

       end select

    end do FILTER_LOOP

  end subroutine get_scoring_bins

!===============================================================================
! SCORE_SURFACE_CURRENT tallies surface crossings in a mesh tally by manually
! determining which mesh surfaces were crossed
!===============================================================================

  subroutine score_surface_current()

    integer :: i                    ! loop indices
    integer :: j                    ! loop indices
    integer :: k                    ! loop indices
    integer :: ijk0(3)              ! indices of starting coordinates
    integer :: ijk1(3)              ! indices of ending coordinates
    integer :: n_cross              ! number of surface crossings
    integer :: n                    ! number of incoming energy bins
    integer :: bins(N_FILTER_TYPES) ! scoring bin combination
    integer :: score_index          ! index of scoring bin
    real(8) :: uvw(3)               ! cosine of angle of particle
    real(8) :: xyz0(3)              ! starting/intermediate coordinates
    real(8) :: xyz1(3)              ! ending coordinates of particle
    real(8) :: xyz_cross(3)         ! coordinates of bounding surfaces
    real(8) :: d(3)                 ! distance to each bounding surface
    real(8) :: distance             ! actual distance traveled
    logical :: start_in_mesh        ! particle's starting xyz in mesh?
    logical :: end_in_mesh          ! particle's ending xyz in mesh?
    logical :: x_same               ! same starting/ending x index (i)
    logical :: y_same               ! same starting/ending y index (j)
    logical :: z_same               ! same starting/ending z index (k)
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()

    bins = 1

    do i = 1, n_current_tallies
       ! Copy starting and ending location of particle
       xyz0 = p % last_xyz
       xyz1 = p % coord0 % xyz

       ! Get pointer to tally
       t => tallies(current_tallies(i))

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
       n = t % n_filter_bins(FILTER_ENERGYIN)
       if (n > 0) then
          ! check if energy of the particle is within energy bins
          if (p % last_E < t % energy_in(1) .or. &
               p % last_E > t % energy_in(n + 1)) cycle

          ! search to find incoming energy bin
          bins(SURF_FILTER_ENERGYIN) = binary_search(t % energy_in, n + 1, p % last_E)
       else
          bins(SURF_FILTER_ENERGYIN) = 1
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
                   bins(SURF_FILTER_SURFACE) = OUT_TOP
                   bins(1:3) = ijk0 + 1
                   score_index = sum((bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(score_index, 1), p % wgt)
                end if
             end do
          else
             do j = ijk0(3) - 1, ijk1(3), -1
                ijk0(3) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(SURF_FILTER_SURFACE) = IN_TOP
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
                   bins(SURF_FILTER_SURFACE) = OUT_FRONT
                   bins(1:3) = ijk0 + 1
                   score_index = sum((bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(score_index, 1), p % wgt)
                end if
             end do
          else
             do j = ijk0(2) - 1, ijk1(2), -1
                ijk0(2) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(SURF_FILTER_SURFACE) = IN_FRONT
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
                   bins(SURF_FILTER_SURFACE) = OUT_RIGHT
                   bins(1:3) = ijk0 + 1
                   score_index = sum((bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(score_index, 1), p % wgt)
                end if
             end do
          else
             do j = ijk0(1) - 1, ijk1(1), -1
                ijk0(1) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(SURF_FILTER_SURFACE) = IN_RIGHT
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
             xyz_cross(j) = m % lower_left(j) + ijk0(j) * m % width(j)
          else
             xyz_cross(j) = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
          end if
       end do

       do k = 1, n_cross
          ! Reset scoring bin index
          bins(SURF_FILTER_SURFACE) = 0

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
                   bins(SURF_FILTER_SURFACE) = OUT_RIGHT
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
                   bins(SURF_FILTER_SURFACE) = IN_RIGHT
                   bins(1:3) = ijk0 + 1
                end if
             end if
          elseif (distance == d(2)) then
             if (uvw(2) > 0) then
                ! Crossing into front mesh cell -- this is treated as outgoing
                ! current in (i,j,k)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(SURF_FILTER_SURFACE) = OUT_FRONT
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
                   bins(SURF_FILTER_SURFACE) = IN_FRONT
                   bins(1:3) = ijk0 + 1
                end if
             end if
          else if (distance == d(3)) then
             if (uvw(3) > 0) then
                ! Crossing into top mesh cell -- this is treated as outgoing
                ! current in (i,j,k)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   bins(SURF_FILTER_SURFACE) = OUT_TOP
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
                   bins(SURF_FILTER_SURFACE) = IN_TOP
                   bins(1:3) = ijk0 + 1
                end if
             end if
          end if

          ! Determine scoring index
          if (bins(SURF_FILTER_SURFACE) > 0) then
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
    
    score % n_events = score % n_events + 1
    score % value    = score % value    + val
    
  end subroutine add_to_score

!===============================================================================
! ACCUMULATE_BATCH_ESTIMATE accumulates scores from many histories (or many
! generations) into a single realization of a random variable.
!===============================================================================

  elemental subroutine accumulate_batch_estimate(score)

    type(TallyScore), intent(inout) :: score

    real(8) :: val

    ! Add the sum and square of the sum of contributions from each cycle
    ! within a cycle to the variables sum and sum_sq. This will later allow us
    ! to calculate a variance on the tallies

    val = score % value/(n_particles*gen_per_batch)
    score % sum    = score % sum    + val
    score % sum_sq = score % sum_sq + val*val

    ! Reset the single batch estimate
    score % value = ZERO

  end subroutine accumulate_batch_estimate

!===============================================================================
! SYNCHRONIZE_TALLIES accumulates the sum of the contributions from each history
! within the batch to a new random variable
!===============================================================================

  subroutine synchronize_tallies()

    integer :: i   ! index in tallies array
    type(TallyObject), pointer :: t => null()

#ifdef MPI
    call reduce_tallies()
    if (.not. master) return
#endif

    do i = 1, n_tallies
       t => tallies(i)

       ! Loop over all filter and scoring bins
       call accumulate_batch_estimate(t % scores)
    end do

  end subroutine synchronize_tallies

!===============================================================================
! REDUCE_TALLIES collects all the results from tallies onto one processor
!===============================================================================

#ifdef MPI
  subroutine reduce_tallies()

    integer :: i      ! loop index for tallies
    integer :: n      ! number of filter bins
    integer :: m      ! number of score bins
    integer :: n_bins ! total number of bins
    real(8), allocatable :: tally_temp(:,:) ! contiguous array of scores
    type(TallyObject), pointer :: t => null()

    do i = 1, n_tallies
       t => tallies(i)

       n = t % n_total_bins
       m = t % n_score_bins
       n_bins = n*m

       allocate(tally_temp(n,m))

       tally_temp = t % scores(:,:) % value

       if (master) then
          ! The MPI_IN_PLACE specifier allows the master to copy values into a
          ! receive buffer without having a temporary variable
          call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, MPI_SUM, &
               0, MPI_COMM_WORLD, mpi_err)

          ! Transfer values to value on master
          t % scores(:,:) % value = tally_temp
       else
          ! Receive buffer not significant at other processors
          call MPI_REDUCE(tally_temp, tally_temp, n_bins, MPI_REAL8, MPI_SUM, &
               0, MPI_COMM_WORLD, mpi_err)

          ! Reset value on other processors
          t % scores(:,:) % value = 0
       end if

       deallocate(tally_temp)

    end do

  end subroutine reduce_tallies
#endif

!===============================================================================
! WRITE_TALLIES creates an output file and writes out the mean values of all
! tallies and their standard deviations
!===============================================================================

  subroutine write_tallies()

    integer :: i                          ! index in tallies array
    integer :: j                          ! level in tally hierarchy
    integer :: k                          ! loop index for scoring bins
    integer :: bins(N_FILTER_TYPES) = 0   ! bins corresponding to each filter
    integer :: indent                     ! number of spaces to preceed output
    integer :: io_error                   ! error in opening/writing file
    integer :: last_filter                ! lowest level filter type
    integer :: score_index                ! index in scores array for filters
    logical :: file_exists                ! does tallies.out file already exists? 
    logical :: has_filter(N_FILTER_TYPES) ! does tally have this filter?
    character(MAX_FILE_LEN) :: filename                    ! name of output file
    character(15)           :: filter_name(N_FILTER_TYPES) ! names of tally filters
    character(27)           :: score_name(N_SCORE_TYPES)   ! names of scoring function
    type(TallyObject), pointer :: t

    ! Skip if there are no tallies
    if (n_tallies == 0) return

    ! Initialize names for tally filter types
    filter_name(FILTER_UNIVERSE)  = "Universe"
    filter_name(FILTER_MATERIAL)  = "Material"
    filter_name(FILTER_CELL)      = "Cell"
    filter_name(FILTER_CELLBORN)  = "Birth Cell"
    filter_name(FILTER_SURFACE)   = "Surface"
    filter_name(FILTER_MESH)      = "Mesh"
    filter_name(FILTER_ENERGYIN)  = "Incoming Energy"
    filter_name(FILTER_ENERGYOUT) = "Outgoing Energy"

    ! Initialize names for scores
    score_name(abs(SCORE_FLUX))       = "Flux"
    score_name(abs(SCORE_TOTAL))      = "Total Reaction Rate"
    score_name(abs(SCORE_SCATTER))    = "Scattering Rate"
    score_name(abs(SCORE_NU_SCATTER)) = "Scattering Production Rate"
    score_name(abs(SCORE_SCATTER_1))  = "First Scattering Moment"
    score_name(abs(SCORE_SCATTER_2))  = "Second Scattering Moment"
    score_name(abs(SCORE_SCATTER_3))  = "Third Scattering Moment"
    score_name(abs(SCORE_DIFFUSION))  = "Diffusion Coefficient"
    score_name(abs(SCORE_N_1N))       = "(n,1n) Rate"
    score_name(abs(SCORE_N_2N))       = "(n,2n) Rate"
    score_name(abs(SCORE_N_3N))       = "(n,3n) Rate"
    score_name(abs(SCORE_N_4N))       = "(n,4n) Rate"
    score_name(abs(SCORE_ABSORPTION)) = "Absorption Rate"
    score_name(abs(SCORE_FISSION))    = "Fission Rate"
    score_name(abs(SCORE_NU_FISSION)) = "Nu-Fission Rate"

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
       call header("TALLY " // trim(to_str(t % id)), unit=UNIT_TALLY, level=3)

       ! Handle surface current tallies separately
       if (t % type == TALLY_SURFACE_CURRENT) then
          call write_surface_current(t)
          cycle
       end if

       ! First determine which filters this tally has
       do j = 1, N_FILTER_TYPES
          if (t % n_filter_bins(j) > 0) then
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

             if ((has_filter(j) .and. bins(j) > t % n_filter_bins(j)) .or. &
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
          do k = 1, t % n_score_bins
             write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') & 
                  repeat(" ", indent), score_name(abs(t % score_bins(k) % scalar)), &
                  to_str(t % scores(score_index,k) % sum), &
                  trim(to_str(t % scores(score_index,k) % sum_sq))
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

    integer :: i                    ! mesh index for x
    integer :: j                    ! mesh index for y
    integer :: k                    ! mesh index for z
    integer :: l                    ! mesh index for energy
    integer :: bins(N_FILTER_TYPES) ! bin combination
    integer :: n                    ! number of incoming energy bins
    integer :: len1                 ! length of string 
    integer :: len2                 ! length of string 
    integer :: score_index          ! index in scores array for filters
    logical :: print_ebin           ! should incoming energy bin be displayed?
    character(MAX_LINE_LEN) :: string
    type(StructuredMesh), pointer :: m => null()

    ! Get pointer to mesh
    m => meshes(t % mesh)

    ! initialize bins array
    bins = 1

    ! determine how many energy in bins there are
    n = t % n_filter_bins(FILTER_ENERGYIN)
    if (n > 0) then
       print_ebin = .true.
    else
       print_ebin = .false.
       n = 1
    end if

    do i = 1, m % dimension(1)
       string = "Mesh Index (" // trim(to_str(i)) // ", "
       len1 = len_trim(string)
       do j = 1, m % dimension(2)
          string = string(1:len1+1) // trim(to_str(j)) // ", "
          len2 = len_trim(string)
          do k = 1, m % dimension(3)
             ! Write mesh cell index
             string = string(1:len2+1) // trim(to_str(k)) // ")"
             write(UNIT=UNIT_TALLY, FMT='(1X,A)') trim(string)

             do l = 1, n
                ! Write incoming energy bin
                if (print_ebin) then
                   write(UNIT=UNIT_TALLY, FMT='(3X,A,1X,A)') &
                        "Incoming Energy", trim(get_label(t, FILTER_ENERGYIN, l))
                end if

                ! Set incoming energy bin
                bins(SURF_FILTER_ENERGYIN) = l

                ! Left Surface
                bins(1:3) = (/ i-1, j, k /) + 1
                bins(SURF_FILTER_SURFACE) = IN_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Left", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                bins(SURF_FILTER_SURFACE) = OUT_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Left", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                ! Right Surface
                bins(1:3) = (/ i, j, k /) + 1
                bins(SURF_FILTER_SURFACE) = IN_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Right", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                bins(SURF_FILTER_SURFACE) = OUT_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Right", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                ! Back Surface
                bins(1:3) = (/ i, j-1, k /) + 1
                bins(SURF_FILTER_SURFACE) = IN_FRONT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Back", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                bins(SURF_FILTER_SURFACE) = OUT_FRONT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Back", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                ! Front Surface
                bins(1:3) = (/ i, j, k /) + 1
                bins(SURF_FILTER_SURFACE) = IN_FRONT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Front", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                bins(SURF_FILTER_SURFACE) = OUT_FRONT
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Front", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                ! Bottom Surface
                bins(1:3) = (/ i, j, k-1 /) + 1
                bins(SURF_FILTER_SURFACE) = IN_TOP
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Bottom", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                bins(SURF_FILTER_SURFACE) = OUT_TOP
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Bottom", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                ! Top Surface
                bins(1:3) = (/ i, j, k /) + 1
                bins(SURF_FILTER_SURFACE) = IN_TOP
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Top", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))

                bins(SURF_FILTER_SURFACE) = OUT_TOP
                score_index = sum((bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Top", &
                     to_str(t % scores(score_index,1) % sum), &
                     trim(to_str(t % scores(score_index,1) % sum_sq))
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
    character(30)              :: label       ! user-specified identifier

    integer              :: i      ! index in cells/surfaces/etc array
    integer, allocatable :: ijk(:) ! indices in mesh
    real(8)              :: E0     ! lower bound for energy bin
    real(8)              :: E1     ! upper bound for energy bin
    type(StructuredMesh), pointer :: m

    select case(filter_type)
    case (FILTER_UNIVERSE)
       i = t % universe_bins(bin) % scalar
       label = to_str(universes(i) % id)
    case (FILTER_MATERIAL)
       i = t % material_bins(bin) % scalar
       label = to_str(materials(i) % id)
    case (FILTER_CELL)
       i = t % cell_bins(bin) % scalar
       label = to_str(cells(i) % id)
    case (FILTER_CELLBORN)
       i = t % cellborn_bins(bin) % scalar
       label = to_str(cells(i) % id)
    case (FILTER_SURFACE)
       i = t % surface_bins(bin) % scalar
       label = to_str(surfaces(i) % id)
    case (FILTER_MESH)
       m => meshes(t % mesh)
       allocate(ijk(m % n_dimension))
       call bin_to_mesh_indices(m, bin, ijk)
       if (m % n_dimension == 2) then
          label = "Index (" // trim(to_str(ijk(1))) // ", " // &
               trim(to_str(ijk(2))) // ")"
       elseif (m % n_dimension == 3) then
          label = "Index (" // trim(to_str(ijk(1))) // ", " // &
               trim(to_str(ijk(2))) // ", " // trim(to_str(ijk(3))) // ")"
       end if
    case (FILTER_ENERGYIN)
       E0 = t % energy_in(bin)
       E1 = t % energy_in(bin + 1)
       label = "[" // trim(to_str(E0)) // ", " // trim(to_str(E1)) // ")"
    case (FILTER_ENERGYOUT)
       E0 = t % energy_out(bin)
       E1 = t % energy_out(bin + 1)
       label = "[" // trim(to_str(E0)) // ", " // trim(to_str(E1)) // ")"
    end select

  end function get_label

!===============================================================================
! CALCULATE_STATISTICS determines the sample mean and the standard deviation of
! the mean for a TallyScore.
!===============================================================================

  elemental subroutine calculate_statistics(score)

    type(TallyScore), intent(inout) :: score

    ! Calculate sample mean and standard deviation of the mean -- note that we
    ! have used Bessel's correction so that the estimator of the variance of the
    ! sample mean is unbiased.

    score % sum    = score % sum/n_active
    score % sum_sq = sqrt((score % sum_sq/n_active - score % sum**2) / &
         (n_active - 1))

  end subroutine calculate_statistics

!===============================================================================
! TALLY_STATISTICS computes the mean and standard deviation of the mean of each
! tally and stores them in the val and val_sq attributes of the TallyScores
! respectively
!===============================================================================

  subroutine tally_statistics()

    integer :: i    ! index in tallies array
    type(TallyObject), pointer :: t => null()

    ! Calculate statistics for user-defined tallies
    do i = 1, n_tallies
       t => tallies(i)

       call calculate_statistics(t % scores)
    end do

    ! Calculate statistics for global tallies
    call calculate_statistics(global_tallies)

  end subroutine tally_statistics

end module tally
