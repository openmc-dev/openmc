module tally

  use constants
  use datatypes_header, only: ListInt
  use error,         only: fatal_error
  use global
  use math,          only: t_percentile
  use mesh,          only: get_mesh_bin, bin_to_mesh_indices, get_mesh_indices, &
                           mesh_indices_to_bin, mesh_intersects
  use mesh_header,   only: StructuredMesh
  use output,        only: header
  use search,        only: binary_search
  use string,        only: to_str
  use tally_header,  only: TallyScore, TallyMapItem, TallyMapElement

#ifdef MPI
  use mpi
#endif

  implicit none

  ! Tally map positioning array
  integer :: position(N_FILTER_TYPES - 3) = 0

contains

!===============================================================================
! CREATE_TALLY_MAP creates a map that allows a quick determination of which
! tallies and bins need to be scored to when a particle makes a collision. This
! subroutine also sets the stride attribute for each tally as well as allocating
! storage for the scores array.
!===============================================================================

  subroutine create_tally_map()

    integer :: i    ! loop index for tallies
    integer :: j    ! loop index for filters
    integer :: k    ! loop index for bins
    integer :: bin  ! filter bin entries
    integer :: type ! type of tally filter
    type(TallyObject),    pointer :: t => null()

    ! allocate tally map array -- note that we don't need a tally map for the
    ! energy_in and energy_out filters
    allocate(tally_maps(N_FILTER_TYPES - 3))

    ! allocate list of items for each different filter type
    allocate(tally_maps(FILTER_UNIVERSE) % items(n_universes))
    allocate(tally_maps(FILTER_MATERIAL) % items(n_materials))
    allocate(tally_maps(FILTER_CELL)     % items(n_cells))
    allocate(tally_maps(FILTER_CELLBORN) % items(n_cells))
    allocate(tally_maps(FILTER_SURFACE)  % items(n_surfaces))

    TALLY_LOOP: do i = 1, n_tallies
       ! Get pointer to tally
       t => tallies(i)

       ! No need to set up tally maps for surface current tallies
       if (t % type == TALLY_SURFACE_CURRENT) cycle

       FILTER_LOOP: do j = 1, t % n_filters
          ! Determine type of filter
          type = t % filters(j) % type

          if (type == FILTER_CELL .or. type == FILTER_SURFACE .or. &
               type == FILTER_MATERIAL .or. type == FILTER_UNIVERSE .or. &
               type == FILTER_CELLBORN) then

             ! Add map elements
             BIN_LOOP: do k = 1, t % filters(j) % n_bins
                bin = t % filters(j) % int_bins(k)
                call add_map_element(tally_maps(type) % items(bin), i, k)
             end do BIN_LOOP
          end if

       end do FILTER_LOOP

    end do TALLY_LOOP

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

    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    integer :: filter_index         ! single index for single bin
    integer :: score_bin            ! scoring bin, e.g. SCORE_FLUX
    integer :: i_nuclide            ! index in nuclides array
    integer :: score_index          ! scoring bin index
    real(8) :: score                ! analog tally score
    real(8) :: last_wgt             ! pre-collision particle weight
    real(8) :: wgt                  ! post-collision particle weight
    real(8) :: mu                   ! cosine of angle of collision
    real(8) :: macro_total          ! material macro total xs
    real(8) :: macro_scatt          ! material macro scatt xs
    logical :: found_bin            ! scoring bin found?
    type(TallyObject), pointer :: t => null()
    type(ListInt), pointer :: curr_ptr => null()

    ! Copy particle's pre- and post-collision weight and angle
    last_wgt = p % last_wgt
    wgt = p % wgt
    mu = p % mu

    ! set current pointer to active analog tallies
    curr_ptr => active_analog_tallies

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do while (associated(curr_ptr))

       t => tallies(analog_tallies(curr_ptr % data))

       ! =======================================================================
       ! DETERMINE SCORING BIN COMBINATION

       call get_scoring_bins(analog_tallies(curr_ptr % data), found_bin)
       if (.not. found_bin) then
         curr_ptr => curr_ptr % next
         cycle
       end if

       ! =======================================================================
       ! CALCULATE SCORES AND ACCUMULATE TALLY

       ! If we have made it here, we have a scoring combination of bins for this
       ! tally -- now we need to determine where in the scores array we should
       ! be accumulating the tally values

       ! Determine scoring index for this filter combination
       filter_index = sum((t % matching_bins - 1) * t % stride) + 1

       ! Check for nuclide bins
       k = 0
       NUCLIDE_LOOP: do while (k < t % n_nuclide_bins)

          ! Increment the index in the list of nuclide bins
          k = k + 1

          if (t % all_nuclides) then
             ! In the case that the user has requested to tally all nuclides, we
             ! can take advantage of the fact that we know exactly how nuclide
             ! bins correspond to nuclide indices.
             if (k == 1) then
                ! If we just entered, set the nuclide bin index to the index in
                ! the nuclides array since this will match the index in the
                ! nuclide bin array.
                k = p % event_nuclide
             elseif (k == p % event_nuclide + 1) then
                ! After we've tallied the individual nuclide bin, we also need
                ! to contribute to the total material bin which is the last bin
                k = n_nuclides_total + 1
             else
                ! After we've tallied in both the individual nuclide bin and the
                ! total material bin, we're done
                exit
             end if

          else
             ! If the user has explicitly specified nuclides (or specified
             ! none), we need to search through the nuclide bin list one by
             ! one. First we need to get the value of the nuclide bin
             i_nuclide = t % nuclide_bins(k)

             ! Now compare the value against that of the colliding nuclide.
             if (i_nuclide /= p % event_nuclide .and. i_nuclide /= -1) cycle
          end if

          ! Determine score for each bin
          SCORE_LOOP: do j = 1, t % n_score_bins
             ! determine what type of score bin
             score_bin = t % score_bins(j)

             ! determine scoring bin index
             score_index = (k - 1)*t % n_score_bins + j

             select case (score_bin)
             case (SCORE_FLUX)
                ! All events score to a flux bin. We actually use a collision
                ! estimator since there is no way to count 'events' exactly for
                ! the flux

                score = last_wgt / material_xs % total

             case (SCORE_TOTAL)
                ! All events will score to the total reaction rate. We can just
                ! use the weight of the particle entering the collision as the
                ! score

                score = last_wgt

             case (SCORE_SCATTER)
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! Since only scattering events make it here, again we can use
                ! the weight entering the collision as the estimator for the
                ! reaction rate

                score = last_wgt

             case (SCORE_NU_SCATTER) 
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! For scattering production, we need to use the post-collision
                ! weight as the estimate for the number of neutrons exiting a
                ! reaction with neutrons in the exit channel

                score = wgt

             case (SCORE_SCATTER_1)
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! The first scattering moment can be determined by using the
                ! rate of scattering reactions multiplied by the cosine of the
                ! change in neutron's angle due to the collision

                score = last_wgt * mu

             case (SCORE_SCATTER_2)
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! The second scattering moment can be determined in a similar
                ! manner to the first scattering moment

                score = last_wgt * 0.5*(3.0*mu*mu - ONE)

             case (SCORE_SCATTER_3)
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! The first scattering moment can be determined by using the
                ! rate of scattering reactions multiplied by the cosine of the
                ! change in neutron's angle due to the collision

                score = last_wgt * 0.5*(5.0*mu*mu*mu - 3.0*mu)

             case (SCORE_TRANSPORT)
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! get material macros
                macro_total = material_xs % total
                macro_scatt = material_xs % total - material_xs % absorption

                ! Score total rate - p1 scatter rate Note estimator needs to be
                ! adjusted since tallying is only occuring when a scatter has
                ! happend. Effectively this means multiplying the estimator by
                ! total/scatter macro
                score = (macro_total - mu*macro_scatt)*(ONE/macro_scatt)

             case (SCORE_DIFFUSION)
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! Temporarily store the scattering cross section
                score = material_xs % total - material_xs % absorption

                ! Since this only gets tallied at every scattering event, the
                ! flux estimator is 1/Sigma_s. Therefore, the diffusion
                ! coefficient times flux is 1/(3*Sigma_s*(Sigma_t -
                ! mu*Sigma_s)).

                score = last_wgt / (3.0_8 * score * (material_xs % total - &
                     mu * score))

             case (SCORE_N_1N)
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! Skip any events where weight of particle changed
                if (wgt /= last_wgt) cycle SCORE_LOOP

                ! All events that reach this point are (n,1n) reactions
                score = last_wgt

             case (SCORE_N_2N)
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! Skip any scattering events where the weight didn't double
                if (nint(wgt/last_wgt) /= 2) cycle SCORE_LOOP

                ! All events that reach this point are (n,2n) reactions
                score = last_wgt

             case (SCORE_N_3N)
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! Skip any scattering events where the weight didn't double
                if (nint(wgt/last_wgt) /= 3) cycle SCORE_LOOP

                ! All events that reach this point are (n,3n) reactions
                score = last_wgt

             case (SCORE_N_4N)
                ! Skip any event where the particle didn't scatter
                if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

                ! Skip any scattering events where the weight didn't double
                if (nint(wgt/last_wgt) /= 4) cycle SCORE_LOOP

                ! All events that reach this point are (n,3n) reactions
                score = last_wgt

             case (SCORE_ABSORPTION)
                ! Skip any event where the particle wasn't absorbed
                if (p % event == EVENT_SCATTER) cycle SCORE_LOOP

                ! All fission and absorption events will contribute here, so we
                ! can just use the particle's weight entering the collision

                score = last_wgt

             case (SCORE_FISSION)
                ! Skip any non-fission events
                if (p % event /= EVENT_FISSION) cycle SCORE_LOOP

                ! All fission events will contribute, so again we can use
                ! particle's weight entering the collision as the estimate for
                ! the fission reaction rate

                score = last_wgt

             case (SCORE_NU_FISSION)
                ! Skip any non-fission events
                if (p % event /= EVENT_FISSION) cycle SCORE_LOOP

                if (t % find_filter(FILTER_ENERGYOUT) > 0) then
                   ! Normally, we only need to make contributions to one scoring
                   ! bin. However, in the case of fission, since multiple
                   ! fission neutrons were emitted with different energies,
                   ! multiple outgoing energy bins may have been scored to. The
                   ! following logic treats this special case and scores to
                   ! multiple bins

                   call score_fission_eout(t, score_index)
                   cycle SCORE_LOOP

                else
                   ! If there is no outgoing energy filter, than we only need to
                   ! score to one bin. For the score to be 'analog', we need to
                   ! score the number of particles that were banked in the
                   ! fission bank. Since this was weighted by 1/keff, we
                   ! multiply by keff to get the proper score.

                   score = keff * p % wgt_bank

                end if

             case default
                message = "Invalid score type on tally " // to_str(t % id) // "."
                call fatal_error()
             end select
             
             ! Add score to tally
             call add_to_score(t % scores(score_index, filter_index), score)

          end do SCORE_LOOP

       end do NUCLIDE_LOOP
             
       ! If the user has specified that we can assume all tallies are spatially
       ! separate, this implies that once a tally has been scored to, we needn't
       ! check the others. This cuts down on overhead when there are many
       ! tallies specified

       if (assume_separate) exit TALLY_LOOP

       ! select next active tally
       curr_ptr => curr_ptr % next

    end do TALLY_LOOP

    ! Reset tally map positioning
    position = 0

    nullify(curr_ptr)

  end subroutine score_analog_tally

!===============================================================================
! SCORE_FISSION_EOUT handles a special case where we need to store neutron
! production rate with an outgoing energy filter (think of a fission matrix). In
! this case, we may need to score to multiple bins if there were multiple
! neutrons produced with different energies.
!===============================================================================

  subroutine score_fission_eout(t, i_score)

    type(TallyObject), pointer :: t
    integer, intent(in)        :: i_score ! index for score

    integer :: i             ! index of outgoing energy filter
    integer :: k             ! loop index for bank sites
    integer :: bin_energyout ! original outgoing energy bin
    integer :: i_filter      ! index for matching filter bin combination
    real(8) :: score         ! actualy score
    real(8) :: E_out         ! energy of fission bank site

    ! save original outgoing energy bin and score index
    i = t % find_filter(FILTER_ENERGYOUT)
    bin_energyout = t % matching_bins(i)

    ! Since the creation of fission sites is weighted such that it is
    ! expected to create n_particles sites, we need to multiply the
    ! score by keff to get the true nu-fission rate. Otherwise, the sum
    ! of all nu-fission rates would be ~1.0.

    ! loop over number of particles banked
    do k = 1, p % n_bank
       ! determine score based on bank site weight and keff
       score = keff * fission_bank(n_bank - p % n_bank + k) % wgt

       ! determine outgoing energy from fission bank
       E_out = fission_bank(n_bank - p % n_bank + k) % E

       ! change outgoing energy bin
       t % matching_bins(i) = binary_search(t % filters(i) % real_bins, &
            size(t % filters(i) % real_bins), E_out)

       ! determine scoring index
       i_filter = sum((t % matching_bins - 1) * t % stride) + 1

       ! Add score to tally
       call add_to_score(t % scores(i_score, i_filter), score)
    end do

    ! reset outgoing energy bin and score index
    t % matching_bins(i) = bin_energyout

  end subroutine score_fission_eout

!===============================================================================
! SCORE_TRACKLENGTH_TALLY calculates fluxes and reaction rates based on the
! track-length estimate of the flux. This is triggered at every event (surface
! crossing, lattice crossing, or collision) and thus cannot be done for tallies
! that require post-collision information.
!===============================================================================

  subroutine score_tracklength_tally(distance)

    real(8), intent(in) :: distance

    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array
    integer :: score_bin            ! scoring type, e.g. SCORE_FLUX
    integer :: score_index          ! scoring bin index
    real(8) :: flux                 ! tracklength estimate of flux
    real(8) :: score                ! actual score (e.g., flux*xs)
    real(8) :: atom_density         ! atom density of single nuclide in atom/b-cm
    logical :: found_bin            ! scoring bin found?
    type(TallyObject), pointer :: t => null()
    type(Material), pointer :: mat => null()
    type(ListInt), pointer :: curr_ptr => null()

    ! Determine track-length estimate of flux
    flux = p % wgt * distance

    ! set current pointer to active tracklength tallies
    curr_ptr => active_tracklength_tallies

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do while (associated(curr_ptr))
       t => tallies(tracklength_tallies(curr_ptr % data))

       ! Check if this tally has a mesh filter -- if so, we treat it separately
       ! since multiple bins can be scored to with a single track
       
       if (t % find_filter(FILTER_MESH) > 0) then
          call score_tl_on_mesh(tracklength_tallies(curr_ptr % data), distance)
          curr_ptr => curr_ptr % next ! select next tally
          cycle
       end if

       ! =======================================================================
       ! DETERMINE SCORING BIN COMBINATION

       call get_scoring_bins(tracklength_tallies(curr_ptr % data), found_bin)
       if (.not. found_bin) then
         curr_ptr => curr_ptr % next
         cycle
       end if

       ! =======================================================================
       ! CALCULATE SCORES AND ACCUMULATE TALLY

       ! If we have made it here, we have a scoring combination of bins for this
       ! tally -- now we need to determine where in the scores array we should
       ! be accumulating the tally values

       ! Determine scoring index for this filter combination
       filter_index = sum((t % matching_bins - 1) * t % stride) + 1

       if (t % all_nuclides) then
          call score_all_nuclides(tracklength_tallies(curr_ptr % data), flux, &
                                  filter_index)
       else

          NUCLIDE_BIN_LOOP: do k = 1, t % n_nuclide_bins
             ! Get index of nuclide in nuclides array
             i_nuclide = t % nuclide_bins(k)

             if (i_nuclide > 0) then
                ! Get pointer to current material
                mat => materials(p % material)

                ! Determine if nuclide is actually in material
                NUCLIDE_MAT_LOOP: do j = 1, mat % n_nuclides
                   ! If index of nuclide matches the j-th nuclide listed in the
                   ! material, break out of the loop
                   if (i_nuclide == mat % nuclide(j)) exit

                   ! If we've reached the last nuclide in the material, it means
                   ! the specified nuclide to be tallied is not in this material
                   if (j == mat % n_nuclides) then
                      cycle NUCLIDE_BIN_LOOP
                   end if
                end do NUCLIDE_MAT_LOOP

                atom_density = mat % atom_density(j)
             end if

             ! Determine score for each bin
             SCORE_LOOP: do j = 1, t % n_score_bins
                ! determine what type of score bin
                score_bin = t % score_bins(j)

                if (i_nuclide > 0) then
                   ! Determine macroscopic nuclide cross section 
                   select case(score_bin)
                   case (SCORE_TOTAL)
                      score = micro_xs(i_nuclide) % total * &
                           atom_density * flux
                   case (SCORE_SCATTER)
                      score = (micro_xs(i_nuclide) % total - &
                           micro_xs(i_nuclide) % absorption) * &
                           atom_density * flux
                   case (SCORE_ABSORPTION)
                      score = micro_xs(i_nuclide) % absorption * &
                           atom_density * flux
                   case (SCORE_FISSION)
                      score = micro_xs(i_nuclide) % fission * &
                           atom_density * flux
                   case (SCORE_NU_FISSION)
                      score = micro_xs(i_nuclide) % nu_fission * &
                           atom_density * flux
                   case default
                      message = "Invalid score type on tally " // to_str(t % id) // "."
                      call fatal_error()
                   end select

                else
                   ! Determine macroscopic material cross section 
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
                end if

                ! Determine scoring bin index
                score_index = (k - 1)*t % n_score_bins + j

                ! Add score to tally
                call add_to_score(t % scores(score_index, filter_index), score)

             end do SCORE_LOOP

          end do NUCLIDE_BIN_LOOP
       end if

       ! If the user has specified that we can assume all tallies are spatially
       ! separate, this implies that once a tally has been scored to, we needn't
       ! check the others. This cuts down on overhead when there are many
       ! tallies specified

       if (assume_separate) exit TALLY_LOOP

       ! select next active tally
       curr_ptr => curr_ptr % next

    end do TALLY_LOOP

    ! Reset tally map positioning
    position = 0

    nullify(curr_ptr)

  end subroutine score_tracklength_tally

!===============================================================================
! SCORE_ALL_NUCLIDES tallies individual nuclide reaction rates specifically when
! the user requests <nuclides>all</nuclides>.
!===============================================================================

  subroutine score_all_nuclides(i_tally, flux, filter_index)

    integer, intent(in) :: i_tally
    real(8), intent(in) :: flux
    integer, intent(in) :: filter_index

    integer :: i             ! loop index for nuclides in material
    integer :: j             ! loop index for scoring bin types
    integer :: i_nuclide     ! index in nuclides array
    integer :: score_bin     ! type of score, e.g. SCORE_FLUX
    integer :: score_index   ! scoring bin index
    real(8) :: score         ! actual scoring tally value
    real(8) :: atom_density  ! atom density of single nuclide in atom/b-cm
    type(TallyObject), pointer :: t => null()
    type(Material),    pointer :: mat => null()

    ! Get pointer to tally
    t => tallies(i_tally)

    ! Get pointer to current material. We need this in order to determine what
    ! nuclides are in the material
    mat => materials(p % material)

    ! ==========================================================================
    ! SCORE ALL INDIVIDUAL NUCLIDE REACTION RATES

    NUCLIDE_LOOP: do i = 1, mat % n_nuclides

       ! Determine index in nuclides array and atom density for i-th nuclide in
       ! current material
       i_nuclide = mat % nuclide(i)
       atom_density = mat % atom_density(i)

       ! Loop over score types for each bin
       SCORE_LOOP: do j = 1, t % n_score_bins
          ! determine what type of score bin
          score_bin = t % score_bins(j)

          ! Determine macroscopic nuclide cross section 
          select case(score_bin)
          case (SCORE_TOTAL)
             score = micro_xs(i_nuclide) % total * atom_density * flux
          case (SCORE_SCATTER)
             score = (micro_xs(i_nuclide) % total - &
                  micro_xs(i_nuclide) % absorption) * atom_density * flux
          case (SCORE_ABSORPTION)
             score = micro_xs(i_nuclide) % absorption * atom_density * flux
          case (SCORE_FISSION)
             score = micro_xs(i_nuclide) % fission * atom_density * flux
          case (SCORE_NU_FISSION)
             score = micro_xs(i_nuclide) % nu_fission * atom_density * flux
          case default
             message = "Invalid score type on tally " // to_str(t % id) // "."
             call fatal_error()
          end select

          ! Determine scoring bin index based on what the index of the nuclide
          ! is in the nuclides array
          score_index = (i_nuclide - 1)*t % n_score_bins + j

          ! Add score to tally
          call add_to_score(t % scores(score_index, filter_index), score)

       end do SCORE_LOOP

    end do NUCLIDE_LOOP

    ! ==========================================================================
    ! SCORE ALL INDIVIDUAL NUCLIDE REACTION RATES

    ! Loop over score types for each bin
    MATERIAL_SCORE_LOOP: do j = 1, t % n_score_bins
       ! determine what type of score bin
       score_bin = t % score_bins(j)

       ! Determine macroscopic material cross section 
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

       ! Determine scoring bin index based on what the index of the nuclide
       ! is in the nuclides array
       score_index = n_nuclides_total*t % n_score_bins + j

       ! Add score to tally
       call add_to_score(t % scores(score_index, filter_index), score)

    end do MATERIAL_SCORE_LOOP

  end subroutine score_all_nuclides

!===============================================================================
! SCORE_TL_ON_MESH calculate fluxes and reaction rates based on the track-length
! estimate of the flux specifically for tallies that have mesh filters. For
! these tallies, it is possible to score to multiple mesh cells for each track.
!===============================================================================

  subroutine score_tl_on_mesh(i_tally, d_track)

    integer, intent(in) :: i_tally
    real(8), intent(in) :: d_track

    integer :: i                    ! loop index for filter/score bins
    integer :: j                    ! loop index for direction
    integer :: k                    ! loop index for mesh cell crossings
    integer :: b                    ! loop index for nuclide bins
    integer :: ijk0(3)              ! indices of starting coordinates
    integer :: ijk1(3)              ! indices of ending coordinates
    integer :: ijk_cross(3)         ! indices of mesh cell crossed
    integer :: n_cross              ! number of surface crossings
    integer :: filter_index         ! single index for single bin
    integer :: score_bin            ! scoring bin, e.g. SCORE_FLUX
    integer :: i_nuclide            ! index in nuclides array
    integer :: score_index          ! scoring bin index
    integer :: i_filter_mesh        ! index of mesh filter in filters array
    real(8) :: atom_density         ! density of individual nuclide in atom/b-cm
    real(8) :: flux                 ! tracklength estimate of flux
    real(8) :: score                ! actual score (e.g., flux*xs)
    real(8) :: uvw(3)               ! cosine of angle of particle
    real(8) :: xyz0(3)              ! starting/intermediate coordinates
    real(8) :: xyz1(3)              ! ending coordinates of particle
    real(8) :: xyz_cross(3)         ! coordinates of next boundary
    real(8) :: d(3)                 ! distance to each bounding surface
    real(8) :: distance             ! distance traveled in mesh cell
    logical :: found_bin            ! was a scoring bin found?
    logical :: start_in_mesh        ! starting coordinates inside mesh?
    logical :: end_in_mesh          ! ending coordinates inside mesh?
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()
    type(Material),       pointer :: mat => null()

    t => tallies(i_tally)
    t % matching_bins = 1

    ! ==========================================================================
    ! CHECK IF THIS TRACK INTERSECTS THE MESH

    ! Copy starting and ending location of particle
    xyz0 = p % coord0 % xyz - (d_track - TINY_BIT) * p % coord0 % uvw
    xyz1 = p % coord0 % xyz  - TINY_BIT * p % coord0 % uvw

    ! Get index for mesh filter
    i_filter_mesh = t % find_filter(FILTER_MESH)

    ! Determine indices for starting and ending location
    m => meshes(t % filters(i_filter_mesh) % int_bins(1))
    call get_mesh_indices(m, xyz0, ijk0, start_in_mesh)
    call get_mesh_indices(m, xyz1, ijk1, end_in_mesh)

    ! Check if start or end is in mesh -- if not, check if track still
    ! intersects with mesh
    if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) then
       if (.not. mesh_intersects(m, xyz0, xyz1)) return
    end if

    ! Reset starting and ending location
    xyz0 = p % coord0 % xyz - d_track * p % coord0 % uvw
    xyz1 = p % coord0 % xyz

    ! =========================================================================
    ! CHECK FOR SCORING COMBINATION FOR FILTERS OTHER THAN MESH

    FILTER_LOOP: do i = 1, t % n_filters

       select case (t % filters(i) % type)
       case (FILTER_UNIVERSE)
          ! determine next universe bin
          ! TODO: Account for multiple universes when performing this filter
          t % matching_bins(i) = get_next_bin(FILTER_UNIVERSE, &
               p % coord % universe, i_tally)

       case (FILTER_MATERIAL)
          t % matching_bins(i) = get_next_bin(FILTER_MATERIAL, &
               p % material, i_tally)

       case (FILTER_CELL)
          ! determine next cell bin
          ! TODO: Account for cells in multiple levels when performing this filter
          t % matching_bins(i) = get_next_bin(FILTER_CELL, &
               p % coord % cell, i_tally)

       case (FILTER_CELLBORN)
          ! determine next cellborn bin
          t % matching_bins(i) = get_next_bin(FILTER_CELLBORN, &
               p % cell_born, i_tally)

       case (FILTER_SURFACE)
          ! determine next surface bin
          t % matching_bins(i) = get_next_bin(FILTER_SURFACE, &
               p % surface, i_tally)

       case (FILTER_ENERGYIN)
          ! determine incoming energy bin
          k = t % filters(i) % n_bins

          ! check if energy of the particle is within energy bins
          if (p % E < t % filters(i) % real_bins(1) .or. &
               p % E > t % filters(i) % real_bins(k + 1)) then
             t % matching_bins(i) = NO_BIN_FOUND
          else
             ! search to find incoming energy bin
             t % matching_bins(i) = binary_search(t % filters(i) % real_bins, &
                  k + 1, p % E)
          end if

       end select
       
       ! Check if no matching bin was found
       if (t % matching_bins(i) == NO_BIN_FOUND) return

    end do FILTER_LOOP

    ! ==========================================================================
    ! DETERMINE WHICH MESH CELLS TO SCORE TO

    ! Calculate number of surface crossings
    n_cross = sum(abs(ijk1(:m % n_dimension) - ijk0(:m % n_dimension))) + 1

    ! Copy particle's direction
    uvw = p % coord0 % uvw

    ! Bounding coordinates
    do j = 1, m % n_dimension
       if (uvw(j) > 0) then
          xyz_cross(j) = m % lower_left(j) + ijk0(j) * m % width(j)
       else
          xyz_cross(j) = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
       end if
    end do

    MESH_LOOP: do k = 1, n_cross
       found_bin = .false.

       ! Calculate distance to each bounding surface. We need to treat special
       ! case where the cosine of the angle is zero since this would result in a
       ! divide-by-zero.

       if (k == n_cross) xyz_cross = xyz1
          
       do j = 1, m % n_dimension
          if (uvw(j) == 0) then
             d(j) = INFINITY
          else
             d(j) = (xyz_cross(j) - xyz0(j))/uvw(j)
          end if
       end do

       ! Determine the closest bounding surface of the mesh cell by calculating
       ! the minimum distance
          
       j = minloc(d(:m % n_dimension), 1)
       distance = d(j)

       ! Now use the minimum distance and diretion of the particle to determine
       ! which surface was crossed

       if (all(ijk0(:m % n_dimension) >= 1) .and. all(ijk0(:m % n_dimension) <= m % dimension)) then
          ijk_cross = ijk0
          found_bin = .true.
       end if

       ! Increment indices and determine new crossing point
       if (uvw(j) > 0) then
          ijk0(j) = ijk0(j) + 1
          xyz_cross(j) = xyz_cross(j) + m % width(j)
       else
          ijk0(j) = ijk0(j) - 1
          xyz_cross(j) = xyz_cross(j) - m % width(j)
       end if

       ! =======================================================================
       ! SCORE TO THIS MESH CELL

       if (found_bin) then
          ! Calculate track-length estimate of flux
          flux = p % wgt * distance

          ! Determine mesh bin
          t % matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, ijk_cross)

          ! Determining scoring index
          filter_index = sum((t % matching_bins - 1) * t % stride) + 1

          if (t % all_nuclides) then
             ! Score reaction rates for each nuclide in material
             call score_all_nuclides(i_tally, flux, filter_index)

          else
             NUCLIDE_BIN_LOOP: do b = 1, t % n_nuclide_bins
                ! Get index of nuclide in nuclides array
                i_nuclide = t % nuclide_bins(b)

                if (i_nuclide > 0) then
                   ! Get pointer to current material
                   mat => materials(p % material)

                   ! Determine if nuclide is actually in material
                   NUCLIDE_MAT_LOOP: do j = 1, mat % n_nuclides
                      ! If index of nuclide matches the j-th nuclide listed in
                      ! the material, break out of the loop
                      if (i_nuclide == mat % nuclide(j)) exit

                      ! If we've reached the last nuclide in the material, it
                      ! means the specified nuclide to be tallied is not in this
                      ! material
                      if (j == mat % n_nuclides) then
                         cycle NUCLIDE_BIN_LOOP
                      end if
                   end do NUCLIDE_MAT_LOOP

                   atom_density = mat % atom_density(j)
                end if

                ! Determine score for each bin
                SCORE_LOOP: do j = 1, t % n_score_bins
                   ! determine what type of score bin
                   score_bin = t % score_bins(j)

                   if (i_nuclide > 0) then
                      ! Determine macroscopic nuclide cross section 
                      select case(score_bin)
                      case (SCORE_TOTAL)
                         score = micro_xs(i_nuclide) % total * &
                              atom_density * flux
                      case (SCORE_SCATTER)
                         score = (micro_xs(i_nuclide) % total - &
                              micro_xs(i_nuclide) % absorption) * &
                              atom_density * flux
                      case (SCORE_ABSORPTION)
                         score = micro_xs(i_nuclide) % absorption * &
                              atom_density * flux
                      case (SCORE_FISSION)
                         score = micro_xs(i_nuclide) % fission * &
                              atom_density * flux
                      case (SCORE_NU_FISSION)
                         score = micro_xs(i_nuclide) % nu_fission * &
                              atom_density * flux
                      case default
                         message = "Invalid score type on tally " // &
                              to_str(t % id) // "."
                         call fatal_error()
                      end select

                   else
                      ! Determine macroscopic material cross section 
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
                         message = "Invalid score type on tally " // &
                              to_str(t % id) // "."
                         call fatal_error()
                      end select
                   end if

                   ! Determine scoring bin index
                   score_index = (b - 1)*t % n_score_bins + j

                   ! Add score to tally
                   call add_to_score(t % scores(score_index, filter_index), score)

                end do SCORE_LOOP

             end do NUCLIDE_BIN_LOOP
          end if
       end if

       ! Calculate new coordinates
       xyz0 = xyz0 + distance * uvw

    end do MESH_LOOP

  end subroutine score_tl_on_mesh

!===============================================================================
! GET_SCORING_BINS determines a combination of filter bins that should be scored
! for a tally based on the particle's current attributes.
!===============================================================================

  subroutine get_scoring_bins(i_tally, found_bin)

    integer, intent(in)     :: i_tally
    logical, intent(out)    :: found_bin

    integer :: i ! loop index for filters
    integer :: n ! number of bins for single filter
    real(8) :: E ! particle energy
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()

    found_bin = .true.
    t => tallies(i_tally)
    t % matching_bins = 1

    FILTER_LOOP: do i = 1, t % n_filters

       select case (t % filters(i) % type)
       case (FILTER_MESH)
          ! determine mesh bin
          m => meshes(t % filters(i) % int_bins(1))

          ! Determine if we're in the mesh first
          call get_mesh_bin(m, p % coord0 % xyz, t % matching_bins(i))

       case (FILTER_UNIVERSE)
          ! determine next universe bin
          ! TODO: Account for multiple universes when performing this filter
          t % matching_bins(i) = get_next_bin(FILTER_UNIVERSE, &
               p % coord % universe, i_tally)

       case (FILTER_MATERIAL)
          t % matching_bins(i) = get_next_bin(FILTER_MATERIAL, &
               p % material, i_tally)

       case (FILTER_CELL)
          ! determine next cell bin
          ! TODO: Account for cells in multiple levels when performing this filter
          t % matching_bins(i) = get_next_bin(FILTER_CELL, &
               p % coord % cell, i_tally)

       case (FILTER_CELLBORN)
          ! determine next cellborn bin
          t % matching_bins(i) = get_next_bin(FILTER_CELLBORN, &
               p % cell_born, i_tally)

       case (FILTER_SURFACE)
          ! determine next surface bin
          t % matching_bins(i) = get_next_bin(FILTER_SURFACE, &
               p % surface, i_tally)

       case (FILTER_ENERGYIN)
          ! determine incoming energy bin
          n = t % filters(i) % n_bins

          ! make sure the correct energy is used
          if (t % estimator == ESTIMATOR_TRACKLENGTH) then
             E = p % E
          else
             E = p % last_E
          end if

          ! check if energy of the particle is within energy bins
          if (E < t % filters(i) % real_bins(1) .or. &
               E > t % filters(i) % real_bins(n + 1)) then
             t % matching_bins(i) = NO_BIN_FOUND
          else
             ! search to find incoming energy bin
             t % matching_bins(i) = binary_search(t % filters(i) % real_bins, &
                  n + 1, E)
          end if

       case (FILTER_ENERGYOUT)
          ! determine outgoing energy bin
          n = t % filters(i) % n_bins

          ! check if energy of the particle is within energy bins
          if (E < t % filters(i) % real_bins(1) .or. &
               E > t % filters(i) % real_bins(n + 1)) then
             t % matching_bins(i) = NO_BIN_FOUND
          else
             ! search to find incoming energy bin
             t % matching_bins(i) = binary_search(t % filters(i) % real_bins, &
                  n + 1, p % E)
          end if

       end select

       ! If the current filter didn't match, exit this subroutine
       if (t % matching_bins(i) == NO_BIN_FOUND) then
          found_bin = .false.
          return
       end if

    end do FILTER_LOOP

  end subroutine get_scoring_bins

!===============================================================================
! SCORE_SURFACE_CURRENT tallies surface crossings in a mesh tally by manually
! determining which mesh surfaces were crossed
!===============================================================================

  subroutine score_surface_current()

    integer :: j                    ! loop indices
    integer :: k                    ! loop indices
    integer :: ijk0(3)              ! indices of starting coordinates
    integer :: ijk1(3)              ! indices of ending coordinates
    integer :: n_cross              ! number of surface crossings
    integer :: n                    ! number of incoming energy bins
    integer :: filter_index         ! index of scoring bin
    integer :: i_filter_mesh        ! index of mesh filter in filters array
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
    type(ListInt), pointer :: curr_ptr => null()

    ! set current pointer to active current tallies
    curr_ptr => active_current_tallies

    TALLY_LOOP: do while (associated(curr_ptr))
       ! Copy starting and ending location of particle
       xyz0 = p % last_xyz
       xyz1 = p % coord0 % xyz

       ! Get pointer to tally
       t => tallies(current_tallies(curr_ptr % data))

       ! Get index for filter mesh
       i_filter_mesh = t % find_filter(FILTER_MESH)

       ! Determine indices for starting and ending location
       m => meshes(t % filters(i_filter_mesh) % int_bins(1))
       call get_mesh_indices(m, xyz0, ijk0, start_in_mesh)
       call get_mesh_indices(m, xyz1, ijk1, end_in_mesh)

       ! Check to if start or end is in mesh -- if not, check if track still
       ! intersects with mesh
       if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) then
          if (.not. mesh_intersects(m, xyz0, xyz1)) then
            curr_ptr => curr_ptr % next ! select next tally
            cycle
           end if
       end if

       ! Calculate number of surface crossings
       n_cross = sum(abs(ijk1 - ijk0))
       if (n_cross == 0) then
         curr_ptr => curr_ptr % next ! select next tally
         cycle
       end if

       ! Copy particle's direction
       uvw = p % coord0 % uvw

       ! determine incoming energy bin
       j = t % find_filter(FILTER_ENERGYIN)
       if (j > 0) then
          n = t % filters(j) % n_bins
          ! check if energy of the particle is within energy bins
          if (p % E < t % filters(j) % real_bins(1) .or. &
               p % E > t % filters(j) % real_bins(n + 1)) then
            curr_ptr => curr_ptr % next ! select next tally
            cycle
          end if

          ! search to find incoming energy bin
          t % matching_bins(j) = binary_search(t % filters(j) % real_bins, &
               n + 1, p % E)
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
                   t % matching_bins(t % n_filters + 1) = OUT_TOP
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                   filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(1, filter_index), p % wgt)
                end if
             end do
          else
             do j = ijk0(3) - 1, ijk1(3), -1
                ijk0(3) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   t % matching_bins(t % n_filters + 1) = IN_TOP
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                   filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(1, filter_index), p % wgt)
                end if
             end do
          end if
          curr_ptr => curr_ptr % next ! select next tally
          cycle
       elseif (x_same .and. z_same) then
          ! Only y crossings
          if (uvw(2) > 0) then
             do j = ijk0(2), ijk1(2) - 1
                ijk0(2) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   t % matching_bins(t % n_filters + 1) = OUT_FRONT
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                   filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(1, filter_index), p % wgt)
                end if
             end do
          else
             do j = ijk0(2) - 1, ijk1(2), -1
                ijk0(2) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   t % matching_bins(t % n_filters + 1) = IN_FRONT
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                   filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(1, filter_index), p % wgt)
                end if
             end do
          end if
          curr_ptr => curr_ptr % next ! select next tally
          cycle
       elseif (y_same .and. z_same) then
          ! Only x crossings
          if (uvw(1) > 0) then
             do j = ijk0(1), ijk1(1) - 1
                ijk0(1) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   t % matching_bins(t % n_filters + 1) = OUT_RIGHT
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                   filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(1, filter_index), p % wgt)
                end if
             end do
          else
             do j = ijk0(1) - 1, ijk1(1), -1
                ijk0(1) = j
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   t % matching_bins(t % n_filters + 1) = IN_RIGHT
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                   filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                   call add_to_score(t % scores(1, filter_index), p % wgt)
                end if
             end do
          end if
          curr_ptr => curr_ptr % next ! select next tally
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
          t % matching_bins(t % n_filters + 1) = 0

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
                   t % matching_bins(t % n_filters + 1) = OUT_RIGHT
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                end if
                ijk0(1) = ijk0(1) + 1
                xyz_cross(1) = xyz_cross(1) + m % width(1)
             else
                ! Crossing into left mesh cell -- this is treated as incoming
                ! current in (i-1,j,k)
                ijk0(1) = ijk0(1) - 1
                xyz_cross(1) = xyz_cross(1) - m % width(1)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   t % matching_bins(t % n_filters + 1) = IN_RIGHT
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                end if
             end if
          elseif (distance == d(2)) then
             if (uvw(2) > 0) then
                ! Crossing into front mesh cell -- this is treated as outgoing
                ! current in (i,j,k)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   t % matching_bins(t % n_filters + 1) = OUT_FRONT
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                end if
                ijk0(2) = ijk0(2) + 1
                xyz_cross(2) = xyz_cross(2) + m % width(2)
             else
                ! Crossing into back mesh cell -- this is treated as incoming
                ! current in (i,j-1,k)
                ijk0(2) = ijk0(2) - 1
                xyz_cross(2) = xyz_cross(2) - m % width(2)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   t % matching_bins(t % n_filters + 1) = IN_FRONT
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                end if
             end if
          else if (distance == d(3)) then
             if (uvw(3) > 0) then
                ! Crossing into top mesh cell -- this is treated as outgoing
                ! current in (i,j,k)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   t % matching_bins(t % n_filters + 1) = OUT_TOP
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                end if
                ijk0(3) = ijk0(3) + 1
                xyz_cross(3) = xyz_cross(3) + m % width(3)
             else
                ! Crossing into bottom mesh cell -- this is treated as incoming
                ! current in (i,j,k-1)
                ijk0(3) = ijk0(3) - 1
                xyz_cross(3) = xyz_cross(3) - m % width(3)
                if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
                   t % matching_bins(t % n_filters + 1) = IN_TOP
                   t % matching_bins(i_filter_mesh) = &
                        mesh_indices_to_bin(m, ijk0 + 1)
                end if
             end if
          end if

          ! Determine scoring index
          if (t % matching_bins(t % n_filters + 1) > 0) then
             filter_index = sum((t % matching_bins - 1) * t % stride) + 1

             ! Check for errors
             if (filter_index <= 0 .or. filter_index > &
                  t % total_filter_bins) then
                message = "Score index outside range."
                call fatal_error()
             end if

             ! Add to surface current tally
             call add_to_score(t % scores(1, filter_index), p % wgt)
          end if

          ! Calculate new coordinates
          xyz0 = xyz0 + distance * uvw
       end do

     ! select next active tally
     curr_ptr => curr_ptr % next

    end do TALLY_LOOP

    nullify(curr_ptr)

  end subroutine score_surface_current

!===============================================================================
! GET_NEXT_BIN determines the next scoring bin for a particular filter variable
!===============================================================================

  function get_next_bin(filter_type, filter_value, i_tally) result(bin)

    integer, intent(in) :: filter_type  ! e.g. FILTER_MATERIAL
    integer, intent(in) :: filter_value ! value of filter, e.g. material 3
    integer, intent(in) :: i_tally      ! index of tally
    integer             :: bin          ! index of filter

    integer :: i_tally_check
    integer :: n

    ! If there are no scoring bins for this item, then return immediately
    if (.not. allocated(tally_maps(filter_type) % items(filter_value) % elements)) then
       bin = NO_BIN_FOUND
       return
    end if

    ! Check how many elements there are for this item
    n = size(tally_maps(filter_type) % items(filter_value) % elements)

    do
       ! Increment position in elements
       position(filter_type) = position(filter_type) + 1

       ! If we've reached the end of the array, there is no more bin to score to
       if (position(filter_type) > n) then
          position(filter_type) = 0
          bin = NO_BIN_FOUND
          return
       end if

       i_tally_check = tally_maps(filter_type) % items(filter_value) % &
            elements(position(filter_type)) % index_tally

       if (i_tally_check > i_tally) then
          ! Since the index being checked against is greater than the index we
          ! need (and the tally indices were added to elements sequentially), we
          ! know that no more bins will be scoring bins for this tally
          position(filter_type) = 0
          bin = NO_BIN_FOUND
          return
       elseif (i_tally_check == i_tally) then
          ! Found a match
          bin = tally_maps(filter_type) % items(filter_value) % &
               elements(position(filter_type)) % index_bin
          return
       end if

    end do

  end function get_next_bin

!===============================================================================
! SYNCHRONIZE_TALLIES accumulates the sum of the contributions from each history
! within the batch to a new random variable
!===============================================================================

  subroutine synchronize_tallies()

    type(ListInt), pointer :: curr_ptr => null()

#ifdef MPI
    if (reduce_tallies) call reduce_tally_values()
#endif

    ! Increase number of realizations (only used for global tallies)
    if (active_batches) then
       if (reduce_tallies) then
          n_realizations = n_realizations + 1
       else
          n_realizations = n_realizations + n_procs
       end if
    end if

    ! Accumulate on master only unless run is not reduced then do it on all
    if (master .or. (.not. reduce_tallies)) then
       ! Accumulate scores for each tally
       curr_ptr => active_tallies
       do while(associated(curr_ptr))
          call accumulate_tally(tallies(curr_ptr % data))
          curr_ptr => curr_ptr % next
       end do

       if (run_mode == MODE_CRITICALITY) then
          ! Before accumulating scores for global_tallies, we need to get the
          ! current batch estimate of k_analog for displaying to output
          if (active_batches) k_batch(current_batch) = global_tallies(K_ANALOG) % value
       end if

       ! Accumulate scores for global tallies
       if (active_batches) call accumulate_score(global_tallies)
    end if

    if (associated(curr_ptr)) nullify(curr_ptr)

  end subroutine synchronize_tallies

!===============================================================================
! REDUCE_TALLY_VALUES collects all the results from tallies onto one processor
!===============================================================================

#ifdef MPI
  subroutine reduce_tally_values()

    integer :: n      ! number of filter bins
    integer :: m      ! number of score bins
    integer :: n_bins ! total number of bins
    real(8), allocatable :: tally_temp(:,:) ! contiguous array of scores
    real(8) :: global_temp(N_GLOBAL_TALLIES)
    real(8) :: dummy  ! temporary receive buffer for non-root reduces
    type(TallyObject), pointer :: t => null()
    type(ListInt), pointer :: curr_ptr => null()

    curr_ptr => active_tallies
    do while(associated(curr_ptr))
       t => tallies(curr_ptr % data)

       m = t % total_score_bins
       n = t % total_filter_bins
       n_bins = m*n

       allocate(tally_temp(m,n))

       tally_temp = t % scores(:,:) % value

       if (master) then
          ! The MPI_IN_PLACE specifier allows the master to copy values into
          ! a receive buffer without having a temporary variable
          call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

          ! Transfer values to value on master
          t % scores(:,:) % value = tally_temp
       else
          ! Receive buffer not significant at other processors
          call MPI_REDUCE(tally_temp, dummy, n_bins, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

          ! Reset value on other processors
          t % scores(:,:) % value = 0
       end if

       deallocate(tally_temp)
       curr_ptr => curr_ptr % next
    end do

    ! Copy global tallies into array to be reduced
    if (active_batches) then
       global_temp = global_tallies(:) % value

       if (master) then
          call MPI_REDUCE(MPI_IN_PLACE, global_temp, N_GLOBAL_TALLIES, &
               MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

          ! Transfer values back to global_tallies on master
          global_tallies(:) % value = global_temp
       else
          ! Receive buffer not significant at other processors
          call MPI_REDUCE(global_temp, dummy, N_GLOBAL_TALLIES, &
               MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

          ! Reset value on other processors
          global_tallies(:) % value = ZERO
       end if
    end if

    ! We also need to determine the total starting weight of particles from the
    ! last realization
    if (master) then
       call MPI_REDUCE(MPI_IN_PLACE, total_weight, 1, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
    else
       ! Receive buffer not significant at other processors
       call MPI_REDUCE(total_weight, dummy, 1, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
    end if

    if (associated(curr_ptr)) nullify(curr_ptr)

  end subroutine reduce_tally_values
#endif

!===============================================================================
! WRITE_TALLIES creates an output file and writes out the mean values of all
! tallies and their standard deviations
!===============================================================================

  subroutine write_tallies()

    integer :: i            ! index in tallies array
    integer :: j            ! level in tally hierarchy
    integer :: k            ! loop index for scoring bins
    integer :: n            ! loop index for nuclides
    integer :: type         ! type of tally filter        
    integer :: indent       ! number of spaces to preceed output
    integer :: filter_index ! index in scores array for filters
    integer :: score_index  ! scoring bin index
    integer :: i_nuclide    ! index in nuclides array
    integer :: i_listing    ! index in xs_listings array
    real(8) :: t_value      ! t-values for confidence intervals
    real(8) :: alpha        ! significance level for CI
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
    score_name(abs(SCORE_TRANSPORT))  = "Transport Rate"
    score_name(abs(SCORE_DIFFUSION))  = "Diffusion Coefficient"
    score_name(abs(SCORE_N_1N))       = "(n,1n) Rate"
    score_name(abs(SCORE_N_2N))       = "(n,2n) Rate"
    score_name(abs(SCORE_N_3N))       = "(n,3n) Rate"
    score_name(abs(SCORE_N_4N))       = "(n,4n) Rate"
    score_name(abs(SCORE_ABSORPTION)) = "Absorption Rate"
    score_name(abs(SCORE_FISSION))    = "Fission Rate"
    score_name(abs(SCORE_NU_FISSION)) = "Nu-Fission Rate"

    ! Create filename for tally output
    if (run_mode == MODE_TALLIES) then
       filename = "tallies." // trim(to_str(restart_batch)) // ".out"
    else
       filename = "tallies.out"
    end if

    ! Open tally file for writing
    open(FILE=filename, UNIT=UNIT_TALLY, STATUS='replace', ACTION='write')

    ! Calculate t-value for confidence intervals
    if (confidence_intervals) then
       alpha = ONE - CONFIDENCE_LEVEL
       t_value = t_percentile(ONE - alpha/TWO, n_realizations - 1)
    end if

    TALLY_LOOP: do i = 1, n_tallies
       t => tallies(i)

       ! Multiply uncertainty by t-value
       if (confidence_intervals) then
          do k = 1, size(t % scores, 2)
             do j = 1, size(t % scores, 1)
                ! Calculate t-value for confidence intervals
                if (confidence_intervals) then
                   alpha = ONE - CONFIDENCE_LEVEL
                   t_value = t_percentile(ONE - alpha/TWO, t % n_realizations - 1)
                end if
                t % scores(j,k) % sum_sq = t_value * t % scores(j,k) % sum_sq
             end do
          end do
       end if

       ! Write header block
       if (t % label == "") then
          call header("TALLY " // trim(to_str(t % id)), unit=UNIT_TALLY, &
             level=3)
       else
          call header("TALLY " // trim(to_str(t % id)) // ": " &
             // trim(t % label), unit=UNIT_TALLY, level=3)
       endif
       
       ! Handle surface current tallies separately
       if (t % type == TALLY_SURFACE_CURRENT) then
          call write_surface_current(t)
          cycle
       end if

       ! WARNING: Admittedly, the logic for moving for printing scores is
       ! extremely confusing and took quite a bit of time to get correct. The
       ! logic is structured this way since it is not practical to have a do
       ! loop for each filter variable (given that only a few filters are likely
       ! to be used for a given tally.

       ! Initialize bins, filter level, and indentation
       t % matching_bins = 0
       j = 1
       indent = 0

       print_bin: do
          find_bin: do
             ! Increment bin combination
             t % matching_bins(j) = t % matching_bins(j) + 1

             ! =================================================================
             ! REACHED END OF BINS FOR THIS FILTER, MOVE TO NEXT FILTER

             if (t % matching_bins(j) > t % filters(j) % n_bins) then
                ! If this is the first filter, then exit
                if (j == 1) exit print_bin

                t % matching_bins(j) = 0
                j = j - 1
                indent = indent - 2

             ! =================================================================
             ! VALID BIN -- WRITE FILTER INFORMATION OR EXIT TO WRITE SCORES

             else
                ! Check for no filters
                if (t % n_filters == 0) exit find_bin

                ! Check if this is last filter
                if (j == t % n_filters) exit find_bin

                ! Print current filter information
                type = t % filters(j) % type
                write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                     trim(filter_name(type)), trim(get_label(t, j))
                indent = indent + 2
                j = j + 1
             end if

          end do find_bin

          ! Print filter information
          if (t % n_filters > 0) then
             type = t % filters(j) % type
             write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                  trim(filter_name(type)), trim(get_label(t, j))
          end if

          ! Determine scoring index for this bin combination -- note that unlike
          ! in the score_tally subroutine, we have to use max(bins,1) since all
          ! bins below the lowest filter level will be zeros

          filter_index = sum((max(t % matching_bins,1) - 1) * t % stride) + 1

          ! Write scores for this filter bin combination
          score_index = 0
          if (t % n_filters > 0) indent = indent + 2
          do n = 1, t % n_nuclide_bins
             ! Write label for nuclide
             i_nuclide = t % nuclide_bins(n)
             if (i_nuclide == -1) then
                write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                     "Total Material"
             else
                i_listing = nuclides(i_nuclide) % listing
                write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                     trim(xs_listings(i_listing) % alias)
             end if

             indent = indent + 2
             do k = 1, t % n_score_bins
                score_index = score_index + 1
                write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') & 
                     repeat(" ", indent), score_name(abs(t % score_bins(k))), &
                     to_str(t % scores(score_index,filter_index) % sum), &
                     trim(to_str(t % scores(score_index,filter_index) % sum_sq))
             end do
             indent = indent - 2

          end do
          indent = indent - 2

       end do print_bin

    end do TALLY_LOOP

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
    integer :: l                    ! index for energy
    integer :: i_filter_mesh
    integer :: i_filter_ein
    integer :: n                    ! number of incoming energy bins
    integer :: len1                 ! length of string 
    integer :: len2                 ! length of string 
    integer :: filter_index         ! index in scores array for filters
    logical :: print_ebin           ! should incoming energy bin be displayed?
    character(MAX_LINE_LEN) :: string
    type(StructuredMesh), pointer :: m => null()

    ! Get pointer to mesh
    i_filter_mesh = t % find_filter(FILTER_MESH)
    m => meshes(t % filters(i_filter_mesh) % int_bins(1))

    ! initialize bins array
    t % matching_bins = 1

    ! determine how many energy in bins there are
    i_filter_ein = t % find_filter(FILTER_ENERGYIN)
    if (i_filter_ein > 0) then
       print_ebin = .true.
       n = t % filters(i_filter_ein) % n_bins
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
                if (print_ebin) then
                   ! Set incoming energy bin
                   t % matching_bins(i_filter_ein) = l

                   ! Write incoming energy bin
                   write(UNIT=UNIT_TALLY, FMT='(3X,A,1X,A)') &
                        "Incoming Energy", trim(get_label(t, i_filter_ein))
                end if

                ! Left Surface
                t % matching_bins(i_filter_mesh) = &
                     mesh_indices_to_bin(m, (/ i-1, j, k /) + 1)
                t % matching_bins(t % n_filters + 1) = IN_RIGHT
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Left", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                t % matching_bins(t % n_filters + 1) = OUT_RIGHT
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Left", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                ! Right Surface
                t % matching_bins(i_filter_mesh) = &
                     mesh_indices_to_bin(m, (/ i, j, k /) + 1)
                t % matching_bins(t % n_filters + 1) = IN_RIGHT
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Right", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                t % matching_bins(t % n_filters + 1) = OUT_RIGHT
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Right", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                ! Back Surface
                t % matching_bins(i_filter_mesh) = &
                     mesh_indices_to_bin(m, (/ i, j-1, k /) + 1)
                t % matching_bins(t % n_filters + 1) = IN_FRONT
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Back", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                t % matching_bins(t % n_filters + 1) = OUT_FRONT
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Back", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                ! Front Surface
                t % matching_bins(i_filter_mesh) = &
                     mesh_indices_to_bin(m, (/ i, j, k /) + 1)
                t % matching_bins(t % n_filters + 1) = IN_FRONT
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Front", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                t % matching_bins(t % n_filters + 1) = OUT_FRONT
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Front", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                ! Bottom Surface
                t % matching_bins(i_filter_mesh) = &
                     mesh_indices_to_bin(m, (/ i, j, k-1 /) + 1)
                t % matching_bins(t % n_filters + 1) = IN_TOP
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Bottom", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                t % matching_bins(t % n_filters + 1) = OUT_TOP
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Bottom", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                ! Top Surface
                t % matching_bins(i_filter_mesh) = &
                     mesh_indices_to_bin(m, (/ i, j, k /) + 1)
                t % matching_bins(t % n_filters + 1) = IN_TOP
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Incoming Current from Top", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))

                t % matching_bins(t % n_filters + 1) = OUT_TOP
                filter_index = sum((t % matching_bins - 1) * t % stride) + 1
                write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') & 
                     "Outgoing Current to Top", &
                     to_str(t % scores(1,filter_index) % sum), &
                     trim(to_str(t % scores(1,filter_index) % sum_sq))
             end do

          end do
       end do
    end do

  end subroutine write_surface_current

!===============================================================================
! GET_LABEL returns a label for a cell/surface/etc given a tally, filter type,
! and corresponding bin
!===============================================================================

  function get_label(t, i_filter) result(label)

    type(TallyObject), pointer :: t        ! tally object
    integer, intent(in)        :: i_filter ! index in filters array
    character(30)              :: label    ! user-specified identifier

    integer :: i      ! index in cells/surfaces/etc array
    integer :: bin
    integer, allocatable :: ijk(:) ! indices in mesh
    real(8)              :: E0     ! lower bound for energy bin
    real(8)              :: E1     ! upper bound for energy bin
    type(StructuredMesh), pointer :: m => null()

    bin = t % matching_bins(i_filter)

    select case(t % filters(i_filter) % type)
    case (FILTER_UNIVERSE)
       i = t % filters(i_filter) % int_bins(bin)
       label = to_str(universes(i) % id)
    case (FILTER_MATERIAL)
       i = t % filters(i_filter) % int_bins(bin)
       label = to_str(materials(i) % id)
    case (FILTER_CELL, FILTER_CELLBORN)
       i = t % filters(i_filter) % int_bins(bin)
       label = to_str(cells(i) % id)
    case (FILTER_SURFACE)
       i = t % filters(i_filter) % int_bins(bin)
       label = to_str(surfaces(i) % id)
    case (FILTER_MESH)
       m => meshes(t % filters(i_filter) % int_bins(1))
       allocate(ijk(m % n_dimension))
       call bin_to_mesh_indices(m, bin, ijk)
       if (m % n_dimension == 2) then
          label = "Index (" // trim(to_str(ijk(1))) // ", " // &
               trim(to_str(ijk(2))) // ")"
       elseif (m % n_dimension == 3) then
          label = "Index (" // trim(to_str(ijk(1))) // ", " // &
               trim(to_str(ijk(2))) // ", " // trim(to_str(ijk(3))) // ")"
       end if
    case (FILTER_ENERGYIN, FILTER_ENERGYOUT)
       E0 = t % filters(i_filter) % real_bins(bin)
       E1 = t % filters(i_filter) % real_bins(bin + 1)
       label = "[" // trim(to_str(E0)) // ", " // trim(to_str(E1)) // ")"
    end select

  end function get_label

!===============================================================================
! ACCUMULATE_TALLY
!===============================================================================

  subroutine accumulate_tally(t)

    type(TallyObject), intent(inout) :: t

    ! Increment number of realizations
    if (reduce_tallies) then
       t % n_realizations = t % n_realizations + 1
    else
       t % n_realizations = t % n_realizations + n_procs
    end if

    ! Accumulate each TallyScore
    call accumulate_score(t % scores)

  end subroutine accumulate_tally

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

       call statistics_score(t % scores, t % n_realizations)
    end do

    ! Calculate statistics for global tallies
    call statistics_score(global_tallies, n_realizations)

  end subroutine tally_statistics

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
! ACCUMULATE_SCORE accumulates scores from many histories (or many generations)
! into a single realization of a random variable.
!===============================================================================

  elemental subroutine accumulate_score(score)

    type(TallyScore), intent(inout) :: score

    real(8) :: val

    ! Add the sum and square of the sum of contributions from each cycle
    ! within a cycle to the variables sum and sum_sq. This will later allow us
    ! to calculate a variance on the tallies

    val = score % value/total_weight
    score % sum    = score % sum    + val
    score % sum_sq = score % sum_sq + val*val

    ! Reset the single batch estimate
    score % value = ZERO

  end subroutine accumulate_score

!===============================================================================
! STATISTICS_SCORE determines the sample mean and the standard deviation of the
! mean for a TallyScore.
!===============================================================================

  elemental subroutine statistics_score(score, n)

    type(TallyScore), intent(inout) :: score
    integer,          intent(in)    :: n

    ! Calculate sample mean and standard deviation of the mean -- note that we
    ! have used Bessel's correction so that the estimator of the variance of the
    ! sample mean is unbiased.

    score % sum    = score % sum/n
    score % sum_sq = sqrt((score % sum_sq/n - score % sum * &
         score % sum) / (n - 1))

  end subroutine statistics_score

!===============================================================================
! RESET_SCORE zeroes out the value and accumulated sum and sum-squared for a
! single TallyScore.
!===============================================================================

  elemental subroutine reset_score(score)

    type(TallyScore), intent(inout) :: score

    score % n_events = 0
    score % value    = ZERO
    score % sum      = ZERO
    score % sum_sq   = ZERO

  end subroutine reset_score

!===============================================================================
! SETUP_ACTIVE_USERTALLIES
!===============================================================================

  subroutine setup_active_usertallies()

    integer                  :: i       ! loop counter
    type(ListInt), pointer :: curr_ptr  ! pointer to current list node
    type(ListInt), pointer :: tall_ptr  ! pointer to active tallies only

    !============================================================
    ! TALLIES

    ! traverse to end of tally list
    tall_ptr => active_tallies
    do while(associated(tall_ptr))
      tall_ptr => tall_ptr % next
    end do

    !============================================================
    ! ANALOG TALLIES

    ! check to see if analog tallies have already been allocated
    curr_ptr => null()
    if (associated(active_analog_tallies)) then

      ! traverse to the end of the linked list
      curr_ptr => active_analog_tallies
      do while(associated(curr_ptr))
        curr_ptr => curr_ptr % next
      end do

    end if

    ! append all analog tallies
    do i = n_user_analog_tallies,1,-1

      ! allocate node
      allocate(curr_ptr)

      ! set the tally index
      curr_ptr % data = i
      curr_ptr % next => active_analog_tallies
      active_analog_tallies => curr_ptr

      ! set indices in active tallies
      allocate(tall_ptr)
      tall_ptr % data = analog_tallies(i)
      tall_ptr % next => active_tallies
      active_tallies => tall_ptr

    end do

    !============================================================
    ! TRACKLENGTH TALLIES

    ! check to see if tracklength tallies have already been allocated
    curr_ptr => null()
    if (associated(active_tracklength_tallies)) then

      ! traverse to the end of the linked list
      curr_ptr => active_tracklength_tallies
      do while(associated(curr_ptr))
        curr_ptr => curr_ptr % next
      end do

    end if


    ! append all tracklength tallies
    do i = n_user_tracklength_tallies,1,-1

      ! allocate node
      allocate(curr_ptr)

      ! set the tally index
      curr_ptr % data = i
      curr_ptr % next => active_tracklength_tallies
      active_tracklength_tallies => curr_ptr

      ! set indices in active tallies
      allocate(tall_ptr)
      tall_ptr % data = tracklength_tallies(i)
      tall_ptr % next => active_tallies
      active_tallies => tall_ptr

    end do

   !============================================================
   ! CURRENT TALLIES

    ! check to see if current tallies have already been allocated
    curr_ptr => null()
    if (associated(active_current_tallies)) then

      ! traverse to the end of the linked list
      curr_ptr => active_current_tallies
      do while(associated(curr_ptr))
        curr_ptr => curr_ptr % next
      end do

    end if

    ! append all current tallies
    do i = n_user_current_tallies,1,-1

      ! allocate node
      allocate(curr_ptr)

      ! set the tally index
      curr_ptr % data = i
      curr_ptr % next => active_current_tallies
      active_current_tallies => curr_ptr

      ! set indices in active tallies
      allocate(tall_ptr)
      tall_ptr % data = current_tallies(i)
      tall_ptr % next => active_tallies
      active_tallies => tall_ptr

    end do

    ! nullify the temporary pointer
    if (associated(curr_ptr)) nullify(curr_ptr)

  end subroutine setup_active_usertallies

!===============================================================================
! SETUP_ACTIVE_CMFDTALLIES
!===============================================================================

  subroutine setup_active_cmfdtallies()

    integer                :: i         ! loop counter
    type(ListInt), pointer :: curr_ptr  ! pointer to current list node
    type(ListInt), pointer :: tall_ptr  ! pointer to active tallies only

    ! check to see if actives tallies has been allocated
    tall_ptr => active_tallies
    if (associated(active_tallies)) then
      message = 'Active tallies should not exist before CMFD tallies!'
      call fatal_error()
    end if

    ! check to see if analog tallies have already been allocated
    curr_ptr => null()
    if (associated(active_analog_tallies)) then
      message = 'Active analog tallies should not exist before CMFD tallies!'
      call fatal_error()
    end if

    do i = n_cmfd_analog_tallies + n_user_analog_tallies, n_user_analog_tallies + 1, -1

      ! allocate node
      allocate(curr_ptr)

      ! set the tally index
      curr_ptr % data = i
      curr_ptr % next => active_analog_tallies
      active_analog_tallies => curr_ptr

      ! set indices in active tallies
      allocate(tall_ptr)
      tall_ptr % data = analog_tallies(i)
      tall_ptr % next => active_tallies
      active_tallies => tall_ptr

    end do

    ! check to see if tracklength tallies have already been allocated
    curr_ptr => null()
    if (associated(active_tracklength_tallies)) then
      message = 'Active tracklength tallies should not exist before CMFD tallies!'
      call fatal_error()
    end if

    do i = n_cmfd_tracklength_tallies + n_user_tracklength_tallies, n_user_tracklength_tallies + 1, -1

      ! allocate node
      allocate(curr_ptr)

      ! set the tally index
      curr_ptr % data = i
      curr_ptr % next => active_tracklength_tallies
      active_tracklength_tallies => curr_ptr

      ! set indices in active tallies
      allocate(tall_ptr)
      tall_ptr % data = tracklength_tallies(i)
      tall_ptr % next => active_tallies
      active_tallies => tall_ptr

    end do

    ! check to see if current tallies have already been allocated
    curr_ptr => null()
    if (associated(active_current_tallies)) then
      message = 'Active current tallies should not exist before CMFD tallies!'
      call fatal_error()
    end if

    do i = n_cmfd_current_tallies + n_user_current_tallies, n_user_current_tallies + 1, -1

      ! allocate node
      allocate(curr_ptr)

      ! set the tally index
      curr_ptr % data = i
      curr_ptr % next => active_current_tallies
      active_current_tallies => curr_ptr

      ! set indices in active tallies
      allocate(tall_ptr)
      tall_ptr % data = current_tallies(i)
      tall_ptr % next => active_tallies
      active_tallies => tall_ptr

    end do

  end subroutine setup_active_cmfdtallies

end module tally
