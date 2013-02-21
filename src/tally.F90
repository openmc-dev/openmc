module tally

  use ace_header,       only: Reaction
  use constants
  use error,            only: fatal_error
  use global
  use math,             only: t_percentile, calc_pn
  use mesh,             only: get_mesh_bin, bin_to_mesh_indices, &
                              get_mesh_indices, mesh_indices_to_bin, &
                              mesh_intersects
  use mesh_header,      only: StructuredMesh
  use output,           only: header
  use particle_header,  only: LocalCoord
  use search,           only: binary_search
  use string,           only: to_str
  use tally_header,     only: TallyResult, TallyMapItem, TallyMapElement

#ifdef MPI
  use mpi
#endif

  implicit none

  ! Tally map positioning array
  integer :: position(N_FILTER_TYPES - 3) = 0

contains

!===============================================================================
! SCORE_ANALOG_TALLY keeps track of how many events occur in a specified cell,
! energy range, etc. Note that since these are "analog" tallies, they are only
! triggered at every collision, not every event
!===============================================================================

  subroutine score_analog_tally()

    integer :: i
    integer :: i_tally
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    integer :: n                    ! loop index for scattering order
    integer :: l                    ! scoring bin loop index, allowing for changing 
                                    ! position during the loop
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

    ! Copy particle's pre- and post-collision weight and angle
    last_wgt = p % last_wgt
    wgt = p % wgt
    mu = p % mu

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do i = 1, active_analog_tallies % size()
      ! Get index of tally and pointer to tally
      i_tally = active_analog_tallies % get_item(i)
      t => tallies(i_tally)

      ! =======================================================================
      ! DETERMINE SCORING BIN COMBINATION

      call get_scoring_bins(i_tally, found_bin)
      if (.not. found_bin) cycle

      ! =======================================================================
      ! CALCULATE RESULTS AND ACCUMULATE TALLY

      ! If we have made it here, we have a scoring combination of bins for this
      ! tally -- now we need to determine where in the results array we should
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
        j = 0
        SCORE_LOOP: do l = 1, t % n_user_score_bins
          j = j + 1
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

            if (survival_biasing) then
              ! We need to account for the fact that some weight was already
              ! absorbed
              score = last_wgt + p % absorb_wgt
            else
              score = last_wgt
            end if

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
            
          case (SCORE_SCATTER_N)
            ! Skip any event where the particle didn't scatter
            if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP

            ! Find the scattering order for a singly requested moment, and
            ! store its moment contribution.

            if (t % scatt_order(j) == 1) then
              score = last_wgt * mu ! avoid function call overhead
            else
              score = last_wgt * calc_pn(t % scatt_order(j), mu)
            endif

          case (SCORE_SCATTER_PN)
            ! Skip any event where the particle didn't scatter
            if (p % event /= EVENT_SCATTER) then
              j = j + t % scatt_order(j)
              cycle SCORE_LOOP
            end if
            score_index = score_index - 1
            ! Find the scattering order for a collection of requested moments
            ! and store the moment contribution of each
            do n = 0, t % scatt_order(j)
              ! determine scoring bin index
              score_index = score_index + 1
              ! get the score and tally it
              score = last_wgt * calc_pn(n, mu)
              
              t % results(score_index, filter_index) % value = &
                t % results(score_index, filter_index) % value + score
            end do
            j = j + t % scatt_order(j)
            cycle SCORE_LOOP

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

          case (SCORE_ABSORPTION)
            if (survival_biasing) then
              ! No absorption events actually occur if survival biasing is on --
              ! just use weight absorbed in survival biasing

              score = p % absorb_wgt

            else
              ! Skip any event where the particle wasn't absorbed
              if (p % event == EVENT_SCATTER) cycle SCORE_LOOP

              ! All fission and absorption events will contribute here, so we
              ! can just use the particle's weight entering the collision

              score = last_wgt
            end if

          case (SCORE_FISSION)
            if (survival_biasing) then
              ! No fission events occur if survival biasing is on -- need to
              ! calculate fraction of absorptions that would have resulted in
              ! fission

              score = p % absorb_wgt * micro_xs(p % event_nuclide) % fission / &
                   micro_xs(p % event_nuclide) % absorption

            else
              ! Skip any non-fission events
              if (p % event /= EVENT_FISSION) cycle SCORE_LOOP

              ! All fission events will contribute, so again we can use
              ! particle's weight entering the collision as the estimate for the
              ! fission reaction rate

              score = last_wgt
            end if

          case (SCORE_NU_FISSION)
            if (survival_biasing) then
              ! No fission events occur if survival biasing is on -- need to
              ! calculate fraction of absorptions that would have resulted in
              ! nu-fission

              score = p % absorb_wgt * micro_xs(p % event_nuclide) % &
                   nu_fission / micro_xs(p % event_nuclide) % absorption

            else
              ! Skip any non-fission events
              if (p % event /= EVENT_FISSION) cycle SCORE_LOOP

              if (t % find_filter(FILTER_ENERGYOUT) > 0) then
                ! Normally, we only need to make contributions to one scoring
                ! bin. However, in the case of fission, since multiple fission
                ! neutrons were emitted with different energies, multiple
                ! outgoing energy bins may have been scored to. The following
                ! logic treats this special case and results to multiple bins

                call score_fission_eout(t, score_index)
                cycle SCORE_LOOP

              else
                ! If there is no outgoing energy filter, than we only need to
                ! score to one bin. For the score to be 'analog', we need to
                ! score the number of particles that were banked in the fission
                ! bank. Since this was weighted by 1/keff, we multiply by keff
                ! to get the proper score.

                score = keff * p % wgt_bank

              end if
            end if
          
          case (SCORE_KAPPA_FISSION)
            if (survival_biasing) then
              ! No fission events occur if survival biasing is on -- need to
              ! calculate fraction of absorptions that would have resulted in
              ! fission and multiply by Q

              score = p % absorb_wgt * &
                      micro_xs(p % event_nuclide) % kappa_fission / &
                      micro_xs(p % event_nuclide) % absorption
              
            else
              ! Skip any non-fission events
              if (p % event /= EVENT_FISSION) cycle SCORE_LOOP

              ! All fission events will contribute, so again we can use
              ! particle's weight entering the collision as the estimate for
              ! the fission energy production rate
              
              n = nuclides(p % event_nuclide) % index_fission(1)
              score = last_wgt * &
                nuclides(p % event_nuclide) % reactions(n) % Q_value
            end if
          case (SCORE_EVENTS)
            ! Simply count number of scoring events
            score = ONE

          case default
            ! Any other score is assumed to be a MT number. Thus, we just need
            ! to check if it matches the MT number of the event
            if (p % event_MT /= score_bin) cycle SCORE_LOOP

            score = last_wgt

          end select

          ! Add score to tally
          t % results(score_index, filter_index) % value = &
               t % results(score_index, filter_index) % value + score

        end do SCORE_LOOP

      end do NUCLIDE_LOOP

      ! If the user has specified that we can assume all tallies are spatially
      ! separate, this implies that once a tally has been scored to, we needn't
      ! check the others. This cuts down on overhead when there are many
      ! tallies specified

      if (assume_separate) exit TALLY_LOOP

    end do TALLY_LOOP

    ! Reset tally map positioning
    position = 0

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
    integer :: n             ! number of energies on filter
    integer :: k             ! loop index for bank sites
    integer :: bin_energyout ! original outgoing energy bin
    integer :: i_filter      ! index for matching filter bin combination
    real(8) :: score         ! actualy score
    real(8) :: E_out         ! energy of fission bank site

    ! save original outgoing energy bin and score index
    i = t % find_filter(FILTER_ENERGYOUT)
    bin_energyout = t % matching_bins(i)

    ! Get number of energies on filter
    n = size(t % filters(i) % real_bins)

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

      ! check if outgoing energy is within specified range on filter
      if (E_out < t % filters(i) % real_bins(1) .or. &
           E_out > t % filters(i) % real_bins(n)) cycle

      ! change outgoing energy bin
      t % matching_bins(i) = binary_search(t % filters(i) % real_bins, n, E_out)

      ! determine scoring index
      i_filter = sum((t % matching_bins - 1) * t % stride) + 1

      ! Add score to tally
      t % results(i_score, i_filter) % value = &
           t % results(i_score, i_filter) % value + score
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

    integer :: i
    integer :: i_tally
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    integer :: l                    ! loop index for nuclides in material
    integer :: m                    ! loop index for reactions
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array (from bins)
    integer :: i_nuc                ! index in nuclides array (from material)
    integer :: i_energy             ! index in nuclide energy grid
    integer :: score_bin            ! scoring type, e.g. SCORE_FLUX
    integer :: score_index          ! scoring bin index
    real(8) :: f                    ! interpolation factor
    real(8) :: flux                 ! tracklength estimate of flux
    real(8) :: score                ! actual score (e.g., flux*xs)
    real(8) :: atom_density         ! atom density of single nuclide in atom/b-cm
    logical :: found_bin            ! scoring bin found?
    type(TallyObject), pointer :: t => null()
    type(Material),    pointer :: mat => null()
    type(Reaction),    pointer :: rxn => null()

    ! Determine track-length estimate of flux
    flux = p % wgt * distance

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do i = 1, active_tracklength_tallies % size()
      ! Get index of tally and pointer to tally
      i_tally = active_tracklength_tallies % get_item(i)
      t => tallies(i_tally)

      ! Check if this tally has a mesh filter -- if so, we treat it separately
      ! since multiple bins can be scored to with a single track

      if (t % find_filter(FILTER_MESH) > 0) then
        call score_tl_on_mesh(i_tally, distance)
        cycle
      end if

      ! =======================================================================
      ! DETERMINE SCORING BIN COMBINATION

      call get_scoring_bins(i_tally, found_bin)
      if (.not. found_bin) cycle

      ! =======================================================================
      ! CALCULATE RESULTS AND ACCUMULATE TALLY

      ! If we have made it here, we have a scoring combination of bins for this
      ! tally -- now we need to determine where in the results array we should
      ! be accumulating the tally values

      ! Determine scoring index for this filter combination
      filter_index = sum((t % matching_bins - 1) * t % stride) + 1

      if (t % all_nuclides) then
        call score_all_nuclides(i_tally, flux, filter_index)
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
              ! ================================================================
              ! DETERMINE NUCLIDE CROSS SECTION

              select case(score_bin)
              case (SCORE_FLUX)
                ! For flux, we need no cross section
                score = flux

              case (SCORE_TOTAL)
                ! Total cross section is pre-calculated
                score = micro_xs(i_nuclide) % total * &
                     atom_density * flux

              case (SCORE_SCATTER)
                ! Scattering cross section is pre-calculated
                score = (micro_xs(i_nuclide) % total - &
                     micro_xs(i_nuclide) % absorption) * &
                     atom_density * flux

              case (SCORE_ABSORPTION)
                ! Absorption cross section is pre-calculated
                score = micro_xs(i_nuclide) % absorption * &
                     atom_density * flux

              case (SCORE_FISSION)
                ! Fission cross section is pre-calculated
                score = micro_xs(i_nuclide) % fission * &
                     atom_density * flux

              case (SCORE_NU_FISSION)
                ! Nu-fission cross section is pre-calculated
                score = micro_xs(i_nuclide) % nu_fission * &
                     atom_density * flux

              case (SCORE_KAPPA_FISSION)
                score = micro_xs(i_nuclide) % kappa_fission * &
                     atom_density * flux

              case (SCORE_EVENTS)
                ! For number of events, just score unity
                score = ONE

              case default
                ! Any other cross section has to be calculated on-the-fly. For
                ! cross sections that are used often (e.g. n2n, ngamma, etc. for
                ! depletion), it might make sense to optimize this section or
                ! pre-calculate cross sections

                if (score_bin > 1) then
                  ! Set default score
                  score = ZERO

                  ! TODO: The following search for the matching reaction could
                  ! be replaced by adding a dictionary on each Nuclide instance
                  ! of the form {MT: i_reaction, ...}

                  REACTION_LOOP: do m = 1, nuclides(i_nuclide) % n_reaction
                    ! Get pointer to reaction
                    rxn => nuclides(i_nuclide) % reactions(m)

                    ! Check if this is the desired MT
                    if (score_bin == rxn % MT) then
                      ! Retrieve index on nuclide energy grid and interpolation
                      ! factor
                      i_energy = micro_xs(i_nuclide) % index_grid
                      f = micro_xs(i_nuclide) % interp_factor

                      if (i_energy >= rxn % threshold) then
                        score = ((ONE - f) * rxn % sigma(i_energy - &
                             rxn%threshold + 1) + f * rxn % sigma(i_energy - &
                             rxn%threshold + 2)) * atom_density * flux
                      end if

                      exit REACTION_LOOP
                    end if
                  end do REACTION_LOOP

                else
                  message = "Invalid score type on tally " // to_str(t % id) // "."
                  call fatal_error()
                end if
              end select

            else
              ! ================================================================
              ! DETERMINE MATERIAL CROSS SECTION

              select case(score_bin)
              case (SCORE_FLUX)
                ! For flux, we need no cross section
                score = flux

              case (SCORE_TOTAL)
                ! Total cross section is pre-calculated
                score = material_xs % total * flux

              case (SCORE_SCATTER)
                ! Scattering cross section is pre-calculated
                score = (material_xs % total - material_xs % absorption) * flux

              case (SCORE_ABSORPTION)
                ! Absorption cross section is pre-calculated
                score = material_xs % absorption * flux

              case (SCORE_FISSION)
                ! Fission cross section is pre-calculated
                score = material_xs % fission * flux

              case (SCORE_NU_FISSION)
                ! Nu-fission cross section is pre-calculated
                score = material_xs % nu_fission * flux

              case (SCORE_KAPPA_FISSION)
                score = material_xs % kappa_fission * flux

              case (SCORE_EVENTS)
                ! For number of events, just score unity
                score = ONE

              case default
                ! Any other cross section has to be calculated on-the-fly. This
                ! is somewhat costly since it requires a loop over each nuclide
                ! in a material and each reaction in the nuclide

                if (score_bin > 1) then
                  ! Set default score
                  score = ZERO

                  ! Get pointer to current material
                  mat => materials(p % material)

                  do l = 1, mat % n_nuclides
                    ! Get atom density
                    atom_density = mat % atom_density(l)
                    
                    ! Get index in nuclides array
                    i_nuc = mat % nuclide(l)

                    ! TODO: The following search for the matching reaction could
                    ! be replaced by adding a dictionary on each Nuclide
                    ! instance of the form {MT: i_reaction, ...}

                    do m = 1, nuclides(i_nuc) % n_reaction
                      ! Get pointer to reaction
                      rxn => nuclides(i_nuc) % reactions(m)

                      ! Check if this is the desired MT
                      if (score_bin == rxn % MT) then
                        ! Retrieve index on nuclide energy grid and interpolation
                        ! factor
                        i_energy = micro_xs(i_nuc) % index_grid
                        f = micro_xs(i_nuc) % interp_factor

                        if (i_energy >= rxn % threshold) then
                          score = score + ((ONE - f) * rxn % sigma(i_energy - &
                               rxn%threshold + 1) + f * rxn % sigma(i_energy - &
                               rxn%threshold + 2)) * atom_density * flux
                        end if

                        exit
                      end if
                    end do

                  end do

                else
                  message = "Invalid score type on tally " // to_str(t % id) // "."
                  call fatal_error()
                end if
              end select
            end if

            ! Determine scoring bin index
            score_index = (k - 1)*t % n_score_bins + j

            ! Add score to tally
            t % results(score_index, filter_index) % value = &
                 t % results(score_index, filter_index) % value + score

          end do SCORE_LOOP

        end do NUCLIDE_BIN_LOOP
      end if

      ! If the user has specified that we can assume all tallies are spatially
      ! separate, this implies that once a tally has been scored to, we needn't
      ! check the others. This cuts down on overhead when there are many
      ! tallies specified

      if (assume_separate) exit TALLY_LOOP

    end do TALLY_LOOP

    ! Reset tally map positioning
    position = 0

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
    integer :: m             ! loop index for reactions in nuclide
    integer :: i_nuclide     ! index in nuclides array
    integer :: score_bin     ! type of score, e.g. SCORE_FLUX
    integer :: score_index   ! scoring bin index
    integer :: i_energy      ! index in nuclide energy grid
    real(8) :: f             ! interpolation factor
    real(8) :: score         ! actual scoring tally value
    real(8) :: atom_density  ! atom density of single nuclide in atom/b-cm
    type(TallyObject), pointer :: t => null()
    type(Material),    pointer :: mat => null()
    type(Reaction),    pointer :: rxn => null()

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
        case (SCORE_FLUX)
          score = flux

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

        case (SCORE_KAPPA_FISSION)
          score = micro_xs(i_nuclide) % kappa_fission * atom_density * flux

        case (SCORE_EVENTS)
          score = ONE

        case default
          ! Any other cross section has to be calculated on-the-fly. For cross
          ! sections that are used often (e.g. n2n, ngamma, etc. for depletion),
          ! it might make sense to optimize this section or pre-calculate cross
          ! sections
          
          if (score_bin > 1) then
            ! Set default score
            score = ZERO

            ! TODO: The following search for the matching reaction could be
            ! replaced by adding a dictionary on each Nuclide instance of the
            ! form {MT: i_reaction, ...}

            REACTION_LOOP: do m = 1, nuclides(i_nuclide) % n_reaction
              ! Get pointer to reaction
              rxn => nuclides(i_nuclide) % reactions(m)

              ! Check if this is the desired MT
              if (score_bin == rxn % MT) then
                ! Retrieve index on nuclide energy grid and interpolation factor
                i_energy = micro_xs(i_nuclide) % index_grid
                f = micro_xs(i_nuclide) % interp_factor

                if (i_energy >= rxn % threshold) then
                  score = ((ONE - f) * rxn % sigma(i_energy - &
                       rxn%threshold + 1) + f * rxn % sigma(i_energy - &
                       rxn%threshold + 2)) * atom_density * flux
                end if

                exit REACTION_LOOP
              end if
            end do REACTION_LOOP

          else
            message = "Invalid score type on tally " // to_str(t % id) // "."
            call fatal_error()
          end if
        end select

        ! Determine scoring bin index based on what the index of the nuclide
        ! is in the nuclides array
        score_index = (i_nuclide - 1)*t % n_score_bins + j

        ! Add score to tally
        t % results(score_index, filter_index) % value = &
             t % results(score_index, filter_index) % value + score

      end do SCORE_LOOP

    end do NUCLIDE_LOOP

    ! ==========================================================================
    ! SCORE TOTAL MATERIAL REACTION RATES

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

      case (SCORE_KAPPA_FISSION)
        score = material_xs % kappa_fission * flux

      case (SCORE_EVENTS)
        score = ONE

      case default
        ! Any other cross section has to be calculated on-the-fly. This is
        ! somewhat costly since it requires a loop over each nuclide in a
        ! material and each reaction in the nuclide
        
        if (score_bin > 1) then
          ! Set default score
          score = ZERO

          ! Get pointer to current material
          mat => materials(p % material)

          do i = 1, mat % n_nuclides
            ! Get atom density
            atom_density = mat % atom_density(i)

            ! Get index in nuclides array
            i_nuclide = mat % nuclide(i)

            ! TODO: The following search for the matching reaction could
            ! be replaced by adding a dictionary on each Nuclide
            ! instance of the form {MT: i_reaction, ...}

            do m = 1, nuclides(i_nuclide) % n_reaction
              ! Get pointer to reaction
              rxn => nuclides(i_nuclide) % reactions(m)

              ! Check if this is the desired MT
              if (score_bin == rxn % MT) then
                ! Retrieve index on nuclide energy grid and interpolation
                ! factor
                i_energy = micro_xs(i_nuclide) % index_grid
                f = micro_xs(i_nuclide) % interp_factor

                if (i_energy >= rxn % threshold) then
                  score = score + ((ONE - f) * rxn % sigma(i_energy - &
                       rxn%threshold + 1) + f * rxn % sigma(i_energy - &
                       rxn%threshold + 2)) * atom_density * flux
                end if

                exit
              end if
            end do

          end do

        else
          message = "Invalid score type on tally " // to_str(t % id) // "."
          call fatal_error()
        end if
      end select

      ! Determine scoring bin index based on what the index of the nuclide
      ! is in the nuclides array
      score_index = n_nuclides_total*t % n_score_bins + j

      ! Add score to tally
      t % results(score_index, filter_index) % value = &
           t % results(score_index, filter_index) % value + score

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
    type(LocalCoord),     pointer :: coord => null()

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
    call get_mesh_indices(m, xyz0, ijk0(:m % n_dimension), start_in_mesh)
    call get_mesh_indices(m, xyz1, ijk1(:m % n_dimension), end_in_mesh)

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
        coord => p % coord0
        do while(associated(coord))
          position(FILTER_CELL) = 0
          t % matching_bins(i) = get_next_bin(FILTER_CELL, &
               coord % cell, i_tally)
          if (t % matching_bins(i) /= NO_BIN_FOUND) exit
          coord => coord % next
        end do
        nullify(coord)

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
                case (SCORE_FLUX)
                  score = flux
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
                case (SCORE_KAPPA_FISSION)
                  score = micro_xs(i_nuclide) % kappa_fission * atom_density * flux
                case (SCORE_EVENTS)
                  score = ONE
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
                case (SCORE_KAPPA_FISSION)
                  score = material_xs % kappa_fission * flux
                case (SCORE_EVENTS)
                  score = ONE
                case default
                  message = "Invalid score type on tally " // &
                       to_str(t % id) // "."
                  call fatal_error()
                end select
              end if

              ! Determine scoring bin index
              score_index = (b - 1)*t % n_score_bins + j

              ! Add score to tally
              t % results(score_index, filter_index) % value = &
                   t % results(score_index, filter_index) % value + score

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
    type(LocalCoord),     pointer :: coord => null()

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
        coord => p % coord0
        do while(associated(coord))
          position(FILTER_CELL) = 0
          t % matching_bins(i) = get_next_bin(FILTER_CELL, &
               coord % cell, i_tally)
          if (t % matching_bins(i) /= NO_BIN_FOUND) exit
          coord => coord % next
        end do
        nullify(coord)

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
        if (p % E < t % filters(i) % real_bins(1) .or. &
             p % E > t % filters(i) % real_bins(n + 1)) then
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

    integer :: i
    integer :: i_tally
    integer :: j                    ! loop indices
    integer :: k                    ! loop indices
    integer :: ijk0(3)              ! indices of starting coordinates
    integer :: ijk1(3)              ! indices of ending coordinates
    integer :: n_cross              ! number of surface crossings
    integer :: n                    ! number of incoming energy bins
    integer :: filter_index         ! index of scoring bin
    integer :: i_filter_mesh        ! index of mesh filter in filters array
    integer :: i_filter_surf        ! index of surface filter in filters
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

    TALLY_LOOP: do i = 1, active_current_tallies % size()
      ! Copy starting and ending location of particle
      xyz0 = p % last_xyz
      xyz1 = p % coord0 % xyz

      ! Get pointer to tally
      i_tally = active_current_tallies % get_item(i)
      t => tallies(i_tally)

      ! Get index for mesh and surface filters
      i_filter_mesh = t % find_filter(FILTER_MESH)
      i_filter_surf = t % find_filter(FILTER_SURFACE)

      ! Determine indices for starting and ending location
      m => meshes(t % filters(i_filter_mesh) % int_bins(1))
      call get_mesh_indices(m, xyz0, ijk0(:m % n_dimension), start_in_mesh)
      call get_mesh_indices(m, xyz1, ijk1(:m % n_dimension), end_in_mesh)

      ! Check to if start or end is in mesh -- if not, check if track still
      ! intersects with mesh
      if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) then
        if (.not. mesh_intersects(m, xyz0, xyz1)) then
          cycle
        end if
      end if

      ! Calculate number of surface crossings
      n_cross = sum(abs(ijk1 - ijk0))
      if (n_cross == 0) then
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
              t % matching_bins(i_filter_surf) = OUT_TOP
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
              filter_index = sum((t % matching_bins - 1) * t % stride) + 1
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if
          end do
        else
          do j = ijk0(3) - 1, ijk1(3), -1
            ijk0(3) = j
            if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
              t % matching_bins(i_filter_surf) = IN_TOP
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
              filter_index = sum((t % matching_bins - 1) * t % stride) + 1
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
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
              t % matching_bins(i_filter_surf) = OUT_FRONT
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
              filter_index = sum((t % matching_bins - 1) * t % stride) + 1
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if
          end do
        else
          do j = ijk0(2) - 1, ijk1(2), -1
            ijk0(2) = j
            if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
              t % matching_bins(i_filter_surf) = IN_FRONT
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
              filter_index = sum((t % matching_bins - 1) * t % stride) + 1
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
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
              t % matching_bins(i_filter_surf) = OUT_RIGHT
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
              filter_index = sum((t % matching_bins - 1) * t % stride) + 1
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if
          end do
        else
          do j = ijk0(1) - 1, ijk1(1), -1
            ijk0(1) = j
            if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
              t % matching_bins(i_filter_surf) = IN_RIGHT
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
              filter_index = sum((t % matching_bins - 1) * t % stride) + 1
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
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
        t % matching_bins(i_filter_surf) = 0

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
              t % matching_bins(i_filter_surf) = OUT_RIGHT
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
            end if
            ijk0(1) = ijk0(1) + 1
            xyz_cross(1) = xyz_cross(1) + m % width(1)
          else
            ! Crossing into left mesh cell -- this is treated as incoming
            ! current in (i-1,j,k)
            ijk0(1) = ijk0(1) - 1
            xyz_cross(1) = xyz_cross(1) - m % width(1)
            if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
              t % matching_bins(i_filter_surf) = IN_RIGHT
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
            end if
          end if
        elseif (distance == d(2)) then
          if (uvw(2) > 0) then
            ! Crossing into front mesh cell -- this is treated as outgoing
            ! current in (i,j,k)
            if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
              t % matching_bins(i_filter_surf) = OUT_FRONT
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
            end if
            ijk0(2) = ijk0(2) + 1
            xyz_cross(2) = xyz_cross(2) + m % width(2)
          else
            ! Crossing into back mesh cell -- this is treated as incoming
            ! current in (i,j-1,k)
            ijk0(2) = ijk0(2) - 1
            xyz_cross(2) = xyz_cross(2) - m % width(2)
            if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
              t % matching_bins(i_filter_surf) = IN_FRONT
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
            end if
          end if
        else if (distance == d(3)) then
          if (uvw(3) > 0) then
            ! Crossing into top mesh cell -- this is treated as outgoing
            ! current in (i,j,k)
            if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
              t % matching_bins(i_filter_surf) = OUT_TOP
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
            end if
            ijk0(3) = ijk0(3) + 1
            xyz_cross(3) = xyz_cross(3) + m % width(3)
          else
            ! Crossing into bottom mesh cell -- this is treated as incoming
            ! current in (i,j,k-1)
            ijk0(3) = ijk0(3) - 1
            xyz_cross(3) = xyz_cross(3) - m % width(3)
            if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) then
              t % matching_bins(i_filter_surf) = IN_TOP
              t % matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0 + 1, .true.)
            end if
          end if
        end if

        ! Determine scoring index
        if (t % matching_bins(i_filter_surf) > 0) then
          filter_index = sum((t % matching_bins - 1) * t % stride) + 1

          ! Check for errors
          if (filter_index <= 0 .or. filter_index > &
               t % total_filter_bins) then
            message = "Score index outside range."
            call fatal_error()
          end if

          ! Add to surface current tally
          t % results(1, filter_index) % value = &
               t % results(1, filter_index) % value + p % wgt
        end if

        ! Calculate new coordinates
        xyz0 = xyz0 + distance * uvw
      end do

    end do TALLY_LOOP

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

    integer :: i
    real(8) :: k_col ! Copy of batch collision estimate of keff
    real(8) :: k_abs ! Copy of batch absorption estimate of keff
    real(8) :: k_tra ! Copy of batch tracklength estimate of keff

#ifdef MPI
    ! Combine tally results onto master process
    if (reduce_tallies) call reduce_tally_results()
#endif

    ! Increase number of realizations (only used for global tallies)
    if (reduce_tallies) then
      n_realizations = n_realizations + 1
    else
      n_realizations = n_realizations + n_procs
    end if

    ! Accumulate on master only unless run is not reduced then do it on all
    if (master .or. (.not. reduce_tallies)) then
      ! Accumulate results for each tally
      do i = 1, active_tallies % size()
        call accumulate_tally(tallies(active_tallies % get_item(i)))
      end do

      if (run_mode == MODE_EIGENVALUE) then
        ! Get the current batch estimate of k_analog for displaying to output
        ! --- this has to be performed after reduce_tally_values and before
        ! accumulate_result

        k_batch(current_batch) = global_tallies(K_TRACKLENGTH) % value

        if (active_batches) then
          ! Accumulate products of different estimators of k
          k_col = global_tallies(K_COLLISION) % value / total_weight
          k_abs = global_tallies(K_ABSORPTION) % value / total_weight
          k_tra = global_tallies(K_TRACKLENGTH) % value / total_weight
          k_col_abs = k_col_abs + k_col * k_abs
          k_col_tra = k_col_tra + k_col * k_tra
          k_abs_tra = k_abs_tra + k_abs * k_tra
        end if
      end if

      ! Accumulate results for global tallies
      call accumulate_result(global_tallies)
    end if

  end subroutine synchronize_tallies

!===============================================================================
! REDUCE_TALLY_RESULTS collects all the results from tallies onto one processor
!===============================================================================

#ifdef MPI
  subroutine reduce_tally_results()

    integer :: i
    integer :: n      ! number of filter bins
    integer :: m      ! number of score bins
    integer :: n_bins ! total number of bins
    real(8), allocatable :: tally_temp(:,:) ! contiguous array of results
    real(8) :: global_temp(N_GLOBAL_TALLIES)
    real(8) :: dummy  ! temporary receive buffer for non-root reduces
    type(TallyObject), pointer :: t => null()

    do i = 1, active_tallies % size()
      t => tallies(active_tallies % get_item(i))

      m = t % total_score_bins
      n = t % total_filter_bins
      n_bins = m*n

      allocate(tally_temp(m,n))

      tally_temp = t % results(:,:) % value

      if (master) then
        ! The MPI_IN_PLACE specifier allows the master to copy values into
        ! a receive buffer without having a temporary variable
        call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

        ! Transfer values to value on master
        t % results(:,:) % value = tally_temp
      else
        ! Receive buffer not significant at other processors
        call MPI_REDUCE(tally_temp, dummy, n_bins, MPI_REAL8, &
             MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

        ! Reset value on other processors
        t % results(:,:) % value = 0
      end if

      deallocate(tally_temp)
    end do

    ! Copy global tallies into array to be reduced
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

  end subroutine reduce_tally_results
#endif

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

    ! Accumulate each TallyResult
    call accumulate_result(t % results)

  end subroutine accumulate_tally

!===============================================================================
! TALLY_STATISTICS computes the mean and standard deviation of the mean of each
! tally and stores them in the val and val_sq attributes of the TallyResults
! respectively
!===============================================================================

  subroutine tally_statistics()

    integer :: i    ! index in tallies array
    type(TallyObject), pointer :: t => null()

    ! Calculate statistics for user-defined tallies
    do i = 1, n_tallies
      t => tallies(i)

      call statistics_result(t % results, t % n_realizations)
    end do

    ! Calculate statistics for global tallies
    call statistics_result(global_tallies, n_realizations)

  end subroutine tally_statistics

!===============================================================================
! ACCUMULATE_RESULT accumulates results from many histories (or many generations)
! into a single realization of a random variable.
!===============================================================================

  elemental subroutine accumulate_result(this)

    type(TallyResult), intent(inout) :: this

    real(8) :: val

    ! Add the sum and square of the sum of contributions from a tally result to
    ! the variables sum and sum_sq. This will later allow us to calculate a
    ! variance on the tallies.

    val = this % value/total_weight
    this % sum    = this % sum    + val
    this % sum_sq = this % sum_sq + val*val

    ! Reset the single batch estimate
    this % value = ZERO

  end subroutine accumulate_result

!===============================================================================
! STATISTICS_RESULT determines the sample mean and the standard deviation of the
! mean for a TallyResult.
!===============================================================================

  elemental subroutine statistics_result(this, n)

    type(TallyResult), intent(inout) :: this
    integer,           intent(in)    :: n

    ! Calculate sample mean and standard deviation of the mean -- note that we
    ! have used Bessel's correction so that the estimator of the variance of the
    ! sample mean is unbiased.

    this % sum    = this % sum/n
    this % sum_sq = sqrt((this % sum_sq/n - this % sum * &
         this % sum) / (n - 1))

  end subroutine statistics_result

!===============================================================================
! RESET_RESULT zeroes out the value and accumulated sum and sum-squared for a
! single TallyResult.
!===============================================================================

  elemental subroutine reset_result(this)

    type(TallyResult), intent(inout) :: this

    this % value    = ZERO
    this % sum      = ZERO
    this % sum_sq   = ZERO

  end subroutine reset_result

!===============================================================================
! SETUP_ACTIVE_USERTALLIES
!===============================================================================

  subroutine setup_active_usertallies()

    integer :: i ! loop counter

    do i = 1, n_user_tallies
      ! Add tally to active tallies
      call active_tallies % add(i_user_tallies + i)

      ! Check what type of tally this is and add it to the appropriate list
      if (user_tallies(i) % type == TALLY_VOLUME) then
        if (user_tallies(i) % estimator == ESTIMATOR_ANALOG) then
          call active_analog_tallies % add(i_user_tallies + i)
        elseif (user_tallies(i) % estimator == ESTIMATOR_TRACKLENGTH) then
          call active_tracklength_tallies % add(i_user_tallies + i)
        end if
      elseif (user_tallies(i) % type == TALLY_SURFACE_CURRENT) then
        call active_current_tallies % add(i_user_tallies + i)
      end if
    end do

  end subroutine setup_active_usertallies

!===============================================================================
! SETUP_ACTIVE_CMFDTALLIES
!===============================================================================

  subroutine setup_active_cmfdtallies()

    integer :: i ! loop counter

    ! check to see if any of the active tally lists has been allocated
    if (active_tallies % size() > 0) then
      message = "Active tallies should not exist before CMFD tallies!"
      call fatal_error()
    else if (active_analog_tallies % size() > 0) then
      message = 'Active analog tallies should not exist before CMFD tallies!'
      call fatal_error()
    else if (active_tracklength_tallies % size() > 0) then
      message = "Active tracklength tallies should not exist before CMFD &
           &tallies!"
      call fatal_error()
    else if (active_current_tallies % size() > 0) then
      message = "Active current tallies should not exist before CMFD tallies!"
      call fatal_error()
    end if

    do i = 1, n_cmfd_tallies
      ! Add CMFD tally to active tallies
      call active_tallies % add(i_cmfd_tallies + i)

      ! Check what type of tally this is and add it to the appropriate list
      if (cmfd_tallies(i) % type == TALLY_VOLUME) then
        if (cmfd_tallies(i) % estimator == ESTIMATOR_ANALOG) then
          call active_analog_tallies % add(i_cmfd_tallies + i)
        elseif (cmfd_tallies(i) % estimator == ESTIMATOR_TRACKLENGTH) then
          call active_tracklength_tallies % add(i_cmfd_tallies + i)
        end if
      elseif (cmfd_tallies(i) % type == TALLY_SURFACE_CURRENT) then
        call active_current_tallies % add(i_cmfd_tallies + i)
      end if
    end do

  end subroutine setup_active_cmfdtallies

end module tally
