module tally

  use, intrinsic :: ISO_C_BINDING

  use algorithm,        only: binary_search
  use constants
  use cross_section,    only: multipole_deriv_eval, calculate_elastic_xs
  use dict_header,      only: EMPTY
  use error,            only: fatal_error
  use geometry_header
  use math,             only: t_percentile, calc_pn, calc_rn
  use mesh_header,      only: RegularMesh, meshes
  use message_passing
  use mgxs_header
  use nuclide_header
  use output,           only: header
  use particle_header,  only: LocalCoord, Particle
  use settings
  use simulation_header
  use string,           only: to_str
  use tally_derivative_header, only: tally_derivs
  use tally_filter
  use tally_header

  implicit none

  procedure(score_general_),      pointer :: score_general => null()
  procedure(score_analog_tally_), pointer :: score_analog_tally => null()

  abstract interface
    subroutine score_general_(p, t, start_index, filter_index, i_nuclide, &
                              atom_density, flux)
      import Particle
      import TallyObject
      type(Particle),    intent(in)    :: p
      type(TallyObject), intent(inout) :: t
      integer,            intent(in)   :: start_index
      integer,            intent(in)   :: i_nuclide
      integer,            intent(in)   :: filter_index   ! for % results
      real(8),            intent(in)   :: flux           ! flux estimate
      real(8),            intent(in)   :: atom_density   ! atom/b-cm
    end subroutine score_general_

    subroutine score_analog_tally_(p)
      import Particle
      type(Particle), intent(in) :: p
    end subroutine score_analog_tally_
  end interface

contains

!===============================================================================
! INIT_TALLY_ROUTINES Sets the procedure pointers needed for minimizing code
! with the CE and MG modes.
!===============================================================================

  subroutine init_tally_routines()
    if (run_CE) then
      score_general      => score_general_ce
      score_analog_tally => score_analog_tally_ce
    else
      score_general      => score_general_mg
      score_analog_tally => score_analog_tally_mg
    end if
  end subroutine init_tally_routines

!===============================================================================
! SCORE_GENERAL* adds scores to the tally array for the given filter and
! nuclide.  This function is called by all volume tallies.  For analog tallies,
! the flux estimate depends on the score type so the flux argument is really
! just used for filter weights.  The atom_density argument is not used for
! analog tallies.
!===============================================================================

  subroutine score_general_ce(p, t, start_index, filter_index, i_nuclide, &
       atom_density, flux)
    type(Particle),    intent(in)    :: p
    type(TallyObject), intent(inout) :: t
    integer,           intent(in)    :: start_index
    integer,           intent(in)    :: i_nuclide
    integer,           intent(in)    :: filter_index   ! for % results
    real(8),           intent(in)    :: flux           ! flux estimate
    real(8),           intent(in)    :: atom_density   ! atom/b-cm

    integer :: i                    ! loop index for scoring bins
    integer :: l                    ! loop index for nuclides in material
    integer :: m                    ! loop index for reactions
    integer :: q                    ! loop index for scoring bins
    integer :: i_temp               ! temperature index
    integer :: i_nuc                ! index in nuclides array (from material)
    integer :: i_energy             ! index in nuclide energy grid
    integer :: score_bin            ! scoring bin, e.g. SCORE_FLUX
    integer :: score_index          ! scoring bin index
    integer :: d                    ! delayed neutron index
    integer :: g                    ! delayed neutron index
    integer :: k                    ! loop index for bank sites
    integer :: d_bin                ! delayed group bin index
    integer :: dg_filter            ! index of delayed group filter
    real(8) :: yield                ! delayed neutron yield
    real(8) :: atom_density_        ! atom/b-cm
    real(8) :: f                    ! interpolation factor
    real(8) :: score                ! analog tally score
    real(8) :: E                    ! particle energy

    ! Pre-collision energy of particle
    E = p % last_E

    i = 0
    SCORE_LOOP: do q = 1, t % n_user_score_bins
      i = i + 1

      ! determine what type of score bin
      score_bin = t % score_bins(i)

      ! determine scoring bin index
      score_index = start_index + i

      !#########################################################################
      ! Determine appropirate scoring value.

      select case(score_bin)


      case (SCORE_FLUX, SCORE_FLUX_YN)
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! All events score to a flux bin. We actually use a collision
          ! estimator in place of an analog one since there is no way to count
          ! 'events' exactly for the flux
          if (survival_biasing) then
            ! We need to account for the fact that some weight was already
            ! absorbed
            score = p % last_wgt + p % absorb_wgt
          else
            score = p % last_wgt
          end if
          score = score / material_xs % total * flux

        else
          ! For flux, we need no cross section
          score = flux
        end if


      case (SCORE_TOTAL, SCORE_TOTAL_YN)
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! All events will score to the total reaction rate. We can just
          ! use the weight of the particle entering the collision as the
          ! score
          if (survival_biasing) then
            ! We need to account for the fact that some weight was already
            ! absorbed
            score = (p % last_wgt + p % absorb_wgt) * flux
          else
            score = p % last_wgt * flux
          end if

        else
          if (i_nuclide > 0) then
            score = micro_xs(i_nuclide) % total * atom_density * flux
          else
            score = material_xs % total * flux
          end if
        end if


      case (SCORE_INVERSE_VELOCITY)
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! All events score to an inverse velocity bin. We actually use a
          ! collision estimator in place of an analog one since there is no way
          ! to count 'events' exactly for the inverse velocity
          if (survival_biasing) then
            ! We need to account for the fact that some weight was already
            ! absorbed
            score = p % last_wgt + p % absorb_wgt
          else
            score = p % last_wgt
          end if

          ! Score the flux weighted inverse velocity with velocity in units of
          ! cm/s
          score = score / material_xs % total &
               / (sqrt(TWO * E / MASS_NEUTRON_EV) * C_LIGHT * 100.0_8) * flux

        else
          ! For inverse velocity, we don't need a cross section. The velocity is
          ! in units of cm/s.
          score = flux / (sqrt(TWO * E / MASS_NEUTRON_EV) * C_LIGHT * 100.0_8)
        end if


      case (SCORE_SCATTER, SCORE_SCATTER_N)
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! Skip any event where the particle didn't scatter
          if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP
          ! Since only scattering events make it here, again we can use
          ! the weight entering the collision as the estimator for the
          ! reaction rate
          score = p % last_wgt * flux

        else
          ! Note SCORE_SCATTER_N not available for tracklength/collision.
          if (i_nuclide > 0) then
            score = (micro_xs(i_nuclide) % total &
                 - micro_xs(i_nuclide) % absorption) * atom_density * flux
          else
            score = (material_xs % total - material_xs % absorption) * flux
          end if
        end if


      case (SCORE_SCATTER_PN)
        ! Only analog estimators are available.
        ! Skip any event where the particle didn't scatter
        if (p % event /= EVENT_SCATTER) then
          i = i + t % moment_order(i)
          cycle SCORE_LOOP
        end if
        ! Since only scattering events make it here, again we can use
        ! the weight entering the collision as the estimator for the
        ! reaction rate
        score = p % last_wgt * flux


      case (SCORE_SCATTER_YN)
        ! Only analog estimators are available.
        ! Skip any event where the particle didn't scatter
        if (p % event /= EVENT_SCATTER) then
          i = i + (t % moment_order(i) + 1)**2 - 1
          cycle SCORE_LOOP
        end if
        ! Since only scattering events make it here, again we can use
        ! the weight entering the collision as the estimator for the
        ! reaction rate
        score = p % last_wgt * flux


      case (SCORE_NU_SCATTER, SCORE_NU_SCATTER_N)
        ! Only analog estimators are available.
        ! Skip any event where the particle didn't scatter
        if (p % event /= EVENT_SCATTER) cycle SCORE_LOOP
        ! For scattering production, we need to use the pre-collision weight
        ! times the yield as the estimate for the number of neutrons exiting a
        ! reaction with neutrons in the exit channel
        if (p % event_MT == ELASTIC .or. p % event_MT == N_LEVEL .or. &
             (p % event_MT >= N_N1 .and. p % event_MT <= N_NC)) then
          ! Don't waste time on very common reactions we know have
          ! multiplicities of one.
          score = p % last_wgt * flux
        else
          m = nuclides(p % event_nuclide) % reaction_index(p % event_MT)

          ! Get yield and apply to score
          associate (rxn => nuclides(p % event_nuclide) % reactions(m))
            score = p % last_wgt * flux &
                 * rxn % products(1) % yield % evaluate(E)
          end associate
        end if


      case (SCORE_NU_SCATTER_PN)
        ! Only analog estimators are available.
        ! Skip any event where the particle didn't scatter
        if (p % event /= EVENT_SCATTER) then
          i = i + t % moment_order(i)
          cycle SCORE_LOOP
        end if
        ! For scattering production, we need to use the pre-collision
        ! weight times the yield as the estimate for the number of
        ! neutrons exiting a reaction with neutrons in the exit channel
        if (p % event_MT == ELASTIC .or. p % event_MT == N_LEVEL .or. &
             (p % event_MT >= N_N1 .and. p % event_MT <= N_NC)) then
          ! Don't waste time on very common reactions we know have
          ! multiplicities of one.
          score = p % last_wgt * flux
        else
          m = nuclides(p % event_nuclide) % reaction_index(p % event_MT)

          ! Get yield and apply to score
          associate (rxn => nuclides(p % event_nuclide) % reactions(m))
            score = p % last_wgt * flux &
                 * rxn % products(1) % yield % evaluate(E)
          end associate
        end if


      case (SCORE_NU_SCATTER_YN)
        ! Only analog estimators are available.
        ! Skip any event where the particle didn't scatter
        if (p % event /= EVENT_SCATTER) then
          i = i + (t % moment_order(i) + 1)**2 - 1
          cycle SCORE_LOOP
        end if
        ! For scattering production, we need to use the pre-collision
        ! weight times the yield as the estimate for the number of
        ! neutrons exiting a reaction with neutrons in the exit channel
        if (p % event_MT == ELASTIC .or. p % event_MT == N_LEVEL .or. &
             (p % event_MT >= N_N1 .and. p % event_MT <= N_NC)) then
          ! Don't waste time on very common reactions we know have
          ! multiplicities of one.
          score = p % last_wgt * flux
        else
          m = nuclides(p % event_nuclide) % reaction_index(p % event_MT)

          ! Get yield and apply to score
          associate (rxn => nuclides(p%event_nuclide)%reactions(m))
            score = p % last_wgt * flux &
                 * rxn % products(1) % yield % evaluate(E)
          end associate
        end if


      case (SCORE_ABSORPTION)
        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No absorption events actually occur if survival biasing is on --
            ! just use weight absorbed in survival biasing
            score = p % absorb_wgt * flux
          else
            ! Skip any event where the particle wasn't absorbed
            if (p % event == EVENT_SCATTER) cycle SCORE_LOOP
            ! All fission and absorption events will contribute here, so we
            ! can just use the particle's weight entering the collision
            score = p % last_wgt * flux
          end if

        else
          if (i_nuclide > 0) then
            score = micro_xs(i_nuclide) % absorption * atom_density * flux
          else
            score = material_xs % absorption * flux
          end if
        end if


      case (SCORE_FISSION)
        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! fission
            if (micro_xs(p % event_nuclide) % absorption > ZERO) then
              score = p % absorb_wgt * micro_xs(p % event_nuclide) % fission &
                   / micro_xs(p % event_nuclide) % absorption * flux
            else
              score = ZERO
            end if
          else
            ! Skip any non-absorption events
            if (p % event == EVENT_SCATTER) cycle SCORE_LOOP
            ! All fission events will contribute, so again we can use
            ! particle's weight entering the collision as the estimate for the
            ! fission reaction rate
            score = p % last_wgt * micro_xs(p % event_nuclide) % fission &
                 / micro_xs(p % event_nuclide) % absorption * flux
          end if

        else
          if (i_nuclide > 0) then
            score = micro_xs(i_nuclide) % fission * atom_density * flux
          else
            score = material_xs % fission * flux
          end if
        end if


      case (SCORE_NU_FISSION)
        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing .or. p % fission) then
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              ! Normally, we only need to make contributions to one scoring
              ! bin. However, in the case of fission, since multiple fission
              ! neutrons were emitted with different energies, multiple
              ! outgoing energy bins may have been scored to. The following
              ! logic treats this special case and results to multiple bins
              call score_fission_eout(p, t, score_index, score_bin)
              cycle SCORE_LOOP
            end if
          end if
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! nu-fission
            if (micro_xs(p % event_nuclide) % absorption > ZERO) then
              score = p % absorb_wgt * micro_xs(p % event_nuclide) % &
                   nu_fission / micro_xs(p % event_nuclide) % absorption * flux
            else
              score = ZERO
            end if
          else
            ! Skip any non-fission events
            if (.not. p % fission) cycle SCORE_LOOP
            ! If there is no outgoing energy filter, than we only need to
            ! score to one bin. For the score to be 'analog', we need to
            ! score the number of particles that were banked in the fission
            ! bank. Since this was weighted by 1/keff, we multiply by keff
            ! to get the proper score.
            score = keff * p % wgt_bank * flux
          end if

        else
          if (i_nuclide > 0) then
            score = micro_xs(i_nuclide) % nu_fission * atom_density * flux
          else
            score = material_xs % nu_fission * flux
          end if
        end if


      case (SCORE_PROMPT_NU_FISSION)
        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing .or. p % fission) then
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              ! Normally, we only need to make contributions to one scoring
              ! bin. However, in the case of fission, since multiple fission
              ! neutrons were emitted with different energies, multiple
              ! outgoing energy bins may have been scored to. The following
              ! logic treats this special case and results to multiple bins
              call score_fission_eout(p, t, score_index, score_bin)
              cycle SCORE_LOOP
            end if
          end if
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! prompt-nu-fission
            if (micro_xs(p % event_nuclide) % absorption > ZERO) then
                score = p % absorb_wgt * micro_xs(p % event_nuclide) % fission &
                     * nuclides(p % event_nuclide) % nu(E, EMISSION_PROMPT) &
                     / micro_xs(p % event_nuclide) % absorption * flux
            else
              score = ZERO
            end if
          else
            ! Skip any non-fission events
            if (.not. p % fission) cycle SCORE_LOOP
            ! If there is no outgoing energy filter, than we only need to
            ! score to one bin. For the score to be 'analog', we need to
            ! score the number of particles that were banked in the fission
            ! bank as prompt neutrons. Since this was weighted by 1/keff, we
            ! multiply by keff to get the proper score.
            score = keff * p % wgt_bank * (ONE - sum(p % n_delayed_bank) &
                 / real(p % n_bank, 8)) * flux
          end if

        else
          if (i_nuclide > 0) then
              score = micro_xs(i_nuclide) % fission * nuclides(i_nuclide) % &
                   nu(E, EMISSION_PROMPT) * atom_density * flux
          else

            score = ZERO

            ! Loop over all nuclides in the current material
            if (p % material /= MATERIAL_VOID) then
              do l = 1, materials(p % material) % n_nuclides

                ! Get atom density
                atom_density_ = materials(p % material) % atom_density(l)

                ! Get index in nuclides array
                i_nuc = materials(p % material) % nuclide(l)

                ! Accumulate the contribution from each nuclide
                score = score + micro_xs(i_nuc) % fission * nuclides(i_nuc) % &
                     nu(E, EMISSION_PROMPT) * atom_density_ * flux
              end do
            end if
          end if
        end if

      case (SCORE_DELAYED_NU_FISSION)
        ! Set the delayedgroup filter index and the number of delayed group bins
        dg_filter = t % find_filter(FILTER_DELAYEDGROUP)

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing .or. p % fission) then
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              ! Normally, we only need to make contributions to one scoring
              ! bin. However, in the case of fission, since multiple fission
              ! neutrons were emitted with different energies, multiple
              ! outgoing energy bins may have been scored to. The following
              ! logic treats this special case and results to multiple bins
              call score_fission_eout(p, t, score_index, score_bin)
              cycle SCORE_LOOP
            end if
          end if
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! delayed-nu-fission
            if (micro_xs(p % event_nuclide) % absorption > ZERO .and. &
                 nuclides(p % event_nuclide) % fissionable) then

              ! Check if the delayed group filter is present
              if (dg_filter > 0) then
                select type(filt => filters(t % filter(dg_filter)) % obj)
                type is (DelayedGroupFilter)

                  ! Loop over all delayed group bins and tally to them
                  ! individually
                  do d_bin = 1, filt % n_bins

                    ! Get the delayed group for this bin
                    d = filt % groups(d_bin)

                    ! Compute the yield for this delayed group
                    yield = nuclides(p % event_nuclide) &
                            % nu(E, EMISSION_DELAYED, d)

                    ! Compute the score and tally to bin
                    score = p % absorb_wgt * yield &
                         * micro_xs(p % event_nuclide) % fission &
                         / micro_xs(p % event_nuclide) % absorption * flux
                    call score_fission_delayed_dg(t, d_bin, score, score_index)
                  end do
                  cycle SCORE_LOOP
                end select
              else
                ! If the delayed group filter is not present, compute the score
                ! by multiplying the absorbed weight by the fraction of the
                ! delayed-nu-fission xs to the absorption xs
                score = p % absorb_wgt * micro_xs(p % event_nuclide) % fission &
                     * nuclides(p % event_nuclide) % nu(E, EMISSION_DELAYED) &
                     / micro_xs(p % event_nuclide) % absorption * flux
              end if
            end if
          else
            ! Skip any non-fission events
            if (.not. p % fission) cycle SCORE_LOOP
            ! If there is no outgoing energy filter, than we only need to
            ! score to one bin. For the score to be 'analog', we need to
            ! score the number of particles that were banked in the fission
            ! bank. Since this was weighted by 1/keff, we multiply by keff
            ! to get the proper score. Loop over the neutrons produced from
            ! fission and check which ones are delayed. If a delayed neutron is
            ! encountered, add its contribution to the fission bank to the
            ! score.

            ! Check if the delayed group filter is present
            if (dg_filter > 0) then
              select type(filt => filters(t % filter(dg_filter)) % obj)
              type is (DelayedGroupFilter)

                ! Loop over all delayed group bins and tally to them
                ! individually
                do d_bin = 1, filt % n_bins

                  ! Get the delayed group for this bin
                  d = filt % groups(d_bin)

                  ! Compute the score and tally to bin
                  score = keff * p % wgt_bank / p % n_bank * &
                       p % n_delayed_bank(d) * flux

                  call score_fission_delayed_dg(t, d_bin, score, score_index)
                end do
                cycle SCORE_LOOP
              end select
            else

              ! Add the contribution from all delayed groups
              score = keff * p % wgt_bank / p % n_bank * &
                   sum(p % n_delayed_bank) * flux
            end if
          end if
        else

          ! Check if tally is on a single nuclide
          if (i_nuclide > 0) then

            ! Check if the delayed group filter is present
            if (dg_filter > 0) then
              select type(filt => filters(t % filter(dg_filter)) % obj)
              type is (DelayedGroupFilter)

                ! Loop over all delayed group bins and tally to them
                ! individually
                do d_bin = 1, filt % n_bins

                  ! Get the delayed group for this bin
                  d = filt % groups(d_bin)

                  ! Compute the yield for this delayed group
                  yield = nuclides(i_nuclide) % nu(E, EMISSION_DELAYED, d)

                  ! Compute the score and tally to bin
                  score = micro_xs(i_nuclide) % fission * yield * &
                       atom_density * flux
                  call score_fission_delayed_dg(t, d_bin, score, score_index)
                end do
                cycle SCORE_LOOP
              end select
            else

              ! If the delayed group filter is not present, compute the score
              ! by multiplying the delayed-nu-fission macro xs by the flux
              score = micro_xs(i_nuclide) % fission * nuclides(i_nuclide) % &
                   nu(E, EMISSION_DELAYED) * atom_density * flux
            end if

          ! Tally is on total nuclides
          else

            ! Check if the delayed group filter is present
            if (dg_filter > 0) then
              select type(filt => filters(t % filter(dg_filter)) % obj)
              type is (DelayedGroupFilter)

                ! Loop over all nuclides in the current material
                if (p % material /= MATERIAL_VOID) then
                  do l = 1, materials(p % material) % n_nuclides

                    ! Get atom density
                    atom_density_ = materials(p % material) % atom_density(l)

                    ! Get index in nuclides array
                    i_nuc = materials(p % material) % nuclide(l)

                    ! Loop over all delayed group bins and tally to them
                    ! individually
                    do d_bin = 1, filt % n_bins

                      ! Get the delayed group for this bin
                      d = filt % groups(d_bin)

                      ! Get the yield for the desired nuclide and delayed group
                      yield = nuclides(i_nuc) % nu(E, EMISSION_DELAYED, d)

                      ! Compute the score and tally to bin
                      score = micro_xs(i_nuc) % fission * yield &
                           * atom_density_ * flux
                      call score_fission_delayed_dg(t, d_bin, score, &
                           score_index)
                    end do
                  end do
                end if
                cycle SCORE_LOOP
              end select
            else

              score = ZERO

              ! Loop over all nuclides in the current material
              if (p % material /= MATERIAL_VOID) then
                do l = 1, materials(p % material) % n_nuclides

                  ! Get atom density
                  atom_density_ = materials(p % material) % atom_density(l)

                  ! Get index in nuclides array
                  i_nuc = materials(p % material) % nuclide(l)

                  ! Accumulate the contribution from each nuclide
                  score = score + micro_xs(i_nuc) % fission * nuclides(i_nuc) %&
                       nu(E, EMISSION_DELAYED) * atom_density_ * flux
                end do
              end if
            end if
          end if
        end if

      case (SCORE_DECAY_RATE)
        ! Set the delayedgroup filter index
        dg_filter = t % find_filter(FILTER_DELAYEDGROUP)

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! delayed-nu-fission
            if (micro_xs(p % event_nuclide) % absorption > ZERO .and. &
                 nuclides(p % event_nuclide) % fissionable) then

              ! Check if the delayed group filter is present
              if (dg_filter > 0) then
                select type(filt => filters(t % filter(dg_filter)) % obj)
                type is (DelayedGroupFilter)

                  ! Loop over all delayed group bins and tally to them
                  ! individually
                  do d_bin = 1, filt % n_bins

                    ! Get the delayed group for this bin
                    d = filt % groups(d_bin)

                    ! Compute the yield for this delayed group
                    yield = nuclides(p % event_nuclide) &
                         % nu(E, EMISSION_DELAYED, d)

                    associate (rxn => nuclides(p % event_nuclide) % &
                         reactions(nuclides(p % event_nuclide) % &
                         index_fission(1)))

                      ! Compute the score
                      score = p % absorb_wgt * yield * &
                           micro_xs(p % event_nuclide) % fission &
                           / micro_xs(p % event_nuclide) % absorption &
                           * rxn % products(1 + d) % decay_rate * flux
                    end associate

                    ! Tally to bin
                    call score_fission_delayed_dg(t, d_bin, score, score_index)
                  end do
                  cycle SCORE_LOOP
                end select
              else

                ! If the delayed group filter is not present, compute the score
                ! by accumulating the absorbed weight times the decay rate times
                ! the fraction of the delayed-nu-fission xs to the absorption xs
                ! for all delayed groups.
                score = ZERO

                associate (rxn => nuclides(p % event_nuclide) % &
                     reactions(nuclides(p % event_nuclide) % index_fission(1)))

                  ! We need to be careful not to overshoot the number of delayed
                  ! groups since this could cause the range of the
                  ! rxn % products array to be exceeded. Hence, we use the size
                  ! of this array and not the MAX_DELAYED_GROUPS constant for
                  ! this loop.
                  do d = 1, size(rxn % products) - 2

                    score = score + rxn % products(1 + d) % decay_rate * &
                         p % absorb_wgt &
                         * micro_xs(p % event_nuclide) % fission &
                         * nuclides(p % event_nuclide) % &
                         nu(E, EMISSION_DELAYED, d) &
                         / micro_xs(p % event_nuclide) % absorption * flux
                  end do
                end associate
              end if
            end if
          else

            ! Skip any non-fission events
            if (.not. p % fission) cycle SCORE_LOOP
            ! If there is no outgoing energy filter, than we only need to
            ! score to one bin. For the score to be 'analog', we need to
            ! score the number of particles that were banked in the fission
            ! bank. Since this was weighted by 1/keff, we multiply by keff
            ! to get the proper score. Loop over the neutrons produced from
            ! fission and check which ones are delayed. If a delayed neutron is
            ! encountered, add its contribution to the fission bank to the
            ! score.

            score = ZERO

            ! loop over number of particles banked
            do k = 1, p % n_bank

              ! get the delayed group
              g = fission_bank(n_bank - p % n_bank + k) % delayed_group

              ! Case for tallying delayed emissions
              if (g /= 0) then

                ! Accumulate the decay rate times delayed nu fission score
                associate (rxn => nuclides(p % event_nuclide) % &
                     reactions(nuclides(p % event_nuclide) % index_fission(1)))

                  ! determine score based on bank site weight and keff.
                  score = score + keff * fission_bank(n_bank - p % n_bank + k) &
                       % wgt * rxn % products(1 + g) % decay_rate * flux
                end associate

                ! if the delayed group filter is present, tally to corresponding
                ! delayed group bin if it exists
                if (dg_filter > 0) then

                  ! declare the delayed group filter type
                  select type(filt => filters(t % filter(dg_filter)) % obj)
                  type is (DelayedGroupFilter)

                    ! loop over delayed group bins until the corresponding bin
                    ! is found
                    do d_bin = 1, filt % n_bins
                      d = filt % groups(d_bin)

                      ! check whether the delayed group of the particle is equal
                      ! to the delayed group of this bin
                      if (d == g) then
                        call score_fission_delayed_dg(t, d_bin, score, &
                                                      score_index)
                      end if
                    end do
                  end select

                  ! Reset the score to zero
                  score = ZERO
                end if
              end if
            end do
          end if
        else

          ! Check if tally is on a single nuclide
          if (i_nuclide > 0) then

            ! Check if the delayed group filter is present
            if (dg_filter > 0) then
              select type(filt => filters(t % filter(dg_filter)) % obj)
              type is (DelayedGroupFilter)

                ! Loop over all delayed group bins and tally to them
                ! individually
                do d_bin = 1, filt % n_bins

                  ! Get the delayed group for this bin
                  d = filt % groups(d_bin)

                  ! Compute the yield for this delayed group
                  yield = nuclides(i_nuclide) % nu(E, EMISSION_DELAYED, d)

                  associate (rxn => nuclides(i_nuclide) % &
                       reactions(nuclides(i_nuclide) % index_fission(1)))

                    ! Compute the score and tally to bin
                    score = micro_xs(i_nuclide) % fission * yield * flux * &
                         atom_density * rxn % products(1 + d) % decay_rate
                  end associate

                  ! Tally to bin
                  call score_fission_delayed_dg(t, d_bin, score, score_index)
                end do
                cycle SCORE_LOOP
              end select
            else

              ! If the delayed group filter is not present, compute the score
              ! by accumulating the absorbed weight times the decay rate times
              ! the fraction of the delayed-nu-fission xs to the absorption xs
              ! for all delayed groups.
              score = ZERO

              associate (rxn => nuclides(p % event_nuclide) % &
                   reactions(nuclides(p % event_nuclide) % index_fission(1)))

                ! We need to be careful not to overshoot the number of delayed
                ! groups since this could cause the range of the rxn % products
                ! array to be exceeded. Hence, we use the size of this array
                ! and not the MAX_DELAYED_GROUPS constant for this loop.
                do d = 1, size(rxn % products) - 2

                  score = score + micro_xs(i_nuclide) % fission * flux * &
                       nuclides(i_nuclide) % nu(E, EMISSION_DELAYED) * &
                       atom_density * rxn % products(1 + d) % decay_rate
                end do
              end associate
            end if

            ! Tally is on total nuclides
          else

            ! Check if the delayed group filter is present
            if (dg_filter > 0) then
              select type(filt => filters(t % filter(dg_filter)) % obj)
              type is (DelayedGroupFilter)

                ! Loop over all nuclides in the current material
                if (p % material /= MATERIAL_VOID) then
                  do l = 1, materials(p % material) % n_nuclides

                    ! Get atom density
                    atom_density_ = materials(p % material) % atom_density(l)

                    ! Get index in nuclides array
                    i_nuc = materials(p % material) % nuclide(l)

                    if (nuclides(i_nuc) % fissionable) then

                      ! Loop over all delayed group bins and tally to them
                      ! individually
                      do d_bin = 1, filt % n_bins

                        ! Get the delayed group for this bin
                        d = filt % groups(d_bin)

                        ! Get the yield for the desired nuclide and delayed
                        ! group
                        yield = nuclides(i_nuc) % nu(E, EMISSION_DELAYED, d)

                        associate (rxn => nuclides(i_nuc) % &
                             reactions(nuclides(i_nuc) % index_fission(1)))

                          ! Compute the score
                          score = micro_xs(i_nuc) % fission * yield * flux * &
                               atom_density_ &
                               * rxn % products(1 + d) % decay_rate
                        end associate

                        ! Tally to bin
                        call score_fission_delayed_dg(t, d_bin, score, &
                                                      score_index)
                      end do
                    end if
                  end do
                end if
                cycle SCORE_LOOP
              end select
            else

              score = ZERO

              ! Loop over all nuclides in the current material
              if (p % material /= MATERIAL_VOID) then
                do l = 1, materials(p % material) % n_nuclides

                  ! Get atom density
                  atom_density_ = materials(p % material) % atom_density(l)

                  ! Get index in nuclides array
                  i_nuc = materials(p % material) % nuclide(l)

                  if (nuclides(i_nuc) % fissionable) then

                    associate (rxn => nuclides(i_nuc) % &
                         reactions(nuclides(i_nuc) % index_fission(1)))

                      ! We need to be careful not to overshoot the number of
                      ! delayed groups since this could cause the range of the
                      ! rxn % products array to be exceeded. Hence, we use the
                      ! size of this array and not the MAX_DELAYED_GROUPS
                      ! constant for this loop.
                      do d = 1, size(rxn % products) - 2

                        ! Accumulate the contribution from each nuclide
                        score = score + micro_xs(i_nuc) % fission &
                             * nuclides(i_nuc) % nu(E, EMISSION_DELAYED) &
                             * atom_density_ * flux &
                             * rxn % products(1 + d) % decay_rate
                      end do
                    end associate
                  end if
                end do
              end if
            end if
          end if
        end if

      case (SCORE_KAPPA_FISSION)
        ! Determine kappa-fission cross section on the fly. The ENDF standard
        ! (ENDF-102) states that MT 18 stores the fission energy as the Q_value
        ! (fission(1))

        score = ZERO

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! fission scaled by kappa-fission
            associate (nuc => nuclides(p % event_nuclide))
              if (micro_xs(p % event_nuclide) % absorption > ZERO .and. &
                   nuc % fissionable) then
                score = p % absorb_wgt * &
                     nuc % reactions(nuc % index_fission(1)) % Q_value * &
                     micro_xs(p % event_nuclide) % fission / &
                     micro_xs(p % event_nuclide) % absorption * flux
              end if
            end associate
          else
            ! Skip any non-absorption events
            if (p % event == EVENT_SCATTER) cycle SCORE_LOOP
            ! All fission events will contribute, so again we can use
            ! particle's weight entering the collision as the estimate for
            ! the fission energy production rate
            associate (nuc => nuclides(p % event_nuclide))
              if (nuc % fissionable) then
                score = p % last_wgt * &
                     nuc % reactions(nuc % index_fission(1)) % Q_value * &
                     micro_xs(p % event_nuclide) % fission / &
                     micro_xs(p % event_nuclide) % absorption * flux
              end if
            end associate
          end if

        else
          if (i_nuclide > 0) then
            associate (nuc => nuclides(i_nuclide))
              if (nuc % fissionable) then
                score = nuc % reactions(nuc % index_fission(1)) % Q_value * &
                     micro_xs(i_nuclide) % fission * atom_density * flux
              end if
            end associate
          else
            if (p % material == MATERIAL_VOID) then
              score = ZERO
            else
              do l = 1, materials(p % material) % n_nuclides
                ! Determine atom density and index of nuclide
                atom_density_ = materials(p % material) % atom_density(l)
                i_nuc = materials(p % material) % nuclide(l)

                ! If nuclide is fissionable, accumulate kappa fission
                associate(nuc => nuclides(i_nuc))
                  if (nuc % fissionable) then
                    score = score + &
                         nuc % reactions(nuc % index_fission(1)) % Q_value * &
                         micro_xs(i_nuc) % fission * atom_density_ * flux
                  end if
                end associate
              end do
            end if
          end if
        end if

      case (SCORE_EVENTS)
        ! Simply count number of scoring events
        score = ONE

      case (ELASTIC)
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! Check if event MT matches
          if (p % event_MT /= ELASTIC) cycle SCORE_LOOP
          score = p % last_wgt * flux

        else
          if (i_nuclide > 0) then
            if (micro_xs(i_nuclide) % elastic == CACHE_INVALID) then
              call calculate_elastic_xs(i_nuclide)
            end if
            score = micro_xs(i_nuclide) % elastic * atom_density * flux
          else
            score = ZERO
            if (p % material /= MATERIAL_VOID) then
              do l = 1, materials(p % material) % n_nuclides
                ! Get atom density
                atom_density_ = materials(p % material) % atom_density(l)

                ! Get index in nuclides array
                i_nuc = materials(p % material) % nuclide(l)
                if (micro_xs(i_nuc) % elastic == CACHE_INVALID) then
                  call calculate_elastic_xs(i_nuc)
                end if

                score = score + micro_xs(i_nuc) % elastic * atom_density_ * flux
              end do
            end if
          end if
        end if

      case (SCORE_FISS_Q_PROMPT)
        score = ZERO

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! fission scaled by Q-value
            associate (nuc => nuclides(p % event_nuclide))
              if (micro_xs(p % event_nuclide) % absorption > ZERO .and. &
                   allocated(nuc % fission_q_prompt)) then
                score = p % absorb_wgt &
                     * nuc % fission_q_prompt % evaluate(E) &
                     * micro_xs(p % event_nuclide) % fission &
                     / micro_xs(p % event_nuclide) % absorption * flux
              end if
            end associate
          else
            ! Skip any non-absorption events
            if (p % event == EVENT_SCATTER) cycle SCORE_LOOP
            ! All fission events will contribute, so again we can use
            ! particle's weight entering the collision as the estimate for
            ! the fission energy production rate
            associate (nuc => nuclides(p % event_nuclide))
              if (allocated(nuc % fission_q_prompt)) then
                score = p % last_wgt &
                     * nuc % fission_q_prompt % evaluate(E) &
                     * micro_xs(p % event_nuclide) % fission &
                     / micro_xs(p % event_nuclide) % absorption * flux
              end if
            end associate
          end if

        else
          if (i_nuclide > 0) then
            if (allocated(nuclides(i_nuclide) % fission_q_prompt)) then
              score = micro_xs(i_nuclide) % fission * atom_density * flux &
                      * nuclides(i_nuclide) % fission_q_prompt % evaluate(E)
            end if
          else
            if (p % material /= MATERIAL_VOID) then
              do l = 1, materials(p % material) % n_nuclides
                atom_density_ = materials(p % material) % atom_density(l)
                i_nuc = materials(p % material) % nuclide(l)
                if (allocated(nuclides(i_nuc) % fission_q_prompt)) then
                  score = score + micro_xs(i_nuc) % fission * atom_density_ &
                          * flux &
                          * nuclides(i_nuc) % fission_q_prompt % evaluate(E)
                end if
              end do
            end if
          end if
        end if

      case (SCORE_FISS_Q_RECOV)
        score = ZERO

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! fission scaled by Q-value
            associate (nuc => nuclides(p % event_nuclide))
              if (micro_xs(p % event_nuclide) % absorption > ZERO .and. &
                   allocated(nuc % fission_q_recov)) then
                score = p % absorb_wgt &
                     * nuc % fission_q_recov % evaluate(E) &
                     * micro_xs(p % event_nuclide) % fission &
                     / micro_xs(p % event_nuclide) % absorption * flux
              end if
            end associate
          else
            ! Skip any non-absorption events
            if (p % event == EVENT_SCATTER) cycle SCORE_LOOP
            ! All fission events will contribute, so again we can use
            ! particle's weight entering the collision as the estimate for
            ! the fission energy production rate
            associate (nuc => nuclides(p % event_nuclide))
              if (allocated(nuc % fission_q_recov)) then
                score = p % last_wgt &
                     * nuc % fission_q_recov % evaluate(E) &
                     * micro_xs(p % event_nuclide) % fission &
                     / micro_xs(p % event_nuclide) % absorption * flux
              end if
            end associate
          end if

        else
          if (i_nuclide > 0) then
            if (allocated(nuclides(i_nuclide) % fission_q_recov)) then
              score = micro_xs(i_nuclide) % fission * atom_density * flux &
                      * nuclides(i_nuclide) % fission_q_recov % evaluate(E)
            end if
          else
            if (p % material /= MATERIAL_VOID) then
              do l = 1, materials(p % material) % n_nuclides
                atom_density_ = materials(p % material) % atom_density(l)
                i_nuc = materials(p % material) % nuclide(l)
                if (allocated(nuclides(i_nuc) % fission_q_recov)) then
                  score = score + micro_xs(i_nuc) % fission * atom_density_ &
                          * flux &
                          * nuclides(i_nuc) % fission_q_recov % evaluate(E)
                end if
              end do
            end if
          end if
        end if

      case (N_2N, N_3N, N_4N, N_GAMMA, N_P, N_A)
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! Check if event MT matches
          if (p % event_MT /= score_bin) cycle SCORE_LOOP
          score = p % last_wgt * flux

        else
          ! Determine index in NuclideMicroXS % reaction array
          select case (score_bin)
          case (N_2N)
            m = 1
          case (N_3N)
            m = 2
          case (N_4N)
            m = 3
          case (N_GAMMA)
            m = 4
          case (N_P)
            m = 5
          case (N_A)
            m = 6
          end select

          if (i_nuclide > 0) then
            score = micro_xs(i_nuclide) % reaction(m) * atom_density * flux
          else
            score = ZERO
            if (p % material /= MATERIAL_VOID) then
              associate (mat => materials(p % material))
                do l = 1, materials(p % material) % n_nuclides
                  i_nuc = mat % nuclide(l)
                  atom_density_ = mat % atom_density(l)
                  score = score + micro_xs(i_nuc) % reaction(m) * atom_density_ * flux
                end do
              end associate
            end if
          end if
        end if

      case default
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! Any other score is assumed to be a MT number. Thus, we just need
          ! to check if it matches the MT number of the event
          if (p % event_MT /= score_bin) cycle SCORE_LOOP
          score = p % last_wgt * flux

        else
          ! Any other cross section has to be calculated on-the-fly. For
          ! cross sections that are used often (e.g. n2n, ngamma, etc. for
          ! depletion), it might make sense to optimize this section or
          ! pre-calculate cross sections
          if (score_bin > 1) then
            ! Set default score
            score = ZERO

            if (i_nuclide > 0) then
              m = nuclides(i_nuclide) % reaction_index(score_bin)
              if (m /= 0) then
                ! Retrieve temperature and energy grid index and interpolation
                ! factor
                i_temp = micro_xs(i_nuclide) % index_temp
                if (i_temp > 0) then
                  i_energy = micro_xs(i_nuclide) % index_grid
                  f = micro_xs(i_nuclide) % interp_factor

                  associate (xs => nuclides(i_nuclide) % reactions(m) &
                             % xs(i_temp))
                    if (i_energy >= xs % threshold) then
                      score = ((ONE - f) * xs % value(i_energy - &
                           xs % threshold + 1) + f * xs % value(i_energy - &
                           xs % threshold + 2)) * atom_density * flux
                    end if
                  end associate
                else
                  ! This block is reached if multipole is turned on and we're in
                  ! the resolved range. Assume xs is zero.
                  score = ZERO
                end if
              end if

            else
              if (p % material /= MATERIAL_VOID) then
                do l = 1, materials(p % material) % n_nuclides
                  ! Get atom density
                  atom_density_ = materials(p % material) % atom_density(l)

                  ! Get index in nuclides array
                  i_nuc = materials(p % material) % nuclide(l)

                  m = nuclides(i_nuc) % reaction_index(score_bin)
                  if (m /= 0) then
                    ! Retrieve temperature and energy grid index and
                    ! interpolation factor
                    i_temp = micro_xs(i_nuc) % index_temp
                    if (i_temp > 0) then
                      i_energy = micro_xs(i_nuc) % index_grid
                      f = micro_xs(i_nuc) % interp_factor

                      associate (xs => nuclides(i_nuc) % reactions(m) &
                                 % xs(i_temp))
                        if (i_energy >= xs % threshold) then
                          score = score + ((ONE - f) * xs % value(i_energy - &
                               xs % threshold + 1) + f * xs % value(i_energy - &
                               xs % threshold + 2)) * atom_density_ * flux
                        end if
                      end associate
                    else
                      ! This block is reached if multipole is turned on and
                      ! we're in the resolved range. Assume xs is zero.
                      score = ZERO
                    end if
                  end if
                end do
              end if
            end if

          else
            call fatal_error("Invalid score type on tally " &
                 // to_str(t % id) // ".")
          end if
        end if

      end select

      !#########################################################################
      ! Add derivative information on score for differential tallies.

      if (t % deriv /= NONE) then
        call apply_derivative_to_score(p, t, i_nuclide, atom_density, &
             score_bin, score)
      end if

      !#########################################################################
      ! Expand score if necessary and add to tally results.
      call expand_and_score(p, t, score_index, filter_index, score_bin, &
                            score, i)

    end do SCORE_LOOP
  end subroutine score_general_ce

  subroutine score_general_mg(p, t, start_index, filter_index, i_nuclide, &
       atom_density, flux)
    type(Particle),    intent(in)    :: p
    type(TallyObject), intent(inout) :: t
    integer,           intent(in)    :: start_index
    integer,           intent(in)    :: i_nuclide
    integer,           intent(in)    :: filter_index   ! for % results
    real(8),           intent(in)    :: flux           ! flux estimate
    real(8),           intent(in)    :: atom_density   ! atom/b-cm

    integer :: i                    ! loop index for scoring bins
    integer :: q                    ! loop index for scoring bins
    integer :: score_bin            ! scoring bin, e.g. SCORE_FLUX
    integer :: score_index          ! scoring bin index
    integer :: d                    ! delayed neutron index
    integer :: g                    ! delayed neutron index
    integer :: k                    ! loop index for bank sites
    integer :: d_bin                ! delayed group bin index
    integer :: dg_filter            ! index of delayed group filter
    real(8) :: score                ! analog tally score
    real(8) :: p_uvw(3)             ! Particle's current uvw
    integer :: p_g                  ! Particle group to use for getting info
                                    ! to tally with.
    class(Mgxs), pointer :: matxs
    class(Mgxs), pointer :: nucxs

    ! Set the direction and group to use with get_xs
    if (t % estimator == ESTIMATOR_ANALOG .or. &
         t % estimator == ESTIMATOR_COLLISION) then

      if (survival_biasing) then

        ! Then we either are alive and had a scatter (and so g changed),
        ! or are dead and g did not change
        if (p % alive) then
          p_uvw = p % last_uvw
          p_g = p % last_g
        else
          p_uvw = p % coord(p % n_coord) % uvw
          p_g = p % g
        end if
      else if (p % event == EVENT_SCATTER) then

        ! Then the energy group has been changed by the scattering routine
        ! meaning gin is now in p % last_g
        p_uvw = p % last_uvw
        p_g = p % last_g
      else

        ! No scatter, no change in g.
        p_uvw = p % coord(p % n_coord) % uvw
        p_g = p % g
      end if
    else

      ! No actual collision so g has not changed.
      p_uvw = p % coord(p % n_coord) % uvw
      p_g = p % g
    end if

    ! To significantly reduce de-referencing, point matxs to the
    ! macroscopic Mgxs for the material of interest
    matxs => macro_xs(p % material) % obj

    ! Do same for nucxs, point it to the microscopic nuclide data of interest
    if (i_nuclide > 0) then
      nucxs => nuclides_MG(i_nuclide) % obj
      ! And since we haven't calculated this temperature index yet, do so now
      call nucxs % find_temperature(p % sqrtkT)
    end if

    i = 0
    SCORE_LOOP: do q = 1, t % n_user_score_bins
      i = i + 1

      ! determine what type of score bin
      score_bin = t % score_bins(i)

      ! determine scoring bin index
      score_index = start_index + i

      !#########################################################################
      ! Determine appropirate scoring value.

      select case(score_bin)


      case (SCORE_FLUX, SCORE_FLUX_YN)
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! All events score to a flux bin. We actually use a collision
          ! estimator in place of an analog one since there is no way to count
          ! 'events' exactly for the flux

          if (survival_biasing) then
            ! We need to account for the fact that some weight was already
            ! absorbed
            score = p % last_wgt + p % absorb_wgt
          else
            score = p % last_wgt
          end if
          score = score / material_xs % total * flux

        else
          ! For flux, we need no cross section
          score = flux
        end if


      case (SCORE_TOTAL, SCORE_TOTAL_YN)
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! All events will score to the total reaction rate. We can just
          ! use the weight of the particle entering the collision as the
          ! score
          if (survival_biasing) then
            ! We need to account for the fact that some weight was already
            ! absorbed
            score = p % last_wgt + p % absorb_wgt
          else
            score = p % last_wgt
          end if

          if (i_nuclide > 0) then
            score = score * atom_density * &
                 nucxs % get_xs('total', p_g, UVW=p_uvw) / &
                 matxs % get_xs('total', p_g, UVW=p_uvw) * flux
          end if

        else
          if (i_nuclide > 0) then
            score = nucxs % get_xs('total', p_g, UVW=p_uvw) * &
                 atom_density * flux
          else
            score = material_xs % total * flux
          end if
        end if


      case (SCORE_INVERSE_VELOCITY)
        if (t % estimator == ESTIMATOR_ANALOG .or. &
             t % estimator == ESTIMATOR_COLLISION) then
          ! All events score to an inverse velocity bin. We actually use a
          ! collision estimator in place of an analog one since there is no way
          ! to count 'events' exactly for the inverse velocity
          if (survival_biasing) then
            ! We need to account for the fact that some weight was already
            ! absorbed
            score = p % last_wgt + p % absorb_wgt
          else
            score = p % last_wgt
          end if

          if (i_nuclide > 0) then
            score = score * nucxs % get_xs('inverse-velocity', p_g, UVW=p_uvw) &
                 / matxs % get_xs('absorption', p_g, UVW=p_uvw) * flux
          else
            score = score * matxs % get_xs('inverse-velocity', p_g, UVW=p_uvw) &
                 / matxs % get_xs('absorption', p_g, UVW=p_uvw) * flux
          end if

        else

          if (i_nuclide > 0) then
            score = flux * nucxs % get_xs('inverse-velocity', p_g, UVW=p_uvw)
          else
            score = flux * matxs % get_xs('inverse-velocity', p_g, UVW=p_uvw)
          end if
        end if


      case (SCORE_SCATTER, SCORE_SCATTER_N, SCORE_SCATTER_PN, SCORE_SCATTER_YN)
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! Skip any event where the particle didn't scatter
          if (p % event /= EVENT_SCATTER) then
            if (score_bin == SCORE_SCATTER_PN) then
              i = i + t % moment_order(i)
            else if (score_bin == SCORE_SCATTER_YN) then
              i = i + (t % moment_order(i) + 1)**2 - 1
            end if
            cycle SCORE_LOOP
          end if

          ! Since only scattering events make it here, again we can use
          ! the weight entering the collision as the estimator for the
          ! reaction rate
          score = p % last_wgt * flux

          ! Since we transport based on material data, the angle selected
          ! was not selected from the f(mu) for the nuclide.  Therefore
          ! adjust the score by the actual probability for that nuclide.
          if (i_nuclide > 0) then
            score = score * atom_density * &
                 nucxs % get_xs('scatter*f_mu/mult', p % last_g, p % g, &
                                UVW=p_uvw, MU=p % mu) / &
                 matxs % get_xs('scatter*f_mu/mult', p % last_g, p % g, &
                                UVW=p_uvw, MU=p % mu)
          end if

        else
          ! Note SCORE_SCATTER_*N not available for tracklength/collision.
          if (i_nuclide > 0) then
            score = atom_density * flux * &
                 nucxs % get_xs('scatter/mult', p_g, UVW=p_uvw)
          else
            ! Get the scattering x/s and take away
            ! the multiplication baked in to sigS
            score = flux * &
                 matxs % get_xs('scatter/mult', p_g, UVW=p_uvw)
          end if
        end if


      case (SCORE_NU_SCATTER, SCORE_NU_SCATTER_N, SCORE_NU_SCATTER_PN, &
            SCORE_NU_SCATTER_YN)
        if (t % estimator == ESTIMATOR_ANALOG) then
          ! Skip any event where the particle didn't scatter
          if (p % event /= EVENT_SCATTER) then
            if (score_bin == SCORE_NU_SCATTER_PN) then
              i = i + t % moment_order(i)
            else if (score_bin == SCORE_NU_SCATTER_YN) then
              i = i + (t % moment_order(i) + 1)**2 - 1
            end if
            cycle SCORE_LOOP
          end if

          ! For scattering production, we need to use the pre-collision
          ! weight times the multiplicity as the estimate for the number of
          ! neutrons exiting a reaction with neutrons in the exit channel
          score = p % wgt * flux

          ! Since we transport based on material data, the angle selected
          ! was not selected from the f(mu) for the nuclide.  Therefore
          ! adjust the score by the actual probability for that nuclide.
          if (i_nuclide > 0) then
            score = score * atom_density * &
                 nucxs % get_xs('scatter*f_mu', p % last_g, p % g, &
                                UVW=p_uvw, MU=p % mu) / &
                 matxs % get_xs('scatter*f_mu', p % last_g, p % g, &
                                UVW=p_uvw, MU=p % mu)
          end if

        else
          ! Note SCORE_NU_SCATTER_*N not available for tracklength/collision.
          if (i_nuclide > 0) then
              score = nucxs % get_xs('scatter', p_g, UVW=p_uvw) * &
                   atom_density * flux
          else
            ! Get the scattering x/s, which includes multiplication
            score = matxs % get_xs('scatter', p_g, UVW=p_uvw) * flux
          end if
        end if


      case (SCORE_ABSORPTION)
        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No absorption events actually occur if survival biasing is on --
            ! just use weight absorbed in survival biasing
            score = p % absorb_wgt * flux
          else
            ! Skip any event where the particle wasn't absorbed
            if (p % event == EVENT_SCATTER) cycle SCORE_LOOP
            ! All fission and absorption events will contribute here, so we
            ! can just use the particle's weight entering the collision
            score = p % last_wgt * flux
          end if
          if (i_nuclide > 0) then
            score = score * atom_density * &
                 nucxs % get_xs('absorption', p_g, UVW=p_uvw) / &
                 matxs % get_xs('absorption', p_g, UVW=p_uvw)
          end if
        else
          if (i_nuclide > 0) then
            score = nucxs % get_xs('absorption', p_g, UVW=p_uvw) * &
                 atom_density * flux
          else
            score = material_xs % absorption * flux
          end if
        end if


      case (SCORE_FISSION)

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! fission
            score = p % absorb_wgt * flux
          else
            ! Skip any non-absorption events
            if (p % event == EVENT_SCATTER) cycle SCORE_LOOP
            ! All fission events will contribute, so again we can use
            ! particle's weight entering the collision as the estimate for the
            ! fission reaction rate
            score = p % last_wgt * flux
          end if
          if (i_nuclide > 0) then
            score = score * atom_density * &
                 nucxs % get_xs('fission', p_g, UVW=p_uvw) / &
                 matxs % get_xs('absorption', p_g, UVW=p_uvw)
          else
            score = score * &
                 matxs % get_xs('fission', p_g, UVW=p_uvw) / &
                 matxs % get_xs('absorption', p_g, UVW=p_uvw)
          end if
        else
          if (i_nuclide > 0) then
            score = nucxs % get_xs('fission', p_g, UVW=p_uvw) * &
                 atom_density * flux
          else
            score = matxs % get_xs('fission', p_g, UVW=p_uvw) * flux
          end if
        end if


      case (SCORE_NU_FISSION)

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing .or. p % fission) then
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              ! Normally, we only need to make contributions to one scoring
              ! bin. However, in the case of fission, since multiple fission
              ! neutrons were emitted with different energies, multiple
              ! outgoing energy bins may have been scored to. The following
              ! logic treats this special case and results to multiple bins
              call score_fission_eout(p, t, score_index, score_bin)
              cycle SCORE_LOOP
            end if
          end if
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! nu-fission
            score = p % absorb_wgt * flux
            if (i_nuclide > 0) then
              score = score * atom_density * &
                   nucxs % get_xs('nu-fission', p_g, UVW=p_uvw) / &
                   matxs % get_xs('absorption', p_g, UVW=p_uvw)
            else
              score = score * &
                   matxs % get_xs('nu-fission', p_g, UVW=p_uvw) / &
                   matxs % get_xs('absorption', p_g, UVW=p_uvw)
            end if
          else
            ! Skip any non-fission events
            if (.not. p % fission) cycle SCORE_LOOP
            ! If there is no outgoing energy filter, than we only need to
            ! score to one bin. For the score to be 'analog', we need to
            ! score the number of particles that were banked in the fission
            ! bank. Since this was weighted by 1/keff, we multiply by keff
            ! to get the proper score.
            score = keff * p % wgt_bank * flux
            if (i_nuclide > 0) then
              score = score * atom_density * &
                   nucxs % get_xs('fission', p_g, UVW=p_uvw) / &
                   matxs % get_xs('fission', p_g, UVW=p_uvw)
            end if
          end if

        else
          if (i_nuclide > 0) then
            score = nucxs % get_xs('nu-fission', p_g, UVW=p_uvw) * &
                 atom_density * flux
          else
            score = matxs % get_xs('nu-fission', p_g, UVW=p_uvw) * flux
          end if
        end if


      case (SCORE_PROMPT_NU_FISSION)

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing .or. p % fission) then
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              ! Normally, we only need to make contributions to one scoring
              ! bin. However, in the case of fission, since multiple fission
              ! neutrons were emitted with different energies, multiple
              ! outgoing energy bins may have been scored to. The following
              ! logic treats this special case and results to multiple bins
              call score_fission_eout(p, t, score_index, score_bin)
              cycle SCORE_LOOP
            end if
          end if
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! nu-fission
            score = p % absorb_wgt * flux
            if (i_nuclide > 0) then
              score = score * atom_density * &
                   nucxs % get_xs('prompt-nu-fission', p_g, UVW=p_uvw) / &
                   matxs % get_xs('absorption', p_g, UVW=p_uvw)
            else
              score = score * &
                   matxs % get_xs('prompt-nu-fission', p_g, UVW=p_uvw) / &
                   matxs % get_xs('absorption', p_g, UVW=p_uvw)
            end if
          else
            ! Skip any non-fission events
            if (.not. p % fission) cycle SCORE_LOOP
            ! If there is no outgoing energy filter, than we only need to
            ! score to one bin. For the score to be 'analog', we need to
            ! score the number of particles that were banked in the fission
            ! bank. Since this was weighted by 1/keff, we multiply by keff
            ! to get the proper score.
            score = keff * p % wgt_bank * (ONE - sum(p % n_delayed_bank) &
                 / real(p % n_bank, 8)) * flux
            if (i_nuclide > 0) then
              score = score * atom_density * &
                   nucxs % get_xs('fission', p_g, UVW=p_uvw) / &
                   matxs % get_xs('fission', p_g, UVW=p_uvw)
            end if
          end if

        else
          if (i_nuclide > 0) then
            score = nucxs % get_xs('prompt-nu-fission', p_g, UVW=p_uvw) * &
                 atom_density * flux
          else
            score = matxs % get_xs('prompt-nu-fission', p_g, UVW=p_uvw) * flux
          end if
        end if


      case (SCORE_DELAYED_NU_FISSION)

        ! Set the delayedgroup filter index and the number of delayed group bins
        dg_filter = t % find_filter(FILTER_DELAYEDGROUP)

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing .or. p % fission) then
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              ! Normally, we only need to make contributions to one scoring
              ! bin. However, in the case of fission, since multiple fission
              ! neutrons were emitted with different energies, multiple
              ! outgoing energy bins may have been scored to. The following
              ! logic treats this special case and results to multiple bins
              call score_fission_eout(p, t, score_index, score_bin)
              cycle SCORE_LOOP
            end if
          end if
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! nu-fission
            if (matxs % get_xs('absorption', p_g, UVW=p_uvw) > ZERO) then

              if (dg_filter > 0) then
                select type(filt => filters(t % filter(dg_filter)) % obj)
                type is (DelayedGroupFilter)

                  ! Loop over all delayed group bins and tally to them
                  ! individually
                  do d_bin = 1, filt % n_bins

                    ! Get the delayed group for this bin
                    d = filt % groups(d_bin)

                    score = p % absorb_wgt * flux
                    if (i_nuclide > 0) then
                      score = score * nucxs % get_xs('delayed-nu-fission', &
                           p_g, UVW=p_uvw, dg=d) / &
                           matxs % get_xs('absorption', p_g, UVW=p_uvw)
                    else
                      score = score * matxs % get_xs('delayed-nu-fission', &
                           p_g, UVW=p_uvw, dg=d) / &
                           matxs % get_xs('absorption', p_g, UVW=p_uvw)
                    end if

                    call score_fission_delayed_dg(t, d_bin, score, score_index)
                  end do
                  cycle SCORE_LOOP
                end select
              else
                score = p % absorb_wgt * flux
                if (i_nuclide > 0) then
                  score = score * nucxs % get_xs('delayed-nu-fission', p_g, &
                       UVW=p_uvw) / matxs % get_xs('absorption', p_g, UVW=p_uvw)
                else
                  score = score * matxs % get_xs('delayed-nu-fission', p_g, &
                       UVW=p_uvw) / matxs % get_xs('absorption', p_g, UVW=p_uvw)
                end if
              end if
            end if
          else
            ! Skip any non-fission events
            if (.not. p % fission) cycle SCORE_LOOP
            ! If there is no outgoing energy filter, than we only need to
            ! score to one bin. For the score to be 'analog', we need to
            ! score the number of particles that were banked in the fission
            ! bank. Since this was weighted by 1/keff, we multiply by keff
            ! to get the proper score.

            ! Check if the delayed group filter is present
            if (dg_filter > 0) then
              select type(filt => filters(t % filter(dg_filter)) % obj)
              type is (DelayedGroupFilter)

                ! Loop over all delayed group bins and tally to them
                ! individually
                do d_bin = 1, filt % n_bins

                  ! Get the delayed group for this bin
                  d = filt % groups(d_bin)

                  score = keff * p % wgt_bank / p % n_bank * &
                       p % n_delayed_bank(d) * flux

                  if (i_nuclide > 0) then
                    score = score * atom_density * &
                         nucxs % get_xs('fission', p_g, UVW=p_uvw) / &
                         matxs % get_xs('fission', p_g, UVW=p_uvw)
                  end if

                  call score_fission_delayed_dg(t, d_bin, score, score_index)
                end do
                cycle SCORE_LOOP
              end select
            else
              score = keff * p % wgt_bank / p % n_bank * sum(p % n_delayed_bank) * flux
              if (i_nuclide > 0) then
                score = score * atom_density * &
                     nucxs % get_xs('fission', p_g, UVW=p_uvw) / &
                     matxs % get_xs('fission', p_g, UVW=p_uvw)
              end if
            end if
          end if
        else

          ! Check if the delayed group filter is present
          if (dg_filter > 0) then
            select type(filt => filters(t % filter(dg_filter)) % obj)
            type is (DelayedGroupFilter)

              ! Loop over all delayed group bins and tally to them
              ! individually
              do d_bin = 1, filt % n_bins

                ! Get the delayed group for this bin
                d = filt % groups(d_bin)

                if (i_nuclide > 0) then
                  score = nucxs % get_xs('delayed-nu-fission', p_g, &
                       UVW=p_uvw, dg=d) * atom_density * flux
                else
                  score = matxs % get_xs('delayed-nu-fission', p_g, &
                       UVW=p_uvw, dg=d) * flux
                end if

                call score_fission_delayed_dg(t, d_bin, score, score_index)
              end do
              cycle SCORE_LOOP
            end select
          else
            if (i_nuclide > 0) then
              score = nucxs % get_xs('delayed-nu-fission', p_g, UVW=p_uvw) &
                   * atom_density * flux
            else
              score = matxs % get_xs('delayed-nu-fission', p_g, UVW=p_uvw) &
                   * flux
            end if
          end if
        end if

      case (SCORE_DECAY_RATE)

        ! Set the delayedgroup filter index and the number of delayed group bins
        dg_filter = t % find_filter(FILTER_DELAYEDGROUP)

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! nu-fission
            if (matxs % get_xs('absorption', p_g, UVW=p_uvw) > ZERO) then

              if (dg_filter > 0) then
                select type(filt => filters(t % filter(dg_filter)) % obj)
                type is (DelayedGroupFilter)

                  ! Loop over all delayed group bins and tally to them
                  ! individually
                  do d_bin = 1, filt % n_bins

                    ! Get the delayed group for this bin
                    d = filt % groups(d_bin)

                    score = p % absorb_wgt * flux
                    if (i_nuclide > 0) then
                      score = score * nucxs % get_xs('decay rate', p_g, &
                           UVW=p_uvw, dg=d) * &
                           nucxs % get_xs('delayed-nu-fission', p_g, &
                           UVW=p_uvw, dg=d) / matxs % get_xs('absorption', &
                           p_g, UVW=p_uvw)
                    else
                      score = score * matxs % get_xs('decay rate', p_g, &
                           UVW=p_uvw, dg=d) * &
                           matxs % get_xs('delayed-nu-fission', p_g, &
                           UVW=p_uvw, dg=d) / matxs % get_xs('absorption', &
                           p_g, UVW=p_uvw)
                    end if

                    call score_fission_delayed_dg(t, d_bin, score, score_index)
                  end do
                  cycle SCORE_LOOP
                end select
              else

                score = ZERO

                ! If the delayed group filter is not present, compute the score
                ! by accumulating the absorbed weight times the decay rate times
                ! the fraction of the delayed-nu-fission xs to the absorption xs
                ! for all delayed groups.
                do d = 1, num_delayed_groups
                  if (i_nuclide > 0) then
                    score = score + p % absorb_wgt * &
                         nucxs % get_xs('decay rate', p_g, UVW=p_uvw, dg=d) * &
                         nucxs % get_xs('delayed-nu-fission', p_g, UVW=p_uvw, &
                         dg=d) / matxs % get_xs('absorption', p_g, UVW=p_uvw) * flux
                  else
                    score = score + p % absorb_wgt * &
                         matxs % get_xs('decay rate', p_g, UVW=p_uvw, dg=d) * &
                         matxs % get_xs('delayed-nu-fission', p_g, UVW=p_uvw, &
                         dg=d) / matxs % get_xs('absorption', p_g, UVW=p_uvw) * flux
                  end if
                end do
              end if
            end if
          else
            ! Skip any non-fission events
            if (.not. p % fission) cycle SCORE_LOOP
            ! If there is no outgoing energy filter, than we only need to
            ! score to one bin. For the score to be 'analog', we need to
            ! score the number of particles that were banked in the fission
            ! bank. Since this was weighted by 1/keff, we multiply by keff
            ! to get the proper score.

            score = ZERO

            ! loop over number of particles banked
            do k = 1, p % n_bank

              ! get the delayed group
              g = fission_bank(n_bank - p % n_bank + k) % delayed_group

              ! Case for tallying delayed emissions
              if (g /= 0) then

                ! determine score based on bank site weight and keff.
                if (i_nuclide > 0) then
                  score = score + keff * atom_density * &
                       fission_bank(n_bank - p % n_bank + k) % wgt * &
                       nucxs % get_xs('decay rate', p_g, UVW=p_uvw, dg=g) * &
                       nucxs % get_xs('fission', p_g, UVW=p_uvw) / &
                       matxs % get_xs('fission', p_g, UVW=p_uvw) * flux
                else
                  score = score + keff * &
                       fission_bank(n_bank - p % n_bank + k) % wgt * &
                       matxs % get_xs('decay rate', p_g, UVW=p_uvw, dg=g) * flux
                end if

                ! if the delayed group filter is present, tally to corresponding
                ! delayed group bin if it exists
                if (dg_filter > 0) then

                  ! declare the delayed group filter type
                  select type(filt => filters(t % filter(dg_filter)) % obj)
                  type is (DelayedGroupFilter)

                    ! loop over delayed group bins until the corresponding bin
                    ! is found
                    do d_bin = 1, filt % n_bins
                      d = filt % groups(d_bin)

                      ! check whether the delayed group of the particle is equal
                      ! to the delayed group of this bin
                      if (d == g) then
                        call score_fission_delayed_dg(t, d_bin, score, &
                             score_index)
                      end if
                    end do
                  end select

                  ! Reset the score to zero
                  score = ZERO
                end if
              end if
            end do

            ! If the delayed group filter is present, cycle because the
            ! score_fission_delayed_dg(...) has already tallied the score
            if (dg_filter > 0) then
              cycle SCORE_LOOP
            end if
          end if
        else

          ! Check if the delayed group filter is present
          if (dg_filter > 0) then
            select type(filt => filters(t % filter(dg_filter)) % obj)
            type is (DelayedGroupFilter)

              ! Loop over all delayed group bins and tally to them
              ! individually
              do d_bin = 1, filt % n_bins

                ! Get the delayed group for this bin
                d = filt % groups(d_bin)

                if (i_nuclide > 0) then
                  score = nucxs % get_xs('decay rate', p_g, UVW=p_uvw, dg=d) * &
                       nucxs % get_xs('delayed-nu-fission', p_g, UVW=p_uvw, &
                       dg=d) * atom_density * flux
                else
                  score = matxs % get_xs('decay rate', p_g, UVW=p_uvw, dg=d) * &
                       matxs % get_xs('delayed-nu-fission', p_g, UVW=p_uvw, &
                       dg=d) * flux
                end if

                call score_fission_delayed_dg(t, d_bin, score, score_index)
              end do
              cycle SCORE_LOOP
            end select
          else
            score = ZERO

            ! If the delayed group filter is not present, compute the score
            ! by accumulating the absorbed weight times the decay rate times
            ! the fraction of the delayed-nu-fission xs to the absorption xs
            ! for all delayed groups.
            do d = 1, num_delayed_groups
              if (i_nuclide > 0) then
                score = score + atom_density * flux * &
                     nucxs % get_xs('decay rate', p_g, UVW=p_uvw, dg=d) * &
                     nucxs % get_xs('delayed-nu-fission', p_g, UVW=p_uvw, dg=d)
              else
                score = score + flux * &
                     matxs % get_xs('decay rate', p_g, UVW=p_uvw, dg=d) * &
                     matxs % get_xs('delayed-nu-fission', p_g, UVW=p_uvw, dg=d)
              end if
            end do
          end if
        end if


      case (SCORE_KAPPA_FISSION)

        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! fission
            score = p % absorb_wgt * flux
          else
            ! Skip any non-absorption events
            if (p % event == EVENT_SCATTER) cycle SCORE_LOOP
            ! All fission events will contribute, so again we can use
            ! particle's weight entering the collision as the estimate for the
            ! fission reaction rate
            score = p % last_wgt * flux
          end if
          if (i_nuclide > 0) then
            score = score * atom_density * &
                 nucxs % get_xs('kappa-fission', p_g, UVW=p_uvw) / &
                 matxs % get_xs('absorption', p_g, UVW=p_uvw)
          else
            score = score * &
                 matxs % get_xs('kappa-fission', p_g, UVW=p_uvw) / &
                 matxs % get_xs('absorption', p_g, UVW=p_uvw)
          end if
        else
          if (i_nuclide > 0) then
            score = nucxs % get_xs('kappa-fission', p_g, UVW=p_uvw) * &
                 atom_density * flux
          else
            score = matxs % get_xs('kappa-fission', p_g, UVW=p_uvw) * flux

          end if
        end if


      case (SCORE_EVENTS)
        ! Simply count number of scoring events
        score = ONE

      end select

      !#########################################################################
      ! Expand score if necessary and add to tally results.
      call expand_and_score(p, t, score_index, filter_index, score_bin, &
                            score, i)

    end do SCORE_LOOP

    nullify(matxs, nucxs)
  end subroutine score_general_mg

!===============================================================================
! EXPAND_AND_SCORE takes a previously determined score value and adjusts it
! if necessary (for functional expansion weighting), and then adds the resultant
! value to the tally results array.
!===============================================================================

  subroutine expand_and_score(p, t, score_index, filter_index, score_bin, &
                              score, i)
    type(Particle),    intent(in)    :: p
    type(TallyObject), intent(inout) :: t
    integer,           intent(inout) :: score_index
    integer,           intent(in)    :: filter_index ! for % results
    integer,           intent(in)    :: score_bin    ! score of concern
    real(8),           intent(inout) :: score        ! data to score
    integer,           intent(inout) :: i            ! Working index

    integer :: num_nm ! Number of N,M orders in harmonic
    integer :: n      ! Moment loop index
    real(8) :: uvw(3)

    select case(score_bin)
    case (SCORE_SCATTER_N, SCORE_NU_SCATTER_N)
      ! Find the scattering order for a singly requested moment, and
      ! store its moment contribution.
      if (t % moment_order(i) == 1) then
        score = score * p % mu ! avoid function call overhead
      else
        score = score * calc_pn(t % moment_order(i), p % mu)
      endif
!$omp atomic
      t % results(RESULT_VALUE, score_index, filter_index) = &
           t % results(RESULT_VALUE, score_index, filter_index) + score


    case(SCORE_SCATTER_YN, SCORE_NU_SCATTER_YN)
      score_index = score_index - 1
      num_nm = 1
      ! Find the order for a collection of requested moments
      ! and store the moment contribution of each
      do n = 0, t % moment_order(i)
        ! determine scoring bin index
        score_index = score_index + num_nm
        ! Update number of total n,m bins for this n (m = [-n: n])
        num_nm = 2 * n + 1

        ! multiply score by the angular flux moments and store
!$omp critical (score_general_scatt_yn)
        t % results(RESULT_VALUE, score_index: score_index + num_nm - 1, &
             filter_index) = t % results(RESULT_VALUE, &
             score_index: score_index + num_nm - 1, filter_index) &
             + score * calc_pn(n, p % mu) * calc_rn(n, p % last_uvw)
!$omp end critical (score_general_scatt_yn)
      end do
      i = i + (t % moment_order(i) + 1)**2 - 1


    case(SCORE_FLUX_YN, SCORE_TOTAL_YN)
      score_index = score_index - 1
      num_nm = 1
      if (t % estimator == ESTIMATOR_ANALOG .or. &
           t % estimator == ESTIMATOR_COLLISION) then
        uvw = p % last_uvw
      else if (t % estimator == ESTIMATOR_TRACKLENGTH) then
        uvw = p % coord(1) % uvw
      end if
      ! Find the order for a collection of requested moments
      ! and store the moment contribution of each
      do n = 0, t % moment_order(i)
        ! determine scoring bin index
        score_index = score_index + num_nm
        ! Update number of total n,m bins for this n (m = [-n: n])
        num_nm = 2 * n + 1

        ! multiply score by the angular flux moments and store
!$omp critical (score_general_flux_tot_yn)
        t % results(RESULT_VALUE, score_index: score_index + num_nm - 1, &
             filter_index) = t % results(RESULT_VALUE, &
             score_index: score_index + num_nm - 1, filter_index) &
             + score * calc_rn(n, uvw)
!$omp end critical (score_general_flux_tot_yn)
      end do
      i = i + (t % moment_order(i) + 1)**2 - 1


    case (SCORE_SCATTER_PN, SCORE_NU_SCATTER_PN)
      score_index = score_index - 1
      ! Find the scattering order for a collection of requested moments
      ! and store the moment contribution of each
      do n = 0, t % moment_order(i)
        ! determine scoring bin index
        score_index = score_index + 1

        ! get the score and tally it
!$omp atomic
        t % results(RESULT_VALUE, score_index, filter_index) = &
             t % results(RESULT_VALUE, score_index, filter_index) &
             + score * calc_pn(n, p % mu)
      end do
      i = i + t % moment_order(i)


    case default
!$omp atomic
      t % results(RESULT_VALUE, score_index, filter_index) = &
           t % results(RESULT_VALUE, score_index, filter_index) + score

    end select

  end subroutine expand_and_score

!===============================================================================
! SCORE_ALL_NUCLIDES tallies individual nuclide reaction rates specifically when
! the user requests <nuclides>all</nuclides>.
!===============================================================================

  subroutine score_all_nuclides(p, t, flux, filter_index)

    type(Particle), intent(in) :: p
    type(TallyObject), intent(inout) :: t
    real(8),        intent(in) :: flux
    integer,        intent(in) :: filter_index

    integer :: i             ! loop index for nuclides in material
    integer :: i_nuclide     ! index in nuclides array
    real(8) :: atom_density  ! atom density of single nuclide in atom/b-cm
    type(Material),    pointer :: mat

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

      ! Determine score for each bin
      call score_general(p, t, (i_nuclide-1)*t % n_score_bins, filter_index, &
           i_nuclide, atom_density, flux)

    end do NUCLIDE_LOOP

    ! ==========================================================================
    ! SCORE TOTAL MATERIAL REACTION RATES

    i_nuclide = -1
    atom_density = ZERO

    ! Determine score for each bin
    call score_general(p, t, n_nuclides*t % n_score_bins, filter_index, &
         i_nuclide, atom_density, flux)

  end subroutine score_all_nuclides

!===============================================================================
! SCORE_ANALOG_TALLY keeps track of how many events occur in a specified cell,
! energy range, etc. Note that since these are "analog" tallies, they are only
! triggered at every collision, not every event
!===============================================================================

  subroutine score_analog_tally_ce(p)

    type(Particle), intent(in) :: p

    integer :: i, j
    integer :: i_tally
    integer :: i_filt
    integer :: i_bin
    integer :: k                    ! loop index for nuclide bins
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array
    real(8) :: filter_weight        ! combined weight of all filters
    logical :: finished             ! found all valid bin combinations

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do i = 1, active_analog_tallies % size()
      ! Get index of tally and pointer to tally
      i_tally = active_analog_tallies % data(i)
      associate (t => tallies(i_tally) % obj)

      ! Find all valid bins in each filter if they have not already been found
      ! for a previous tally.
      do j = 1, size(t % filter)
        i_filt = t % filter(j)
        if (.not. filter_matches(i_filt) % bins_present) then
          call filter_matches(i_filt) % bins % clear()
          call filter_matches(i_filt) % weights % clear()
          call filters(i_filt) % obj % get_all_bins(p, t % estimator, &
               filter_matches(i_filt))
          filter_matches(i_filt) % bins_present = .true.
        end if
        ! If there are no valid bins for this filter, then there is nothing to
        ! score and we can move on to the next tally.
        if (filter_matches(i_filt) % bins % size() == 0) cycle TALLY_LOOP

        ! Set the index of the bin used in the first filter combination
        filter_matches(i_filt) % i_bin = 1
      end do

      ! ========================================================================
      ! Loop until we've covered all valid bins on each of the filters.

      FILTER_LOOP: do

        ! Reset scoring index and weight
        filter_index = 1
        filter_weight = ONE

        ! Determine scoring index and weight for this filter combination
        do j = 1, size(t % filter)
          i_filt = t % filter(j)
          i_bin = filter_matches(i_filt) % i_bin
          filter_index = filter_index + (filter_matches(i_filt) % bins % &
               data(i_bin) - 1) * t % stride(j)
          filter_weight = filter_weight * filter_matches(i_filt) % weights % &
               data(i_bin)
        end do

        ! ======================================================================
        ! Nuclide logic

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
              k = n_nuclides + 1
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
          call score_general(p, t, (k-1)*t % n_score_bins, filter_index, &
               i_nuclide, ZERO, filter_weight)

        end do NUCLIDE_LOOP

        ! ======================================================================
        ! Filter logic

        ! Increment the filter bins, starting with the last filter to find the
        ! next valid bin combination
        finished = .true.
        do j = size(t % filter), 1, -1
          i_filt = t % filter(j)
          if (filter_matches(i_filt) % i_bin < filter_matches(i_filt) % &
               bins % size()) then
            filter_matches(i_filt) % i_bin = filter_matches(i_filt) % i_bin + 1
            finished = .false.
            exit
          else
            filter_matches(i_filt) % i_bin = 1
          end if
        end do

        ! Once we have finished all valid bins for each of the filters, exit
        ! the loop.
        if (finished) exit FILTER_LOOP

      end do FILTER_LOOP

      ! If the user has specified that we can assume all tallies are spatially
      ! separate, this implies that once a tally has been scored to, we needn't
      ! check the others. This cuts down on overhead when there are many
      ! tallies specified

      if (assume_separate) exit TALLY_LOOP

      end associate
    end do TALLY_LOOP

    ! Reset filter matches flag
    filter_matches(:) % bins_present = .false.

  end subroutine score_analog_tally_ce

  subroutine score_analog_tally_mg(p)

    type(Particle), intent(in) :: p

    integer :: i, j
    integer :: i_tally
    integer :: i_filt
    integer :: i_bin
    integer :: k                    ! loop index for nuclide bins
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array
    real(8) :: filter_weight        ! combined weight of all filters
    real(8) :: atom_density
    logical :: finished             ! found all valid bin combinations
    type(Material),    pointer :: mat

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do i = 1, active_analog_tallies % size()
      ! Get index of tally and pointer to tally
      i_tally = active_analog_tallies % data(i)
      associate (t => tallies(i_tally) % obj)

      ! Find all valid bins in each filter if they have not already been found
      ! for a previous tally.
      do j = 1, size(t % filter)
        i_filt = t % filter(j)
        if (.not. filter_matches(i_filt) % bins_present) then
          call filter_matches(i_filt) % bins % clear()
          call filter_matches(i_filt) % weights % clear()
          call filters(i_filt) % obj % get_all_bins(p, t % estimator, &
               filter_matches(i_filt))
          filter_matches(i_filt) % bins_present = .true.
        end if
        ! If there are no valid bins for this filter, then there is nothing to
        ! score and we can move on to the next tally.
        if (filter_matches(i_filt) % bins % size() == 0) cycle TALLY_LOOP

        ! Set the index of the bin used in the first filter combination
        filter_matches(i_filt) % i_bin = 1
      end do

      ! ========================================================================
      ! Loop until we've covered all valid bins on each of the filters.

      FILTER_LOOP: do

        ! Reset scoring index and weight
        filter_index = 1
        filter_weight = ONE

        ! Determine scoring index and weight for this filter combination
        do j = 1, size(t % filter)
          i_filt = t % filter(j)
          i_bin = filter_matches(i_filt) % i_bin
          filter_index = filter_index + (filter_matches(i_filt) % bins % &
               data(i_bin) - 1) * t % stride(j)
          filter_weight = filter_weight * filter_matches(i_filt) % weights % &
               data(i_bin)
        end do

        ! ======================================================================
        ! Nuclide logic

        ! Check for nuclide bins
        NUCLIDE_LOOP: do k = 1, t % n_nuclide_bins
          ! Get index of nuclide in nuclides array
          i_nuclide = t % nuclide_bins(k)

          if (i_nuclide > 0) then
            if (p % material /= MATERIAL_VOID) then
              ! Get pointer to current material
              mat => materials(p % material)

              ! Determine index of nuclide in Material % atom_density array
              j = mat % mat_nuclide_index(i_nuclide)
              if (j == 0) cycle NUCLIDE_LOOP

              ! Copy corresponding atom density
              atom_density = mat % atom_density(j)
            end if
          else
            atom_density = ZERO
          end if

          call score_general(p, t, (k-1)*t % n_score_bins, filter_index, &
               i_nuclide, atom_density, filter_weight)
        end do NUCLIDE_LOOP

        ! ======================================================================
        ! Filter logic

        ! Increment the filter bins, starting with the last filter to find the
        ! next valid bin combination
        finished = .true.
        do j = size(t % filter), 1, -1
          i_filt = t % filter(j)
          if (filter_matches(i_filt) % i_bin < filter_matches(i_filt) % &
               bins % size()) then
            filter_matches(i_filt) % i_bin = filter_matches(i_filt) % i_bin + 1
            finished = .false.
            exit
          else
            filter_matches(i_filt) % i_bin = 1
          end if
        end do

        ! Once we have finished all valid bins for each of the filters, exit
        ! the loop.
        if (finished) exit FILTER_LOOP

      end do FILTER_LOOP

      ! If the user has specified that we can assume all tallies are spatially
      ! separate, this implies that once a tally has been scored to, we needn't
      ! check the others. This cuts down on overhead when there are many
      ! tallies specified

      if (assume_separate) exit TALLY_LOOP

      end associate
    end do TALLY_LOOP

    ! Reset filter matches flag
    filter_matches(:) % bins_present = .false.

  end subroutine score_analog_tally_mg

!===============================================================================
! SCORE_FISSION_EOUT handles a special case where we need to store neutron
! production rate with an outgoing energy filter (think of a fission matrix). In
! this case, we may need to score to multiple bins if there were multiple
! neutrons produced with different energies.
!===============================================================================

  subroutine score_fission_eout(p, t, i_score, score_bin)

    type(Particle), intent(in)       :: p
    type(TallyObject), intent(inout) :: t
    integer, intent(in)              :: i_score ! index for score
    integer, intent(in)              :: score_bin

    integer :: i             ! index of outgoing energy filter
    integer :: j             ! index of delayedgroup filter
    integer :: d             ! delayed group
    integer :: g             ! another delayed group
    integer :: d_bin         ! delayed group bin index
    integer :: n             ! number of energies on filter
    integer :: k             ! loop index for bank sites
    integer :: l             ! loop index for tally filters
    integer :: f             ! index in filters array
    integer :: b             ! index of filter bin
    integer :: i_bin         ! index of matching filter bin
    integer :: bin_energyout ! original outgoing energy bin
    integer :: i_filter      ! index for matching filter bin combination
    real(8) :: filter_weight ! combined weight of all filters
    real(8) :: score         ! actual score
    real(8) :: E_out         ! energy of fission bank site
    integer :: g_out         ! energy group of fission bank site

    ! save original outgoing energy bin and score index
    i = t % filter(t % find_filter(FILTER_ENERGYOUT))
    i_bin = filter_matches(i) % i_bin
    bin_energyout = filter_matches(i) % bins % data(i_bin)

    ! declare the energyout filter type
    select type(eo_filt => filters(i) % obj)
    type is (EnergyoutFilter)

      ! Get number of energies on filter
      n = size(eo_filt % bins)

      ! Since the creation of fission sites is weighted such that it is
      ! expected to create n_particles sites, we need to multiply the
      ! score by keff to get the true nu-fission rate. Otherwise, the sum
      ! of all nu-fission rates would be ~1.0.

      ! loop over number of particles banked
      do k = 1, p % n_bank

        ! get the delayed group
        g = fission_bank(n_bank - p % n_bank + k) % delayed_group

        ! determine score based on bank site weight and keff
        score = keff * fission_bank(n_bank - p % n_bank + k) % wgt

        ! Add derivative information for differential tallies.  Note that the
        ! i_nuclide and atom_density arguments do not matter since this is an
        ! analog estimator.
        if (t % deriv /= NONE) then
          call apply_derivative_to_score(p, t, 0, ZERO, SCORE_NU_FISSION, score)
        end if

        if (.not. run_CE .and. eo_filt % matches_transport_groups) then

          ! determine outgoing energy group from fission bank
          g_out = int(fission_bank(n_bank - p % n_bank + k) % E)

          ! modify the value so that g_out = 1 corresponds to the highest
          ! energy bin
          g_out = size(eo_filt % bins) - g_out

          ! change outgoing energy bin
          filter_matches(i) % bins % data(i_bin) = g_out

        else

          ! determine outgoing energy from fission bank
          if (run_CE) then
            E_out = fission_bank(n_bank - p % n_bank + k) % E
          else
            E_out = energy_bin_avg(int(fission_bank(n_bank - p % n_bank + k) &
                 % E))
          end if

          ! check if outgoing energy is within specified range on filter
          if (E_out < eo_filt % bins(1) .or. E_out > eo_filt % bins(n)) cycle

          ! change outgoing energy bin
          filter_matches(i) % bins % data(i_bin) = &
               binary_search(eo_filt % bins, n, E_out)

        end if

        ! Case for tallying prompt neutrons
        if (score_bin == SCORE_NU_FISSION .or. &
             (score_bin == SCORE_PROMPT_NU_FISSION .and. g == 0)) then

          ! determine scoring index and weight for this filter combination
          i_filter = 1
          do l = 1, size(t % filter)
            i_filter = i_filter + (filter_matches(t % filter(l)) % bins % &
                 data(filter_matches(t % filter(l)) % i_bin) - 1) * &
                 t % stride(l)
          end do

          ! Add score to tally
!$omp atomic
          t % results(RESULT_VALUE, i_score, i_filter) = &
               t % results(RESULT_VALUE, i_score, i_filter) + score

        ! Case for tallying delayed emissions
        else if (score_bin == SCORE_DELAYED_NU_FISSION .and. g /= 0) then

          ! Get the index of delayed group filter
          j = t % find_filter(FILTER_DELAYEDGROUP)

          ! if the delayed group filter is present, tally to corresponding
          ! delayed group bin if it exists
          if (j > 0) then

            ! declare the delayed group filter type
            select type(dg_filt => filters(t % filter(j)) % obj)
            type is (DelayedGroupFilter)

              ! loop over delayed group bins until the corresponding bin is
              ! found
              do d_bin = 1, dg_filt % n_bins
                d = dg_filt % groups(d_bin)

                ! check whether the delayed group of the particle is equal to
                ! the delayed group of this bin
                if (d == g) then

                  ! Reset scoring index and filter weight
                  i_filter = 1
                  filter_weight = ONE

                  ! determine scoring index and weight for this filter
                  ! combination
                  do l = 1, size(t % filter)
                    f = t % filter(l)
                    b = filter_matches(f) % i_bin
                    i_filter = i_filter + (filter_matches(f) % bins % &
                         data(b) - 1) * t % stride(l)
                    filter_weight = filter_weight * filter_matches(f) % &
                         weights % data(b)
                  end do

                  call score_fission_delayed_dg(t, d_bin, &
                       score * filter_weight, i_score)
                end if
              end do
            end select

          ! if the delayed group filter is not present, add score to tally
          else

            ! Reset scoring index and filter weight
            i_filter = 1
            filter_weight = ONE

            ! determine scoring index and weight for this filter combination
            do l = 1, size(t % filter)
              f = t % filter(l)
              b = filter_matches(f) % i_bin
              i_filter = i_filter + (filter_matches(f) % bins % data(b) - 1) &
                   * t % stride(l)
              filter_weight = filter_weight * filter_matches(f) % weights % &
                   data(b)
            end do

            ! Add score to tally
!$omp atomic
            t % results(RESULT_VALUE, i_score, i_filter) = &
                 t % results(RESULT_VALUE, i_score, i_filter) + score * filter_weight
          end if
        end if
      end do
    end select

    ! reset outgoing energy bin and score index
    filter_matches(i) % bins % data(i_bin) = bin_energyout

  end subroutine score_fission_eout

!===============================================================================
! SCORE_FISSION_DELAYED_DG helper function used to increment the tally when a
! delayed group filter is present.
!===============================================================================

  subroutine score_fission_delayed_dg(t, d_bin, score, score_index)

    type(TallyObject), intent(inout) :: t
    integer, intent(in)              :: d_bin       ! delayed group bin index
    real(8), intent(in)              :: score       ! actual score
    integer, intent(in)              :: score_index ! index for score

    integer :: i             ! loop over tally filters
    integer :: i_filt        ! index in filters array
    integer :: i_bin         ! index of matching filter bin
    integer :: bin_original  ! original bin index
    integer :: filter_index  ! index for matching filter bin combination

    ! save original delayed group bin
    i_filt = t % filter(t % find_filter(FILTER_DELAYEDGROUP))
    i_bin = filter_matches(i_filt) % i_bin
    bin_original = filter_matches(i_filt) % bins % data(i_bin)
    filter_matches(i_filt) % bins % data(i_bin) = d_bin

    ! determine scoring index and weight on the modified matching bins
    filter_index = 1
    do i = 1, size(t % filter)
      filter_index = filter_index + (filter_matches(t % filter(i)) % bins % &
           data(filter_matches(t % filter(i)) % i_bin) - 1) * t % stride(i)
    end do

!$omp atomic
    t % results(RESULT_VALUE, score_index, filter_index) = &
         t % results(RESULT_VALUE, score_index, filter_index) + score

    ! reset original delayed group bin
    filter_matches(i_filt) % bins % data(i_bin) = bin_original

  end subroutine score_fission_delayed_dg

!===============================================================================
! SCORE_TRACKLENGTH_TALLY calculates fluxes and reaction rates based on the
! track-length estimate of the flux. This is triggered at every event (surface
! crossing, lattice crossing, or collision) and thus cannot be done for tallies
! that require post-collision information.
!===============================================================================

  subroutine score_tracklength_tally(p, distance)

    type(Particle), intent(in) :: p
    real(8),        intent(in) :: distance

    integer :: i
    integer :: i_tally
    integer :: i_filt
    integer :: i_bin
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array (from bins)
    real(8) :: flux                 ! tracklength estimate of flux
    real(8) :: atom_density         ! atom density of single nuclide in atom/b-cm
    real(8) :: filter_weight        ! combined weight of all filters
    logical :: finished             ! found all valid bin combinations
    type(Material),    pointer :: mat

    ! Determine track-length estimate of flux
    flux = p % wgt * distance

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do i = 1, active_tracklength_tallies % size()
      ! Get index of tally and pointer to tally
      i_tally = active_tracklength_tallies % data(i)
      associate (t => tallies(i_tally) % obj)

      ! Find all valid bins in each filter if they have not already been found
      ! for a previous tally.
      do j = 1, size(t % filter)
        i_filt = t % filter(j)
        if (.not. filter_matches(i_filt) % bins_present) then
          call filter_matches(i_filt) % bins % clear()
          call filter_matches(i_filt) % weights % clear()
          call filters(i_filt) % obj % get_all_bins(p, t % estimator, &
               filter_matches(i_filt))
          filter_matches(i_filt) % bins_present = .true.
        end if
        ! If there are no valid bins for this filter, then there is nothing to
        ! score and we can move on to the next tally.
        if (filter_matches(i_filt) % bins % size() == 0) cycle TALLY_LOOP

        ! Set the index of the bin used in the first filter combination
        filter_matches(i_filt) % i_bin = 1
      end do

      ! ========================================================================
      ! Loop until we've covered all valid bins on each of the filters.

      FILTER_LOOP: do

        ! Reset scoring index and weight
        filter_index = 1
        filter_weight = ONE

        ! Determine scoring index and weight for this filter combination
        do j = 1, size(t % filter)
          i_filt = t % filter(j)
          i_bin = filter_matches(i_filt) % i_bin
          filter_index = filter_index + (filter_matches(i_filt) % bins % &
               data(i_bin) - 1) * t % stride(j)
          filter_weight = filter_weight * filter_matches(i_filt) % weights % &
               data(i_bin)
        end do

        ! ======================================================================
        ! Nuclide logic

        if (t % all_nuclides) then
          if (p % material /= MATERIAL_VOID) then
            call score_all_nuclides(p, t, flux * filter_weight, filter_index)
          end if
        else

          NUCLIDE_BIN_LOOP: do k = 1, t % n_nuclide_bins
            ! Get index of nuclide in nuclides array
            i_nuclide = t % nuclide_bins(k)

            if (i_nuclide > 0) then
              if (p % material /= MATERIAL_VOID) then
                ! Get pointer to current material
                mat => materials(p % material)

                ! Determine index of nuclide in Material % atom_density array
                j = mat % mat_nuclide_index(i_nuclide)
                if (j == 0) cycle NUCLIDE_BIN_LOOP

                ! Copy corresponding atom density
                atom_density = mat % atom_density(j)
              else
                atom_density = ZERO
              end if
            end if

            ! Determine score for each bin
            call score_general(p, t, (k-1)*t % n_score_bins, filter_index, &
                 i_nuclide, atom_density, flux * filter_weight)

          end do NUCLIDE_BIN_LOOP

        end if

        ! ======================================================================
        ! Filter logic

        ! Increment the filter bins, starting with the last filter to find the
        ! next valid bin combination
        finished = .true.
        do j = size(t % filter), 1, -1
          i_filt = t % filter(j)
          if (filter_matches(i_filt) % i_bin < filter_matches(i_filt) % &
               bins % size()) then
            filter_matches(i_filt) % i_bin = filter_matches(i_filt) % i_bin + 1
            finished = .false.
            exit
          else
            filter_matches(i_filt) % i_bin = 1
          end if
        end do

        ! Once we have finished all valid bins for each of the filters, exit
        ! the loop.
        if (finished) exit FILTER_LOOP

      end do FILTER_LOOP

      ! If the user has specified that we can assume all tallies are spatially
      ! separate, this implies that once a tally has been scored to, we needn't
      ! check the others. This cuts down on overhead when there are many
      ! tallies specified

      if (assume_separate) exit TALLY_LOOP

      end associate
    end do TALLY_LOOP

    ! Reset filter matches flag
    filter_matches(:) % bins_present = .false.

  end subroutine score_tracklength_tally

!===============================================================================
! SCORE_COLLISION_TALLY calculates fluxes and reaction rates based on the
! 1/Sigma_t estimate of the flux.  This is triggered after every collision.  It
! is invalid for tallies that require post-collison information because it can
! score reactions that didn't actually occur, and we don't a priori know what
! the outcome will be for reactions that we didn't sample.
!===============================================================================

  subroutine score_collision_tally(p)

    type(Particle), intent(in) :: p

    integer :: i
    integer :: i_tally
    integer :: i_filt
    integer :: i_bin
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array (from bins)
    real(8) :: flux                 ! collision estimate of flux
    real(8) :: atom_density         ! atom density of single nuclide
                                    !   in atom/b-cm
    real(8) :: filter_weight        ! combined weight of all filters
    logical :: finished             ! found all valid bin combinations
    type(Material),    pointer :: mat

    ! Determine collision estimate of flux
    if (survival_biasing) then
      ! We need to account for the fact that some weight was already absorbed
      flux = (p % last_wgt + p % absorb_wgt) / material_xs % total
    else
      flux = p % last_wgt / material_xs % total
    end if

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do i = 1, active_collision_tallies % size()
      ! Get index of tally and pointer to tally
      i_tally = active_collision_tallies % data(i)
      associate (t => tallies(i_tally) % obj)

      ! Find all valid bins in each filter if they have not already been found
      ! for a previous tally.
      do j = 1, size(t % filter)
        i_filt = t % filter(j)
        if (.not. filter_matches(i_filt) % bins_present) then
          call filter_matches(i_filt) % bins % clear()
          call filter_matches(i_filt) % weights % clear()
          call filters(i_filt) % obj % get_all_bins(p, t % estimator, &
               filter_matches(i_filt))
          filter_matches(i_filt) % bins_present = .true.
        end if
        ! If there are no valid bins for this filter, then there is nothing to
        ! score and we can move on to the next tally.
        if (filter_matches(i_filt) % bins % size() == 0) cycle TALLY_LOOP

        ! Set the index of the bin used in the first filter combination
        filter_matches(i_filt) % i_bin = 1
      end do

      ! ========================================================================
      ! Loop until we've covered all valid bins on each of the filters.

      FILTER_LOOP: do

        ! Reset scoring index and weight
        filter_index = 1
        filter_weight = ONE

        ! Determine scoring index and weight for this filter combination
        do j = 1, size(t % filter)
          i_filt = t % filter(j)
          i_bin = filter_matches(i_filt) % i_bin
          filter_index = filter_index + (filter_matches(i_filt) % bins % &
               data(i_bin) - 1) * t % stride(j)
          filter_weight = filter_weight * filter_matches(i_filt) % weights % &
               data(i_bin)
        end do

        ! ======================================================================
        ! Nuclide logic

        if (t % all_nuclides) then
          if (p % material /= MATERIAL_VOID) then
            call score_all_nuclides(p, t, flux * filter_weight, filter_index)
          end if
        else

          NUCLIDE_BIN_LOOP: do k = 1, t % n_nuclide_bins
            ! Get index of nuclide in nuclides array
            i_nuclide = t % nuclide_bins(k)

            if (i_nuclide > 0) then
              if (p % material /= MATERIAL_VOID) then
                ! Get pointer to current material
                mat => materials(p % material)

                ! Determine index of nuclide in Material % atom_density array
                j = mat % mat_nuclide_index(i_nuclide)
                if (j == 0) cycle NUCLIDE_BIN_LOOP

                ! Copy corresponding atom density
                atom_density = mat % atom_density(j)
              else
                atom_density = ZERO
              end if
            end if

            ! Determine score for each bin
            call score_general(p, t, (k-1)*t % n_score_bins, filter_index, &
                 i_nuclide, atom_density, flux * filter_weight)

          end do NUCLIDE_BIN_LOOP

        end if

        ! ======================================================================
        ! Filter logic

        ! Increment the filter bins, starting with the last filter to find the
        ! next valid bin combination
        finished = .true.
        do j = size(t % filter), 1, -1
          i_filt = t % filter(j)
          if (filter_matches(i_filt) % i_bin < filter_matches(i_filt) % &
               bins % size()) then
            filter_matches(i_filt) % i_bin = filter_matches(i_filt) % i_bin + 1
            finished = .false.
            exit
          else
            filter_matches(i_filt) % i_bin = 1
          end if
        end do

        ! Once we have finished all valid bins for each of the filters, exit
        ! the loop.
        if (finished) exit FILTER_LOOP

      end do FILTER_LOOP

      ! If the user has specified that we can assume all tallies are spatially
      ! separate, this implies that once a tally has been scored to, we needn't
      ! check the others. This cuts down on overhead when there are many
      ! tallies specified

      if (assume_separate) exit TALLY_LOOP

      end associate
    end do TALLY_LOOP

    ! Reset filter matches flag
    filter_matches(:) % bins_present = .false.

  end subroutine score_collision_tally

!===============================================================================
! score_surface_tally is called at every surface crossing and can be used to
! tally total or partial currents between two cells
!===============================================================================

    subroutine  score_surface_tally(p)

    type(Particle), intent(in)    :: p

    integer :: i
    integer :: i_tally
    integer :: i_filt
    integer :: i_bin
    integer :: q                    ! loop index for scoring bins
    integer :: k                    ! working index for expand and score
    integer :: score_bin            ! scoring bin, e.g. SCORE_FLUX
    integer :: score_index          ! scoring bin index
    integer :: j                    ! loop index for scoring bins
    integer :: filter_index         ! single index for single bin
    real(8) :: flux                 ! collision estimate of flux
    real(8) :: filter_weight        ! combined weight of all filters
    real(8) :: score                ! analog tally score
    logical :: finished             ! found all valid bin combinations

    ! No collision, so no weight change when survival biasing
    flux = p % wgt

    TALLY_LOOP: do i = 1, active_surface_tallies % size()
      ! Get index of tally and pointer to tally
      i_tally = active_surface_tallies % data(i)
      associate (t => tallies(i_tally) % obj)

      ! Find all valid bins in each filter if they have not already been found
      ! for a previous tally.
      do j = 1, size(t % filter)
        i_filt = t % filter(j)
        if (.not. filter_matches(i_filt) % bins_present) then
          call filter_matches(i_filt) % bins % clear()
          call filter_matches(i_filt) % weights % clear()
          call filters(i_filt) % obj % get_all_bins(p, t % estimator, &
               filter_matches(i_filt))
          filter_matches(i_filt) % bins_present = .true.
        end if
        ! If there are no valid bins for this filter, then there is nothing to
        ! score and we can move on to the next tally.
        if (filter_matches(i_filt) % bins % size() == 0) cycle TALLY_LOOP

        ! Set the index of the bin used in the first filter combination
        filter_matches(i_filt) % i_bin = 1
      end do

      ! ========================================================================
      ! Loop until we've covered all valid bins on each of the filters.

      FILTER_LOOP: do

        ! Reset scoring index and weight
        filter_index = 1
        filter_weight = ONE

        ! Determine scoring index and weight for this filter combination
        do j = 1, size(t % filter)
          i_filt = t % filter(j)
          i_bin = filter_matches(i_filt) % i_bin
          filter_index = filter_index + (filter_matches(i_filt) % bins % &
               data(i_bin) - 1) * t % stride(j)
          filter_weight = filter_weight * filter_matches(i_filt) % weights % &
               data(i_bin)
        end do

        ! Determine score
        score = flux * filter_weight

        ! Currently only one score type
        k = 0
        SCORE_LOOP: do q = 1, t % n_user_score_bins
          k = k + 1

          ! determine what type of score bin
          score_bin = t % score_bins(q)

          ! determine scoring bin index, no offset from nuclide bins
          score_index = q

          ! Expand score if necessary and add to tally results.
          call expand_and_score(p, t, score_index, filter_index, score_bin, &
               score, k)
        end do SCORE_LOOP

        ! ======================================================================
        ! Filter logic

        ! Increment the filter bins, starting with the last filter to find the
        ! next valid bin combination
        finished = .true.
        do j = size(t % filter), 1, -1
          i_filt = t % filter(j)
          if (filter_matches(i_filt) % i_bin < filter_matches(i_filt) % &
               bins % size()) then
            filter_matches(i_filt) % i_bin = filter_matches(i_filt) % i_bin + 1
            finished = .false.
            exit
          else
            filter_matches(i_filt) % i_bin = 1
          end if
        end do

        ! Once we have finished all valid bins for each of the filters, exit
        ! the loop.
        if (finished) exit FILTER_LOOP

      end do FILTER_LOOP

      end associate

      ! If the user has specified that we can assume all tallies are spatially
      ! separate, this implies that once a tally has been scored to, we needn't
      ! check the others. This cuts down on overhead when there are many
      ! tallies specified

      if (assume_separate) exit TALLY_LOOP

    end do TALLY_LOOP

    ! Reset filter matches flag
    filter_matches(:) % bins_present = .false.

  end subroutine score_surface_tally

!===============================================================================
! SCORE_SURFACE_CURRENT tallies surface crossings in a mesh tally by manually
! determining which mesh surfaces were crossed
!===============================================================================

  subroutine score_surface_current(p)

    type(Particle), intent(in) :: p

    integer :: i
    integer :: i_tally
    integer :: j, k                 ! loop indices
    integer :: n_dim                ! num dimensions of the mesh
    integer :: d1                   ! dimension index
    integer :: d2                   ! dimension index
    integer :: d3                   ! dimension index
    integer :: ijk0(3)              ! indices of starting coordinates
    integer :: ijk1(3)              ! indices of ending coordinates
    integer :: n_cross              ! number of surface crossings
    integer :: filter_index         ! index of scoring bin
    integer :: i_filter_mesh        ! index of mesh filter in filters array
    integer :: i_filter_surf        ! index of surface filter in filters
    integer :: i_filter_energy      ! index of energy filter in filters
    real(8) :: uvw(3)               ! cosine of angle of particle
    real(8) :: xyz0(3)              ! starting/intermediate coordinates
    real(8) :: xyz1(3)              ! ending coordinates of particle
    real(8) :: xyz_cross(3)         ! coordinates of bounding surfaces
    real(8) :: d(3)                 ! distance to each bounding surface
    real(8) :: distance             ! actual distance traveled
    integer :: matching_bin         ! next valid filter bin
    logical :: start_in_mesh        ! particle's starting xyz in mesh?
    logical :: end_in_mesh          ! particle's ending xyz in mesh?
    logical :: cross_surface        ! whether the particle crosses a surface
    logical :: energy_filter        ! energy filter present
    type(RegularMesh), pointer :: m

    TALLY_LOOP: do i = 1, active_current_tallies % size()
      ! Copy starting and ending location of particle
      xyz0 = p % last_xyz_current
      xyz1 = p % coord(1) % xyz

      ! Get pointer to tally
      i_tally = active_current_tallies % data(i)
      associate (t => tallies(i_tally) % obj)

      ! Check for energy filter
      energy_filter = (t % find_filter(FILTER_ENERGYIN) > 0)

      ! Get index for mesh, surface, and energy filters
      i_filter_mesh = t % filter(t % find_filter(FILTER_MESH))
      i_filter_surf = t % filter(t % find_filter(FILTER_SURFACE))
      if (energy_filter) then
        i_filter_energy = t % filter(t % find_filter(FILTER_ENERGYIN))
      end if

      ! Reset the matching bins arrays
      call filter_matches(i_filter_mesh) % bins % resize(1)
      call filter_matches(i_filter_surf) % bins % resize(1)
      if (energy_filter) then
        call filter_matches(i_filter_energy) % bins % resize(1)
      end if

      ! Get pointer to mesh
      select type(filt => filters(i_filter_mesh) % obj)
      type is (MeshFilter)
        m => meshes(filt % mesh)
      end select

      n_dim = m % n_dimension

      ! Determine indices for starting and ending location
      call m % get_indices(xyz0, ijk0, start_in_mesh)
      call m % get_indices(xyz1, ijk1, end_in_mesh)

      ! Check to see if start or end is in mesh -- if not, check if track still
      ! intersects with mesh
      if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) then
        if (.not. m % intersects(xyz0, xyz1)) cycle
      end if

      ! Calculate number of surface crossings
      n_cross = sum(abs(ijk1(:n_dim) - ijk0(:n_dim)))
      if (n_cross == 0) then
        cycle
      end if

      ! Copy particle's direction
      uvw = p % coord(1) % uvw

      ! Determine incoming energy bin.  We need to tell the energy filter this
      ! is a tracklength tally so it uses the pre-collision energy.
      if (energy_filter) then
        call filter_matches(i_filter_energy) % bins % clear()
        call filter_matches(i_filter_energy) % weights % clear()
        call filters(i_filter_energy) % obj % get_all_bins(p, &
             ESTIMATOR_TRACKLENGTH, filter_matches(i_filter_energy))
        if (filter_matches(i_filter_energy) % bins % size() == 0) cycle
        matching_bin = filter_matches(i_filter_energy) % bins % data(1)
        filter_matches(i_filter_energy) % bins % data(1) = matching_bin
      end if

      ! Bounding coordinates
      do d1 = 1, n_dim
        if (uvw(d1) > 0) then
          xyz_cross(d1) = m % lower_left(d1) + ijk0(d1) * m % width(d1)
        else
          xyz_cross(d1) = m % lower_left(d1) + (ijk0(d1) - 1) * m % width(d1)
        end if
      end do

      do j = 1, n_cross
        ! Reset scoring bin index
        filter_matches(i_filter_surf) % bins % data(1) = 0

        ! Set the distances to infinity
        d = INFINITY

        ! Calculate distance to each bounding surface. We need to treat
        ! special case where the cosine of the angle is zero since this would
        ! result in a divide-by-zero.
        do d1 = 1, n_dim
          if (uvw(d1) == 0) then
            d(d1) = INFINITY
          else
            d(d1) = (xyz_cross(d1) - xyz0(d1))/uvw(d1)
          end if
        end do

        ! Determine the closest bounding surface of the mesh cell by
        ! calculating the minimum distance. Then use the minimum distance and
        ! direction of the particle to determine which surface was crossed.
        distance = minval(d)

        ! Loop over the dimensions
        do d1 = 1, n_dim

          ! Get the other dimensions.
          if (d1 == 1) then
            d2 = mod(d1, 3) + 1
            d3 = mod(d1 + 1, 3) + 1
          else
            d2 = mod(d1 + 1, 3) + 1
            d3 = mod(d1, 3) + 1
          end if

          ! Check whether distance is the shortest distance
          if (distance == d(d1)) then

            ! Check whether particle is moving in positive d1 direction
            if (uvw(d1) > 0) then

              ! Outward current on d1 max surface
              if (all(ijk0(:n_dim) >= 1) .and. &
                   all(ijk0(:n_dim) <= m % dimension)) then
                filter_matches(i_filter_surf) % bins % data(1) = d1 * 4 - 1
                filter_matches(i_filter_mesh) % bins % data(1) = &
                     m % get_bin_from_indices(ijk0)
                filter_index = 1
                do k = 1, size(t % filter)
                  filter_index = filter_index + (filter_matches(t % &
                       filter(k)) % bins % data(1) - 1) * t % stride(k)
                end do
!$omp atomic
                t % results(RESULT_VALUE, 1, filter_index) = &
                     t % results(RESULT_VALUE, 1, filter_index) + p % wgt
              end if

              ! Inward current on d1 min surface
              cross_surface = .false.
              select case(n_dim)

              case (1)
                if (ijk0(d1) >= 0 .and. ijk0(d1) <  m % dimension(d1)) then
                  cross_surface = .true.
                end if

              case (2)
                if (ijk0(d1) >= 0 .and. ijk0(d1) <  m % dimension(d1) .and. &
                     ijk0(d2) >= 1 .and. ijk0(d2) <= m % dimension(d2)) then
                  cross_surface = .true.
                end if

              case (3)
                if (ijk0(d1) >= 0 .and. ijk0(d1) <  m % dimension(d1) .and. &
                     ijk0(d2) >= 1 .and. ijk0(d2) <= m % dimension(d2) .and. &
                     ijk0(d3) >= 1 .and. ijk0(d3) <= m % dimension(d3)) then
                  cross_surface = .true.
                end if
              end select

              ! If the particle crossed the surface, tally the current
              if (cross_surface) then
                ijk0(d1) = ijk0(d1) + 1
                filter_matches(i_filter_surf) % bins % data(1) = d1 * 4 - 2
                filter_matches(i_filter_mesh) % bins % data(1) = &
                     m % get_bin_from_indices(ijk0)
                filter_index = 1
                do k = 1, size(t % filter)
                  filter_index = filter_index + (filter_matches(t % &
                       filter(k)) % bins % data(1) - 1) * t % stride(k)
                end do
!$omp atomic
                t % results(RESULT_VALUE, 1, filter_index) = &
                     t % results(RESULT_VALUE, 1, filter_index) + p % wgt
                ijk0(d1) = ijk0(d1) - 1
              end if

              ijk0(d1) = ijk0(d1) + 1
              xyz_cross(d1) = xyz_cross(d1) + m % width(d1)

              ! The particle is moving in the negative d1 direction
            else

              ! Outward current on d1 min surface
              if (all(ijk0(:n_dim) >= 1) .and. &
                   all(ijk0(:n_dim) <= m % dimension)) then
                filter_matches(i_filter_surf) % bins % data(1) = d1 * 4 - 3
                filter_matches(i_filter_mesh) % bins % data(1) = &
                     m % get_bin_from_indices(ijk0)
                filter_index = 1
                do k = 1, size(t % filter)
                  filter_index = filter_index + (filter_matches(t % &
                       filter(k)) % bins % data(1) - 1) * t % stride(k)
                end do
!$omp atomic
                t % results(RESULT_VALUE, 1, filter_index) = &
                     t % results(RESULT_VALUE, 1, filter_index) + p % wgt
              end if

              ! Inward current on d1 max surface
              cross_surface = .false.
              select case(n_dim)

              case (1)
                if (ijk0(d1) >  1 .and. ijk0(d1) <= m % dimension(d1) + 1) then
                  cross_surface = .true.
                end if

              case (2)
                if (ijk0(d1) >  1 .and. ijk0(d1) <= m % dimension(d1) + 1 .and.&
                     ijk0(d2) >= 1 .and. ijk0(d2) <= m % dimension(d2)) then
                  cross_surface = .true.
                end if

              case (3)
                if (ijk0(d1) >  1 .and. ijk0(d1) <= m % dimension(d1) + 1 .and.&
                     ijk0(d2) >= 1 .and. ijk0(d2) <= m % dimension(d2) .and. &
                     ijk0(d3) >= 1 .and. ijk0(d3) <= m % dimension(d3)) then
                  cross_surface = .true.
                end if
              end select

              ! If the particle crossed the surface, tally the current
              if (cross_surface) then
                ijk0(d1) = ijk0(d1) - 1
                filter_matches(i_filter_surf) % bins % data(1) = d1 * 4
                filter_matches(i_filter_mesh) % bins % data(1) = &
                     m % get_bin_from_indices(ijk0)
                filter_index = 1
                do k = 1, size(t % filter)
                  filter_index = filter_index + (filter_matches(t % &
                       filter(k)) % bins % data(1) - 1) * t % stride(k)
                end do
!$omp atomic
                t % results(RESULT_VALUE, 1, filter_index) = &
                     t % results(RESULT_VALUE, 1, filter_index) + p % wgt
                ijk0(d1) = ijk0(d1) + 1
              end if

              ijk0(d1) = ijk0(d1) - 1
              xyz_cross(d1) = xyz_cross(d1) - m % width(d1)
            end if
          end if
        end do

        ! Calculate new coordinates
        xyz0 = xyz0 + distance * uvw
      end do

      end associate
    end do TALLY_LOOP

  end subroutine score_surface_current

!===============================================================================
! APPLY_DERIVATIVE_TO_SCORE multiply the given score by its relative derivative
!===============================================================================

  subroutine apply_derivative_to_score(p, t, i_nuclide, atom_density, &
                                       score_bin, score)
    type(Particle),    intent(in)    :: p
    type(TallyObject), intent(in)    :: t
    integer,           intent(in)    :: i_nuclide
    real(8),           intent(in)    :: atom_density   ! atom/b-cm
    integer,           intent(in)    :: score_bin
    real(8),           intent(inout) :: score

    integer :: l
    logical :: scoring_diff_nuclide
    real(8) :: flux_deriv
    real(8) :: dsig_t, dsig_a, dsig_f, cum_dsig

    if (score == ZERO) return

    ! If our score was previously c then the new score is
    ! c * (1/f * d_f/d_p + 1/c * d_c/d_p)
    ! where (1/f * d_f/d_p) is the (logarithmic) flux derivative and p is the
    ! perturbated variable.

    associate(deriv => tally_derivs(t % deriv))
      flux_deriv = deriv % flux_deriv

      select case (tally_derivs(t % deriv) % variable)

      !=========================================================================
      ! Density derivative:
      ! c = Sigma_MT
      ! c = sigma_MT * N
      ! c = sigma_MT * rho * const
      ! d_c / d_rho = sigma_MT * const
      ! (1 / c) * (d_c / d_rho) = 1 / rho

      case (DIFF_DENSITY)
        select case (t % estimator)

        case (ESTIMATOR_ANALOG)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL, SCORE_SCATTER, SCORE_ABSORPTION, SCORE_FISSION, &
                SCORE_NU_FISSION)
            if (materials(p % material) % id == deriv % diff_material) then
              score = score * (flux_deriv + ONE &
                   / materials(p % material) % density_gpcc)
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id)))
          end select

        case (ESTIMATOR_COLLISION)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL, SCORE_SCATTER, SCORE_ABSORPTION, SCORE_FISSION, &
                SCORE_NU_FISSION)
            if (materials(p % material) % id == deriv % diff_material) then
              score = score * (flux_deriv + ONE &
                   / materials(p % material) % density_gpcc)
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id)))
          end select

        case default
            call fatal_error("Differential tallies are only implemented for &
                 &analog and collision estimators.")
        end select

      !=========================================================================
      ! Nuclide density derivative:
      ! If we are scoring a reaction rate for a single nuclide then
      ! c = Sigma_MT_i
      ! c = sigma_MT_i * N_i
      ! d_c / d_N_i = sigma_MT_i
      ! (1 / c) * (d_c / d_N_i) = 1 / N_i
      ! If the score is for the total material (i_nuclide = -1)
      ! c = Sum_i(Sigma_MT_i)
      ! d_c / d_N_i = sigma_MT_i
      ! (1 / c) * (d_c / d_N) = sigma_MT_i / Sigma_MT
      ! where i is the perturbed nuclide.

      case (DIFF_NUCLIDE_DENSITY)
        select case (t % estimator)

        case (ESTIMATOR_ANALOG)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL, SCORE_SCATTER, SCORE_ABSORPTION, SCORE_FISSION, &
                SCORE_NU_FISSION)
            if (materials(p % material) % id == deriv % diff_material &
                 .and. p % event_nuclide == deriv % diff_nuclide) then
              associate(mat => materials(p % material))
                ! Search for the index of the perturbed nuclide.
                do l = 1, mat % n_nuclides
                  if (mat % nuclide(l) == deriv % diff_nuclide) exit
                end do

                score = score * (flux_deriv &
                     + ONE / mat % atom_density(l))
              end associate
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id)))
          end select

        case (ESTIMATOR_COLLISION)
          scoring_diff_nuclide = &
               (materials(p % material) % id == deriv % diff_material) &
               .and. (i_nuclide == deriv % diff_nuclide)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL)
            if (i_nuclide == -1 .and. &
                 materials(p % material) % id == deriv % diff_material .and. &
                 material_xs % total /= ZERO) then
              score = score * (flux_deriv &
                   + micro_xs(deriv % diff_nuclide) % total &
                   / material_xs % total)
            else if (scoring_diff_nuclide .and. &
                 micro_xs(deriv % diff_nuclide) % total /= ZERO) then
              score = score * (flux_deriv + ONE / atom_density)
            else
              score = score * flux_deriv
            end if

          case (SCORE_SCATTER)
            if (i_nuclide == -1 .and. &
                 materials(p % material) % id == deriv % diff_material .and. &
                 material_xs % total - material_xs % absorption /= ZERO) then
              score = score * (flux_deriv &
                   + (micro_xs(deriv % diff_nuclide) % total &
                   - micro_xs(deriv % diff_nuclide) % absorption) &
                   / (material_xs % total - material_xs % absorption))
            else if (scoring_diff_nuclide .and. &
                 (micro_xs(deriv % diff_nuclide) % total &
                 - micro_xs(deriv % diff_nuclide) % absorption) /= ZERO) then
              score = score * (flux_deriv + ONE / atom_density)
            else
              score = score * flux_deriv
            end if

          case (SCORE_ABSORPTION)
            if (i_nuclide == -1 .and. &
                 materials(p % material) % id == deriv % diff_material .and. &
                 material_xs % absorption /= ZERO) then
              score = score * (flux_deriv &
                   + micro_xs(deriv % diff_nuclide) % absorption &
                   / material_xs % absorption )
            else if (scoring_diff_nuclide .and. &
                 micro_xs(deriv % diff_nuclide) % absorption /= ZERO) then
              score = score * (flux_deriv + ONE / atom_density)
            else
              score = score * flux_deriv
            end if

          case (SCORE_FISSION)
            if (i_nuclide == -1 .and. &
                 materials(p % material) % id == deriv % diff_material .and. &
                 material_xs % fission /= ZERO) then
              score = score * (flux_deriv &
                   + micro_xs(deriv % diff_nuclide) % fission &
                   / material_xs % fission)
            else if (scoring_diff_nuclide .and. &
                 micro_xs(deriv % diff_nuclide) % fission /= ZERO) then
              score = score * (flux_deriv + ONE / atom_density)
            else
              score = score * flux_deriv
            end if

          case (SCORE_NU_FISSION)
            if (i_nuclide == -1 .and. &
                 materials(p % material) % id == deriv % diff_material .and. &
                 material_xs % nu_fission /= ZERO) then
              score = score * (flux_deriv &
                   + micro_xs(deriv % diff_nuclide) % nu_fission &
                   / material_xs % nu_fission)
            else if (scoring_diff_nuclide .and. &
                 micro_xs(deriv % diff_nuclide) % nu_fission /= ZERO) then
              score = score * (flux_deriv + ONE / atom_density)
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id)))
          end select

        case default
            call fatal_error("Differential tallies are only implemented for &
                 &analog and collision estimators.")
        end select

      !=========================================================================
      ! Temperature derivative:
      ! If we are scoring a reaction rate for a single nuclide then
      ! c = Sigma_MT_i
      ! c = sigma_MT_i * N_i
      ! d_c / d_T = (d_sigma_Mt_i / d_T) * N_i
      ! (1 / c) * (d_c / d_T) = (d_sigma_MT_i / d_T) / sigma_MT_i
      ! If the score is for the total material (i_nuclide = -1)
      ! (1 / c) * (d_c / d_T) = Sum_i((d_sigma_MT_i / d_T) * N_i) / Sigma_MT_i
      ! where i is the perturbed nuclide.  The d_sigma_MT_i / d_T term is
      ! computed by multipole_deriv_eval.  It only works for the resolved
      ! resonance range and requires multipole data.

      case (DIFF_TEMPERATURE)
        select case (t % estimator)

        case (ESTIMATOR_ANALOG)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL)
            if (materials(p % material) % id == deriv % diff_material .and. &
                 micro_xs(p % event_nuclide) % total > ZERO) then
              associate(mat => materials(p % material))
                ! Search for the index of the perturbed nuclide.
                do l = 1, mat % n_nuclides
                  if (mat % nuclide(l) == p % event_nuclide) exit
                end do

                dsig_t = ZERO
                associate (nuc => nuclides(p % event_nuclide))
                  if (nuc % mp_present .and. &
                       p % last_E >= nuc % multipole % start_E .and. &
                       p % last_E <= nuc % multipole % end_E) then
                    call multipole_deriv_eval(nuc % multipole, p % last_E, &
                         p % sqrtkT, dsig_t, dsig_a, dsig_f)
                  end if
                end associate
                score = score * (flux_deriv &
                     + dsig_t * mat % atom_density(l) / material_xs % total)
              end associate
            else
              score = score * flux_deriv
            end if

          case (SCORE_SCATTER)
            if (materials(p % material) % id == deriv % diff_material .and. &
                 (micro_xs(p % event_nuclide) % total &
                 - micro_xs(p % event_nuclide) % absorption) > ZERO) then
              associate(mat => materials(p % material))
                ! Search for the index of the perturbed nuclide.
                do l = 1, mat % n_nuclides
                  if (mat % nuclide(l) == p % event_nuclide) exit
                end do

                dsig_t = ZERO
                dsig_a = ZERO
                associate (nuc => nuclides(p % event_nuclide))
                  if (nuc % mp_present .and. &
                       p % last_E >= nuc % multipole % start_E .and. &
                       p % last_E <= nuc % multipole % end_E) then
                    call multipole_deriv_eval(nuc % multipole, p % last_E, &
                         p % sqrtkT, dsig_t, dsig_a, dsig_f)
                  end if
                end associate
                score = score * (flux_deriv + (dsig_t - dsig_a) &
                     * mat % atom_density(l) / &
                     (material_xs % total - material_xs % absorption))
              end associate
            else
              score = score * flux_deriv
            end if

          case (SCORE_ABSORPTION)
            if (materials(p % material) % id == deriv % diff_material .and. &
                 micro_xs(p % event_nuclide) % absorption > ZERO) then
              associate(mat => materials(p % material))
                ! Search for the index of the perturbed nuclide.
                do l = 1, mat % n_nuclides
                  if (mat % nuclide(l) == p % event_nuclide) exit
                end do

                dsig_a = ZERO
                associate (nuc => nuclides(p % event_nuclide))
                  if (nuc % mp_present .and. &
                       p % last_E >= nuc % multipole % start_E .and. &
                       p % last_E <= nuc % multipole % end_E) then
                    call multipole_deriv_eval(nuc % multipole, p % last_E, &
                         p % sqrtkT, dsig_t, dsig_a, dsig_f)
                  end if
                end associate
                score = score * (flux_deriv + dsig_a * mat % atom_density(l) &
                                              / material_xs % absorption)
              end associate
            else
              score = score * flux_deriv
            end if

          case (SCORE_FISSION)
            if (materials(p % material) % id == deriv % diff_material .and. &
                 micro_xs(p % event_nuclide) % fission > ZERO) then
              associate(mat => materials(p % material))
                ! Search for the index of the perturbed nuclide.
                do l = 1, mat % n_nuclides
                  if (mat % nuclide(l) == p % event_nuclide) exit
                end do

                dsig_f = ZERO
                associate (nuc => nuclides(p % event_nuclide))
                  if (nuc % mp_present .and. &
                       p % last_E >= nuc % multipole % start_E .and. &
                       p % last_E <= nuc % multipole % end_E) then
                    call multipole_deriv_eval(nuc % multipole, p % last_E, &
                         p % sqrtkT, dsig_t, dsig_a, dsig_f)
                  end if
                end associate
                score = score * (flux_deriv &
                     + dsig_f * mat % atom_density(l) / material_xs % fission)
              end associate
            else
              score = score * flux_deriv
            end if

          case (SCORE_NU_FISSION)
            if (materials(p % material) % id == deriv % diff_material .and. &
                 micro_xs(p % event_nuclide) % nu_fission > ZERO) then
              associate(mat => materials(p % material))
                ! Search for the index of the perturbed nuclide.
                do l = 1, mat % n_nuclides
                  if (mat % nuclide(l) == p % event_nuclide) exit
                end do

                dsig_f = ZERO
                associate (nuc => nuclides(p % event_nuclide))
                  if (nuc % mp_present .and. &
                       p % last_E >= nuc % multipole % start_E .and. &
                       p % last_E <= nuc % multipole % end_E) then
                    call multipole_deriv_eval(nuc % multipole, p % last_E, &
                         p % sqrtkT, dsig_t, dsig_a, dsig_f)
                  end if
                end associate
                score = score * (flux_deriv &
                     + dsig_f * mat % atom_density(l) / material_xs % nu_fission&
                     * micro_xs(p % event_nuclide) % nu_fission &
                     / micro_xs(p % event_nuclide) % fission)
              end associate
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id)))
          end select

        case (ESTIMATOR_COLLISION)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL)
            if (i_nuclide == -1 .and. &
                 materials(p % material) % id == deriv % diff_material .and. &
                 material_xs % total > ZERO) then
              cum_dsig = ZERO
              associate(mat => materials(p % material))
                do l = 1, mat % n_nuclides
                  associate (nuc => nuclides(mat % nuclide(l)))
                    if (nuc % mp_present .and. &
                         p % last_E >= nuc % multipole % start_E .and. &
                         p % last_E <= nuc % multipole % end_E .and. &
                         micro_xs(mat % nuclide(l)) % total > ZERO) then
                      call multipole_deriv_eval(nuc % multipole, p % last_E, &
                           p % sqrtkT, dsig_t, dsig_a, dsig_f)
                      cum_dsig = cum_dsig + dsig_t * mat % atom_density(l)
                    end if
                  end associate
                end do
              end associate
              score = score * (flux_deriv &
                   + cum_dsig / material_xs % total)
            else if (materials(p % material) % id == deriv % diff_material &
                 .and. material_xs % total > ZERO) then
              dsig_t = ZERO
              associate (nuc => nuclides(i_nuclide))
                if (nuc % mp_present .and. &
                     p % last_E >= nuc % multipole % start_E .and. &
                     p % last_E <= nuc % multipole % end_E) then
                  call multipole_deriv_eval(nuc % multipole, p % last_E, &
                       p % sqrtkT, dsig_t, dsig_a, dsig_f)
                end if
              end associate
              score = score * (flux_deriv &
                   + dsig_t / micro_xs(i_nuclide) % total)
            else
              score = score * flux_deriv
            end if

          case (SCORE_SCATTER)
            if (i_nuclide == -1 .and. &
                 materials(p % material) % id == deriv % diff_material .and. &
                 (material_xs % total - material_xs % absorption) > ZERO) then
              cum_dsig = ZERO
              associate(mat => materials(p % material))
                do l = 1, mat % n_nuclides
                  associate (nuc => nuclides(mat % nuclide(l)))
                    if (nuc % mp_present .and. &
                         p % last_E >= nuc % multipole % start_E .and. &
                         p % last_E <= nuc % multipole % end_E .and. &
                         (micro_xs(mat % nuclide(l)) % total &
                         - micro_xs(mat % nuclide(l)) % absorption) > ZERO) then
                      call multipole_deriv_eval(nuc % multipole, p % last_E, &
                           p % sqrtkT, dsig_t, dsig_a, dsig_f)
                      cum_dsig = cum_dsig &
                           + (dsig_t - dsig_a) * mat % atom_density(l)
                    end if
                  end associate
                end do
              end associate
              score = score * (flux_deriv + cum_dsig &
                   / (material_xs % total - material_xs % absorption))
            else if ( materials(p % material) % id == deriv % diff_material &
                 .and. (material_xs % total - material_xs % absorption) > ZERO)&
                 then
              dsig_t = ZERO
              dsig_a = ZERO
              associate (nuc => nuclides(i_nuclide))
                if (nuc % mp_present .and. &
                     p % last_E >= nuc % multipole % start_E .and. &
                     p % last_E <= nuc % multipole % end_E) then
                  call multipole_deriv_eval(nuc % multipole, p % last_E, &
                       p % sqrtkT, dsig_t, dsig_a, dsig_f)
                end if
              end associate
              score = score * (flux_deriv + (dsig_t - dsig_a) &
                   / (micro_xs(i_nuclide) % total &
                   - micro_xs(i_nuclide) % absorption))
            else
              score = score * flux_deriv
            end if

          case (SCORE_ABSORPTION)
            if (i_nuclide == -1 .and. &
                 materials(p % material) % id == deriv % diff_material .and. &
                 material_xs % absorption > ZERO) then
              cum_dsig = ZERO
              associate(mat => materials(p % material))
                do l = 1, mat % n_nuclides
                  associate (nuc => nuclides(mat % nuclide(l)))
                    if (nuc % mp_present .and. &
                         p % last_E >= nuc % multipole % start_E .and. &
                         p % last_E <= nuc % multipole % end_E .and. &
                         micro_xs(mat % nuclide(l)) % absorption > ZERO) then
                      call multipole_deriv_eval(nuc % multipole, p % last_E, &
                           p % sqrtkT, dsig_t, dsig_a, dsig_f)
                      cum_dsig = cum_dsig + dsig_a * mat % atom_density(l)
                    end if
                  end associate
                end do
              end associate
              score = score * (flux_deriv &
                   + cum_dsig / material_xs % absorption)
            else if (materials(p % material) % id == deriv % diff_material &
                 .and. material_xs % absorption > ZERO) then
              dsig_a = ZERO
              associate (nuc => nuclides(i_nuclide))
                if (nuc % mp_present .and. &
                     p % last_E >= nuc % multipole % start_E .and. &
                     p % last_E <= nuc % multipole % end_E) then
                  call multipole_deriv_eval(nuc % multipole, p % last_E, &
                       p % sqrtkT, dsig_t, dsig_a, dsig_f)
                end if
              end associate
              score = score * (flux_deriv &
                   + dsig_a / micro_xs(i_nuclide) % absorption)
            else
              score = score * flux_deriv
            end if

          case (SCORE_FISSION)
            if (i_nuclide == -1 .and. &
                 materials(p % material) % id == deriv % diff_material .and. &
                 material_xs % fission > ZERO) then
              cum_dsig = ZERO
              associate(mat => materials(p % material))
                do l = 1, mat % n_nuclides
                  associate (nuc => nuclides(mat % nuclide(l)))
                    if (nuc % mp_present .and. &
                         p % last_E >= nuc % multipole % start_E .and. &
                         p % last_E <= nuc % multipole % end_E .and. &
                         micro_xs(mat % nuclide(l)) % fission > ZERO) then
                      call multipole_deriv_eval(nuc % multipole, p % last_E, &
                           p % sqrtkT, dsig_t, dsig_a, dsig_f)
                      cum_dsig = cum_dsig + dsig_f * mat % atom_density(l)
                    end if
                  end associate
                end do
              end associate
              score = score * (flux_deriv &
                   + cum_dsig / material_xs % fission)
            else if (materials(p % material) % id == deriv % diff_material &
                 .and. material_xs % fission > ZERO) then
              dsig_f = ZERO
              associate (nuc => nuclides(i_nuclide))
                if (nuc % mp_present .and. &
                     p % last_E >= nuc % multipole % start_E .and. &
                     p % last_E <= nuc % multipole % end_E) then
                  call multipole_deriv_eval(nuc % multipole, p % last_E, &
                       p % sqrtkT, dsig_t, dsig_a, dsig_f)
                end if
              end associate
              score = score * (flux_deriv &
                   + dsig_f / micro_xs(i_nuclide) % fission)
            else
              score = score * flux_deriv
            end if

          case (SCORE_NU_FISSION)
            if (i_nuclide == -1 .and. &
                 materials(p % material) % id == deriv % diff_material .and. &
                 material_xs % nu_fission > ZERO) then
              cum_dsig = ZERO
              associate(mat => materials(p % material))
                do l = 1, mat % n_nuclides
                  associate (nuc => nuclides(mat % nuclide(l)))
                    if (nuc % mp_present .and. &
                         p % last_E >= nuc % multipole % start_E .and. &
                         p % last_E <= nuc % multipole % end_E .and. &
                         micro_xs(mat % nuclide(l)) % nu_fission > ZERO) then
                      call multipole_deriv_eval(nuc % multipole, p % last_E, &
                           p % sqrtkT, dsig_t, dsig_a, dsig_f)
                      cum_dsig = cum_dsig + dsig_f * mat % atom_density(l) &
                           * micro_xs(mat % nuclide(l)) % nu_fission &
                           / micro_xs(mat % nuclide(l)) % fission
                    end if
                  end associate
                end do
              end associate
              score = score * (flux_deriv &
                   + cum_dsig / material_xs % nu_fission)
            else if (materials(p % material) % id == deriv % diff_material &
                 .and. material_xs % nu_fission > ZERO) then
              dsig_f = ZERO
              associate (nuc => nuclides(i_nuclide))
                if (nuc % mp_present .and. &
                     p % last_E >= nuc % multipole % start_E .and. &
                     p % last_E <= nuc % multipole % end_E) then
                  call multipole_deriv_eval(nuc % multipole, p % last_E, &
                       p % sqrtkT, dsig_t, dsig_a, dsig_f)
                end if
              end associate
              score = score * (flux_deriv &
                   + dsig_f / micro_xs(i_nuclide) % fission)
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id)))
          end select

        case default
            call fatal_error("Differential tallies are only implemented for &
                 &analog and collision estimators.")
        end select
      end select
    end associate
  end subroutine apply_derivative_to_score

!===============================================================================
! SCORE_TRACK_DERIVATIVE Adjust flux derivatives on differential tallies to
! account for a neutron travelling through a perturbed material.
!===============================================================================

  subroutine score_track_derivative(p, distance)
    type(Particle), intent(in) :: p
    real(8),        intent(in) :: distance ! Neutron flight distance

    integer :: i, l
    real(8) :: dsig_t, dsig_a, dsig_f

    ! A void material cannot be perturbed so it will not affect flux derivatives
    if (p % material == MATERIAL_VOID) return

    do i = 1, size(tally_derivs)
      associate(deriv => tally_derivs(i))
        select case (deriv % variable)

        case (DIFF_DENSITY)
          associate (mat => materials(p % material))
            if (mat % id == deriv % diff_material) then
              ! phi is proportional to e^(-Sigma_tot * dist)
              ! (1 / phi) * (d_phi / d_rho) = - (d_Sigma_tot / d_rho) * dist
              ! (1 / phi) * (d_phi / d_rho) = - Sigma_tot / rho * dist
              deriv % flux_deriv = deriv % flux_deriv &
                   - distance * material_xs % total / mat % density_gpcc
            end if
          end associate

        case (DIFF_NUCLIDE_DENSITY)
          associate (mat => materials(p % material))
            if (mat % id == deriv % diff_material) then
              ! phi is proportional to e^(-Sigma_tot * dist)
              ! (1 / phi) * (d_phi / d_N) = - (d_Sigma_tot / d_N) * dist
              ! (1 / phi) * (d_phi / d_N) = - sigma_tot * dist
              deriv % flux_deriv = deriv % flux_deriv &
                   - distance * micro_xs(deriv % diff_nuclide) % total
            end if
          end associate

        case (DIFF_TEMPERATURE)
          associate (mat => materials(p % material))
            if (mat % id == deriv % diff_material) then
              do l=1, mat % n_nuclides
                associate (nuc => nuclides(mat % nuclide(l)))
                  if (nuc % mp_present .and. &
                       p % E >= nuc % multipole % start_E .and. &
                       p % E <= nuc % multipole % end_E) then
                    ! phi is proportional to e^(-Sigma_tot * dist)
                    ! (1 / phi) * (d_phi / d_T) = - (d_Sigma_tot / d_T) * dist
                    ! (1 / phi) * (d_phi / d_T) = - N (d_sigma_tot / d_T) * dist
                    call multipole_deriv_eval(nuc % multipole, p % E, &
                         p % sqrtkT, dsig_t, dsig_a, dsig_f)
                    deriv % flux_deriv = deriv % flux_deriv &
                         - distance * dsig_t * mat % atom_density(l)
                  end if
                end associate
              end do
            end if
          end associate
        end select
      end associate
    end do
  end subroutine score_track_derivative

!===============================================================================
! SCORE_COLLISION_DERIVATIVE Adjust flux derivatives on differential tallies to
! account for a neutron scattering in the perturbed material.  Note that this
! subroutine will be called after absorption events in addition to scattering
! events, but any flux derivatives scored after an absorption will never be
! tallied.  This is because the order of operations is
! 1. Particle is moved.
! 2. score_track_derivative is called.
! 3. Collision physics are computed, and the particle is labeled absorbed.
! 4. Analog- and collision-estimated tallies are scored.
! 5. This subroutine is called.
! 6. Particle is killed and no more tallies are scored.
! Hence, it is safe to assume that only derivative of the scattering cross
! section need to be computed here.
!===============================================================================

  subroutine score_collision_derivative(p)
    type(Particle), intent(in) :: p

    integer :: i, j, l
    real(8) :: dsig_t, dsig_a, dsig_f

    ! A void material cannot be perturbed so it will not affect flux derivatives
    if (p % material == MATERIAL_VOID) return

    do i = 1, size(tally_derivs)
      associate(deriv => tally_derivs(i))
        select case (deriv % variable)

        case (DIFF_DENSITY)
          associate (mat => materials(p % material))
            if (mat % id == deriv % diff_material) then
              ! phi is proportional to Sigma_s
              ! (1 / phi) * (d_phi / d_rho) = (d_Sigma_s / d_rho) / Sigma_s
              ! (1 / phi) * (d_phi / d_rho) = 1 / rho
              deriv % flux_deriv = deriv % flux_deriv &
                   + ONE / mat % density_gpcc
            end if
          end associate

        case (DIFF_NUCLIDE_DENSITY)
          associate (mat => materials(p % material))
            if (mat % id == deriv % diff_material &
                 .and. p % event_nuclide == deriv % diff_nuclide) then
              ! Find the index in this material for the diff_nuclide.
              do j = 1, mat % n_nuclides
                if (mat % nuclide(j) == deriv % diff_nuclide) exit
              end do
              ! Make sure we found the nuclide.
              if (mat % nuclide(j) /= deriv % diff_nuclide) then
                call fatal_error("Couldn't find the right nuclide.")
              end if
              ! phi is proportional to Sigma_s
              ! (1 / phi) * (d_phi / d_N) = (d_Sigma_s / d_N) / Sigma_s
              ! (1 / phi) * (d_phi / d_N) = sigma_s / Sigma_s
              ! (1 / phi) * (d_phi / d_N) = 1 / N
              deriv % flux_deriv = deriv % flux_deriv &
                   + ONE / mat % atom_density(j)
            end if
          end associate

        case (DIFF_TEMPERATURE)
          associate (mat => materials(p % material))
            if (mat % id == deriv % diff_material) then
              do l=1, mat % n_nuclides
                associate (nuc => nuclides(mat % nuclide(l)))
                  if (mat % nuclide(l) == p % event_nuclide .and. &
                       nuc % mp_present .and. &
                       p % last_E >= nuc % multipole % start_E .and. &
                       p % last_E <= nuc % multipole % end_E) then
                    ! phi is proportional to Sigma_s
                    ! (1 / phi) * (d_phi / d_T) = (d_Sigma_s / d_T) / Sigma_s
                    ! (1 / phi) * (d_phi / d_T) = (d_sigma_s / d_T) / sigma_s
                    call multipole_deriv_eval(nuc % multipole, p % last_E, &
                         p % sqrtkT, dsig_t, dsig_a, dsig_f)
                    deriv % flux_deriv = deriv % flux_deriv + (dsig_t - dsig_a)&
                         / (micro_xs(mat % nuclide(l)) % total &
                         - micro_xs(mat % nuclide(l)) % absorption)
                    ! Note that this is an approximation!  The real scattering
                    ! cross section is Sigma_s(E'->E, uvw'->uvw) =
                    ! Sigma_s(E') * P(E'->E, uvw'->uvw).  We are assuming that
                    ! d_P(E'->E, uvw'->uvw) / d_T = 0 and only computing
                    ! d_S(E') / d_T.  Using this approximation in the vicinity
                    ! of low-energy resonances causes errors (~2-5% for PWR
                    ! pincell eigenvalue derivatives).
                  end if
                end associate
              end do
            end if
          end associate
        end select
      end associate
    end do
  end subroutine score_collision_derivative

!===============================================================================
! ZERO_FLUX_DERIVS Set the flux derivatives on differential tallies to zero.
!===============================================================================

  subroutine zero_flux_derivs()
    integer :: i
    do i = 1, size(tally_derivs)
      tally_derivs(i) % flux_deriv = ZERO
    end do
  end subroutine zero_flux_derivs

!===============================================================================
! ACCUMULATE_TALLIES accumulates the sum of the contributions from each history
! within the batch to a new random variable
!===============================================================================

  subroutine accumulate_tallies()

    integer :: i
    real(C_DOUBLE) :: k_col ! Copy of batch collision estimate of keff
    real(C_DOUBLE) :: k_abs ! Copy of batch absorption estimate of keff
    real(C_DOUBLE) :: k_tra ! Copy of batch tracklength estimate of keff
    real(C_DOUBLE) :: val

#ifdef OPENMC_MPI
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
        call tallies(active_tallies % data(i)) % obj % accumulate()
      end do

      if (run_mode == MODE_EIGENVALUE) then
        if (current_batch > n_inactive) then
          ! Accumulate products of different estimators of k
          k_col = global_tallies(RESULT_VALUE, K_COLLISION) / total_weight
          k_abs = global_tallies(RESULT_VALUE, K_ABSORPTION) / total_weight
          k_tra = global_tallies(RESULT_VALUE, K_TRACKLENGTH) / total_weight
          k_col_abs = k_col_abs + k_col * k_abs
          k_col_tra = k_col_tra + k_col * k_tra
          k_abs_tra = k_abs_tra + k_abs * k_tra
        end if
      end if

      ! Accumulate results for global tallies
      do i = 1, size(global_tallies, 2)
        val = global_tallies(RESULT_VALUE, i)/total_weight
        global_tallies(RESULT_VALUE, i) = ZERO

        global_tallies(RESULT_SUM, i) = global_tallies(RESULT_SUM, i) + val
        global_tallies(RESULT_SUM_SQ, i) = &
             global_tallies(RESULT_SUM_SQ, i) + val*val
      end do
    end if

  end subroutine accumulate_tallies

!===============================================================================
! REDUCE_TALLY_RESULTS collects all the results from tallies onto one processor
!===============================================================================

#ifdef OPENMC_MPI
  subroutine reduce_tally_results()

    integer :: i
    integer :: n       ! number of filter bins
    integer :: m       ! number of score bins
    integer :: n_bins  ! total number of bins
    integer :: mpi_err ! MPI error code
    real(C_DOUBLE), allocatable :: tally_temp(:,:)  ! contiguous array of results
    real(C_DOUBLE), allocatable :: tally_temp2(:,:) ! reduced contiguous results
    real(C_DOUBLE) :: temp(N_GLOBAL_TALLIES), temp2(N_GLOBAL_TALLIES)

    do i = 1, active_tallies % size()
      associate (t => tallies(active_tallies % data(i)) % obj)

        m = size(t % results, 2)
        n = size(t % results, 3)
        n_bins = m*n

        allocate(tally_temp(m,n), tally_temp2(m,n))

        ! Reduce contiguous set of tally results
        tally_temp = t % results(RESULT_VALUE,:,:)
        call MPI_REDUCE(tally_temp, tally_temp2, n_bins, MPI_DOUBLE, &
             MPI_SUM, 0, mpi_intracomm, mpi_err)

        if (master) then
          ! Transfer values to value on master
          t % results(RESULT_VALUE,:,:) = tally_temp2
        else
          ! Reset value on other processors
          t % results(RESULT_VALUE,:,:) = ZERO
        end if

        deallocate(tally_temp, tally_temp2)
      end associate
    end do

    ! Reduce global tallies onto master
    temp = global_tallies(RESULT_VALUE, :)
    call MPI_REDUCE(temp, temp2, N_GLOBAL_TALLIES, MPI_DOUBLE, MPI_SUM, &
         0, mpi_intracomm, mpi_err)
    if (master) then
      global_tallies(RESULT_VALUE, :) = temp2
    else
      global_tallies(RESULT_VALUE, :) = ZERO
    end if

    ! We also need to determine the total starting weight of particles from the
    ! last realization
    temp(1) = total_weight
    call MPI_REDUCE(temp, total_weight, 1, MPI_REAL8, MPI_SUM, &
         0, mpi_intracomm, mpi_err)

  end subroutine reduce_tally_results
#endif

!===============================================================================
! SETUP_ACTIVE_TALLIES
!===============================================================================

  subroutine setup_active_tallies()

    integer :: i ! loop counter

    call active_tallies % clear()
    call active_analog_tallies % clear()
    call active_collision_tallies % clear()
    call active_tracklength_tallies % clear()
    call active_surface_tallies % clear()
    call active_current_tallies % clear()

    do i = 1, n_tallies
      associate (t => tallies(i) % obj)
        if (t % active) then
          ! Add tally to active tallies
          call active_tallies % push_back(i)

          ! Check what type of tally this is and add it to the appropriate list
          if (t % type == TALLY_VOLUME) then
            if (t % estimator == ESTIMATOR_ANALOG) then
              call active_analog_tallies % push_back(i)
            elseif (t % estimator == ESTIMATOR_TRACKLENGTH) then
              call active_tracklength_tallies % push_back(i)
            elseif (t % estimator == ESTIMATOR_COLLISION) then
              call active_collision_tallies % push_back(i)
            end if
          elseif (t % type == TALLY_MESH_CURRENT) then
            call active_current_tallies % push_back(i)
          elseif (t % type == TALLY_SURFACE) then
            call active_surface_tallies % push_back(i)
          end if

          ! Check if tally contains depletion reactions and if so, set flag
          if (t % depletion_rx) need_depletion_rx = .true.
        end if
      end associate
    end do

  end subroutine setup_active_tallies

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_tally_set_type(index, type) result(err) bind(C)
    ! Set the type of the tally
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR), intent(in) :: type(*)
    integer(C_INT) :: err

    integer(C_INT32_T) :: empty(0)
    character(:), allocatable :: type_

    ! Convert C string to Fortran string
    type_ = to_f_string(type)

    err = 0
    if (index >= 1 .and. index <= n_tallies) then
      if (allocated(tallies(index) % obj)) then
        err = E_ALLOCATE
        call set_errmsg("Tally type has already been set.")
      else
        select case (type_)
        case ('generic')
          allocate(TallyObject :: tallies(index) % obj)
        case default
          err = E_UNASSIGNED
          call set_errmsg("Unknown tally type: " // trim(type_))
        end select

        ! When a tally is allocated, set it to have 0 filters
        err = tallies(index) % obj % set_filters(empty)
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in tallies array is out of bounds.")
    end if
  end function openmc_tally_set_type

end module tally
