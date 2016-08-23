module tally

  use constants
  use error,            only: fatal_error
  use geometry_header
  use global
  use math,             only: t_percentile, calc_pn, calc_rn
  use mesh,             only: get_mesh_bin, bin_to_mesh_indices, &
                              get_mesh_indices, mesh_indices_to_bin, &
                              mesh_intersects_2d, mesh_intersects_3d
  use mesh_header,      only: RegularMesh
  use output,           only: header
  use particle_header,  only: LocalCoord, Particle
  use search,           only: binary_search
  use string,           only: to_str
  use tally_header,     only: TallyResult
  use tally_filter

#ifdef MPI
  use message_passing
#endif

  implicit none

  integer :: position(N_FILTER_TYPES - 3) = 0 ! Tally map positioning array

!$omp threadprivate(position)

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
    integer :: i_nuc                ! index in nuclides array (from material)
    integer :: i_energy             ! index in nuclide energy grid
    integer :: score_bin            ! scoring bin, e.g. SCORE_FLUX
    integer :: score_index          ! scoring bin index
    integer :: d                    ! delayed neutron index
    integer :: d_bin                ! delayed group bin index
    integer :: dg_filter            ! index of delayed group filter
    real(8) :: yield                ! delayed neutron yield
    real(8) :: atom_density_        ! atom/b-cm
    real(8) :: f                    ! interpolation factor
    real(8) :: score                ! analog tally score
    real(8) :: E                    ! particle energy

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
            score = p % last_wgt + p % absorb_wgt * flux
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
        ! make sure the correct energy is used
        if (t % estimator == ESTIMATOR_TRACKLENGTH) then
          E = p % E
        else
          E = p % last_E
        end if

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
               / (sqrt(TWO * E / (MASS_NEUTRON_MEV)) * C_LIGHT * 100.0_8) * flux

        else
          ! For inverse velocity, we don't need a cross section. The velocity is
          ! in units of cm/s.
          score = flux / (sqrt(TWO * E / (MASS_NEUTRON_MEV)) * C_LIGHT * 100.0_8)
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
          m = nuclides(p % event_nuclide) % reaction_index % &
               get_key(p % event_MT)

          ! Get yield and apply to score
          associate (rxn => nuclides(p % event_nuclide) % reactions(m))
            score = p % last_wgt * flux &
                 * rxn % products(1) % yield % evaluate(p % last_E)
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
          ! Don't waste time on very common reactions we know have multiplicities
          ! of one.
          score = p % last_wgt * flux
        else
          m = nuclides(p%event_nuclide)%reaction_index% &
               get_key(p % event_MT)

          ! Get yield and apply to score
          associate (rxn => nuclides(p % event_nuclide) % reactions(m))
            score = p % last_wgt * flux &
                 * rxn % products(1) % yield % evaluate(p % last_E)
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
          ! Don't waste time on very common reactions we know have multiplicities
          ! of one.
          score = p % last_wgt * flux
        else
          m = nuclides(p%event_nuclide)%reaction_index% &
               get_key(p % event_MT)

          ! Get yield and apply to score
          associate (rxn => nuclides(p%event_nuclide)%reactions(m))
            score = p % last_wgt * flux &
                 * rxn % products(1) % yield % evaluate(p % last_E)
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
              call score_fission_eout_ce(p, t, score_index, score_bin)
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
              call score_fission_eout_ce(p, t, score_index, score_bin)
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
                     / micro_xs(p % event_nuclide) % absorption
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
                 / real(p % n_bank, 8))
          end if

        else
          ! make sure the correct energy is used
          if (t % estimator == ESTIMATOR_TRACKLENGTH) then
            E = p % E
          else
            E = p % last_E
          end if

          if (i_nuclide > 0) then
              score = micro_xs(i_nuclide) % fission * nuclides(i_nuclide) % &
                   nu(E, EMISSION_PROMPT) * atom_density * flux
          else

            score = ZERO

            ! Loop over all nuclides in the current material
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


      case (SCORE_DELAYED_NU_FISSION)

        ! make sure the correct energy is used
        if (t % estimator == ESTIMATOR_TRACKLENGTH) then
          E = p % E
        else
          E = p % last_E
        end if

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
              call score_fission_eout_ce(p, t, score_index, score_bin)
              cycle SCORE_LOOP
            end if
          end if
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! delayed-nu-fission
            if (micro_xs(p % event_nuclide) % absorption > ZERO) then

              ! Check if the delayed group filter is present
              if (dg_filter > 0) then
                select type(filt => t % filters(dg_filter) % obj)
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
                         / micro_xs(p % event_nuclide) % absorption
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
                     / micro_xs(p % event_nuclide) % absorption
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
              select type(filt => t % filters(dg_filter) % obj)
              type is (DelayedGroupFilter)

                ! Loop over all delayed group bins and tally to them
                ! individually
                do d_bin = 1, filt % n_bins

                  ! Get the delayed group for this bin
                  d = filt % groups(d_bin)

                  ! Compute the score and tally to bin
                  score = keff * p % wgt_bank / p % n_bank * p % n_delayed_bank(d)
                  call score_fission_delayed_dg(t, d_bin, score, score_index)
                end do
                cycle SCORE_LOOP
              end select
            else

              ! Add the contribution from all delayed groups
              score = keff * p % wgt_bank / p % n_bank * sum(p % n_delayed_bank)
            end if
          end if
        else

          ! Check if tally is on a single nuclide
          if (i_nuclide > 0) then

            ! Check if the delayed group filter is present
            if (dg_filter > 0) then
              select type(filt => t % filters(dg_filter) % obj)
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
              select type(filt => t % filters(dg_filter) % obj)
              type is (DelayedGroupFilter)

                ! Loop over all nuclides in the current material
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
                    score = micro_xs(i_nuc) % fission * yield * atom_density_ * flux
                    call score_fission_delayed_dg(t, d_bin, score, score_index)
                  end do
                end do
                cycle SCORE_LOOP
              end select
            else

              score = ZERO

              ! Loop over all nuclides in the current material
              do l = 1, materials(p % material) % n_nuclides

                ! Get atom density
                atom_density_ = materials(p % material) % atom_density(l)

                ! Get index in nuclides array
                i_nuc = materials(p % material) % nuclide(l)

                ! Accumulate the contribution from each nuclide
                score = score + micro_xs(i_nuc) % fission * nuclides(i_nuc) % &
                     nu(E, EMISSION_DELAYED) * atom_density_ * flux
              end do
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
            score = micro_xs(i_nuclide) % elastic * atom_density * flux
          else
            score = material_xs % elastic * flux
          end if
        end if

      case (SCORE_FISS_Q_PROMPT)
        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! fission scaled by Q-value
            associate (nuc => nuclides(p % event_nuclide))
              if (micro_xs(p % event_nuclide) % absorption > ZERO .and. &
                   allocated(nuc % fission_q_prompt)) then
                score = p % absorb_wgt &
                     * nuc % fission_q_prompt % evaluate(p % last_E) &
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
                     * nuc % fission_q_prompt % evaluate(p % last_E) &
                     * micro_xs(p % event_nuclide) % fission &
                     / micro_xs(p % event_nuclide) % absorption * flux
              end if
            end associate
          end if

        else
          if (t % estimator == ESTIMATOR_COLLISION) then
            E = p % last_E
          else
            E = p % E
          end if

          if (i_nuclide > 0) then
            if (allocated(nuclides(i_nuclide) % fission_q_prompt)) then
              score = micro_xs(i_nuclide) % fission * atom_density * flux &
                      * nuclides(i_nuclide) % fission_q_prompt % evaluate(E)
            else
              score = ZERO
            end if
          else
            score = ZERO
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

      case (SCORE_FISS_Q_RECOV)
        if (t % estimator == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! fission scaled by Q-value
            associate (nuc => nuclides(p % event_nuclide))
              if (micro_xs(p % event_nuclide) % absorption > ZERO .and. &
                   allocated(nuc % fission_q_recov)) then
                score = p % absorb_wgt &
                     * nuc % fission_q_recov % evaluate(p % last_E) &
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
                     * nuc % fission_q_recov % evaluate(p % last_E) &
                     * micro_xs(p % event_nuclide) % fission &
                     / micro_xs(p % event_nuclide) % absorption * flux
              end if
            end associate
          end if

        else
          if (t % estimator == ESTIMATOR_COLLISION) then
            E = p % last_E
          else
            E = p % E
          end if

          if (i_nuclide > 0) then
            if (allocated(nuclides(i_nuclide) % fission_q_recov)) then
              score = micro_xs(i_nuclide) % fission * atom_density * flux &
                      * nuclides(i_nuclide) % fission_q_recov % evaluate(E)
            else
              score = ZERO
            end if
          else
            score = ZERO
            do l = 1, materials(p % material) % n_nuclides
              atom_density_ = materials(p % material) % atom_density(l)
              i_nuc = materials(p % material) % nuclide(l)
              if (allocated(nuclides(i_nuc) % fission_q_recov)) then
                score = score + micro_xs(i_nuc) % fission * atom_density_ &
                        * flux * nuclides(i_nuc) % fission_q_recov % evaluate(E)
              end if
            end do
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
              if (nuclides(i_nuclide)%reaction_index%has_key(score_bin)) then
                m = nuclides(i_nuclide)%reaction_index%get_key(score_bin)
                associate (rxn => nuclides(i_nuclide) % reactions(m))

                  ! Retrieve index on nuclide energy grid and interpolation
                  ! factor
                  i_energy = micro_xs(i_nuclide) % index_grid
                  f = micro_xs(i_nuclide) % interp_factor
                  if (i_energy >= rxn % threshold) then
                    score = ((ONE - f) * rxn % sigma(i_energy - &
                         rxn%threshold + 1) + f * rxn % sigma(i_energy - &
                         rxn%threshold + 2)) * atom_density * flux
                  end if
                end associate
              end if

            else
              do l = 1, materials(p % material) % n_nuclides
                ! Get atom density
                atom_density_ = materials(p % material) % atom_density(l)

                ! Get index in nuclides array
                i_nuc = materials(p % material) % nuclide(l)

                if (nuclides(i_nuc)%reaction_index%has_key(score_bin)) then
                  m = nuclides(i_nuc)%reaction_index%get_key(score_bin)
                  associate (rxn => nuclides(i_nuc) % reactions(m))
                    ! Retrieve index on nuclide energy grid and interpolation
                    ! factor
                    i_energy = micro_xs(i_nuc) % index_grid
                    f = micro_xs(i_nuc) % interp_factor
                    if (i_energy >= rxn % threshold) then
                      score = score + ((ONE - f) * rxn % sigma(i_energy - &
                           rxn%threshold + 1) + f * rxn % sigma(i_energy - &
                           rxn%threshold + 2)) * atom_density_ * flux
                    end if
                  end associate
                end if
              end do
            end if

          else
            call fatal_error("Invalid score type on tally " &
                 // to_str(t % id) // ".")
          end if
        end if

      end select

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
    real(8) :: score                ! analog tally score
    real(8) :: p_uvw(3)             ! Particle's current uvw
    integer :: p_g                  ! Particle group to use for getting info
                                    ! to tally with.
    class(Mgxs), pointer :: matxs
    class(Mgxs), pointer :: nucxs

    ! Set the direction and group to use with get_xs
    ! this only depends on if we
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
          score = score * inverse_velocities(p_g) / material_xs % total * flux

        else
          ! For inverse velocity, we need no cross section
          score = flux * inverse_velocities(p_g)
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
            score = flux * material_xs % fission

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
              call score_fission_eout_mg(p, t, score_index, i_nuclide, &
                                         atom_density)
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
                   nucxs % get_xs('nu_fission', p_g, UVW=p_uvw) / &
                   matxs % get_xs('absorption', p_g, UVW=p_uvw)
            else
              score = score * &
                   matxs % get_xs('nu_fission', p_g, UVW=p_uvw) / &
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
            score = nucxs % get_xs('nu_fission', p_g, UVW=p_uvw) * &
                 atom_density * flux
          else
            score = material_xs % nu_fission * flux
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
                 nucxs % get_xs('kappa_fission', p_g, UVW=p_uvw) / &
                 matxs % get_xs('absorption', p_g, UVW=p_uvw)
          else
            score = score * &
                 matxs % get_xs('kappa_fission', p_g, UVW=p_uvw) / &
                 matxs % get_xs('absorption', p_g, UVW=p_uvw)
          end if
        else
          if (i_nuclide > 0) then
            score = flux * atom_density * &
                 nucxs % get_xs('kappa_fission', p_g, UVW=p_uvw)
          else
            score = flux * matxs % get_xs('kappa_fission', p_g, UVW=p_uvw)

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

    nullify(matxs,nucxs)
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
      t % results(score_index, filter_index) % value = &
           t % results(score_index, filter_index) % value + score


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
        t % results(score_index: score_index + num_nm - 1, filter_index) &
             % value = t &
             % results(score_index: score_index + num_nm - 1, filter_index)&
             % value &
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
        t % results(score_index: score_index + num_nm - 1, filter_index) &
             % value = t &
             % results(score_index: score_index + num_nm - 1, filter_index)&
             % value &
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
        t % results(score_index, filter_index) % value = &
             t % results(score_index, filter_index) % value &
             + score * calc_pn(n, p % mu)
      end do
      i = i + t % moment_order(i)


    case default
!$omp atomic
      t % results(score_index, filter_index) % value = &
           t % results(score_index, filter_index) % value + score


    end select

  end subroutine expand_and_score

!===============================================================================
! SCORE_ALL_NUCLIDES tallies individual nuclide reaction rates specifically when
! the user requests <nuclides>all</nuclides>.
!===============================================================================

  subroutine score_all_nuclides(p, i_tally, flux, filter_index)

    type(Particle), intent(in) :: p
    integer,        intent(in) :: i_tally
    real(8),        intent(in) :: flux
    integer,        intent(in) :: filter_index

    integer :: i             ! loop index for nuclides in material
    integer :: i_nuclide     ! index in nuclides array
    real(8) :: atom_density  ! atom density of single nuclide in atom/b-cm
    type(TallyObject), pointer :: t
    type(Material),    pointer :: mat

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

      ! Determine score for each bin
      call score_general(p, t, (i_nuclide-1)*t % n_score_bins, filter_index, &
           i_nuclide, atom_density, flux)

    end do NUCLIDE_LOOP

    ! ==========================================================================
    ! SCORE TOTAL MATERIAL REACTION RATES

    i_nuclide = -1
    atom_density = ZERO

    ! Determine score for each bin
    call score_general(p, t, n_nuclides_total*t % n_score_bins, filter_index, &
         i_nuclide, atom_density, flux)

  end subroutine score_all_nuclides

!===============================================================================
! SCORE_ANALOG_TALLY keeps track of how many events occur in a specified cell,
! energy range, etc. Note that since these are "analog" tallies, they are only
! triggered at every collision, not every event
!===============================================================================

  subroutine score_analog_tally_ce(p)

    type(Particle), intent(in) :: p

    integer :: i
    integer :: i_tally
    integer :: i_filt
    integer :: k                    ! loop index for nuclide bins
                                    ! position during the loop
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array
    real(8) :: filter_weight        ! combined weight of all filters
    type(TallyObject), pointer :: t

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do i = 1, active_analog_tallies % size()
      ! Get index of tally and pointer to tally
      i_tally = active_analog_tallies % get_item(i)
      t => tallies(i_tally)

      ! Find the first bin in each filter. There may be more than one matching
      ! bin per filter, but we'll deal with those later.
      do i_filt = 1, size(t % filters)
        call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
             NO_BIN_FOUND, matching_bins(i_filt), filter_weights(i_filt))
        ! If there are no valid bins for this filter, then there is nothing to
        ! score and we can move on to the next tally.
        if (matching_bins(i_filt) == NO_BIN_FOUND) cycle TALLY_LOOP
      end do

      ! ========================================================================
      ! Loop until we've covered all valid bins on each of the filters.

      FILTER_LOOP: do

        ! Determine scoring index and weight for this filter combination
        filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
             * t % stride) + 1
        filter_weight = product(filter_weights(:size(t % filters)))

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
          call score_general(p, t, (k-1)*t % n_score_bins, filter_index, &
               i_nuclide, ZERO, filter_weight)

        end do NUCLIDE_LOOP

        ! ======================================================================
        ! Filter logic

        ! If there are no filters, then we are done.
        if (size(t % filters) == 0) exit FILTER_LOOP

        ! Increment the filter bins, starting with the last filter. If we get a
        ! NO_BIN_FOUND for the last filter, it means we finished all valid bins
        ! for that filter, but next-to-last filter might have more than one
        ! valid bin so we need to increment that one as well, and so on.
        do i_filt = size(t % filters), 1, -1
          call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
               matching_bins(i_filt), matching_bins(i_filt), &
               filter_weights(i_filt))
          if (matching_bins(i_filt) /= NO_BIN_FOUND) exit
        end do

        ! If we got all NO_BIN_FOUNDs, then we have finished all valid bins for
        ! each of the filters. Exit the loop.
        if (all(matching_bins(:size(t % filters)) == NO_BIN_FOUND)) &
             exit FILTER_LOOP

        ! Reset all the filters with NO_BIN_FOUND. This will set them back to
        ! their first valid bin.
        do i_filt = 1, size(t % filters)
          if (matching_bins(i_filt) == NO_BIN_FOUND) then
            call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
                 matching_bins(i_filt), matching_bins(i_filt), &
                 filter_weights(i_filt))
          end if
        end do
      end do FILTER_LOOP

      ! If the user has specified that we can assume all tallies are spatially
      ! separate, this implies that once a tally has been scored to, we needn't
      ! check the others. This cuts down on overhead when there are many
      ! tallies specified

      if (assume_separate) exit TALLY_LOOP

    end do TALLY_LOOP

    ! Reset tally map positioning
    position = 0

  end subroutine score_analog_tally_ce

  subroutine score_analog_tally_mg(p)

    type(Particle), intent(in) :: p

    integer :: i, m
    integer :: i_tally
    integer :: i_filt
    integer :: k                    ! loop index for nuclide bins
                                    ! position during the loop
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array
    real(8) :: filter_weight        ! combined weight of all filters
    real(8) :: atom_density
    type(TallyObject), pointer :: t
    type(Material),    pointer :: mat

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do i = 1, active_analog_tallies % size()
      ! Get index of tally and pointer to tally
      i_tally = active_analog_tallies % get_item(i)
      t => tallies(i_tally)

      ! Get pointer to current material. We need this in order to determine what
      ! nuclides are in the material
      mat => materials(p % material)

      ! Find the first bin in each filter. There may be more than one matching
      ! bin per filter, but we'll deal with those later.
      do i_filt = 1, size(t % filters)
        call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
             NO_BIN_FOUND, matching_bins(i_filt), filter_weights(i_filt))
        ! If there are no valid bins for this filter, then there is nothing to
        ! score and we can move on to the next tally.
        if (matching_bins(i_filt) == NO_BIN_FOUND) cycle TALLY_LOOP
      end do

      ! ========================================================================
      ! Loop until we've covered all valid bins on each of the filters.

      FILTER_LOOP: do

        ! Determine scoring index and weight for this filter combination
        filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
             * t % stride) + 1
        filter_weight = product(filter_weights(:size(t % filters)))

        ! ======================================================================
        ! Nuclide logic

        ! Check for nuclide bins
        k = 0
        NUCLIDE_LOOP: do while (k < t % n_nuclide_bins)

          ! Increment the index in the list of nuclide bins
          k = k + 1

          i_nuclide = t % nuclide_bins(k)

          ! Check to see if this nuclide was in the material of our collision.
          do m = 1, mat % n_nuclides
            if (mat % nuclide(m) == i_nuclide) then
              atom_density = mat % atom_density(m)
              exit
            end if
          end do

          ! Determine score for each bin
          call score_general(p, t, (k-1)*t % n_score_bins, filter_index, &
               i_nuclide, atom_density, filter_weight)

        end do NUCLIDE_LOOP

        ! ======================================================================
        ! Filter logic

        ! If there are no filters, then we are done.
        if (size(t % filters) == 0) exit FILTER_LOOP

        ! Increment the filter bins, starting with the last filter. If we get a
        ! NO_BIN_FOUND for the last filter, it means we finished all valid bins
        ! for that filter, but next-to-last filter might have more than one
        ! valid bin so we need to increment that one as well, and so on.
        do i_filt = size(t % filters), 1, -1
          call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
               matching_bins(i_filt), matching_bins(i_filt), &
               filter_weights(i_filt))
          if (matching_bins(i_filt) /= NO_BIN_FOUND) exit
        end do

        ! If we got all NO_BIN_FOUNDs, then we have finished all valid bins for
        ! each of the filters. Exit the loop.
        if (all(matching_bins(:size(t % filters)) == NO_BIN_FOUND)) &
             exit FILTER_LOOP

        ! Reset all the filters with NO_BIN_FOUND. This will set them back to
        ! their first valid bin.
        do i_filt = 1, size(t % filters)
          if (matching_bins(i_filt) == NO_BIN_FOUND) then
            call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
                 matching_bins(i_filt), matching_bins(i_filt), &
                 filter_weights(i_filt))
          end if
        end do
      end do FILTER_LOOP

      ! If the user has specified that we can assume all tallies are spatially
      ! separate, this implies that once a tally has been scored to, we needn't
      ! check the others. This cuts down on overhead when there are many
      ! tallies specified

      if (assume_separate) exit TALLY_LOOP

    end do TALLY_LOOP

    ! Reset tally map positioning
    position = 0

  end subroutine score_analog_tally_mg

!===============================================================================
! SCORE_FISSION_EOUT handles a special case where we need to store neutron
! production rate with an outgoing energy filter (think of a fission matrix). In
! this case, we may need to score to multiple bins if there were multiple
! neutrons produced with different energies.
!===============================================================================

  subroutine score_fission_eout_ce(p, t, i_score, score_bin)

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
    integer :: bin_energyout ! original outgoing energy bin
    integer :: i_filter      ! index for matching filter bin combination
    real(8) :: filter_weight ! combined weight of all filters
    real(8) :: score         ! actual score
    real(8) :: E_out         ! energy of fission bank site

    ! save original outgoing energy bin and score index
    i = t % find_filter(FILTER_ENERGYOUT)
    bin_energyout = matching_bins(i)

    ! declare the energyout filter type
    select type(eo_filt => t % filters(i) % obj)
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

        ! determine outgoing energy from fission bank
        E_out = fission_bank(n_bank - p % n_bank + k) % E

        ! check if outgoing energy is within specified range on filter
        if (E_out < eo_filt % bins(1) .or. E_out > eo_filt % bins(n)) cycle

        ! change outgoing energy bin
        Matching_bins(i) = binary_search(eo_filt % bins, n, E_out)

        ! Case for tallying prompt neutrons
        if (score_bin == SCORE_NU_FISSION .or. &
             (score_bin == SCORE_PROMPT_NU_FISSION .and. g == 0)) then

          ! determine scoring index and weight for this filter combination
          i_filter = sum((matching_bins(1:size(t%filters)) - 1) * t % stride) &
               + 1
          filter_weight = product(filter_weights(:size(t % filters)))

          ! Add score to tally
!$omp atomic
          t % results(i_score, i_filter) % value = &
               t % results(i_score, i_filter) % value + score * filter_weight

        ! Case for tallying delayed emissions
        else if (score_bin == SCORE_DELAYED_NU_FISSION .and. g /= 0) then

          ! Get the index of delayed group filter
          j = t % find_filter(FILTER_DELAYEDGROUP)

          ! if the delayed group filter is present, tally to corresponding
          ! delayed group bin if it exists
          if (j > 0) then

            ! declare the delayed group filter type
            select type(dg_filt => t % filters(j) % obj)
            type is (DelayedGroupFilter)

              ! loop over delayed group bins until the corresponding bin is
              ! found
              do d_bin = 1, dg_filt % n_bins
                d = dg_filt % groups(d_bin)

                ! check whether the delayed group of the particle is equal to
                ! the delayed group of this bin
                if (d == g) then
                  call score_fission_delayed_dg(t, d_bin, score, i_score)
                end if
              end do
            end select

          ! if the delayed group filter is not present, add score to tally
          else

            ! determine scoring index and weight for this filter combination
            i_filter = sum((matching_bins(1:size(t%filters)) - 1) * t % stride)&
                 + 1
            filter_weight = product(filter_weights(:size(t % filters)))

            ! Add score to tally
!$omp atomic
            t % results(i_score, i_filter) % value = &
                 t % results(i_score, i_filter) % value + score * filter_weight
          end if
        end if
      end do
    end select

    ! reset outgoing energy bin and score index
    matching_bins(i) = bin_energyout

  end subroutine score_fission_eout_ce

  subroutine score_fission_eout_mg(p, t, i_score, i_nuclide, atom_density)
    type(Particle),    intent(in)    :: p
    type(TallyObject), intent(inout) :: t
    integer,           intent(in)    :: i_score   ! index for score
    integer,           intent(in)    :: i_nuclide ! index for nuclide
    real(8),           intent(in)    :: atom_density

    integer :: i             ! index of outgoing energy filter
    integer :: n             ! number of energies on filter
    integer :: k             ! loop index for bank sites
    integer :: bin_energyout ! original outgoing energy bin
    integer :: i_filter      ! index for matching filter bin combination
    real(8) :: filter_weight ! combined weight of all filters
    real(8) :: score         ! actual score
    integer :: gout          ! energy group of fission bank site
    integer :: gin           ! energy group of incident particle
    real(8) :: E_out

    ! save original outgoing energy bin and score index
    i = t % find_filter(FILTER_ENERGYOUT)
    bin_energyout = matching_bins(i)

    ! Declare the filter type
    select type(filt => t % filters(i) % obj)
    type is (EnergyoutFilter)

      ! Get number of energies on filter
      n = size(filt % bins)

      ! Since the creation of fission sites is weighted such that it is
      ! expected to create n_particles sites, we need to multiply the
      ! score by keff to get the true nu-fission rate. Otherwise, the sum
      ! of all nu-fission rates would be ~1.0.

      ! loop over number of particles banked
      do k = 1, p % n_bank
        ! determine score based on bank site weight and keff
        score = keff * fission_bank(n_bank - p % n_bank + k) % wgt
        if (i_nuclide > 0) then
          if (survival_biasing) then
            gin = p % g
          else
            gin = p % last_g
          end if
          score = score * atom_density * &
               nuclides_MG(i_nuclide) % obj % get_xs('fission', gin, &
                                                     UVW=p % last_uvw) / &
               macro_xs(p % material) % obj % get_xs('fission', gin, &
                                                     UVW=p % last_uvw)
        end if

        if (filt % matches_transport_groups) then
          ! determine outgoing energy from fission bank
          gout = int(fission_bank(n_bank - p % n_bank + k) % E)

          ! change outgoing energy bin
          matching_bins(i) = gout
        else
          ! determine outgoing energy from fission bank
          E_out = energy_bin_avg(int(fission_bank(n_bank - p % n_bank + k) % E))

          ! check if outgoing energy is within specified range on filter
          if (E_out < filt % bins(1) .or. E_out > filt % bins(n)) cycle

          ! change outgoing energy bin
          matching_bins(i) = binary_search(filt % bins, n, E_out)
        end if

        ! determine scoring index and weight for this filter combination
        i_filter = sum((matching_bins(1:size(t%filters)) - 1) * t % stride) + 1
        filter_weight = product(filter_weights(:size(t % filters)))

        ! Add score to tally
!$omp atomic
        t % results(i_score, i_filter) % value = &
             t % results(i_score, i_filter) % value + score * filter_weight
      end do

    ! reset outgoing energy bin and score index
    matching_bins(i) = bin_energyout
  end select

  end subroutine score_fission_eout_mg

!===============================================================================
! SCORE_FISSION_DELAYED_DG helper function used to increment the tally when a
! delayed group filter is present.
!===============================================================================

  subroutine score_fission_delayed_dg(t, d_bin, score, score_index)

    type(TallyObject), intent(inout) :: t
    integer, intent(in)              :: d_bin       ! delayed group bin index
    real(8), intent(in)              :: score       ! actual score
    integer, intent(in)              :: score_index ! index for score

    integer :: bin_original  ! original bin index
    integer :: filter_index  ! index for matching filter bin combination
    real(8) :: filter_weight ! combined weight of all filters

    ! save original delayed group bin
    bin_original = matching_bins(t % find_filter(FILTER_DELAYEDGROUP))
    matching_bins(t % find_filter(FILTER_DELAYEDGROUP)) = d_bin

    ! determine scoring index and weight on the modified matching_bins
    filter_index = sum((matching_bins(1:size(t % filters)) - 1) * t % stride) &
         + 1
    filter_weight = product(filter_weights(:size(t % filters)))

!$omp atomic
    t % results(score_index, filter_index) % value = &
         t % results(score_index, filter_index) % value + score * filter_weight

    ! reset original delayed group bin
    matching_bins(t % find_filter(FILTER_DELAYEDGROUP)) = bin_original

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
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array (from bins)
    real(8) :: flux                 ! tracklength estimate of flux
    real(8) :: atom_density         ! atom density of single nuclide in atom/b-cm
    real(8) :: filter_weight        ! combined weight of all filters
    type(TallyObject), pointer :: t
    type(Material),    pointer :: mat

    ! Determine track-length estimate of flux
    flux = p % wgt * distance

    ! A loop over all tallies is necessary because we need to simultaneously
    ! determine different filter bins for the same tally in order to score to it

    TALLY_LOOP: do i = 1, active_tracklength_tallies % size()
      ! Get index of tally and pointer to tally
      i_tally = active_tracklength_tallies % get_item(i)
      t => tallies(i_tally)

      ! Find the first bin in each filter. There may be more than one matching
      ! bin per filter, but we'll deal with those later.
      do i_filt = 1, size(t % filters)
        call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
             NO_BIN_FOUND, matching_bins(i_filt), filter_weights(i_filt))
        ! If there are no valid bins for this filter, then there is nothing to
        ! score and we can move on to the next tally.
        if (matching_bins(i_filt) == NO_BIN_FOUND) cycle TALLY_LOOP
      end do

      ! ========================================================================
      ! Loop until we've covered all valid bins on each of the filters.

      FILTER_LOOP: do

        ! Determine scoring index and weight for this filter combination
        filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
             * t % stride) + 1
        filter_weight = product(filter_weights(:size(t % filters)))

        ! ======================================================================
        ! Nuclide logic

        if (t % all_nuclides) then
          if (p % material /= MATERIAL_VOID) then
            call score_all_nuclides(p, i_tally, flux, filter_index)
          end if
        else

          NUCLIDE_BIN_LOOP: do k = 1, t % n_nuclide_bins
            ! Get index of nuclide in nuclides array
            i_nuclide = t % nuclide_bins(k)

            if (i_nuclide > 0) then
              if (p % material /= MATERIAL_VOID) then
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

        ! If there are no filters, then we are done.
        if (size(t % filters) == 0) exit FILTER_LOOP

        ! Increment the filter bins, starting with the last filter. If we get a
        ! NO_BIN_FOUND for the last filter, it means we finished all valid bins
        ! for that filter, but next-to-last filter might have more than one
        ! valid bin so we need to increment that one as well, and so on.
        do i_filt = size(t % filters), 1, -1
          call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
               matching_bins(i_filt), matching_bins(i_filt), &
               filter_weights(i_filt))
          if (matching_bins(i_filt) /= NO_BIN_FOUND) exit
        end do

        ! If we got all NO_BIN_FOUNDs, then we have finished all valid bins for
        ! each of the filters. Exit the loop.
        if (all(matching_bins(:size(t % filters)) == NO_BIN_FOUND)) &
             exit FILTER_LOOP

        ! Reset all the filters with NO_BIN_FOUND. This will set them back to
        ! their first valid bin.
        do i_filt = 1, size(t % filters)
          if (matching_bins(i_filt) == NO_BIN_FOUND) then
            call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
                 matching_bins(i_filt), matching_bins(i_filt), &
                 filter_weights(i_filt))
          end if
        end do
      end do FILTER_LOOP

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
    integer :: j                    ! loop index for scoring bins
    integer :: k                    ! loop index for nuclide bins
    integer :: filter_index         ! single index for single bin
    integer :: i_nuclide            ! index in nuclides array (from bins)
    real(8) :: flux                 ! collision estimate of flux
    real(8) :: atom_density         ! atom density of single nuclide
                                    !   in atom/b-cm
    real(8) :: filter_weight        ! combined weight of all filters
    type(TallyObject), pointer :: t
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
      i_tally = active_collision_tallies % get_item(i)
      t => tallies(i_tally)

      ! Find the first bin in each filter. There may be more than one matching
      ! bin per filter, but we'll deal with those later.
      do i_filt = 1, size(t % filters)
        call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
             NO_BIN_FOUND, matching_bins(i_filt), filter_weights(i_filt))
        ! If there are no valid bins for this filter, then there is nothing to
        ! score and we can move on to the next tally.
        if (matching_bins(i_filt) == NO_BIN_FOUND) cycle TALLY_LOOP
      end do

      ! ========================================================================
      ! Loop until we've covered all valid bins on each of the filters.

      FILTER_LOOP: do

        ! Determine scoring index and weight for this filter combination
        filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
             * t % stride) + 1
        filter_weight = product(filter_weights(:size(t % filters)))

        ! ======================================================================
        ! Nuclide logic

        if (t % all_nuclides) then
          if (p % material /= MATERIAL_VOID) then
            call score_all_nuclides(p, i_tally, flux, filter_index)
          end if
        else

          NUCLIDE_BIN_LOOP: do k = 1, t % n_nuclide_bins
            ! Get index of nuclide in nuclides array
            i_nuclide = t % nuclide_bins(k)

            if (i_nuclide > 0) then
              if (p % material /= MATERIAL_VOID) then
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

        ! If there are no filters, then we are done.
        if (size(t % filters) == 0) exit FILTER_LOOP

        ! Increment the filter bins, starting with the last filter. If we get a
        ! NO_BIN_FOUND for the last filter, it means we finished all valid bins
        ! for that filter, but next-to-last filter might have more than one
        ! valid bin so we need to increment that one as well, and so on.
        do i_filt = size(t % filters), 1, -1
          call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
               matching_bins(i_filt), matching_bins(i_filt), &
               filter_weights(i_filt))
          if (matching_bins(i_filt) /= NO_BIN_FOUND) exit
        end do

        ! If we got all NO_BIN_FOUNDs, then we have finished all valid bins for
        ! each of the filters. Exit the loop.
        if (all(matching_bins(:size(t % filters)) == NO_BIN_FOUND)) &
             exit FILTER_LOOP

        ! Reset all the filters with NO_BIN_FOUND. This will set them back to
        ! their first valid bin.
        do i_filt = 1, size(t % filters)
          if (matching_bins(i_filt) == NO_BIN_FOUND) then
            call t % filters(i_filt) % obj % get_next_bin(p, t % estimator, &
                 matching_bins(i_filt), matching_bins(i_filt), &
                 filter_weights(i_filt))
          end if
        end do
      end do FILTER_LOOP

      ! If the user has specified that we can assume all tallies are spatially
      ! separate, this implies that once a tally has been scored to, we needn't
      ! check the others. This cuts down on overhead when there are many
      ! tallies specified

      if (assume_separate) exit TALLY_LOOP

    end do TALLY_LOOP

    ! Reset tally map positioning
    position = 0

  end subroutine score_collision_tally

!===============================================================================
! SCORE_SURFACE_CURRENT tallies surface crossings in a mesh tally by manually
! determining which mesh surfaces were crossed
!===============================================================================

  subroutine score_surface_current(p)

    type(Particle), intent(in) :: p

    integer :: i
    integer :: i_tally
    integer :: j                    ! loop indices
    integer :: k                    ! loop indices
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
    real(8) :: filt_score           ! score applied by filters
    logical :: start_in_mesh        ! particle's starting xyz in mesh?
    logical :: end_in_mesh          ! particle's ending xyz in mesh?
    logical :: x_same               ! same starting/ending x index (i)
    logical :: y_same               ! same starting/ending y index (j)
    logical :: z_same               ! same starting/ending z index (k)
    type(TallyObject), pointer :: t
    type(RegularMesh), pointer :: m

    TALLY_LOOP: do i = 1, active_current_tallies % size()
      ! Copy starting and ending location of particle
      xyz0 = p % last_xyz_current
      xyz1 = p % coord(1) % xyz

      ! Get pointer to tally
      i_tally = active_current_tallies % get_item(i)
      t => tallies(i_tally)

      ! Get index for mesh, surface, and energy filters
      i_filter_mesh = t % find_filter(FILTER_MESH)
      i_filter_surf = t % find_filter(FILTER_SURFACE)
      i_filter_energy = t % find_filter(FILTER_ENERGYIN)

      ! Get pointer to mesh
      select type(filt => t % filters(i_filter_mesh) % obj)
      type is (MeshFilter)
        m => meshes(filt % mesh)
      end select

      ! Determine indices for starting and ending location
      call get_mesh_indices(m, xyz0, ijk0(:m % n_dimension), start_in_mesh)
      call get_mesh_indices(m, xyz1, ijk1(:m % n_dimension), end_in_mesh)

      ! Check to see if start or end is in mesh -- if not, check if track still
      ! intersects with mesh
      if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) then
        if (m % n_dimension == 2) then
          if (.not. mesh_intersects_2d(m, xyz0, xyz1)) cycle
        else
          if (.not. mesh_intersects_3d(m, xyz0, xyz1)) cycle
        end if
      end if

      ! Calculate number of surface crossings
      n_cross = sum(abs(ijk1 - ijk0))
      if (n_cross == 0) then
        cycle
      end if

      ! Copy particle's direction
      uvw = p % coord(1) % uvw

      ! Determine incoming energy bin.  We need to tell the energy filter this
      ! is a tracklength tally so it uses the pre-collision energy.
      if (i_filter_energy > 0) then
        call t % filters(i_filter_energy) % obj % get_next_bin(p, &
             ESTIMATOR_TRACKLENGTH, NO_BIN_FOUND, &
             matching_bins(i_filter_energy), filt_score)
        if (matching_bins(i_filter_energy) == NO_BIN_FOUND) cycle
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

            ! OUT_TOP
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_TOP
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_BOTTOM
            if (ijk0(1) >= 1 .and. ijk0(2) >= 1 .and. ijk0(3) >= 0 .and. &
                 ijk0(1) <= m % dimension(1) .and. &
                 ijk0(2) <= m % dimension(2) .and. &
                 ijk0(3) <  m % dimension(3)) then
              ijk0(3) = ijk0(3) + 1
              matching_bins(i_filter_surf) = IN_BOTTOM
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(3) = ijk0(3) - 1
            end if
          end do
        else
          do j = ijk0(3), ijk1(3) + 1, -1
            ijk0(3) = j

            ! OUT_BOTTOM
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_BOTTOM
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_TOP
            if (ijk0(1) >= 1 .and. ijk0(2) >= 1 .and. ijk0(3) > 1 .and. &
                 ijk0(1) <= m % dimension(1) .and. &
                 ijk0(2) <= m % dimension(2) .and. &
                 ijk0(3) <= m % dimension(3) + 1) then
              ijk0(3) = ijk0(3) - 1
              matching_bins(i_filter_surf) = IN_TOP
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(3) = ijk0(3) + 1
            end if
          end do
        end if
        cycle
      elseif (x_same .and. z_same) then
        ! Only y crossings
        if (uvw(2) > 0) then
          do j = ijk0(2), ijk1(2) - 1
            ijk0(2) = j

            ! OUT_FRONT
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_FRONT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_BACK
            if (ijk0(1) >= 1 .and. ijk0(2) >= 0 .and. ijk0(3) >= 1 .and. &
                 ijk0(1) <= m % dimension(1) .and. &
                 ijk0(2) <  m % dimension(2) .and. &
                 ijk0(3) <= m % dimension(3)) then
              ijk0(2) = ijk0(2) + 1
              matching_bins(i_filter_surf) = IN_BACK
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(2) = ijk0(2) - 1
            end if
          end do
        else
          do j = ijk0(2), ijk1(2) + 1, -1
            ijk0(2) = j

            ! OUT_BACK
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_BACK
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_FRONT
            if (ijk0(1) >= 1 .and. ijk0(2) > 1 .and. ijk0(3) >= 1 .and. &
                 ijk0(1) <= m % dimension(1) .and. &
                 ijk0(2) <= m % dimension(2) + 1 .and. &
                 ijk0(3) <= m % dimension(3)) then
              ijk0(2) = ijk0(2) - 1
              matching_bins(i_filter_surf) = IN_FRONT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(2) = ijk0(2) + 1
            end if
          end do
        end if
        cycle
      elseif (y_same .and. z_same) then
        ! Only x crossings
        if (uvw(1) > 0) then
          do j = ijk0(1), ijk1(1) - 1
            ijk0(1) = j

            ! OUT_RIGHT
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_RIGHT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_LEFT
            if (ijk0(1) >= 0 .and. ijk0(2) >= 1 .and. ijk0(3) >= 1 .and. &
                 ijk0(1) <  m % dimension(1) .and. &
                 ijk0(2) <= m % dimension(2) .and. &
                 ijk0(3) <= m % dimension(3)) then
              ijk0(1) = ijk0(1) + 1
              matching_bins(i_filter_surf) = IN_LEFT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(1) = ijk0(1) - 1
            end if
          end do
        else
          do j = ijk0(1), ijk1(1) + 1, -1
            ijk0(1) = j

            ! OUT_LEFT
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_LEFT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_RIGHT
            if (ijk0(1) > 1 .and. ijk0(2) >= 1 .and. ijk0(3) >= 1 .and. &
                 ijk0(1) <= m % dimension(1) + 1 .and. &
                 ijk0(2) <= m % dimension(2) .and. &
                 ijk0(3) <= m % dimension(3)) then
              ijk0(1) = ijk0(1) - 1
              matching_bins(i_filter_surf) = IN_RIGHT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(1) = ijk0(1) + 1
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
        matching_bins(i_filter_surf) = 0

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

        ! Now use the minimum distance and direction of the particle to
        ! determine which surface was crossed

        if (distance == d(1)) then
          if (uvw(1) > 0) then

            ! OUT_RIGHT
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_RIGHT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_LEFT
            if (ijk0(1) >= 0 .and. ijk0(2) >= 1 .and. ijk0(3) >= 1 .and. &
                 ijk0(1) <  m % dimension(1) .and. &
                 ijk0(2) <= m % dimension(2) .and. &
                 ijk0(3) <= m % dimension(3)) then
              ijk0(1) = ijk0(1) + 1
              matching_bins(i_filter_surf) = IN_LEFT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(1) = ijk0(1) - 1
            end if

            ijk0(1) = ijk0(1) + 1
            xyz_cross(1) = xyz_cross(1) + m % width(1)
          else

            ! OUT_LEFT
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_LEFT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_RIGHT
            if (ijk0(1) > 1 .and. ijk0(2) >= 1 .and. ijk0(3) >= 1 .and. &
                 ijk0(1) <= m % dimension(1) + 1 .and. &
                 ijk0(2) <= m % dimension(2) .and. &
                 ijk0(3) <= m % dimension(3)) then
              ijk0(1) = ijk0(1) - 1
              matching_bins(i_filter_surf) = IN_RIGHT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(1) = ijk0(1) + 1
            end if

            ijk0(1) = ijk0(1) - 1
            xyz_cross(1) = xyz_cross(1) - m % width(1)
          end if
        elseif (distance == d(2)) then
          if (uvw(2) > 0) then

            ! OUT_FRONT
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_FRONT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_BACK
            if (ijk0(1) >= 1 .and. ijk0(2) >= 0 .and. ijk0(3) >= 1 .and. &
                 ijk0(1) <= m % dimension(1) .and. &
                 ijk0(2) <  m % dimension(2) .and. &
                 ijk0(3) <= m % dimension(3)) then
              ijk0(2) = ijk0(2) + 1
              matching_bins(i_filter_surf) = IN_BACK
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(2) = ijk0(2) - 1
            end if

            ijk0(2) = ijk0(2) + 1
            xyz_cross(2) = xyz_cross(2) + m % width(2)
          else

            ! OUT_BACK
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_BACK
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_FRONT
            if (ijk0(1) >= 1 .and. ijk0(2) > 1 .and. ijk0(3) >= 1 .and. &
                 ijk0(1) <= m % dimension(1) .and. &
                 ijk0(2) <= m % dimension(2) + 1 .and. &
                 ijk0(3) <= m % dimension(3)) then
              ijk0(2) = ijk0(2) - 1
              matching_bins(i_filter_surf) = IN_FRONT
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(2) = ijk0(2) + 1
            end if

            ijk0(2) = ijk0(2) - 1
            xyz_cross(2) = xyz_cross(2) - m % width(2)
          end if
        else if (distance == d(3)) then
          if (uvw(3) > 0) then

            ! OUT_TOP
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_TOP
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_BOTTOM
            if (ijk0(1) >= 1 .and. ijk0(2) >= 1 .and. ijk0(3) >= 0 .and. &
                 ijk0(1) <= m % dimension(1) .and. &
                 ijk0(2) <= m % dimension(2) .and. &
                 ijk0(3) <  m % dimension(3)) then
              ijk0(3) = ijk0(3) + 1
              matching_bins(i_filter_surf) = IN_BOTTOM
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(3) = ijk0(3) - 1
            end if

            ijk0(3) = ijk0(3) + 1
            xyz_cross(3) = xyz_cross(3) + m % width(3)
          else

            ! OUT_BOTTOM
            if (all(ijk0 >= 1) .and. all(ijk0 <= m % dimension)) then
              matching_bins(i_filter_surf) = OUT_BOTTOM
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
            end if

            ! IN_TOP
            if (ijk0(1) >= 1 .and. ijk0(2) >= 1 .and. ijk0(3) > 1 .and. &
                 ijk0(1) <= m % dimension(1) .and. &
                 ijk0(2) <= m % dimension(2) .and. &
                 ijk0(3) <= m % dimension(3) + 1) then
              ijk0(3) = ijk0(3) - 1
              matching_bins(i_filter_surf) = IN_TOP
              matching_bins(i_filter_mesh) = &
                   mesh_indices_to_bin(m, ijk0)
              filter_index = sum((matching_bins(1:size(t % filters)) - 1) &
                   * t % stride) + 1
!$omp atomic
              t % results(1, filter_index) % value = &
                   t % results(1, filter_index) % value + p % wgt
              ijk0(3) = ijk0(3) + 1
            end if

            ijk0(3) = ijk0(3) - 1
            xyz_cross(3) = xyz_cross(3) - m % width(3)
          end if
        end if

        ! Calculate new coordinates
        xyz0 = xyz0 + distance * uvw
      end do

    end do TALLY_LOOP

  end subroutine score_surface_current

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
    type(TallyObject), pointer :: t

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
    type(TallyObject), pointer :: t

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
        elseif (user_tallies(i) % estimator == ESTIMATOR_COLLISION) then
          call active_collision_tallies % add(i_user_tallies + i)
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
      call fatal_error("Active tallies should not exist before CMFD tallies!")
    else if (active_analog_tallies % size() > 0) then
      call fatal_error('Active analog tallies should not exist before CMFD &
           &tallies!')
    else if (active_tracklength_tallies % size() > 0) then
      call fatal_error("Active tracklength tallies should not exist before &
           &CMFD tallies!")
    else if (active_current_tallies % size() > 0) then
      call fatal_error("Active current tallies should not exist before CMFD &
           &tallies!")
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
