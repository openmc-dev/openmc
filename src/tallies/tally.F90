module tally

  use, intrinsic :: ISO_C_BINDING

  use algorithm,        only: binary_search
  use bank_header
  use constants
  use dict_header,      only: EMPTY
  use error,            only: fatal_error
  use geometry_header
  use material_header
  use math,             only: t_percentile
  use message_passing
  use mgxs_interface
  use nuclide_header
  use output,           only: header
  use particle_header,  only: LocalCoord, Particle
  use settings
  use simulation_header
  use string,           only: to_str
  use tally_derivative_header
  use tally_filter
  use tally_header

  implicit none

  procedure(score_general_),      pointer :: score_general => null()
  procedure(score_analog_tally_), pointer :: score_analog_tally => null()

  abstract interface
    subroutine score_general_(p, i_tally, start_index, filter_index, i_nuclide, &
                              atom_density, flux)
      import Particle, C_INT, C_DOUBLE
      type(Particle), intent(in)        :: p
      integer(C_INT), intent(in), value :: i_tally
      integer(C_INT), intent(in), value :: start_index
      integer(C_INT), intent(in), value :: i_nuclide
      integer(C_INT), intent(in), value :: filter_index   ! for % results
      real(C_DOUBLE), intent(in), value :: flux           ! flux estimate
      real(C_DOUBLE), intent(in), value :: atom_density   ! atom/b-cm
    end subroutine score_general_

    subroutine score_analog_tally_(p)
      import Particle
      type(Particle), intent(in) :: p
    end subroutine score_analog_tally_
  end interface

  interface
    subroutine score_general_ce(p, i_tally, start_index, filter_index, &
         i_nuclide, atom_density, flux) bind(C)
      import Particle, C_INT, C_DOUBLE
      type(Particle), intent(in) :: p
      integer(C_INT), intent(in), value :: i_tally
      integer(C_INT), intent(in), value :: start_index
      integer(C_INT), intent(in), value :: filter_index
      integer(C_INT), intent(in), value :: i_nuclide
      real(C_DOUBLE), intent(in), value :: atom_density
      real(C_DOUBLE), intent(in), value :: flux
    end subroutine

    subroutine score_general_mg_c(p, i_tally, start_index, filter_index, &
         i_nuclide, atom_density, flux) bind(C)
      import Particle, C_INT, C_DOUBLE
      type(Particle), intent(in) :: p
      integer(C_INT), intent(in), value :: i_tally
      integer(C_INT), intent(in), value :: start_index
      integer(C_INT), intent(in), value :: filter_index
      integer(C_INT), intent(in), value :: i_nuclide
      real(C_DOUBLE), intent(in), value :: atom_density
      real(C_DOUBLE), intent(in), value :: flux
    end subroutine

    subroutine score_fission_delayed_dg(i_tally, d_bin, score, score_index) bind(C)
      import C_INT, C_DOUBLE
      integer(C_INT), value :: i_tally
      integer(C_INT), value :: d_bin
      real(C_DOUBLE), value :: score
      integer(C_INT), value :: score_index
    end subroutine

    subroutine score_fission_eout(p, i_tally, i_score, score_bin) bind(C)
      import Particle, C_INT
      type(Particle)         :: p
      integer(C_INT), value :: i_tally
      integer(C_INT), value :: i_score
      integer(C_INT), value :: score_bin
    end subroutine

    subroutine score_analog_tally_ce(p) bind(C)
      import Particle
      type(Particle), intent(in) :: p
    end subroutine

    subroutine score_analog_tally_mg(p) bind(C)
      import Particle
      type(Particle), intent(in) :: p
    end subroutine

    subroutine score_tracklength_tally(p, distance) bind(C)
      import Particle, C_DOUBLE
      type(Particle) :: p
      real(C_DOUBLE), value :: distance
    end subroutine

    subroutine score_collision_tally(p) bind(C)
      import Particle
      type(Particle) :: p
    end subroutine

    subroutine score_meshsurface_tally(p) bind(C)
      import Particle
      type(Particle) :: p
    end subroutine

    subroutine score_surface_tally(p) bind(C)
      import Particle
      type(Particle) :: p
    end subroutine

    subroutine zero_flux_derivs() bind(C)
    end subroutine
  end interface

contains

!===============================================================================
! INIT_TALLY_ROUTINES Sets the procedure pointers needed for minimizing code
! with the CE and MG modes.
!===============================================================================

  subroutine init_tally_routines() bind(C)
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

  subroutine score_general_mg(p, i_tally, start_index, filter_index, i_nuclide, &
       atom_density, flux) bind(C)
    type(Particle), intent(in)        :: p
    integer(C_INT), intent(in), value :: i_tally
    integer(C_INT), intent(in), value :: start_index
    integer(C_INT), intent(in), value :: i_nuclide
    integer(C_INT), intent(in), value :: filter_index   ! for % results
    real(C_DOUBLE), intent(in), value :: flux           ! flux estimate
    real(C_DOUBLE), intent(in), value :: atom_density   ! atom/b-cm

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

    associate (t => tallies(i_tally) % obj)

    ! Set the direction and group to use with get_xs
    if (t % estimator() == ESTIMATOR_ANALOG .or. &
         t % estimator() == ESTIMATOR_COLLISION) then

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
    call set_macro_angle_index_c(p % material, p_uvw)

    ! Do same for nucxs, point it to the microscopic nuclide data of interest
    if (i_nuclide >= 0) then
      ! And since we haven't calculated this temperature index yet, do so now
      call set_nuclide_temperature_index_c(i_nuclide+1, p % sqrtkT)
      call set_nuclide_angle_index_c(i_nuclide+1, p_uvw)
    end if

    call score_general_mg_c(p, i_tally, start_index, filter_index, i_nuclide, &
         atom_density, flux)

    i = 0
    SCORE_LOOP: do q = 1, t % n_score_bins()
      i = i + 1

      ! determine what type of score bin
      score_bin = t % score_bins(i)

      ! determine scoring bin index
      score_index = start_index + i

      !#########################################################################
      ! Determine appropirate scoring value.

      select case(score_bin)


      case (SCORE_FLUX)
        cycle SCORE_LOOP


      case (SCORE_TOTAL)
        cycle SCORE_LOOP


      case (SCORE_INVERSE_VELOCITY)
        cycle SCORE_LOOP


      case (SCORE_SCATTER)
        cycle SCORE_LOOP


      case (SCORE_NU_SCATTER)
        cycle SCORE_LOOP


      case (SCORE_ABSORPTION)
        cycle SCORE_LOOP


      case (SCORE_FISSION)
        cycle SCORE_LOOP


      case (SCORE_NU_FISSION)
        cycle SCORE_LOOP


      case (SCORE_PROMPT_NU_FISSION)

        if (t % estimator() == ESTIMATOR_ANALOG) then
          if (survival_biasing .or. p % fission) then
            if (t % energyout_filter() > 0) then
              ! Normally, we only need to make contributions to one scoring
              ! bin. However, in the case of fission, since multiple fission
              ! neutrons were emitted with different energies, multiple
              ! outgoing energy bins may have been scored to. The following
              ! logic treats this special case and results to multiple bins
              call score_fission_eout(p, i_tally, score_index, score_bin)
              cycle SCORE_LOOP
            end if
          end if
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! nu-fission
            score = p % absorb_wgt * flux
            if (i_nuclide >= 0) then
              score = score * atom_density * &
                   get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_PROMPT_NU_FISSION, p_g) / &
                   get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
            else
              score = score * &
                   get_macro_xs_c(p % material, MG_GET_XS_PROMPT_NU_FISSION, p_g) / &
                   get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
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
            if (i_nuclide >= 0) then
              score = score * atom_density * &
                   get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_FISSION, p_g) / &
                   get_macro_xs_c(p % material, MG_GET_XS_FISSION, p_g)
            end if
          end if

        else
          if (i_nuclide >= 0) then
            score = get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_PROMPT_NU_FISSION, p_g) * &
                 atom_density * flux
          else
            score = get_macro_xs_c(p % material, MG_GET_XS_PROMPT_NU_FISSION, p_g) * flux
          end if
        end if


      case (SCORE_DELAYED_NU_FISSION)

        ! Set the delayedgroup filter index and the number of delayed group bins
        dg_filter = t % delayedgroup_filter()

        if (t % estimator() == ESTIMATOR_ANALOG) then
          if (survival_biasing .or. p % fission) then
            if (t % energyout_filter() > 0) then
              ! Normally, we only need to make contributions to one scoring
              ! bin. However, in the case of fission, since multiple fission
              ! neutrons were emitted with different energies, multiple
              ! outgoing energy bins may have been scored to. The following
              ! logic treats this special case and results to multiple bins
              call score_fission_eout(p, i_tally, score_index, score_bin)
              cycle SCORE_LOOP
            end if
          end if
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! nu-fission
            if (get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g) > ZERO) then

              if (dg_filter > 0) then
                select type(filt => filters(t % filter(dg_filter) + 1) % obj)
                type is (DelayedGroupFilter)

                  ! Loop over all delayed group bins and tally to them
                  ! individually
                  do d_bin = 1, filt % n_bins

                    ! Get the delayed group for this bin
                    d = filt % groups(d_bin)

                    score = p % absorb_wgt * flux
                    if (i_nuclide >= 0) then
                      score = score * &
                           get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d) / &
                           get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
                    else
                      score = score * &
                           get_macro_xs_c(p % material, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d) / &
                           get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
                    end if

                    call score_fission_delayed_dg(i_tally, d_bin, score, score_index)
                  end do
                  cycle SCORE_LOOP
                end select
              else
                score = p % absorb_wgt * flux
                if (i_nuclide >= 0) then
                  score = score * &
                       get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DELAYED_NU_FISSION, p_g) / &
                       get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
                else
                  score = score * &
                       get_macro_xs_c(p % material, MG_GET_XS_DELAYED_NU_FISSION, p_g) / &
                       get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
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
              select type(filt => filters(t % filter(dg_filter) + 1) % obj)
              type is (DelayedGroupFilter)

                ! Loop over all delayed group bins and tally to them
                ! individually
                do d_bin = 1, filt % n_bins

                  ! Get the delayed group for this bin
                  d = filt % groups(d_bin)

                  score = keff * p % wgt_bank / p % n_bank * &
                       p % n_delayed_bank(d) * flux

                  if (i_nuclide >= 0) then
                    score = score * atom_density * &
                         get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_FISSION, p_g) / &
                         get_macro_xs_c(p % material, MG_GET_XS_FISSION, p_g)
                  end if

                  call score_fission_delayed_dg(i_tally, d_bin, score, score_index)
                end do
                cycle SCORE_LOOP
              end select
            else
              score = keff * p % wgt_bank / p % n_bank * sum(p % n_delayed_bank) * flux
              if (i_nuclide >= 0) then
                score = score * atom_density * &
                     get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_FISSION, p_g) / &
                     get_macro_xs_c(p % material, MG_GET_XS_FISSION, p_g)
              end if
            end if
          end if
        else

          ! Check if the delayed group filter is present
          if (dg_filter > 0) then
            select type(filt => filters(t % filter(dg_filter) + 1) % obj)
            type is (DelayedGroupFilter)

              ! Loop over all delayed group bins and tally to them
              ! individually
              do d_bin = 1, filt % n_bins

                ! Get the delayed group for this bin
                d = filt % groups(d_bin)

                if (i_nuclide >= 0) then
                  score = atom_density * flux * &
                       get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d)
                else
                  score = flux * &
                       get_macro_xs_c(p % material, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d)
                end if

                call score_fission_delayed_dg(i_tally, d_bin, score, score_index)
              end do
              cycle SCORE_LOOP
            end select
          else
            if (i_nuclide >= 0) then
              score = atom_density * flux * &
                   get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DELAYED_NU_FISSION, p_g)

            else
              score = flux * &
                   get_macro_xs_c(p % material, MG_GET_XS_DELAYED_NU_FISSION, p_g)
            end if
          end if
        end if

      case (SCORE_DECAY_RATE)

        ! Set the delayedgroup filter index and the number of delayed group bins
        dg_filter = t % delayedgroup_filter()

        if (t % estimator() == ESTIMATOR_ANALOG) then
          if (survival_biasing) then
            ! No fission events occur if survival biasing is on -- need to
            ! calculate fraction of absorptions that would have resulted in
            ! nu-fission
            if (get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g) > ZERO) then

              if (dg_filter > 0) then
                select type(filt => filters(t % filter(dg_filter) + 1) % obj)
                type is (DelayedGroupFilter)

                  ! Loop over all delayed group bins and tally to them
                  ! individually
                  do d_bin = 1, filt % n_bins

                    ! Get the delayed group for this bin
                    d = filt % groups(d_bin)

                    score = p % absorb_wgt * flux
                    if (i_nuclide >= 0) then
                      score = score * &
                           get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DECAY_RATE, p_g, DG=d) * &
                           get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d) / &
                           get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
                    else
                      score = score * &
                           get_macro_xs_c(p % material, MG_GET_XS_DECAY_RATE, p_g, DG=d) * &
                           get_macro_xs_c(p % material, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d) / &
                           get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
                    end if

                    call score_fission_delayed_dg(i_tally, d_bin, score, score_index)
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
                  if (i_nuclide >= 0) then
                    score = score + p % absorb_wgt * flux * &
                         get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DECAY_RATE, p_g, DG=d) * &
                         get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d) / &
                         get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
                  else
                    score = score + p % absorb_wgt * flux * &
                         get_macro_xs_c(p % material, MG_GET_XS_DECAY_RATE, p_g, DG=d) * &
                         get_macro_xs_c(p % material, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d) / &
                         get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
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
              g = fission_bank_delayed_group(n_bank - p % n_bank + k)

              ! Case for tallying delayed emissions
              if (g /= 0) then

                ! determine score based on bank site weight and keff.
                if (i_nuclide >= 0) then
                  score = score + keff * atom_density * &
                       fission_bank_wgt(n_bank - p % n_bank + k) * &
                       get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DECAY_RATE, p_g, DG=g) * &
                       get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_FISSION, p_g) / &
                       get_macro_xs_c(p % material, MG_GET_XS_FISSION, p_g) * flux
                else
                  score = score + keff * &
                       fission_bank_wgt(n_bank - p % n_bank + k) * &
                       get_macro_xs_c(p % material, MG_GET_XS_DECAY_RATE, p_g, DG=g) * flux
                end if

                ! if the delayed group filter is present, tally to corresponding
                ! delayed group bin if it exists
                if (dg_filter > 0) then

                  ! declare the delayed group filter type
                  select type(filt => filters(t % filter(dg_filter) + 1) % obj)
                  type is (DelayedGroupFilter)

                    ! loop over delayed group bins until the corresponding bin
                    ! is found
                    do d_bin = 1, filt % n_bins
                      d = filt % groups(d_bin)

                      ! check whether the delayed group of the particle is equal
                      ! to the delayed group of this bin
                      if (d == g) then
                        call score_fission_delayed_dg(i_tally, d_bin, score, &
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
            select type(filt => filters(t % filter(dg_filter) + 1) % obj)
            type is (DelayedGroupFilter)

              ! Loop over all delayed group bins and tally to them
              ! individually
              do d_bin = 1, filt % n_bins

                ! Get the delayed group for this bin
                d = filt % groups(d_bin)

                if (i_nuclide >= 0) then
                  score = atom_density * flux * &
                       get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DECAY_RATE, p_g, DG=d) * &
                       get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d)
                else
                  score = flux * &
                       get_macro_xs_c(p % material, MG_GET_XS_DECAY_RATE, p_g, DG=d) * &
                       get_macro_xs_c(p % material, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d)
                end if

                call score_fission_delayed_dg(i_tally, d_bin, score, score_index)
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
              if (i_nuclide >= 0) then
                score = score + atom_density * flux * &
                     get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DECAY_RATE, p_g, DG=d) * &
                     get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d)
              else
                score = score + flux * &
                     get_macro_xs_c(p % material, MG_GET_XS_DECAY_RATE, p_g, DG=d) * &
                     get_macro_xs_c(p % material, MG_GET_XS_DELAYED_NU_FISSION, p_g, DG=d)
              end if
            end do
          end if
        end if


      case (SCORE_KAPPA_FISSION)

        if (t % estimator() == ESTIMATOR_ANALOG) then
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
          if (i_nuclide >= 0) then
            score = score * atom_density * &
                 get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_KAPPA_FISSION, p_g) / &
                 get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
          else
            score = score * &
                 get_macro_xs_c(p % material, MG_GET_XS_KAPPA_FISSION, p_g) / &
                 get_macro_xs_c(p % material, MG_GET_XS_ABSORPTION, p_g)
          end if
        else
          if (i_nuclide >= 0) then
            score = get_nuclide_xs_c(i_nuclide+1, MG_GET_XS_KAPPA_FISSION, p_g) * &
                 atom_density * flux
          else
            score = flux * &
                 get_macro_xs_c(p % material, MG_GET_XS_KAPPA_FISSION, p_g)

          end if
        end if


      case (SCORE_EVENTS)
        ! Simply count number of scoring events
        score = ONE

      end select

      !#########################################################################
      ! Expand score if necessary and add to tally results.
!$omp atomic
      t % results(RESULT_VALUE, score_index, filter_index) = &
           t % results(RESULT_VALUE, score_index, filter_index) + score

    end do SCORE_LOOP
    end associate
  end subroutine score_general_mg

!===============================================================================
! APPLY_DERIVATIVE_TO_SCORE multiply the given score by its relative derivative
!===============================================================================

  subroutine apply_derivative_to_score(p, i_tally, i_nuclide, atom_density, &
                                       score_bin, score) bind(C)
    type(Particle),        intent(in)    :: p
    integer(C_INT), value, intent(in)    :: i_tally
    integer(C_INT), value, intent(in)    :: i_nuclide
    real(C_DOUBLE), value, intent(in)    :: atom_density   ! atom/b-cm
    integer(C_INT), value, intent(in)    :: score_bin
    real(C_DOUBLE),        intent(inout) :: score

    type(TallyDerivative), pointer :: deriv
    integer :: l
    integer :: i_nuc
    logical :: scoring_diff_nuclide
    real(8) :: flux_deriv
    real(8) :: dsig_s, dsig_a, dsig_f, cum_dsig

    associate (t => tallies(i_tally) % obj)

    if (score == ZERO) return

    ! If our score was previously c then the new score is
    ! c * (1/f * d_f/d_p + 1/c * d_c/d_p)
    ! where (1/f * d_f/d_p) is the (logarithmic) flux derivative and p is the
    ! perturbated variable.

    !associate(deriv => tally_derivs(t % deriv()))
    deriv => tally_deriv_c(t % deriv())
      flux_deriv = deriv % flux_deriv

      !select case (tally_derivs(t % deriv()) % variable)
      select case (deriv % variable)

      !=========================================================================
      ! Density derivative:
      ! c = Sigma_MT
      ! c = sigma_MT * N
      ! c = sigma_MT * rho * const
      ! d_c / d_rho = sigma_MT * const
      ! (1 / c) * (d_c / d_rho) = 1 / rho

      case (DIFF_DENSITY)
        select case (t % estimator())

        case (ESTIMATOR_ANALOG)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL, SCORE_SCATTER, SCORE_ABSORPTION, SCORE_FISSION, &
                SCORE_NU_FISSION)
            if (material_id(p % material) == deriv % diff_material) then
              score = score * (flux_deriv + ONE &
                   / material_density_gpcc(p % material))
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id())))
          end select

        case (ESTIMATOR_COLLISION)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL, SCORE_SCATTER, SCORE_ABSORPTION, SCORE_FISSION, &
                SCORE_NU_FISSION)
            if (material_id(p % material) == deriv % diff_material) then
              score = score * (flux_deriv + ONE &
                   / material_density_gpcc(p % material))
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id())))
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
        select case (t % estimator())

        case (ESTIMATOR_ANALOG)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL, SCORE_SCATTER, SCORE_ABSORPTION, SCORE_FISSION, &
                SCORE_NU_FISSION)
            if (material_id(p % material) == deriv % diff_material &
                 .and. p % event_nuclide == deriv % diff_nuclide) then
                ! Search for the index of the perturbed nuclide.
                do l = 1, material_nuclide_size(p % material)
                  if (material_nuclide(p % material, l) == deriv % diff_nuclide) exit
                end do

                score = score * (flux_deriv &
                     + ONE / material_atom_density(p % material, l))
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id())))
          end select

        case (ESTIMATOR_COLLISION)
          scoring_diff_nuclide = &
               (material_id(p % material) == deriv % diff_material) &
               .and. (i_nuclide+1 == deriv % diff_nuclide)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL)
            if (i_nuclide == -1 .and. &
                 material_id(p % material) == deriv % diff_material .and. &
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
                 material_id(p % material) == deriv % diff_material .and. &
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
                 material_id(p % material) == deriv % diff_material .and. &
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
                 material_id(p % material) == deriv % diff_material .and. &
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
                 material_id(p % material) == deriv % diff_material .and. &
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
                 &tally ' // trim(to_str(t % id())))
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
        select case (t % estimator())

        case (ESTIMATOR_ANALOG)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL)
            if (material_id(p % material) == deriv % diff_material .and. &
                 micro_xs(p % event_nuclide) % total > ZERO) then
                ! Search for the index of the perturbed nuclide.
                do l = 1, material_nuclide_size(p % material)
                  if (material_nuclide(p % material, l) == p % event_nuclide) exit
                end do

                dsig_s = ZERO
                dsig_a = ZERO
                associate (nuc => nuclides(p % event_nuclide))
                  if (multipole_in_range(nuc % ptr, p % last_E)) then
                    call multipole_deriv_eval(nuc % ptr, p % last_E, &
                         p % sqrtkT, dsig_s, dsig_a, dsig_f)
                  end if
                end associate
                score = score * (flux_deriv &
                     + (dsig_s + dsig_a) * material_atom_density(p % material, l) &
                     / material_xs % total)
            else
              score = score * flux_deriv
            end if

          case (SCORE_SCATTER)
            if (material_id(p % material) == deriv % diff_material .and. &
                 (micro_xs(p % event_nuclide) % total &
                 - micro_xs(p % event_nuclide) % absorption) > ZERO) then
                ! Search for the index of the perturbed nuclide.
                do l = 1, material_nuclide_size(p % material)
                  if (material_nuclide(p % material, l) == p % event_nuclide) exit
                end do

                dsig_s = ZERO
                associate (nuc => nuclides(p % event_nuclide))
                  if (multipole_in_range(nuc % ptr, p % last_E)) then
                    call multipole_deriv_eval(nuc % ptr, p % last_E, &
                         p % sqrtkT, dsig_s, dsig_a, dsig_f)
                  end if
                end associate
                score = score * (flux_deriv + dsig_s * material_atom_density(p % material, l) / &
                     (material_xs % total - material_xs % absorption))
            else
              score = score * flux_deriv
            end if

          case (SCORE_ABSORPTION)
            if (material_id(p % material) == deriv % diff_material .and. &
                 micro_xs(p % event_nuclide) % absorption > ZERO) then
                ! Search for the index of the perturbed nuclide.
                do l = 1, material_nuclide_size(p % material)
                  if (material_nuclide(p % material, l) == p % event_nuclide) exit
                end do

                dsig_a = ZERO
                associate (nuc => nuclides(p % event_nuclide))
                  if (multipole_in_range(nuc % ptr, p % last_E)) then
                    call multipole_deriv_eval(nuc % ptr, p % last_E, &
                         p % sqrtkT, dsig_s, dsig_a, dsig_f)
                  end if
                end associate
                score = score * (flux_deriv + dsig_a * material_atom_density(p % material, l) &
                                              / material_xs % absorption)
            else
              score = score * flux_deriv
            end if

          case (SCORE_FISSION)
            if (material_id(p % material) == deriv % diff_material .and. &
                 micro_xs(p % event_nuclide) % fission > ZERO) then
                ! Search for the index of the perturbed nuclide.
                do l = 1, material_nuclide_size(p % material)
                  if (material_nuclide(p % material, l) == p % event_nuclide) exit
                end do

                dsig_f = ZERO
                associate (nuc => nuclides(p % event_nuclide))
                  if (multipole_in_range(nuc % ptr, p % last_E)) then
                    call multipole_deriv_eval(nuc % ptr, p % last_E, &
                         p % sqrtkT, dsig_s, dsig_a, dsig_f)
                  end if
                end associate
                score = score * (flux_deriv &
                     + dsig_f * material_atom_density(p % material, l) / material_xs % fission)
            else
              score = score * flux_deriv
            end if

          case (SCORE_NU_FISSION)
            if (material_id(p % material) == deriv % diff_material .and. &
                 micro_xs(p % event_nuclide) % nu_fission > ZERO) then
                ! Search for the index of the perturbed nuclide.
                do l = 1, material_nuclide_size(p % material)
                  if (material_nuclide(p % material, l) == p % event_nuclide) exit
                end do

                dsig_f = ZERO
                associate (nuc => nuclides(p % event_nuclide))
                  if (multipole_in_range(nuc % ptr, p % last_E)) then
                    call multipole_deriv_eval(nuc % ptr, p % last_E, &
                         p % sqrtkT, dsig_s, dsig_a, dsig_f)
                  end if
                end associate
                score = score * (flux_deriv &
                     + dsig_f * material_atom_density(p % material, l) / material_xs % nu_fission&
                     * micro_xs(p % event_nuclide) % nu_fission &
                     / micro_xs(p % event_nuclide) % fission)
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id())))
          end select

        case (ESTIMATOR_COLLISION)

          select case (score_bin)

          case (SCORE_FLUX)
            score = score * flux_deriv

          case (SCORE_TOTAL)
            if (i_nuclide == -1 .and. &
                 material_id(p % material) == deriv % diff_material .and. &
                 material_xs % total > ZERO) then
              cum_dsig = ZERO
                do l = 1, material_nuclide_size(p % material)
                  i_nuc = material_nuclide(p % material, l)
                  associate (nuc => nuclides(i_nuc))
                    if (multipole_in_range(nuc % ptr, p % last_E) .and. &
                         micro_xs(i_nuc) % total > ZERO) then
                      call multipole_deriv_eval(nuc % ptr, p % last_E, &
                           p % sqrtkT, dsig_s, dsig_a, dsig_f)
                      cum_dsig = cum_dsig + (dsig_s + dsig_a) &
                           * material_atom_density(p % material, l)
                    end if
                  end associate
                end do
              score = score * (flux_deriv &
                   + cum_dsig / material_xs % total)
            else if (material_id(p % material) == deriv % diff_material &
                 .and. material_xs % total > ZERO) then
              dsig_s = ZERO
              dsig_a = ZERO
              associate (nuc => nuclides(i_nuclide+1))
                if (multipole_in_range(nuc % ptr, p % last_E)) then
                  call multipole_deriv_eval(nuc % ptr, p % last_E, &
                       p % sqrtkT, dsig_s, dsig_a, dsig_f)
                end if
              end associate
              score = score * (flux_deriv &
                   + (dsig_s + dsig_a) / micro_xs(i_nuclide+1) % total)
            else
              score = score * flux_deriv
            end if

          case (SCORE_SCATTER)
            if (i_nuclide == -1 .and. &
                 material_id(p % material) == deriv % diff_material .and. &
                 (material_xs % total - material_xs % absorption) > ZERO) then
              cum_dsig = ZERO
                do l = 1, material_nuclide_size(p % material)
                  i_nuc = material_nuclide(p % material, l)
                  associate (nuc => nuclides(i_nuc))
                    if (multipole_in_range(nuc % ptr, p % last_E) .and. &
                         (micro_xs(i_nuc) % total &
                         - micro_xs(i_nuc) % absorption) > ZERO) then
                      call multipole_deriv_eval(nuc % ptr, p % last_E, &
                           p % sqrtkT, dsig_s, dsig_a, dsig_f)
                      cum_dsig = cum_dsig + dsig_s * material_atom_density(p % material, l)
                    end if
                  end associate
                end do
              score = score * (flux_deriv + cum_dsig &
                   / (material_xs % total - material_xs % absorption))
            else if ( material_id(p % material) == deriv % diff_material &
                 .and. (material_xs % total - material_xs % absorption) > ZERO)&
                 then
              dsig_s = ZERO
              associate (nuc => nuclides(i_nuclide+1))
                if (multipole_in_range(nuc % ptr, p % last_E)) then
                  call multipole_deriv_eval(nuc % ptr, p % last_E, &
                       p % sqrtkT, dsig_s, dsig_a, dsig_f)
                end if
              end associate
              score = score * (flux_deriv + dsig_s &
                   / (micro_xs(i_nuclide+1) % total &
                   - micro_xs(i_nuclide+1) % absorption))
            else
              score = score * flux_deriv
            end if

          case (SCORE_ABSORPTION)
            if (i_nuclide == -1 .and. &
                 material_id(p % material) == deriv % diff_material .and. &
                 material_xs % absorption > ZERO) then
              cum_dsig = ZERO
                do l = 1, material_nuclide_size(p % material)
                  i_nuc = material_nuclide(p % material, l)
                  associate (nuc => nuclides(i_nuc))
                    if (multipole_in_range(nuc % ptr, p % last_E) .and. &
                         micro_xs(i_nuc) % absorption > ZERO) then
                      call multipole_deriv_eval(nuc % ptr, p % last_E, &
                           p % sqrtkT, dsig_s, dsig_a, dsig_f)
                      cum_dsig = cum_dsig + dsig_a * material_atom_density(p % material, l)
                    end if
                  end associate
                end do
              score = score * (flux_deriv &
                   + cum_dsig / material_xs % absorption)
            else if (material_id(p % material) == deriv % diff_material &
                 .and. material_xs % absorption > ZERO) then
              dsig_a = ZERO
              associate (nuc => nuclides(i_nuclide+1))
                if (multipole_in_range(nuc % ptr, p % last_E)) then
                  call multipole_deriv_eval(nuc % ptr, p % last_E, &
                       p % sqrtkT, dsig_s, dsig_a, dsig_f)
                end if
              end associate
              score = score * (flux_deriv &
                   + dsig_a / micro_xs(i_nuclide+1) % absorption)
            else
              score = score * flux_deriv
            end if

          case (SCORE_FISSION)
            if (i_nuclide == -1 .and. &
                 material_id(p % material) == deriv % diff_material .and. &
                 material_xs % fission > ZERO) then
              cum_dsig = ZERO
                do l = 1, material_nuclide_size(p % material)
                  i_nuc = material_nuclide(p % material, l)
                  associate (nuc => nuclides(i_nuc))
                    if (multipole_in_range(nuc % ptr, p % last_E) .and. &
                         micro_xs(i_nuc) % fission > ZERO) then
                      call multipole_deriv_eval(nuc % ptr, p % last_E, &
                           p % sqrtkT, dsig_s, dsig_a, dsig_f)
                      cum_dsig = cum_dsig + dsig_f * material_atom_density(p % material, l)
                    end if
                  end associate
                end do
              score = score * (flux_deriv &
                   + cum_dsig / material_xs % fission)
            else if (material_id(p % material) == deriv % diff_material &
                 .and. material_xs % fission > ZERO) then
              dsig_f = ZERO
              associate (nuc => nuclides(i_nuclide+1))
                if (multipole_in_range(nuc % ptr, p % last_E)) then
                  call multipole_deriv_eval(nuc % ptr, p % last_E, &
                       p % sqrtkT, dsig_s, dsig_a, dsig_f)
                end if
              end associate
              score = score * (flux_deriv &
                   + dsig_f / micro_xs(i_nuclide+1) % fission)
            else
              score = score * flux_deriv
            end if

          case (SCORE_NU_FISSION)
            if (i_nuclide == -1 .and. &
                 material_id(p % material) == deriv % diff_material .and. &
                 material_xs % nu_fission > ZERO) then
              cum_dsig = ZERO
                do l = 1, material_nuclide_size(p % material)
                  i_nuc = material_nuclide(p % material, l)
                  associate (nuc => nuclides(i_nuc))
                    if (multipole_in_range(nuc % ptr, p % last_E) .and. &
                         micro_xs(i_nuc) % nu_fission > ZERO) then
                      call multipole_deriv_eval(nuc % ptr, p % last_E, &
                           p % sqrtkT, dsig_s, dsig_a, dsig_f)
                      cum_dsig = cum_dsig + dsig_f * material_atom_density(p % material, l) &
                           * micro_xs(i_nuc) % nu_fission &
                           / micro_xs(i_nuc) % fission
                    end if
                  end associate
                end do
              score = score * (flux_deriv &
                   + cum_dsig / material_xs % nu_fission)
            else if (material_id(p % material) == deriv % diff_material &
                 .and. material_xs % nu_fission > ZERO) then
              dsig_f = ZERO
              associate (nuc => nuclides(i_nuclide+1))
                if (multipole_in_range(nuc % ptr, p % last_E)) then
                  call multipole_deriv_eval(nuc % ptr, p % last_E, &
                       p % sqrtkT, dsig_s, dsig_a, dsig_f)
                end if
              end associate
              score = score * (flux_deriv &
                   + dsig_f / micro_xs(i_nuclide+1) % fission)
            else
              score = score * flux_deriv
            end if

          case default
            call fatal_error('Tally derivative not defined for a score on &
                 &tally ' // trim(to_str(t % id())))
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

    type(TallyDerivative), pointer :: deriv
    integer :: i, l
    real(8) :: dsig_s, dsig_a, dsig_f

    ! A void material cannot be perturbed so it will not affect flux derivatives
    if (p % material == MATERIAL_VOID) return

    do i = 0, n_tally_derivs() - 1
      !associate(deriv => tally_derivs(i))
      deriv => tally_deriv_c(i)
        select case (deriv % variable)

        case (DIFF_DENSITY)
            if (material_id(p % material) == deriv % diff_material) then
              ! phi is proportional to e^(-Sigma_tot * dist)
              ! (1 / phi) * (d_phi / d_rho) = - (d_Sigma_tot / d_rho) * dist
              ! (1 / phi) * (d_phi / d_rho) = - Sigma_tot / rho * dist
              deriv % flux_deriv = deriv % flux_deriv &
                    - distance * material_xs % total / material_density_gpcc(p % material)
            end if

        case (DIFF_NUCLIDE_DENSITY)
            if (material_id(p % material) == deriv % diff_material) then
              ! phi is proportional to e^(-Sigma_tot * dist)
              ! (1 / phi) * (d_phi / d_N) = - (d_Sigma_tot / d_N) * dist
              ! (1 / phi) * (d_phi / d_N) = - sigma_tot * dist
              deriv % flux_deriv = deriv % flux_deriv &
                   - distance * micro_xs(deriv % diff_nuclide) % total
            end if

        case (DIFF_TEMPERATURE)
            if (material_id(p % material) == deriv % diff_material) then
              do l=1, material_nuclide_size(p % material)
                associate (nuc => nuclides(material_nuclide(p % material, l)))
                  if (multipole_in_range(nuc % ptr, p % E)) then
                    ! phi is proportional to e^(-Sigma_tot * dist)
                    ! (1 / phi) * (d_phi / d_T) = - (d_Sigma_tot / d_T) * dist
                    ! (1 / phi) * (d_phi / d_T) = - N (d_sigma_tot / d_T) * dist
                    call multipole_deriv_eval(nuc % ptr, p % E, &
                         p % sqrtkT, dsig_s, dsig_a, dsig_f)
                    deriv % flux_deriv = deriv % flux_deriv &
                         - distance * (dsig_s + dsig_a) * material_atom_density(p % material, l)
                  end if
                end associate
              end do
            end if
        end select
      !end associate
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

    type(TallyDerivative), pointer :: deriv
    integer :: i, j, l, i_nuc
    real(8) :: dsig_s, dsig_a, dsig_f

    ! A void material cannot be perturbed so it will not affect flux derivatives
    if (p % material == MATERIAL_VOID) return

    do i = 0, n_tally_derivs() - 1
      !associate(deriv => tally_derivs(i))
      deriv => tally_deriv_c(i)
        select case (deriv % variable)

        case (DIFF_DENSITY)
            if (material_id(p % material) == deriv % diff_material) then
              ! phi is proportional to Sigma_s
              ! (1 / phi) * (d_phi / d_rho) = (d_Sigma_s / d_rho) / Sigma_s
              ! (1 / phi) * (d_phi / d_rho) = 1 / rho
              deriv % flux_deriv = deriv % flux_deriv &
                   + ONE / material_density_gpcc(p % material)
            end if

        case (DIFF_NUCLIDE_DENSITY)
            if (material_id(p % material) == deriv % diff_material &
                 .and. p % event_nuclide == deriv % diff_nuclide) then
              ! Find the index in this material for the diff_nuclide.
              do j = 1, material_nuclide_size(p % material)
                if (material_nuclide(p % material, j) == deriv % diff_nuclide) exit
              end do
              ! Make sure we found the nuclide.
              if (material_nuclide(p % material, j) /= deriv % diff_nuclide) then
                call fatal_error("Couldn't find the right nuclide.")
              end if
              ! phi is proportional to Sigma_s
              ! (1 / phi) * (d_phi / d_N) = (d_Sigma_s / d_N) / Sigma_s
              ! (1 / phi) * (d_phi / d_N) = sigma_s / Sigma_s
              ! (1 / phi) * (d_phi / d_N) = 1 / N
              deriv % flux_deriv = deriv % flux_deriv &
                   + ONE / material_atom_density(p % material, j)
            end if

        case (DIFF_TEMPERATURE)
            if (material_id(p % material) == deriv % diff_material) then
              do l=1, material_nuclide_size(p % material)
                i_nuc = material_nuclide(p % material, l)
                associate (nuc => nuclides(i_nuc))
                  if (i_nuc == p % event_nuclide .and. &
                       multipole_in_range(nuc % ptr, p % last_E)) then
                    ! phi is proportional to Sigma_s
                    ! (1 / phi) * (d_phi / d_T) = (d_Sigma_s / d_T) / Sigma_s
                    ! (1 / phi) * (d_phi / d_T) = (d_sigma_s / d_T) / sigma_s
                    call multipole_deriv_eval(nuc % ptr, p % last_E, &
                         p % sqrtkT, dsig_s, dsig_a, dsig_f)
                    deriv % flux_deriv = deriv % flux_deriv + dsig_s&
                         / (micro_xs(i_nuc) % total &
                         - micro_xs(i_nuc) % absorption)
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
        end select
      !end associate
    end do
  end subroutine score_collision_derivative

!===============================================================================
! ACCUMULATE_TALLIES accumulates the sum of the contributions from each history
! within the batch to a new random variable
!===============================================================================

  subroutine accumulate_tallies() bind(C)

    integer :: i
    real(C_DOUBLE) :: k_col ! Copy of batch collision estimate of keff
    real(C_DOUBLE) :: k_abs ! Copy of batch absorption estimate of keff
    real(C_DOUBLE) :: k_tra ! Copy of batch tracklength estimate of keff
    real(C_DOUBLE) :: val

#ifdef OPENMC_MPI
    interface
      subroutine reduce_tally_results() bind(C)
      end subroutine
    end interface

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

    ! Accumulate results for each tally
    do i = 1, active_tallies_size()
      call tallies(active_tallies_data(i)) % obj % accumulate()
    end do

  end subroutine accumulate_tallies

!===============================================================================
! SETUP_ACTIVE_TALLIES
!===============================================================================

  subroutine setup_active_tallies() bind(C)

    integer :: i
    integer(C_INT) :: err
    logical(C_BOOL) :: active

    interface
      subroutine setup_active_tallies_c() bind(C)
      end subroutine
    end interface

    call setup_active_tallies_c()

    do i = 1, n_tallies
      associate (t => tallies(i) % obj)
        err = openmc_tally_get_active(i, active)
        if (active) then
          ! Check if tally contains depletion reactions and if so, set flag
          if (t % depletion_rx()) need_depletion_rx = .true.
        end if
      end associate
    end do

  end subroutine setup_active_tallies

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_tally_allocate(index, type) result(err) bind(C)
    ! Set the type of the tally
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR), intent(in) :: type(*)
    integer(C_INT) :: err

    integer(C_INT32_T) :: empty(0)
    character(:), allocatable :: type_

    interface
      function tally_pointer(indx) bind(C) result(ptr)
        import C_INT, C_PTR
        integer(C_INT), value :: indx
        type(C_PTR)           :: ptr
      end function
    end interface

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

        ! Assign the pointer to the C++ tally
        tallies(index) % obj % ptr = tally_pointer(index - 1)

        ! When a tally is allocated, set it to have 0 filters
        err = openmc_tally_set_filters(index, 0, empty)
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in tallies array is out of bounds.")
    end if
  end function openmc_tally_allocate

end module tally
