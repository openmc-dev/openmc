module trigger

  use, intrinsic :: ISO_C_BINDING

  use constants
  use eigenvalue,     only: openmc_get_keff
  use endf,           only: reaction_name
  use error,          only: warning, write_message
  use string,         only: to_str
  use message_passing, only: master
  use settings
  use simulation_header
  use trigger_header
  use tally,          only: TallyObject
  use tally_header, only: tallies, n_tallies

  implicit none

  interface
    subroutine get_tally_uncertainty(i_tally, score_index, filter_index, &
         std_dev, rel_err) bind(C)
      import C_INT, C_DOUBLE
      integer(C_INT), value :: i_tally
      integer(C_INT), value :: score_index
      integer(C_INT), value :: filter_index
      real(C_DOUBLE) :: std_dev
      real(C_DOUBLE) :: rel_err
    end subroutine
  end interface

contains

!===============================================================================
! CHECK_TRIGGERS checks any user-specified precision triggers' for convergence
! and predicts the number of remainining batches to convergence.
!===============================================================================

  subroutine check_triggers() bind(C)

    implicit none

    real(8) :: keff_ratio ! uncertainty/threshold ratio for keff
    real(8) :: tally_ratio ! max tally uncertainty/threshold ratio
    integer :: tally_id ! id for tally with max ratio
    integer :: score ! tally score with max ratio
    integer :: n_pred_batches  ! predicted # batches to satisfy all triggers

    if (.not. master) return

    ! Checks if current_batch is one for which the triggers must be checked
    if (current_batch < n_batches .or. (.not. trigger_on)) return
    if (mod((current_batch - n_batches), n_batch_interval) /= 0 .and. &
         current_batch /= n_max_batches) return

    ! By default, assume all triggers are satisfied
    satisfy_triggers = .true.

    ! Check the eigenvalue and tally triggers
    call check_keff_trigger(keff_ratio)
    call check_tally_triggers(tally_ratio, tally_id, score)

    ! Alert the user if the triggers are satisfied
    if (satisfy_triggers) then
      call write_message("Triggers satisfied for batch " // &
           trim(to_str(current_batch)), 7)
      return
    end if

    ! Specify which trigger is unsatisfied
    if (keff_ratio >= tally_ratio) then
      call write_message("Triggers unsatisfied, max unc./thresh. is " // &
           trim(to_str(keff_ratio)) // " for eigenvalue", 7)
    else
      call write_message("Triggers unsatisfied, max unc./thresh. is " // &
           trim(to_str(tally_ratio)) // " for " // trim(reaction_name(score)) &
           // " in tally " // trim(to_str(tally_id)), 7)
    end if

    ! Estimate batches till triggers are satisfied
    if (pred_batches) then
      ! Estimate the number of remaining batches to convergence
      ! The prediction uses the fact that tally variances are proportional
      ! to 1/N where N is the number of the batches/particles
      n_batch_interval = int((current_batch-n_inactive) * &
           (max(keff_ratio, tally_ratio) ** 2)) + n_inactive-n_batches + 1
      n_pred_batches = n_batch_interval + n_batches

      ! Write the predicted number of batches for the user
      if (n_pred_batches > n_max_batches) then
        call warning("The estimated number of batches is " // &
             trim(to_str(n_pred_batches)) // &
             " --  greater than max batches. ")
      else
        call write_message("The estimated number of batches is " // &
             trim(to_str(n_pred_batches)), 7)
      end if
    end if
  end subroutine check_triggers

!===============================================================================
! CHECK_KEFF_TRIGGER computes the uncertainty/threshold ratio for the eigenvalue
! trigger and updates the global satisfy_tiggers variable if the trigger is
! unsatisfied.
!===============================================================================

  subroutine check_keff_trigger(ratio)
    real(8), intent(out) :: ratio
    integer(C_INT) :: err
    real(C_DOUBLE) :: k_combined(2)
    real(8) :: uncertainty

    ratio = 0
    if (run_mode == MODE_EIGENVALUE) then
      if (keff_trigger % trigger_type /= 0) then
        err = openmc_get_keff(k_combined)
        select case (keff_trigger % trigger_type)
        case(VARIANCE)
          uncertainty = k_combined(2) ** 2
        case(STANDARD_DEVIATION)
          uncertainty = k_combined(2)
        case default
          uncertainty = k_combined(2) / k_combined(1)
        end select

        ! If uncertainty is above threshold, store uncertainty ratio
        if (uncertainty > keff_trigger % threshold) then
          satisfy_triggers = .false.
          if (keff_trigger % trigger_type == VARIANCE) then
            ratio = sqrt(uncertainty / keff_trigger % threshold)
          else
            ratio = uncertainty / keff_trigger % threshold
          end if
        end if
      end if
    end if
  end subroutine check_keff_trigger

!===============================================================================
! CHECK_TALLY_TRIGGERS computes the uncertainty/threshold ratio for all tally
! triggers and updates the global satisfy_tiggers variable if any trigger is
! unsatisfied.
!===============================================================================

  subroutine check_tally_triggers(max_ratio, tally_id, score)

    ! Variables to reflect distance to trigger convergence criteria
    real(8), intent(out) :: max_ratio       ! max uncertainty/thresh ratio
    integer, intent(out) :: tally_id        ! id for tally with max ratio
    integer, intent(out) :: score

    integer :: i              ! index in tallies array
    integer :: j              ! index in tally filters
    integer :: n              ! loop index for nuclides
    integer :: s              ! loop index for triggers
    integer :: filter_index   ! index in results array for filters
    integer :: score_index    ! scoring bin index
    integer(C_INT) :: err
    real(8) :: uncertainty    ! trigger uncertainty
    real(8) :: std_dev = ZERO ! trigger standard deviation
    real(8) :: rel_err = ZERO ! trigger relative error
    real(8) :: ratio          ! ratio of the uncertainty/trigger threshold
    real(C_DOUBLE) :: k_combined(2)

    ! Initialize tally trigger maximum uncertainty ratio to zero
    max_ratio = 0

    ! Compute uncertainties for all tallies, scores with triggers
    TALLY_LOOP: do i = 1, n_tallies
      associate (t => tallies(i) % obj)

      ! Cycle through if only one batch has been simumlate
      if (t % n_realizations == 1) then
        cycle TALLY_LOOP
      end if

      TRIGGER_LOOP: do s = 1, t % n_triggers
        associate (trigger => t % triggers(s))

          ! Initialize trigger uncertainties to zero
          trigger % std_dev = ZERO
          trigger % rel_err = ZERO
          trigger % variance = ZERO

          FILTER_LOOP: do filter_index = 1, t % n_filter_bins()

            ! Initialize score index
            score_index = trigger % score_index

            ! Initialize score bin index
            NUCLIDE_LOOP: do n = 1, t % n_nuclide_bins()

              call get_tally_uncertainty(i, score_index, filter_index, &
                   std_dev, rel_err)

              if (trigger % variance < variance) then
                trigger % variance = std_dev ** 2
              end if
              if (trigger % std_dev < std_dev) then
                trigger % std_dev = std_dev
              end if
              if (trigger % rel_err < rel_err) then
                trigger % rel_err = rel_err
              end if

              select case (trigger % type)
              case(VARIANCE)
                uncertainty = trigger % variance
              case(STANDARD_DEVIATION)
                uncertainty = trigger % std_dev
              case default
                uncertainty = trigger % rel_err
              end select

              if (uncertainty > trigger % threshold) then
                satisfy_triggers = .false.

                if (trigger % type == VARIANCE) then
                  ratio = sqrt(uncertainty / trigger % threshold)
                else
                  ratio = uncertainty / trigger % threshold
                end if

                if (max_ratio < ratio) then
                  max_ratio = ratio
                  score = t % score_bins(trigger % score_index)
                  tally_id = t % id
                end if
              end if
            end do NUCLIDE_LOOP
            if (t % n_filters() == 0) exit FILTER_LOOP
          end do FILTER_LOOP
        end associate
      end do TRIGGER_LOOP
      end associate
    end do TALLY_LOOP
  end subroutine check_tally_triggers

end module trigger
