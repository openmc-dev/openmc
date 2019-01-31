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
    function check_keff_trigger() bind(C) result(ratio)
      import C_DOUBLE
      real(C_DOUBLE) :: ratio
    end function

    subroutine check_tally_triggers(ratio, tally_id, score) bind(C)
      import C_DOUBLE, C_INT
      real(C_DOUBLE) :: ratio
      integer(C_INT) :: tally_id
      integer(C_INT) :: score
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

    ! Check the eigenvalue and tally triggers
    keff_ratio = check_keff_trigger()
    call check_tally_triggers(tally_ratio, tally_id, score)

    if (max(keff_ratio, tally_ratio) <= ONE) then
      satisfy_triggers = .true.
    else
      satisfy_triggers = .false.
    end if

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

  function n_tally_triggers(i_tally) bind(C) result(n)
    integer(C_INT), value :: i_tally
    integer(C_INT) :: n
    n = tallies(i_tally) % obj % n_triggers
  end function

  function get_tally_trigger(i_tally, i_trig) bind(C) result(trigger)
    integer(C_INT), value :: i_tally
    integer(C_INT), value :: i_trig
    type(C_PTR) :: trigger
    trigger = C_LOC(tallies(i_tally) % obj % triggers(i_trig))
  end function

end module trigger
