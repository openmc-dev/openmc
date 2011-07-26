module timing

  use global
  use types, only: TimerObj
  use error, only: warning

  implicit none

contains

!===============================================================================
! TIMER_START
!===============================================================================

  subroutine timer_start(timer)

    type(TimerObj), intent(inout) :: timer

    character(max_line_len) :: msg ! warning message

    ! Check if timer is already running
    if (timer % running) then
       msg = "Tried to start a timer that was already running!"
       call warning(msg)
       return
    end if

    ! Turn timer on and measure starting time
    timer % running = .true.
    call system_clock(timer % start_counts)

  end subroutine timer_start

!===============================================================================
! TIMER_GET_VALUE
!===============================================================================

  function timer_get_value(timer) result(elapsed)

    type(TimerObj), intent(in) :: timer   ! the timer
    real(8)                    :: elapsed ! total elapsed time

    integer :: end_counts   ! current number of counts
    integer :: count_rate   ! system-dependent counting rate
    real    :: elapsed_time ! elapsed time since last start

    if (timer % running) then
       call system_clock(end_counts, count_rate)
       elapsed_time = real(end_counts - timer % start_counts)/real(count_rate)
       elapsed = timer % elapsed + elapsed_time
    else
       elapsed = timer % elapsed
    end if

  end function timer_get_value

!===============================================================================
! TIMER_STOP
!===============================================================================

  subroutine timer_stop(timer)

    type(TimerObj), intent(inout) :: timer

    character(max_line_len) :: msg

    ! Check to make sure timer was running
    if (.not. timer % running) then
       msg = "Tried to stop a timer that was not running!"
       call warning(msg)
       return
    end if

    ! Stop timer and add time
    timer % elapsed = timer_get_value(timer)
    timer % running = .false.

  end subroutine timer_stop

!===============================================================================
! TIMER_RESET
!===============================================================================

  subroutine timer_reset(timer)

    type(TimerObj), intent(inout) :: timer

    timer % running      = .false.
    timer % start_counts = 0
    timer % elapsed      = ZERO

  end subroutine timer_reset

end module timing
