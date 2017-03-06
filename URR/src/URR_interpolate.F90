module URR_interpolate

  use URR_openmc_wrapper, only: prn
  use URR_constants, only: ZERO,&
                           ONE,&
                           HISTOGRAM,&
                           LINEAR_LINEAR,&
                           LINEAR_LOG,&
                           LOG_LINEAR,&
                           LOG_LOG,&
                           SQRT_LINEAR,&
                           SQRT_LOG,&
                           STATISTICAL,&
                           LOW_NEIGHBOR
  use URR_error, only: exit_status,&
                       EXIT_FAILURE

  implicit none
  private
  public :: interp_factor,&
            interpolate

contains


!> Computes an interpolation factor
  function interp_factor(val, val_low, val_up, scheme) result(factor)

    integer :: scheme ! interpolation scheme
    real(8) :: val     ! value we're interpolating to
    real(8) :: val_low ! lower bounding value
    real(8) :: val_up  ! upper bounding value
    real(8) :: factor  ! interpolation factor
    real(8) :: xi      ! prn for statistical interpolation
    real(8) :: prob    ! probability for statistical interpolation

    select case(scheme)
    case(HISTOGRAM)
      factor = ONE

    case(LINEAR_LINEAR)
      factor = (val - val_low) / (val_up - val_low)

    case(LINEAR_LOG)
      factor = (val - val_low) / (val_up - val_low)

    case(LOG_LINEAR)
      factor = log(val / val_low) / log(val_up / val_low)

    case(LOG_LOG)
      factor = log(val / val_low) / log(val_up / val_low)

    case(SQRT_LINEAR)
      factor = (sqrt(val) - sqrt(val_low)) / (sqrt(val_up) - sqrt(val_low))

    case(SQRT_LOG)
      factor = (sqrt(val) - sqrt(val_low)) / (sqrt(val_up) - sqrt(val_low))

    case(STATISTICAL)
      xi = prn()
      prob = (val - val_low) / (val_up - val_low)
      if (xi > prob) then
        factor = ZERO
      else
        factor = ONE
      end if

    case(LOW_NEIGHBOR)
      factor = ZERO

    case default
      call exit_status(EXIT_FAILURE, 'Interpolation scheme not recognized')

    end select

  end function interp_factor


!> Computes an interpolated (or extrapolated) value
  function interpolate(factor, val_low, val_up, scheme) result(val)

    integer :: scheme ! interpolation scheme
    real(8) :: factor  ! interpolation factor
    real(8) :: val_low ! lower bounding value
    real(8) :: val_up  ! upper bounding value
    real(8) :: val     ! interpolated value

    select case(scheme)
    case(HISTOGRAM)
      val = val_low

    case(LINEAR_LINEAR)
      val = val_low + factor * (val_up - val_low)

    case(LINEAR_LOG)
      val = val_low * exp(factor * log(val_up / val_low))

    case(LOG_LINEAR)
      val = val_low + factor * (val_up - val_low)

    case(LOG_LOG)
      val = val_low * exp(factor * log(val_up / val_low))

    case(SQRT_LINEAR)
      val = val_low + factor * (val_up - val_low)

    case(SQRT_LOG)
      val = val_low * exp(factor * log(val_up / val_low))

    case(STATISTICAL)
      val = val_low + factor * (val_up - val_low)

    case(LOW_NEIGHBOR)
      val = val_low

    case default
      call exit_status(EXIT_FAILURE, 'Interpolation scheme not recognized')

    end select

  end function interpolate
  

end module URR_interpolate
