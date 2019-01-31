module trigger_header

  use, intrinsic :: ISO_C_BINDING

  use constants, only: NONE, ZERO

  implicit none
  private

!===============================================================================
! TRIGGEROBJECT stores the variance, relative error and standard deviation
! for some user-specified trigger.
!===============================================================================
  type, public, bind(C) :: TriggerObject
    integer(C_INT) :: type          ! "variance", "std_dev" or "rel_err"
    real(C_DOUBLE) :: threshold     ! a convergence threshold
    integer(C_INT) :: score_index   ! the index of the score
    real(C_DOUBLE) :: variance = ZERO ! temp variance container
    real(C_DOUBLE) :: std_dev  = ZERO ! temp std. dev. container
    real(C_DOUBLE) :: rel_err  = ZERO ! temp rel. err. container
  end type TriggerObject

!===============================================================================
! KTRIGGER describes a user-specified precision trigger for k-effective
!===============================================================================
  type, public, bind(C) :: KTrigger
    integer(C_INT)    :: trigger_type = 0
    real(C_DOUBLE)    :: threshold    = ZERO
  end type KTrigger

  type(KTrigger), public, bind(C) :: keff_trigger  ! trigger for k-effective

end module trigger_header
