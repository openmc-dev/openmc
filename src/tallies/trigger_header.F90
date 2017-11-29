module trigger_header

  use constants, only: NONE, N_FILTER_TYPES, ZERO

  implicit none
  private

!===============================================================================
! TRIGGEROBJECT stores the variance, relative error and standard deviation
! for some user-specified trigger.
!===============================================================================
  type, public :: TriggerObject
    integer            :: type          ! "variance", "std_dev" or "rel_err"
    real(8)            :: threshold     ! a convergence threshold
    character(len=52)  :: score_name    ! the name of the score
    integer            :: score_index   ! the index of the score
    real(8)            :: variance = ZERO ! temp variance container
    real(8)            :: std_dev  = ZERO ! temp std. dev. container
    real(8)            :: rel_err  = ZERO ! temp rel. err. container
  end type TriggerObject

!===============================================================================
! KTRIGGER describes a user-specified precision trigger for k-effective
!===============================================================================
  type, public :: KTrigger
    integer    :: trigger_type = 0
    real(8)    :: threshold    = ZERO
  end type KTrigger

  type(KTrigger), public :: keff_trigger  ! trigger for k-effective

end module trigger_header
