module trigger_header

  use constants, only: NONE, N_FILTER_TYPES

  implicit none

!===============================================================================
! TRIGGEROBJECT stores the variance, relative error and standard deviation
! for some user-specified trigger.
!===============================================================================
  type TriggerObject
    integer            :: type          ! "variance", "std_dev" or "rel_err"
    real(8)            :: threshold     ! a convergence threshold
    character(len=52)  :: score_name    ! the name of the score
    integer            :: score_index   ! the index of the score
    real(8)            :: variance=0.0  ! temp variance container
    real(8)            :: std_dev =0.0  ! temp std. dev. container
    real(8)            :: rel_err =0.0  ! temp rel. err. container
  end type TriggerObject

!===============================================================================
! KTRIGGER describes a user-specified precision trigger for k-effective
!===============================================================================
  type KTrigger
    integer    :: trigger_type = 0
    real(8)    :: threshold    = 0
  end type KTrigger

end module trigger_header
