module dd_header 

  implicit none
  private
  public :: allocate_dd, deallocate_dd

  type, public :: dd_type


  end type dd_type

contains

!==============================================================================
! ALLOCATE_DD allocates all data in dd type
!==============================================================================

  subroutine allocate_dd(this)

    type(dd_type), intent(inout) :: this      ! dd instance


  end subroutine allocate_dd

!===============================================================================
! DEALLOCATE_DD frees all memory of dd type 
!===============================================================================

  subroutine deallocate_dd(this)

    type(dd_type), intent(inout) :: this ! dd instance


  end subroutine deallocate_dd

end module dd_header
