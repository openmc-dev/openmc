module lists

  implicit none
  private
  public :: RealList

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! REALELEM contains one element of a linked list of reals
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type RealElem

    real(8) :: val ! value of this element
    type(RealElem), pointer :: next => null() ! next element in list

  end type RealElem

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! REALLIST is a linked list object containing reals
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type RealList

    type(RealElem) :: first ! first element
    type(RealElem), pointer :: head => null() ! current element

  ! type-bound procedures
  contains

    ! deallocates a linked list of reals
    procedure :: clear => clear_real_list

  end type RealList

contains

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CLEAR_REAL_LIST deallocates a linked list of reals
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine clear_real_list(this)

    class(RealList), target, intent(inout) :: this ! real linked list object
    type(RealElem), pointer :: next => null() ! placeholder for next element

    ! start at beginning of list
    this % head => this % first

    ! while the current element is associated
    do while(associated(this % head))

      ! stash the next element
      next => this % head % next

      ! if the current element is the first, don't need to deallocated
      if (.not. associated(this % head, target=this % first))&
           deallocate(this % head)

      ! point to the next element
      this % head => next

    end do

    ! finish it off
    nullify(next)

  end subroutine clear_real_list

end module lists
