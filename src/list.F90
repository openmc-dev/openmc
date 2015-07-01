module list

  implicit none

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! REALLISTELEM contains one element of a linked list of reals
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type RealListElem

    real(8) :: val ! value of this element
    type(RealListElem), pointer :: next => null() ! next element in list

  end type RealListElem

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! REALLIST is a linked list object containing reals
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type RealList

    type(RealListElem) :: first ! first element
    type(RealListElem), pointer :: head => null() ! current element

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
    type(RealListElem), pointer :: next => null() ! placeholder for next element

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

end module list
