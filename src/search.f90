module search

  use constants, only: ONE, MAX_LINE_LEN
  use error,     only: fatal_error

contains

!===============================================================================
! BINARY_SEARCH performs a binary search of an array to find where a specific
! value lies in the array. This is used extensively for energy grid searching
!===============================================================================

  function binary_search(array, n, val) result(index)

    real(8), intent(in) :: array(n)
    integer, intent(in) :: n
    real(8), intent(in) :: val
    integer :: index

    integer :: L
    integer :: R
    real(8) :: testval
    character(MAX_LINE_LEN) :: msg

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
       msg = "Value outside of array during binary search"
       call fatal_error(msg)
    end if
    
    do while (R - L > 1)
       
       ! Check boundaries
       if (val > array(L) .and. val < array(L+1)) then
          index = L
          return
       elseif (val > array(R-1) .and. val < array(R)) then
          index = R-1
          return
       end if

       ! Find values at midpoint
       index = L + (R - L)/2
       testval = array(index)
       if (val > testval) then
          L = index
       elseif (val < testval) then
          R = index
       end if
    end do

    index = L

  end function binary_search

!===============================================================================
! INTERPOLATE
!===============================================================================

  function interpolate(array, n, index, f) result(val)

    real(8), intent(in) :: array(n)
    integer, intent(in) :: n
    integer, intent(in) :: index
    real(8), intent(in) :: f
    real(8)             :: val

    val = (ONE-f) * array(index) + f * array(index+1)
    
  end function interpolate

end module search
