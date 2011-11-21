module search

  use constants, only: ONE, MAX_LINE_LEN
  use error,     only: fatal_error
  use global,    only: message

contains

!===============================================================================
! BINARY_SEARCH performs a binary search of an array to find where a specific
! value lies in the array. This is used extensively for energy grid searching
!===============================================================================

  function binary_search(array, n, val) result(array_index)

    integer, intent(in) :: n
    real(8), intent(in) :: array(n)
    real(8), intent(in) :: val
    integer             :: array_index

    integer :: L
    integer :: R
    real(8) :: testval

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
       message = "Value outside of array during binary search"
       call fatal_error()
    end if
    
    do while (R - L > 1)
       
       ! Check boundaries
       if (val > array(L) .and. val < array(L+1)) then
          array_index = L
          return
       elseif (val > array(R-1) .and. val < array(R)) then
          array_index = R - 1
          return
       end if

       ! Find values at midpoint
       array_index = L + (R - L)/2
       testval = array(array_index)
       if (val > testval) then
          L = array_index
       elseif (val < testval) then
          R = array_index
       end if
    end do

    array_index = L

  end function binary_search

end module search
