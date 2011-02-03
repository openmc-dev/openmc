module search

contains

!=====================================================================
! BINARY_SEARCH performs a binary search of an array to find where a
! specific value lies in the array. This is used extensively for
! energy grid searching
!=====================================================================

  function binary_search( array, n, val ) result( index )

    real(8), intent(in) :: array(n)
    integer, intent(in) :: n
    real(8), intent(in) :: val
    integer :: index

    integer :: L
    integer :: R
    real(8) :: testval

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
       ! error
    end if
    
    do while (R - L > 1)
       
       ! Check boundaries
       if (val > array(L) .and. val < array(L+1)) then
          index = L
          return
       elseif (val > array(R+1) .and. val < array(R)) then
          index = R
          return
       end if

       ! Find values at midpoint
       index = L + (R - L)/2
       testval = array(index)
       if (val > testval) then
          L = index + 1
       elseif (val < testval) then
          R = index - 1
       end if
    end do

    index = L

  end function binary_search

end module search
