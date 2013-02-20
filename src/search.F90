module search

  use error,     only: fatal_error
  use global,    only: message

  integer, parameter :: MAX_ITERATION = 64

  interface binary_search
    module procedure binary_search_real, binary_search_int4, binary_search_int8
  end interface binary_search

contains

!===============================================================================
! BINARY_SEARCH performs a binary search of an array to find where a specific
! value lies in the array. This is used extensively for energy grid searching
!===============================================================================

  function binary_search_real(array, n, val) result(array_index)

    integer, intent(in) :: n
    real(8), intent(in) :: array(n)
    real(8), intent(in) :: val
    integer             :: array_index

    integer :: L
    integer :: R

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
      message = "Value outside of array during binary search"
      call fatal_error()
    end if

    do while (R > L)
      ! Find values at midpoint
      array_index = (R + L)/2
      if (val > array(array_index)) then
        L = array_index + 1
      else
        R = array_index
      end if
    end do

    array_index = L - 1

  end function binary_search_real

  function binary_search_int4(array, n, val) result(array_index)

    integer, intent(in) :: n
    integer, intent(in) :: array(n)
    integer, intent(in) :: val
    integer             :: array_index

    integer :: L
    integer :: R

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
      message = "Value outside of array during binary search"
      call fatal_error()
    end if

    do while (R > L)
      ! Find values at midpoint
      array_index = (R + L)/2
      if (val >= array(array_index)) then
        L = array_index + 1
      else
        R = array_index
      end if
    end do

    array_index = L - 1

  end function binary_search_int4

  function binary_search_int8(array, n, val) result(array_index)

    integer,    intent(in) :: n
    integer(8), intent(in) :: array(n)
    integer(8), intent(in) :: val
    integer                :: array_index

    integer :: L
    integer :: R

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
      message = "Value outside of array during binary search"
      call fatal_error()
    end if

    do while (R > L)
      ! Find values at midpoint
      array_index = L + (R + L)/2
      testval = 
      if (val > array(array_index)) then
        L = array_index + 1
      else
        R = array_index
      end if
    end do

    array_index = L - 1

  end function binary_search_int8

end module search
