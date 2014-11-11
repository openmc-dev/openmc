module search

  use constants
  use error,     only: fatal_error

  implicit none

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
    integer :: n_iteration
    real(8) :: testval

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
      call fatal_error("Value outside of array during binary search")
    end if

    n_iteration = 0
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
      if (val >= testval) then
        L = array_index
      elseif (val < testval) then
        R = array_index
      end if

      ! check for large number of iterations
      n_iteration = n_iteration + 1
      if (n_iteration == MAX_ITERATION) then
        call fatal_error("Reached maximum number of iterations on binary &
            &search.")
      end if
    end do

    array_index = L

  end function binary_search_real

  function binary_search_int4(array, n, val) result(array_index)

    integer, intent(in) :: n
    integer, intent(in) :: array(n)
    integer, intent(in) :: val
    integer             :: array_index

    integer :: L
    integer :: R
    integer :: n_iteration
    real(8) :: testval

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
      call fatal_error("Value outside of array during binary search")
    end if

    n_iteration = 0
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
      if (val >= testval) then
        L = array_index
      elseif (val < testval) then
        R = array_index
      end if

      ! check for large number of iterations
      n_iteration = n_iteration + 1
      if (n_iteration == MAX_ITERATION) then
        call fatal_error("Reached maximum number of iterations on binary &
            &search.")
      end if
    end do

    array_index = L

  end function binary_search_int4

  function binary_search_int8(array, n, val) result(array_index)

    integer,    intent(in) :: n
    integer(8), intent(in) :: array(n)
    integer(8), intent(in) :: val
    integer                :: array_index

    integer :: L
    integer :: R
    integer :: n_iteration
    real(8) :: testval

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
      call fatal_error("Value outside of array during binary search")
    end if

    n_iteration = 0
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
      if (val >= testval) then
        L = array_index
      elseif (val < testval) then
        R = array_index
      end if

      ! check for large number of iterations
      n_iteration = n_iteration + 1
      if (n_iteration == MAX_ITERATION) then
        call fatal_error("Reached maximum number of iterations on binary &
            &search.")
      end if
    end do

    array_index = L

  end function binary_search_int8

end module search
