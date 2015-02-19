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
! NEAR_NEIGHB finds the nearest-neighbor index given bounds and a value
!===============================================================================

  function near_neighb(i_low, lower, upper, val) result(i_near)

    real(8), intent(in) :: lower  ! lower bound value
    real(8), intent(in) :: upper  ! upper bound value
    real(8), intent(in) :: val    ! value
    integer             :: i_low  ! lower bound index
    integer             :: i_near ! nearest-neighbor index
    real(8)             :: dlow   ! distance to lower bound
    real(8)             :: dup    ! distance to upper bound

    dlow = val - lower
    dup  = upper - val

    i_near = i_low

    if (dup < dlow) i_near = i_near + 1

  end function near_neighb

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
      write(*,'(ES16.9,ES16.9,ES16.9)') val, array(L), array(R)
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
      write(*,'(ES16.9,ES16.9,ES16.9)') val, array(L), array(R)
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
      write(*,'(ES16.9,ES16.9,ES16.9)') val, array(L), array(R)
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
