module algorithm

  use constants
  use stl_vector, only: VectorInt, VectorReal

  implicit none

  integer, parameter :: MAX_ITERATION = 64

  interface binary_search
    module procedure binary_search_real, binary_search_int4, binary_search_int8
  end interface binary_search

  interface sort
    module procedure sort_int, sort_real, sort_vector_int, sort_vector_real
  end interface sort

  interface find
    module procedure find_int, find_real, find_vector_int, find_vector_real
  end interface find

contains

!===============================================================================
! BINARY_SEARCH performs a binary search of an array to find where a specific
! value lies in the array. This is used extensively for energy grid searching
!===============================================================================

  pure function binary_search_real(array, n, val) result(array_index)

    integer, intent(in) :: n
    real(8), intent(in) :: array(n)
    real(8), intent(in) :: val
    integer             :: array_index

    integer :: L
    integer :: R
    integer :: n_iteration

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
      array_index = -1
      return
    end if

    n_iteration = 0
    do while (R - L > 1)
      ! Find values at midpoint
      array_index = L + (R - L)/2
      if (val >= array(array_index)) then
        L = array_index
      else
        R = array_index
      end if

      ! check for large number of iterations
      n_iteration = n_iteration + 1
      if (n_iteration == MAX_ITERATION) then
        array_index = -2
        return
      end if
    end do

    array_index = L

  end function binary_search_real

  pure function binary_search_int4(array, n, val) result(array_index)

    integer, intent(in) :: n
    integer, intent(in) :: array(n)
    integer, intent(in) :: val
    integer             :: array_index

    integer :: L
    integer :: R
    integer :: n_iteration

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
      array_index = -1
      return
    end if

    n_iteration = 0
    do while (R - L > 1)
      ! Find values at midpoint
      array_index = L + (R - L)/2
      if (val >= array(array_index)) then
        L = array_index
      else
        R = array_index
      end if

      ! check for large number of iterations
      n_iteration = n_iteration + 1
      if (n_iteration == MAX_ITERATION) then
        array_index = -2
        return
      end if
    end do

    array_index = L

  end function binary_search_int4

  pure function binary_search_int8(array, n, val) result(array_index)

    integer,    intent(in) :: n
    integer(8), intent(in) :: array(n)
    integer(8), intent(in) :: val
    integer                :: array_index

    integer :: L
    integer :: R
    integer :: n_iteration

    L = 1
    R = n

    if (val < array(L) .or. val > array(R)) then
      array_index = -1
      return
    end if

    n_iteration = 0
    do while (R - L > 1)
      ! Find values at midpoint
      array_index = L + (R - L)/2
      if (val >= array(array_index)) then
        L = array_index
      else
        R = array_index
      end if

      ! check for large number of iterations
      n_iteration = n_iteration + 1
      if (n_iteration == MAX_ITERATION) then
        array_index = -2
        return
      end if
    end do

    array_index = L

  end function binary_search_int8

!===============================================================================
! SORT sorts an array in place using an insertion sort.
!===============================================================================

  pure subroutine sort_int(array)
    integer, intent(inout) :: array(:)

    integer :: k, m
    integer :: temp

    if (size(array) > 1) then
      SORT: do k = 2, size(array)
        ! Save value to move
        m = k
        temp = array(k)

        MOVE_OVER: do while (m > 1)
          ! Check if insertion value is greater than (m-1)th value
          if (temp >= array(m - 1)) exit

          ! Move values over until hitting one that's not larger
          array(m) = array(m - 1)
          m = m - 1
        end do MOVE_OVER

        ! Put the original value into its new position
        array(m) = temp
      end do SORT
    end if
  end subroutine sort_int

  pure subroutine sort_real(array)
    real(8), intent(inout) :: array(:)

    integer :: k, m
    real(8) :: temp

    if (size(array) > 1) then
      SORT: do k = 2, size(array)
        ! Save value to move
        m = k
        temp = array(k)

        MOVE_OVER: do while (m > 1)
          ! Check if insertion value is greater than (m-1)th value
          if (temp >= array(m - 1)) exit

          ! Move values over until hitting one that's not larger
          array(m) = array(m - 1)
          m = m - 1
        end do MOVE_OVER

        ! Put the original value into its new position
        array(m) = temp
      end do SORT
    end if
  end subroutine sort_real

  pure subroutine sort_vector_int(vec)
    type(VectorInt), intent(inout) :: vec

    call sort_int(vec % data(1:vec%size()))
  end subroutine sort_vector_int

  pure subroutine sort_vector_real(vec)
    type(VectorReal), intent(inout) :: vec

    call sort_real(vec % data(1:vec%size()))
  end subroutine sort_vector_real

!===============================================================================
! FIND determines the index of the first occurrence of a value in an array. If
! the value does not appear in the array, -1 is returned.
!===============================================================================

  pure function find_int(array, val) result(index)
    integer, intent(in) :: array(:)
    integer, intent(in) :: val
    integer             :: index

    integer :: i

    index = -1
    do i = 1, size(array)
      if (array(i) == val) then
        index = i
        exit
      end if
    end do
  end function find_int

  pure function find_real(array, val) result(index)
    real(8), intent(in) :: array(:)
    real(8), intent(in) :: val
    integer             :: index

    integer :: i

    index = -1
    do i = 1, size(array)
      if (array(i) == val) then
        index = i
        exit
      end if
    end do
  end function find_real

  pure function find_vector_int(vec, val) result(index)
    type(VectorInt), intent(in) :: vec
    integer,         intent(in) :: val
    integer                     :: index

    index = find_int(vec % data(1:vec % size()), val)
  end function find_vector_int

  pure function find_vector_real(vec, val) result(index)
    type(VectorReal), intent(in) :: vec
    real(8),          intent(in) :: val
    integer                      :: index

    index = find_real(vec % data(1:vec % size()), val)
  end function find_vector_real

end module algorithm
