module algorithm

  use constants
  use stl_vector, only: VectorInt, VectorReal

  implicit none

  interface sort
    module procedure sort_int, sort_real, sort_vector_int, sort_vector_real
  end interface sort

  interface find
    module procedure find_int, find_real, find_vector_int, find_vector_real
  end interface find

contains

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
