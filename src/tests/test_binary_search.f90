program test_binary_search

  use search, only: binary_search

  implicit none

  integer :: i
  integer :: index
  integer, parameter :: n = 100
  real(8) :: array(n)
  real(8) :: lbound
  real(8) :: ubound
  real(8) :: val

  array = (/ (log(real(i+2)), i = 1, n) /)

  val = 4.134
  index = binary_search(array, n, val)

  write(UNIT=*, FMT='(A,I2,A,ES11.4)') "array(", index, ") = ", array(index)
  write(UNIT=*, FMT='(A,ES11.4)') "value     = ", val
  write(UNIT=*, FMT='(A,I2,A,ES11.4)') "array(", index + 1, ") = ", array(index + 1)

end program test_binary_search
