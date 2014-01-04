module fox_m_fsys_string

#ifndef DUMMYLIB
  ! Assorted generally useful string manipulation functions

  implicit none
  private

  character(len=26), parameter :: lowerAlphabet = &
    "abcdefghijklmnopqrstuvwxyz"
  character(len=26), parameter :: UPPERAlphabet = &
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

  public :: toLower

contains

  function toLower(s_in) result(s)
    character(len=*), intent(in) :: s_in
    character(len=len(s_in)) :: s

    integer :: i, n
    do i = 1, len(s)
      n = index(UPPERAlphabet, s_in(i:i))
      if (n>0) then
        s(i:i) = lowerAlphabet(n:n)
      else
        s(i:i) = s_in(i:i)
      endif
    enddo

  end function toLower

#endif
end module fox_m_fsys_string
