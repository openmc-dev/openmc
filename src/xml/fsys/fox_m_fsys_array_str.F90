module fox_m_fsys_array_str
#ifndef DUMMYLIB

  implicit none
  private

  interface destroy
    module procedure destroy_vs
  end interface

  interface concat
    module procedure vs_s_concat
  end interface

  interface alloc
    module procedure vs_str_alloc
  end interface

  public :: str_vs
  public :: vs_str
  public :: vs_str_alloc
  public :: vs_vs_alloc

  public :: concat
  public :: alloc

  public :: destroy

contains

  pure function str_vs(vs) result(s)
    character, dimension(:), intent(in) :: vs
    character(len=size(vs)) :: s
#ifdef PGF90
!PGI crashes on this use of transfer. Knob-ends.
    integer :: i
    do i = 1, size(vs)
      s(i:i) = vs(i)
    enddo
#else
    s = transfer(vs, s)
#endif
  end function str_vs


  pure function vs_str(s) result(vs)
    character(len=*), intent(in) :: s
    character, dimension(len(s)) :: vs

    vs = transfer(s, vs)
  end function vs_str

  pure function vs_str_alloc(s) result(vs)
    character(len=*), intent(in) :: s
    character, dimension(:), pointer :: vs

    allocate(vs(len(s)))
    vs = vs_str(s)
  end function vs_str_alloc

  pure function vs_vs_alloc(s) result(vs)
    character, dimension(:), pointer :: s
    character, dimension(:), pointer :: vs

    if (associated(s)) then
      allocate(vs(size(s)))
      vs = s
    else
      vs => null()
    endif
  end function vs_vs_alloc

  pure function vs_s_concat(vs, s) result(vs2)
    character, dimension(:), intent(in) :: vs
    character(len=*), intent(in) :: s
    character, dimension(:), pointer :: vs2

    allocate(vs2(size(vs)+len(s)))
    vs2(:size(vs)) = vs
    vs2(size(vs)+1:) = vs_str(s)
  end function vs_s_concat

  subroutine destroy_vs(vs)
    character, dimension(:), pointer :: vs

    deallocate(vs)
  end subroutine destroy_vs

#endif
end module fox_m_fsys_array_str
