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
    integer :: i
    ! Note that we could use s = transfer(vs, s)
    ! here and (in principle) this could be an O(1)
    ! operation. However, this leakes memory on all
    ! recent gfortrans and crashes on some versions of
    ! PGF90. So the loop would need to be included with 
    ! an #ifdef PGF90. In any case, the looping version
    ! seems to be quickers (probably due to the character
    ! copying being vectorized in the loop). See:
    ! http://gcc.gnu.org/bugzilla/show_bug.cgi?id=51175
    do i = 1, size(vs)
      s(i:i) = vs(i)
    enddo
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
