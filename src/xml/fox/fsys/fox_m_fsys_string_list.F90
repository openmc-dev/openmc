module fox_m_fsys_string_list
#ifndef DUMMYLIB

  use fox_m_fsys_array_str, only: str_vs, vs_str_alloc
  implicit none
  private

  type string_t
    character, pointer :: s(:) => null()
  end type string_t

  type string_list
    type(string_t), pointer :: list(:) => null()
  end type string_list

  public :: string_t
  public :: string_list

  public :: init_string_list
  public :: destroy_string_list
  public :: add_string
  public :: remove_last_string
  public :: get_last_string
  public :: tokenize_to_string_list
  public :: tokenize_and_add_strings
  public :: registered_string

  interface destroy
    module procedure destroy_string_list
  end interface

  public :: destroy

contains

  subroutine init_string_list(s_list)
    type(string_list), intent(inout) :: s_list

    allocate(s_list%list(0))
  end subroutine init_string_list

  subroutine destroy_string_list(s_list)
    type(string_list), intent(inout) :: s_list

    integer :: i

    if (associated(s_list%list)) then
      do i = 1, ubound(s_list%list, 1)
        deallocate(s_list%list(i)%s)
      enddo
      deallocate(s_list%list)
    endif
  end subroutine destroy_string_list

  subroutine add_string(s_list, s)
    type(string_list), intent(inout) :: s_list
    character(len=*), intent(in) :: s

    integer :: i
    type(string_t), pointer :: temp(:)

    temp => s_list%list
    allocate(s_list%list(size(temp)+1))
    do i = 1, size(temp)
      s_list%list(i)%s => temp(i)%s
    enddo
    deallocate(temp)
    s_list%list(i)%s => vs_str_alloc(s)
  end subroutine add_string

  subroutine remove_last_string(s_list)
    type(string_list), intent(inout) :: s_list

    integer :: i
    type(string_t), pointer :: temp(:)

    temp => s_list%list
    allocate(s_list%list(size(temp)-1))
    do i = 1, size(temp)-1
      s_list%list(i)%s => temp(i)%s
    enddo
    deallocate(temp)

  end subroutine remove_last_string

  function get_last_string(s_list) result(s)
    type(string_list), intent(in) :: s_list
    character(len=size(s_list%list(size(s_list%list))%s)) :: s

    s = str_vs(s_list%list(size(s_list%list))%s)
  end function get_last_string

  function tokenize_to_string_list(s) result(s_list)
    character(len=*), intent(in) :: s
    type(string_list) :: s_list

    ! tokenize a whitespace-separated list of strings
    ! and place results in a string list

    character(len=*), parameter :: &
      WHITESPACE = achar(9)//achar(10)//achar(13)//achar(32)
    integer :: i, j

    call init_string_list(s_list)

    i = verify(s, WHITESPACE)
    if (i==0) return
    j = scan(s(i:), WHITESPACE)
    if (j==0) then
      j = len(s)
    else
      j = i + j - 2
    endif
    do
      call add_string(s_list, s(i:j))
      i = j + 1
      j = verify(s(i:), WHITESPACE)
      if (j==0) exit
      i = i + j - 1
      j = scan(s(i:), WHITESPACE)
      if (j==0) then
        j = len(s)
      else
        j = i + j - 2
      endif
    enddo

  end function tokenize_to_string_list

  function registered_string(s_list, s) result(p)
    type(string_list), intent(in) :: s_list
    character(len=*), intent(in) :: s
    logical :: p

    integer :: i

    p = .false.
    do i = 1, size(s_list%list)
      if (str_vs(s_list%list(i)%s)//"x"==s//"x") then
        p = .true.
        exit
      endif
    enddo
  end function registered_string

  subroutine tokenize_and_add_strings(s_list, s, uniquify)
    type(string_list), intent(inout) :: s_list
    character(len=*), intent(in) :: s
    logical, intent(in), optional :: uniquify

    ! tokenize a whitespace-separated list of strings
    ! and place results in the given string list

    character(len=*), parameter :: & 
      WHITESPACE = achar(9)//achar(10)//achar(13)//achar(32)
    integer :: i, j
    logical :: uniquify_
    
    if (present(uniquify)) then
      uniquify_ = uniquify
    else
      uniquify_ = .false.
    endif

    i = verify(s, WHITESPACE)
    if (i==0) return
    j = scan(s(i:), WHITESPACE)
    if (j==0) then
      j = len(s)
    else
      j = i + j - 2
    endif
    do
      if (uniquify_.and..not.registered_string(s_list, s(i:j))) &
        call add_string(s_list, s(i:j))
      i = j + 1
      j = verify(s(i:), WHITESPACE)
      if (j==0) exit
      i = i + j - 1
      j = scan(s(i:), WHITESPACE)
      if (j==0) then
        j = len(s)
      else
        j = i + j - 2
      endif
    enddo

  end subroutine tokenize_and_add_strings

#endif
end module fox_m_fsys_string_list
