module fox_m_fsys_varstr
  implicit none

  private

  public :: varstr
  public :: init_varstr
  public :: destroy_varstr
  public :: varstr_len
  public :: is_varstr_empty
  public :: set_varstr_empty
  public :: is_varstr_null
  public :: set_varstr_null
  public :: vs_varstr_alloc
  public :: move_varstr_vs
  public :: move_varstr_varstr
  public :: str_varstr
  public :: append_varstr
  public :: varstr_str
  public :: varstr_vs
  public :: equal_varstr_str
  public :: equal_varstr_varstr

  ! Allocation step in which the data within varstr will be allocated
  integer, parameter :: VARSTR_INIT_SIZE=1024
  integer, parameter :: VARSTR_ALLOC_SIZE=1024

  ! Variable size string type.
  ! It is used for token field, so that it can be grown
  ! by one character at a time without incurring constant
  ! penalty on allocating/deallocating vs data.
  ! See functions and subroutines at the end of this module.
  type varstr
    private
    character, dimension(:), pointer :: data
    integer :: length
  end type varstr

contains

! Initialise varstr type. The string is initialised as null (i.e. invalid)
subroutine init_varstr(vstr)
  type(varstr), intent(inout) :: vstr
  allocate(vstr%data(VARSTR_INIT_SIZE))
  vstr%length = -1
end subroutine init_varstr

! Clean up memory (leaves varstr null, and drops the data field)
subroutine destroy_varstr(vstr)
  type(varstr), intent(inout) :: vstr
  if (associated(vstr%data)) deallocate(vstr%data)
  call set_varstr_null(vstr)
end subroutine destroy_varstr

! Return real length of varstr
function varstr_len(vstr) result(l)
  type(varstr), intent(in) :: vstr
  integer :: l
  if (vstr%length<0) print *, "WARNING: asking for length of null varstr"
  l = vstr%length
end function varstr_len


! Make sure that varstr is at least size-n.
! The data will be kept (copied) if keep is true (default)
! This can be called on a null varstr, but should noty be called on
! one which was destroyed.
subroutine ensure_varstr_size(vstr,n,keep)
  type(varstr), intent(inout) :: vstr
  integer, intent(in) :: n
  logical, optional, intent(in) :: keep

  character, pointer, dimension(:) :: new_data
  integer :: new_size, old_size
  logical :: keep_flag

  if (present(keep)) then
    keep_flag = keep
  else
    keep_flag = .true.
  end if

  old_size = size(vstr%data)
  if (n <= old_size ) return

  new_size = old_size + ((n-old_size)/VARSTR_ALLOC_SIZE+1) * VARSTR_ALLOC_SIZE
  allocate(new_data(new_size))

  if (keep_flag) new_data(1:old_size) = vstr%data(1:old_size)

  deallocate( vstr%data )
  vstr%data => new_data
end subroutine ensure_varstr_size

! Returns whether varstr is empty: ""
function is_varstr_empty(vstr)
  type(varstr), intent(in) :: vstr
  logical is_varstr_empty
  is_varstr_empty = (vstr%length == 0)
end function is_varstr_empty

! Set vstr to empty string
subroutine set_varstr_empty(vstr)
  type(varstr), intent(inout) :: vstr
  vstr%length = 0
end subroutine set_varstr_empty

! Returns whether varstr is null (i.e. invalid)
function is_varstr_null(vstr)
  type(varstr), intent(in) :: vstr
  logical is_varstr_null
  is_varstr_null = (vstr%length < 0)
end function is_varstr_null

! Set vstr to null
subroutine set_varstr_null(vstr)
  type(varstr), intent(inout) :: vstr
  vstr%length = -1
end subroutine set_varstr_null

! Convert varstr to newly allocated array of characters
function vs_varstr_alloc(vstr) result(vs)
  type(varstr) :: vstr
  character, dimension(:), pointer :: vs

  if (is_varstr_null(vstr)) then
    print *, "WARNING: Converting null varstr to string... making it empty first"
    call set_varstr_empty(vstr)
  end if

  allocate(vs(vstr%length))
  vs = vstr%data(1:vstr%length)
end function vs_varstr_alloc

! This call moves data from varstr to vs (i.e. vs is overwritten and vstr is made null)
subroutine move_varstr_vs(vstr,vs)
  type(varstr), intent(inout) :: vstr
  character, dimension(:), pointer, intent(inout) :: vs

  if (associated(vs)) deallocate(vs)
  vs => vs_varstr_alloc(vstr)
  call set_varstr_null(vstr)
end subroutine move_varstr_vs

! This call moves data from varstr to varstr (src becomes null)
subroutine move_varstr_varstr(src,dst)
  type(varstr), intent(inout) :: src
  type(varstr), intent(inout) :: dst
  character, dimension(:), pointer :: tmpdata

  tmpdata => dst%data
  dst%data => src%data
  src%data => tmpdata
  dst%length = src%length

  call set_varstr_null(src)
end subroutine move_varstr_varstr


! Convert varstr to string type
function str_varstr(vstr) result(s)
  type(varstr), intent(in) :: vstr
  character(len=vstr%length) :: s
  integer :: i

  if (is_varstr_null(vstr)) then
    ! Can we really end-up here? Or will it blow on allocation with len=-1 ?
    print *, "WARNING: Trying to convert null varstr to str... returning empty string"
    s = ""
  end if

  do i = 1, vstr%length
    s(i:i) = vstr%data(i)
  enddo
end function str_varstr

! Append string to varstr
subroutine append_varstr(vstr,str)
  type(varstr), intent(inout) :: vstr
  character(len=*), intent(in) :: str
  character, dimension(:), pointer :: tmp
  integer :: i

  if (is_varstr_null(vstr)) then
    print *, "WARNING: Trying to append to null varstr... making it empty first"
    call set_varstr_empty(vstr)
  end if

  call ensure_varstr_size(vstr,vstr%length+len(str))

  ! Note: on a XML file with very large tokens, this loop
  ! is consistently faster than equivalent 'transfer' intrinsic
  do i=1,len(str)
    vstr%data(vstr%length+i) = str(i:i)
  end do
  vstr%length = vstr%length + len(str)
end subroutine append_varstr

! Convert string to varstr in place
subroutine varstr_str(vstr,str)
  type(varstr), intent(inout) :: vstr
  character(len=*), intent(in) :: str
  integer :: i

  call ensure_varstr_size(vstr,len(str),.false.)
  do i=1,len(str)
    vstr%data(i) = str(i:i)
  end do
  vstr%length = len(str)
end subroutine varstr_str

! Convert character array to varstr in place
subroutine varstr_vs(vstr,vs)
  type(varstr), intent(inout) :: vstr
  character, dimension(:), intent(in) :: vs

  call ensure_varstr_size(vstr,size(vs),.false.)

  vstr%length = size(vs)
  vstr%data(1:size(vs)) = vs
end subroutine varstr_vs

! Compare varstr to str (true if equal)
function equal_varstr_str(vstr,str) result(r)
  type(varstr), intent(in) :: vstr
  character(len=*) :: str
  logical :: r

  integer :: i

  r = .false.
  if ( len(str) /= varstr_len(vstr) ) return
  do i=1,len(str)
    if ( str(i:i) /= vstr%data(i) ) return
  end do
  r = .true.
end function equal_varstr_str

! Compare varstr to varstr (true if equal)
function equal_varstr_varstr(vstr1,vstr2) result(r)
  type(varstr), intent(in) :: vstr1, vstr2
  logical :: r

  integer :: i

  r = .false.
  if ( varstr_len(vstr1) /= varstr_len(vstr2) ) return
  do i=1,varstr_len(vstr1)
    if ( vstr1%data(i) /= vstr2%data(i) ) return
  end do
  r = .true.
end function equal_varstr_varstr

end module fox_m_fsys_varstr
