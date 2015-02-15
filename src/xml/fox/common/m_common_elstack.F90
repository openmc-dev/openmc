module m_common_elstack

#ifndef DUMMYLIB
  use fox_m_fsys_array_str, only: str_vs, vs_str
  use m_common_error, only: FoX_fatal
  use m_common_content_model, only: content_particle_t, checkCP, &
    elementContentCP, emptyContentCP, checkCPToEnd

  implicit none
  private

  ! Element stack during parsing. Keeps track of element names
  ! and optionally tracks validity of content model

  ! Initial stack size:
  integer, parameter :: STACK_SIZE_INIT = 10
  ! Multiplier when stack is exceeded:
  real, parameter :: STACK_SIZE_MULT = 1.5

  type :: elstack_item
    character, dimension(:), pointer          :: name => null()
    type(content_particle_t), pointer         :: cp => null()
  end type elstack_item

  type :: elstack_t
    private
    integer                                   :: n_items
    type(elstack_item), pointer, dimension(:) :: stack => null()
  end type elstack_t

  public :: elstack_t

  public  :: push_elstack, pop_elstack, init_elstack, destroy_elstack, reset_elstack, print_elstack
  public  :: get_top_elstack, is_empty
  public :: checkContentModel
  public :: checkContentModelToEnd
  public :: elementContent
  public :: emptyContent
  public  :: len

  interface len
    module procedure number_of_items
  end interface

  interface is_empty
    module procedure is_empty_elstack
  end interface

contains

  subroutine init_elstack(elstack)
    type(elstack_t), intent(inout)  :: elstack

    ! We go from 0 (and initialize the 0th string to "")
    ! in order that we can safely check the top of an
    ! empty stack
    allocate(elstack%stack(0:STACK_SIZE_INIT))
    elstack%n_items = 0
    allocate(elstack%stack(0)%name(0))

  end subroutine init_elstack

  subroutine destroy_elstack(elstack)
    type(elstack_t), intent(inout)  :: elstack
    integer :: i
    do i = 0, elstack % n_items
      deallocate(elstack%stack(i)%name)
    enddo
    deallocate(elstack%stack)
  end subroutine destroy_elstack

  subroutine reset_elstack(elstack)
    type(elstack_t), intent(inout)  :: elstack

    call destroy_elstack(elstack)
    call init_elstack(elstack)

  end subroutine reset_elstack

  subroutine resize_elstack(elstack)
    type(elstack_t), intent(inout)  :: elstack
    type(elstack_item), dimension(0:ubound(elstack%stack,1)) :: temp
    integer :: i, s

    s = ubound(elstack%stack, 1)

    do i = 0, s
      temp(i)%name => elstack%stack(i)%name
      temp(i)%cp => elstack%stack(i)%cp
    enddo
    deallocate(elstack%stack)
    allocate(elstack%stack(0:nint(s*STACK_SIZE_MULT)))
    do i = 0, s
      elstack%stack(i)%name => temp(i)%name
      elstack%stack(i)%cp => temp(i)%cp
    enddo

  end subroutine resize_elstack

  pure function is_empty_elstack(elstack) result(answer)
    type(elstack_t), intent(in)  :: elstack
    logical                    :: answer

    answer = (elstack%n_items == 0)
  end function is_empty_elstack

  function number_of_items(elstack) result(n)
    type(elstack_t), intent(in)  :: elstack
    integer                      :: n

    n = elstack%n_items
  end function number_of_items

  subroutine push_elstack(elstack, name, cp)
    type(elstack_t), intent(inout)              :: elstack
    character(len=*), intent(in)                :: name
    type(content_particle_t), pointer, optional :: cp

    integer :: n

    n = elstack%n_items
    n = n + 1
    if (n == size(elstack%stack)) then
      call resize_elstack(elstack)
    endif
    allocate(elstack%stack(n)%name(len(name)))
    elstack%stack(n)%name = vs_str(name)
    if (present(cp)) elstack%stack(n)%cp => cp
    elstack%n_items = n

  end subroutine push_elstack

  function pop_elstack(elstack) result(item)
    type(elstack_t), intent(inout)     :: elstack
    character(len=merge(size(elstack%stack(elstack%n_items)%name), 0, elstack%n_items > 0)) :: item

    integer :: n

    n = elstack%n_items
    if (n == 0) then
      call FoX_fatal("Element stack empty")
    endif
    item = str_vs(elstack%stack(n)%name)
    deallocate(elstack%stack(n)%name)
    elstack%n_items = n - 1

  end function pop_elstack

  pure function get_top_elstack(elstack) result(item)
    ! Get the top element of the stack, *without popping it*.
    type(elstack_t), intent(in)        :: elstack
    character(len=merge(size(elstack%stack(elstack%n_items)%name), 0, elstack%n_items > 0)) :: item 

    integer :: n

    n = elstack%n_items

    if (n==0) then
      item = ""
    else
      item = str_vs(elstack%stack(n)%name)
    endif

  end function get_top_elstack

  function checkContentModel(elstack, name) result(p)
    type(elstack_t), intent(inout) :: elstack
    character(len=*), intent(in) :: name
    logical :: p

    type(content_particle_t), pointer :: cp

    integer :: n
    n = elstack%n_items
    if (n==0) then
      p = .true.
    else
      cp => elstack%stack(n)%cp
      p = checkCP(cp, name)
      elstack%stack(n)%cp => cp
    endif
  end function checkContentModel

  function checkContentModelToEnd(elstack) result(p)
    type(elstack_t), intent(inout) :: elstack
    logical :: p

    type(content_particle_t), pointer :: cp

    integer :: n
    n = elstack%n_items

    cp => elstack%stack(n)%cp
    p = checkCPToEnd(cp)

  end function checkContentModelToEnd

  function elementContent(elstack) result(p)
    type(elstack_t), intent(in) :: elstack
    logical :: p

    integer :: n
    n = elstack%n_items
    if (n==0) then
      p = .false.
    else
      p = elementContentCP(elstack%stack(n)%cp)
    endif
  end function elementContent

  function emptyContent(elstack) result(p)
    type(elstack_t), intent(in) :: elstack
    logical :: p

    integer :: n
    n = elstack%n_items
    if (n==0) then
      p = .false.
    else
      p = emptyContentCP(elstack%stack(n)%cp)
    endif
  end function emptyContent

  subroutine print_elstack(elstack,unit)
    type(elstack_t), intent(in)   :: elstack
    integer, intent(in)           :: unit
    integer   :: i

    do i = elstack%n_items, 1, -1
      write(unit=unit,fmt=*) elstack%stack(i)%name
    enddo

  end subroutine print_elstack

#endif
end module m_common_elstack
