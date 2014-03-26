module m_common_entities

#ifndef DUMMYLIB

  use fox_m_fsys_array_str, only: str_vs, vs_str_alloc
  use fox_m_fsys_format, only: str_to_int_10, str_to_int_16
  use fox_m_utils_uri, only: URI, destroyURI
  use m_common_charset, only: digits, hexdigits
  use m_common_error, only: FoX_error

  implicit none
  private

  type entity_t
    logical :: external
    logical :: wfc ! Was this entity declared externally or in a PE, where
                   ! a non-validating processor might not see it?
    character(len=1), dimension(:), pointer :: name => null()
    character(len=1), dimension(:), pointer :: text => null()
    character(len=1), dimension(:), pointer :: publicId => null()
    character(len=1), dimension(:), pointer :: systemId => null()
    character(len=1), dimension(:), pointer :: notation => null()
    type(URI), pointer :: baseURI => null()
  end type entity_t

  type entity_list
    private
    type(entity_t), dimension(:), pointer :: list => null()
  end type entity_list

  public :: is_unparsed_entity
  public :: is_external_entity

  public :: expand_entity_text
  public :: expand_entity_text_len
  public :: existing_entity

  public :: expand_char_entity

  public :: expand_entity
  public :: expand_entity_len

  public :: entity_t
  public :: entity_list
  public :: init_entity_list
  public :: reset_entity_list
  public :: destroy_entity_list
  public :: print_entity_list
  public :: add_internal_entity
  public :: add_external_entity
  public :: pop_entity_list

  interface size
    module procedure size_el
  end interface

  interface is_unparsed_entity
    module procedure is_unparsed_entity_
    module procedure is_unparsed_entity_from_list
  end interface

  public :: getEntityByIndex
  public :: getEntityByName

  public :: size

contains

  function size_el(el) result(n)
    type(entity_list), intent(in) :: el
    integer :: n

    n = ubound(el%list, 1)
  end function size_el

  function shallow_copy_entity(ent1) result(ent2)
    type(entity_t), intent(in) :: ent1
    type(entity_t) :: ent2
    
    ent2%external = ent1%external
    ent2%wfc = ent1%wfc
    ent2%name => ent1%name
    ent2%text => ent1%text
    ent2%publicId => ent1%publicId
    ent2%systemId => ent1%systemId
    ent2%notation => ent1%notation
    ent2%baseURI => ent1%baseURI

  end function shallow_copy_entity

  function getEntityByIndex(el, i) result(e)
    type(entity_list), intent(in) :: el
    integer, intent(in) :: i
    type(entity_t), pointer :: e

    e => el%list(i)
  end function getEntityByIndex

  function getEntityNameByIndex(el, i) result(c)
    type(entity_list), intent(in) :: el
    integer, intent(in) :: i
    character(len=size(el%list(i)%name)) :: c

    c = str_vs(el%list(i)%name)
  end function getEntityNameByIndex

  function getEntityByName(el, name) result(e)
    type(entity_list), intent(in) :: el
    character(len=*), intent(in) :: name
    type(entity_t), pointer :: e

    integer :: i

    e => null()
    do i = 1, size(el%list)
      if (str_vs(el%list(i)%name)==name) then
        e => el%list(i)
        exit
      endif
    enddo
  end function getEntityByName


  subroutine destroy_entity(ent)
    type(entity_t), intent(inout) :: ent
    
    deallocate(ent%name)
    deallocate(ent%text)
    deallocate(ent%publicId)
    deallocate(ent%systemId)
    deallocate(ent%notation)

    if (associated(ent%baseURI)) call destroyURI(ent%baseURI)

  end subroutine destroy_entity


  subroutine init_entity_list(ents)
    type(entity_list), intent(inout) :: ents

    if (associated(ents%list)) deallocate(ents%list)
    allocate(ents%list(0))

  end subroutine init_entity_list


  subroutine reset_entity_list(ents)
    type(entity_list), intent(inout) :: ents

    call destroy_entity_list(ents)
    call init_entity_list(ents)

  end subroutine reset_entity_list


  subroutine destroy_entity_list(ents)
    type(entity_list), intent(inout) :: ents

    integer :: i, n

    n = size(ents%list)
    do i = 1, n
      call destroy_entity(ents%list(i))
    enddo
    deallocate(ents%list)
  end subroutine destroy_entity_list

  function pop_entity_list(ents) result(name)
    type(entity_list), intent(inout) :: ents
    character(len=size(ents%list(size(ents%list))%name)) :: name
    
    type(entity_t), pointer :: ents_tmp(:)
    integer :: i, n
    n = size(ents%list)

    ents_tmp => ents%list
    allocate(ents%list(n-1))
    do i = 1, n - 1
      ents%list(i) = shallow_copy_entity(ents_tmp(i))
    enddo
    name = str_vs(ents_tmp(i)%name)

    call destroy_entity(ents_tmp(i))
    deallocate(ents_tmp)
  end function pop_entity_list

  subroutine print_entity_list(ents)
    type(entity_list), intent(in) :: ents

    integer :: i, n

    n = size(ents%list)
    write(*,'(a)') '>ENTITYLIST'
    do i = 1, n
      write(*,'(a)') str_vs(ents%list(i)%name)
      write(*,'(a)') str_vs(ents%list(i)%text)
      write(*,'(a)') str_vs(ents%list(i)%publicId)
      write(*,'(a)') str_vs(ents%list(i)%systemId)
      write(*,'(a)') str_vs(ents%list(i)%notation)
    enddo
    write(*,'(a)') '<ENTITYLIST'
  end subroutine print_entity_list


  subroutine add_entity(ents, name, text, publicId, systemId, notation, baseURI, wfc)
    type(entity_list), intent(inout) :: ents
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: text
    character(len=*), intent(in) :: publicId
    character(len=*), intent(in) :: systemId
    character(len=*), intent(in) :: notation
    type(URI), pointer :: baseURI
    logical, intent(in) :: wfc

    type(entity_t), pointer :: ents_tmp(:)
    integer :: i, n

    ! This should only ever be called by add_internal_entity or add_external_entity
    ! below, so we don't bother sanity-checking input. Note especially we don't 
    ! check for duplication of entities, so this will happily add another entity
    ! of the same name if you ask it to. This should't matter though, since the
    ! first defined will always be picked up first, which is what the XML spec
    ! requires.

    n = size(ents%list)
    ents_tmp => ents%list
    allocate(ents%list(n+1))
    do i = 1, n
      ents%list(i) = shallow_copy_entity(ents_tmp(i))
    enddo
    deallocate(ents_tmp)
    ents%list(i)%external = len(systemId)>0
    ents%list(i)%wfc = wfc
    ents%list(i)%name => vs_str_alloc(name)
    ents%list(i)%text => vs_str_alloc(text)
    ents%list(i)%publicId => vs_str_alloc(publicId)
    ents%list(i)%systemId => vs_str_alloc(systemId)
    ents%list(i)%notation => vs_str_alloc(notation)
    ents%list(i)%baseURI => baseURI
  end subroutine add_entity


  subroutine add_internal_entity(ents, name, text, baseURI, wfc)
    type(entity_list), intent(inout) :: ents
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: text
    type(URI), pointer :: baseURI
    logical, intent(in) :: wfc

    call add_entity(ents, name=name, text=text, &
      publicId="", systemId="", notation="", baseURI=baseURI, wfc=wfc)
  end subroutine add_internal_entity

  
  subroutine add_external_entity(ents, name, systemId, baseURI, wfc, publicId, notation)
    type(entity_list), intent(inout) :: ents
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: systemId
    character(len=*), intent(in), optional :: publicId
    character(len=*), intent(in), optional :: notation
    type(URI), pointer :: baseURI
    logical, intent(in) :: wfc

    if (present(publicId) .and. present(notation)) then
      call add_entity(ents, name=name, text="", &
        publicId=publicId, systemId=systemId, notation=notation, &
        wfc=wfc, baseURI=baseURI)
    elseif (present(publicId)) then
      call add_entity(ents, name=name, text="", &
        publicId=publicId, systemId=systemId, notation="", &
        wfc=wfc, baseURI=baseURI)
    elseif (present(notation)) then
      call add_entity(ents, name=name, text="", &
        publicId="", systemId=systemId, notation=notation, &
        wfc=wfc, baseURI=baseURI)
    else
      call add_entity(ents, name=name, text="", &
        publicId="", systemId=systemId, notation="", &
        wfc=wfc, baseURI=baseURI)
    endif
  end subroutine add_external_entity


  function is_unparsed_entity_from_list(ents, name) result(p)
    type(entity_list), intent(in) :: ents
    character(len=*), intent(in) :: name
    logical :: p

    integer :: i

    p = .false.

    do i = 1, size(ents%list)
      if (name == str_vs(ents%list(i)%name)) then
        p = (size(ents%list(i)%notation)>0)
        exit
      endif
    enddo
  end function is_unparsed_entity_from_list

  function is_unparsed_entity_(ent) result(p)
    type(entity_t), intent(in) :: ent
    logical :: p

    p = (size(ent%notation)>0)

  end function is_unparsed_entity_


  function is_external_entity(ents, name) result(p)
    type(entity_list), intent(in) :: ents
    character(len=*), intent(in) :: name
    logical :: p

    integer :: i

    p = .false.

    do i = 1, size(ents%list)
      if (name == str_vs(ents%list(i)%name)) then
        p = ents%list(i)%external
        exit
      endif
    enddo
  end function is_external_entity
 
  pure function expand_char_entity_len(name) result(n)
    character(len=*), intent(in) :: name
    integer :: n

    integer :: number

    if (name(1:1) == "#") then
      if (name(2:2) == "x") then       ! hex character reference
        if (verify(name(3:), hexdigits) == 0) then
          number = str_to_int_16(name(3:))   
          if (0 <= number .and. number <= 128) then
            n = 1
          else
            n = len(name) + 2
          endif
        else 
           n = 0
        endif
      else                             ! decimal character reference
        if (verify(name(3:), digits) == 0) then
          number = str_to_int_10(name(2:))
          if (0 <= number .and. number <= 128) then
            n = 1
          else
            n = len(name) + 2
          endif
        else 
          n = 0
        endif
      endif
    else
      n = 0
    endif
  end function expand_char_entity_len


  function expand_char_entity(name) result(text)
    character(len=*), intent(in) :: name
    character(len=expand_char_entity_len(name)) :: text

    integer :: number

    select case (len(text))
    case (0)
      call FoX_error("Invalid character entity reference")
    case (1)  
      if (name(2:2) == "x") then       ! hex character reference
        number = str_to_int_16(name(3:))   
      else                             ! decimal character reference
        number = str_to_int_10(name(2:))
      endif
      text = achar(number)
      ! FIXME what about > 127 ...
    case default
      text = "&"//name//";"
    end select

  end function expand_char_entity


  pure function existing_entity(ents, name) result(p)
    type(entity_list), intent(in) :: ents
    character(len=*), intent(in)  :: name
    logical :: p

    integer :: i

    p = .false.
    
!FIXME the following test is not entirely in accordance with the valid chars check we do elsewhere...

    do i = 1, size(ents%list)
      if (name == str_vs(ents%list(i)%name)) then
        p = .true.
        return
      endif
    enddo

  end function existing_entity


  pure function expand_entity_text_len(ents, name) result(n)
    type(entity_list), intent(in) :: ents
    character(len=*), intent(in)  :: name
    integer :: n

    integer :: i

    do i = 1, size(ents%list)
      if (name == str_vs(ents%list(i)%name)) then
        n = size(ents%list(i)%text)
      endif
    enddo

  end function expand_entity_text_len


  function expand_entity_text(ents, name) result(text)
    type(entity_list), intent(in) :: ents
    character(len=*), intent(in)  :: name
    character(len=expand_entity_text_len(ents, name)) :: text

    integer :: i

    ! No error checking - make sure entity exists before calling it.

    do i = 1, size(ents%list)
      if (name == str_vs(ents%list(i)%name)) then
        text = str_vs(ents%list(i)%text)
        exit
      endif
    enddo

  end function expand_entity_text


  pure function expand_entity_len(ents, name) result(n)
    type(entity_list), intent(in) :: ents
    character(len=*), intent(in)  :: name
    integer :: n

    integer :: i

    do i = 1, size(ents%list)
      if (name == str_vs(ents%list(i)%name)) then
        n = size(ents%list(i)%text)
      endif
    enddo

  end function expand_entity_len


  function expand_entity(ents, name) result(text)
    type(entity_list), intent(in) :: ents
    character(len=*), intent(in)  :: name
    character(len=expand_entity_len(ents, name)) :: text
    
    integer :: i
    
    do i = 1, size(ents%list)
      if (name == str_vs(ents%list(i)%name)) then
        text = str_vs(ents%list(i)%text)
      endif
    enddo

  end function expand_entity

#endif
end module m_common_entities
