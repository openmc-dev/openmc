module m_common_attrs

#ifndef DUMMYLIB
  use fox_m_fsys_array_str, only : str_vs, vs_str_alloc
  use m_common_element, only: get_att_type_enum, ATT_CDATA, ATT_CDAMB, ATT_TYPES
  use m_common_error, only : FoX_error, FoX_fatal

  implicit none
  private

  !Initial length of dictionary
  integer, parameter :: DICT_INIT_LEN = 10 
  !Multiplier if we need to extend it.
  real, parameter :: DICT_LEN_MULT = 1.5

  type dict_item
    character(len=1), pointer, dimension(:) :: nsURI => null()
    character(len=1), pointer, dimension(:) :: localName => null()
    character(len=1), pointer, dimension(:) :: prefix => null()
    character(len=1), pointer, dimension(:) :: key => null()
    character(len=1), pointer, dimension(:) :: value => null()
    logical :: specified = .true.
    logical :: declared = .false.
    logical :: isId = .false.
    integer :: type = 11
  end type dict_item

  type dict_item_ptr
    type(dict_item), pointer :: d => null()
  end type dict_item_ptr
#endif

  type dictionary_t
    private
#ifndef DUMMYLIB
    type(dict_item_ptr), dimension(:), pointer :: list => null()
    character, dimension(:), pointer           :: base => null()
#else
    integer :: i
#endif
  end type dictionary_t

  public :: dictionary_t

  ! Building procedures
#ifndef DUMMYLIB
  public :: init_dict
  public :: reset_dict
  public :: add_item_to_dict
  public :: destroy_dict
#endif
  ! Query and extraction procedures

  ! SAX names:
  public :: getIndex
  public :: getLength
  public :: getLocalName
  public :: getQName
  public :: getURI
  public :: getValue
  public :: getType
  public :: isSpecified
  public :: setSpecified
  public :: isDeclared
  public :: setDeclared
  public :: hasKey

#ifndef DUMMYLIB
  public :: len
  public :: get_key 
  public :: get_value
  public :: remove_key
  public :: has_key
  public :: print_dict

  ! Namespaces
  public :: get_prefix
  public :: get_localName
  public :: set_nsURI
  public :: set_prefix
  public :: set_localName
#endif

  ! For internal FoX use only:
  public :: get_att_index_pointer
  public :: getWhitespaceHandling
  public :: setIsId
  public :: getIsId

  public :: setBase
  public :: getBase

  public :: sortAttrs

  interface len
    module procedure getLength
  end interface

  interface hasKey
    module procedure has_key, has_key_ns
  end interface

  interface getIndex
    module procedure get_key_index, get_key_index_ns
  end interface

  interface getQName
    module procedure get_key
  end interface

  interface getValue
    module procedure get_value_by_key, get_value_by_index, get_value_by_key_ns
  end interface
#ifndef DUMMYLIB
  interface get_value
    module procedure get_value_by_key, get_value_by_index
  end interface
  interface remove_key
    module procedure remove_key_by_index
  end interface
#endif

  interface getURI
    module procedure get_nsURI_by_index
  end interface
#ifndef DUMMYLIB
  interface get_prefix
    module procedure get_prefix_by_index
  end interface
#endif
  interface getLocalName
    module procedure get_localName_by_index
  end interface
#ifndef DUMMYLIB
  interface get_localName
    module procedure get_localName_by_index
  end interface
  interface set_nsURI
    module procedure set_nsURI_by_index
  end interface
  interface set_prefix
    module procedure set_prefix_by_index
  end interface
  interface set_localName
    module procedure set_localName_by_index_s
    module procedure set_localName_by_index_vs
  end interface
#endif

  interface getType
    module procedure getType_by_index
    module procedure getType_by_keyname
  end interface

  interface isSpecified
    module procedure isSpecified_by_index
    module procedure isSpecified_by_key
    module procedure isSpecified_by_keyNS
  end interface

  interface isDeclared
    module procedure isDeclared_by_index
    module procedure isDeclared_by_key
    module procedure isDeclared_by_keyNS
  end interface

#ifndef DUMMYLIB
  interface getIsId
    module procedure getIsId_by_index
  end interface

  interface setIsId
    module procedure setIsId_by_index
  end interface

  interface destroy
    module procedure destroy_dict_item
    module procedure destroy_dict
  end interface

#endif

contains

  pure function getLength(dict) result(n)
    type(dictionary_t), intent(in) :: dict
    integer :: n

#ifndef DUMMYLIB
    n = ubound(dict%list, 1)
#else
    n = 0
#endif
  end function getLength


  function has_key(dict, key) result(found)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: key
    logical :: found

#ifndef DUMMYLIB
    integer :: i

    do i = 1, ubound(dict%list, 1)
      if (key==str_vs(dict%list(i)%d%key)) then
        found = .true.
        return
      endif
    enddo
    found = .false.
#else
    found = .false.
#endif
  end function has_key

  function has_key_ns(dict, uri, localname) result(found)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: uri, localname
    logical :: found

#ifndef DUMMYLIB
    integer  ::  i
    found = .false.
    do i = 1, ubound(dict%list, 1)
      ! FIXME xlf 10.01 segfaults if the below is done as
      ! an AND rather than two separate ifs. This is
      ! probably due to the Heisenbug
      if (uri==str_vs(dict%list(i)%d%nsURI)) then
        if (localname==str_vs(dict%list(i)%d%localname)) then
          found = .true.
          exit
        endif
      endif
    enddo
#else
    found = .false.
#endif
  end function has_key_ns

  pure function get_key_index(dict, key) result(ind)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: key
    integer :: ind

#ifndef DUMMYLIB
    integer  ::  i
    do i = 1, ubound(dict%list, 1)
      if (key == str_vs(dict%list(i)%d%key)) then
        ind = i
        return
      endif
    enddo
    ind = 0
#else
    ind = 0
#endif
  end function get_key_index

  pure function get_key_index_ns(dict, uri, localname) result(ind)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: uri, localname
    integer :: ind

#ifndef DUMMYLIB
    integer  ::  i
    ind = -1
    do  i = 1, ubound(dict%list, 1)
      if (uri==str_vs(dict%list(i)%d%nsURI) &
        .and. localname==str_vs(dict%list(i)%d%localname)) then
        ind = i
        exit
      endif
    enddo
#else
    ind = 0
#endif
  end function get_key_index_ns

#ifndef DUMMYLIB
  pure function get_value_by_key_len(dict, key) result(n)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: key
    integer :: n

    integer :: i

    do i = 1, ubound(dict%list, 1)
      if (key == str_vs(dict%list(i)%d%key)) then
        n = size(dict%list(i)%d%value)
        return
      endif
    enddo
    n = 0
  end function get_value_by_key_len
#endif

  function get_value_by_key(dict, key) result(value)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: key
#ifndef DUMMYLIB
    character(len=get_value_by_key_len(dict,key)) :: value

    integer :: i

    do i = 1, ubound(dict%list, 1)
      if (key == str_vs(dict%list(i)%d%key)) then
        value = str_vs(dict%list(i)%d%value)
        return
      endif
    enddo
    value = ""
#else
    character(len=1) :: value
    value = ""
#endif
  end function get_value_by_key

#ifndef DUMMYLIB
  pure function get_value_by_key_ns_len(dict, nsUri, localname) result(n)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: nsUri
    character(len=*), intent(in) :: localname
    integer :: n

    integer :: i

    do i = 1, ubound(dict%list, 1)
      if (nsUri==str_vs(dict%list(i)%d%nsURI) &
        .and.localname==str_vs(dict%list(i)%d%localname)) then
        n = size(dict%list(i)%d%value)
        return
      endif
    enddo
    n = 0
  end function get_value_by_key_ns_len
#endif

  function get_value_by_key_ns(dict, nsUri, localname) result(value)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: nsUri
    character(len=*), intent(in) :: localname
#ifndef DUMMYLIB
    character(len=get_value_by_key_ns_len(dict, nsURI, localname)) :: value

    integer :: i

    do i = 1, ubound(dict%list, 1)
      if (nsUri==str_vs(dict%list(i)%d%nsURI) &
        .and.localname==str_vs(dict%list(i)%d%localname)) then
        value = str_vs(dict%list(i)%d%value)
        return
      endif
    enddo
    value = ""
#else
    character(len=1) :: value
    value = ""
#endif
  end function get_value_by_key_ns

#ifndef DUMMYLIB
  subroutine get_att_index_pointer(dict, key, i, value)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: key
    integer, intent(out) :: i
    character, pointer :: value(:)

    value => null()
    do i = 1, ubound(dict%list, 1)
      if (key == str_vs(dict%list(i)%d%key)) then
        value => dict%list(i)%d%value
        return
      endif
    enddo
    i = 0
  end subroutine get_att_index_pointer

  subroutine remove_key_by_index(dict, ind)
    type(dictionary_t), intent(inout) :: dict
    integer, intent(in) :: ind

    integer :: i, n
    type(dict_item_ptr), pointer :: tempList(:)

    n = ubound(dict%list, 1)

    if (ind<=0.or.ind>n) return

    allocate(tempList(0:n-1))
    do i = 0, ind-1
      tempList(i)%d => dict%list(i)%d
    enddo
    call destroy(dict%list(ind)%d)
    do i = ind+1, n
      tempList(i-1)%d => dict%list(i)%d
    enddo
    deallocate(dict%list)
    dict%list => tempList
  end subroutine remove_key_by_index
#endif

  pure function get_value_by_index_len(dict, i) result(n)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
    integer :: n

#ifdef DUMMYLIB
    n = 1
#else
    if (i>0.and.i<=ubound(dict%list, 1)) then
      n = size(dict%list(i)%d%value)
    else
      n = 0
    endif
#endif
  end function get_value_by_index_len

  function get_value_by_index(dict, i) result(value)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
#ifndef DUMMYLIB
    character(len=get_value_by_index_len(dict, i)) :: value

    if (i>0.and.i<=ubound(dict%list, 1)) then
      value = str_vs(dict%list(i)%d%value)
    else
      value = ""
    endif
#else
    character(len=1) :: value
    value = ""
#endif
  end function get_value_by_index

#ifndef DUMMYLIB
  pure function get_key_len(dict, i) result(n)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
    integer :: n

    if (i>0.and.i<=ubound(dict%list, 1)) then
      n = size(dict%list(i)%d%key)
    else
      n = 0
    endif
  end function get_key_len
#endif

  function get_key(dict, i) result(key)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
#ifndef DUMMYLIB
    character(len=get_key_len(dict,i)) :: key

    if (i>0.and.i<=ubound(dict%list, 1)) then
      key = str_vs(dict%list(i)%d%key)
    else
      key = ""
    endif
#else
    character(len=1) :: key
    key = ""
#endif
  end function get_key

#ifndef DUMMYLIB
  subroutine add_item_to_dict(dict, key, value, prefix, nsURI, type, itype, specified, declared)
    type(dictionary_t), intent(inout) :: dict
    character(len=*), intent(in)           :: key
    character(len=*), intent(in)           :: value
    character(len=*), intent(in), optional :: prefix
    character(len=*), intent(in), optional :: nsURI
    character(len=*), intent(in), optional :: type
    integer, intent(in), optional :: itype
    logical, intent(in), optional :: specified
    logical, intent(in), optional :: declared

    type(dict_item_ptr), pointer :: tempList(:)
    integer :: i, n

    if (present(prefix) .eqv. .not.present(nsURI)) &
      call FoX_Error('Namespace improperly specified')

    n = ubound(dict%list, 1)
    allocate(tempList(0:n+1))
    do i = 0, n
      tempList(i)%d => dict%list(i)%d
    enddo
    n = n + 1

    allocate(tempList(n)%d)
    tempList(n)%d%value => vs_str_alloc(value)
    if (present(prefix)) then
      tempList(n)%d%key => vs_str_alloc(prefix//":"//key)
      tempList(n)%d%localname => vs_str_alloc(key)
      tempList(n)%d%prefix => vs_str_alloc(prefix)
      tempList(n)%d%nsURI => vs_str_alloc(nsURI)
    else
      tempList(n)%d%key => vs_str_alloc(key)
      tempList(n)%d%localname => vs_str_alloc(key)
      allocate(tempList(n)%d%prefix(0))
      allocate(tempList(n)%d%nsURI(0))
    endif
    if (present(type)) then
      if (present(itype)) &
        call FoX_fatal("internal library error in add_item_to_dict")
      tempList(n)%d%type = get_att_type_enum(type)
    elseif (present(itype)) then
      tempList(n)%d%type = itype
    else
      tempList(n)%d%type = ATT_CDAMB
    endif
    if (present(specified)) then
      tempList(n)%d%specified = specified
    else
      tempList(n)%d%specified = .true.
    endif
    if (present(declared)) then
      tempList(n)%d%declared = declared
    else
      tempList(n)%d%declared = .false.
    endif

    deallocate(dict%list)
    dict%list => tempList

  end subroutine add_item_to_dict

  subroutine set_nsURI_by_index(dict, i, nsURI)
    type(dictionary_t), intent(inout) :: dict
    integer, intent(in) :: i
    character(len=*), intent(in) :: nsURI

    if (associated(dict%list(i)%d%nsURI)) &
      deallocate(dict%list(i)%d%nsURI)
    dict%list(i)%d%nsURI => vs_str_alloc(nsURI)
  end subroutine set_nsURI_by_index

  subroutine set_prefix_by_index(dict, i, prefix)
    type(dictionary_t), intent(inout) :: dict
    integer, intent(in) :: i
    character(len=*), intent(in) :: prefix

    if (associated(dict%list(i)%d%prefix)) &
      deallocate(dict%list(i)%d%prefix)
    dict%list(i)%d%prefix => vs_str_alloc(prefix)
  end subroutine set_prefix_by_index

  subroutine set_localName_by_index_s(dict, i, localName)
    type(dictionary_t), intent(inout) :: dict
    integer, intent(in) :: i
    character(len=*), intent(in) :: localName

    if (associated(dict%list(i)%d%localName)) &
      deallocate(dict%list(i)%d%localName)
    dict%list(i)%d%localName => vs_str_alloc(localName)
  end subroutine set_localName_by_index_s

  subroutine set_localName_by_index_vs(dict, i, localName)
    type(dictionary_t), intent(inout) :: dict
    integer, intent(in) :: i
    character(len=1), dimension(:), intent(in) :: localName

    if (associated(dict%list(i)%d%localName)) &
      deallocate(dict%list(i)%d%localName)
    allocate(dict%list(i)%d%localName(size(localName)))
    dict%list(i)%d%localName = localName
  end subroutine set_localName_by_index_vs
#endif

  pure function get_nsURI_by_index(dict, i) result(nsURI)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
#ifndef DUMMYLIB
    character(len=size(dict%list(i)%d%nsURI)) :: nsURI
    nsURI = str_vs(dict%list(i)%d%nsURI)
#else
    character(len=1) :: nsURI
    nsURI = ""
#endif
  end function get_nsURI_by_index

#ifndef DUMMYLIB
  pure function get_prefix_by_index(dict, i) result(prefix)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
    character(len=size(dict%list(i)%d%prefix)) :: prefix

    prefix = str_vs(dict%list(i)%d%prefix)
  end function get_prefix_by_index
#endif

  pure function get_localName_by_index(dict, i) result(localName)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
#ifndef DUMMYLIB
    character(len=size(dict%list(i)%d%localName)) :: localName
    localName = str_vs(dict%list(i)%d%localName)
#else
    character(len=1) :: localName
    localName = ""
#endif
  end function get_localName_by_index

#ifndef DUMMYLIB
  pure function get_nsURI_by_keyname_len(dict, keyname) result(n)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: keyname
    integer :: n

    integer :: i

    i = get_key_index(dict, keyname)
    n = size(dict%list(i)%d%nsURI)
  end function get_nsURI_by_keyname_len
#endif

#ifndef DUMMYLIB
  pure function get_nsURI_by_keyname(dict, keyname) result(nsURI)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: keyname
    character(len=get_nsURI_by_keyname_len(dict, keyname)) :: nsURI
    integer :: i

    i = get_key_index(dict, keyname)
    nsURI = str_vs(dict%list(i)%d%nsURI)
  end function get_nsURI_by_keyname

  pure function get_prefix_by_keyname_len(dict, keyname) result(n)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: keyname
    integer :: n

    integer :: i

    i = get_key_index(dict, keyname)
    n = size(dict%list(i)%d%prefix)

  end function get_prefix_by_keyname_len

  pure function get_prefix_by_keyname(dict, keyname) result(prefix)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: keyname
    character(len=get_prefix_by_keyname_len(dict,keyname)) :: prefix
    integer :: i

    i = get_key_index(dict, keyname)
    prefix = str_vs(dict%list(i)%d%prefix)

  end function get_prefix_by_keyname

  pure function get_localname_by_keyname_len(dict, keyname) result(n)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: keyname
    integer :: n

    integer :: i

    i = get_key_index(dict, keyname)
    n = size(dict%list(i)%d%localName)

  end function get_localname_by_keyname_len

  pure function get_localName_by_keyname(dict, keyname) result(localName)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: keyname
    character(get_localname_by_keyname_len(dict, keyname)) :: localname
    integer :: i

    i=get_key_index(dict, keyname)
    localName = str_vs(dict%list(i)%d%localName)

  end function get_localName_by_keyname
#endif

#ifndef DUMMYLIB
  pure function getType_by_index_len(dict, i) result(n)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
    integer :: n

    if (i>0.and.i<=ubound(dict%list, 1)) then
      n = len_trim(ATT_TYPES(dict%list(i)%d%type))
    else
      n = 0
    endif
  end function getType_by_index_len
#endif

  function getType_by_index(dict, i) result(type)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
#ifndef DUMMYLIB
    character(len=getType_by_index_len(dict, i)) :: type

    if (i>0.and.i<=ubound(dict%list, 1)) then
      type = ATT_TYPES(dict%list(i)%d%type)
    else
      type = ""
    endif
#else
    character(len=1) :: type
    type = ""
#endif
  end function getType_by_index

#ifndef DUMMYLIB
  pure function getType_by_keyname_len(dict, keyname) result(n)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: keyname
    integer :: n

    integer :: i
    i = get_key_index(dict, keyname)

    if (i>0) then
      n = len_trim(ATT_TYPES(i))
    else
      n = 0
    endif
  end function getType_by_keyname_len
#endif

  function getType_by_keyname(dict, keyname) result(type)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: keyname
#ifndef DUMMYLIB
    character(len=getType_by_keyname_len(dict, keyname)) :: type

    integer :: i
    i = get_key_index(dict, keyname)
    if (i>0) then
      type = ATT_TYPES(dict%list(i)%d%type)
    else
      type = ""
    endif
#else
    character(len=1) :: type
    type = ""
#endif
  end function getType_by_keyname

  function isSpecified_by_index(dict, i) result(p)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
    logical :: p

#ifndef DUMMYLIB
    if (i>0.and.i<=ubound(dict%list, 1)) then
      p = dict%list(i)%d%specified
    else
      p = .false.
    endif
#else
    p = .false.
#endif
  end function isSpecified_by_index

  function isSpecified_by_key(dict, qName) result(p)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: qName
    logical :: p

#ifndef DUMMYLIB
    integer :: i
    i = getIndex(dict, qName)
    if (i>0.and.i<=ubound(dict%list, 1)) then
      p = dict%list(i)%d%specified
    else
      p = .false.
    endif
#else
    p = .false.
#endif
  end function isSpecified_by_key

  subroutine setSpecified(dict, i, p)
    type(dictionary_t), intent(inout) :: dict
    integer, intent(in) :: i
    logical, intent(in) :: p

#ifndef DUMMYLIB
    if (i>0.and.i<=ubound(dict%list, 1)) then
      dict%list(i)%d%specified = p
    endif
#endif
  end subroutine setSpecified

  function isSpecified_by_keyNS(dict, uri, localName) result(p)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: uri
    character(len=*), intent(in) :: localName
    logical :: p

#ifndef DUMMYLIB
    integer :: i
    i = getIndex(dict, uri, localName)
    if (i>0.and.i<=ubound(dict%list, 1)) then
      p = dict%list(i)%d%specified
    else
      p = .false.
    endif
#else
    p = .false.
#endif
  end function isSpecified_by_keyNS

  function isDeclared_by_index(dict, i) result(p)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
    logical :: p

#ifndef DUMMYLIB
    if (i>0.and.i<=ubound(dict%list, 1)) then
      p = dict%list(i)%d%declared
    else
      p = .false.
    endif
#else
    p = .false.
#endif
  end function isDeclared_by_index

  function isDeclared_by_key(dict, qName) result(p)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: qName
    logical :: p

#ifndef DUMMYLIB
    integer :: i
    i = getIndex(dict, qName)
    if (i>0.and.i<=ubound(dict%list, 1)) then
      p = dict%list(i)%d%declared
    else
      p = .false.
    endif
#else
    p = .false.
#endif
  end function isDeclared_by_key

  function isDeclared_by_keyNS(dict, uri, localName) result(p)
    type(dictionary_t), intent(in) :: dict
    character(len=*), intent(in) :: uri
    character(len=*), intent(in) :: localName
    logical :: p

#ifndef DUMMYLIB
    integer :: i
    i = getIndex(dict, uri, localName)
    if (i>0.and.i<=ubound(dict%list, 1)) then
      p = dict%list(i)%d%declared
    else
      p = .false.
    endif
#else
    p = .false.
#endif
  end function isDeclared_by_keyNS

  subroutine setDeclared(dict, i, p)
    type(dictionary_t), intent(inout) :: dict
    integer, intent(in) :: i
    logical, intent(in) :: p

#ifndef DUMMYLIB
    if (i>0.and.i<=ubound(dict%list, 1)) then
      dict%list(i)%d%declared = p
    endif
#endif
  end subroutine setDeclared

#ifndef DUMMYLIB
  function getIsId_by_index(dict, i) result(p)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
    logical :: p

    if (i>0.and.i<=ubound(dict%list, 1)) then
      p = dict%list(i)%d%isId
    else
      p = .false.
    endif

  end function getIsId_by_index

  subroutine setIsId_by_index(dict, i, p)
    type(dictionary_t), intent(inout) :: dict
    integer, intent(in) :: i
    logical, intent(in) :: p

    if (i>0.and.i<=ubound(dict%list, 1)) then
      dict%list(i)%d%isId = p
    endif
  end subroutine setIsId_by_index

  function getWhitespaceHandling(dict, i) result(j)
    type(dictionary_t), intent(in) :: dict
    integer, intent(in) :: i
    integer :: j

    if (i<=ubound(dict%list, 1)) then
      select case(dict%list(i)%d%type)
      case (ATT_CDATA)
        j = 0 !
      case (ATT_CDAMB)
        j = 1
      case default
        j = 2
      end select
    else
      j = 2
    endif

  end function getWhitespaceHandling

  subroutine setBase(dict, base)
    type(dictionary_t), intent(inout) :: dict
    character(len=*), intent(in) :: base

    if (associated(dict%base)) deallocate(dict%base)
    dict%base => vs_str_alloc(base)
  end subroutine setBase

  pure function getBase_len(dict) result(n)
    type(dictionary_t), intent(in) :: dict
    integer :: n

    if (associated(dict%base)) then
      n = size(dict%base)
    else
      n = 0
    endif
  end function getBase_len

  function getBase(dict) result(base)
    type(dictionary_t), intent(in) :: dict
    character(len=getBase_len(dict)) :: base

    if (associated(dict%base)) then
      base = str_vs(dict%base)
    else
      base = ""
    endif
  end function getBase

  subroutine destroy_dict_item(d)
    type(dict_item), pointer :: d

    if (associated(d)) then
      deallocate(d%key)
      deallocate(d%value)
      deallocate(d%nsURI)
      deallocate(d%prefix)
      deallocate(d%localName)
      deallocate(d)
    endif
  end subroutine destroy_dict_item

  subroutine init_dict(dict)
    type(dictionary_t), intent(out) :: dict


    allocate(dict%list(0:0))
    allocate(dict%list(0)%d)
    allocate(dict%list(0)%d%key(0))

  end subroutine init_dict

  subroutine destroy_dict(dict)
    type(dictionary_t), intent(inout) :: dict
    integer :: i

    if (associated(dict%list)) then
      deallocate(dict%list(0)%d%key)
      deallocate(dict%list(0)%d)
      do i = 1, ubound(dict%list, 1)
        call destroy(dict%list(i)%d)
      enddo
      deallocate(dict%list)
    endif
    if (associated(dict%base)) deallocate(dict%base)

  end subroutine destroy_dict


  subroutine reset_dict(dict)
    type(dictionary_t), intent(inout) :: dict

    call destroy_dict(dict)
    call init_dict(dict)

  end subroutine reset_dict

  subroutine sortAttrs(dict)
    type(dictionary_t), intent(inout) :: dict

    logical :: done(ubound(dict%list, 1))
    type(dict_item_ptr), dimension(:), pointer :: list => null()
    integer :: i, j, n, firstIndex
    character, pointer :: firstKey(:)

    !Ridiculously naive sort algorithm. We are unlikely
    !to ever be sorting more then ten or so attributes though

    n = ubound(dict%list, 1)

    allocate(list(0:n))
    list(0)%d => dict%list(0)%d

    j = 1
    done = .false.

    firstIndex = 1
    do while (firstIndex/=0)
      firstIndex = 0
      firstKey => null()
      do i = 1, n
        if (.not.done(i).and.str_vs(dict%list(i)%d%key)=="xmlns" &
          .or. str_vs(dict%list(i)%d%prefix)=="xmlns") then
          firstIndex = i
          if (associated(firstKey)) then
            if (llt(str_vs(dict%list(i)%d%key),str_vs(firstKey))) then
              firstIndex = i
              firstKey => dict%list(i)%d%key
            endif
          else
            firstIndex = i
            firstKey => dict%list(i)%d%key
          endif
        endif
      enddo
      if (firstIndex/=0) then
        done(firstIndex) = .true.
        list(j)%d => dict%list(firstIndex)%d
        j = j + 1
      endif
    enddo

    do while (any(.not.done))
      firstIndex = 0
      firstKey => null()
      do i = 1, n
        if (.not.done(i)) then
          if (associated(firstKey)) then
            if (llt(str_vs(dict%list(i)%d%key),str_vs(firstKey))) then
              firstIndex = i
              firstKey => dict%list(i)%d%key
            endif
          else
            firstIndex = i
            firstKey => dict%list(i)%d%key
          endif
        endif
      enddo
      done(firstIndex) = .true.
      list(j)%d => dict%list(firstIndex)%d
      j = j + 1
    enddo

    deallocate(dict%list)
    dict%list => list

  end subroutine sortAttrs


  subroutine print_dict(dict)
    type(dictionary_t), intent(in) :: dict

    integer  :: i

    do i = 1, ubound(dict%list, 1)
      write(*,'(7a)') str_vs(dict%list(i)%d%key), " [ {", str_vs(dict%list(i)%d%nsURI), &
        "}", str_vs(dict%list(i)%d%localName), " ]  = ", str_vs(dict%list(i)%d%value)
    enddo

  end subroutine print_dict

#endif
end module m_common_attrs
