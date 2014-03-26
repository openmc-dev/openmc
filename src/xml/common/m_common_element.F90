module m_common_element

#ifndef DUMMYLIB
  ! Structure and manipulation of element specification

  use fox_m_fsys_array_str, only: str_vs, vs_str_alloc, vs_vs_alloc
  use fox_m_fsys_string_list, only: string_list, init_string_list, &
    destroy_string_list, add_string, tokenize_to_string_list, &
    registered_string
  use m_common_charset, only: isInitialNameChar, isNameChar, &
    upperCase, XML_WHITESPACE
  use m_common_content_model, only: content_particle_t, newCP, destroyCPtree, &
    OP_MIXED, OP_CHOICE, OP_SEQ, OP_NAME, &
    REP_QUESTION_MARK, REP_ASTERISK, &
    transformCPPlus ! , dumpCPtree ! For debugging - see below.
  use m_common_error, only: error_stack, add_error, in_error
  use m_common_namecheck, only: checkName, checkNames, checkNCName, &
    checkNCNames, checkQName, checkNmtoken, checkNmtokens

  implicit none
  private

  integer, parameter :: ST_START               = 0
  integer, parameter :: ST_EMPTYANY            = 1
  integer, parameter :: ST_FIRSTCHILD          = 2
  integer, parameter :: ST_END                 = 3
  integer, parameter :: ST_PCDATA              = 4
  integer, parameter :: ST_NAME                = 5
  integer, parameter :: ST_CHILD               = 6
  integer, parameter :: ST_AFTERBRACKET        = 7
  integer, parameter :: ST_AFTERLASTBRACKET    = 8
  integer, parameter :: ST_SEPARATOR           = 9
  integer, parameter :: ST_AFTERNAME           = 10
  integer, parameter :: ST_ATTTYPE             = 11
  integer, parameter :: ST_AFTER_NOTATION      = 12
  integer, parameter :: ST_NOTATION_LIST       = 13
  integer, parameter :: ST_ENUMERATION         = 14
  integer, parameter :: ST_ENUM_NAME           = 15
  integer, parameter :: ST_AFTER_ATTTYPE_SPACE = 16
  integer, parameter :: ST_AFTER_ATTTYPE       = 17
  integer, parameter :: ST_DEFAULT_DECL        = 18
  integer, parameter :: ST_AFTERDEFAULTDECL    = 19
  integer, parameter :: ST_DEFAULTVALUE        = 20

  integer, parameter :: ATT_NULL = 0

  integer, parameter :: ATT_CDATA = 1
  integer, parameter :: ATT_ID = 2
  integer, parameter :: ATT_IDREF = 3
  integer, parameter :: ATT_IDREFS = 4
  integer, parameter :: ATT_ENTITY = 5
  integer, parameter :: ATT_ENTITIES = 6
  integer, parameter :: ATT_NMTOKEN = 7
  integer, parameter :: ATT_NMTOKENS = 8
  integer, parameter :: ATT_NOTATION = 9
  integer, parameter :: ATT_ENUM = 10
  integer, parameter :: ATT_CDANO = 11
  integer, parameter :: ATT_CDAMB = 12

  character(len=8), parameter :: ATT_TYPES(12) = (/ &
    "CDATA   ", &
    "ID      ", &
    "IDREF   ", &
    "IDREFS  ", &
    "ENTITY  ", &
    "ENTITIES", &
    "NMTOKEN ", &
    "NMTOKENS", &
    "NOTATION", &
    "ENUM    ", &
    "CDANO   ", &
    "CDAMB   "/)

  integer, parameter :: ATT_REQUIRED = 1
  integer, parameter :: ATT_IMPLIED = 2
  integer, parameter :: ATT_DEFAULT = 4
  integer, parameter :: ATT_FIXED = 3


  type attribute_t
    character, pointer :: name(:) => null()
    integer :: attType = ATT_NULL 
    integer :: attDefault = ATT_NULL
    type(string_list) :: enumerations
    character, pointer :: default(:) => null()
    logical :: internal = .true.
  end type attribute_t

  type attribute_list
    type(attribute_t), pointer :: list(:) => null()
  end type attribute_list

  type element_t
    character, pointer :: name(:) => null()
    logical :: empty = .false.
    logical :: any = .false.
    logical :: mixed = .false.
    logical :: id_declared = .false.
    logical :: internal = .true.
    type (content_particle_t), pointer :: cp => null()
    character, pointer :: model(:) => null()
    type(attribute_list) :: attlist
  end type element_t

  type element_list
    type(element_t), pointer :: list(:) => null()
  end type element_list


  public :: element_t
  public :: element_list

  public :: attribute_t
  public :: attribute_list

  public :: init_element_list
  public :: destroy_element_list
  public :: existing_element
  public :: declared_element
  public :: get_element
  public :: add_element

  public :: parse_dtd_element

  public :: init_attribute_list
  public :: destroy_attribute_list


  public :: parse_dtd_attlist

  public :: report_declarations

  public :: attribute_has_default
  public :: get_attlist_size
  public :: get_attribute_declaration
  public :: express_attribute_declaration

  public :: att_value_normalize

  public :: get_att_type_enum

  public :: ATT_NULL
  public :: ATT_CDATA
  public :: ATT_ID 
  public :: ATT_IDREF
  public :: ATT_IDREFS
  public :: ATT_ENTITY
  public :: ATT_ENTITIES
  public :: ATT_NMTOKEN
  public :: ATT_NMTOKENS
  public :: ATT_NOTATION
  public :: ATT_ENUM

  public :: ATT_CDANO
  public :: ATT_CDAMB

  public :: ATT_REQUIRED
  public :: ATT_IMPLIED
  public :: ATT_DEFAULT
  public :: ATT_FIXED

  public :: ATT_TYPES

  interface get_attribute_declaration
    module procedure get_attdecl_by_index
    module procedure get_attdecl_by_name
  end interface

contains

  subroutine init_element_list(e_list)
    type(element_list), intent(inout) :: e_list

    allocate(e_list%list(0))
  end subroutine init_element_list

  subroutine destroy_element_list(e_list)
    type(element_list), intent(inout) :: e_list

    integer :: i

    do i = 1, size(e_list%list)
      deallocate(e_list%list(i)%name)
      if (associated(e_list%list(i)%cp)) call destroyCPtree(e_list%list(i)%cp)
      if (associated(e_list%list(i)%model)) deallocate(e_list%list(i)%model)
      call destroy_attribute_list(e_list%list(i)%attlist)
    enddo
    deallocate(e_list%list)
  end subroutine destroy_element_list

  function existing_element(e_list, name) result(p)
    type(element_list), intent(in) :: e_list
    character(len=*), intent(in) :: name
    logical :: p

    integer :: i

    p = .false.
    do i = 1, size(e_list%list)
      if (str_vs(e_list%list(i)%name)==name) then
        p = .true.
        exit
      endif
    enddo
  end function existing_element

  function declared_element(e_list, name) result(p)
    type(element_list), intent(in) :: e_list
    character(len=*), intent(in) :: name
    logical :: p

    integer :: i

    p = .false.
    do i = 1, size(e_list%list)
      if (str_vs(e_list%list(i)%name)==name) then
        p = associated(e_list%list(i)%model)
        exit
      endif
    enddo
  end function declared_element

  function get_element(e_list, name) result(e)
    type(element_list), intent(in) :: e_list
    character(len=*), intent(in) :: name
    type(element_t), pointer :: e

    integer :: i

    do i = 1, size(e_list%list)
      if (str_vs(e_list%list(i)%name)==name) then
        e => e_list%list(i)
        return
      endif
    enddo
    e => null()
  end function get_element

  function add_element(e_list, name) result(e)
    type(element_list), intent(inout) :: e_list
    character(len=*), intent(in) :: name
    type(element_t), pointer :: e

    type(element_t), pointer :: temp(:)
    integer :: i

    temp => e_list%list

    allocate(e_list%list(size(temp)+1))
    do i = 1, size(temp)
      e_list%list(i)%name => temp(i)%name
      e_list%list(i)%model => temp(i)%model
      e_list%list(i)%empty = temp(i)%empty
      e_list%list(i)%any = temp(i)%any
      e_list%list(i)%mixed = temp(i)%mixed
      e_list%list(i)%cp => temp(i)%cp
      e_list%list(i)%id_declared = temp(i)%id_declared
      e_list%list(i)%internal = temp(i)%internal
      e_list%list(i)%attlist%list => temp(i)%attlist%list
    enddo
    deallocate(temp)
    e => e_list%list(i)
    e%name => vs_str_alloc(name)
    call init_attribute_list(e%attlist)

  end function add_element

  subroutine parse_dtd_element(contents, xv, stack, element, internal)
    character(len=*), intent(in) :: contents
    integer, intent(in) :: xv
    type(error_stack), intent(inout) :: stack
    type(element_t), pointer :: element
    logical, intent(in) :: internal

    integer :: state
    integer :: i, nbrackets
    logical :: mixed, empty, any
    character :: c
    character, pointer :: order(:), name(:), temp(:)
    type(content_particle_t), pointer :: top, current, tcp
    logical :: mixed_additional, firstChild

    ! FIXME should we check namespaces here (for element names)
    ! checking duplicates - valid or wf? - and only for MIXED?

    order => null()
    name => null()
    temp => null()

    any = .false.
    empty = .false.
    mixed = .false.
    nbrackets = 0
    mixed_additional = .false.
    firstChild = .true.
    state = ST_START

    top => null()

    do i = 1, len(contents) + 1
      if (i<=len(contents)) then
        c = contents(i:i)
      else
        c = ' '
      endif

      if (state==ST_START) then
        !write(*,*)'ST_START'
        if (verify(c, XML_WHITESPACE)==0) then
          continue
        elseif (verify(c, 'EMPTYANY')==0) then
          name => vs_str_alloc(c)
          state = ST_EMPTYANY
        elseif (c=='(') then
          order => vs_str_alloc(" ")
          nbrackets = 1
          top => newCP()
          current => top
          state = ST_FIRSTCHILD
        else
          call add_error(stack, &
            'Unexpected character "'//c//'" at start of ELEMENT specification')
          goto 100
        endif

      elseif (state==ST_EMPTYANY) then
        !write(*,*)'ST_EMPTYANY'
        if (verify(c, upperCase)==0) then
          temp => name
          name => vs_str_alloc(str_vs(temp)//c)
          deallocate(temp)
        elseif (verify(c, XML_WHITESPACE)==0) then
          if (str_vs(name)=='EMPTY') then
            empty = .true.
            top => newCP(empty=.true.)
            current => top
          elseif (str_vs(name)=='ANY') then
            any = .true.
            top => newCP(any=.true.)
            current => top
          else
            call add_error(stack, &
              'Unexpected ELEMENT specification; expecting EMPTY or ANY')
            goto 100
          endif
          deallocate(name)
          state = ST_END
        else
          call add_error(stack, &
            'Unexpected ELEMENT specification; expecting EMPTY or ANY')
          goto 100
        endif

      elseif (state==ST_FIRSTCHILD) then
        !write(*,*)'ST_FIRSTCHILD'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (c=='#') then
          mixed = .true.
          state = ST_PCDATA
          name => vs_str_alloc("")
        elseif (isInitialNameChar(c, xv)) then
          allocate(name(1))
          name(1) = c
          state = ST_NAME
        elseif (c=='(') then
          nbrackets = nbrackets + 1
          deallocate(order)
          tcp => newCP()
          current%firstChild => tcp
          tcp%parent => current
          current => tcp
          order => vs_str_alloc("  ")
          state = ST_CHILD
        else
          call add_error(stack, &
            'Unexpected character in ELEMENT specification')
          goto 100
        endif

      elseif (state==ST_PCDATA) then
        !write(*,*)'ST_PCDATA'
        if (verify(c, 'PCDATA')==0) then
          temp => name
          name => vs_str_alloc(str_vs(temp)//c)
          deallocate(temp)
        elseif (verify(c, XML_WHITESPACE)==0) then
          if (str_vs(name)=='PCDATA') then
            deallocate(name)
          else
            call add_error(stack, &
              'Unexpected token after #')
            goto 100
          endif
          ! Must be first child
          current%operator = OP_MIXED
          tcp => newCP(name="#PCDATA")
          current%firstChild => tcp
          tcp%parent => current
          current => tcp
          firstChild = .false.
          state = ST_SEPARATOR
        elseif (c==')') then
          if (str_vs(name)=='PCDATA') then
            deallocate(name)
            nbrackets = 0
            state = ST_AFTERLASTBRACKET
            deallocate(order)
          else
            call add_error(stack, &
              'Unexpected token after #')
            goto 100
          endif
          ! Must be first child
          current%operator = OP_MIXED
          tcp => newCP(name="#PCDATA")
          current%firstChild => tcp
          tcp%parent => current
          firstChild = .false.
        elseif (c=='|') then
          if (str_vs(name)=='PCDATA') then
            firstChild = .false.
            deallocate(name)
          else
            call add_error(stack, &
              'Unexpected token after #')
            goto 100
          endif
          ! Must be first child
          current%operator = OP_MIXED
          tcp => newCP(name="#PCDATA")
          current%firstChild => tcp
          tcp%parent => current
          current => tcp
          firstChild = .false.
          order(1) = '|'
          state = ST_CHILD
        elseif (c==',') then
          call add_error(stack, &
            'Ordered specification not allowed for Mixed elements')
          goto 100
        else
          call add_error(stack, &
            'Unexpected character in ELEMENT specification')
          goto 100
        endif

      elseif (state==ST_NAME) then
        !write(*,*)'ST_NAME'
        if (isNameChar(c, xv)) then
          temp => name
          name => vs_str_alloc(str_vs(temp)//c)
          deallocate(temp)
        elseif (scan(c, "?+*")>0) then
          if (mixed) then
            call add_error(stack, &
              'Repeat operators forbidden for Mixed elements')
            goto 100
          endif
          tcp => newCP(name=str_vs(name), repeat=c)
          deallocate(name)
          if (firstChild) then
            current%firstChild => tcp
            tcp%parent => current
            firstChild = .false.
          else
            current%nextSibling => tcp
            tcp%parent => current%parent
          endif
          current => tcp
          if (c=="+") call transformCPPlus(current)
          state = ST_SEPARATOR
        elseif (verify(c, XML_WHITESPACE)==0) then
          if (mixed) mixed_additional = .true.
          tcp => newCP(name=str_vs(name))
          deallocate(name)
          if (firstChild) then
            current%firstChild => tcp
            tcp%parent => current
            firstChild = .false.
          else
            current%nextSibling => tcp
            tcp%parent => current%parent
          endif
          current => tcp
          state = ST_SEPARATOR
        elseif (scan(c,',|')>0) then
          if (order(nbrackets)=='') then
            order(nbrackets)=c
          elseif (order(nbrackets)/=c) then
            call add_error(stack, &
              'Cannot mix ordered and unordered elements')
            goto 100
          endif
          if (mixed) mixed_additional = .true.
          tcp => newCP(name=str_vs(name))
          deallocate(name)
          if (firstChild) then
            current%firstChild => tcp
            tcp%parent => current
            firstChild = .false.
          else
            current%nextSibling => tcp
            tcp%parent => current%parent
          endif
          current => tcp
          if (c=="|".and.current%parent%operator/=OP_MIXED) &
            current%parent%operator = OP_CHOICE
          state = ST_CHILD
        elseif (c==')') then
          if (mixed) mixed_additional = .true.
          nbrackets = nbrackets - 1
          if (nbrackets==0) then
            state = ST_AFTERLASTBRACKET
            deallocate(order)
          else
            temp => order
            allocate(order(nbrackets))
            order = temp(:size(order))
            deallocate(temp)
            state = ST_AFTERBRACKET
          endif
          tcp => newCP(name=str_vs(name))
          deallocate(name)
          if (firstChild) then
            current%firstChild => tcp
            tcp%parent => current
            firstChild = .false.
          else
            current%nextSibling => tcp
            tcp%parent => current%parent
            current => current%parent
            if (.not.check_duplicates(current)) &
              goto 100
          endif
        else
          call add_error(stack, &
            'Unexpected character found after element name')
          goto 100
        endif

      elseif (state==ST_CHILD) then
        !write(*,*)'ST_CHILD'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (c=='#') then
          call add_error(stack, &
            '# forbidden except as first child element')
          goto 100
        elseif (isInitialNameChar(c, xv)) then
          name => vs_str_alloc(c)
          state = ST_NAME
        elseif (c=='(') then
          if (mixed) then
            call add_error(stack, &
              'Nested brackets forbidden for Mixed content')
            goto 100
          endif
          tcp => newCP()
          if (firstChild) then
            current%firstChild => tcp
            tcp%parent => current
          else
            current%nextSibling => tcp
            tcp%parent => current%parent
            firstChild = .true.
          endif
          current => tcp
          nbrackets = nbrackets + 1
          temp => order
          order => vs_str_alloc(str_vs(temp)//" ")
          deallocate(temp)
        else
          call add_error(stack, &
            'Unexpected character "'//c//'" found after (')
          goto 100
        endif

      elseif (state==ST_SEPARATOR) then
        !write(*,*)'ST_SEPARATOR'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (c=='#') then
          call add_error(stack, &
            '#PCDATA must be first in list')
          goto 100
        elseif (scan(c,'|,')>0) then
          if (order(nbrackets)=='') then
            order(nbrackets) = c
          elseif (order(nbrackets)/=c) then
            call add_error(stack, &
              'Cannot mix ordered and unordered elements')
            goto 100
          endif
          if (c=="|".and.current%parent%operator/=OP_MIXED) &
            current%parent%operator = OP_CHOICE
          state = ST_CHILD
        elseif (c==')') then
          nbrackets = nbrackets - 1
          if (nbrackets==0) then
            state = ST_AFTERLASTBRACKET
            deallocate(order)
          else
            temp => order
            allocate(order(nbrackets))
            order = temp(:size(order))
            deallocate(temp)
            state = ST_AFTERBRACKET
          endif
          current => current%parent
          if (.not.check_duplicates(current)) &
            goto 100
        else
          call add_error(stack, &
            'Unexpected character found in element declaration.')
          goto 100
        endif

      elseif (state==ST_AFTERBRACKET) then
        !write(*,*)'ST_AFTERBRACKET'
        if (c=='*') then
          current%repeater = REP_ASTERISK
          state = ST_SEPARATOR
        elseif (c=='+') then
          call transformCPPlus(current)
          state = ST_SEPARATOR
        elseif (c=='?') then
          current%repeater = REP_QUESTION_MARK
          state = ST_SEPARATOR
        elseif (verify(c, XML_WHITESPACE)==0) then
          state = ST_SEPARATOR
        elseif (scan(c,'|,')>0) then
          if (order(nbrackets)=='') then
            order(nbrackets) = c
          elseif (order(nbrackets)/=c) then
            call add_error(stack, &
              'Cannot mix ordered and unordered elements')
            goto 100
          endif
          if (c=="|".and.current%parent%operator/=OP_MIXED) &
            current%parent%operator = OP_CHOICE
          state = ST_CHILD
        elseif (c==')') then
          nbrackets = nbrackets - 1
          if (nbrackets==0) then
            deallocate(order)
            state = ST_AFTERLASTBRACKET
          else
            temp => order
            allocate(order(nbrackets))
            order = temp(:size(order))
            deallocate(temp)
            state = ST_AFTERBRACKET
          endif
          current => current%parent
          if (.not.check_duplicates(current)) &
            goto 100
        else
          call add_error(stack, &
            'Unexpected character "'//c//'"found after ")"')
          goto 100
        endif

      elseif (state==ST_AFTERLASTBRACKET) then
        !write(*,*)'ST_AFTERLASTBRACKET'
        if (c=='*') then
          state = ST_END
          current%repeater = REP_ASTERISK
        elseif (c=='+') then
          if (mixed) then
            call add_error(stack, &
              '+ operator disallowed for Mixed elements')
            goto 100
          endif
          call transformCPPlus(current)
          state = ST_END
        elseif (c=='?') then
          if (mixed) then
            call add_error(stack, &
              '? operator disallowed for Mixed elements')
            goto 100
          endif
          current%repeater = REP_QUESTION_MARK
          state = ST_END
        elseif (verify(c, XML_WHITESPACE)==0) then
          if (mixed) then
            if (mixed_additional) then
              call add_error(stack, &
                'Missing "*" at end of Mixed element specification')
              goto 100
            endif
          endif
          state = ST_END
        else
          call add_error(stack, &
            'Unexpected character "'//c//'" found after final ")"')
          goto 100
        endif

      elseif (state==ST_END) then
        !write(*,*)'ST_END'
        if (verify(c, XML_WHITESPACE)==0) then
          continue
        else
          call add_error(stack, &
            'Unexpected token found after end of element specification')
          goto 100
        endif

      endif

    enddo

    if (state/=ST_END) then
      call add_error(stack, "Error in parsing contents of element declaration")
      goto 100
    endif

    if (associated(element)) then
      element%any = any
      element%empty = empty
      element%mixed = mixed
      element%model => vs_str_alloc(trim(strip_spaces(contents)))
      element%cp => top
      element%internal = internal
! For debugging it may be useful to dump the result here...
! Also need to use the subroutine.
!      call dumpCPtree(top)
    else
      if (associated(top)) call destroyCPtree(top)
    endif
    return

100 if (associated(order)) deallocate(order)
    if (associated(name)) deallocate(name)
    if (associated(top)) call destroyCPtree(top)

    contains
      function strip_spaces(s1) result(s2)
        character(len=*) :: s1
        character(len=len(s1)) :: s2
        integer :: i, i2
        i2 = 1
        do i = 1, len(s1)
          if (verify(s1(i:i), XML_WHITESPACE)==0) cycle
          s2(i2:i2) = s1(i:i)
          i2 = i2 + 1
        end do
        s2(i2:) = ''
      end function strip_spaces

      function check_duplicates(cp) result(p)
        type(content_particle_t), pointer :: cp
        logical :: p

        type(string_list) :: sl
        type(content_particle_t), pointer :: tcp

        if (cp%operator==OP_SEQ) then
          p = .true.
          return
        endif

        call init_string_list(sl)
        tcp => cp%firstChild
        p = .false.
        do while (associated(tcp))
          if (tcp%operator==OP_NAME) then
            if (registered_string(sl, str_vs(tcp%name))) then
              call destroy_string_list(sl)
              if (cp%operator==OP_MIXED) then
                call add_error(stack, &
                  "Duplicate element names found in MIXED")
              elseif (cp%operator==OP_CHOICE) then
                call add_error(stack, &
                  "Duplicate element names found in CHOICE")
              endif
              return
            else
              call add_string(sl, str_vs(tcp%name))
            endif
          endif
          tcp => tcp%nextSibling
        enddo
        p = .true.
        call destroy_string_list(sl)
      end function check_duplicates
  end subroutine parse_dtd_element


  subroutine init_attribute_list(a_list)
    type(attribute_list), intent(inout) :: a_list

    allocate(a_list%list(0))
  end subroutine init_attribute_list

  subroutine destroy_attribute_t(a)
    type(attribute_t), pointer :: a 

    if (associated(a%name)) deallocate(a%name)
    if (associated(a%default)) deallocate(a%default)
    call destroy_string_list(a%enumerations)

    deallocate(a)
  end subroutine destroy_attribute_t

  subroutine destroy_attribute_list(a_list)
    type(attribute_list), intent(inout) :: a_list

    integer :: i

    do i = 1, size(a_list%list)
      deallocate(a_list%list(i)%name)
      if (associated(a_list%list(i)%default)) deallocate(a_list%list(i)%default)
      call destroy_string_list(a_list%list(i)%enumerations)
    enddo
    deallocate(a_list%list)

  end subroutine destroy_attribute_list

  function existing_attribute(a_list, name) result(p)
    type(attribute_list), intent(inout) :: a_list
    character(len=*), intent(in) :: name
    logical :: p

    integer :: i
    p = .false.
    do i = 1, size(a_list%list)
      p = (str_vs(a_list%list(i)%name)==name)
      if (p) exit
    enddo
  end function existing_attribute

  function add_attribute(a_list, name, internal) result(a)
    type(attribute_list), intent(inout) :: a_list
    character(len=*), intent(in) :: name
    logical, intent(in) :: internal
    type(attribute_t), pointer :: a

    integer :: i
    type(attribute_t), pointer :: temp(:)

    temp => a_list%list
    allocate(a_list%list(size(temp)+1))
    do i = 1, size(temp)
      a_list%list(i)%name => temp(i)%name
      a_list%list(i)%atttype = temp(i)%atttype
      a_list%list(i)%attdefault = temp(i)%attdefault
      a_list%list(i)%default => temp(i)%default
      a_list%list(i)%enumerations%list => temp(i)%enumerations%list
      a_list%list(i)%internal = temp(i)%internal
    enddo
    deallocate(temp)
    a => a_list%list(i)

    a%name => vs_str_alloc(name)
    call init_string_list(a%enumerations)
    a%internal = internal

  end function add_attribute
  
  function get_attribute(a_list, name) result(a)
    type(attribute_list), intent(inout) :: a_list
    character(len=*), intent(in) :: name
    type(attribute_t), pointer :: a

    integer :: i
    do i = 1, size(a_list%list)
      if (str_vs(a_list%list(i)%name)==name) then
        a => a_list%list(i)
        exit
      endif
    enddo
  end function get_attribute

  subroutine parse_dtd_attlist(contents, xv, namespaces, validCheck, stack, elem, internal)
    character(len=*), intent(in) :: contents
    integer, intent(in) :: xv
    logical, intent(in) :: validCheck
    logical, intent(in) :: namespaces
    type(error_stack), intent(inout) :: stack
    type(element_t), pointer :: elem
    logical, intent(in) :: internal

    integer :: i
    integer :: state
    character :: c, q
    character, pointer :: name(:), attType(:), default(:), value(:), temp(:)

    type(attribute_t), pointer :: ca
    type(attribute_t), pointer :: ignore_att

    ignore_att => null()
    ! We need ignore_att to process but not take account of duplicate attributes
    ! elem is optional so we can not record declarations if necessary.
    ca => null()
    name => null()
    attType => null()
    default => null()
    value => null()
    temp => null()

    state = ST_START

    do i = 1, len(contents) + 1
      if (in_error(stack)) exit
      if (i<=len(contents)) then
        c = contents(i:i)
      else
        c = " "
      endif

      if (state==ST_START) then
        !write(*,*)'ST_START'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (isInitialNameChar(c, xv)) then
          name => vs_str_alloc(c)
          state = ST_NAME
        else
          call add_error(stack, &
            'Unexpected character in Attlist')
        endif

      elseif (state==ST_NAME) then
        !write(*,*)'ST_NAME'
        if (isNameChar(c, xv)) then
          temp => vs_str_alloc(str_vs(name)//c)
          deallocate(name)
          name => temp
        elseif (verify(c, XML_WHITESPACE)==0) then
          if (namespaces.and..not.checkQName(str_vs(name), xv)) then
            call add_error(stack, &
              "Attribute name in ATTLIST must be QName")
          elseif (associated(elem)) then
            if (existing_attribute(elem%attlist, str_vs(name))) then
              if (associated(ignore_att)) call destroy_attribute_t(ignore_att)
              allocate(ignore_att)
              call init_string_list(ignore_att%enumerations)
              ignore_att%name => vs_vs_alloc(name)
              ca => ignore_att
            else
              ca => add_attribute(elem%attlist, str_vs(name), internal)
            endif
          else
            if (associated(ignore_att)) call destroy_attribute_t(ignore_att)
            allocate(ignore_att)
            call init_string_list(ignore_att%enumerations)
            ignore_att%name => vs_vs_alloc(name)
            ca => ignore_att
          endif
          deallocate(name)
          state = ST_AFTERNAME
        else
          call add_error(stack, &
            'Unexpected character in Attlist Name')
        endif

      elseif (state==ST_AFTERNAME) then
        !write(*,*)'ST_AFTERNAME'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (verify(c, upperCase)==0) then
          attType => vs_str_alloc(c)
          state = ST_ATTTYPE
        elseif (c=='(') then
          allocate(value(0))
          ca%attType = ATT_ENUM
          state = ST_ENUMERATION
        else
          call add_error(stack, &
            'Unexpected error after Attlist Name')
        endif

      elseif (state==ST_ATTTYPE) then
        !write(*,*)'ST_ATTTYPE'
        if (verify(c, upperCase)==0) then
          temp => attType
          attType => vs_str_alloc(str_vs(temp)//c)
          deallocate(temp)
        elseif (verify(c, XML_WHITESPACE)==0) then
          ! xml:id constraint
          if (str_vs(ca%name)=="xml:id" &
            .and..not.str_vs(attType)=="ID") then
            call add_error(stack, &
              "xml:id attribute must be declared as type ID")
          elseif (str_vs(attType)=='CDATA') then
            ca%attType = ATT_CDATA
            state = ST_AFTER_ATTTYPE
          elseif (str_vs(attType)=='ID') then
            if (validCheck) then
              ! Validity Constraint: One ID per Element Type
              if (associated(elem)) then
                if (elem%id_declared) then
                  call add_error(stack, &
                    "Cannot have two declared attributes of type ID on one element type.")
                else
                  elem%id_declared = .true.
                endif
              endif
            endif
            ca%attType = ATT_ID
            state = ST_AFTER_ATTTYPE
          elseif (str_vs(attType)=='IDREF') then
            ca%attType = ATT_IDREF
            state = ST_AFTER_ATTTYPE
          elseif (str_vs(attType)=='IDREFS') then
            ca%attType = ATT_IDREFS
            state = ST_AFTER_ATTTYPE
          elseif (str_vs(attType)=='ENTITY') then
            ca%attType = ATT_ENTITY
            state = ST_AFTER_ATTTYPE
          elseif (str_vs(attType)=='ENTITIES') then
            ca%attType = ATT_ENTITIES
            state = ST_AFTER_ATTTYPE
          elseif (str_vs(attType)=='NMTOKEN') then
            ca%attType = ATT_NMTOKEN
            state = ST_AFTER_ATTTYPE
          elseif (str_vs(attType)=='NMTOKENS') then
            ca%attType = ATT_NMTOKENS
            state = ST_AFTER_ATTTYPE
          elseif (str_vs(attType)=='NOTATION') then
            ca%attType = ATT_NOTATION
            state = ST_AFTER_NOTATION
          else
            call add_error(stack, &
              'Unknown AttType')
          endif
          deallocate(attType)
        else
          call add_error(stack, &
            'Unexpected character in AttType')
        endif

      elseif (state==ST_AFTER_NOTATION) then
        !write(*,*)'ST_AFTER_NOTATION'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (c=='(') then
          state = ST_NOTATION_LIST
        else
          call add_error(stack, &
            'Unexpected character after Notation')
        endif

      elseif (state==ST_NOTATION_LIST) then
        !write(*,*)'ST_NOTATION_LIST'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (isInitialNameChar(c, xv)) then
          value => vs_str_alloc(c)
          state = ST_ENUM_NAME
        else
          call add_error(stack, &
            'Unexpected character in Notation list')
        endif

      elseif (state==ST_ENUMERATION) then
        !write(*,*)'ST_ENUMERATION'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (isNameChar(c, xv)) then
          temp => vs_str_alloc(str_vs(value)//c)
          deallocate(value)
          value => temp
          state = ST_ENUM_NAME
        elseif (c=='|') then
          call add_error(stack, &
            "Missing token in Enumeration")
        elseif (c==')') then
          call add_error(stack, &
            "Missing tokens in Enumeration")
        else
          call add_error(stack, &
            'Unexpected character in attlist enumeration')
        endif

      elseif (state==ST_ENUM_NAME) then
        !write(*,*)'ST_ENUM_NAME'
        if (isNameChar(c, xv)) then
          temp => vs_str_alloc(str_vs(value)//c)
          deallocate(value)
          value => temp
        elseif (verify(c, XML_WHITESPACE)==0) then
          if (validCheck.and.registered_string(ca%enumerations, str_vs(value))) then
            call add_error(stack, &
              "Duplicate enumeration value in ATTLIST")
          elseif (namespaces.and.ca%attType==ATT_NOTATION &
            .and..not.checkNCName(str_vs(value), xv)) then
            call add_error(stack, &
              "Notation name must be NCName")
          else
            call add_string(ca%enumerations, str_vs(value))
          endif
          deallocate(value)
          state = ST_SEPARATOR
        elseif (c=='|') then
          if (validCheck.and.registered_string(ca%enumerations, str_vs(value))) then
            call add_error(stack, &
              "Duplicate enumeration value in ATTLIST")
          elseif (namespaces.and.ca%attType==ATT_NOTATION &
            .and..not.checkNCName(str_vs(value), xv)) then
            call add_error(stack, &
              "Notation name must be NCName")
          else
            call add_string(ca%enumerations, str_vs(value))
          endif
          deallocate(value)
          if (ca%attType==ATT_NOTATION) then
            state = ST_NOTATION_LIST
          else
            allocate(value(0))
            state = ST_ENUMERATION
          endif
        elseif (c==')') then
          if (size(value)==0) then
            call add_error(stack, &
              'Missing token in Enumeration list')
          endif
          if (validCheck.and.registered_string(ca%enumerations, str_vs(value))) then
            call add_error(stack, &
              "Duplicate enumeration value in ATTLIST")
          elseif (namespaces.and.ca%attType==ATT_NOTATION &
            .and..not.checkNCName(str_vs(value), xv)) then
            call add_error(stack, &
              "Notation name must be NCName")
          else
            call add_string(ca%enumerations, str_vs(value))
          endif
          deallocate(value)
          state = ST_AFTER_ATTTYPE_SPACE
        else
          call add_error(stack, &
            'Unexpected character in attlist enumeration')
        endif

      elseif (state==ST_SEPARATOR) then
        !write(*,*)'ST_SEPARATOR'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (c=='|') then
          if (ca%attType==ATT_NOTATION) then
            state = ST_NOTATION_LIST
          else
            allocate(value(0))
            state = ST_ENUMERATION
          endif
        elseif (c==')') then
          state = ST_AFTER_ATTTYPE_SPACE
        else
          call add_error(stack, &
            'Unexpected character in attlist enumeration')
        endif

      elseif (state==ST_AFTER_ATTTYPE_SPACE) then
        if (verify(c, XML_WHITESPACE)/=0) then
          call add_error(stack, &
            'Missing whitespace in attlist enumeration')
        endif
        state = ST_AFTER_ATTTYPE

      elseif (state==ST_AFTER_ATTTYPE) then
        !write(*,*)'ST_AFTER_ATTTYPE'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (c=='#') then
          allocate(default(0))
          state = ST_DEFAULT_DECL
        elseif (c=='"'.or.c=="'") then
          if (validCheck) then
            ! Validity Constraint: ID Attribute Default
            if (ca%attType==ATT_ID) &
              call add_error(stack, &
              "Attribute of type ID may not have default value")
          endif
          ca%attDefault = ATT_DEFAULT
          q = c
          allocate(value(0))
          state = ST_DEFAULTVALUE
        else
          call add_error(stack, &
            'Unexpected character after AttType')
        endif

      elseif (state==ST_DEFAULT_DECL) then
        !write(*,*)'ST_DEFAULT_DECL'
        if (verify(c, upperCase)==0) then
          temp => vs_str_alloc(str_vs(default)//c)
          deallocate(default)
          default => temp
        elseif (verify(c, XML_WHITESPACE)==0) then
          if (str_vs(default)=='REQUIRED') then
            ca%attdefault = ATT_REQUIRED
            deallocate(default)
            state = ST_START
          elseif (str_vs(default)=='IMPLIED') then
            ca%attdefault = ATT_IMPLIED
            deallocate(default)
            state = ST_START
          elseif (str_vs(default)=='FIXED') then
            if (validCheck) then
              ! Validity Constraint: ID Attribute Default
              if (ca%attType==ATT_ID) &
                call add_error(stack, &
                "Attribute of type ID may not have FIXED value")
            endif
            ca%attdefault = ATT_FIXED
            deallocate(default)
            state = ST_AFTERDEFAULTDECL
          else
            call add_error(stack, &
              'Unknown Default declaration')
          endif
        else
          call add_error(stack, &
            'Unexpected character in Default declaration')
        endif

      elseif (state==ST_AFTERDEFAULTDECL) then
        !write(*,*)'ST_AFTERDEFAULTDECL'
        if (verify(c, XML_WHITESPACE)==0) cycle
        if (c=='"') then
          q = c
          allocate(value(0))
          state = ST_DEFAULTVALUE
        elseif (c=="'") then
          q = c
          allocate(value(0))
          state = ST_DEFAULTVALUE
        else
          call add_error(stack, &
            'Unexpected character after Default declaration')
        endif

      elseif (state==ST_DEFAULTVALUE) then
        !write(*,*)'ST_DEFAULTVALUE'
        if (c==q) then
          if (ca%attType/=ATT_CDATA) then
            temp => vs_str_alloc(att_value_normalize(str_vs(value)))
            deallocate(value)
            value => temp
          endif
          if (validCheck) then
            select case(ca%attType)
              ! Can't have ID with defaults
            case (ATT_IDREF)
              ! VC: IDREF
              if (namespaces) then
                if (.not.checkNCName(str_vs(value), xv)) &
                  call add_error(stack, &
                  "Attributes of type IDREF must have a value which is an XML NCName")
              else
                if (.not.checkName(str_vs(value), xv)) &
                  call add_error(stack, &
                  "Attributes of type IDREF must have a value which is an XML Name")
              endif
            case (ATT_IDREFS)
              ! VC: IDREF
              if (namespaces) then
                if (.not.checkNCNames(str_vs(value), xv)) &
                  call add_error(stack, &
                  "Attributes of type IDREFS must have a value which contains only XML NCNames")
              else
                if (.not.checkNames(str_vs(value), xv)) &
                  call add_error(stack, &
                  "Attributes of type IDREFS must have a value which contains only XML Names")
              endif
            case (ATT_ENTITY)
              ! VC: Entity Name
              if (namespaces) then
                if (.not.checkNCName(str_vs(value), xv)) &
                  call add_error(stack, &
                  "Attributes of type ENTITY must have a value which is an XML NCName")
              else
                if (.not.checkName(str_vs(value), xv)) &
                  call add_error(stack, &
                  "Attributes of type ENTITY must have a value which is an XML Name")
              endif
            case (ATT_ENTITIES)
              ! VC: Entity Name
              if (namespaces) then
                if (.not.checkNames(str_vs(value), xv)) &
                  call add_error(stack, &
                  "Attributes of type ENTITIES must have a value which contains only XML NCNames")
              else
                if (.not.checkNames(str_vs(value), xv)) &
                  call add_error(stack, &
                  "Attributes of type ENTITIES must have a value which contains only XML Names")
              endif
            case (ATT_NMTOKEN)
              ! VC Name Token
              if (.not.checkNmtoken(str_vs(value), xv)) &
                call add_error(stack, &
                "Attributes of type NMTOKEN must have a value which is a NMTOKEN")
            case (ATT_NMTOKENS)
              ! VC: Name Token
              if (.not.checkNmtokens(str_vs(value), xv)) &
                call add_error(stack, &
                "Attributes of type NMTOKENS must have a value which contain only NMTOKENs")
            case (ATT_NOTATION)
              ! VC: Notation Attributes
              if (namespaces) then
                if (.not.checkNCName(str_vs(value), xv)) &
                  call add_error(stack, &
                  "Attributes of type NOTATION must have a value which is an XMLNCName")
              else
                if (.not.checkName(str_vs(value), xv)) &
                  call add_error(stack, &
                  "Attributes of type NOTATION must have a value which is an XML Name")
              endif
            case (ATT_ENUM)
              ! VC: Enumeration
              if (.not.checkNmtoken(str_vs(value), xv)) &
                call add_error(stack, &
                "Attributes of type ENUM must have a value which is an NMTOKENs")
              if (.not.registered_string(ca%enumerations, str_vs(value))) &
                call add_error(stack, &
                "Default value of ENUM does not match permitted values")
            end select
          endif
          if (.not.in_error(stack)) then
            if (ca%attType==ATT_ENTITIES) then
              call destroy_string_list(ca%enumerations)
              ca%enumerations = tokenize_to_string_list(str_vs(value))
            endif
            ca%default => value
            value => null()
            state = ST_START
          endif
        else
          temp => vs_str_alloc(str_vs(value)//c)
          deallocate(value)
          value => temp
        endif

      endif

    enddo

    if (associated(ignore_att)) call destroy_attribute_t(ignore_att)

    if (.not.in_error(stack)) then
      if (state==ST_START) then
        return
      else
        call add_error(stack, &
          'Incomplete Attlist declaration')
      endif
    endif
    
    if (associated(name)) deallocate(name)
    if (associated(attType)) deallocate(attType)
    if (associated(default)) deallocate(default)
    if (associated(value)) deallocate(value)

  end subroutine parse_dtd_attlist

  subroutine report_declarations(elem, attributeDecl_handler)
    type(element_t), intent(in) :: elem
    interface
      subroutine attributeDecl_handler(eName, aName, type, mode, value)
        character(len=*), intent(in) :: eName
        character(len=*), intent(in) :: aName
        character(len=*), intent(in) :: type
        character(len=*), intent(in), optional :: mode
        character(len=*), intent(in), optional :: value
      end subroutine attributeDecl_handler
    end interface

    integer :: i
    character(len=8) :: type
    character(len=8) :: mode
    type(attribute_t), pointer :: a

    do i = 1, size(elem%attlist%list)
      a => elem%attlist%list(i)
      type = ATT_TYPES(a%attType)
      select case (a%attDefault)
      case (ATT_REQUIRED)
        mode = "REQUIRED"
      case (ATT_IMPLIED)
        mode = "IMPLIED"
      case (ATT_FIXED)
        mode = "FIXED"
      end select

      if (a%attType==ATT_NOTATION) then
        if (a%attDefault==ATT_DEFAULT) then
          if (associated(a%default)) then
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              'NOTATION '//make_token_group(a%enumerations), value=str_vs(a%default))
          else
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              'NOTATION '//make_token_group(a%enumerations))
          endif
        else
          if (associated(a%default)) then
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              'NOTATION '//make_token_group(a%enumerations), mode=trim(mode), &
              value=str_vs(a%default))
          else
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              'NOTATION '//make_token_group(a%enumerations), mode=trim(mode))
          endif
        endif
      elseif (a%attType==ATT_ENUM) then
        if (a%attDefault==ATT_DEFAULT) then
          if (associated(a%default)) then
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              make_token_group(a%enumerations), value=str_vs(a%default))
          else
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              make_token_group(a%enumerations))
          endif
        else
          if (associated(a%default)) then
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              make_token_group(a%enumerations), mode=trim(mode), &
              value=str_vs(a%default))
          else
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              make_token_group(a%enumerations), mode=trim(mode))
          endif
        endif
      else
        if (a%attDefault==ATT_DEFAULT) then
          if (associated(a%default)) then
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              trim(type), value=str_vs(a%default))
          else
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              trim(type))
          endif
        else
          if (associated(a%default)) then
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              trim(type), mode=trim(mode), value=str_vs(a%default))
          else
            call attributeDecl_handler(str_vs(elem%name), str_vs(a%name), &
              trim(type), mode=trim(mode))
          endif
        endif
      endif
    enddo


  end subroutine report_declarations

  pure function make_token_group_len(s_list) result(n)
    type(string_list), intent(in) :: s_list
    integer :: n

    integer :: i
    n = size(s_list%list) + 1
    do i = 1, size(s_list%list)
      n = n + size(s_list%list(i)%s)
    enddo
  end function make_token_group_len

  function make_token_group(s_list) result(s)
    type(string_list), intent(in) :: s_list
    character(len=make_token_group_len(s_list)) :: s
    
    integer :: i, m, n
    s(1:1) = '('
    n = 2
    do i = 1, size(s_list%list)-1
      m = size(s_list%list(i)%s)
      s(n:n+m) = str_vs(s_list%list(i)%s)//'|'
      n = n + m + 1
    enddo
    s(n:) = str_vs(s_list%list(i)%s)//')'
  end function make_token_group

  function attribute_has_default(att) result(p)
    type(attribute_t), pointer :: att
    logical :: p

    if (associated(att)) then
      p = att%attDefault==ATT_DEFAULT.or.att%attDefault==ATT_FIXED
    else
      p = .false.
    endif
  end function attribute_has_default

  function get_attlist_size(elem) result(n)
    type(element_t), pointer :: elem
    integer :: n

    if (associated(elem)) then
      n = size(elem%attlist%list)
    else
      n = 0
    endif
  end function get_attlist_size

  function get_attdecl_by_index(elem, n) result(att)
    type(element_t), pointer :: elem
    integer, intent(in) :: n
    type(attribute_t), pointer :: att

    att => null()
    if (associated(elem)) then
      if (n>0.and.n<=size(elem%attlist%list)) then
        att => elem%attlist%list(n)
      endif
    endif
  end function get_attdecl_by_index

  function get_attdecl_by_name(elem, name) result(att)
    type(element_t), pointer :: elem
    character(len=*), intent(in) :: name
    type(attribute_t), pointer :: att

    integer :: i
    att => null()
    if (associated(elem)) then
      do i = 1, size(elem%attlist%list)
        if (str_vs(elem%attlist%list(i)%name)==name) then
          att => elem%attlist%list(i)
          return
        endif
      enddo
    endif
  end function get_attdecl_by_name

  pure function express_att_decl_len(a) result(n)
    type(attribute_t), intent(in) :: a
    integer :: n

    if (a%attType==ATT_ENUM) then
      n = size(a%name)
    else
      n = size(a%name)+1+len_trim(ATT_TYPES(a%attType))
    endif

    if (a%attType==ATT_NOTATION &
      .or.a%attType==ATT_ENUM) &
      n = n + 1 + make_token_group_len(a%enumerations)

    select case(a%attDefault)
    case (ATT_REQUIRED)
      n = n + len(" #REQUIRED")
    case (ATT_IMPLIED)
      n = n + len(" #IMPLIED")
    case (ATT_DEFAULT)
      n = n + len(" ")
    case (ATT_FIXED)
      n = n + len(" #FIXED")
    end select
    
    if (associated(a%default)) &
      n = n + 3 + size(a%default)
  end function express_att_decl_len

  function express_attribute_declaration(a) result(s)
    type(attribute_t), intent(in) :: a
    character(len=express_att_decl_len(a)) :: s

    if (a%attType==ATT_ENUM) then
      s = str_vs(a%name)
    else
      s = str_vs(a%name)//" "//ATT_TYPES(a%attType)
    endif
    if (a%attType==ATT_NOTATION &
      .or.a%attType==ATT_ENUM) &
      s = trim(s)//" "//make_token_group(a%enumerations)

    select case(a%attDefault)
    case (ATT_REQUIRED)
      s = trim(s)//" #REQUIRED"
    case (ATT_IMPLIED)
      s = trim(s)//" #IMPLIED"
    case (ATT_DEFAULT)
      s = trim(s)//" "
    case (ATT_FIXED)
      s = trim(s)//" #FIXED"
    end select
    
    if (associated(a%default)) &
      s = trim(s)//" """//str_vs(a%default)//""""
  end function express_attribute_declaration

  function get_att_type_enum(s) result(n)
    character(len=*), intent(in) :: s
    integer :: n

    select case(s)
    case ('CDATA')
      n = ATT_CDATA
    case ('ID')
      n = ATT_ID
    case ('IDREF')
      n = ATT_IDREF
    case ('IDREFS')
      n = ATT_IDREFS
    case ('NMTOKEN')
      n = ATT_NMTOKEN
    case ('NMTOKENS')
      n = ATT_NMTOKENS
    case ('ENTITY')
      n = ATT_ENTITY
    case ('ENTITIES')
      n = ATT_ENTITIES
    case ('NOTATION')
      n = ATT_NOTATION
    case ('CDANO')
      n= ATT_CDANO
    case ('CDAMB')
      n = ATT_CDAMB
    end select
  end function get_att_type_enum

  pure function att_value_normalize_len(s1) result(n)
    character(len=*), intent(in) :: s1
    integer :: n

    integer :: i
    logical :: w

    n = 0
    w = .true.
    do i = 1, len(s1)
      if (w.and.(verify(s1(i:i),XML_WHITESPACE)==0)) cycle
      w = .false.
      n = n + 1
      if (verify(s1(i:i),XML_WHITESPACE)==0) w = .true.
    enddo
    if (w) n = n - 1 ! Discard final space
    
  end function att_value_normalize_len

  function att_value_normalize(s1) result(s2)
    character(len=*), intent(in) :: s1
    character(len=att_value_normalize_len(s1)) :: s2
    
    integer :: i, i2
    logical :: w
    
    i = 0
    i2 = 1
    w = .true.
    do while (i2<=len(s2))
      i = i + 1
      if (w.and.(verify(s1(i:i),XML_WHITESPACE)==0)) cycle
      w = .false.
      s2(i2:i2) = s1(i:i)
      i2 = i2 + 1
      if (verify(s1(i:i),XML_WHITESPACE)==0) w = .true.
    enddo

  end function att_value_normalize

#endif
end module m_common_element
