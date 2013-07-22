module m_sax_types

#ifndef DUMMYLIB
  use m_common_attrs, only: dictionary_t
  use m_common_elstack, only: elstack_t
  use m_common_entities, only: entity_list
  use m_common_error, only: error_stack
  use m_common_namespaces, only: namespacedictionary
  use m_common_notations, only: notation_list
  use m_common_struct, only: xml_doc_state

  use m_sax_reader, only: file_buffer_t

  implicit none

  ! Context

  integer, parameter :: CTXT_NULL = -1
  integer, parameter :: CTXT_INIT = 0
  integer, parameter :: CTXT_BEFORE_DTD = 1
  integer, parameter :: CTXT_IN_DTD = 2
  integer, parameter :: CTXT_IGNORE = 3
  integer, parameter :: CTXT_BEFORE_CONTENT = 4
  integer, parameter :: CTXT_IN_CONTENT = 5
  integer, parameter :: CTXT_AFTER_CONTENT = 6

  ! State

  integer, parameter :: ST_STOP                     = -1
  integer, parameter :: ST_NULL                     = 0
  integer, parameter :: ST_MISC                     = 1
  integer, parameter :: ST_BANG_TAG                 = 2
  integer, parameter :: ST_START_PI                 = 3
  integer, parameter :: ST_PI_CONTENTS              = 4
  integer, parameter :: ST_PI_END                   = 5
  integer, parameter :: ST_START_COMMENT            = 6
  integer, parameter :: ST_COMMENT_END              = 7
  integer, parameter :: ST_START_TAG                = 8
  integer, parameter :: ST_START_CDATA_DECLARATION  = 9
  integer, parameter :: ST_FINISH_CDATA_DECLARATION = 10
  integer, parameter :: ST_IN_TAG                   = 11
  integer, parameter :: ST_ATT_NAME                 = 12
  integer, parameter :: ST_ATT_EQUALS               = 13
  integer, parameter :: ST_CHAR_IN_CONTENT          = 14
  integer, parameter :: ST_CLOSING_TAG              = 15
  integer, parameter :: ST_CDATA_CONTENTS           = 16
  integer, parameter :: ST_IN_CLOSING_TAG           = 17
  integer, parameter :: ST_TAG_IN_CONTENT           = 18
  integer, parameter :: ST_CDATA_END                = 19
  integer, parameter :: ST_IN_DOCTYPE               = 20
  integer, parameter :: ST_DOC_NAME                 = 21
  integer, parameter :: ST_DOC_SYSTEM               = 22
  integer, parameter :: ST_DOC_PUBLIC               = 23
  integer, parameter :: ST_DOC_DECL                 = 24
  integer, parameter :: ST_CLOSE_DOCTYPE            = 25
  integer, parameter :: ST_START_ENTITY             = 26
  integer, parameter :: ST_START_PE                 = 27
  integer, parameter :: ST_IN_SUBSET                = 28

! DTD states
  integer, parameter :: ST_DTD_NULL                = 50
  integer, parameter :: ST_DTD_SUBSET              = 51
  integer, parameter :: ST_DTD_START_SECTION_DECL  = 52
  integer, parameter :: ST_DTD_FINISH_SECTION_DECL = 53
  integer, parameter :: ST_DTD_IN_IGNORE_SECTION   = 54
  integer, parameter :: ST_DTD_BANG_TAG            = 55
  integer, parameter :: ST_DTD_START_PI            = 56
  integer, parameter :: ST_DTD_PI_CONTENTS         = 57
  integer, parameter :: ST_DTD_PI_END              = 58
  integer, parameter :: ST_DTD_COMMENT_END         = 59
  integer, parameter :: ST_DTD_START_COMMENT       = 60
  integer, parameter :: ST_DTD_ATTLIST             = 61
  integer, parameter :: ST_DTD_ELEMENT             = 62
  integer, parameter :: ST_DTD_ENTITY              = 63
  integer, parameter :: ST_DTD_NOTATION            = 64
  integer, parameter :: ST_DTD_NOTATION_ID         = 65
  integer, parameter :: ST_DTD_NOTATION_SYSTEM     = 66
  integer, parameter :: ST_DTD_NOTATION_PUBLIC     = 67
  integer, parameter :: ST_DTD_NOTATION_PUBLIC_2   = 68
  integer, parameter :: ST_DTD_NOTATION_END        = 69
  integer, parameter :: ST_DTD_ENTITY_PE           = 70
  integer, parameter :: ST_DTD_ENTITY_ID           = 71
  integer, parameter :: ST_DTD_ENTITY_PUBLIC       = 72
  integer, parameter :: ST_DTD_ENTITY_SYSTEM       = 73
  integer, parameter :: ST_DTD_ENTITY_NDATA        = 74
  integer, parameter :: ST_DTD_ENTITY_NDATA_VALUE  = 75
  integer, parameter :: ST_DTD_ENTITY_END          = 76
  integer, parameter :: ST_DTD_ATTLIST_CONTENTS    = 77
  integer, parameter :: ST_DTD_ATTLIST_END         = 78
  integer, parameter :: ST_DTD_ELEMENT_CONTENTS    = 79
  integer, parameter :: ST_DTD_ELEMENT_END         = 80
  integer, parameter :: ST_DTD_DONE                = 81

! token types

  integer, parameter :: TOK_NULL = 0
  integer, parameter :: TOK_PI_TAG = 1 ! <?
  integer, parameter :: TOK_BANG_TAG = 2 ! <!
  integer, parameter :: TOK_OPEN_TAG = 3 ! <
  integer, parameter :: TOK_OPEN_SB = 4 ! [
  integer, parameter :: TOK_CLOSE_SB = 5 ! [
  integer, parameter :: TOK_OPEN_COMMENT = 6 ! --
  integer, parameter :: TOK_NAME = 7 ! name (+token)
  integer, parameter :: TOK_CHAR = 8 ! character data (+token)
  integer, parameter :: TOK_PI_END = 9 ! ?>
  integer, parameter :: TOK_COMMENT_END = 10 ! -->
  integer, parameter :: TOK_SECTION_START = 11 ! <![
  integer, parameter :: TOK_SECTION_END = 12 ! ]]>
  integer, parameter :: TOK_END_TAG = 13 ! >
  integer, parameter :: TOK_END_TAG_CLOSE = 14 ! />
  integer, parameter :: TOK_CLOSE_TAG = 15 ! </
  integer, parameter :: TOK_ENTITY = 16 ! % or &
  integer, parameter :: TOK_EQUALS = 17 ! =
  integer, parameter :: TOK_DTD_CONTENTS = 18 ! for element and attlist
  integer, parameter :: TOK_OPEN_PAR = 19 ! (
  integer, parameter :: TOK_CLOSE_PAR = 20 ! )

  type sax_parser_t
    type(xml_doc_state), pointer :: xds
    logical :: xds_used = .false. ! is the xds used by DOM? If so, we must
                                  ! not destroy it once we are finished
    integer :: context 
    integer :: state = ST_NULL
    integer :: state_dtd = ST_DTD_SUBSET
    logical :: well_formed = .false.
    logical :: skippedExternal = .false.
    character, dimension(:), pointer :: token => null()
    character, dimension(:), pointer :: content => null()
    integer :: tokenType = TOK_NULL
    integer :: nextTokenType = TOK_NULL
    character, dimension(:), pointer :: name => null()
    character, dimension(:), pointer :: attname => null()
    logical :: error = .false.
    type(error_stack) :: error_stack
    ! Aspects of document structure
    character, dimension(:), pointer :: root_element => null()
    type(elstack_t) :: elstack
    type(dictionary_t) :: attributes
    type(namespacedictionary) :: nsdict
    type(notation_list) :: nlist
    type(entity_list) :: predefined_e_list
    type(entity_list) :: forbidden_pe_list
    type(entity_list) :: forbidden_ge_list
    character(len=1), dimension(:), pointer :: PublicId => null()
    character(len=1), dimension(:), pointer :: SystemId => null()
    character(len=1), dimension(:), pointer :: Ndata => null()
    logical :: inIntSubset = .false.
    logical :: spaceBeforeEntity = .false.
  end type sax_parser_t
#endif

  type xml_t
#ifndef DUMMYLIB
    type(file_buffer_t) :: fb
    type(sax_parser_t) :: fx
#else
    integer :: i = 0
#endif
  end type xml_t

end module m_sax_types
