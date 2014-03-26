module m_dom_error

  use fox_m_fsys_abort_flush, only: pxfabort

  use m_common_error, only: error_stack, add_error, in_error, destroy_error_stack

  implicit none
  private

  type DOMException
    private
    type(error_stack) :: stack
  end type DOMException

  integer, parameter, public :: INDEX_SIZE_ERR              = 1
  integer, parameter, public :: DOMSTRING_SIZE_ERR          = 2
  integer, parameter, public :: HIERARCHY_REQUEST_ERR       = 3
  integer, parameter, public :: WRONG_DOCUMENT_ERR          = 4
  integer, parameter, public :: INVALID_CHARACTER_ERR       = 5
  integer, parameter, public :: NO_DATA_ALLOWED_ERR         = 6
  integer, parameter, public :: NO_MODIFICATION_ALLOWED_ERR = 7
  integer, parameter, public :: NOT_FOUND_ERR               = 8
  integer, parameter, public :: NOT_SUPPORTED_ERR           = 9
  integer, parameter, public :: INUSE_ATTRIBUTE_ERR         = 10
  integer, parameter, public :: INVALID_STATE_ERR           = 11
  integer, parameter, public :: SYNTAX_ERR                  = 12
  integer, parameter, public :: INVALID_MODIFICATION_ERR    = 13
  integer, parameter, public :: NAMESPACE_ERR               = 14
  integer, parameter, public :: INVALID_ACCESS_ERR          = 15
  integer, parameter, public :: VALIDATION_ERR              = 16
  integer, parameter, public :: TYPE_MISMATCH_ERR           = 17

  integer, parameter, public :: INVALID_EXPRESSION_ERR      = 51
  integer, parameter, public :: TYPE_ERR                    = 52

  integer, parameter, public :: PARSE_ERR                   = 81
  integer, parameter, public :: SERIALIZE_ERR               = 82

  integer, parameter, public :: FoX_INVALID_NODE            = 201
  integer, parameter, public :: FoX_INVALID_CHARACTER       = 202
  integer, parameter, public :: FoX_NO_SUCH_ENTITY          = 203
  integer, parameter, public :: FoX_INVALID_PI_DATA         = 204
  integer, parameter, public :: FoX_INVALID_CDATA_SECTION   = 205
  integer, parameter, public :: FoX_HIERARCHY_REQUEST_ERR   = 206
  integer, parameter, public :: FoX_INVALID_PUBLIC_ID       = 207
  integer, parameter, public :: FoX_INVALID_SYSTEM_ID       = 208
  integer, parameter, public :: FoX_INVALID_COMMENT         = 209
  integer, parameter, public :: FoX_NODE_IS_NULL            = 210
  integer, parameter, public :: FoX_INVALID_ENTITY          = 211
  integer, parameter, public :: FoX_INVALID_URI             = 212
  integer, parameter, public :: FoX_IMPL_IS_NULL            = 213
  integer, parameter, public :: FoX_MAP_IS_NULL             = 214
  integer, parameter, public :: FoX_LIST_IS_NULL            = 215

  integer, parameter, public :: FoX_INTERNAL_ERROR          = 999


  public :: DOMException
  public :: getExceptionCode
  public :: throw_exception
  public :: inException
  public :: dom_error
  public :: internal_error

contains

  pure function errorString(code) result(s)
    integer, intent(in) :: code
    character(len=27) :: s

    select case(code)
    case(1)
      s = "INDEX_SIZE_ERR"
    case(2)
      s = "DOMSTRING_SIZE_ERR"
    case(3)
      s = "HIERARCHY_REQUEST_ERR"
    case(4)
      s = "WRONG_DOCUMENT_ERR"
    case(5)
      s = "INVALID_CHARACTER_ERR"
    case(6)
      s = "NO_DATA_ALLOWED_ERR"
    case(7)
      s = "NO_MODIFICATION_ALLOWED_ERR"
    case(8)
      s = "NOT_FOUND_ERR"
    case(9)
      s = "NOT_SUPPORTED_ERR"
    case(10)
      s = "INUSE_ATTRIBUTE_ERR"
    case(11)
      s = "INVALID_STATE_ERR"
    case(12)
      s = "SYNTAX_ERR"
    case(13)
      s = "INVALID_MODIFICATION_ERR"
    case(14)
      s = "NAMESPACE_ERR"
    case(15)
      s = "INVALID_ACCESS_ERR"
    case(16)
      s = "VALIDATION_ERR"
    case(18)
      s = "TYPE_MISMATCH_ERR"
    case(51)
      s = "INVALID_EXPRESSION_ERR"
    case(52)
      s = "TYPE_ERR"
    case(81)
      s = "PARSE_ERR"
    case(82)
      s = "SERIALIZE_ERR"
    case(201)
      s = "FoX_INVALID_NODE"
    case(202)
      s = "FoX_INVALID_CHARACTER"
    case(203)
      s = "FoX_NO_SUCH_ENTITY"
    case(204)
      s = "FoX_INVALID_PI_DATA"
    case(205)
      s = "FoX_INVALID_CDATA_SECTION"
    case(206)
      s = "FoX_HIERARCHY_REQUEST_ERR"
    case(207)
      s = "FoX_INVALID_PUBLIC_ID"
    case(208)
      s = "FoX_INVALID_SYSTEM_ID"
    case(209)
      s = "FoX_INVALID_COMMENT"
    case(210)
      s = "FoX_NODE_IS_NULL"
    case(211)
      s = "FoX_INVALID_ENTITY"
    case(212)
      s = "FoX_NO_DOCTYPE"
    case(213)
      s = "FoX_IMPL_IS_NULL"
    case(214)
      s = "FoX_MAP_IS_NULL"
    case(215)
      s = "FoX_LIST_IS_NULL"
    case default
      s = "INTERNAL ERROR!!!!"
    end select

  end function errorString

  function getExceptionCode(ex) result(n)
    type(DOMException), intent(inout) :: ex
    integer :: n

    if (in_error(ex%stack)) then
      n = ex%stack%stack(size(ex%stack%stack))%error_code
      call destroy_error_stack(ex%stack)
    else
      n = 0
    endif
  end function getExceptionCode

  subroutine throw_exception(code, msg, ex)
    integer, intent(in) :: code
    character(len=*), intent(in) :: msg
    type(DOMException), intent(inout), optional :: ex

    if (present(ex)) then
      call add_error(ex%stack, msg, error_code=code) ! FIXME
    else
      write(0,'(a)') errorString(code)
      write(0,'(i0,a)') code, " "//msg
      call pxfabort()
    endif
      
  end subroutine throw_exception


  subroutine dom_error(name,code,msg)
    character(len=*), intent(in) :: name, msg
    integer, intent(in)          :: code

    write(0,'(4a)') "Routine ", name, ":", msg
    write(0,'(a)') errorString(code)
    call pxfabort()

  end subroutine dom_error

  
  subroutine destroyDOMException(ex)
    type(DOMException), intent(inout) :: ex

    call destroy_error_stack(ex%stack)
  end subroutine destroyDOMException


  subroutine internal_error(name,msg)
    character(len=*), intent(in) :: name, msg

    write(0,'(4a)') "Internal error in ", name, ":", msg
    call pxfabort()

  end subroutine internal_error

  function inException(ex) result(p)
    type(DOMException), intent(in) :: ex
    logical :: p

    p = in_error(ex%stack)
  end function inException

end module m_dom_error
