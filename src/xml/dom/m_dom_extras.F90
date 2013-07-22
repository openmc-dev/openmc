module m_dom_extras

  use fox_m_fsys_realtypes, only: sp, dp
  use fox_m_fsys_parse_input, only: rts

  use m_dom_error, only: DOMException, inException, throw_exception,           &
    FoX_NODE_IS_NULL, FoX_INVALID_NODE
  use m_dom_dom, only: Node, ELEMENT_NODE,                                     &
    getAttribute, getAttributeNS, getTextContent, getNodeType, getFoX_checks

  implicit none
  private

  public :: extractDataContent
  public :: extractDataAttribute
  public :: extractDataAttributeNS

  interface extractDataContent
     module procedure extractDataContentCmplxDpSca
     module procedure extractDataContentCmplxDpArr
     module procedure extractDataContentCmplxDpMat
     module procedure extractDataContentCmplxSpSca
     module procedure extractDataContentCmplxSpArr
     module procedure extractDataContentCmplxSpMat
     module procedure extractDataContentRealDpSca
     module procedure extractDataContentRealDpArr
     module procedure extractDataContentRealDpMat
     module procedure extractDataContentRealSpSca
     module procedure extractDataContentRealSpArr
     module procedure extractDataContentRealSpMat
     module procedure extractDataContentIntSca
     module procedure extractDataContentIntArr
     module procedure extractDataContentIntMat
     module procedure extractDataContentLgSca
     module procedure extractDataContentLgArr
     module procedure extractDataContentLgMat
     module procedure extractDataContentChSca
     module procedure extractDataContentChArr
     module procedure extractDataContentChMat

  end interface extractDataContent

  interface extractDataAttribute
     module procedure extractDataAttributeCmplxDpSca
     module procedure extractDataAttributeCmplxDpArr
     module procedure extractDataAttributeCmplxDpMat
     module procedure extractDataAttributeCmplxSpSca
     module procedure extractDataAttributeCmplxSpArr
     module procedure extractDataAttributeCmplxSpMat
     module procedure extractDataAttributeRealDpSca
     module procedure extractDataAttributeRealDpArr
     module procedure extractDataAttributeRealDpMat
     module procedure extractDataAttributeRealSpSca
     module procedure extractDataAttributeRealSpArr
     module procedure extractDataAttributeRealSpMat
     module procedure extractDataAttributeIntSca
     module procedure extractDataAttributeIntArr
     module procedure extractDataAttributeIntMat
     module procedure extractDataAttributeLgSca
     module procedure extractDataAttributeLgArr
     module procedure extractDataAttributeLgMat
     module procedure extractDataAttributeChSca
     module procedure extractDataAttributeChArr
     module procedure extractDataAttributeChMat

  end interface extractDataAttribute

  interface extractDataAttributeNS
     module procedure extractDataAttNSCmplxDpSca
     module procedure extractDataAttNSCmplxDpArr
     module procedure extractDataAttNSCmplxDpMat
     module procedure extractDataAttNSCmplxSpSca
     module procedure extractDataAttNSCmplxSpArr
     module procedure extractDataAttNSCmplxSpMat
     module procedure extractDataAttNSRealDpSca
     module procedure extractDataAttNSRealDpArr
     module procedure extractDataAttNSRealDpMat
     module procedure extractDataAttNSRealSpSca
     module procedure extractDataAttNSRealSpArr
     module procedure extractDataAttNSRealSpMat
     module procedure extractDataAttNSIntSca
     module procedure extractDataAttNSIntArr
     module procedure extractDataAttNSIntMat
     module procedure extractDataAttNSLgSca
     module procedure extractDataAttNSLgArr
     module procedure extractDataAttNSLgMat
     module procedure extractDataAttNSChSca
     module procedure extractDataAttNSChArr
     module procedure extractDataAttNSChMat

  end interface extractDataAttributeNS

contains

subroutine extractDataContentCmplxDpSca(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    complex(dp), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentCmplxDpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentCmplxDpSca

subroutine extractDataContentCmplxSpSca(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    complex(sp), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentCmplxSpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentCmplxSpSca

subroutine extractDataContentRealDpSca(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    real(dp), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentRealDpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentRealDpSca

subroutine extractDataContentRealSpSca(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    real(sp), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentRealSpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentRealSpSca

subroutine extractDataContentIntSca(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    integer, intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentIntSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentIntSca

subroutine extractDataContentLgSca(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    logical, intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentLgSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentLgSca

subroutine extractDataContentChSca(arg, data, separator, csv, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(out) :: data

    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator
    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentChSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    endif
    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, separator, csv, num, iostat)
    else
      call rts(getTextContent(arg), data, separator, csv, num, iostat)
    endif

  end subroutine extractDataContentChSca


subroutine extractDataContentCmplxDpArr(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    complex(dp), dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentCmplxDpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentCmplxDpArr

subroutine extractDataContentCmplxSpArr(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    complex(sp), dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentCmplxSpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentCmplxSpArr

subroutine extractDataContentRealDpArr(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    real(dp), dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentRealDpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentRealDpArr

subroutine extractDataContentRealSpArr(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    real(sp), dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentRealSpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentRealSpArr

subroutine extractDataContentIntArr(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    integer, dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentIntArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentIntArr

subroutine extractDataContentLgArr(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    logical, dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentLgArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentLgArr

subroutine extractDataContentChArr(arg, data, separator, csv, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), dimension(:), intent(out) :: data

    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator
    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentChArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    endif
    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, separator, csv, num, iostat)
    else
      call rts(getTextContent(arg), data, separator, csv, num, iostat)
    endif

  end subroutine extractDataContentChArr


subroutine extractDataContentCmplxDpMat(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    complex(dp), dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentCmplxDpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentCmplxDpMat

subroutine extractDataContentCmplxSpMat(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    complex(sp), dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentCmplxSpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentCmplxSpMat

subroutine extractDataContentRealDpMat(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    real(dp), dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentRealDpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentRealDpMat

subroutine extractDataContentRealSpMat(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    real(sp), dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentRealSpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentRealSpMat

subroutine extractDataContentIntMat(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    integer, dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentIntMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentIntMat

subroutine extractDataContentLgMat(arg, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    logical, dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentLgMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif

  end subroutine extractDataContentLgMat

subroutine extractDataContentChMat(arg, data, separator, csv, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), dimension(:,:), intent(out) :: data

    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator
    integer, intent(out), optional :: num, iostat

    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataContentChMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    endif
    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, separator, csv, num, iostat)
    else
      call rts(getTextContent(arg), data, separator, csv, num, iostat)
    endif

  end subroutine extractDataContentChMat


subroutine extractDataAttributeCmplxDpSca(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    complex(dp), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeCmplxDpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeCmplxDpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeCmplxDpSca

subroutine extractDataAttributeCmplxSpSca(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    complex(sp), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeCmplxSpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeCmplxSpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeCmplxSpSca

subroutine extractDataAttributeRealDpSca(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    real(dp), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeRealDpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeRealDpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeRealDpSca

subroutine extractDataAttributeRealSpSca(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    real(sp), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeRealSpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeRealSpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeRealSpSca

subroutine extractDataAttributeIntSca(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    integer, intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeIntSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeIntSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeIntSca

subroutine extractDataAttributeLgSca(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    logical, intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeLgSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeLgSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeLgSca

subroutine extractDataAttributeChSca(arg, name, data, separator, csv, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator
    character(len=*), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeChSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeChSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, separator, csv, num, iostat)
    else
      call rts(getAttribute(arg, name), data, separator, csv, num, iostat)
    endif


  end subroutine extractDataAttributeChSca


subroutine extractDataAttributeCmplxDpArr(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    complex(dp), dimension(:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeCmplxDpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeCmplxDpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeCmplxDpArr

subroutine extractDataAttributeCmplxSpArr(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    complex(sp), dimension(:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeCmplxSpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeCmplxSpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeCmplxSpArr

subroutine extractDataAttributeRealDpArr(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    real(dp), dimension(:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeRealDpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeRealDpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeRealDpArr

subroutine extractDataAttributeRealSpArr(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    real(sp), dimension(:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeRealSpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeRealSpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeRealSpArr

subroutine extractDataAttributeIntArr(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    integer, dimension(:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeIntArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeIntArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeIntArr

subroutine extractDataAttributeLgArr(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    logical, dimension(:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeLgArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeLgArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeLgArr

subroutine extractDataAttributeChArr(arg, name, data, separator, csv, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator
    character(len=*), dimension(:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeChArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeChArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, separator, csv, num, iostat)
    else
      call rts(getAttribute(arg, name), data, separator, csv, num, iostat)
    endif


  end subroutine extractDataAttributeChArr


subroutine extractDataAttributeCmplxDpMat(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    complex(dp), dimension(:,:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeCmplxDpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeCmplxDpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeCmplxDpMat

subroutine extractDataAttributeCmplxSpMat(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    complex(sp), dimension(:,:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeCmplxSpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeCmplxSpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeCmplxSpMat

subroutine extractDataAttributeRealDpMat(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    real(dp), dimension(:,:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeRealDpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeRealDpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeRealDpMat

subroutine extractDataAttributeRealSpMat(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    real(sp), dimension(:,:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeRealSpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeRealSpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeRealSpMat

subroutine extractDataAttributeIntMat(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    integer, dimension(:,:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeIntMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeIntMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeIntMat

subroutine extractDataAttributeLgMat(arg, name, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    logical, dimension(:,:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeLgMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeLgMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif


  end subroutine extractDataAttributeLgMat

subroutine extractDataAttributeChMat(arg, name, data, separator, csv, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator
    character(len=*), dimension(:,:), intent(out) :: data
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttributeChMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttributeChMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, separator, csv, num, iostat)
    else
      call rts(getAttribute(arg, name), data, separator, csv, num, iostat)
    endif


  end subroutine extractDataAttributeChMat


subroutine extractDataAttNSCmplxDpSca(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    complex(dp), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSCmplxDpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSCmplxDpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSCmplxDpSca

subroutine extractDataAttNSCmplxSpSca(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    complex(sp), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSCmplxSpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSCmplxSpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSCmplxSpSca

subroutine extractDataAttNSRealDpSca(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    real(dp), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSRealDpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSRealDpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSRealDpSca

subroutine extractDataAttNSRealSpSca(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    real(sp), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSRealSpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSRealSpSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSRealSpSca

subroutine extractDataAttNSIntSca(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    integer, intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSIntSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSIntSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSIntSca

subroutine extractDataAttNSLgSca(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    logical, intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSLgSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSLgSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSLgSca

subroutine extractDataAttNSChSca(arg, namespaceURI, localName, data, separator, csv, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    character(len=*), intent(out) :: data

    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSChSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSChSca", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, separator, csv, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, separator, csv, num, iostat)
    endif


  end subroutine extractDataAttNSChSca


subroutine extractDataAttNSCmplxDpArr(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    complex(dp), dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSCmplxDpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSCmplxDpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSCmplxDpArr

subroutine extractDataAttNSCmplxSpArr(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    complex(sp), dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSCmplxSpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSCmplxSpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSCmplxSpArr

subroutine extractDataAttNSRealDpArr(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    real(dp), dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSRealDpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSRealDpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSRealDpArr

subroutine extractDataAttNSRealSpArr(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    real(sp), dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSRealSpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSRealSpArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSRealSpArr

subroutine extractDataAttNSIntArr(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    integer, dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSIntArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSIntArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSIntArr

subroutine extractDataAttNSLgArr(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    logical, dimension(:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSLgArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSLgArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSLgArr

subroutine extractDataAttNSChArr(arg, namespaceURI, localName, data, separator, csv, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    character(len=*), dimension(:), intent(out) :: data

    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSChArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSChArr", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, separator, csv, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, separator, csv, num, iostat)
    endif


  end subroutine extractDataAttNSChArr


subroutine extractDataAttNSCmplxDpMat(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    complex(dp), dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSCmplxDpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSCmplxDpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSCmplxDpMat

subroutine extractDataAttNSCmplxSpMat(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    complex(sp), dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSCmplxSpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSCmplxSpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSCmplxSpMat

subroutine extractDataAttNSRealDpMat(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    real(dp), dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSRealDpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSRealDpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSRealDpMat

subroutine extractDataAttNSRealSpMat(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    real(sp), dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSRealSpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSRealSpMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSRealSpMat

subroutine extractDataAttNSIntMat(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    integer, dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSIntMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSIntMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSIntMat

subroutine extractDataAttNSLgMat(arg, namespaceURI, localName, data, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    logical, dimension(:,:), intent(out) :: data

    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSLgMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSLgMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif


    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif


  end subroutine extractDataAttNSLgMat

subroutine extractDataAttNSChMat(arg, namespaceURI, localName, data, separator, csv, num, iostat, ex)
    type(DOMException), intent(out), optional :: ex
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    character(len=*), dimension(:,:), intent(out) :: data

    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator
    integer, intent(out), optional :: num, iostat
    if (.not.associated(arg)) then
      if (getFoX_checks().or.FoX_NODE_IS_NULL<200) then
  call throw_exception(FoX_NODE_IS_NULL, "extractDataAttNSChMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    elseif (getNodeType(arg)/=ELEMENT_NODE) then
      if (getFoX_checks().or.FoX_INVALID_NODE<200) then
  call throw_exception(FoX_INVALID_NODE, "extractDataAttNSChMat", ex)
  if (present(ex)) then
    if (inException(ex)) then
        data = ""
       return
    endif
  endif
endif

    endif

    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, separator, csv, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, separator, csv, num, iostat)
    endif


  end subroutine extractDataAttNSChMat



end module m_dom_extras
