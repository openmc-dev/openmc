module xml_interface

  use fox_dom

  implicit none
  private
  public :: open_xmldoc
  public :: close_xmldoc

contains

!===============================================================================
! OPEN_XMLDOC
!===============================================================================

  subroutine open_xmldoc(ptr, filename)

    character(len=*) :: filename
    type(Node), pointer :: ptr

    ptr => parseFile(trim(filename))

  end subroutine open_xmldoc

!===============================================================================
! CLOSE_XMLDOC
!===============================================================================

  subroutine close_xmldoc(ptr)

    type(Node), pointer :: ptr

    call destroy(ptr)

  end subroutine close_xmldoc

!===============================================================================
! GET_ELEMENT
!===============================================================================

  subroutine get_element(inptr, elename, outptr)

    character(len=*) :: elename
    type(Node), pointer :: inptr
    type(Node), pointer :: outptr
    type(NodeList), pointer :: tempptr

    tempptr => getElementsByTagname(inptr, elename)
    outptr => item(tempptr,0)

  end subroutine get_element

end module xml_interface
