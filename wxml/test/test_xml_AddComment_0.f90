program test

  use FoX_wxml, only : xmlf_t
  use FoX_wxml, only : xml_AddComment
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_AddComment(xf, "comment")

end program test
