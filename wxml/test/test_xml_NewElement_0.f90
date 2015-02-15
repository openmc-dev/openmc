program test

  use FoX_wxml, only : xmlf_t
  use FoX_wxml, only : xml_NewElement
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_NewElement(xf, "root")

end program test
