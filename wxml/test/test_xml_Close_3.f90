program test

  use FoX_wxml, only : xmlf_t, xml_Close
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_Close(xf)

end program test
