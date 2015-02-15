program test

  use FoX_wxml, only : xmlf_t
  use FoX_wxml, only : xml_DeclareNamespace
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_DeclareNamespace(xf, "http://www.xml-cml.org/schema")

end program test
