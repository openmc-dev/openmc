program test

  use FoX_wxml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_AddXMLDeclaration(xf)
  call xml_NewElement(xf, "a")

end program test
