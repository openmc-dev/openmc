program test

  use FoX_wxml, only : xmlf_t, xml_OpenFile, xml_Close
  use FoX_wxml, only : xml_NewElement, xml_DeclareNamespace, xml_UndeclareNamespace
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_OpenFile(filename, xf)
  call xml_DeclareNamespace(xf, "http://www.xml-cml.org/schema", "cml")
  call xml_NewElement(xf, "cml:cml")
  call xml_DeclareNamespace(xf, "http://www.w3.org/1999/svg", "svg")
  call xml_NewElement(xf, "svg:svg")
  call xml_UndeclareNamespace(xf, "cml")
  call xml_NewElement(xf, "svg:toby")
  call xml_NewElement(xf, "cml:toby")

  call xml_Close(xf)

end program test
