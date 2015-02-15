program test

  use FoX_wxml, only : xmlf_t, xml_OpenFile, xml_Close
  use FoX_wxml, only : xml_NewElement, xml_DeclareNamespace
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_OpenFile(filename, xf)
  call xml_DeclareNamespace(xf, "http://www.xml-cml.org/schema", "cml")
  call xml_NewElement(xf, "cml:cml")
  call xml_Close(xf)

end program test
