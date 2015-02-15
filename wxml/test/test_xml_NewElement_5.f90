program test

  use FoX_wxml, only : xmlf_t, xml_OpenFile, xml_Close
  use FoX_wxml, only : xml_NewElement, xml_AddDOCTYPE
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_OpenFile(filename, xf, validate=.true.)
  call xml_AddDOCTYPE(xf, "html", "http://www.w3.org/1999/xhtml")
  call xml_NewElement(xf, "root")
  call xml_Close(xf)

end program test
