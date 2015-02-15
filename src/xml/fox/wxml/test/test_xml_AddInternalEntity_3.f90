program test

  use FoX_wxml, only : xmlf_t, xml_OpenFile, xml_Close
  use FoX_wxml, only : xml_NewElement, xml_AddDOCTYPE, xml_AddInternalEntity
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_OpenFile(filename, xf)
  call xml_AddDOCTYPE(xf, 'html')
  call xml_NewElement(xf, 'html')
  call xml_AddInternalEntity(xf, 'copy', 'Copyright Toby White 2006')
  call xml_Close(xf)

end program test
