program test

  use FoX_wxml, only : xmlf_t, xml_OpenFile, xml_Close
  use FoX_wxml, only : xml_NewElement, xml_AddAttribute
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_OpenFile(filename, xf)
  call xml_NewElement(xf, "root")
  call xml_AddAttribute(xf, "att1", "value")
  call xml_AddAttribute(xf, "att2", "value")
  call xml_AddAttribute(xf, "att3", "value")
  call xml_AddAttribute(xf, "att4", "value")
  call xml_AddAttribute(xf, "att5", "value")
  call xml_AddAttribute(xf, "att6", "value")
  call xml_AddAttribute(xf, "att7", "value")
  call xml_AddAttribute(xf, "att8", "value")
  call xml_AddAttribute(xf, "att9", "value")
  call xml_AddAttribute(xf, "att10", "value")
  call xml_AddAttribute(xf, "att11", "value")
  call xml_AddAttribute(xf, "att12", "value")
  call xml_AddAttribute(xf, "att13", "value")
  call xml_Close(xf)

end program test
