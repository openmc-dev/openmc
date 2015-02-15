program test

  use FoX_wxml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_OpenFile(filename, xf)
  call xml_AddXMLStylesheet(xf, "display.xsl", "text/stylesheet", title="Title")
  call xml_NewElement(xf, "a")
  call xml_Close(xf)

end program test
