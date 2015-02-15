program test

  use FoX_wxml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_AddXMLStylesheet(xf, "display.xsl", "text/stylesheet")

end program test
