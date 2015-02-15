program test_xml_Openfile

  use FoX_wxml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'

  type(xmlf_t) :: xf1

  call xml_OpenFile(filename, xf1)
  call xml_Close(xf1, empty=.true.)

end program test_xml_Openfile

