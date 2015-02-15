program test_xml_Openfile

  use FoX_wxml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'

  type(xmlf_t) :: xf1

  call xml_OpenFile(filename, xf1, replace=.false.)

end program test_xml_Openfile

