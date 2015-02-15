program test_xml_Openfile

  use FoX_wxml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'

  type(xmlf_t) :: xf1

  call xml_OpenFile(filename, xf1, unit=20)

  !TOHW FIXME use check 20 is open and attached to filename

end program test_xml_Openfile
