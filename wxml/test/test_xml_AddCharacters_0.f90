program test

  use FoX_wxml, only : xmlf_t
  use FoX_wxml, only : xml_AddCharacters
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_AddCharacters(xf, "characters")

end program test
