program test

  use FoX_wxml, only : xmlf_t, xml_OpenFile, xml_Close
  use FoX_wxml, only : xml_NewElement, xml_AddCharacters
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_OpenFile(filename, xf)
  call xml_NewElement(xf, "root")
  call xml_AddCharacters(xf, (reshape(    &
     (/(10.0D0, -1.2D1),(10.0D0,-1.2D1),  &
       (10.0D0, -1.2D1),(10.0D0,-1.2D1)/),&
     (/2,2/))))
  call xml_Close(xf)

end program test
