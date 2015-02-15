program test

  use FoX_common, only: str
  use FoX_wcml
  use FoX_wxml, only: xml_AddCharacters
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  call cmlAddParameter(xf, name="name", &
    value=reshape((/(1.0d0,2.0d0), (3.0d0,-4.0d0), (-1.0d0,2.5d0),(3.4d0,-10.9d0)/),(/2,2/)), fmt="s3", units="siUnits:m")

  call cmlFinishFile(xf)

end program test
