program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  call cmlAddParameter(xf, name="name", value=(/(1.0d0,2.0d0), (3.0d0,-4.0d0)/), units="siUnits:m", fmt="s3")

  call cmlFinishFile(xf)

end program test
