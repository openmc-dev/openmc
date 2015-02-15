program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  call cmlStartBand(xf)
  call cmlAddEigenValue(xf, 5.43d0, units="units:eV")
  call cmlEndBand(xf)
  call cmlFinishFile(xf)

end program test
