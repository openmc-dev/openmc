program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  call cmlStartBand(xf)
  call cmlAddEigenValueVector(xf, eigval=5.43, &
    eigvec=reshape((/(1.0,0.0), (2.0,0.0), (3.0,0.0)/), (/1,3/)), units="units:eV")
  call cmlEndBand(xf)
  call cmlFinishFile(xf)

end program test
