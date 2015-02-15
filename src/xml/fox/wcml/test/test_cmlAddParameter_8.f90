program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  call cmlAddParameter(xf, name="name", value=(/"value1", "value2"/))

  call cmlFinishFile(xf)

end program test
