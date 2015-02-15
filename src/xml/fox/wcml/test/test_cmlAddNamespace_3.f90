program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  call cmlAddNamespace(xf, 'xhtml', 'http://www.w3.org/1999/xhtml')
  call cmlFinishFile(xf)

end program test
