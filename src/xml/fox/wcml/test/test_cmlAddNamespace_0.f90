program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call cmlAddNamespace(xf, 'xhtml', 'http://www.w3.org/1999/xhtml')

end program test
