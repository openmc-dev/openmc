program test_cmlBeginFile

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'

  type(xmlf_t) :: xf1

  call cmlBeginFile(xf1, filename, unit=-1)
  call cmlFinishFile(xf1)

end program test_cmlBeginFile

