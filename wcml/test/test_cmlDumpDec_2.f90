program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf, compchem=.true.)
  call wcmlDumpDec(xf, '_no_file_to_find.input', 80, .false.)

  call cmlFinishFile(xf)

end program test


