program test

  use FoX_wcml
  use FoX_common, only: str
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  call cmlAddProperty(xf, title="name", &
    value=reshape((/"value1", "value2", "value3", "value4"/), (/2,2/)))

  call cmlFinishFile(xf)

end program test
