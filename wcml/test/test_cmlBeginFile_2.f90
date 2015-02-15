program test_xml_Openfile

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'

  type(xmlf_t) :: xf1

  integer :: n

  call cmlBeginFile(xf1, filename, unit=20)

  inquire(file=filename, number=n)

  if (n==20) then
    write(*,'(a)') "file is attached"
  else
    write(*,'(a)') "file is not attached"
  endif

end program test_xml_Openfile
