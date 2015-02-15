program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1, .true., 'testdoc')
  call kmlFinishFile(myfile)

end program test
