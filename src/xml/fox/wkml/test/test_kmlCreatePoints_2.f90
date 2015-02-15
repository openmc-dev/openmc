program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreatePoints(myfile, 1.0, 1.0, 500.0)
  call kmlFinishFile(myfile)

end program test
