program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreatePoints(myfile, (/1.0, 2.0/), (/1.0, 2.0/), (/500.0, 1000.0/), extrude=.true.)
  call kmlFinishFile(myfile)

end program test
