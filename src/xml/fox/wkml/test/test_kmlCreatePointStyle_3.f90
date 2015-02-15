program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreatePointStyle(myfile, id='myid', scale=10.0, colorname='red')
  call kmlFinishFile(myfile)

end program test
