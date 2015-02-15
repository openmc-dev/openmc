program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreateLineStyle(myfile, id='myid', width=100)
  call kmlFinishFile(myfile)

end program test
