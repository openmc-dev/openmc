program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreatePointStyle(myfile, id='myid', iconhref='http://www.e x a m p l e.com/')
  call kmlFinishFile(myfile)

end program test
