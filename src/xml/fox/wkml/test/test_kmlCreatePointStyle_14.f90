program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreatePointStyle(myfile, id='myid', colorhex='F90000FF', colormode='modnar')
  call kmlFinishFile(myfile)

end program test
