program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlStartRegion(myfile, reshape((/10.0, 10.0, 250.0, 11.0, 11.0, 350.0, 12.0, 12.0, 300.0/),(/3,3/)))
  call kmlEndRegion(myfile)
  call kmlFinishFile(myfile)

end program test
