program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlStartRegion(myfile, reshape((/10.0D0, 10.0D0, 250.0D0, 11.0D0, 11.0D0, 350.0D0, 12.0D0, 12.0D0, 300.0D0/),(/3,3/)))
  call kmlEndRegion(myfile)
  call kmlFinishFile(myfile)

end program test
