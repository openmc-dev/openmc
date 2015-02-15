program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlStartRegion(myfile, (/10.0, 10.0, 11.0, 11.0/), & 
                              (/10.0, 11.0, 11.0, 10.0/))
  call kmlAddInnerBoundary(myfile, reshape((/10.25, 10.25, 1.1,& 
                 10.25, 10.75, 1.2, 10.75, 10.75, 1.3, 10.75, 10.25, 1.4/),(/3,4/)))

  call kmlEndRegion(myfile)
  call kmlFinishFile(myfile)

end program test
