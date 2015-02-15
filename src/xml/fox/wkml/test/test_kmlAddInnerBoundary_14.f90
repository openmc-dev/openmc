program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlStartRegion(myfile, (/10.0, 10.0, 11.0, 11.0/), & 
                              (/10.0, 11.0, 11.0, 10.0/))
  call kmlAddInnerBoundary(myfile, reshape((/10.25D0, 10.25D0,& 
                 10.25D0, 10.75D0, 10.75D0, 10.75D0, 10.75D0, 10.25D0/),(/4,2/)))
  call kmlEndRegion(myfile)
  call kmlFinishFile(myfile)

end program test
