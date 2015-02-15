program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlStartRegion(myfile, (/10.0, 10.0, 11.0, 11.0/), & 
                              (/10.0, 11.0, 11.0, 10.0/))
  call kmlAddInnerBoundary(myfile, reshape((/10.25, 10.25,& 
                 10.25, 10.75, 10.75, 10.75, 10.75, 10.25/),(/2,4/)), &
                           altitude=(/1.1, 1.2, 1.3, 1.4/))

  call kmlEndRegion(myfile)
  call kmlFinishFile(myfile)

end program test
