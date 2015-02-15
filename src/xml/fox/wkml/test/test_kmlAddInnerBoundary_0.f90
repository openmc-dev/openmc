program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlAddInnerBoundary(myfile, (/10.25, 10.25, 10.75, 10.75/), &
                                   (/10.25, 10.75, 10.75, 10.25/))
  call kmlEndRegion(myfile)
  call kmlFinishFile(myfile)

end program test
