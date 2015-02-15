program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlStartRegion(myfile, (/10.0, 11.0, 12.0/), (/10.0, 11.0, 12.0/))
  call kmlEndRegion(myfile)
  call kmlFinishFile(myfile)

end program test
