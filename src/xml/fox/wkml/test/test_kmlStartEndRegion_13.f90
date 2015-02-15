program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlStartRegion(myfile, (/10.0, 11.0, 12.0/), (/10.0, 11.0, 12.0/), &
           altitude=(/350.0, 200.0/))
  call kmlEndRegion(myfile)
  call kmlFinishFile(myfile)

end program test
