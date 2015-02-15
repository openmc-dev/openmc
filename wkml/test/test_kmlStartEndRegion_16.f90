program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlStartRegion(myfile, reshape((/10.0D0, 10.0D0, 11.0D0, 11.0D0, 12.0D0, 12.0D0/),(/2,3/)), &
           altitude=(/350.0D0, 200.0D0/))

  call kmlEndRegion(myfile)
  call kmlFinishFile(myfile)

end program test
