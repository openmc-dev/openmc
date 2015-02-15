program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreateLine(myfile, (/10.0D0, 11.0D0, 12.0D0, 13.0D0/), (/10.0D0, 11.0D0, 12.0D0, 13.0D0/))
  call kmlFinishFile(myfile)

end program test
