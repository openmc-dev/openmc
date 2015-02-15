program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreateLine(myfile, (/10.0, 11.0/), (/10.0, 11.0/))
  call kmlFinishFile(myfile)

end program test
