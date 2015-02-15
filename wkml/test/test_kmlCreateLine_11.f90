program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreateLine(myfile, reshape((/10.0, 10.0, 1.0, 11.0, 11.0, 2.0, 12.0, 12.0, 3.0/),(/3,3/)), &
                       altitude=(/1.0, 2.0, 3.0/))
  call kmlFinishFile(myfile)

end program test
