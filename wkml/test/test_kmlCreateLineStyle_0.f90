program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlCreateLineStyle(myfile, id='myid')
  call kmlFinishFile(myfile)

end program test
