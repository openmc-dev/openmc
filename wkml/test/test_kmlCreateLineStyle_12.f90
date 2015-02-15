program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile
  type(color_t) :: mycolor

  call kmlSetCustomColor(mycolor, 'F90000FF')

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreateLineStyle(myfile, id='myid', width=1, colorhex='f90000ff')
  call kmlFinishFile(myfile)

end program test
