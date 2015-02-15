program test

  use FoX_wkml

  implicit none

  type(color_t) :: colourmap(1)
  type(xmlf_t) :: myfile

  colourmap(1) = kmlGetCustomColor(193)


  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreatePointStyle(myfile, id='myid', color=colourmap(1))
  call kmlFinishFile(myfile)

end program test
