program test

  use FoX_wkml

  implicit none

  type(color_t) :: colourmap(1)
  type(xmlf_t) :: myfile

  call kmlSetCustomColor(colourmap(1), 'j90000ff')


  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlCreatePointStyle(myfile, id='myid', color=colourmap(1))
  call kmlFinishFile(myfile)

end program test
