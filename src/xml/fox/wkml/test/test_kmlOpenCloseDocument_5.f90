program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1, .true., 'testdoc')
  call kmlOpenDocument(myfile, "NewName")
  call kmlCloseDocument(myfile)
  call kmlFinishFile(myfile)

end program test
