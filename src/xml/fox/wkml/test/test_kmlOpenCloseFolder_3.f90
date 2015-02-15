program test

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  call kmlBeginFile(myfile, "test.xml", -1)
  call kmlOpenFolder(myfile, id='folderID')
  call kmlCloseFolder(myfile)
  call kmlFinishFile(myfile)

end program test
