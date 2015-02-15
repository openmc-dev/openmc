program test_sax_reader

  use FoX_sax
  type(xml_t) :: xp
  integer :: iostat

  ! Can we parse the same file twice?

  call open_xml_file(xp, "testin.xml", iostat)

  write(*,'(i1)') iostat

  call parse(xp)

  call parse(xp)

  call close_xml_t(xp) 
end program test_sax_reader
