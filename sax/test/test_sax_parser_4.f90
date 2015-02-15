program test_sax_reader

  use FoX_sax
  type(xml_t) :: xp
  integer :: iostat

  ! Can we parse the same file twice?

  call open_xml_file(xp, "test_sax_fsm_1.in", iostat)

  write(*,'(i1)') iostat

  call parse(xp)

  call close_xml_t(xp) 

  call open_xml_file(xp, "test_sax_fsm_1.in", iostat)

  write(*,'(i1)') iostat

  call parse(xp)

  call close_xml_t(xp) 
end program test_sax_reader
