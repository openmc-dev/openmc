program test_sax_reader

  use m_sax_reader
  use m_common_error,  only: error_stack

  type(file_buffer_t) :: fb

  type(error_stack) :: es
  integer :: iostat
  logical :: eof
  character(len=1) :: c

  call open_file(fb, file="test_sax_reader_1.in",  iostat=iostat, es=es)

  c = get_character(fb, eof, es)
  write(*,'(2a)') 'char:', c
  write(*,'(a,i0)') 'iost:', iostat

  call close_file(fb) 

end program test_sax_reader
