program parse_nullURI

  ! Check that URIs can be parsed and expressed, even if all 
  ! elements are not defined. This is intended to catch 
  ! cases where the parser returns undefined variables.
  use fox_m_utils_uri

  type(URI), pointer :: u

  print*, "Testing empty URI"
  u => parseURI("")
  call check

  print*, "All parts defined"
  u => parseURI("http://uname:pass@www.example.org:80/dir1/dir2/test.html?a=b;c=d#frag")
  call check

  print*, "No #frag"
  u => parseURI("http://uname:pass@www.example.org:80/dir1/dir2/test.html?a=b;c=d")
  call check

  print*, "No ?query"
  u => parseURI("http://uname:pass@www.example.org:80/dir1/dir2/test.html#frag")
  call check

  print*, "No /path"
  u => parseURI("http://uname:pass@www.example.org:80?a=b;c=d#frag")
  call check

  print*, "No :port"
  u => parseURI("http://uname:pass@www.example.org/dir1/dir2/test.html?a=b;c=d#frag")
  call check

  print*, "No user@"
  u => parseURI("http://www.example.org:80/dir1/dir2/test.html?a=b;c=d#frag")
  call check

  print*, "No pass"
  u => parseURI("http://uname@www.example.org:80/dir1/dir2/test.html?a=b;c=d#frag")
  call check

  print*, "No scheme"
  u => parseURI("uname@www.example.org:80/dir1/dir2/test.html?a=b;c=d#frag")
  call check

  contains
    subroutine check
      character(len=100) :: us
      if (associated(u))  then
        call dumpURI(u)
        us = expressURI(u)
        print*, us
        print*
        call destroyURI(u)
      else
        print*, "parsing failed"
      endif
    end subroutine check

end program parse_nullURI
  
