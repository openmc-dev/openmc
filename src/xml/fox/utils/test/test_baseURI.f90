program test_baseURI

  use fox_m_utils_uri

  type(URI), pointer :: u, base, u3

  base => parseURI("/here/we/are")

  u => parseURI("http://tow:bum@www.example.org:80/dir1/dir2/test.html?a=b;c=d#frag")
  call check
  u => parseURI("http://tow:bum@www.example.org:80/dir1/dir2/test.html?a=b;c=d")
  call check
  u => parseURI("http://tow:bum@www.example.org:80/dir1/dir2/test.html#frag")
  call check
  u => parseURI("http://tow:bum@www.example.org:80#frag")
  call check
  u => parseURI("http://tow:bum@www.example.org:80/%7E%20tow/?frag")
  call check
  u => parseURI("file:///Users/tow/devel/FoX/")
  call check
  u => parseURI("file:/Users/tow/devel/FoX/")
  call check

  u => parseURI("//p/Users/tow/devel/FoX/")
  call check

  u => parseURI("/p/Users/tow/devel/FoX/")
  call check

  u => parseURI("Users/tow/devel/FoX/")
  call check

  u => parseURI("./Users/../tow/devel/FoX/")
  call check

  u3 => parseURI("../../../../../tow/devel/FoX/")
  u => rebaseURI(base, u3)
  call check
  call destroyURI(u3)

  call destroyURI(base)

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

end program test_baseURI
  
