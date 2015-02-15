program testURI

  use fox_m_utils_uri

  type(URI), pointer :: u

  u => parseURI("file:///C:/a%20ctivity.xml")
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

end program testURI
  
