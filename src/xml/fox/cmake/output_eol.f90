program output_eol
  
  character(len=1)  :: single_char = "a"
  integer           :: ios
  integer(kind=1)   :: n1,n2,n3
  
  n1 = 0
  n2 = 0
  n3 = 0
  ! open a file, write a single character and close it again
  open(unit=10,file="eol.out",status="unknown")
  write(10,"(a)") single_char
  close(10)
  
  ! open the file unformatted, stream
  open(unit=11,file="eol.out",status="old",form="unformatted",access="stream",iostat=ios)
  
  ! read first character (should be an a)
  read(11,iostat=ios) n1 ! should be 97 (ascii a)
  if ( ios < 0 ) then
    !print *, "end of file reached after 0 characters"
  else
    ! read second character, should be lf (ascii 10) on dos and unix,
    ! and cr on mac (ascii 13)
    read(11,iostat=ios) n2
    if ( ios < 0 ) then
      ! this should not happen now, there are at least 2 characters in the file
      !print *, "end of file reached after 1 character"
    else
      ! read third character if we're on dos
      read(11,iostat=ios) n3
      if ( ios < 0 ) then
        !print *, "end of file reached after 2 characters"
      end if
    end if
  end if
  
  ! analyse n2 and n3
  if ( n2 == 10 .and. n3 == 0 ) then
    print *, "LF"
  elseif ( n2 == 13 .and. n3 == 10 ) then
    print *, "CRLF"
  elseif ( n2 == 13 .and. n3 == 0 ) then
    print *, "CR"
  else
    print *, "could not determine line-ending convention"
  end if
  close(11)
  
  
end program

