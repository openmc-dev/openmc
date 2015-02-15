program testUUID

  use fox_m_utils_uuid

  ! This just checks that we can generate the various 
  ! types of UUID (without crashing) and checks that 
  ! they have the correct syntax. We could also check 
  ! that the UUID changes for each call and I think 
  ! there is an additional check we could make within
  ! the UUID itself. But for now this is enough.

  character(len=36) :: uuid

  print*, 'Version 0'
  uuid = generate_uuid(0) 
  if ((check_uuid(uuid)).and.(uuid =='00000000-0000-0000-0000-000000000000')) then
    print*, "OK (is null)"
  else
    print*, "Error: ", uuid
  endif

  print*, 'Version 1'
  uuid = generate_uuid(1) 
  if (check_uuid(uuid)) then
     print*, "OK"
  else 
     print*, "Error: ", uuid
  endif

  print*, 'Version 2'
  uuid = generate_uuid(2) 
  if (uuid == '' ) then
     print*, "OK (not implemented)"
  else 
     print*, "Error: ", uuid
  endif

  print*, 'Version 3'
  uuid = generate_uuid(3) 
  if (uuid == '' ) then
     print*, "OK (not implemented)"
  else 
     print*, "Error: ", uuid
  endif

  print*, 'Version 4'
  uuid = generate_uuid(4) 
  if (check_uuid(uuid)) then
     print*, "OK"
  else 
     print*, "Error: ", uuid
  endif

  print*, 'Version 5'
  uuid = generate_uuid(5) 
  if (uuid == '' ) then
     print*, "OK (not implemented)"
  else 
     print*, "Error: ", uuid
  endif

contains

  function check_uuid(chars) result(lout)

    ! Return true if the string is permitted by the UUID
    ! BFN in RFC

    character(len=*) :: chars
    character(len=22), parameter :: hex = '0123456789abcdefABCDEF'
    logical :: lout

    lout = (len_trim(chars) == 36)
    if (lout) then
        lout = lout.and.(verify(chars(1:8), hex) == 0) 
        lout = lout.and.(verify(chars(9:9), '-') == 0) 
        lout = lout.and.(verify(chars(10:13), hex) == 0) 
        lout = lout.and.(verify(chars(14:14), '-') == 0) 
        lout = lout.and.(verify(chars(15:18), hex) == 0) 
        lout = lout.and.(verify(chars(19:19), '-') == 0) 
        lout = lout.and.(verify(chars(20:23), hex) == 0) 
        lout = lout.and.(verify(chars(24:24), '-') == 0) 
        lout = lout.and.(verify(chars(25:36), hex) == 0) 
    endif

  end function check_uuid

end program testUUID
  
