program flush_nag
  use f90_unix_io, only: flush
  print*
  call flush(5)
end program flush_nag
   
