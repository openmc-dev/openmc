program test_fsys_array_str
  ! Some versions of the Intel Fortran Compiler have 
  ! trouble with the creation of long strings using 
  ! fox_m_fsys_array_str - see:
  ! http://software.intel.com/en-us/forums/showthread.php?t=77847&o=a&s=lr
  ! We cannot do much about this compiler bug, but we 
  ! can test for it.
  use fox_m_fsys_array_str, only: str_vs, vs_str_alloc  
  character, pointer :: tempString1(:) => null()  
  character, pointer :: tempString2(:) => null()  
  
  tempString2 => vs_str_alloc( "" )  
  
  ! No stack overflow
  do i = 1, 998  
     tempString1 => tempString2  
     tempString2 => vs_str_alloc( str_vs( tempString1 ) // 'x' )  
     deallocate( tempString1 )  
  end do  
  
  write(*,*) size(tempString2), "'", tempString2, "'"  
  
  deallocate( tempString2 )  
  tempString2 => vs_str_alloc( "" )  
  
  ! Fails - stack overflow.
  do i = 1, 999  
     tempString1 => tempString2  
     tempString2 => vs_str_alloc( str_vs( tempString1 ) // 'x' )  
     deallocate( tempString1 )  
  end do  

  write(*,*) size(tempString2), "'", tempString2, "'"  

  ! Should fail.
  do i = 1, 9999
     tempString1 => tempString2  
     tempString2 => vs_str_alloc( str_vs( tempString1 ) // 'x' )  
     deallocate( tempString1 )  
  end do  
  
  write(*,*) size(tempString2), "'", tempString2, "'"  
  
end program  
