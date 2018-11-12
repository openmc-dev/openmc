module initialize

  use, intrinsic :: ISO_C_BINDING

  use settings
  use string, only: to_f_string

  implicit none

  interface
    function openmc_path_input() result(ptr) bind(C)
      import C_PTR
      type(C_PTR) :: ptr
    end function
    function openmc_path_output() result(ptr) bind(C)
      import C_PTR
      type(C_PTR) :: ptr
    end function
    function openmc_path_particle_restart() result(ptr) bind(C)
      import C_PTR
      type(C_PTR) :: ptr
    end function
    function openmc_path_statepoint() result(ptr) bind(C)
      import C_PTR
      type(C_PTR) :: ptr
    end function
    function openmc_path_sourcepoint() result(ptr) bind(C)
      import C_PTR
      type(C_PTR) :: ptr
    end function
  end interface

contains

!===============================================================================
! READ_COMMAND_LINE reads all parameters from the command line
!===============================================================================

  subroutine read_command_line() bind(C)
    ! Arguments were already read on C++ side (initialize.cpp). Here we just
    ! convert the C-style strings to Fortran style

    character(kind=C_CHAR), pointer :: string(:)
    interface
      function is_null(ptr) result(x) bind(C)
        import C_PTR, C_BOOL
        type(C_PTR), value :: ptr
        logical(C_BOOL) :: x
      end function is_null
    end interface

    if (.not. is_null(openmc_path_input())) then
      call c_f_pointer(openmc_path_input(), string, [255])
      path_input = to_f_string(string)
    else
      path_input = ''
    end if
    if (.not. is_null(openmc_path_statepoint())) then
      call c_f_pointer(openmc_path_statepoint(), string, [255])
      path_state_point = to_f_string(string)
    end if
    if (.not. is_null(openmc_path_sourcepoint())) then
      call c_f_pointer(openmc_path_sourcepoint(), string, [255])
      path_source_point = to_f_string(string)
    end if
    if (.not. is_null(openmc_path_particle_restart())) then
      call c_f_pointer(openmc_path_particle_restart(), string, [255])
      path_particle_restart = to_f_string(string)
    end if
  end subroutine read_command_line

end module initialize
