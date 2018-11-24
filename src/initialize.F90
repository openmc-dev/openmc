module initialize

  use, intrinsic :: ISO_C_BINDING

  use settings
  use string, only: to_f_string

  implicit none

  interface
    function path_input_c() result(ptr) bind(C)
      import C_PTR
      type(C_PTR) :: ptr
    end function
    function path_output_c() result(ptr) bind(C)
      import C_PTR
      type(C_PTR) :: ptr
    end function
    function path_particle_restart_c() result(ptr) bind(C)
      import C_PTR
      type(C_PTR) :: ptr
    end function
    function path_statepoint_c() result(ptr) bind(C)
      import C_PTR
      type(C_PTR) :: ptr
    end function
    function path_sourcepoint_c() result(ptr) bind(C)
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

    if (.not. is_null(path_input_c())) then
      call c_f_pointer(path_input_c(), string, [255])
      path_input = to_f_string(string)
    else
      path_input = ''
    end if
    if (.not. is_null(path_statepoint_c())) then
      call c_f_pointer(path_statepoint_c(), string, [255])
      path_state_point = to_f_string(string)
    end if
    if (.not. is_null(path_sourcepoint_c())) then
      call c_f_pointer(path_sourcepoint_c(), string, [255])
      path_source_point = to_f_string(string)
    end if
    if (.not. is_null(path_particle_restart_c())) then
      call c_f_pointer(path_particle_restart_c(), string, [255])
      path_particle_restart = to_f_string(string)
    end if
  end subroutine read_command_line

end module initialize
