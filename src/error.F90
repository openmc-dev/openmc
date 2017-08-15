module error

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use constants
  use message_passing

  implicit none

  private
  public :: fatal_error
  public :: warning

  ! Error codes
  integer(C_INT), public, bind(C) :: E_UNASSIGNED = -1
  integer(C_INT), public, bind(C) :: E_OUT_OF_BOUNDS = -2
  integer(C_INT), public, bind(C) :: E_CELL_NOT_ALLOCATED = -3
  integer(C_INT), public, bind(C) :: E_CELL_INVALID_ID = -4
  integer(C_INT), public, bind(C) :: E_CELL_NOT_FOUND = -5
  integer(C_INT), public, bind(C) :: E_NUCLIDE_NOT_ALLOCATED = -6
  integer(C_INT), public, bind(C) :: E_NUCLIDE_NOT_LOADED = -7
  integer(C_INT), public, bind(C) :: E_NUCLIDE_NOT_IN_LIBRARY = -8
  integer(C_INT), public, bind(C) :: E_MATERIAL_NOT_ALLOCATED = -9
  integer(C_INT), public, bind(C) :: E_MATERIAL_INVALID_ID = -10
  integer(C_INT), public, bind(C) :: E_TALLY_NOT_ALLOCATED = -11
  integer(C_INT), public, bind(C) :: E_TALLY_INVALID_ID = -12
  integer(C_INT), public, bind(C) :: E_INVALID_SIZE = -13
  integer(C_INT), public, bind(C) :: E_CELL_NO_MATERIAL = -14
  integer(C_INT), public, bind(C) :: E_ALREADY_ALLOCATED = -15
  integer(C_INT), public, bind(C) :: E_ARGUMENT_INVALID = -16
  integer(C_INT), public, bind(C) :: E_WRONG_TYPE = -17
  integer(C_INT), public, bind(C) :: E_FILTER_NOT_ALLOCATED = -18
  integer(C_INT), public, bind(C) :: E_FILTER_INVALID_ID = -19

  ! Warning codes
  integer(C_INT), public, bind(C) :: W_BELOW_MIN_BOUND = 1
  integer(C_INT), public, bind(C) :: W_ABOVE_MAX_BOUND = 2

contains

!===============================================================================
! WARNING issues a warning to the user in the log file and the standard output
! stream.
!===============================================================================

  subroutine warning(message)

    character(*) :: message

    integer :: i_start   ! starting position
    integer :: i_end     ! ending position
    integer :: line_wrap ! length of line
    integer :: length    ! length of message
    integer :: indent    ! length of indentation

    ! Write warning at beginning
    write(ERROR_UNIT, fmt='(1X,A)', advance='no') 'WARNING: '

    ! Set line wrapping and indentation
    line_wrap = 80
    indent = 10

    ! Determine length of message
    length = len_trim(message)

    i_start = 0
    do
      if (length - i_start < line_wrap - indent + 1) then
        ! Remainder of message will fit on line
        write(ERROR_UNIT, fmt='(A)') message(i_start+1:length)
        exit

      else
        ! Determine last space in current line
        i_end = i_start + index(message(i_start+1:i_start+line_wrap-indent+1), &
             ' ', BACK=.true.)

        if (i_end == i_start) then
          ! This is a special case where there is no space
          i_end = i_start + line_wrap - indent + 1
          write(ERROR_UNIT, fmt='(A/A)', advance='no') &
               message(i_start+1:i_end-1), repeat(' ', indent)
          i_end = i_end - 1
        else
          ! Write up to last space
          write(ERROR_UNIT, fmt='(A/A)', advance='no') &
               message(i_start+1:i_end-1), repeat(' ', indent)
        end if

        ! Advance starting position
        i_start = i_end
        if (i_start > length) exit
      end if
    end do

  end subroutine warning

!===============================================================================
! FATAL_ERROR alerts the user that an error has been encountered and displays a
! message about the particular problem. Errors are considered 'fatal' and hence
! the program is aborted.
!===============================================================================

  subroutine fatal_error(message, error_code)

    character(*) :: message
    integer, optional :: error_code ! error code

    integer :: code      ! error code
    integer :: i_start   ! starting position
    integer :: i_end     ! ending position
    integer :: line_wrap ! length of line
    integer :: length    ! length of message
    integer :: indent    ! length of indentation
#ifdef MPI
    integer :: mpi_err
#endif


    ! set default error code
    if (present(error_code)) then
      code = error_code
    else
      code = -1
    end if

    ! Write error at beginning
    write(ERROR_UNIT, fmt='(1X,A)', advance='no') 'ERROR: '

    ! Set line wrapping and indentation
    line_wrap = 80
    indent = 8

    ! Determine length of message
    length = len_trim(message)

    i_start = 0
    do
      if (length - i_start < line_wrap - indent + 1) then
        ! Remainder of message will fit on line
        write(ERROR_UNIT, fmt='(A)') message(i_start+1:length)
        exit

      else
        ! Determine last space in current line
        i_end = i_start + index(message(i_start+1:i_start+line_wrap-indent+1), &
             ' ', BACK=.true.)

        if (i_end == i_start) then
          ! This is a special case where there is no space
          i_end = i_start + line_wrap - indent + 1
          write(ERROR_UNIT, fmt='(A/A)', advance='no') &
               message(i_start+1:i_end-1), repeat(' ', indent)
          i_end = i_end - 1
        else
          ! Write up to last space
          write(ERROR_UNIT, fmt='(A/A)', advance='no') &
               message(i_start+1:i_end-1), repeat(' ', indent)
        end if

        ! Advance starting position
        i_start = i_end
        if (i_start > length) exit
      end if
    end do

#ifdef MPI
    ! Abort MPI
    call MPI_ABORT(mpi_intracomm, code, mpi_err)
#endif

    ! Abort program
#ifdef NO_F2008
    stop
#else
    error stop
#endif

  end subroutine fatal_error

end module error
