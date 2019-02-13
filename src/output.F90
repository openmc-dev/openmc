module output

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use constants
  use eigenvalue,      only: openmc_get_keff
  use endf,            only: reaction_name
  use error,           only: fatal_error, warning
  use geometry_header
  use math,            only: t_percentile
  use message_passing, only: master, n_procs
  use mgxs_interface
  use nuclide_header
  use particle_header, only: LocalCoord, Particle
  use settings
  use simulation_header
  use surface_header,  only: surfaces
  use string,          only: to_upper, to_str
  use tally_header
  use tally_derivative_header
  use tally_filter

  implicit none

  ! Short names for output and error units
  integer :: ou = OUTPUT_UNIT
  integer :: eu = ERROR_UNIT

  interface
    subroutine print_particle(p) bind(C)
      import Particle
      type(Particle), intent(in) :: p
    end subroutine

    subroutine write_tallies() bind(C)
    end subroutine
  end interface

contains

!===============================================================================
! TIME_STAMP returns the current date and time in a formatted string
!===============================================================================

  function time_stamp() result(current_time)

    character(19) :: current_time ! ccyy-mm-dd hh:mm:ss
    character(8)  :: date_        ! ccyymmdd
    character(10) :: time_        ! hhmmss.sss

    call date_and_time(DATE=date_, TIME=time_)
    current_time = date_(1:4) // "-" // date_(5:6) // "-" // date_(7:8) // &
         " " // time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)

  end function time_stamp

!===============================================================================
! HEADER displays a header block according to a specified level. If no level is
! specified, it is assumed to be a minor header block.
!===============================================================================

  subroutine header(msg, level, unit)
    character(*), intent(in)      :: msg   ! header message
    integer, intent(in)           :: level
    integer, intent(in), optional :: unit  ! unit to write to

    integer :: n            ! number of = signs on left
    integer :: m            ! number of = signs on right
    integer :: unit_        ! unit to write to
    character(MAX_LINE_LEN) :: line

    ! set default unit
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! determine how many times to repeat '=' character
    n = (63 - len_trim(msg))/2
    m = n
    if (mod(len_trim(msg),2) == 0) m = m + 1

    ! convert line to upper case
    line = to_upper(msg)

    ! print header based on level
    if (verbosity >= level) then
      write(UNIT=unit_, FMT='(/1X,A/)') repeat('=', n) // '>     ' // &
           trim(line) // '     <' // repeat('=', m)
    end if

  end subroutine header

end module output
