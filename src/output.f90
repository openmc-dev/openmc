module output

  use ISO_FORTRAN_ENV
  use global,  only: PATHLENGTH, COLLISION, ABSORPTION, FISSION, &
       &             NUFISSION, nhist, ncycle, localmesh, ijk_to_n, &
       &             nprocs, verbosity, VERSION_MAJOR, VERSION_MINOR, &
       &             VERSION_RELEASE, free_memory

  implicit none

  contains

!-----------------------------------------------------------------------

    subroutine title()

      character(10) :: date
      character(8)  :: time

      write(OUTPUT_UNIT,*)
      write(OUTPUT_UNIT,*) ' .d88888b.                             888b     d888  .d8888b. '
      write(OUTPUT_UNIT,*) 'd88P" "Y88b                            8888b   d8888 d88P  Y88b'
      write(OUTPUT_UNIT,*) '888     888                            88888b.d88888 888    888'
      write(OUTPUT_UNIT,*) '888     888 88888b.   .d88b.  88888b.  888Y88888P888 888       '
      write(OUTPUT_UNIT,*) '888     888 888 "88b d8P  Y8b 888 "88b 888 Y888P 888 888       '
      write(OUTPUT_UNIT,*) '888     888 888  888 88888888 888  888 888  Y8P  888 888    888'
      write(OUTPUT_UNIT,*) 'Y88b. .d88P 888 d88P Y8b.     888  888 888   "   888 Y88b  d88P'
      write(OUTPUT_UNIT,*) ' "Y88888P"  88888P"   "Y8888  888  888 888       888  "Y8888P" '
      write(OUTPUT_UNIT,*) '____________888________________________________________________'
      write(OUTPUT_UNIT,*) '            888                                                '
      write(OUTPUT_UNIT,*) '            888                                                '
      write(OUTPUT_UNIT,*)

      ! Write version information
!      write(OUTPUT_UNIT,*) '     ______________________________________________________'
      write(OUTPUT_UNIT,*) '     Developed At:    Massachusetts Institute of Technology'
      write(OUTPUT_UNIT,*) '     Lead Developer:  Paul K. Romano'
      write(OUTPUT_UNIT,100) VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
100   format (6X,"Version:",9X,I1,".",I1,".",I1)

      ! Write the date and time
      call get_today( date, time )
      write(OUTPUT_UNIT,101) trim(date), trim(time)
101   format (6X,"Date/Time:",7X,A,1X,A)
!      write(OUTPUT_UNIT,*) '     ______________________________________________________'
      write(OUTPUT_UNIT,*)

    end subroutine title

!-----------------------------------------------------------------------

    subroutine message( msg, level )

      character(*), intent(in) :: msg
      integer, intent(in) :: level

      if ( level <= verbosity ) then
         write (OUTPUT_UNIT,*) trim(msg)
      end if

    end subroutine message

!-----------------------------------------------------------------------

    subroutine warning( msg )

      character(*), intent(in) :: msg

      integer :: n_lines
      integer :: i

      write(OUTPUT_UNIT, fmt='(1X,A9)', advance='no') 'WARNING: '

      n_lines = (len_trim(msg)-1)/70 + 1
      do i = 1, n_lines
         if ( i == 1 ) then
            write(OUTPUT_UNIT, fmt='(A70)') msg(70*(i-1)+1:70*i)
         else
            write(OUTPUT_UNIT, fmt='(10X,A70)') msg(70*(i-1)+1:70*i)
         end if
      end do
      
    end subroutine warning

!-----------------------------------------------------------------------

    subroutine error( msg )

      character(*), intent(in) :: msg

      integer :: n_lines
      integer :: i

      write(ERROR_UNIT, fmt='(1X,A7)', advance='no') 'ERROR: '

      n_lines = (len_trim(msg)-1)/72 + 1
      do i = 1, n_lines
         if ( i == 1 ) then
            write(ERROR_UNIT, fmt='(A72)') msg(72*(i-1)+1:72*i)
         else
            write(ERROR_UNIT, fmt='(7X,A72)') msg(72*(i-1)+1:72*i)
         end if
      end do
      write(ERROR_UNIT,*)

      call free_memory()
      
    end subroutine error

!-----------------------------------------------------------------------

    subroutine get_today( today_date, today_time )

      character(10) :: today_date
      character(8)  :: today_time

      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer       :: val(8)

      call date_and_time(date, time, zone, val)
      ! val(1) = year (YYYY)
      ! val(2) = month (MM)
      ! val(3) = day (DD)
      ! val(4) = timezone
      ! val(5) = hours (HH)
      ! val(6) = minutes (MM)
      ! val(7) = seconds (SS)
      ! val(8) = milliseconds
      
      if ( val(2) < 10 ) then
         if ( val(3) < 10 ) then
            today_date = date(6:6) // "/" // date(7:7) // "/" // date(1:4)
         else
            today_date = date(6:6) // "/" // date(7:8) // "/" // date(1:4)
         end if
      else
         if ( val(3) < 10 ) then
            today_date = date(5:6) // "/" // date(7:7) // "/" // date(1:4)
         else
            today_date = date(5:6) // "/" // date(7:8) // "/" // date(1:4)
         end if
      end if
      today_time = time(1:2) // ":" // time(3:4) // ":" // time(5:6)

    end subroutine get_today

!-----------------------------------------------------------------------

end module output




