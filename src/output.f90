module output

  use ISO_FORTRAN_ENV
  use global

  implicit none

  contains

!=====================================================================
! TITLE prints the main title banner as well as information about the
! program developers, version, and date/time which the problem was
! run.
!=====================================================================

    subroutine title()

      character(10) :: date
      character(8)  :: time
      integer       :: ou

      ou = OUTPUT_UNIT

      write(ou,*)
      write(ou,*) ' .d88888b.                             888b     d888  .d8888b. '
      write(ou,*) 'd88P" "Y88b                            8888b   d8888 d88P  Y88b'
      write(ou,*) '888     888                            88888b.d88888 888    888'
      write(ou,*) '888     888 88888b.   .d88b.  88888b.  888Y88888P888 888       '
      write(ou,*) '888     888 888 "88b d8P  Y8b 888 "88b 888 Y888P 888 888       '
      write(ou,*) '888     888 888  888 88888888 888  888 888  Y8P  888 888    888'
      write(ou,*) 'Y88b. .d88P 888 d88P Y8b.     888  888 888   "   888 Y88b  d88P'
      write(ou,*) ' "Y88888P"  88888P"   "Y8888  888  888 888       888  "Y8888P" '
      write(ou,*) '____________888________________________________________________'
      write(ou,*) '            888                                                '
      write(ou,*) '            888                                                '
      write(ou,*)

      ! Write version information
!      write(ou,*) '     ______________________________________________________'
      write(ou,*) '     Developed At:    Massachusetts Institute of Technology'
      write(ou,*) '     Lead Developer:  Paul K. Romano'
      write(ou,100) VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
100   format (6X,"Version:",9X,I1,".",I1,".",I1)

      ! Write the date and time
      call get_today( date, time )
      write(ou,101) trim(date), trim(time)
101   format (6X,"Date/Time:",7X,A,1X,A)
!      write(ou,*) '     ______________________________________________________'
      write(ou,*)

    end subroutine title

!=====================================================================
! ECHO_INPUT displays summary information about the problem about to
! be run after reading all the input.
!=====================================================================

    subroutine echo_input()

      character(32) :: string
      integer :: ou

      ou = OUTPUT_UNIT

      ! Display problem summary
      write(ou,*) '============================================='
      write(ou,*) '=>             PROBLEM SUMMARY             <='
      write(ou,*) '============================================='
      write(ou,*)
      if ( problem_type == PROB_CRITICALITY ) then
         write(ou,100) 'Problem type:', 'Criticality'
         write(ou,100) 'Number of Cycles:', int_to_str(n_cycles)
         write(ou,100) 'Number of Inactive Cycles:', int_to_str(n_inactive) 
      elseif ( problem_type == PROB_SOURCE ) then
         write(ou,100) 'Problem type:', 'External Source'
      end if
      write(ou,100) 'Number of Particles:', int_to_str(n_particles)
      write(ou,*)      
      ! Display geometry summary
      write(ou,*) '============================================='
      write(ou,*) '=>              GEOMETRY SUMMARY           <='
      write(ou,*) '============================================='
      write(ou,*)
      write(ou,100) 'Number of Cells:', int_to_str(n_cells)
      write(ou,100) 'Number of Surfaces:', int_to_str(n_surfaces)
      write(ou,100) 'Number of Materials:', int_to_str(n_materials)
      write(ou,*)

      ! Format descriptor for columns
100   format (1X,A,T35,A)

    end subroutine echo_input

!=====================================================================
! MESSAGE displays an informational message to the log file and the
! standard output stream.
!=====================================================================

    subroutine message( msg, level )

      character(*), intent(in) :: msg
      integer, intent(in) :: level

      if ( level <= verbosity ) then
         write (OUTPUT_UNIT,*) trim(msg)
      end if

    end subroutine message

!=====================================================================
! WARNING issues a warning to the user in the log file and the
! standard output stream.
!=====================================================================

    subroutine warning( msg )

      character(*), intent(in) :: msg

      integer :: n_lines
      integer :: i, ou

      ou = OUTPUT_UNIT

      write(ou, fmt='(1X,A9)', advance='no') 'WARNING: '

      n_lines = (len_trim(msg)-1)/70 + 1
      do i = 1, n_lines
         if ( i == 1 ) then
            write(ou, fmt='(A70)') msg(70*(i-1)+1:70*i)
         else
            write(ou, fmt='(10X,A70)') msg(70*(i-1)+1:70*i)
         end if
      end do
      
    end subroutine warning

!=====================================================================
! ERROR alerts the user that an error has been encountered and
! displays a message about the particular problem. Errors are
! considered 'fatal' and hence the program is aborted.
!=====================================================================

    subroutine error( msg )

      character(*), intent(in) :: msg

      integer :: n_lines
      integer :: i, eu

      eu = ERROR_UNIT

      write(eu, fmt='(1X,A7)', advance='no') 'ERROR: '

      n_lines = (len_trim(msg)-1)/72 + 1
      do i = 1, n_lines
         if ( i == 1 ) then
            write(eu, fmt='(A72)') msg(72*(i-1)+1:72*i)
         else
            write(eu, fmt='(7X,A72)') msg(72*(i-1)+1:72*i)
         end if
      end do
      write(eu,*)

      call free_memory()
      
    end subroutine error

!=====================================================================
! GET_TODAY determines the date and time at which the program began
! execution and returns it in a readable format
!=====================================================================

    subroutine get_today( today_date, today_time )

      character(10), intent(out) :: today_date
      character(8),  intent(out) :: today_time

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

end module output
