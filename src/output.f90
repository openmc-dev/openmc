module output

  use ISO_FORTRAN_ENV
  use global
  use types, only: Cell, Universe, Surface
  use data_structures, only: dict_get_key

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

      integer :: ou
      integer :: n_lines
      integer :: i

      if ( level <= verbosity ) then
         ou = OUTPUT_UNIT
         
         n_lines = (len_trim(msg)-1)/79 + 1
         do i = 1, n_lines
            write(ou, fmt='(1X,A79)') msg(79*(i-1)+1:79*i)
         end do

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
            today_date = date(6:6) // "/" // date(8:8) // "/" // date(1:4)
         else
            today_date = date(6:6) // "/" // date(7:8) // "/" // date(1:4)
         end if
      else
         if ( val(3) < 10 ) then
            today_date = date(5:6) // "/" // date(8:8) // "/" // date(1:4)
         else
            today_date = date(5:6) // "/" // date(7:8) // "/" // date(1:4)
         end if
      end if
      today_time = time(1:2) // ":" // time(3:4) // ":" // time(5:6)

    end subroutine get_today

!=====================================================================
! PRINT_CELL displays the attributes of a cell
!=====================================================================

    subroutine print_cell(c)

      type(Cell), pointer :: c

      integer :: ou
      integer :: temp
      integer :: i
      type(Universe), pointer :: u => null()
      type(Material), pointer :: m => null()
      character(250) :: string

      ou = OUTPUT_UNIT

      write(ou,*) 'Cell ' // int_to_str(c % uid)
      temp = dict_get_key(cell_dict, c % uid)
      write(ou,*) '    Array Index = ' // int_to_str(temp)
      u => universes(c % universe)
      write(ou,*) '    Universe = ' // int_to_str(u % uid)
      if (c % fill == 0) then
         write(ou,*) '    Fill = NONE'
      else
         u => universes(c % fill)
         write(ou,*) '    Fill = ' // int_to_str(u % uid)
      end if
      if (c % material == 0) then
         write(ou,*) '    Material = NONE'
      else
         m => materials(c % material)
         write(ou,*) '    Material = ' // int_to_str(m % uid)
      end if
      write(ou,*) '    Parent Cell = ' // int_to_str(c % parent)
      string = ""
      do i = 1, c % n_items
         select case (c % boundary_List(i))
         case (OP_LEFT_PAREN)
            string = trim(string) // ' ('
         case (OP_RIGHT_PAREN)
            string = trim(string) // ' )'
         case (OP_UNION)
            string = trim(string) // ' :'
         case (OP_DIFFERENCE)
            string = trim(string) // ' !'
         case default
            string = trim(string) // ' ' // int_to_str(c % boundary_list(i))
         end select
      end do
      write(ou,*) '    Surface Specification:' // trim(string)
      write(ou,*)

      ! nullify associated pointers
      nullify(u)
      nullify(m)

    end subroutine print_cell

!=====================================================================
! PRINT_UNIVERSE displays the attributes of a universe
!=====================================================================

    subroutine print_universe(univ)

      type(Universe), pointer :: univ

      integer :: ou
      integer :: i
      character(250) :: string
      type(Cell), pointer :: c

      ou = OUTPUT_UNIT

      write(ou,*) 'Universe ' // int_to_str(univ % uid)
      write(ou,*) '    Level = ' // int_to_str(univ % level)
      string = ""
      do i = 1, univ % n_cells
         c => cells(univ % cells(i))
         string = trim(string) // ' ' // int_to_str(c % uid)
      end do
      write(ou,*) '    Cells =' // trim(string)
      write(ou,*)

    end subroutine print_universe

!=====================================================================
! PRINT_SURFACE displays the attributes of a surface
!=====================================================================

    subroutine print_surface(surf)

      type(Surface), pointer :: surf

      integer :: ou
      integer :: i
      character(80) :: string

      ou = OUTPUT_UNIT

      write(ou,*) 'Surface ' // int_to_str(surf % uid)
      select case (surf % type)
      case (SURF_PX)
         string = "X Plane"
      case (SURF_PY)
         string = "Y Plane"
      case (SURF_PZ)
         string = "Z Plane"
      case (SURF_PLANE)
         string = "Plane"
      case (SURF_CYL_X)
         string = "X Cylinder"
      case (SURF_CYL_Y)
         string = "Y Cylinder"
      case (SURF_CYL_Z)
         string = "Z Cylinder"
      case (SURF_SPHERE)
         string = "Sphere"
      case (SURF_BOX_X)
      case (SURF_BOX_Y)
      case (SURF_BOX_Z)
      case (SURF_BOX)
      case (SURF_GQ)
         string = "General Quadratic"
      end select
      write(ou,*) '    Type = ' // trim(string)
      write(ou,*) '    Coefficients = ', surf % coeffs

      string = ""
      if (allocated(surf % neighbor_pos)) then
         do i = 1, size(surf % neighbor_pos)
            string = trim(string) // ' ' // int_to_str(surf % neighbor_pos(i))
         end do
      end if
      write(ou,*) '    Positive Neighbors = ', string

      string = ""
      if (allocated(surf % neighbor_neg)) then
         do i = 1, size(surf % neighbor_neg)
            string = trim(string) // ' ' // int_to_str(surf % neighbor_neg(i))
         end do
      end if
      write(ou,*) '    Negative Neighbors =', string

    end subroutine print_surface

!=====================================================================
! PRINT_MATERIAL displays the attributes of a material
!=====================================================================

    subroutine print_material(mat)

      type(Material), pointer :: mat

      integer :: ou
      integer :: i
      integer :: n_lines
      type(AceContinuous), pointer :: table
      character(250) :: string

      ou = OUTPUT_UNIT

      write(ou,*) 'Material ' // int_to_str(mat % uid)
      ! Make string of all isotopes
      string = ""
      do i = 1, mat % n_isotopes
         table => xs_continuous(mat % table(i))
         string = trim(string) // ' ' // table % name
      end do
      ! Print isotopes with word wrap
      ! TODO: Change this to generic word wrap subroutine?
      n_lines = (len_trim(string)-1)/75 + 1
      do i = 1, n_lines
         if ( i == 1 ) then
            write(ou, fmt='(5X,A75)') 'Isotopes =' // string(75*(i-1)+1:75*i)
         else
            write(ou, fmt='(5X,A75)') string(75*(i-1)+1:75*i)
         end if
      end do
      write(ou,'(5X,A,G12.4,A)') 'Atom Density = ', mat % atom_density, & 
           & ' atom/b-cm'
      write(ou,*)

    end subroutine print_material

!=====================================================================
! PRINT_SUMMARY displays the attributes of all cells, universes,
! surfaces and materials read in the input file. Very useful for
! debugging!
!=====================================================================

    subroutine print_summary()

      type(Surface),  pointer :: s
      type(Cell),     pointer :: c
      type(Universe), pointer :: u
      type(Material), pointer :: m
      integer :: i
      integer :: ou

      ou = OUTPUT_UNIT

      ! print summary of cells
      write(ou,*) '============================================='
      write(ou,*) '=>              CELL SUMMARY               <='
      write(ou,*) '============================================='
      write(ou,*)
      do i = 1, n_cells
         c => cells(i)
         call print_cell(c)
      end do

      ! print summary of universes
      write(ou,*) '============================================='
      write(ou,*) '=>             UNIVERSE SUMMARY            <='
      write(ou,*) '============================================='
      write(ou,*)
      do i = 1, n_universes
         u => universes(i)
         call print_universe(u)
      end do

      ! print summary of surfaces
      write(ou,*) '============================================='
      write(ou,*) '=>              SURFACE SUMMARY            <='
      write(ou,*) '============================================='
      write(ou,*)
      do i = 1, n_surfaces
         s => surfaces(i)
         call print_surface(s)
      end do

      ! print summary of materials
      write(ou,*) '============================================='
      write(ou,*) '=>             MATERIAL SUMMARY            <='
      write(ou,*) '============================================='
      write(ou,*)
      do i = 1, n_materials
         m => materials(i)
         call print_material(m)
      end do

    end subroutine print_summary

end module output
