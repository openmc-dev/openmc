module output

  use ISO_FORTRAN_ENV
  use global
  use types,           only: Cell, Universe, Surface
  use data_structures, only: dict_get_key
  use endf,            only: reaction_name

  implicit none

  integer :: ou = OUTPUT_UNIT
  integer :: eu = ERROR_UNIT

contains

!===============================================================================
! TITLE prints the main title banner as well as information about the program
! developers, version, and date/time which the problem was run.
!===============================================================================

  subroutine title()

    character(10) :: date
    character(8)  :: time

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
    write(ou,*) '     Developed At:    Massachusetts Institute of Technology'
    write(ou,*) '     Lead Developer:  Paul K. Romano'
    write(ou,100) VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
100 format (6X,"Version:",9X,I1,".",I1,".",I1)

    ! Write the date and time
    call get_today(date, time)
    write(ou,101) trim(date), trim(time)
101 format (6X,"Date/Time:",7X,A,1X,A)
    write(ou,*)

  end subroutine title

!===============================================================================
! ECHO_INPUT displays summary information about the problem about to be run
! after reading all the input.
!===============================================================================

  subroutine echo_input()

    character(32) :: string

    ! Display problem summary
    write(ou,*) '============================================='
    write(ou,*) '=>             PROBLEM SUMMARY             <='
    write(ou,*) '============================================='
    write(ou,*)
    if (problem_type == PROB_CRITICALITY) then
       write(ou,100) 'Problem type:', 'Criticality'
       write(ou,100) 'Number of Cycles:', int_to_str(n_cycles)
       write(ou,100) 'Number of Inactive Cycles:', int_to_str(n_inactive) 
    elseif (problem_type == PROB_SOURCE) then
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
100 format (1X,A,T35,A)

  end subroutine echo_input

!===============================================================================
! MESSAGE displays an informational message to the log file and the standard
! output stream.
!===============================================================================

  subroutine message(msg, level)

    character(*), intent(in) :: msg
    integer,      optional   :: level

    integer :: n_lines
    integer :: i

    ! Only allow master to print to screen
    if (.not. master .and. present(level)) return

    if (.not. present(level) .or. level <= verbosity) then
       n_lines = (len_trim(msg)-1)/79 + 1
       do i = 1, n_lines
          write(ou, fmt='(1X,A)') trim(msg(79*(i-1)+1:79*i))
       end do
    end if

  end subroutine message

!===============================================================================
! WARNING issues a warning to the user in the log file and the standard output
! stream.
!===============================================================================

  subroutine warning(msg)

    character(*), intent(in) :: msg

    integer :: n_lines
    integer :: i

    ! Only allow master to print to screen
    if (.not. master) return

    write(ou, fmt='(1X,A9)', advance='no') 'WARNING: '

    n_lines = (len_trim(msg)-1)/70 + 1
    do i = 1, n_lines
       if (i == 1) then
          write(ou, fmt='(A70)') msg(70*(i-1)+1:70*i)
       else
          write(ou, fmt='(10X,A70)') msg(70*(i-1)+1:70*i)
       end if
    end do

  end subroutine warning

!===============================================================================
! ERROR alerts the user that an error has been encountered and displays a
! message about the particular problem. Errors are considered 'fatal' and hence
! the program is aborted.
!===============================================================================

  subroutine error(msg)

    character(*), intent(in) :: msg

    integer :: n_lines
    integer :: i

    ! Only allow master to print to screen
    if (master) then
       write(eu, fmt='(1X,A7)', advance='no') 'ERROR: '

       n_lines = (len_trim(msg)-1)/72 + 1
       do i = 1, n_lines
          if (i == 1) then
             write(eu, fmt='(A72)') msg(72*(i-1)+1:72*i)
          else
             write(eu, fmt='(7X,A72)') msg(72*(i-1)+1:72*i)
          end if
       end do
       write(eu,*)
    end if

    ! All processors abort
    call free_memory()

  end subroutine error

!===============================================================================
! GET_TODAY determines the date and time at which the program began execution
! and returns it in a readable format
!===============================================================================

  subroutine get_today(today_date, today_time)

    character(10), intent(out) :: today_date
    character(8),  intent(out) :: today_time

    integer       :: val(8)
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone

    call date_and_time(date, time, zone, val)
    ! val(1) = year (YYYY)
    ! val(2) = month (MM)
    ! val(3) = day (DD)
    ! val(4) = timezone
    ! val(5) = hours (HH)
    ! val(6) = minutes (MM)
    ! val(7) = seconds (SS)
    ! val(8) = milliseconds

    if (val(2) < 10) then
       if (val(3) < 10) then
          today_date = date(6:6) // "/" // date(8:8) // "/" // date(1:4)
       else
          today_date = date(6:6) // "/" // date(7:8) // "/" // date(1:4)
       end if
    else
       if (val(3) < 10) then
          today_date = date(5:6) // "/" // date(8:8) // "/" // date(1:4)
       else
          today_date = date(5:6) // "/" // date(7:8) // "/" // date(1:4)
       end if
    end if
    today_time = time(1:2) // ":" // time(3:4) // ":" // time(5:6)

  end subroutine get_today

!===============================================================================
! PRINT_PARTICLE displays the attributes of a particle
!===============================================================================

  subroutine print_particle(p)

    type(Particle), pointer :: p

    integer :: i
    character(250) :: string
    type(Cell),     pointer :: c => null()
    type(Surface),  pointer :: s => null()
    type(Universe), pointer :: u => null()

    select case (p % type)
    case (NEUTRON)
       write(ou,*) 'Neutron ' // int_to_str(p % uid)
    case (PHOTON)
       write(ou,*) 'Photon ' // int_to_str(p % uid)
    case (ELECTRON)
       write(ou,*) 'Electron ' // int_to_str(p % uid)
    case default
       write(ou,*) 'Unknown Particle ' // int_to_str(p % uid)
    end select
    write(ou,*) '    x = ' // real_to_str(p % xyz(1))
    write(ou,*) '    y = ' // real_to_str(p % xyz(2))
    write(ou,*) '    z = ' // real_to_str(p % xyz(3))
    write(ou,*) '    x local = ' // real_to_str(p % xyz_local(1))
    write(ou,*) '    y local = ' // real_to_str(p % xyz_local(2))
    write(ou,*) '    z local = ' // real_to_str(p % xyz_local(3))
    write(ou,*) '    u = ' // real_to_str(p % uvw(1))
    write(ou,*) '    v = ' // real_to_str(p % uvw(2))
    write(ou,*) '    w = ' // real_to_str(p % uvw(3))
    write(ou,*) '    Weight = ' // real_to_str(p % wgt)
    write(ou,*) '    Energy = ' // real_to_str(p % E)
    write(ou,*) '    x index = ' // int_to_str(p % index_x)
    write(ou,*) '    y index = ' // int_to_str(p % index_y)
    write(ou,*) '    IE = ' // int_to_str(p % IE)
    write(ou,*) '    Interpolation factor = ' // real_to_str(p % interp)

    if (p % cell > 0) then
       c => cells(p % cell)
       write(ou,*) '    Cell = ' // int_to_str(c % uid)
    else
       write(ou,*) '    Cell not determined'
    end if

    if (p % surface > 0) then
       s => surfaces(p % surface)
       write(ou,*) '    Surface = ' // int_to_str(s % uid)
    else
       write(ou,*) '    Surface = None'
    end if

    u => universes(p % universe)
    write(ou,*) '    Universe = ' // int_to_str(u % uid)
    write(ou,*)

    nullify(c)
    nullify(s)
    nullify(u)

  end subroutine print_particle

!===============================================================================
! PRINT_REACTION displays the attributes of a reaction
!===============================================================================

  subroutine print_reaction(rxn)

    type(AceReaction), pointer :: rxn

    write(ou,*) 'Reaction ' // reaction_name(rxn % MT)
    write(ou,*) '    MT = ' // int_to_str(rxn % MT)
    write(ou,*) '    Q-value = ' // real_to_str(rxn % Q_value)
    write(ou,*) '    TY = ' // int_to_str(rxn % TY)
    write(ou,*) '    Starting index = ' // int_to_str(rxn % IE)
    if (rxn % has_energy_dist) then
       write(ou,*) '    Energy: Law ' // int_to_str(rxn % edist % law)
    end if
    write(ou,*)

  end subroutine print_reaction

!===============================================================================
! PRINT_CELL displays the attributes of a cell
!===============================================================================

  subroutine print_cell(c)

    type(Cell), pointer :: c

    integer :: temp
    integer :: i
    character(250) :: string
    type(Universe), pointer :: u => null()
    type(Lattice),  pointer :: l => null()
    type(Material), pointer :: m => null()

    write(ou,*) 'Cell ' // int_to_str(c % uid)
    temp = dict_get_key(cell_dict, c % uid)
    write(ou,*) '    Array Index = ' // int_to_str(temp)
    u => universes(c % universe)
    write(ou,*) '    Universe = ' // int_to_str(u % uid)
    select case (c % type)
    case (CELL_NORMAL)
       write(ou,*) '    Fill = NONE'
    case (CELL_FILL)
       u => universes(c % fill)
       write(ou,*) '    Fill = Universe ' // int_to_str(u % uid)
    case (CELL_LATTICE)
       l => lattices(c % fill)
       write(ou,*) '    Fill = Lattice ' // int_to_str(l % uid)
    end select
    if (c % material == 0) then
       write(ou,*) '    Material = NONE'
    else
       m => materials(c % material)
       write(ou,*) '    Material = ' // int_to_str(m % uid)
    end if
    write(ou,*) '    Parent Cell = ' // int_to_str(c % parent)
    string = ""
    do i = 1, c % n_surfaces
       select case (c % surfaces(i))
       case (OP_LEFT_PAREN)
          string = trim(string) // ' ('
       case (OP_RIGHT_PAREN)
          string = trim(string) // ' )'
       case (OP_UNION)
          string = trim(string) // ' :'
       case (OP_DIFFERENCE)
          string = trim(string) // ' !'
       case default
          string = trim(string) // ' ' // int_to_str(c % surfaces(i))
       end select
    end do
    write(ou,*) '    Surface Specification:' // trim(string)
    write(ou,*)

    ! nullify associated pointers
    nullify(u)
    nullify(m)

  end subroutine print_cell

!===============================================================================
! PRINT_UNIVERSE displays the attributes of a universe
!===============================================================================

  subroutine print_universe(univ)

    type(Universe), pointer :: univ

    integer :: i
    character(250) :: string
    type(Cell), pointer :: c => null()

    write(ou,*) 'Universe ' // int_to_str(univ % uid)
    write(ou,*) '    Level = ' // int_to_str(univ % level)
    string = ""
    do i = 1, univ % n_cells
       c => cells(univ % cells(i))
       string = trim(string) // ' ' // int_to_str(c % uid)
    end do
    write(ou,*) '    Cells =' // trim(string)
    write(ou,*)

    nullify(c)

  end subroutine print_universe

!===============================================================================
! PRINT_LATTICE displays the attributes of a lattice
!===============================================================================

  subroutine print_lattice(lat)

    type(Lattice), pointer :: lat

    write(ou,*) 'Lattice ' // int_to_str(lat % uid)
    write(ou,*) '    n_x = ' // int_to_str(lat % n_x)
    write(ou,*) '    n_y = ' // int_to_str(lat % n_y)
    write(ou,*) '    x0 = ' // real_to_str(lat % x0)
    write(ou,*) '    y0 = ' // real_to_str(lat % y0)
    write(ou,*) '    width_x = ' // real_to_str(lat % width_x)
    write(ou,*) '    width_y = ' // real_to_str(lat % width_y)
    write(ou,*)

  end subroutine print_lattice

!===============================================================================
! PRINT_SURFACE displays the attributes of a surface
!===============================================================================

  subroutine print_surface(surf)

    type(Surface), pointer :: surf

    integer :: i
    character(80) :: string

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

    string = ""
    do i = 1, size(surf % coeffs)
       string = trim(string) // ' ' // real_to_str(surf % coeffs(i), 4)
    end do
    write(ou,*) '    Coefficients = ' // trim(string)

    string = ""
    if (allocated(surf % neighbor_pos)) then
       do i = 1, size(surf % neighbor_pos)
          string = trim(string) // ' ' // int_to_str(surf % neighbor_pos(i))
       end do
    end if
    write(ou,*) '    Positive Neighbors = ' // trim(string)

    string = ""
    if (allocated(surf % neighbor_neg)) then
       do i = 1, size(surf % neighbor_neg)
          string = trim(string) // ' ' // int_to_str(surf % neighbor_neg(i))
       end do
    end if
    write(ou,*) '    Negative Neighbors =' // trim(string)
    select case (surf % bc)
    case (BC_TRANSMIT)
       write(ou,*) '    Boundary Condition = Transmission'
    case (BC_VACUUM)
       write(ou,*) '    Boundary Condition = Vacuum'
    case (BC_REFLECT)
       write(ou,*) '    Boundary Condition = Reflective'
    case (BC_PERIODIC)
       write(ou,*) '    Boundary Condition = Periodic'
    end select
    write(ou,*)

  end subroutine print_surface

!===============================================================================
! PRINT_MATERIAL displays the attributes of a material
!===============================================================================

  subroutine print_material(mat)

    type(Material), pointer :: mat

    integer        :: i
    integer        :: n_lines
    real(8)        :: density
    character(250) :: string
    type(AceContinuous), pointer :: table => null()

    write(ou,*) 'Material ' // int_to_str(mat % uid)
    write(ou,*) '    Atom Density = ' // trim(real_to_str(mat % atom_density)) &
         & // ' atom/b-cm'
    do i = 1, mat % n_isotopes
       table => xs_continuous(mat % table(i))
       density = mat % atom_density * mat % atom_percent(i)
       string = '    ' // trim(table % name) // ' = ' // &
            & trim(real_to_str(density)) // ' atom/b-cm'
       write(ou,*) trim(string)
    end do
    write(ou,*)

    nullify(table)

  end subroutine print_material

!===============================================================================
! PRINT_TALLY displays the attributes of a tally
!===============================================================================

  subroutine print_tally(tal)

    type(Tally), pointer :: tal

    integer :: i
    integer :: MT
    character(250) :: string

    write(ou,*) 'Tally ' // int_to_str(tal % uid)

    select case (tal % reaction_type)
    case (TALLY_FLUX)
       write(ou,*) '    Type: Flux'
    case (TALLY_ALL)
       write(ou,*) '    Type: Total Collision Rate'
    case (TALLY_BINS)
       write(ou,*) '    Type: Partial reactions'
    case (TALLY_SUM)
       write(ou,*) '    Type: Partial reactions (summed)'
    end select

    select case (tal % cell_type)
    case (TALLY_BINS)
       write(ou,*) '    Cell Type: Separate bins'
    case (TALLY_SUM)
       write(ou,*) '    Cell Type: Sum over cells'
    end select

    if (allocated(tal % reactions)) then
       string = ""
       do i = 1, size(tal % reactions)
          MT = tal % reactions(i)
          string = trim(string) // ' ' // trim(reaction_name(MT))
       end do
       write(ou,*) '    Reactions:' // trim(string)
    end if

    if (allocated(cells)) then
       string = ""
       do i = 1, size(tal % cells)
          string = trim(string) // ' ' // trim(int_to_str(tal % cells(i)))
       end do
       write(ou,*) '    Cells:' // trim(string)
    end if
    
    if (allocated(tal % energies)) then
       string = ""
       do i = 1, size(tal % energies)
          string = trim(string) // ' ' // trim(real_to_str(tal % energies(i)))
       end do
       write(ou,*) '    Energies:' // trim(string)
    end if
    write(ou,*)

  end subroutine print_tally

!===============================================================================
! PRINT_SUMMARY displays the attributes of all cells, universes,
! surfaces and materials read in the input file. Very useful for
! debugging!
!===============================================================================

  subroutine print_summary()

    type(Surface),  pointer :: s => null()
    type(Cell),     pointer :: c => null()
    type(Universe), pointer :: u => null()
    type(Lattice),  pointer :: l => null()
    type(Material), pointer :: m => null()
    type(Tally),    pointer :: t => null()
    integer :: i

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

    ! print summary of lattices
    if (n_lattices > 0) then
       write(ou,*) '============================================='
       write(ou,*) '=>              LATTICE SUMMARY            <='
       write(ou,*) '============================================='
       write(ou,*)
       do i = 1, n_lattices
          l => lattices(i)
          call print_lattice(l)
       end do
    end if

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

    ! print summary of tallies
    if (n_tallies > 0) then
       write(ou,*) '============================================='
       write(ou,*) '=>               TALLY SUMMARY             <='
       write(ou,*) '============================================='
       write(ou,*)
       do i = 1, n_tallies
          t=> tallies(i)
          call print_tally(t)
       end do
    end if

    nullify(s)
    nullify(c)
    nullify(u)
    nullify(l)
    nullify(m)
    nullify(t)

  end subroutine print_summary

end module output
