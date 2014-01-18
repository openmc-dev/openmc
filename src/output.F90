module output

  use, intrinsic :: ISO_FORTRAN_ENV

  use ace_header,      only: Nuclide, Reaction, UrrData
  use constants
  use endf,            only: reaction_name
  use error,           only: warning
  use geometry_header, only: Cell, Universe, Surface, BASE_UNIVERSE
  use global
  use math,            only: t_percentile
  use mesh_header,     only: StructuredMesh
  use mesh,            only: mesh_indices_to_bin, bin_to_mesh_indices
  use particle_header, only: LocalCoord, Particle
  use plot_header
  use string,          only: upper_case, to_str
  use tally_header,    only: TallyObject

  implicit none

  ! Short names for output and error units
  integer :: ou = OUTPUT_UNIT
  integer :: eu = ERROR_UNIT

contains

!===============================================================================
! TITLE prints the main title banner as well as information about the program
! developers, version, and date/time which the problem was run.
!===============================================================================

  subroutine title()

#ifdef _OPENMP
    use omp_lib
#endif

    write(UNIT=OUTPUT_UNIT, FMT='(/11(A/))') &
         '       .d88888b.                             888b     d888  .d8888b.', &
         '      d88P" "Y88b                            8888b   d8888 d88P  Y88b', &
         '      888     888                            88888b.d88888 888    888', &
         '      888     888 88888b.   .d88b.  88888b.  888Y88888P888 888       ', &
         '      888     888 888 "88b d8P  Y8b 888 "88b 888 Y888P 888 888       ', &
         '      888     888 888  888 88888888 888  888 888  Y8P  888 888    888', &
         '      Y88b. .d88P 888 d88P Y8b.     888  888 888   "   888 Y88b  d88P', &
         '       "Y88888P"  88888P"   "Y8888  888  888 888       888  "Y8888P"', &
         '__________________888______________________________________________________', &
         '                  888', &
         '                  888'

    ! Write version information
    write(UNIT=OUTPUT_UNIT, FMT=*) &
         '     Copyright:      2011-2014 Massachusetts Institute of Technology'
    write(UNIT=OUTPUT_UNIT, FMT=*) &
         '     License:        http://mit-crpg.github.io/openmc/license.html'
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Version:",8X,I1,".",I1,".",I1)') &
         VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
#ifdef GIT_SHA1
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Git SHA1:",7X,A)') GIT_SHA1
#endif

    ! Write the date and time
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Date/Time:",6X,A)') &
         time_stamp()

#ifdef MPI
    ! Write number of processors
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"MPI Processes:",2X,A)') &
         trim(to_str(n_procs))
#endif

#ifdef _OPENMP
    ! Write number of OpenMP threads
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"OpenMP Threads:",1X,A)') &
         trim(to_str(omp_get_max_threads()))
#endif

  end subroutine title

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
! specified, it is assumed to be a minor header block (H3).
!===============================================================================

  subroutine header(msg, unit, level)

    character(*), intent(in) :: msg ! header message
    integer, optional :: unit       ! unit to write to
    integer, optional :: level      ! specified header level

    integer :: n            ! number of = signs on left
    integer :: m            ! number of = signs on right
    integer :: unit_        ! unit to write to
    integer :: header_level ! actual header level
    character(MAX_LINE_LEN) :: line

    ! set default level
    if (present(level)) then
      header_level = level
    else
      header_level = 3
    end if

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
    line = msg
    call upper_case(line)

    ! print header based on level
    select case (header_level)
    case (1)
      write(UNIT=unit_, FMT='(/3(1X,A/))') repeat('=', 75), &
           repeat('=', n) // '>     ' // trim(line) // '     <' // &
           repeat('=', m), repeat('=', 75)
    case (2)
      write(UNIT=unit_, FMT='(/2(1X,A/))') trim(line), repeat('-', 75)
    case (3)
      write(UNIT=unit_, FMT='(/1X,A/)') repeat('=', n) // '>     ' // &
           trim(line) // '     <' // repeat('=', m)
    end select

  end subroutine header

!===============================================================================
! PRINT_VERSION shows the current version as well as copright and license
! information
!===============================================================================

  subroutine print_version()

    if (master) then
      write(UNIT=OUTPUT_UNIT, FMT='(1X,A,1X,I1,".",I1,".",I1)') &
           "OpenMC version", VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
      write(UNIT=OUTPUT_UNIT, FMT=*) "Copyright (c) 2011-2013 &
           &Massachusetts Institute of Technology"
      write(UNIT=OUTPUT_UNIT, FMT=*) "MIT/X license at &
           &<http://mit-crpg.github.io/openmc/license.html>"
    end if

  end subroutine print_version

!===============================================================================
! PRINT_USAGE displays information about command line usage of OpenMC
!===============================================================================

  subroutine print_usage()

    if (master) then
      write(OUTPUT_UNIT,*) 'Usage: openmc [options] [directory]'
      write(OUTPUT_UNIT,*)
      write(OUTPUT_UNIT,*) 'Options:'
      write(OUTPUT_UNIT,*) '  -g, --geometry-debug   Run in geometry debugging mode'
      write(OUTPUT_UNIT,*) '  -n, --particles        Number of particles per generation'
      write(OUTPUT_UNIT,*) '  -p, --plot             Run in plotting mode'
      write(OUTPUT_UNIT,*) '  -r, --restart          Restart a previous run from a state point'
      write(OUTPUT_UNIT,*) '                         or a particle restart file'
      write(OUTPUT_UNIT,*) '  -s, --threads          Number of OpenMP threads'
      write(OUTPUT_UNIT,*) '  -t, --track            Write tracks for all particles'
      write(OUTPUT_UNIT,*) '  -v, --version          Show version information'
      write(OUTPUT_UNIT,*) '  -?, --help             Show this message'
    end if

  end subroutine print_usage

!===============================================================================
! WRITE_MESSAGE displays an informational message to the log file and the
! standard output stream.
!===============================================================================

  subroutine write_message(level)

    integer, optional :: level ! verbosity level

    integer :: i_start    ! starting position
    integer :: i_end      ! ending position
    integer :: line_wrap  ! length of line
    integer :: length     ! length of message
    integer :: last_space ! index of last space (relative to start)

    ! Set length of line
    line_wrap = 80

    ! Only allow master to print to screen
    if (.not. master .and. present(level)) return

    if (.not. present(level) .or. level <= verbosity) then
      ! Determine length of message
      length = len_trim(message)

      i_start = 0
      do
        if (length - i_start < line_wrap - 1) then
          ! Remainder of message will fit on line
          write(ou, fmt='(1X,A)') message(i_start+1:length)
          exit

        else
          ! Determine last space in current line
          last_space = index(message(i_start+1:i_start+line_wrap), &
               ' ', BACK=.true.)
          if (last_space == 0) then
            i_end = min(length + 1, i_start+line_wrap) - 1
            write(ou, fmt='(1X,A)') message(i_start+1:i_end)
          else
            i_end = i_start + last_space
            write(ou, fmt='(1X,A)') message(i_start+1:i_end-1)
          end if

          ! Write up to last space

          ! Advance starting position
          i_start = i_end
          if (i_start > length) exit
        end if
      end do
    end if

  end subroutine write_message

!===============================================================================
! PRINT_PARTICLE displays the attributes of a particle
!===============================================================================

  subroutine print_particle(p)

    type(Particle), intent(in) :: p

    integer :: i ! index for coordinate levels
    type(Cell),       pointer :: c => null()
    type(Surface),    pointer :: s => null()
    type(Universe),   pointer :: u => null()
    type(Lattice),    pointer :: l => null()
    type(LocalCoord), pointer :: coord => null()

    ! display type of particle
    select case (p % type)
    case (NEUTRON)
      write(ou,*) 'Neutron ' // to_str(p % id)
    case (PHOTON)
      write(ou,*) 'Photon ' // to_str(p % id)
    case (ELECTRON)
      write(ou,*) 'Electron ' // to_str(p % id)
    case default
      write(ou,*) 'Unknown Particle ' // to_str(p % id)
    end select

    ! loop through each level of universes
    coord => p % coord0
    i = 0
    do while(associated(coord))
      ! Print level
      write(ou,*) '  Level ' // trim(to_str(i))

      ! Print cell for this level
      if (coord % cell /= NONE) then
        c => cells(coord % cell)
        write(ou,*) '    Cell             = ' // trim(to_str(c % id))
      end if

      ! Print universe for this level
      if (coord % universe /= NONE) then
        u => universes(coord % universe)
        write(ou,*) '    Universe         = ' // trim(to_str(u % id))
      end if

      ! Print information on lattice
      if (coord % lattice /= NONE) then
        l => lattices(coord % lattice)
        write(ou,*) '    Lattice          = ' // trim(to_str(l % id))
        write(ou,*) '    Lattice position = (' // trim(to_str(&
             p % coord % lattice_x)) // ',' // trim(to_str(&
             p % coord % lattice_y)) // ')'
      end if

      ! Print local coordinates
      write(ou,'(1X,A,3ES12.4)') '    xyz = ', coord % xyz
      write(ou,'(1X,A,3ES12.4)') '    uvw = ', coord % uvw

      coord => coord % next
      i = i + 1
    end do

    ! Print surface
    if (p % surface /= NONE) then
      s => surfaces(abs(p % surface))
      write(ou,*) '  Surface = ' // to_str(sign(s % id, p % surface))
    end if

    ! Display weight, energy, grid index, and interpolation factor
    write(ou,*) '  Weight = ' // to_str(p % wgt)
    write(ou,*) '  Energy = ' // to_str(p % E)
    write(ou,*)

  end subroutine print_particle

!===============================================================================
! PRINT_REACTION displays the attributes of a reaction
!===============================================================================

  subroutine print_reaction(rxn)

    type(Reaction), pointer :: rxn

    write(ou,*) 'Reaction ' // reaction_name(rxn % MT)
    write(ou,*) '    MT = ' // to_str(rxn % MT)
    write(ou,*) '    Q-value = ' // to_str(rxn % Q_value)
    write(ou,*) '    Multiplicity = ' // to_str(rxn % multiplicity)
    write(ou,*) '    Threshold = ' // to_str(rxn % threshold)
    if (rxn % has_energy_dist) then
      write(ou,*) '    Energy: Law ' // to_str(rxn % edist % law)
    end if
    write(ou,*)

  end subroutine print_reaction

!===============================================================================
! PRINT_CELL displays the attributes of a cell
!===============================================================================

  subroutine print_cell(c, unit)

    type(Cell), pointer :: c
    integer,   optional :: unit ! specified unit to write to

    integer :: index_cell ! index in cells array
    integer :: i          ! loop index for surfaces
    integer :: index_surf ! index in surfaces array
    integer :: unit_      ! unit to write to
    character(MAX_LINE_LEN) :: string
    type(Universe), pointer :: u => null()
    type(Lattice),  pointer :: l => null()
    type(Material), pointer :: m => null()

    ! Set unit to stdout if not already set
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! Write user-specified id for cell
    write(unit_,*) 'Cell ' // to_str(c % id)

    ! Find index in cells array and write
    index_cell = cell_dict % get_key(c % id)
    write(unit_,*) '    Array Index = ' // to_str(index_cell)

    ! Write what universe this cell is in
    u => universes(c % universe)
    write(unit_,*) '    Universe = ' // to_str(u % id)

    ! Write information on fill for cell
    select case (c % type)
    case (CELL_NORMAL)
      write(unit_,*) '    Fill = NONE'
    case (CELL_FILL)
      u => universes(c % fill)
      write(unit_,*) '    Fill = Universe ' // to_str(u % id)
    case (CELL_LATTICE)
      l => lattices(c % fill)
      write(unit_,*) '    Fill = Lattice ' // to_str(l % id)
    end select

    ! Write information on material
    if (c % material == 0) then
      write(unit_,*) '    Material = NONE'
    elseif (c % material == MATERIAL_VOID) then
      write(unit_,*) '    Material = Void'
    else
      m => materials(c % material)
      write(unit_,*) '    Material = ' // to_str(m % id)
    end if

    ! Write surface specification
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
        index_surf = abs(c % surfaces(i))
        string = trim(string) // ' ' // to_str(sign(&
             surfaces(index_surf) % id, c % surfaces(i)))
      end select
    end do
    write(unit_,*) '    Surface Specification:' // trim(string)
    write(unit_,*)

  end subroutine print_cell

!===============================================================================
! PRINT_UNIVERSE displays the attributes of a universe
!===============================================================================

  subroutine print_universe(univ, unit)

    type(Universe), pointer :: univ
    integer,       optional :: unit

    integer :: i     ! loop index for cells in this universe
    integer :: unit_ ! unit to write to
    character(MAX_LINE_LEN) :: string
    type(Cell), pointer     :: c => null()
    type(Universe), pointer :: base_u => null()

    ! Set default unit to stdout if not specified
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! Get a pointer to the base universe
    base_u => universes(BASE_UNIVERSE)

    ! Write user-specified id for this universe
    write(unit_,*) 'Universe ' // to_str(univ % id)

    ! If this is the base universe, indicate so
    if (associated(univ, base_u)) then
      write(unit_,*) '    Base Universe'
    end if

    ! Write list of cells in this universe
    string = ""
    do i = 1, univ % n_cells
      c => cells(univ % cells(i))
      string = trim(string) // ' ' // to_str(c % id)
    end do
    write(unit_,*) '    Cells =' // trim(string)
    write(unit_,*)

  end subroutine print_universe

!===============================================================================
! PRINT_LATTICE displays the attributes of a lattice
!===============================================================================

  subroutine print_lattice(lat, unit)

    type(Lattice), pointer :: lat
    integer,      optional :: unit

    integer :: i     ! loop index
    integer :: unit_ ! unit to write to
    character(MAX_LINE_LEN) :: string


    ! set default unit if not specified
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! Write information about lattice
    write(unit_,*) 'Lattice ' // to_str(lat % id)

    ! Write dimension of lattice
    string = ""
    do i = 1, lat % n_dimension
      string = trim(string) // ' ' // to_str(lat % dimension(i))
    end do
    write(unit_,*) '    Dimension =' // string

    ! Write lower-left coordinates of lattice
    string = ""
    do i = 1, lat % n_dimension
      string = trim(string) // ' ' // to_str(lat % lower_left(i))
    end do
    write(unit_,*) '    Lower-left =' // string

    ! Write width of each lattice cell
    string = ""
    do i = 1, lat % n_dimension
      string = trim(string) // ' ' // to_str(lat % width(i))
    end do
    write(unit_,*) '    Width =' // string
    write(unit_,*)

  end subroutine print_lattice

!===============================================================================
! PRINT_SURFACE displays the attributes of a surface
!===============================================================================

  subroutine print_surface(surf, unit)

    type(Surface), pointer :: surf
    integer,      optional :: unit ! specified unit to write to

    integer :: i     ! loop index for coefficients
    integer :: unit_ ! unit to write to
    character(MAX_LINE_LEN) :: string
    type(Cell), pointer :: c => null()

    ! set default unit if not specified
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! Write user-specified id of surface
    write(unit_,*) 'Surface ' // to_str(surf % id)

    ! Write type of surface
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
    case (SURF_CONE_X)
      string = "X Cone"
    case (SURF_CONE_Y)
      string = "Y Cone"
    case (SURF_CONE_Z)
      string = "Z Cone"
    end select
    write(unit_,*) '    Type = ' // trim(string)

    ! Write coefficients for this surface
    string = ""
    do i = 1, size(surf % coeffs)
      string = trim(string) // ' ' // to_str(surf % coeffs(i), 4)
    end do
    write(unit_,*) '    Coefficients = ' // trim(string)

    ! Write neighboring cells on positive side of this surface
    string = ""
    if (allocated(surf % neighbor_pos)) then
      do i = 1, size(surf % neighbor_pos)
        c => cells(abs(surf % neighbor_pos(i)))
        string = trim(string) // ' ' // to_str(&
             sign(c % id, surf % neighbor_pos(i)))
      end do
    end if
    write(unit_,*) '    Positive Neighbors = ' // trim(string)

    ! Write neighboring cells on negative side of this surface
    string = ""
    if (allocated(surf % neighbor_neg)) then
      do i = 1, size(surf % neighbor_neg)
        c => cells(abs(surf % neighbor_neg(i)))
        string = trim(string) // ' ' // to_str(&
             sign(c % id, surf % neighbor_neg(i)))
      end do
    end if
    write(unit_,*) '    Negative Neighbors =' // trim(string)

    ! Write boundary condition for this surface
    select case (surf % bc)
    case (BC_TRANSMIT)
      write(unit_,*) '    Boundary Condition = Transmission'
    case (BC_VACUUM)
      write(unit_,*) '    Boundary Condition = Vacuum'
    case (BC_REFLECT)
      write(unit_,*) '    Boundary Condition = Reflective'
    case (BC_PERIODIC)
      write(unit_,*) '    Boundary Condition = Periodic'
    end select
    write(unit_,*)

  end subroutine print_surface

!===============================================================================
! PRINT_MATERIAL displays the attributes of a material
!===============================================================================

  subroutine print_material(mat, unit)

    type(Material), pointer :: mat
    integer,       optional :: unit

    integer :: i       ! loop index for nuclides
    integer :: unit_   ! unit to write to
    real(8) :: density ! density in atom/b-cm
    character(MAX_LINE_LEN) :: string
    type(Nuclide),  pointer :: nuc => null()

    ! set default unit to stdout if not specified
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! Write identifier for material
    write(unit_,*) 'Material ' // to_str(mat % id)

    ! Write total atom density in atom/b-cm
    write(unit_,*) '    Atom Density = ' // trim(to_str(mat % density)) &
         // ' atom/b-cm'

    ! Write atom density for each nuclide in material
    write(unit_,*) '    Nuclides:'
    do i = 1, mat % n_nuclides
      nuc => nuclides(mat % nuclide(i))
      density = mat % atom_density(i)
      string = '        ' // trim(nuc % name) // ' = ' // &
           trim(to_str(density)) // ' atom/b-cm'
      write(unit_,*) trim(string)
    end do

    ! Write information on S(a,b) table
    if (mat % n_sab > 0) then
        write(unit_,*) '    S(a,b) tables:'
      do i = 1, mat % n_sab
        write(unit_,*) '      ' // trim(&
             sab_tables(mat % i_sab_tables(i)) % name)
      end do
    end if
    write(unit_,*)

  end subroutine print_material

!===============================================================================
! PRINT_TALLY displays the attributes of a tally
!===============================================================================

  subroutine print_tally(t, unit)

    type(TallyObject), pointer :: t
    integer,          optional :: unit

    integer :: i     ! index for filter or score bins
    integer :: j     ! index in filters array
    integer :: id    ! user-specified id
    integer :: unit_ ! unit to write to
    integer :: n     ! scattering order to include in name
    character(MAX_LINE_LEN) :: string
    character(MAX_WORD_LEN) :: pn_string
    type(Cell),           pointer :: c => null()
    type(Surface),        pointer :: s => null()
    type(Universe),       pointer :: u => null()
    type(Material),       pointer :: m => null()
    type(StructuredMesh), pointer :: sm => null()

    ! set default unit to stdout if not specified
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! Write user-specified id of tally
    write(unit_,*) 'Tally ' // to_str(t % id)

    ! Write the type of tally
    select case(t % type)
    case (TALLY_VOLUME)
      write(unit_,*) '    Type: Volume'
    case (TALLY_SURFACE_CURRENT)
      write(unit_,*) '    Type: Surface Current'
    end select

    ! Write the estimator used
    select case(t % estimator)
    case(ESTIMATOR_ANALOG)
      write(unit_,*) '    Estimator: Analog'
    case(ESTIMATOR_TRACKLENGTH)
      write(unit_,*) '    Estimator: Track-length'
    end select

    ! Write any cells bins if present
    j = t % find_filter(FILTER_CELL)
    if (j > 0) then
      string = ""
      do i = 1, t % filters(j) % n_bins
        id = t % filters(j) % int_bins(i)
        c => cells(id)
        string = trim(string) // ' ' // trim(to_str(c % id))
      end do
      write(unit_, *) '    Cell Bins:' // trim(string)
    end if

    ! Write any surface bins if present
    j = t % find_filter(FILTER_SURFACE)
    if (j > 0) then
      string = ""
      do i = 1, t % filters(j) % n_bins
        id = t % filters(j) % int_bins(i)
        s => surfaces(id)
        string = trim(string) // ' ' // trim(to_str(s % id))
      end do
      write(unit_, *) '    Surface Bins:' // trim(string)
    end if

    ! Write any universe bins if present
    j = t % find_filter(FILTER_UNIVERSE)
    if (j > 0) then
      string = ""
      do i = 1, t % filters(j) % n_bins
        id = t % filters(j) % int_bins(i)
        u => universes(id)
        string = trim(string) // ' ' // trim(to_str(u % id))
      end do
      write(unit_, *) '    Universe Bins:' // trim(string)
    end if

    ! Write any material bins if present
    j = t % find_filter(FILTER_MATERIAL)
    if (j > 0) then
      string = ""
      do i = 1, t % filters(j) % n_bins
        id = t % filters(j) % int_bins(i)
        m => materials(id)
        string = trim(string) // ' ' // trim(to_str(m % id))
      end do
      write(unit_, *) '    Material Bins:' // trim(string)
    end if

    ! Write any mesh bins if present
    j = t % find_filter(FILTER_MESH)
    if (j > 0) then
      string = ""
      id = t % filters(j) % int_bins(1)
      sm => meshes(id)
      string = trim(string) // ' ' // trim(to_str(sm % dimension(1)))
      do i = 2, sm % n_dimension
        string = trim(string) // ' x ' // trim(to_str(sm % dimension(i)))
      end do
      write(unit_, *) '    Mesh Bins:' // trim(string)
    end if

    ! Write any birth region bins if present
    j = t % find_filter(FILTER_CELLBORN)
    if (j > 0) then
      string = ""
      do i = 1, t % filters(j) % n_bins
        id = t % filters(j) % int_bins(i)
        c => cells(id)
        string = trim(string) // ' ' // trim(to_str(c % id))
      end do
      write(unit_, *) '    Birth Region Bins:' // trim(string)
    end if

    ! Write any incoming energy bins if present
    j = t % find_filter(FILTER_ENERGYIN)
    if (j > 0) then
      string = ""
      do i = 1, t % filters(j) % n_bins + 1
        string = trim(string) // ' ' // trim(to_str(&
             t % filters(j) % real_bins(i)))
      end do
      write(unit_,*) '    Incoming Energy Bins:' // trim(string)
    end if

    ! Write any outgoing energy bins if present
    j = t % find_filter(FILTER_ENERGYOUT)
    if (j > 0) then
      string = ""
      do i = 1, t % filters(j) % n_bins + 1
        string = trim(string) // ' ' // trim(to_str(&
             t % filters(j) % real_bins(i)))
      end do
      write(unit_,*) '    Outgoing Energy Bins:' // trim(string)
    end if

    ! Write nuclides bins
    write(unit_,fmt='(1X,A)',advance='no') '    Nuclide Bins:'
    do i = 1, t % n_nuclide_bins
      if (t % nuclide_bins(i) == -1) then
        write(unit_,fmt='(A)',advance='no') ' total'
      else
        write(unit_,fmt='(A)',advance='no') ' ' // trim(adjustl(&
             nuclides(t % nuclide_bins(i)) % name))
      end if
      if (mod(i,4) == 0 .and. i /= t % n_nuclide_bins) &
           write(unit_,'(/18X)',advance='no')
    end do
    write(unit_,*)


    ! Write score bins
    string   = ""
    j = 0
    do i = 1, t % n_user_score_bins
      j = j + 1
      select case (t % score_bins(j))
      case (SCORE_FLUX)
        string = trim(string) // ' flux'
      case (SCORE_TOTAL)
        string = trim(string) // ' total'
      case (SCORE_SCATTER)
        string = trim(string) // ' scatter'
      case (SCORE_NU_SCATTER)
        string = trim(string) // ' nu-scatter'
      case (SCORE_SCATTER_N)
        pn_string = ' scatter-' // trim(to_str(t % scatt_order(j)))
        string = trim(string) // pn_string
      case (SCORE_SCATTER_PN)
        pn_string = ' scatter'
        string = trim(string) // pn_string
        do n = 1, t % scatt_order(j)
          pn_string = ' scatter-' // trim(to_str(n))
          string = trim(string) // pn_string
        end do
        j = j + n - 1
      case (SCORE_TRANSPORT)
        string = trim(string) // ' transport'
      case (SCORE_N_1N)
        string = trim(string) // ' n1n'
      case (SCORE_ABSORPTION)
        string = trim(string) // ' absorption'
      case (SCORE_FISSION)
        string = trim(string) // ' fission'
      case (SCORE_NU_FISSION)
        string = trim(string) // ' nu-fission'
      case (SCORE_KAPPA_FISSION)
        string = trim(string) // ' kappa-fission'
      case (SCORE_CURRENT)
        string = trim(string) // ' current'
      case default
        string = trim(string) // ' ' // reaction_name(t % score_bins(j))
      end select
    end do
    write(unit_,*) '    Scores:' // trim(string)
    write(unit_,*)

  end subroutine print_tally

!===============================================================================
! PRINT_GEOMETRY displays the attributes of all cells, surfaces, universes,
! surfaces, and lattices read in the input files.
!===============================================================================

  subroutine print_geometry()

    integer :: i ! loop index for various arrays
    type(Surface),     pointer :: s => null()
    type(Cell),        pointer :: c => null()
    type(Universe),    pointer :: u => null()
    type(Lattice),     pointer :: l => null()

    ! print summary of surfaces
    call header("SURFACE SUMMARY", unit=UNIT_SUMMARY)
    do i = 1, n_surfaces
      s => surfaces(i)
      call print_surface(s, unit=UNIT_SUMMARY)
    end do

    ! print summary of cells
    call header("CELL SUMMARY", unit=UNIT_SUMMARY)
    do i = 1, n_cells
      c => cells(i)
      call print_cell(c, unit=UNIT_SUMMARY)
    end do

    ! print summary of universes
    call header("UNIVERSE SUMMARY", unit=UNIT_SUMMARY)
    do i = 1, n_universes
      u => universes(i)
      call print_universe(u, unit=UNIT_SUMMARY)
    end do

    ! print summary of lattices
    if (n_lattices > 0) then
      call header("LATTICE SUMMARY", unit=UNIT_SUMMARY)
      do i = 1, n_lattices
        l => lattices(i)
        call print_lattice(l, unit=UNIT_SUMMARY)
      end do
    end if

  end subroutine print_geometry

!===============================================================================
! PRINT_NUCLIDE displays information about a continuous-energy neutron
! cross_section table and its reactions and secondary angle/energy distributions
!===============================================================================

  subroutine print_nuclide(nuc, unit)

    type(Nuclide), pointer :: nuc
    integer,      optional :: unit

    integer :: i                 ! loop index over nuclides
    integer :: unit_             ! unit to write to
    integer :: size_total        ! memory used by nuclide (bytes)
    integer :: size_angle_total  ! total memory used for angle dist. (bytes)
    integer :: size_energy_total ! total memory used for energy dist. (bytes)
    integer :: size_xs           ! memory used for cross-sections (bytes)
    integer :: size_angle        ! memory used for an angle distribution (bytes)
    integer :: size_energy       ! memory used for a  energy distributions (bytes)
    integer :: size_urr          ! memory used for probability tables (bytes)
    character(11) :: law         ! secondary energy distribution law
    type(Reaction), pointer :: rxn => null()
    type(UrrData),  pointer :: urr => null()

    ! set default unit for writing information
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! Initialize totals
    size_angle_total = 0
    size_energy_total = 0
    size_urr = 0
    size_xs = 0

    ! Basic nuclide information
    write(unit_,*) 'Nuclide ' // trim(nuc % name)
    write(unit_,*) '  zaid = ' // trim(to_str(nuc % zaid))
    write(unit_,*) '  awr = ' // trim(to_str(nuc % awr))
    write(unit_,*) '  kT = ' // trim(to_str(nuc % kT))
    write(unit_,*) '  # of grid points = ' // trim(to_str(nuc % n_grid))
    write(unit_,*) '  Fissionable = ', nuc % fissionable
    write(unit_,*) '  # of fission reactions = ' // trim(to_str(nuc % n_fission))
    write(unit_,*) '  # of reactions = ' // trim(to_str(nuc % n_reaction))

    ! Information on each reaction
    write(unit_,*) '  Reaction     Q-value  COM  Law    IE    size(angle) size(energy)'
    do i = 1, nuc % n_reaction
      rxn => nuc % reactions(i)

      ! Determine size of angle distribution
      if (rxn % has_angle_dist) then
        size_angle = rxn % adist % n_energy * 16 + size(rxn % adist % data) * 8
      else
        size_angle = 0
      end if

      ! Determine size of energy distribution and law
      if (rxn % has_energy_dist) then
        size_energy = size(rxn % edist % data) * 8
        law = to_str(rxn % edist % law)
      else
        size_energy = 0
        law = 'None'
      end if

      write(unit_,'(3X,A11,1X,F8.3,3X,L1,3X,A4,1X,I6,1X,I11,1X,I11)') &
           reaction_name(rxn % MT), rxn % Q_value, rxn % scatter_in_cm, &
           law(1:4), rxn % threshold, size_angle, size_energy

      ! Accumulate data size
      size_xs = size_xs + (nuc % n_grid - rxn%threshold + 1) * 8
      size_angle_total = size_angle_total + size_angle
      size_energy_total = size_energy_total + size_energy
    end do

    ! Add memory required for summary reactions (total, absorption, fission,
    ! nu-fission)
    size_xs = 8 * nuc % n_grid * 4

    ! Write information about URR probability tables
    size_urr = 0
    if (nuc % urr_present) then
      urr => nuc % urr_data
      write(unit_,*) '  Unresolved resonance probability table:'
      write(unit_,*) '    # of energies = ' // trim(to_str(urr % n_energy))
      write(unit_,*) '    # of probabilities = ' // trim(to_str(urr % n_prob))
      write(unit_,*) '    Interpolation =  ' // trim(to_str(urr % interp))
      write(unit_,*) '    Inelastic flag = ' // trim(to_str(urr % inelastic_flag))
      write(unit_,*) '    Absorption flag = ' // trim(to_str(urr % absorption_flag))
      write(unit_,*) '    Multiply by smooth? ', urr % multiply_smooth
      write(unit_,*) '    Min energy = ', trim(to_str(urr % energy(1)))
      write(unit_,*) '    Max energy = ', trim(to_str(urr % energy(urr % n_energy)))

      ! Calculate memory used by probability tables and add to total
      size_urr = urr % n_energy * (urr % n_prob * 6 + 1) * 8
    end if

    ! Calculate total memory
    size_total = size_xs + size_angle_total + size_energy_total + size_urr

    ! Write memory used
    write(unit_,*) '  Memory Requirements'
    write(unit_,*) '    Cross sections = ' // trim(to_str(size_xs)) // ' bytes'
    write(unit_,*) '    Secondary angle distributions = ' // &
         trim(to_str(size_angle_total)) // ' bytes'
    write(unit_,*) '    Secondary energy distributions = ' // &
         trim(to_str(size_energy_total)) // ' bytes'
    write(unit_,*) '    Probability Tables = ' // &
         trim(to_str(size_urr)) // ' bytes'
    write(unit_,*) '    Total = ' // trim(to_str(size_total)) // ' bytes'

    ! Blank line at end of nuclide
    write(unit_,*)

  end subroutine print_nuclide

!===============================================================================
! PRINT_SAB_TABLE displays information about a S(a,b) table containing data
! describing thermal scattering from bound materials such as hydrogen in water.
!===============================================================================

  subroutine print_sab_table(sab, unit)

    type(SAlphaBeta), pointer :: sab
    integer,         optional :: unit

    integer :: size_sab   ! memory used by S(a,b) table
    integer :: unit_      ! unit to write to
    integer :: i          ! Loop counter for parsing through sab % zaid
    integer :: char_count ! Counter for the number of characters on a line

    ! set default unit for writing information
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! Basic S(a,b) table information
    write(unit_,*) 'S(a,b) Table ' // trim(sab % name)
    write(unit_,'(A)',advance="no") '   zaids = '
    ! Initialize the counter based on the above string
    char_count = 11
    do i = 1, sab % n_zaid
      ! Deal with a line thats too long
      if (char_count >= 73) then  ! 73 = 80 - (5 ZAID chars + 1 space + 1 comma)
        ! End the line
        write(unit_,*) ""
        ! Add 11 leading blanks
        write(unit_,'(A)', advance="no") "           "
        ! reset the counter to 11
        char_count = 11
      end if
      if (i < sab % n_zaid) then
        ! Include a comma
        write(unit_,'(A)',advance="no") trim(to_str(sab % zaid(i))) // ", "
        char_count = char_count + len(trim(to_str(sab % zaid(i)))) + 2
      else
        ! Don't include a comma, since we are all done
        write(unit_,'(A)',advance="no") trim(to_str(sab % zaid(i)))
      end if

    end do
    write(unit_,*) "" ! Move to next line
    write(unit_,*) '  awr = ' // trim(to_str(sab % awr))
    write(unit_,*) '  kT = ' // trim(to_str(sab % kT))

    ! Inelastic data
    write(unit_,*) '  # of Incoming Energies (Inelastic) = ' // &
         trim(to_str(sab % n_inelastic_e_in))
    write(unit_,*) '  # of Outgoing Energies (Inelastic) = ' // &
         trim(to_str(sab % n_inelastic_e_out))
    write(unit_,*) '  # of Outgoing Angles (Inelastic) = ' // &
         trim(to_str(sab % n_inelastic_mu))
    write(unit_,*) '  Threshold for Inelastic = ' // &
         trim(to_str(sab % threshold_inelastic))

    ! Elastic data
    if (sab % n_elastic_e_in > 0) then
      write(unit_,*) '  # of Incoming Energies (Elastic) = ' // &
           trim(to_str(sab % n_elastic_e_in))
      write(unit_,*) '  # of Outgoing Angles (Elastic) = ' // &
           trim(to_str(sab % n_elastic_mu))
      write(unit_,*) '  Threshold for Elastic = ' // &
           trim(to_str(sab % threshold_elastic))
    end if

    ! Determine memory used by S(a,b) table and write out
    size_sab = 8 * (sab % n_inelastic_e_in * (2 + sab % n_inelastic_e_out * &
         (1 + sab % n_inelastic_mu)) + sab % n_elastic_e_in * &
         (2 + sab % n_elastic_mu))
    write(unit_,*) '  Memory Used = ' // trim(to_str(size_sab)) // ' bytes'

    ! Blank line at end
    write(unit_,*)

  end subroutine print_sab_table

!===============================================================================
! WRITE_SUMMARY displays summary information about the problem about to be run
! after reading all input files
!===============================================================================

  subroutine write_summary()

    integer                 :: i      ! loop index
    character(MAX_FILE_LEN) :: path   ! path of summary file
    type(Material),    pointer :: m => null()
    type(TallyObject), pointer :: t => null()

    ! Create filename for log file
    path = trim(path_output) // "summary.out"

    ! Open log file for writing
    open(UNIT=UNIT_SUMMARY, FILE=path, STATUS='replace', ACTION='write')

    call header("OpenMC Monte Carlo Code", unit=UNIT_SUMMARY, level=1)
    write(UNIT=UNIT_SUMMARY, FMT=*) &
         "Copyright:     2011-2013 Massachusetts Institute of Technology"
    write(UNIT=UNIT_SUMMARY, FMT='(1X,A,7X,2(I1,"."),I1)') &
         "Version:", VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
#ifdef GIT_SHA1
    write(UNIT=UNIT_SUMMARY, FMT='(1X,"Git SHA1:",6X,A)') GIT_SHA1
#endif
    write(UNIT=UNIT_SUMMARY, FMT='(1X,"Date/Time:",5X,A)') &
         time_stamp()

    ! Write information on number of processors
#ifdef MPI
    write(UNIT=UNIT_SUMMARY, FMT='(1X,"MPI Processes:",1X,A)') &
         trim(to_str(n_procs))
#endif

    ! Display problem summary
    call header("PROBLEM SUMMARY", unit=UNIT_SUMMARY)
    select case(run_mode)
    case (MODE_EIGENVALUE)
      write(UNIT_SUMMARY,100) 'Problem type:', 'k eigenvalue'
      write(UNIT_SUMMARY,101) 'Number of Batches:', n_batches
      write(UNIT_SUMMARY,101) 'Number of Inactive Batches:', n_inactive
      write(UNIT_SUMMARY,101) 'Generations per Batch:', gen_per_batch
    case (MODE_FIXEDSOURCE)
      write(UNIT_SUMMARY,100) 'Problem type:', 'fixed source'
    end select
    write(UNIT_SUMMARY,101) 'Number of Particles:', n_particles

    ! Display geometry summary
    call header("GEOMETRY SUMMARY", unit=UNIT_SUMMARY)
    write(UNIT_SUMMARY,101) 'Number of Cells:', n_cells
    write(UNIT_SUMMARY,101) 'Number of Surfaces:', n_surfaces
    write(UNIT_SUMMARY,101) 'Number of Materials:', n_materials

    ! print summary of all geometry
    call print_geometry()

    ! print summary of materials
    call header("MATERIAL SUMMARY", unit=UNIT_SUMMARY)
    do i = 1, n_materials
      m => materials(i)
      call print_material(m, unit=UNIT_SUMMARY)
    end do

    ! print summary of tallies
    if (n_tallies > 0) then
      call header("TALLY SUMMARY", unit=UNIT_SUMMARY)
      do i = 1, n_tallies
        t=> tallies(i)
        call print_tally(t, unit=UNIT_SUMMARY)
      end do
    end if

    ! print summary of unionized energy grid
    call header("UNIONIZED ENERGY GRID", unit=UNIT_SUMMARY)
    write(UNIT_SUMMARY,*) "Points on energy grid:  " // trim(to_str(n_grid))
    write(UNIT_SUMMARY,*) "Extra storage required: " // trim(to_str(&
         n_grid*n_nuclides_total*4)) // " bytes"

    ! print summary of variance reduction
    call header("VARIANCE REDUCTION", unit=UNIT_SUMMARY)
    if (survival_biasing) then
      write(UNIT_SUMMARY,100) "Survival Biasing:", "on"
    else
      write(UNIT_SUMMARY,100) "Survival Biasing:", "off"
    end if
    write(UNIT_SUMMARY,100) "Weight Cutoff:", trim(to_str(weight_cutoff))
    write(UNIT_SUMMARY,100) "Survival weight:", trim(to_str(weight_survive))

    ! Close summary file
    close(UNIT_SUMMARY)

    ! Format descriptor for columns
100 format (1X,A,T35,A)
101 format (1X,A,T35,I11)

  end subroutine write_summary

!===============================================================================
! WRITE_XS_SUMMARY writes information about each nuclide and S(a,b) table to a
! file called cross_sections.out. This file shows the list of reactions as well
! as information about their secondary angle/energy distributions, how much
! memory is consumed, thresholds, etc.
!===============================================================================

  subroutine write_xs_summary()

    integer                  :: i    ! loop index
    character(MAX_FILE_LEN)  :: path ! path of summary file
    type(Nuclide),    pointer :: nuc => null()
    type(SAlphaBeta), pointer :: sab => null()

    ! Create filename for log file
    path = trim(path_output) // "cross_sections.out"

    ! Open log file for writing
    open(UNIT=UNIT_XS, FILE=path, STATUS='replace', ACTION='write')

    ! Write header
    call header("CROSS SECTION TABLES", unit=UNIT_XS)

    NUCLIDE_LOOP: do i = 1, n_nuclides_total
      ! Get pointer to nuclide
      nuc => nuclides(i)

      ! Print information about nuclide
      call print_nuclide(nuc, unit=UNIT_XS)
    end do NUCLIDE_LOOP

    SAB_TABLES_LOOP: do i = 1, n_sab_tables
      ! Get pointer to S(a,b) table
      sab => sab_tables(i)

      ! Print information about S(a,b) table
      call print_sab_table(sab, unit=UNIT_XS)
    end do SAB_TABLES_LOOP

    ! Close cross section summary file
    close(UNIT_XS)

  end subroutine write_xs_summary

!===============================================================================
! PRINT_COLUMNS displays a header listing what physical values will displayed
! below them
!===============================================================================

  subroutine print_columns()

    write(UNIT=ou, FMT='(2X,A9,3X)', ADVANCE='NO') "Bat./Gen."
    write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "   k    "
    if (entropy_on) write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "Entropy "
    write(UNIT=ou, FMT='(A20,3X)', ADVANCE='NO') "     Average k      "
    if (cmfd_run) then
      write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') " CMFD k "
      select case(trim(cmfd_display))
        case('entropy')
          write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "CMFD Ent"
        case('balance')
          write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "RMS Bal "
        case('source')
          write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "RMS Src "
        case('dominance')
          write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "Dom Rat "
      end select
    end if
    write(UNIT=ou, FMT=*)

    write(UNIT=ou, FMT='(2X,A9,3X)', ADVANCE='NO') "========="
    write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "========"
    if (entropy_on) write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "========"
    write(UNIT=ou, FMT='(A20,3X)', ADVANCE='NO') "===================="
    if (cmfd_run) then
      write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "========"
      if (cmfd_display /= '') &
        write(UNIT=ou, FMT='(A8,3X)', ADVANCE='NO') "========"
    end if
    write(UNIT=ou, FMT=*)

  end subroutine print_columns

!===============================================================================
! PRINT_GENERATION displays information for a generation of neutrons.
!===============================================================================

  subroutine print_generation()

    ! write out information about batch and generation
    write(UNIT=OUTPUT_UNIT, FMT='(2X,A9)', ADVANCE='NO') &
         trim(to_str(current_batch)) // "/" // trim(to_str(current_gen))
    write(UNIT=OUTPUT_UNIT, FMT='(3X,F8.5)', ADVANCE='NO') &
         k_generation(overall_gen)

    ! write out entropy info
    if (entropy_on) write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5)', ADVANCE='NO') &
         entropy(overall_gen)

    if (overall_gen - n_inactive*gen_per_batch > 1) then
      write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5," +/-",F8.5)', ADVANCE='NO') &
           keff, keff_std
    end if

    ! next line
    write(UNIT=OUTPUT_UNIT, FMT=*)

  end subroutine print_generation

!===============================================================================
! PRINT_BATCH_KEFF displays the last batch's tallied value of the neutron
! multiplication factor as well as the average value if we're in active batches
!===============================================================================

  subroutine print_batch_keff()

    ! write out information batch and option independent output
    write(UNIT=OUTPUT_UNIT, FMT='(2X,A9)', ADVANCE='NO') &
         trim(to_str(current_batch)) // "/" // trim(to_str(gen_per_batch))
    write(UNIT=OUTPUT_UNIT, FMT='(3X,F8.5)', ADVANCE='NO') &
         k_generation(overall_gen)

    ! write out entropy info
    if (entropy_on) write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5)', ADVANCE='NO') &
         entropy(current_batch*gen_per_batch)

    ! write out accumulated k-effective if after first active batch
    if (overall_gen - n_inactive*gen_per_batch > 1) then
      write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5," +/-",F8.5)', ADVANCE='NO') &
           keff, keff_std
    else
      write(UNIT=OUTPUT_UNIT, FMT='(23X)', ADVANCE='NO')
    end if

    ! write out cmfd keff if it is active and other display info
    if (cmfd_on) then
      write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5)', ADVANCE='NO') &
         cmfd % k_cmfd(current_batch)
      select case(trim(cmfd_display))
        case('entropy')
          write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5)', ADVANCE='NO') &
            cmfd % entropy(current_batch)
        case('balance')
          write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5)', ADVANCE='NO') &
            cmfd % balance(current_batch)
        case('source')
          write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5)', ADVANCE='NO') &
            cmfd % src_cmp(current_batch)
        case('dominance')
          write(UNIT=OUTPUT_UNIT, FMT='(3X, F8.5)', ADVANCE='NO') &
            cmfd % dom(current_batch)
      end select
    end if

    ! next line
    write(UNIT=OUTPUT_UNIT, FMT=*)

  end subroutine print_batch_keff

!===============================================================================
! PRINT_PLOT displays selected options for plotting
!===============================================================================

  subroutine print_plot()

    integer :: i ! loop index for plots
    type(ObjectPlot), pointer :: pl => null()

    ! Display header for plotting
    call header("PLOTTING SUMMARY")

    do i = 1, n_plots
      pl => plots(i)

      ! Plot id
      write(ou,100) "Plot ID:", trim(to_str(pl % id))

      ! Plot type
      if (pl % type == PLOT_TYPE_SLICE) then
        write(ou,100) "Plot Type:", "Slice"
      else if (pl % type == PLOT_TYPE_VOXEL) then
        write(ou,100) "Plot Type:", "Voxel"
      end if

      ! Plot parameters
      write(ou,100) "Origin:", trim(to_str(pl % origin(1))) // &
           " " // trim(to_str(pl % origin(2))) // " " // &
           trim(to_str(pl % origin(3)))
      if (pl % type == PLOT_TYPE_SLICE) then
        write(ou,100) "Width:", trim(to_str(pl % width(1))) // &
             " " // trim(to_str(pl % width(2)))
      else if (pl % type == PLOT_TYPE_VOXEL) then
        write(ou,100) "Width:", trim(to_str(pl % width(1))) // &
             " " // trim(to_str(pl % width(2))) // &
             " " // trim(to_str(pl % width(3)))
      end if
      if (pl % color_by == PLOT_COLOR_CELLS) then
        write(ou,100) "Coloring:", "Cells"
      else if (pl % color_by == PLOT_COLOR_MATS) then
        write(ou,100) "Coloring:", "Materials"
      end if
      if (pl % type == PLOT_TYPE_SLICE) then
        select case (pl % basis)
        case (PLOT_BASIS_XY)
          write(ou,100) "Basis:", "xy"
        case (PLOT_BASIS_XZ)
          write(ou,100) "Basis:", "xz"
        case (PLOT_BASIS_YZ)
          write(ou,100) "Basis:", "yz"
        end select
        write(ou,100) "Pixels:", trim(to_str(pl % pixels(1))) // " " // &
             trim(to_str(pl % pixels(2)))
      else if (pl % type == PLOT_TYPE_VOXEL) then
        write(ou,100) "Voxels:", trim(to_str(pl % pixels(1))) // " " // &
             trim(to_str(pl % pixels(2))) // " " // trim(to_str(pl % pixels(3)))
      end if

      write(ou,*)

    end do

    ! Format descriptor for columns
100 format (1X,A,T25,A)

  end subroutine print_plot

!===============================================================================
! PRINT_RUNTIME displays the total time elapsed for the entire run, for
! initialization, for computation, and for intergeneration synchronization.
!===============================================================================

  subroutine print_runtime()

    real(8)       :: speed_inactive  ! # of neutrons/second in inactive batches
    real(8)       :: speed_active    ! # of neutrons/second in active batches
    character(15) :: string

    ! display header block
    call header("Timing Statistics")

    ! display time elapsed for various sections
    write(ou,100) "Total time for initialization", time_initialize % elapsed
    write(ou,100) "  Reading cross sections", time_read_xs % elapsed
    write(ou,100) "  Unionizing energy grid", time_unionize % elapsed
    write(ou,100) "Total time in simulation", time_inactive % elapsed + &
         time_active % elapsed
    write(ou,100) "  Time in transport only", time_transport % elapsed
    write(ou,100) "  Time in inactive batches", time_inactive % elapsed
    write(ou,100) "  Time in active batches", time_active % elapsed
    write(ou,100) "  Time synchronizing fission bank", time_bank % elapsed
    write(ou,100) "    Sampling source sites", time_bank_sample % elapsed
    write(ou,100) "    SEND/RECV source sites", time_bank_sendrecv % elapsed
    write(ou,100) "  Time accumulating tallies", time_tallies % elapsed
    if (cmfd_run) write(ou,100) "  Time in CMFD", time_cmfd % elapsed
    if (cmfd_run) write(ou,100) "    Building matrices", &
                  time_cmfdbuild % elapsed
    if (cmfd_run) write(ou,100) "    Solving matrices", &
                  time_cmfdsolve % elapsed
    write(ou,100) "Total time for finalization", time_finalize % elapsed
    write(ou,100) "Total time elapsed", time_total % elapsed

    ! Calculate particle rate in active/inactive batches
    if (restart_run) then
      if (restart_batch < n_inactive) then
        speed_inactive = real(n_particles * (n_inactive - restart_batch) * &
             gen_per_batch) / time_inactive % elapsed
        speed_active = real(n_particles * n_active * gen_per_batch) / &
           time_active % elapsed
      else
        speed_inactive = ZERO
        speed_active = real(n_particles * (n_batches - restart_batch) * &
             gen_per_batch) / time_active % elapsed
      end if
    else
      if (n_inactive > 0) then
        speed_inactive = real(n_particles * n_inactive * gen_per_batch) / &
             time_inactive % elapsed
      end if
      speed_active = real(n_particles * n_active * gen_per_batch) / &
           time_active % elapsed
    end if

    ! display calculation rate
    if (.not. (restart_run .and. (restart_batch >= n_inactive)) &
         .and. n_inactive > 0) then
      string = to_str(speed_inactive)
      write(ou,101) "Calculation Rate (inactive)", trim(string)
    end if
    string = to_str(speed_active)
    write(ou,101) "Calculation Rate (active)", trim(string)

    ! format for write statements
100 format (1X,A,T36,"= ",ES11.4," seconds")
101 format (1X,A,T36,"=  ",A," neutrons/second")

  end subroutine print_runtime

!===============================================================================
! PRINT_RESULTS displays various estimates of k-effective as well as the global
! leakage rate.
!===============================================================================

  subroutine print_results()

    real(8) :: alpha   ! significance level for CI
    real(8) :: t_value ! t-value for confidence intervals

    ! display header block for results
    call header("Results")

    if (confidence_intervals) then
      ! Calculate t-value for confidence intervals
      alpha = ONE - CONFIDENCE_LEVEL
      t_value = t_percentile(ONE - alpha/TWO, n_realizations - 1)

      ! Adjust sum_sq
      global_tallies(:) % sum_sq = t_value * global_tallies(:) % sum_sq

      ! Adjust combined estimator
      if (n_realizations > 3) then
        t_value = t_percentile(ONE - alpha/TWO, n_realizations - 3)
        k_combined(2) = t_value * k_combined(2)
      end if
    end if

    ! write global tallies
    if (n_realizations > 1) then
      write(ou,102) "k-effective (Collision)", global_tallies(K_COLLISION) &
           % sum, global_tallies(K_COLLISION) % sum_sq
      write(ou,102) "k-effective (Track-length)", global_tallies(K_TRACKLENGTH) &
           % sum, global_tallies(K_TRACKLENGTH) % sum_sq
      write(ou,102) "k-effective (Absorption)", global_tallies(K_ABSORPTION) &
           % sum, global_tallies(K_ABSORPTION) % sum_sq
      if (n_realizations > 3) write(ou,102) "Combined k-effective", k_combined
      write(ou,102) "Leakage Fraction", global_tallies(LEAKAGE) % sum, &
           global_tallies(LEAKAGE) % sum_sq
    else
      message = "Could not compute uncertainties -- only one active batch simulated!"
      call warning()

      write(ou,103) "k-effective (Collision)", global_tallies(K_COLLISION) % sum
      write(ou,103) "k-effective (Track-length)", global_tallies(K_TRACKLENGTH)  % sum
      write(ou,103) "k-effective (Absorption)", global_tallies(K_ABSORPTION) % sum
      write(ou,103) "Leakage Fraction", global_tallies(LEAKAGE) % sum
    end if
    write(ou,*)

102 format (1X,A,T30,"= ",F8.5," +/- ",F8.5)
103 format (1X,A,T30,"= ",F8.5)

  end subroutine print_results

!===============================================================================
! PRINT_OVERLAP_DEBUG displays information regarding overlap checking results
!===============================================================================

  subroutine print_overlap_check

    integer :: i, j
    integer :: num_sparse = 0

    ! display header block for geometry debugging section
    call header("Cell Overlap Check Summary")

    write(ou,100) 'Cell ID','No. Overlap Checks'

    do i = 1, n_cells
      write(ou,101) cells(i) % id, overlap_check_cnt(i)
      if (overlap_check_cnt(i) < 10) num_sparse = num_sparse + 1
    end do
    write(ou,*)
    write(ou,'(1X,A)') 'There were ' // trim(to_str(num_sparse)) // &
                       ' cells with less than 10 overlap checks'
    j = 0
    do i = 1, n_cells
      if (overlap_check_cnt(i) < 10) then
        j = j + 1
        write(ou,'(1X,A8)', advance='no') trim(to_str(cells(i) % id))
        if (modulo(j,8) == 0) write(ou,*)
      end if
    end do
    write(ou,*)

100 format (1X,A,T15,A)
101 format (1X,I8,T15,I12)

  end subroutine print_overlap_check

!===============================================================================
! WRITE_TALLIES creates an output file and writes out the mean values of all
! tallies and their standard deviations
!===============================================================================

  subroutine write_tallies()

    integer :: i            ! index in tallies array
    integer :: j            ! level in tally hierarchy
    integer :: k            ! loop index for scoring bins
    integer :: n            ! loop index for nuclides
    integer :: l            ! loop index for user scores
    integer :: type         ! type of tally filter
    integer :: indent       ! number of spaces to preceed output
    integer :: filter_index ! index in results array for filters
    integer :: score_index  ! scoring bin index
    integer :: i_nuclide    ! index in nuclides array
    integer :: i_listing    ! index in xs_listings array
    integer :: n_order      ! loop index for scattering orders
    real(8) :: t_value      ! t-values for confidence intervals
    real(8) :: alpha        ! significance level for CI
    character(MAX_FILE_LEN) :: filename                    ! name of output file
    character(15)           :: filter_name(N_FILTER_TYPES) ! names of tally filters
    character(27)           :: score_names(N_SCORE_TYPES)  ! names of scoring function
    character(27)           :: score_name                  ! names of scoring function
                                                           ! to be applied at write-time
    type(TallyObject), pointer :: t

    ! Skip if there are no tallies
    if (n_tallies == 0) return

    ! Initialize names for tally filter types
    filter_name(FILTER_UNIVERSE)  = "Universe"
    filter_name(FILTER_MATERIAL)  = "Material"
    filter_name(FILTER_CELL)      = "Cell"
    filter_name(FILTER_CELLBORN)  = "Birth Cell"
    filter_name(FILTER_SURFACE)   = "Surface"
    filter_name(FILTER_MESH)      = "Mesh"
    filter_name(FILTER_ENERGYIN)  = "Incoming Energy"
    filter_name(FILTER_ENERGYOUT) = "Outgoing Energy"

    ! Initialize names for scores
    score_names(abs(SCORE_FLUX))          = "Flux"
    score_names(abs(SCORE_TOTAL))         = "Total Reaction Rate"
    score_names(abs(SCORE_SCATTER))       = "Scattering Rate"
    score_names(abs(SCORE_NU_SCATTER))    = "Scattering Production Rate"
    score_names(abs(SCORE_SCATTER_N))     = ""
    score_names(abs(SCORE_SCATTER_PN))    = ""
    score_names(abs(SCORE_TRANSPORT))     = "Transport Rate"
    score_names(abs(SCORE_N_1N))          = "(n,1n) Rate"
    score_names(abs(SCORE_ABSORPTION))    = "Absorption Rate"
    score_names(abs(SCORE_FISSION))       = "Fission Rate"
    score_names(abs(SCORE_NU_FISSION))    = "Nu-Fission Rate"
    score_names(abs(SCORE_KAPPA_FISSION)) = "Kappa-Fission Rate"
    score_names(abs(SCORE_EVENTS))        = "Events"

    ! Create filename for tally output
    filename = trim(path_output) // "tallies.out"

    ! Open tally file for writing
    open(FILE=filename, UNIT=UNIT_TALLY, STATUS='replace', ACTION='write')

    ! Calculate t-value for confidence intervals
    if (confidence_intervals) then
      alpha = ONE - CONFIDENCE_LEVEL
      t_value = t_percentile(ONE - alpha/TWO, n_realizations - 1)
    end if

    TALLY_LOOP: do i = 1, n_tallies
      t => tallies(i)

      if (confidence_intervals) then
        ! Calculate t-value for confidence intervals
        if (confidence_intervals) then
          alpha = ONE - CONFIDENCE_LEVEL
          t_value = t_percentile(ONE - alpha/TWO, t % n_realizations - 1)
        end if

        ! Multiply uncertainty by t-value
        t % results % sum_sq = t_value * t % results % sum_sq
      end if

      ! Write header block
      if (t % label == "") then
        call header("TALLY " // trim(to_str(t % id)), unit=UNIT_TALLY, &
             level=3)
      else
        call header("TALLY " // trim(to_str(t % id)) // ": " &
             // trim(t % label), unit=UNIT_TALLY, level=3)
      endif

      ! Handle surface current tallies separately
      if (t % type == TALLY_SURFACE_CURRENT) then
        call write_surface_current(t)
        cycle
      end if

      ! WARNING: Admittedly, the logic for moving for printing results is
      ! extremely confusing and took quite a bit of time to get correct. The
      ! logic is structured this way since it is not practical to have a do
      ! loop for each filter variable (given that only a few filters are likely
      ! to be used for a given tally.

      ! Initialize bins, filter level, and indentation
      matching_bins(1:t%n_filters) = 0
      j = 1
      indent = 0

      print_bin: do
        find_bin: do
          ! Check for no filters
          if (t % n_filters == 0) exit find_bin

          ! Increment bin combination
          matching_bins(j) = matching_bins(j) + 1

          ! =================================================================
          ! REACHED END OF BINS FOR THIS FILTER, MOVE TO NEXT FILTER

          if (matching_bins(j) > t % filters(j) % n_bins) then
            ! If this is the first filter, then exit
            if (j == 1) exit print_bin

            matching_bins(j) = 0
            j = j - 1
            indent = indent - 2

            ! =================================================================
            ! VALID BIN -- WRITE FILTER INFORMATION OR EXIT TO WRITE RESULTS

          else
            ! Check if this is last filter
            if (j == t % n_filters) exit find_bin

            ! Print current filter information
            type = t % filters(j) % type
            write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                 trim(filter_name(type)), trim(get_label(t, j))
            indent = indent + 2
            j = j + 1
          end if

        end do find_bin

        ! Print filter information
        if (t % n_filters > 0) then
          type = t % filters(j) % type
          write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
               trim(filter_name(type)), trim(get_label(t, j))
        end if

        ! Determine scoring index for this bin combination -- note that unlike
        ! in the score_tally subroutine, we have to use max(bins,1) since all
        ! bins below the lowest filter level will be zeros

        if (t % n_filters > 0) then
          filter_index = sum((max(matching_bins(1:t%n_filters),1) - 1) * t % stride) + 1
        else
          filter_index = 1
        end if

        ! Write results for this filter bin combination
        score_index = 0
        if (t % n_filters > 0) indent = indent + 2
        do n = 1, t % n_nuclide_bins
          ! Write label for nuclide
          i_nuclide = t % nuclide_bins(n)
          if (i_nuclide == -1) then
            write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                 "Total Material"
          else
            i_listing = nuclides(i_nuclide) % listing
            write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                 trim(xs_listings(i_listing) % alias)
          end if

          indent = indent + 2
          k = 0
          do l = 1, t % n_user_score_bins
            k = k + 1
            score_index = score_index + 1
            if (t % score_bins(k) == SCORE_SCATTER_N) then
              if (t % scatt_order(k) == 0) then
                score_name = "Scattering Rate"
              else
                score_name = 'P' // trim(to_str(t % scatt_order(k))) // &
                  ' Scattering Moment'
              end if
              write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') &
                repeat(" ", indent), score_name, &
                to_str(t % results(score_index,filter_index) % sum), &
                trim(to_str(t % results(score_index,filter_index) % sum_sq))
            else if (t % score_bins(k) == SCORE_SCATTER_PN) then
              score_name = "Scattering Rate"
              write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') &
                repeat(" ", indent), score_name, &
                to_str(t % results(score_index,filter_index) % sum), &
                trim(to_str(t % results(score_index,filter_index) % sum_sq))
              do n_order = 1, t % scatt_order(k)
                score_index = score_index + 1
                score_name = 'P' // trim(to_str(n_order)) // &
                  ' Scattering Moment'
                write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') &
                  repeat(" ", indent), score_name, &
                  to_str(t % results(score_index,filter_index) % sum), &
                  trim(to_str(t % results(score_index,filter_index) % sum_sq))
              end do
              k = k + n_order - 1
            else
              if (t % score_bins(k) > 0) then
                score_name = reaction_name(t % score_bins(k))
              else
                score_name = score_names(abs(t % score_bins(k)))
              end if
              write(UNIT=UNIT_TALLY, FMT='(1X,2A,1X,A,"+/- ",A)') &
                repeat(" ", indent), score_name, &
                to_str(t % results(score_index,filter_index) % sum), &
                trim(to_str(t % results(score_index,filter_index) % sum_sq))
            end if
          end do
          indent = indent - 2

        end do
        indent = indent - 2

        if (t % n_filters == 0) exit print_bin

      end do print_bin

    end do TALLY_LOOP

    close(UNIT=UNIT_TALLY)

  end subroutine write_tallies

!===============================================================================
! WRITE_SURFACE_CURRENT writes out surface current tallies over a mesh to the
! tallies.out file.
!===============================================================================

  subroutine write_surface_current(t)

    type(TallyObject), pointer :: t

    integer :: i                    ! mesh index for x
    integer :: j                    ! mesh index for y
    integer :: k                    ! mesh index for z
    integer :: l                    ! index for energy
    integer :: i_filter_mesh        ! index for mesh filter
    integer :: i_filter_ein         ! index for incoming energy filter
    integer :: i_filter_surf        ! index for surface filter
    integer :: n                    ! number of incoming energy bins
    integer :: len1                 ! length of string
    integer :: len2                 ! length of string
    integer :: filter_index         ! index in results array for filters
    logical :: print_ebin           ! should incoming energy bin be displayed?
    character(MAX_LINE_LEN) :: string
    type(StructuredMesh), pointer :: m => null()

    ! Get pointer to mesh
    i_filter_mesh = t % find_filter(FILTER_MESH)
    i_filter_surf = t % find_filter(FILTER_SURFACE)
    m => meshes(t % filters(i_filter_mesh) % int_bins(1))

    ! initialize bins array
    matching_bins(1:t%n_filters) = 1

    ! determine how many energy in bins there are
    i_filter_ein = t % find_filter(FILTER_ENERGYIN)
    if (i_filter_ein > 0) then
      print_ebin = .true.
      n = t % filters(i_filter_ein) % n_bins
    else
      print_ebin = .false.
      n = 1
    end if

    do i = 1, m % dimension(1)
      string = "Mesh Index (" // trim(to_str(i)) // ", "
      len1 = len_trim(string)
      do j = 1, m % dimension(2)
        string = string(1:len1+1) // trim(to_str(j)) // ", "
        len2 = len_trim(string)
        do k = 1, m % dimension(3)
          ! Write mesh cell index
          string = string(1:len2+1) // trim(to_str(k)) // ")"
          write(UNIT=UNIT_TALLY, FMT='(1X,A)') trim(string)

          do l = 1, n
            if (print_ebin) then
              ! Set incoming energy bin
              matching_bins(i_filter_ein) = l

              ! Write incoming energy bin
              write(UNIT=UNIT_TALLY, FMT='(3X,A,1X,A)') &
                   "Incoming Energy", trim(get_label(t, i_filter_ein))
            end if

            ! Left Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i-1, j, k /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_RIGHT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Left", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_RIGHT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Left", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            ! Right Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i, j, k /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_RIGHT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Right", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_RIGHT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Right", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            ! Back Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i, j-1, k /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_FRONT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Back", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_FRONT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Back", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            ! Front Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i, j, k /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_FRONT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Front", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_FRONT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Front", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            ! Bottom Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i, j, k-1 /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_TOP
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Bottom", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_TOP
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Bottom", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            ! Top Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i, j, k /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_TOP
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Top", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_TOP
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=UNIT_TALLY, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Top", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))
          end do

        end do
      end do
    end do

  end subroutine write_surface_current

!===============================================================================
! GET_LABEL returns a label for a cell/surface/etc given a tally, filter type,
! and corresponding bin
!===============================================================================

  function get_label(t, i_filter) result(label)

    type(TallyObject), pointer :: t        ! tally object
    integer, intent(in)        :: i_filter ! index in filters array
    character(30)              :: label    ! user-specified identifier

    integer :: i      ! index in cells/surfaces/etc array
    integer :: bin
    integer, allocatable :: ijk(:) ! indices in mesh
    real(8)              :: E0     ! lower bound for energy bin
    real(8)              :: E1     ! upper bound for energy bin
    type(StructuredMesh), pointer :: m => null()

    bin = matching_bins(i_filter)

    select case(t % filters(i_filter) % type)
    case (FILTER_UNIVERSE)
      i = t % filters(i_filter) % int_bins(bin)
      label = to_str(universes(i) % id)
    case (FILTER_MATERIAL)
      i = t % filters(i_filter) % int_bins(bin)
      label = to_str(materials(i) % id)
    case (FILTER_CELL, FILTER_CELLBORN)
      i = t % filters(i_filter) % int_bins(bin)
      label = to_str(cells(i) % id)
    case (FILTER_SURFACE)
      i = t % filters(i_filter) % int_bins(bin)
      label = to_str(surfaces(i) % id)
    case (FILTER_MESH)
      m => meshes(t % filters(i_filter) % int_bins(1))
      allocate(ijk(m % n_dimension))
      call bin_to_mesh_indices(m, bin, ijk)
      if (m % n_dimension == 2) then
        label = "Index (" // trim(to_str(ijk(1))) // ", " // &
             trim(to_str(ijk(2))) // ")"
      elseif (m % n_dimension == 3) then
        label = "Index (" // trim(to_str(ijk(1))) // ", " // &
             trim(to_str(ijk(2))) // ", " // trim(to_str(ijk(3))) // ")"
      end if
    case (FILTER_ENERGYIN, FILTER_ENERGYOUT)
      E0 = t % filters(i_filter) % real_bins(bin)
      E1 = t % filters(i_filter) % real_bins(bin + 1)
      label = "[" // trim(to_str(E0)) // ", " // trim(to_str(E1)) // ")"
    end select

  end function get_label

end module output
