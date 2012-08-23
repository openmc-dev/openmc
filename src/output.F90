module output

  use, intrinsic :: ISO_FORTRAN_ENV

  use ace_header,      only: Nuclide, Reaction, UrrData
  use constants
  use datatypes,       only: dict_get_key
  use endf,            only: reaction_name
  use geometry_header, only: Cell, Universe, Surface, BASE_UNIVERSE
  use global
  use math,            only: t_percentile
  use mesh_header,     only: StructuredMesh
  use particle_header, only: LocalCoord
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

    character(10) :: date_
    character(10) :: time_

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
         '     Developed At:  Massachusetts Institute of Technology'
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Version:",7X,I1,".",I1,".",I1)') &
         VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
#ifdef GIT_SHA1
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Git SHA1:",6X,A)') GIT_SHA1
#endif

    ! Write the date and time
    call date_and_time(DATE=date_, TIME=time_)
    date_ = date_(1:4) // "-" // date_(5:6) // "-" // date_(7:8)
    time_ = time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)
    write(UNIT=OUTPUT_UNIT, FMT='(6X,"Date/Time:",5X,A,1X,A)') &
         trim(date_), trim(time_)

    ! Write information to summary file
    if (output_summary) then
       call header("OpenMC Monte Carlo Code", unit=UNIT_SUMMARY, level=1)
       write(UNIT=UNIT_SUMMARY, FMT=*) &
            "Copyright:     2011-2012 Massachusetts Institute of Technology"
       write(UNIT=UNIT_SUMMARY, FMT='(1X,A,7X,2(I1,"."),I1)') &
            "Version:", VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
#ifdef GIT_SHA1
       write(UNIT=UNIT_SUMMARY, FMT='(1X,"Git SHA1:",6X,A)') GIT_SHA1
#endif
       write(UNIT=UNIT_SUMMARY, FMT='(1X,"Date/Time:",5X,A,1X,A)') &
            trim(date_), trim(time_)

       ! Write information on number of processors
#ifdef MPI
       write(UNIT=OUTPUT_UNIT, FMT='(1X,A)') '     MPI Processes: ' // &
            trim(to_str(n_procs))
       write(UNIT=UNIT_SUMMARY, FMT='(1X,"MPI Processes:",1X,A)') &
            trim(to_str(n_procs))
#endif
    end if

  end subroutine title

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
       write(UNIT=OUTPUT_UNIT, FMT=*) "Copyright (c) 2011-2012 &
            &Massachusetts Institute of Technology"
       write(UNIT=OUTPUT_UNIT, FMT=*) "MIT/X license at &
            &<http://mit-crpg.github.com/openmc/license.html>"
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
       write(OUTPUT_UNIT,*) '  -p, --plot      Run in plotting mode'
       write(OUTPUT_UNIT,*) '  -r, --restart   Restart a previous run'
       write(OUTPUT_UNIT,*) '  -t, --tallies   Write tally results from state point'
       write(OUTPUT_UNIT,*) '  -v, --version   Show version information'
       write(OUTPUT_UNIT,*) '  -?, --help      Show this message'
    end if

  end subroutine print_usage

!===============================================================================
! WRITE_MESSAGE displays an informational message to the log file and the 
! standard output stream.
!===============================================================================

  subroutine write_message(level)

    integer, optional :: level ! verbosity level

    integer :: n_lines ! number of lines needed
    integer :: i       ! index for lines

    ! Only allow master to print to screen
    if (.not. master .and. present(level)) return

    ! TODO: Take care of line wrapping so words don't get cut off
    if (.not. present(level) .or. level <= verbosity) then
       n_lines = (len_trim(message)-1)/79 + 1
       do i = 1, n_lines
          write(ou, fmt='(1X,A)') trim(message(79*(i-1)+1:79*i))
       end do
    end if

  end subroutine write_message

!===============================================================================
! PRINT_PARTICLE displays the attributes of a particle
!===============================================================================

  subroutine print_particle()

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
    write(ou,*) '  Energy grid index = ' // to_str(p % index_grid)
    write(ou,*) '  Interpolation factor = ' // to_str(p % interp)
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
    index_cell = dict_get_key(cell_dict, c % id)
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
          string = trim(string) // ' ' // to_str(c % surfaces(i))
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

    integer :: unit_ ! unit to write to

    ! set default unit if not specified
    if (present(unit)) then
       unit_ = unit
    else
       unit_ = OUTPUT_UNIT
    end if

    ! Write information about lattice
    write(unit_,*) 'Lattice ' // to_str(lat % id)
    write(unit_,*) '    n_x = ' // to_str(lat % n_x)
    write(unit_,*) '    n_y = ' // to_str(lat % n_y)
    write(unit_,*) '    x0 = ' // to_str(lat % x0)
    write(unit_,*) '    y0 = ' // to_str(lat % y0)
    write(unit_,*) '    width_x = ' // to_str(lat % width_x)
    write(unit_,*) '    width_y = ' // to_str(lat % width_y)
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
    case (SURF_BOX_X)
    case (SURF_BOX_Y)
    case (SURF_BOX_Z)
    case (SURF_BOX)
    case (SURF_GQ)
       string = "General Quadratic"
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
          string = trim(string) // ' ' // to_str(surf % neighbor_pos(i))
       end do
    end if
    write(unit_,*) '    Positive Neighbors = ' // trim(string)

    ! Write neighboring cells on negative side of this surface
    string = ""
    if (allocated(surf % neighbor_neg)) then
       do i = 1, size(surf % neighbor_neg)
          string = trim(string) // ' ' // to_str(surf % neighbor_neg(i))
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
    if (mat % has_sab_table) then
       write(unit_,*) '    S(a,b) table = ' // trim(mat % sab_name)
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
    integer :: id    ! user-specified id
    integer :: unit_ ! unit to write to
    character(MAX_LINE_LEN) :: string
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
    if (t % n_filter_bins(FILTER_CELL) > 0) then
       string = ""
       do i = 1, t % n_filter_bins(FILTER_CELL)
          id = t % cell_bins(i)
          c => cells(id)
          string = trim(string) // ' ' // trim(to_str(c % id))
       end do
       write(unit_, *) '    Cell Bins:' // trim(string)
    end if

    ! Write any surface bins if present
    if (t % n_filter_bins(FILTER_SURFACE) > 0) then
       string = ""
       do i = 1, t % n_filter_bins(FILTER_SURFACE)
          id = t % surface_bins(i)
          s => surfaces(id)
          string = trim(string) // ' ' // trim(to_str(s % id))
       end do
       write(unit_, *) '    Surface Bins:' // trim(string)
    end if

    ! Write any universe bins if present
    if (t % n_filter_bins(FILTER_UNIVERSE) > 0) then
       string = ""
       do i = 1, t % n_filter_bins(FILTER_UNIVERSE)
          id = t % universe_bins(i)
          u => universes(id)
          string = trim(string) // ' ' // trim(to_str(u % id))
       end do
       write(unit_, *) '    Universe Bins:' // trim(string)
    end if

    ! Write any material bins if present
    if (t % n_filter_bins(FILTER_MATERIAL) > 0) then
       string = ""
       do i = 1, t % n_filter_bins(FILTER_MATERIAL)
          id = t % material_bins(i)
          m => materials(id)
          string = trim(string) // ' ' // trim(to_str(m % id))
       end do
       write(unit_, *) '    Material Bins:' // trim(string)
    end if

    ! Write any mesh bins if present
    if (t % n_filter_bins(FILTER_MESH) > 0) then
       string = ""
       id = t % mesh
       sm => meshes(id)
       string = trim(string) // ' ' // trim(to_str(sm % dimension(1)))
       do i = 2, sm % n_dimension
          string = trim(string) // ' x ' // trim(to_str(sm % dimension(i)))
       end do
       write(unit_, *) '    Mesh Bins:' // trim(string)
    end if

    ! Write any birth region bins if present
    if (t % n_filter_bins(FILTER_CELLBORN) > 0) then
       string = ""
       do i = 1, t % n_filter_bins(FILTER_CELLBORN)
          id = t % cellborn_bins(i)
          c => cells(id)
          string = trim(string) // ' ' // trim(to_str(c % id))
       end do
       write(unit_, *) '    Birth Region Bins:' // trim(string)
    end if

    ! Write any incoming energy bins if present
    if (t % n_filter_bins(FILTER_ENERGYIN) > 0) then
       string = ""
       do i = 1, t % n_filter_bins(FILTER_ENERGYIN) + 1
          string = trim(string) // ' ' // trim(to_str(&
               t % energy_in(i)))
       end do
       write(unit_,*) '    Incoming Energy Bins:' // trim(string)
    end if

    ! Write any outgoing energy bins if present
    if (t % n_filter_bins(FILTER_ENERGYOUT) > 0) then
       string = ""
       do i = 1, t % n_filter_bins(FILTER_ENERGYOUT) + 1
          string = trim(string) // ' ' // trim(to_str(&
               t % energy_out(i)))
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
    string = ""
    do i = 1, t % n_score_bins
       select case (t % score_bins(i))
       case (SCORE_FLUX)
          string = trim(string) // ' flux'
       case (SCORE_TOTAL)
          string = trim(string) // ' total'
       case (SCORE_SCATTER)
          string = trim(string) // ' scatter'
       case (SCORE_NU_SCATTER)
          string = trim(string) // ' nu-scatter'
       case (SCORE_SCATTER_1)
          string = trim(string) // ' scatter-1'
       case (SCORE_SCATTER_2)
          string = trim(string) // ' scatter-2'
       case (SCORE_SCATTER_3)
          string = trim(string) // ' scatter-3'
       case (SCORE_TRANSPORT)
          string = trim(string) // ' transport'
       case (SCORE_DIFFUSION)
          string = trim(string) // ' diffusion'
       case (SCORE_N_1N)
          string = trim(string) // ' n1n'
       case (SCORE_N_2N)
          string = trim(string) // ' n2n'
       case (SCORE_N_3N)
          string = trim(string) // ' n3n'
       case (SCORE_N_4N)
          string = trim(string) // ' n4n'
       case (SCORE_ABSORPTION)
          string = trim(string) // ' absorption'
       case (SCORE_FISSION)
          string = trim(string) // ' fission'
       case (SCORE_NU_FISSION)
          string = trim(string) // ' nu-fission'
       case (SCORE_CURRENT)
          string = trim(string) // ' current'
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

    type(SAB_Table), pointer :: sab
    integer,        optional :: unit

    integer :: size_sab ! memory used by S(a,b) table
    integer :: unit_    ! unit to write to

    ! set default unit for writing information
    if (present(unit)) then
       unit_ = unit
    else
       unit_ = OUTPUT_UNIT
    end if

    ! Basic S(a,b) table information
    write(unit_,*) 'S(a,b) Table ' // trim(sab % name)
    write(unit_,*) '  zaid = ' // trim(to_str(sab % zaid))
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
! PRINT_SUMMARY displays summary information about the problem about to be run
! after reading all input files
!===============================================================================

  subroutine print_summary()

    integer :: i ! loop index
    character(15) :: string
    type(Material),    pointer :: m => null()
    type(TallyObject), pointer :: t => null()

    ! Display problem summary
    call header("PROBLEM SUMMARY", unit=UNIT_SUMMARY)
    select case(run_mode)
    case (MODE_CRITICALITY)
       write(UNIT_SUMMARY,100) 'Problem type:', 'Criticality'
       write(UNIT_SUMMARY,101) 'Number of Batches:', n_batches
       write(UNIT_SUMMARY,101) 'Number of Inactive Batches:', n_inactive
       write(UNIT_SUMMARY,101) 'Generations per Batch:', gen_per_batch
    case (MODE_FIXEDSOURCE)
       write(UNIT_SUMMARY,100) 'Problem type:', 'External Source'
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
    string = to_str(weight_cutoff)
    write(UNIT_SUMMARY,100) "Weight Cutoff:", trim(string)
    string = to_str(weight_survive)
    write(UNIT_SUMMARY,100) "Survival weight:", trim(string)

    ! Format descriptor for columns
100 format (1X,A,T35,A)
101 format (1X,A,T35,I11)

  end subroutine print_summary

!===============================================================================
! PRINT_COLUMNS displays a header listing what physical values will displayed
! below them
!===============================================================================

  subroutine print_columns()

    if (entropy_on) then
       message = " Batch   k(batch)   Entropy         Average k"
       call write_message(1)
       message = " =====   ========   =======    ==================="
       call write_message(1)
    else
       message = " Batch   k(batch)          Average k"
       call write_message(1)
       message = " =====   ========     ==================="
       call write_message(1)
    end if

  end subroutine print_columns

!===============================================================================
! PRINT_BATCH_KEFF displays the last batch's tallied value of the neutron
! multiplication factor as well as the average value if we're in active batches
!===============================================================================

  subroutine print_batch_keff()

    if (current_batch <= n_inactive) then
       ! ======================================================================
       ! INACTIVE BATCHES

       if (entropy_on) then
          write(UNIT=OUTPUT_UNIT, FMT=102) current_batch, &
               k_batch(current_batch), entropy(current_batch)
       else
          write(UNIT=OUTPUT_UNIT, FMT=100) current_batch, &
               k_batch(current_batch)
       end if

    elseif (current_batch == n_inactive + 1) then
       ! ======================================================================
       ! ACTIVE BATCHES

       if (entropy_on) then
          write(UNIT=OUTPUT_UNIT, FMT=102) current_batch, &
               k_batch(current_batch), entropy(current_batch)
       else
          write(UNIT=OUTPUT_UNIT, FMT=100) current_batch, &
               k_batch(current_batch)
       end if

    elseif (current_batch > n_inactive + 1) then

       if (entropy_on) then
          write(UNIT=OUTPUT_UNIT, FMT=103) current_batch, &
               k_batch(current_batch), entropy(current_batch), keff, keff_std
       else
          write(UNIT=OUTPUT_UNIT, FMT=101) current_batch, &
               k_batch(current_batch), keff, keff_std
       end if

    end if

100 format (2X,I5,2X,F8.5)
101 format (2X,I5,2X,F8.5,5X,F8.5," +/-",F8.5)
102 format (2X,I5,2X,F8.5,3X,F8.5)
103 format (2X,I5,2X,F8.5,3X,F8.5,3X,F8.5," +/-",F8.5)

  end subroutine print_batch_keff

!===============================================================================
! PRINT_PLOT displays selected options for plotting
!===============================================================================

  subroutine print_plot()

    integer :: i ! loop index for plots
    type(Plot), pointer :: pl => null()

    ! Display header for plotting
    call header("PLOTTING SUMMARY")

    do i = 1, n_plots
      pl => plots(i)

      ! Write plot id
      write(ou,100) "Plot ID:", trim(to_str(pl % id))

      ! Write plotting origin
      write(ou,100) "Origin:", trim(to_str(pl % origin(1))) // &
           " " // trim(to_str(pl % origin(2))) // " " // &
           trim(to_str(pl % origin(3)))

      ! Write plotting width
      if (pl % type == PLOT_TYPE_SLICE) then

        write(ou,100) "Width:", trim(to_str(pl % width(1))) // &
             " " // trim(to_str(pl % width(2)))
        write(ou,100) "Coloring:", trim(to_str(pl % color_by))
        write(ou,100) "Basis:", trim(to_str(pl % basis))
        write(ou,100) "Pixels:", trim(to_str(pl % pixels(1))) // " " // &
                                 trim(to_str(pl % pixels(2)))
      end if

      write(ou,*)

    end do

    ! Format descriptor for columns
100 format (1X,A,T25,A)

  end subroutine print_plot

!===============================================================================
! PRINT_RUNTIME displays the total time elapsed for the entire run, for
! initialization, for computation, and for intercycle synchronization.
!===============================================================================

  subroutine print_runtime()

    integer(8)    :: total_particles ! total # of particles simulated
    real(8)       :: speed           ! # of neutrons/second
    real(8)       :: alpha           ! significance level for CI
    real(8)       :: t_value         ! t-value for confidence intervals
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
    write(ou,100) "  Time between generations", time_intercycle % elapsed
    write(ou,100) "    Accumulating tallies", time_ic_tallies % elapsed
    write(ou,100) "    Sampling source sites", time_ic_sample % elapsed
    write(ou,100) "    SEND/RECV source sites", time_ic_sendrecv % elapsed
    write(ou,100) "Total time for finalization", time_finalize % elapsed
    write(ou,100) "Total time elapsed", time_total % elapsed

    if (restart_run) then
       total_particles = n_particles * (n_batches - &
            restart_batch) * gen_per_batch
    else
       total_particles = n_particles * n_batches * gen_per_batch
    end if

    ! display calculate rate
    speed = real(total_particles) / (time_inactive % elapsed + &
         time_active % elapsed)
    string = to_str(speed)
    write(ou,101) "Calculation Rate", trim(string)

    ! display header block for results
    call header("Results")

    if (confidence_intervals) then
       ! Calculate t-value for confidence intervals
       alpha = ONE - CONFIDENCE_LEVEL
       t_value = t_percentile(ONE - alpha/TWO, n_realizations - 1)

       ! Adjust sum_sq
       global_tallies(:) % sum_sq = t_value * global_tallies(:) % sum_sq
    end if
    
    ! write global tallies
    write(ou,102) "k-effective (Analog)", global_tallies(K_ANALOG) % sum, &
         global_tallies(K_ANALOG) % sum_sq
    write(ou,102) "k-effective (Collision)", global_tallies(K_COLLISION) % sum, &
         global_tallies(K_COLLISION) % sum_sq
    write(ou,102) "k-effective (Track-length)", global_tallies(K_TRACKLENGTH) % sum, &
         global_tallies(K_TRACKLENGTH) % sum_sq
    write(ou,102) "Leakage Fraction", global_tallies(LEAKAGE) % sum, &
         global_tallies(LEAKAGE) % sum_sq
    write(ou,*)

    ! format for write statements
100 format (1X,A,T35,"= ",ES11.4," seconds")
101 format (1X,A,T35,"=  ",A," neutrons/second")
102 format (1X,A,T30,"= ",F8.5," +/- ",F8.5)
 
  end subroutine print_runtime

!===============================================================================
! CREATE_SUMMARY_FILE opens the summary.out file for logging information about
! the simulation
!===============================================================================

  subroutine create_summary_file()

    logical :: file_exists  ! does log file already exist?
    character(MAX_FILE_LEN) :: path ! path of summary file

    ! Create filename for log file
    path = "summary.out"

    ! Check if log file already exists
    inquire(FILE=path, EXIST=file_exists)
    if (file_exists) then
       ! Possibly copy old log file
    end if

    ! Open log file for writing
    open(UNIT=UNIT_SUMMARY, FILE=path, STATUS='replace', &
         ACTION='write')

  end subroutine create_summary_file

!===============================================================================
! CREATE_XS_SUMMARY_FILE creates an output file to write information about the
! cross section tables used in the simulation.
!===============================================================================

  subroutine create_xs_summary_file()

    logical :: file_exists  ! does log file already exist?
    character(MAX_FILE_LEN) :: path ! path of summary file

    ! Create filename for log file
    path = "cross_sections.out"

    ! Check if log file already exists
    inquire(FILE=path, EXIST=file_exists)
    if (file_exists) then
       ! Possibly copy old log file
    end if

    ! Open log file for writing
    open(UNIT=UNIT_XS, FILE=path, STATUS='replace', &
         ACTION='write')

  end subroutine create_xs_summary_file

end module output
