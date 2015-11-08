module output

  use, intrinsic :: ISO_FORTRAN_ENV

  use ace_header,      only: Nuclide, Reaction, UrrData
  use constants
  use endf,            only: reaction_name
  use error,           only: fatal_error, warning
  use geometry_header, only: Cell, Universe, Lattice, RectLattice, &
                             HexLattice, BASE_UNIVERSE
  use global
  use math,            only: t_percentile
  use mesh_header,     only: RegularMesh
  use mesh,            only: mesh_indices_to_bin, bin_to_mesh_indices
  use particle_header, only: LocalCoord, Particle
  use plot_header
  use string,          only: to_upper, to_str
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
         '     Copyright:      2011-2015 Massachusetts Institute of Technology'
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
    line = to_upper(msg)

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
      write(UNIT=OUTPUT_UNIT, FMT=*) "Copyright (c) 2011-2015 &
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
      write(OUTPUT_UNIT,*) '  -h, --help             Show this message'
    end if

  end subroutine print_usage

!===============================================================================
! WRITE_MESSAGE displays an informational message to the log file and the
! standard output stream.
!===============================================================================

  subroutine write_message(message, level)

    character(*) :: message
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
        if (length - i_start < line_wrap + 1) then
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
    type(Cell),       pointer :: c
    type(Universe),   pointer :: u
    class(Lattice),   pointer :: l

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
    do i = 1, p % n_coord
      ! Print level
      write(ou,*) '  Level ' // trim(to_str(i - 1))

      ! Print cell for this level
      if (p % coord(i) % cell /= NONE) then
        c => cells(p % coord(i) % cell)
        write(ou,*) '    Cell             = ' // trim(to_str(c % id))
      end if

      ! Print universe for this level
      if (p % coord(i) % universe /= NONE) then
        u => universes(p % coord(i) % universe)
        write(ou,*) '    Universe         = ' // trim(to_str(u % id))
      end if

      ! Print information on lattice
      if (p % coord(i) % lattice /= NONE) then
        l => lattices(p % coord(i) % lattice) % obj
        write(ou,*) '    Lattice          = ' // trim(to_str(l % id))
        write(ou,*) '    Lattice position = (' // trim(to_str(&
             p % coord(i) % lattice_x)) // ',' // trim(to_str(&
             p % coord(i) % lattice_y)) // ')'
      end if

      ! Print local coordinates
      write(ou,'(1X,A,3ES12.4)') '    xyz = ', p % coord(i) % xyz
      write(ou,'(1X,A,3ES12.4)') '    uvw = ', p % coord(i) % uvw
    end do

    ! Print surface
    if (p % surface /= NONE) then
      write(ou,*) '  Surface = ' // to_str(sign(surfaces(i)%obj%id, p % surface))
    end if

    ! Display weight, energy, grid index, and interpolation factor
    write(ou,*) '  Weight = ' // to_str(p % wgt)
    write(ou,*) '  Energy = ' // to_str(p % E)
    write(ou,*) '  Delayed Group = ' // to_str(p % delayed_group)
    write(ou,*)

  end subroutine print_particle

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
! WRITE_XS_SUMMARY writes information about each nuclide and S(a,b) table to a
! file called cross_sections.out. This file shows the list of reactions as well
! as information about their secondary angle/energy distributions, how much
! memory is consumed, thresholds, etc.
!===============================================================================

  subroutine write_xs_summary()

    integer :: i       ! loop index
    integer :: unit_xs ! cross_sections.out file unit
    character(MAX_FILE_LEN)  :: path ! path of summary file
    type(Nuclide),    pointer :: nuc => null()
    type(SAlphaBeta), pointer :: sab => null()

    ! Create filename for log file
    path = trim(path_output) // "cross_sections.out"

    ! Open log file for writing
    open(NEWUNIT=unit_xs, FILE=path, STATUS='replace', ACTION='write')

    ! Write header
    call header("CROSS SECTION TABLES", unit=unit_xs)

    NUCLIDE_LOOP: do i = 1, n_nuclides_total
      ! Get pointer to nuclide
      nuc => nuclides(i)

      ! Print information about nuclide
      call print_nuclide(nuc, unit=unit_xs)
    end do NUCLIDE_LOOP

    SAB_TABLES_LOOP: do i = 1, n_sab_tables
      ! Get pointer to S(a,b) table
      sab => sab_tables(i)

      ! Print information about S(a,b) table
      call print_sab_table(sab, unit=unit_xs)
    end do SAB_TABLES_LOOP

    ! Close cross section summary file
    close(unit_xs)

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

      ! Plot filename
      write(ou,100) "Plot file:", trim(pl % path_plot)

      ! Plot level
      write(ou,100) "Universe depth:", trim(to_str(pl % level))

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
    write(ou,100) "Total time in simulation", time_inactive % elapsed + &
         time_active % elapsed
    write(ou,100) "  Time in transport only", time_transport % elapsed
    if (run_mode == MODE_EIGENVALUE) then
      write(ou,100) "  Time in inactive batches", time_inactive % elapsed
    end if
    write(ou,100) "  Time in active batches", time_active % elapsed
    if (run_mode == MODE_EIGENVALUE) then
      write(ou,100) "  Time synchronizing fission bank", time_bank % elapsed
      write(ou,100) "    Sampling source sites", time_bank_sample % elapsed
      write(ou,100) "    SEND/RECV source sites", time_bank_sendrecv % elapsed
    end if
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
      if (run_mode == MODE_EIGENVALUE) then
        write(ou,102) "k-effective (Collision)", global_tallies(K_COLLISION) &
             % sum, global_tallies(K_COLLISION) % sum_sq
        write(ou,102) "k-effective (Track-length)", global_tallies(K_TRACKLENGTH) &
             % sum, global_tallies(K_TRACKLENGTH) % sum_sq
        write(ou,102) "k-effective (Absorption)", global_tallies(K_ABSORPTION) &
             % sum, global_tallies(K_ABSORPTION) % sum_sq
        if (n_realizations > 3) write(ou,102) "Combined k-effective", k_combined
      end if
      write(ou,102) "Leakage Fraction", global_tallies(LEAKAGE) % sum, &
           global_tallies(LEAKAGE) % sum_sq
    else
      if (master) call warning("Could not compute uncertainties -- only one &
           &active batch simulated!")

      if (run_mode == MODE_EIGENVALUE) then
        write(ou,103) "k-effective (Collision)", global_tallies(K_COLLISION) % sum
        write(ou,103) "k-effective (Track-length)", global_tallies(K_TRACKLENGTH)  % sum
        write(ou,103) "k-effective (Absorption)", global_tallies(K_ABSORPTION) % sum
      end if
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
    integer :: n_order      ! loop index for moment orders
    integer :: nm_order     ! loop index for Ynm moment orders
    integer :: unit_tally   ! tallies.out file unit
    real(8) :: t_value      ! t-values for confidence intervals
    real(8) :: alpha        ! significance level for CI
    character(MAX_FILE_LEN) :: filename                    ! name of output file
    character(16)           :: filter_name(N_FILTER_TYPES) ! names of tally filters
    character(36)           :: score_names(N_SCORE_TYPES)  ! names of scoring function
    character(36)           :: score_name                  ! names of scoring function
                                                           ! to be applied at write-time
    type(TallyObject), pointer :: t

    ! Skip if there are no tallies
    if (n_tallies == 0) return

    ! Initialize names for tally filter types
    filter_name(FILTER_UNIVERSE)     = "Universe"
    filter_name(FILTER_MATERIAL)     = "Material"
    filter_name(FILTER_DISTRIBCELL)  = "Distributed Cell"
    filter_name(FILTER_CELL)         = "Cell"
    filter_name(FILTER_CELLBORN)     = "Birth Cell"
    filter_name(FILTER_SURFACE)      = "Surface"
    filter_name(FILTER_MESH)         = "Mesh"
    filter_name(FILTER_ENERGYIN)     = "Incoming Energy"
    filter_name(FILTER_ENERGYOUT)    = "Outgoing Energy"
    filter_name(FILTER_MU)           = "Change-in-Angle"
    filter_name(FILTER_POLAR)        = "Polar Angle"
    filter_name(FILTER_AZIMUTHAL)    = "Azimuthal Angle"
    filter_name(FILTER_DELAYEDGROUP) = "Delayed Group"

    ! Initialize names for scores
    score_names(abs(SCORE_FLUX))               = "Flux"
    score_names(abs(SCORE_TOTAL))              = "Total Reaction Rate"
    score_names(abs(SCORE_SCATTER))            = "Scattering Rate"
    score_names(abs(SCORE_NU_SCATTER))         = "Scattering Production Rate"
    score_names(abs(SCORE_TRANSPORT))          = "Transport Rate"
    score_names(abs(SCORE_N_1N))               = "(n,1n) Rate"
    score_names(abs(SCORE_ABSORPTION))         = "Absorption Rate"
    score_names(abs(SCORE_FISSION))            = "Fission Rate"
    score_names(abs(SCORE_NU_FISSION))         = "Nu-Fission Rate"
    score_names(abs(SCORE_KAPPA_FISSION))      = "Kappa-Fission Rate"
    score_names(abs(SCORE_EVENTS))             = "Events"
    score_names(abs(SCORE_FLUX_YN))            = "Flux Moment"
    score_names(abs(SCORE_TOTAL_YN))           = "Total Reaction Rate Moment"
    score_names(abs(SCORE_SCATTER_N))          = "Scattering Rate Moment"
    score_names(abs(SCORE_SCATTER_PN))         = "Scattering Rate Moment"
    score_names(abs(SCORE_SCATTER_YN))         = "Scattering Rate Moment"
    score_names(abs(SCORE_NU_SCATTER_N))       = "Scattering Prod. Rate Moment"
    score_names(abs(SCORE_NU_SCATTER_PN))      = "Scattering Prod. Rate Moment"
    score_names(abs(SCORE_NU_SCATTER_YN))      = "Scattering Prod. Rate Moment"
    score_names(abs(SCORE_DELAYED_NU_FISSION)) = "Delayed-Nu-Fission Rate"
    score_names(abs(SCORE_INVERSE_VELOCITY))   = "Flux-Weighted Inverse Velocity"

    ! Create filename for tally output
    filename = trim(path_output) // "tallies.out"

    ! Open tally file for writing
    open(FILE=filename, NEWUNIT=unit_tally, STATUS='replace', ACTION='write')

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
      if (t % name == "") then
        call header("TALLY " // trim(to_str(t % id)), unit=unit_tally, &
             level=3)
      else
        call header("TALLY " // trim(to_str(t % id)) // ": " &
             // trim(t % name), unit=unit_tally, level=3)
      endif

      ! Handle surface current tallies separately
      if (t % type == TALLY_SURFACE_CURRENT) then
        call write_surface_current(t, unit_tally)
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
            write(UNIT=unit_tally, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                 trim(filter_name(type)), trim(get_label(t, j))
            indent = indent + 2
            j = j + 1
          end if

        end do find_bin

        ! Print filter information
        if (t % n_filters > 0) then
          type = t % filters(j) % type
          write(UNIT=unit_tally, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
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
            write(UNIT=unit_tally, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                 "Total Material"
          else
            i_listing = nuclides(i_nuclide) % listing
            write(UNIT=unit_tally, FMT='(1X,2A,1X,A)') repeat(" ", indent), &
                 trim(xs_listings(i_listing) % alias)
          end if

          indent = indent + 2
          k = 0
          do l = 1, t % n_user_score_bins
            k = k + 1
            score_index = score_index + 1
            select case(t % score_bins(k))
            case (SCORE_SCATTER_N, SCORE_NU_SCATTER_N)
              score_name = 'P' // trim(to_str(t % moment_order(k))) // " " // &
                   score_names(abs(t % score_bins(k)))
              write(UNIT=unit_tally, FMT='(1X,2A,1X,A,"+/- ",A)') &
                   repeat(" ", indent), score_name, &
                   to_str(t % results(score_index,filter_index) % sum), &
                   trim(to_str(t % results(score_index,filter_index) % sum_sq))
            case (SCORE_SCATTER_PN, SCORE_NU_SCATTER_PN)
              score_index = score_index - 1
              do n_order = 0, t % moment_order(k)
                score_index = score_index + 1
                score_name = 'P' // trim(to_str(n_order)) //  " " //&
                     score_names(abs(t % score_bins(k)))
                write(UNIT=unit_tally, FMT='(1X,2A,1X,A,"+/- ",A)') &
                     repeat(" ", indent), score_name, &
                     to_str(t % results(score_index,filter_index) % sum), &
                     trim(to_str(t % results(score_index,filter_index) &
                     % sum_sq))
              end do
              k = k + t % moment_order(k)
            case (SCORE_SCATTER_YN, SCORE_NU_SCATTER_YN, SCORE_FLUX_YN, &
                  SCORE_TOTAL_YN)
              score_index = score_index - 1
              do n_order = 0, t % moment_order(k)
                do nm_order = -n_order, n_order
                  score_index = score_index + 1
                  score_name = 'Y' // trim(to_str(n_order)) // ',' // &
                       trim(to_str(nm_order)) // " " &
                       // score_names(abs(t % score_bins(k)))
                  write(UNIT=unit_tally, FMT='(1X,2A,1X,A,"+/- ",A)') &
                       repeat(" ", indent), score_name, &
                       to_str(t % results(score_index,filter_index) % sum), &
                       trim(to_str(t % results(score_index,filter_index)&
                       % sum_sq))
                end do
              end do
              k = k + (t % moment_order(k) + 1)**2 - 1
            case default
              if (t % score_bins(k) > 0) then
                score_name = reaction_name(t % score_bins(k))
              else
                score_name = score_names(abs(t % score_bins(k)))
              end if
              write(UNIT=unit_tally, FMT='(1X,2A,1X,A,"+/- ",A)') &
                   repeat(" ", indent), score_name, &
                   to_str(t % results(score_index,filter_index) % sum), &
                   trim(to_str(t % results(score_index,filter_index) % sum_sq))
            end select
          end do
          indent = indent - 2

        end do
        indent = indent - 2

        if (t % n_filters == 0) exit print_bin

      end do print_bin

    end do TALLY_LOOP

    close(UNIT=unit_tally)

  end subroutine write_tallies

!===============================================================================
! WRITE_SURFACE_CURRENT writes out surface current tallies over a mesh to the
! tallies.out file.
!===============================================================================

  subroutine write_surface_current(t, unit_tally)
    type(TallyObject), pointer :: t
    integer, intent(in) :: unit_tally

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
    type(RegularMesh), pointer :: m

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
          write(UNIT=unit_tally, FMT='(1X,A)') trim(string)

          do l = 1, n
            if (print_ebin) then
              ! Set incoming energy bin
              matching_bins(i_filter_ein) = l

              ! Write incoming energy bin
              write(UNIT=unit_tally, FMT='(3X,A,1X,A)') &
                   "Incoming Energy", trim(get_label(t, i_filter_ein))
            end if

            ! Left Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i-1, j, k /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_RIGHT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Left", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_RIGHT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Left", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            ! Right Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i, j, k /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_RIGHT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Right", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_RIGHT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Right", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            ! Back Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i, j-1, k /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_FRONT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Back", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_FRONT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Back", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            ! Front Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i, j, k /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_FRONT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Front", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_FRONT
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Front", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            ! Bottom Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i, j, k-1 /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_TOP
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Outgoing Current to Bottom", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_TOP
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Bottom", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            ! Top Surface
            matching_bins(i_filter_mesh) = &
                 mesh_indices_to_bin(m, (/ i, j, k /) + 1, .true.)
            matching_bins(i_filter_surf) = IN_TOP
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
                 "Incoming Current from Top", &
                 to_str(t % results(1,filter_index) % sum), &
                 trim(to_str(t % results(1,filter_index) % sum_sq))

            matching_bins(i_filter_surf) = OUT_TOP
            filter_index = sum((matching_bins(1:t%n_filters) - 1) * t % stride) + 1
            write(UNIT=unit_tally, FMT='(5X,A,T35,A,"+/- ",A)') &
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
    character(100)              :: label    ! user-specified identifier

    integer :: i      ! index in cells/surfaces/etc array
    integer :: bin
    integer :: offset
    integer, allocatable :: ijk(:) ! indices in mesh
    real(8)              :: E0     ! lower bound for energy bin
    real(8)              :: E1     ! upper bound for energy bin
    type(RegularMesh), pointer :: m
    type(Universe), pointer :: univ

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
    case (FILTER_DISTRIBCELL)
      label = ''
      univ => universes(BASE_UNIVERSE)
      offset = 0
      call find_offset(t % filters(i_filter) % offset, &
           t % filters(i_filter) % int_bins(1), &
           univ, bin-1, offset, label)
    case (FILTER_SURFACE)
      i = t % filters(i_filter) % int_bins(bin)
      label = to_str(surfaces(i)%obj%id)
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
    case (FILTER_ENERGYIN, FILTER_ENERGYOUT, FILTER_MU, FILTER_POLAR, &
          FILTER_AZIMUTHAL)
      E0 = t % filters(i_filter) % real_bins(bin)
      E1 = t % filters(i_filter) % real_bins(bin + 1)
      label = "[" // trim(to_str(E0)) // ", " // trim(to_str(E1)) // ")"
    case (FILTER_DELAYEDGROUP)
      i = t % filters(i_filter) % int_bins(bin)
      label = to_str(i)
    end select

  end function get_label

!===============================================================================
! FIND_OFFSET uses a given map number, a target cell ID, and a target offset
! to build a string which is the path from the base universe to the target cell
! with the given offset
!===============================================================================

  recursive subroutine find_offset(map, goal, univ, final, offset, path)

    integer, intent(in) :: map                   ! Index in maps vector
    integer, intent(in) :: goal                  ! The target cell ID
    type(Universe), pointer, intent(in) :: univ  ! Universe to begin search
    integer, intent(in) :: final                 ! Target offset
    integer, intent(inout) :: offset             ! Current offset
    character(100) :: path                       ! Path to offset

    integer :: i, j                 ! Index over cells
    integer :: k, l, m              ! Indices in lattice
    integer :: old_k, old_l, old_m  ! Previous indices in lattice
    integer :: n_x, n_y, n_z        ! Lattice cell array dimensions
    integer :: n                    ! Number of cells to search
    integer :: cell_index           ! Index in cells array
    integer :: lat_offset           ! Offset from lattice
    integer :: temp_offset          ! Looped sum of offsets
    logical :: this_cell = .false.  ! Advance in this cell?
    logical :: later_cell = .false. ! Fill cells after this one?
    type(Cell),     pointer:: c           ! Pointer to current cell
    type(Universe), pointer :: next_univ  ! Next universe to loop through
    class(Lattice), pointer :: lat        ! Pointer to current lattice

    n = univ % n_cells

    ! Write to the geometry stack
    if (univ%id == 0) then
      path = trim(path) // to_str(univ%id)
    else
      path = trim(path) // "->" // to_str(univ%id)
    end if

    ! Look through all cells in this universe
    do i = 1, n

      cell_index = univ % cells(i)
      c => cells(cell_index)

      ! If the cell ID matches the goal and the offset matches final,
      ! write to the geometry stack
      if (cell_dict % get_key(c % id) == goal .AND. offset == final) then
        path = trim(path) // "->" // to_str(c%id)
        return
      end if

    end do

    ! Find the fill cell or lattice cell that we need to enter
    do i = 1, n

      later_cell = .false.

      cell_index = univ % cells(i)
      c => cells(cell_index)

      this_cell = .false.

      ! If we got here, we still think the target is in this universe
      ! or further down, but it's not this exact cell.
      ! Compare offset to next cell to see if we should enter this cell
      if (i /= n) then

        do j = i+1, n

          cell_index = univ % cells(j)
          c => cells(cell_index)

          ! Skip normal cells which do not have offsets
          if (c % type == CELL_NORMAL) then
            cycle
          end if

          ! Break loop once we've found the next cell with an offset
          exit
        end do

        ! Ensure we didn't just end the loop by iteration
        if (c % type /= CELL_NORMAL) then

          ! There are more cells in this universe that it could be in
          later_cell = .true.

          ! Two cases, lattice or fill cell
          if (c % type == CELL_FILL) then
            temp_offset = c % offset(map)

          ! Get the offset of the first lattice location
          else
            lat => lattices(c % fill) % obj
            temp_offset = lat % offset(map, 1, 1, 1)
          end if

          ! If the final offset is in the range of offset - temp_offset+offset
          ! then the goal is in this cell
          if (final < temp_offset + offset) then
            this_cell = .true.
          end if
        end if
      end if

      if (n == 1 .and. c % type /= CELL_NORMAL) then
        this_cell = .true.
      end if

      if (.not. later_cell) then
        this_cell = .true.
      end if

      ! Get pointer to THIS cell because target must be in this cell
      if (this_cell) then

        cell_index = univ % cells(i)
        c => cells(cell_index)

        path = trim(path) // "->" // to_str(c%id)

        ! ====================================================================
        ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
        if (c % type == CELL_FILL) then

          ! Enter this cell to update the current offset
          offset = c % offset(map) + offset

          next_univ => universes(c % fill)
          call find_offset(map, goal, next_univ, final, offset, path)
          return

        ! ====================================================================
        ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
        elseif (c % type == CELL_LATTICE) then

          ! Set current lattice
          lat => lattices(c % fill) % obj

          select type (lat)

          ! ==================================================================
          ! RECTANGULAR LATTICES
          type is (RectLattice)

            ! Write to the geometry stack
            path = trim(path) // "->" // to_str(lat%id)

            n_x = lat % n_cells(1)
            n_y = lat % n_cells(2)
            n_z = lat % n_cells(3)
            old_m = 1
            old_l = 1
            old_k = 1

            ! Loop over lattice coordinates
            do k = 1, n_x
              do l = 1, n_y
                do m = 1, n_z

                  if (final >= lat % offset(map, k, l, m) + offset) then
                    if (k == n_x .and. l == n_y .and. m == n_z) then
                      ! This is last lattice cell, so target must be here
                      lat_offset = lat % offset(map, k, l, m)
                      offset = offset + lat_offset
                      next_univ => universes(lat % universes(k, l, m))
                      path = trim(path) // "(" // trim(to_str(k)) // &
                           "," // trim(to_str(l)) // "," // &
                           trim(to_str(m)) // ")"
                      call find_offset(map, goal, next_univ, final, offset, path)
                      return
                    else
                      old_m = m
                      old_l = l
                      old_k = k
                      cycle
                    end if
                  else
                    ! Target is at this lattice position
                    lat_offset = lat % offset(map, old_k, old_l, old_m)
                    offset = offset + lat_offset
                    next_univ => universes(lat % universes(old_k, old_l, old_m))
                    path = trim(path) // "(" // trim(to_str(old_k)) // &
                         "," // trim(to_str(old_l)) // "," // &
                         trim(to_str(old_m)) // ")"
                    call find_offset(map, goal, next_univ, final, offset, path)
                    return
                  end if

                end do
              end do
            end do

          ! ==================================================================
          ! HEXAGONAL LATTICES
          type is (HexLattice)

            ! Write to the geometry stack
            path = trim(path) // "->" // to_str(lat%id)

            n_z = lat % n_axial
            n_y = 2 * lat % n_rings - 1
            n_x = 2 * lat % n_rings - 1
            old_m = 1
            old_l = 1
            old_k = 1

            ! Loop over lattice coordinates
            do m = 1, n_z
              do l = 1, n_y
                do k = 1, n_x

                  ! This array position is never used
                  if (k + l < lat % n_rings + 1) then
                    cycle
                  ! This array position is never used
                  else if (k + l > 3*lat % n_rings - 1) then
                    cycle
                  end if

                  if (final >= lat % offset(map, k, l, m) + offset) then
                    if (k == lat % n_rings .and. l == n_y .and. m == n_z) then
                      ! This is last lattice cell, so target must be here
                      lat_offset = lat % offset(map, k, l, m)
                      offset = offset + lat_offset
                      next_univ => universes(lat % universes(k, l, m))
                      path = trim(path) // "(" // &
                           trim(to_str(k - lat % n_rings)) // "," // &
                           trim(to_str(l - lat % n_rings)) // "," // &
                           trim(to_str(m)) // ")"
                      call find_offset(map, goal, next_univ, final, offset, &
                                       path)
                      return
                    else
                      old_m = m
                      old_l = l
                      old_k = k
                      cycle
                    end if
                  else
                    ! Target is at this lattice position
                    lat_offset = lat % offset(map, old_k, old_l, old_m)
                    offset = offset + lat_offset
                    next_univ => universes(lat % universes(old_k, old_l, old_m))
                    path = trim(path) // "(" // &
                         trim(to_str(old_k - lat % n_rings)) // "," // &
                         trim(to_str(old_l - lat % n_rings)) // "," // &
                         trim(to_str(old_m)) // ")"
                    call find_offset(map, goal, next_univ, final, offset, path)
                    return
                  end if

                end do
              end do
            end do

          end select

        end if
      end if
    end do
  end subroutine find_offset

end module output
