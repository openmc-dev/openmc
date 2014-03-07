module ace

  use ace_header,       only: Nuclide, Reaction, SAlphaBeta, XsListing, &
                              DistEnergy
  use constants
  use endf,             only: reaction_name
  use error,            only: fatal_error, warning
  use fission,          only: nu_total
  use global
  use material_header,  only: Material
  use output,           only: write_message
  use set_header,       only: SetChar
  use string,           only: to_str

  implicit none

  integer :: NXS(16)             ! Descriptors for ACE XSS tables
  integer :: JXS(32)             ! Pointers into ACE XSS tables
  real(8), allocatable :: XSS(:) ! Cross section data
  integer :: XSS_index           ! current index in XSS data

  private :: NXS
  private :: JXS
  private :: XSS

contains

!===============================================================================
! READ_XS reads all the cross sections for the problem and stores them in
! nuclides and sab_tables arrays
!===============================================================================

  subroutine read_xs()

    integer :: i            ! index in materials array
    integer :: j            ! index over nuclides in material
    integer :: k            ! index over S(a,b) tables in material
    integer :: i_listing    ! index in xs_listings array
    integer :: i_nuclide    ! index in nuclides
    integer :: i_sab        ! index in sab_tables
    integer :: m            ! position for sorting
    integer :: temp_nuclide ! temporary value for sorting
    integer :: temp_table   ! temporary value for sorting
    character(12)  :: name  ! name of isotope, e.g. 92235.03c
    character(12)  :: alias ! alias of nuclide, e.g. U-235.03c
    type(Material),   pointer :: mat => null()
    type(Nuclide),    pointer :: nuc => null()
    type(SAlphaBeta), pointer :: sab => null()
    type(SetChar) :: already_read

    ! allocate arrays for ACE table storage and cross section cache
    allocate(nuclides(n_nuclides_total))
    allocate(sab_tables(n_sab_tables))
!$omp parallel
    allocate(micro_xs(n_nuclides_total))
!$omp end parallel

    ! ==========================================================================
    ! READ ALL ACE CROSS SECTION TABLES

    ! Loop over all files
    MATERIAL_LOOP: do i = 1, n_materials
      mat => materials(i)

      NUCLIDE_LOOP: do j = 1, mat % n_nuclides
        name = mat % names(j)

        if (.not. already_read % contains(name)) then
          i_listing = xs_listing_dict % get_key(name)
          i_nuclide = nuclide_dict % get_key(name)
          name  = xs_listings(i_listing) % name
          alias = xs_listings(i_listing) % alias

          ! Keep track of what listing is associated with this nuclide
          nuc => nuclides(i_nuclide)
          nuc % listing = i_listing

          ! Read the ACE table into the appropriate entry on the nuclides
          ! array
          call read_ace_table(i_nuclide, i_listing)

          ! Add name and alias to dictionary
          call already_read % add(name)
          call already_read % add(alias)
        end if
      end do NUCLIDE_LOOP

      SAB_LOOP: do k = 1, mat % n_sab
        ! Get name of S(a,b) table
        name = mat % sab_names(k)

        if (.not. already_read % contains(name)) then
          i_listing = xs_listing_dict % get_key(name)
          i_sab  = sab_dict % get_key(name)

          ! Read the ACE table into the appropriate entry on the sab_tables
          ! array
          call read_ace_table(i_sab, i_listing)

          ! Add name to dictionary
          call already_read % add(name)
        end if
      end do SAB_LOOP
    end do MATERIAL_LOOP

    ! ==========================================================================
    ! ASSIGN S(A,B) TABLES TO SPECIFIC NUCLIDES WITHIN MATERIALS

    MATERIAL_LOOP2: do i = 1, n_materials
      ! Get pointer to material
      mat => materials(i)

      ASSIGN_SAB: do k = 1, mat % n_sab
        ! In order to know which nuclide the S(a,b) table applies to, we need to
        ! search through the list of nuclides for one which has a matching zaid
        sab => sab_tables(mat % i_sab_tables(k))

        ! Loop through nuclides and find match
        FIND_NUCLIDE: do j = 1, mat % n_nuclides
          if (any(sab % zaid == nuclides(mat % nuclide(j)) % zaid)) then
            mat % i_sab_nuclides(k) = j
            exit FIND_NUCLIDE
          end if
        end do FIND_NUCLIDE

        ! Check to make sure S(a,b) table matched a nuclide
        if (mat % i_sab_nuclides(k) == NONE) then
          message = "S(a,b) table " // trim(mat % sab_names(k)) // " did not &
               &match any nuclide on material " // trim(to_str(mat % id))
          call fatal_error()
        end if
      end do ASSIGN_SAB

      ! If there are multiple S(a,b) tables, we need to make sure that the
      ! entries in i_sab_nuclides are sorted or else they won't be applied
      ! correctly in the cross_section module. The algorithm here is a simple
      ! insertion sort -- don't need anything fancy!

      if (mat % n_sab > 1) then
        SORT_SAB: do k = 2, mat % n_sab
          ! Save value to move
          m = k
          temp_nuclide = mat % i_sab_nuclides(k)
          temp_table   = mat % i_sab_tables(k)

          MOVE_OVER: do
            ! Check if insertion value is greater than (m-1)th value
            if (temp_nuclide >= mat % i_sab_nuclides(m-1)) exit

            ! Move values over until hitting one that's not larger
            mat % i_sab_nuclides(m) = mat % i_sab_nuclides(m-1)
            mat % i_sab_tables(m)   = mat % i_sab_tables(m-1)
            m = m - 1

            ! Exit if we've reached the beginning of the list
            if (m == 1) exit
          end do MOVE_OVER

          ! Put the original value into its new position
          mat % i_sab_nuclides(m) = temp_nuclide
          mat % i_sab_tables(m)   = temp_table
        end do SORT_SAB
      end if

      ! Deallocate temporary arrays for names of nuclides and S(a,b) tables
      if (allocated(mat % names)) deallocate(mat % names)
      if (allocated(mat % sab_names)) deallocate(mat % sab_names)

    end do MATERIAL_LOOP2

    ! Avoid some valgrind leak errors
    call already_read % clear()

  end subroutine read_xs

!===============================================================================
! READ_ACE_TABLE reads a single cross section table in either ASCII or binary
! format. This routine reads the header data for each table and then calls
! appropriate subroutines to parse the actual data.
!===============================================================================

  subroutine read_ace_table(i_table, i_listing)

    integer, intent(in) :: i_table   ! index in nuclides/sab_tables
    integer, intent(in) :: i_listing ! index in xs_listings

    integer       :: i             ! loop index for XSS records
    integer       :: j, j1, j2     ! indices in XSS
    integer       :: record_length ! Fortran record length
    integer       :: location      ! location of ACE table
    integer       :: entries       ! number of entries on each record
    integer       :: length        ! length of ACE table
    integer       :: in = 7        ! file unit
    integer       :: zaids(16)     ! list of ZAIDs (only used for S(a,b))
    integer       :: filetype      ! filetype (ASCII or BINARY)
    real(8)       :: kT            ! temperature of table
    real(8)       :: awrs(16)      ! list of atomic weight ratios (not used)
    real(8)       :: awr           ! atomic weight ratio for table
    logical       :: file_exists   ! does ACE library exist?
    character(7)  :: readable      ! is ACE library readable?
    character(10) :: name          ! name of ACE table
    character(10) :: date_         ! date ACE library was processed
    character(10) :: mat           ! material identifier
    character(70) :: comment       ! comment for ACE table
    character(MAX_FILE_LEN) :: filename ! path to ACE cross section library
    type(Nuclide),   pointer :: nuc => null()
    type(SAlphaBeta), pointer :: sab => null()
    type(XsListing), pointer :: listing => null()

    ! determine path, record length, and location of table
    listing => xs_listings(i_listing)
    filename      = listing % path
    record_length = listing % recl
    location      = listing % location
    entries       = listing % entries
    filetype      = listing % filetype

    ! Check if ACE library exists and is readable
    inquire(FILE=filename, EXIST=file_exists, READ=readable)
    if (.not. file_exists) then
      message = "ACE library '" // trim(filename) // "' does not exist!"
      call fatal_error()
    elseif (readable(1:3) == 'NO') then
      message = "ACE library '" // trim(filename) // "' is not readable! &
           &Change file permissions with chmod command."
      call fatal_error()
    end if

    ! display message
    message = "Loading ACE cross section table: " // listing % name
    call write_message(6)

    if (filetype == ASCII) then
      ! =======================================================================
      ! READ ACE TABLE IN ASCII FORMAT

      ! Find location of table
      open(UNIT=in, FILE=filename, STATUS='old', ACTION='read')
      rewind(UNIT=in)
      do i = 1, location - 1
        read(UNIT=in, FMT=*)
      end do

      ! Read first line of header
      read(UNIT=in, FMT='(A10,2G12.0,1X,A10)') name, awr, kT, date_

      ! Check that correct xs was found -- if cross_sections.xml is broken, the
      ! location of the table may be wrong
      if(adjustl(name) /= adjustl(listing % name)) then
        message = "XS listing entry " // trim(listing % name) // " did not &
             &match ACE data, " // trim(name) // " found instead."
        call fatal_error()
      end if

      ! Read more header and NXS and JXS
      read(UNIT=in, FMT=100) comment, mat, &
           (zaids(i), awrs(i), i=1,16), NXS, JXS
100   format(A70,A10/4(I7,F11.0)/4(I7,F11.0)/4(I7,F11.0)/4(I7,F11.0)/&
           ,8I9/8I9/8I9/8I9/8I9/8I9)

      ! determine table length
      length = NXS(1)
      allocate(XSS(length))

      ! Read XSS array
      read(UNIT=in, FMT='(4G20.0)') XSS

      ! Close ACE file
      close(UNIT=in)

    elseif (filetype == BINARY) then
      ! =======================================================================
      ! READ ACE TABLE IN BINARY FORMAT

      ! Open ACE file
      open(UNIT=in, FILE=filename, STATUS='old', ACTION='read', &
           ACCESS='direct', RECL=record_length)

      ! Read all header information
      read(UNIT=in, REC=location) name, awr, kT, date_, &
           comment, mat, (zaids(i), awrs(i), i=1,16), NXS, JXS

      ! determine table length
      length = NXS(1)
      allocate(XSS(length))

      ! Read remaining records with XSS
      do i = 1, (length + entries - 1)/entries
        j1 = 1 + (i-1)*entries
        j2 = min(length, j1 + entries - 1)
        read(UNIT=IN, REC=location + i) (XSS(j), j=j1,j2)
      end do

      ! Close ACE file
      close(UNIT=in)
    end if

    ! ==========================================================================
    ! PARSE DATA BASED ON NXS, JXS, AND XSS ARRAYS

    select case(listing % type)
    case (ACE_NEUTRON)
      nuc => nuclides(i_table)
      nuc % name = name
      nuc % awr  = awr
      nuc % kT   = kT
      nuc % zaid = NXS(2)

      ! read all blocks
      call read_esz(nuc)
      call read_nu_data(nuc)
      call read_reactions(nuc)
      call read_angular_dist(nuc)
      call read_energy_dist(nuc)
      call read_unr_res(nuc)

      ! Currently subcritical fixed source calculations are not allowed. Thus,
      ! if any fissionable material is found in a fixed source calculation,
      ! abort the run.
      if (run_mode == MODE_FIXEDSOURCE .and. nuc % fissionable) then
        message = "Cannot have fissionable material in a fixed source run."
        call fatal_error()
      end if

      ! for fissionable nuclides, precalculate microscopic nu-fission cross
      ! sections so that we don't need to call the nu_total function during
      ! cross section lookups

      if (nuc % fissionable) call generate_nu_fission(nuc)

    case (ACE_THERMAL)
      sab => sab_tables(i_table)
      sab % name = name
      sab % awr  = awr
      sab % kT   = kT
      ! Find sab % n_zaid
      do i = 1, 16
        if (zaids(i) == 0) then
          sab % n_zaid = i - 1
          exit
        end if
      end do
      allocate(sab % zaid(sab % n_zaid))
      sab % zaid = zaids(1: sab % n_zaid)

      call read_thermal_data(sab)
    end select

    deallocate(XSS)
    if(associated(nuc)) nullify(nuc)
    if(associated(sab)) nullify(sab)

  end subroutine read_ace_table

!===============================================================================
! READ_ESZ - reads through the ESZ block. This block contains the energy grid,
! total xs, absorption xs, elastic scattering xs, and heating numbers.
!===============================================================================

  subroutine read_esz(nuc)

    type(Nuclide), pointer :: nuc

    integer :: NE ! number of energy points for total and elastic cross sections

    ! determine number of energy points
    NE = NXS(3)
    nuc % n_grid = NE

    ! allocate storage for energy grid and cross section arrays
    allocate(nuc % energy(NE))
    allocate(nuc % total(NE))
    allocate(nuc % elastic(NE))
    allocate(nuc % fission(NE))
    allocate(nuc % nu_fission(NE))
    allocate(nuc % absorption(NE))

    ! initialize cross sections
    nuc % total      = ZERO
    nuc % elastic    = ZERO
    nuc % fission    = ZERO
    nuc % nu_fission = ZERO
    nuc % absorption = ZERO

    ! Read data from XSS -- only the energy grid, elastic scattering and heating
    ! cross section values are actually read from here. The total and absorption
    ! cross sections are reconstructed from the partial reaction data.

    XSS_index = 1
    nuc % energy = get_real(NE)

    ! Skip total and absorption
    XSS_index = XSS_index + 2*NE

    ! Continue reading elastic scattering and heating
    nuc % elastic = get_real(NE)

  end subroutine read_esz

!===============================================================================
! READ_NU_DATA reads data given on the number of neutrons emitted from fission
! as a function of the incoming energy of a neutron. This data may be broken
! down into prompt and delayed neutrons emitted as well.
!===============================================================================

  subroutine read_nu_data(nuc)

    type(Nuclide), pointer :: nuc

    integer :: i      ! loop index
    integer :: JXS2   ! location for fission nu data
    integer :: JXS24  ! location for delayed neutron data
    integer :: KNU    ! location for nu data
    integer :: LNU    ! type of nu data (polynomial or tabular)
    integer :: NC     ! number of polynomial coefficients
    integer :: NR     ! number of interpolation regions
    integer :: NE     ! number of energies
    integer :: NPCR   ! number of delayed neutron precursor groups
    integer :: LED    ! location of energy distribution locators
    integer :: LDIS   ! location of all energy distributions
    integer :: LOCC   ! location of energy distributions for given MT
    integer :: lc     ! locator
    integer :: length ! length of data to allocate
    type(DistEnergy), pointer :: edist => null()

    JXS2  = JXS(2)
    JXS24 = JXS(24)

    if (JXS2 == 0) then
      ! =======================================================================
      ! NO PROMPT/TOTAL NU DATA
      nuc % nu_t_type = NU_NONE
      nuc % nu_p_type = NU_NONE

    elseif (XSS(JXS2) > 0) then
      ! =======================================================================
      ! PROMPT OR TOTAL NU DATA
      KNU = JXS2
      LNU = int(XSS(KNU))
      if (LNU == 1) then
        ! Polynomial data
        nuc % nu_t_type = NU_POLYNOMIAL
        nuc % nu_p_type = NU_NONE

        ! allocate determine how many coefficients for polynomial
        NC = int(XSS(KNU+1))
        length = NC + 1
      elseif (LNU == 2) then
        ! Tabular data
        nuc % nu_t_type = NU_TABULAR
        nuc % nu_p_type = NU_NONE

        ! determine number of interpolation regions and number of energies
        NR = int(XSS(KNU+1))
        NE = int(XSS(KNU+2+2*NR))
        length = 2 + 2*NR + 2*NE
      end if

      ! allocate space for nu data storage
      allocate(nuc % nu_t_data(length))

      ! read data -- for polynomial, this is the number of coefficients and the
      ! coefficients themselves, and for tabular, this is interpolation data
      ! and tabular E/nu
      XSS_index = KNU + 1
      nuc % nu_t_data = get_real(length)

    elseif (XSS(JXS2) < 0) then
      ! =======================================================================
      ! PROMPT AND TOTAL NU DATA -- read prompt data first
      KNU = JXS2 + 1
      LNU = int(XSS(KNU))
      if (LNU == 1) then
        ! Polynomial data
        nuc % nu_p_type = NU_POLYNOMIAL

        ! allocate determine how many coefficients for polynomial
        NC = int(XSS(KNU+1))
        length = NC + 1
      elseif (LNU == 2) then
        ! Tabular data
        nuc % nu_p_type = NU_TABULAR

        ! determine number of interpolation regions and number of energies
        NR = int(XSS(KNU+1))
        NE = int(XSS(KNU+2+2*NR))
        length = 2 + 2*NR + 2*NE
      end if

      ! allocate space for nu data storage
      allocate(nuc % nu_p_data(length))

      ! read data
      XSS_index = KNU + 1
      nuc % nu_p_data = get_real(length)

      ! Now read total nu data
      KNU = JXS2 + int(abs(XSS(JXS2))) + 1
      LNU = int(XSS(KNU))
      if (LNU == 1) then
        ! Polynomial data
        nuc % nu_t_type = NU_POLYNOMIAL

        ! allocate determine how many coefficients for polynomial
        NC = int(XSS(KNU+1))
        length = NC + 1
      elseif (LNU == 2) then
        ! Tabular data
        nuc % nu_t_type = NU_TABULAR

        ! determine number of interpolation regions and number of energies
        NR = int(XSS(KNU+1))
        NE = int(XSS(KNU+2+2*NR))
        length = 2 + 2*NR + 2*NE
      end if

      ! allocate space for nu data storage
      allocate(nuc % nu_t_data(length))

      ! read data
      XSS_index = KNU + 1
      nuc % nu_t_data = get_real(length)
    end if

    if (JXS24 > 0) then
      ! =======================================================================
      ! DELAYED NU DATA

      nuc % nu_d_type = NU_TABULAR
      KNU = JXS24

      ! determine size of tabular delayed nu data
      NR = int(XSS(KNU+1))
      NE = int(XSS(KNU+2+2*NR))
      length = 2 + 2*NR + 2*NE

      ! allocate space for delayed nu data
      allocate(nuc % nu_d_data(length))

      ! read delayed nu data
      XSS_index = KNU + 1
      nuc % nu_d_data = get_real(length)

      ! =======================================================================
      ! DELAYED NEUTRON ENERGY DISTRIBUTION

      ! Allocate space for secondary energy distribution
      NPCR = NXS(8)
      nuc % n_precursor = NPCR
      allocate(nuc % nu_d_edist(NPCR))

      LED  = JXS(26)
      LDIS = JXS(27)

      ! Loop over all delayed neutron precursor groups
      do i = 1, NPCR
        ! find location of energy distribution data
        LOCC = int(XSS(LED + i - 1))

        ! read energy distribution data
        edist => nuc % nu_d_edist(i)
        call get_energy_dist(edist, LOCC, .true.)
      end do

      ! =======================================================================
      ! DELAYED NEUTRON PRECUSOR YIELDS AND CONSTANTS

      ! determine length of all precursor constants/yields/interp data
      length = 0
      lc = JXS(25)
      do i = 1, NPCR
        NR = int(XSS(lc + length + 1))
        NE = int(XSS(lc + length + 2 + 2*NR))
        length = length + 3 + 2*NR + 2*NE
      end do

      ! allocate space for precusor data
      allocate(nuc % nu_d_precursor_data(length))

      ! read delayed neutron precursor data
      XSS_index = lc
      nuc % nu_d_precursor_data = get_real(length)

    else
      nuc % nu_d_type = NU_NONE
    end if

  end subroutine read_nu_data

!===============================================================================
! READ_REACTIONS - Get the list of reaction MTs for this cross-section
! table. The MT values are somewhat arbitrary. Also read in Q-values, neutron
! multiplicities, and cross-sections.
!===============================================================================

  subroutine read_reactions(nuc)

    type(Nuclide), pointer :: nuc

    integer :: i         ! loop indices
    integer :: i_fission ! index in nuc % index_fission
    integer :: LMT       ! index of MT list in XSS
    integer :: NMT       ! Number of reactions
    integer :: JXS4      ! index of Q values in XSS
    integer :: JXS5      ! index of neutron multiplicities in XSS
    integer :: JXS7      ! index of reactions cross-sections in XSS
    integer :: LXS       ! location of cross-section locators
    integer :: LOCA      ! location of cross-section for given MT
    integer :: IE        ! reaction's starting index on energy grid
    integer :: NE        ! number of energies for reaction
    type(Reaction), pointer :: rxn => null()

    LMT  = JXS(3)
    JXS4 = JXS(4)
    JXS5 = JXS(5)
    LXS  = JXS(6)
    JXS7 = JXS(7)
    NMT  = NXS(4)

    ! allocate array of reactions. Add one since we need to include an elastic
    ! scattering channel
    nuc % n_reaction = NMT + 1
    allocate(nuc % reactions(NMT+1))

    ! Store elastic scattering cross-section on reaction one -- note that the
    ! sigma array is not allocated or stored for elastic scattering since it is
    ! already stored in nuc % elastic
    rxn => nuc % reactions(1)
    rxn % MT            = 2
    rxn % Q_value       = ZERO
    rxn % multiplicity  = 1
    rxn % threshold     = 1
    rxn % scatter_in_cm = .true.
    rxn % has_angle_dist = .false.
    rxn % has_energy_dist = .false.

    ! Add contribution of elastic scattering to total cross section
    nuc % total = nuc % total + nuc % elastic

    ! By default, set nuclide to not fissionable and then change if fission
    ! reactions are encountered
    nuc % fissionable = .false.
    nuc % has_partial_fission = .false.
    nuc % n_fission = 0
    i_fission = 0

    do i = 1, NMT
      rxn => nuc % reactions(i+1)

      ! set defaults
      rxn % has_angle_dist  = .false.
      rxn % has_energy_dist = .false.

      ! read MT number, Q-value, and neutrons produced
      rxn % MT            = int(XSS(LMT + i - 1))
      rxn % Q_value       = XSS(JXS4 + i - 1)
      rxn % multiplicity  = abs(nint(XSS(JXS5 + i - 1)))
      rxn % scatter_in_cm = (nint(XSS(JXS5 + i - 1)) < 0)

      ! read starting energy index
      LOCA = int(XSS(LXS + i - 1))
      IE   = int(XSS(JXS7 + LOCA - 1))
      rxn % threshold = IE

      ! read number of energies cross section values
      NE = int(XSS(JXS7 + LOCA))
      allocate(rxn % sigma(NE))
      XSS_index = JXS7 + LOCA + 1
      rxn % sigma = get_real(NE)

      ! Skip redundant reactions -- this includes total inelastic level
      ! scattering, gas production cross sections (MT=200+), and (n,p), (n,d),
      ! etc. reactions leaving the nucleus in an excited state
      if (rxn % MT == N_LEVEL .or. rxn % MT > N_DA) cycle

      ! Add contribution to total cross section
      nuc % total(IE:IE+NE-1) = nuc % total(IE:IE+NE-1) + rxn % sigma

      ! Add contribution to absorption cross section
      if (rxn % MT >= N_GAMMA .and. rxn % MT <= N_DA) then
        nuc % absorption(IE:IE+NE-1) = nuc % absorption(IE:IE+NE-1) + rxn % sigma
      end if

      ! Information about fission reactions
      if (rxn % MT == N_FISSION) then
        allocate(nuc % index_fission(1))
      elseif (rxn % MT == N_F) then
        allocate(nuc % index_fission(PARTIAL_FISSION_MAX))
        nuc % has_partial_fission = .true.
      end if

      ! Add contribution to fission cross section
      if (rxn % MT == N_FISSION .or. rxn % MT == N_F .or. rxn % MT == N_NF &
           .or. rxn % MT == N_2NF .or. rxn % MT == N_3NF) then
        nuc % fissionable = .true.
        nuc % fission(IE:IE+NE-1) = nuc % fission(IE:IE+NE-1) + rxn % sigma

        ! Also need to add fission cross sections to absorption
        nuc % absorption(IE:IE+NE-1) = nuc % absorption(IE:IE+NE-1) + rxn % sigma

        ! If total fission reaction is present, there's no need to store the
        ! reaction cross-section since it was copied to nuc % fission
        if (rxn % MT == N_FISSION) deallocate(rxn % sigma)

        ! Keep track of this reaction for easy searching later
        i_fission = i_fission + 1
        nuc % index_fission(i_fission) = i + 1
        nuc % n_fission = nuc % n_fission + 1
      end if
    end do

  end subroutine read_reactions

!===============================================================================
! READ_ANGULAR_DIST parses the angular distribution for each reaction with
! secondary neutrons
!===============================================================================

  subroutine read_angular_dist(nuc)

    type(Nuclide), pointer :: nuc

    integer :: JXS8   ! location of angular distribution locators
    integer :: JXS9   ! location of angular distributions
    integer :: LOCB   ! location of angular distribution for given MT
    integer :: NE     ! number of incoming energies
    integer :: NP     ! number of points for cosine distribution
    integer :: LC     ! locator
    integer :: i      ! index in reactions array
    integer :: j      ! index over incoming energies
    integer :: length ! length of data array to allocate
    type(Reaction), pointer :: rxn => null()

    JXS8 = JXS(8)
    JXS9 = JXS(9)

    ! loop over all reactions with secondary neutrons -- NXS(5) does not include
    ! elastic scattering
    do i = 1, NXS(5) + 1
      rxn => nuc%reactions(i)

      ! find location of angular distribution
      LOCB = int(XSS(JXS8 + i - 1))
      if (LOCB == -1) then
        ! Angular distribution data are specified through LAWi = 44 in the DLW
        ! block
        cycle
      elseif (LOCB == 0) then
        ! No angular distribution data are given for this reaction, isotropic
        ! scattering is asssumed (in CM if TY < 0 and in LAB if TY > 0)
        cycle
      end if
      rxn % has_angle_dist = .true.

      ! allocate space for incoming energies and locations
      NE = int(XSS(JXS9 + LOCB - 1))
      rxn % adist % n_energy = NE
      allocate(rxn % adist % energy(NE))
      allocate(rxn % adist % type(NE))
      allocate(rxn % adist % location(NE))

      ! read incoming energy grid and location of nucs
      XSS_index = JXS9 + LOCB
      rxn % adist % energy   = get_real(NE)
      rxn % adist % location = get_int(NE)

      ! determine dize of data block
      length = 0
      do j = 1, NE
        LC = rxn % adist % location(j)
        if (LC == 0) then
          ! isotropic
          rxn % adist % type(j) = ANGLE_ISOTROPIC
        elseif (LC > 0) then
          ! 32 equiprobable bins
          rxn % adist % type(j) = ANGLE_32_EQUI
          length = length + 33
        elseif (LC < 0) then
          ! tabular distribution
          rxn % adist % type(j) = ANGLE_TABULAR
          NP = int(XSS(JXS9 + abs(LC)))
          length = length + 2 + 3*NP
        end if
      end do

      ! allocate angular distribution data and read
      allocate(rxn % adist % data(length))

      ! read angular distribution -- currently this does not actually parse the
      ! angular distribution tables for each incoming energy, that must be done
      ! on-the-fly
      LC = rxn % adist % location(1)
      XSS_index = JXS9 + abs(LC) - 1
      rxn % adist % data = get_real(length)

      ! change location pointers since they are currently relative to JXS(9)
      LC = abs(rxn % adist % location(1))
      rxn % adist % location = abs(rxn % adist % location) - LC

    end do

  end subroutine read_angular_dist

!===============================================================================
! READ_ENERGY_DIST parses the secondary energy distribution for each reaction
! with seconary neutrons (except elastic scattering)
!===============================================================================

  subroutine read_energy_dist(nuc)

    type(Nuclide), pointer :: nuc

    integer :: LED   ! location of energy distribution locators
    integer :: LOCC  ! location of energy distributions for given MT
    integer :: i     ! loop index
    type(Reaction), pointer :: rxn => null()

    LED  = JXS(10)

    ! Loop over all reactions
    do i = 1, NXS(5)
      rxn => nuc % reactions(i+1) ! skip over elastic scattering
      rxn % has_energy_dist = .true.

      ! find location of energy distribution data
      LOCC = int(XSS(LED + i - 1))

      ! allocate energy distribution
      allocate(rxn % edist)

      ! read data for energy distribution
      call get_energy_dist(rxn % edist, LOCC)
    end do

  end subroutine read_energy_dist

!===============================================================================
! GET_ENERGY_DIST reads in data for a single law for an energy distribution and
! calls itself recursively if there are multiple energy distributions for a
! single reaction
!===============================================================================

  recursive subroutine get_energy_dist(edist, loc_law, delayed_n)

    type(DistEnergy), pointer :: edist     ! energy distribution
    integer, intent(in)       :: loc_law   ! locator for data
    logical, optional         :: delayed_n ! is this for delayed neutrons?

    integer :: LDIS   ! location of all energy distributions
    integer :: LNW    ! location of next energy distribution if multiple
    integer :: LAW    ! secondary energy distribution law
    integer :: NR     ! number of interpolation regions
    integer :: NE     ! number of incoming energies
    integer :: IDAT   ! location of first energy distribution for given MT
    integer :: lc     ! locator
    integer :: length ! length of data to allocate
    integer :: length_interp_data ! length of interpolation data

    ! determine location of energy distribution
    if (present(delayed_n)) then
      LDIS = JXS(27)
    else
      LDIS = JXS(11)
    end if

    ! locator for next law and information on this law
    LNW  = int(XSS(LDIS + loc_law - 1))
    LAW  = int(XSS(LDIS + loc_law))
    IDAT = int(XSS(LDIS + loc_law + 1))
    NR   = int(XSS(LDIS + loc_law + 2))
    edist % law = LAW
    edist % p_valid % n_regions = NR

    ! allocate space for ENDF interpolation parameters
    if (NR > 0) then
      allocate(edist % p_valid % nbt(NR))
      allocate(edist % p_valid % int(NR))
    end if

    ! read ENDF interpolation parameters
    XSS_index = LDIS + loc_law + 3
    if (NR > 0) then
      edist % p_valid % nbt = int(get_real(NR))
      edist % p_valid % int = int(get_real(NR))
    end if

    ! allocate space for law validity data
    NE = int(XSS(LDIS + loc_law + 3 + 2*NR))
    edist % p_valid % n_pairs = NE
    allocate(edist % p_valid % x(NE))
    allocate(edist % p_valid % y(NE))

    length_interp_data = 5 + 2*(NR + NE)

    ! read law validity data
    XSS_index = LDIS + loc_law + 4 + 2*NR
    edist % p_valid % x = get_real(NE)
    edist % p_valid % y = get_real(NE)

    ! Set index to beginning of IDAT array
    lc = LDIS + IDAT - 2

    ! determine length of energy distribution
    length = length_energy_dist(lc, LAW, loc_law, length_interp_data)

    ! allocate secondary energy distribution array
    allocate(edist % data(length))

    ! read secondary energy distribution
    XSS_index = lc + 1
    edist % data = get_real(length)

    ! read next energy distribution if present
    if (LNW > 0) then
      allocate(edist % next)
      call get_energy_dist(edist % next, LNW)
    end if

  end subroutine get_energy_dist

!===============================================================================
! LENGTH_ENERGY_DIST determines how many values are contained in an LDAT energy
! distribution array based on the secondary energy law and location in XSS
!===============================================================================

  function length_energy_dist(lc, law, LOCC, lid) result(length)

    integer, intent(in) :: lc     ! location in XSS array
    integer, intent(in) :: law    ! energy distribution law
    integer, intent(in) :: LOCC   ! location of energy distribution
    integer, intent(in) :: lid    ! length of interpolation data
    integer             :: length ! length of energy distribution (LDAT)

    integer :: i     ! loop index for incoming energies
    integer :: j     ! loop index for outgoing energies
    integer :: k     ! dummy index in XSS
    integer :: NR    ! number of interpolation regions
    integer :: NE    ! number of incoming energies
    integer :: NP    ! number of points in outgoing energy distribution
    integer :: NMU   ! number of points in outgoing cosine distribution
    integer :: NRa   ! number of interpolation regions for Watt 'a'
    integer :: NEa   ! number of energies for Watt 'a'
    integer :: NRb   ! number of interpolation regions for Watt 'b'
    integer :: NEb   ! number of energies for Watt 'b'
    real(8), allocatable :: L(:)  ! locations of distributions for each Ein

    ! initialize length
    length = 0

    select case (law)
    case (1)
      ! Tabular equiprobable energy bins
      NR = int(XSS(lc + 1))
      NE = int(XSS(lc + 2 + 2*NR))
      NP = int(XSS(lc + 3 + 2*NR + NE))
      length = 3 + 2*NR + NE + 3*NP*NE

    case (2)
      ! Discrete photon energy
      length = 2

    case (3)
      ! Level scattering
      length = 2

    case (4)
      ! Continuous tabular distribution
      NR = int(XSS(lc + 1))
      NE = int(XSS(lc + 2 + 2*NR))
      ! Before progressing, check to see if data set uses L(I) values
      ! in a way inconsistent with the current form of the ACE Format Guide
      ! (MCNP5 Manual, Vol 3)
      allocate(L(NE))
      L = int(XSS(lc + 3 + 2*NR + NE: lc + 3 + 2*NR + 2*NE - 1))
      do i = 1,NE
        ! Now check to see if L(i) is equal to any other entries
        ! If so, then we must exit
        if (count(L == L(i)) > 1) then
          message = "Invalid usage of L(I) in ACE data; &
                    &Consider using more recent data set."
          call fatal_error()
        end if
      end do
      deallocate(L)
      ! Continue with finding data length
      length = length + 2 + 2*NR + 2*NE
      do i = 1,NE
        ! determine length
        NP = int(XSS(lc + length + 2))
        length = length + 2 + 3*NP

        ! adjust location for this block
        j = lc + 2 + 2*NR + NE + i
        XSS(j) = XSS(j) - LOCC - lid
      end do

    case (5)
      ! General evaporation spectrum
      NR = int(XSS(lc + 1))
      NE = int(XSS(lc + 2 + 2*NR))
      NP = int(XSS(lc + 3 + 2*NR + 2*NE))
      length = 3 + 2*NR + 2*NE + NP

    case (7)
      ! Maxwell fission spectrum
      NR = int(XSS(lc + 1))
      NE = int(XSS(lc + 2 + 2*NR))
      length = 3 + 2*NR + 2*NE

    case (9)
      ! Evaporation spectrum
      NR = int(XSS(lc + 1))
      NE = int(XSS(lc + 2 + 2*NR))
      length = 3 + 2*NR + 2*NE

    case (11)
      ! Watt spectrum
      NRa = int(XSS(lc + 1))
      NEa = int(XSS(lc + 2 + 2*NRa))
      NRb = int(XSS(lc + 3 + 2*(NRa+NEa)))
      NEb = int(XSS(lc + 4 + 2*(NRa+NEa+NRb)))
      length = 5 + 2*(NRa + NEa + NRb + NEb)

    case (44)
      ! Kalbach-Mann correlated scattering
      NR = int(XSS(lc + 1))
      NE = int(XSS(lc + 2 + 2*NR))
      ! Before progressing, check to see if data set uses L(I) values
      ! in a way inconsistent with the current form of the ACE Format Guide
      ! (MCNP5 Manual, Vol 3)
      allocate(L(NE))
      L = int(XSS(lc + 3 + 2*NR + NE: lc + 3 + 2*NR + 2*NE - 1))
      do i = 1,NE
        ! Now check to see if L(i) is equal to any other entries
        ! If so, then we must exit
        if (count(L == L(i)) > 1) then
          message = "Invalid usage of L(I) in ACE data; &
                    &Consider using more recent data set."
          call fatal_error()
        end if
      end do
      deallocate(L)
      ! Continue with finding data length
      length = length + 2 + 2*NR + 2*NE
      do i = 1,NE
        NP = int(XSS(lc + length + 2))
        length = length + 2 + 5*NP

        ! adjust location for this block
        j = lc + 2 + 2*NR + NE + i
        XSS(j) = XSS(j) - LOCC - lid
      end do

    case (61)
      ! Correlated energy and angle distribution
      NR = int(XSS(lc + 1))
      NE = int(XSS(lc + 2 + 2*NR))
      ! Before progressing, check to see if data set uses L(I) values
      ! in a way inconsistent with the current form of the ACE Format Guide
      ! (MCNP5 Manual, Vol 3)
      allocate(L(NE))
      L = int(XSS(lc + 3 + 2*NR + NE: lc + 3 + 2*NR + 2*NE - 1))
      do i = 1,NE
        ! Now check to see if L(i) is equal to any other entries
        ! If so, then we must exit
        if (count(L == L(i)) > 1) then
          message = "Invalid usage of L(I) in ACE data; &
                    &Consider using more recent data set."
          call fatal_error()
        end if
      end do
      deallocate(L)
      ! Continue with finding data length
      length = length + 2 + 2*NR + 2*NE
      do i = 1,NE
        ! outgoing energy distribution
        NP = int(XSS(lc + length + 2))

        ! adjust locators for angular distribution
        do j = 1, NP
          k = lc + length + 2 + 3*NP + j
          if (XSS(k) /= 0) XSS(k) = XSS(k) - LOCC - lid
        end do

        length = length + 2 + 4*NP
        do j = 1, NP
          ! outgoing angle distribution -- NMU here is actually
          ! referred to as NP in the MCNP documentation
          NMU = int(XSS(lc + length + 2))
          length = length + 2 + 3*NMU
        end do

        ! adjust locators for energy distribution
        j = lc + 2 + 2*NR + NE + i
        XSS(j) = XSS(j) - LOCC - lid
      end do

    case (66)
      ! N-body phase space distribution
      length = 2

    case (67)
      ! Laboratory energy-angle law
      NR  = int(XSS(lc + 1))
      NE  = int(XSS(lc + 2 + 2*NR))
      ! Before progressing, check to see if data set uses L(I) values
      ! in a way inconsistent with the current form of the ACE Format Guide
      ! (MCNP5 Manual, Vol 3)
      allocate(L(NE))
      L = int(XSS(lc + 3 + 2*NR + NE: lc + 3 + 2*NR + 2*NE - 1))
      do i = 1,NE
        ! Now check to see if L(i) is equal to any other entries
        ! If so, then we must exit
        if (count(L == L(i)) > 1) then
          message = "Invalid usage of L(I) in ACE data; &
                    &Consider using more recent data set."
          call fatal_error()
        end if
      end do
      deallocate(L)
      ! Continue with finding data length
      NMU = int(XSS(lc + 4 + 2*NR + 2*NE))
      length = 4 + 2*(NR + NE + NMU)

    end select

  end function length_energy_dist

!===============================================================================
! READ_UNR_RES reads in unresolved resonance probability tables if present.
!===============================================================================

  subroutine read_unr_res(nuc)

    type(Nuclide), pointer :: nuc

    integer :: JXS23 ! location of URR data
    integer :: lc    ! locator
    integer :: N     ! # of incident energies
    integer :: M     ! # of probabilities
    integer :: i     ! index over incoming energies
    integer :: j     ! index over values
    integer :: k     ! index over probabilities

    ! determine locator for URR data
    JXS23 = JXS(23)

    ! check if URR data is present
    if (JXS23 /= 0) then
      nuc % urr_present = .true.
      allocate(nuc % urr_data)
      lc = JXS23
    else
      nuc % urr_present = .false.
      return
    end if

    ! read parameters
    nuc % urr_data % n_energy = int(XSS(lc))
    nuc % urr_data % n_prob = int(XSS(lc + 1))
    nuc % urr_data % interp = int(XSS(lc + 2))
    nuc % urr_data % inelastic_flag = int(XSS(lc + 3))
    nuc % urr_data % absorption_flag = int(XSS(lc + 4))
    if (int(XSS(lc + 5)) == 0) then
      nuc % urr_data % multiply_smooth = .false.
    else
      nuc % urr_data % multiply_smooth = .true.
    end if

    ! if the inelastic competition flag indicates that the inelastic cross
    ! section should be determined from a normal reaction cross section, we need
    ! to set up a pointer to that reaction
    nuc % urr_inelastic = NONE
    if (nuc % urr_data % inelastic_flag > 0) then
      do i = 1, nuc % n_reaction
        if (nuc % reactions(i) % MT == nuc % urr_data % inelastic_flag) then
          nuc % urr_inelastic = i
        end if
      end do

      ! Abort if no corresponding inelastic reaction was found
      if (nuc % urr_inelastic == NONE) then
        message = "Could not find inelastic reaction specified on " &
             // "unresolved resonance probability table."
        call fatal_error()
      end if
    end if

    ! allocate incident energies and probability tables
    N = nuc % urr_data % n_energy
    M = nuc % urr_data % n_prob
    allocate(nuc % urr_data % energy(N))
    allocate(nuc % urr_data % prob(N,6,M))

    ! read incident energies
    XSS_index = lc + 6
    nuc % urr_data % energy = get_real(N)

    ! read probability tables
    do i = 1, N
      do j = 1, 6
        do k = 1, M
          nuc % urr_data % prob(i,j,k) = XSS(XSS_index)
          XSS_index = XSS_index + 1
        end do
      end do
    end do

    ! Check for negative values
    if (any(nuc % urr_data % prob < ZERO)) then
      message = "Negative value(s) found on probability table for nuclide " &
           // nuc % name
      call warning()
    end if

  end subroutine read_unr_res

!===============================================================================
! GENERATE_NU_FISSION precalculates the microscopic nu-fission cross section for
! a given nuclide. This is done so that the nu_total function does not need to
! be called during cross section lookups.
!===============================================================================

  subroutine generate_nu_fission(nuc)

    type(Nuclide), pointer :: nuc

    integer :: i  ! index on nuclide energy grid
    real(8) :: E  ! energy
    real(8) :: nu ! # of neutrons per fission

    do i = 1, nuc % n_grid
      ! determine energy
      E = nuc % energy(i)

      ! determine total nu at given energy
      nu = nu_total(nuc, E)

      ! determine nu-fission microscopic cross section
      nuc % nu_fission(i) = nu * nuc % fission(i)
    end do

  end subroutine generate_nu_fission

!===============================================================================
! READ_THERMAL_DATA reads elastic and inelastic cross sections and corresponding
! secondary energy/angle distributions derived from experimental S(a,b)
! data. Namely, this routine reads the ITIE, ITCE, ITXE, and ITCA blocks.
!===============================================================================

  subroutine read_thermal_data(table)

    type(SAlphaBeta), pointer :: table

    integer :: i      ! index for incoming energies
    integer :: j      ! index for outgoing energies
    integer :: k      ! index for outoging angles
    integer :: lc     ! location in XSS array
    integer :: NE_in  ! number of incoming energies
    integer :: NE_out ! number of outgoing energies
    integer :: NMU    ! number of outgoing angles
    integer :: JXS4   ! location of elastic energy table
    integer(8), allocatable :: LOCC(:) ! Location of inelastic data

    ! read secondary energy mode for inelastic scattering
    table % secondary_mode = NXS(7)

    ! read number of inelastic energies and allocate arrays
    NE_in = int(XSS(JXS(1)))
    table % n_inelastic_e_in = NE_in
    allocate(table % inelastic_e_in(NE_in))
    allocate(table % inelastic_sigma(NE_in))

    ! read inelastic energies and cross-sections
    XSS_index = JXS(1) + 1
    table % inelastic_e_in = get_real(NE_in)
    table % inelastic_sigma = get_real(NE_in)

    ! set threshold value
    table % threshold_inelastic = table % inelastic_e_in(NE_in)

    ! allocate space for outgoing energy/angle for inelastic
    ! scattering
    if (table % secondary_mode == SAB_SECONDARY_EQUAL .or. &
         table % secondary_mode == SAB_SECONDARY_SKEWED) then
      NMU = NXS(3) + 1
      table % n_inelastic_mu = NMU
      NE_out = NXS(4)
      table % n_inelastic_e_out = NE_out
      allocate(table % inelastic_e_out(NE_out, NE_in))
      allocate(table % inelastic_mu(NMU, NE_out, NE_in))
    else if (table % secondary_mode == SAB_SECONDARY_CONT) then
      NMU = NXS(3) - 1
      table % n_inelastic_mu = NMU
      allocate(table % inelastic_data(NE_in))
      allocate(LOCC(NE_in))
      ! NE_out will be determined later
    end if

    ! read outgoing energy/angle distribution for inelastic scattering
    if (table % secondary_mode == SAB_SECONDARY_EQUAL .or. &
         table % secondary_mode == SAB_SECONDARY_SKEWED) then
      lc = JXS(3) - 1
      do i = 1, NE_in
        do j = 1, NE_out
          ! read outgoing energy
          table % inelastic_e_out(j,i) = XSS(lc + 1)

          ! read outgoing angles for this outgoing energy
          do k = 1, NMU
            table % inelastic_mu(k,j,i) = XSS(lc + 1 + k)
          end do

          ! advance pointer
          lc = lc + 1 + NMU
        end do
      end do
    else if (table % secondary_mode == SAB_SECONDARY_CONT) then
      ! Get the location pointers to each Ein's DistEnergySAB data
      LOCC = get_int(NE_in)
      ! Get the number of outgoing energies and allocate space accordingly
      do i = 1, NE_in
        NE_out = int(XSS(XSS_index + i - 1))
        table % inelastic_data(i) % n_e_out = NE_out
        allocate(table % inelastic_data(i) % e_out (NE_out))
        allocate(table % inelastic_data(i) % e_out_pdf (NE_out))
        allocate(table % inelastic_data(i) % e_out_cdf (NE_out))
        allocate(table % inelastic_data(i) % mu (NMU, NE_out))
      end do

      ! Now we can fill the inelastic_data(i) attributes
      do i = 1, NE_in
        XSS_index = LOCC(i)
        NE_out = table % inelastic_data(i) % n_e_out
        do j = 1, NE_out
          table % inelastic_data(i) % e_out(j) = XSS(XSS_index + 1)
          table % inelastic_data(i) % e_out_pdf(j) = XSS(XSS_index + 2)
          table % inelastic_data(i) % e_out_cdf(j) = XSS(XSS_index + 3)
          table % inelastic_data(i) % mu(:, j) = &
            XSS(XSS_index + 4: XSS_index + 4 + NMU - 1)
          XSS_index = XSS_index + 4 + NMU - 1
        end do
      end do
    end if

    ! read number of elastic energies and allocate arrays
    JXS4 = JXS(4)
    if (JXS4 /= 0) then
      NE_in = int(XSS(JXS4))
      table % n_elastic_e_in = NE_in
      allocate(table % elastic_e_in(NE_in))
      allocate(table % elastic_P(NE_in))

      ! read elastic energies and P
      XSS_index = JXS4 + 1
      table % elastic_e_in = get_real(NE_in)
      table % elastic_P    = get_real(NE_in)

      ! set threshold
      table % threshold_elastic = table % elastic_e_in(NE_in)

      ! determine whether sigma=P or sigma = P/E
      table % elastic_mode = NXS(5)
    else
      table % threshold_elastic = ZERO
      table % n_elastic_e_in = 0
    end if

    ! allocate space for outgoing energy/angle for elastic scattering
    NMU = NXS(6) + 1
    table % n_elastic_mu = NMU
    if (NMU > 0) then
      allocate(table % elastic_mu(NMU, NE_in))
    end if

    ! read equiprobable outgoing cosines for elastic scattering each
    ! incoming energy
    if (JXS4 /= 0 .and. NMU /= 0) then
      lc = JXS(6) - 1
      do i = 1, NE_in
        do j = 1, NMU
          table % elastic_mu(j,i) = XSS(lc + j)
        end do
        lc = lc + NMU
      end do
    end if

  end subroutine read_thermal_data

!===============================================================================
! GET_INT returns an array of integers read from the current position in the XSS
! array
!===============================================================================

  function get_int(n_values) result(array)

    integer, intent(in) :: n_values        ! number of values to read
    integer             :: array(n_values) ! array of values

    array = int(XSS(XSS_index:XSS_index + n_values - 1))
    XSS_index = XSS_index + n_values

  end function get_int

!===============================================================================
! GET_REAL returns an array of real(8)s read from the current position in the
! XSS array
!===============================================================================

  function get_real(n_values) result(array)

    integer, intent(in) :: n_values        ! number of values to read
    real(8)             :: array(n_values) ! array of values

    array = XSS(XSS_index:XSS_index + n_values - 1)
    XSS_index = XSS_index + n_values

  end function get_real

end module ace
