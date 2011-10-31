module cross_section

  use constants
  use cross_section_header, only: Nuclide, Reaction, SAB_Table, xsData
  use datatypes,            only: dict_create, dict_add_key, dict_get_key, &
                                  dict_has_key, dict_delete
  use datatypes_header,     only: DictionaryCI
  use endf,                 only: reaction_name
  use error,                only: fatal_error
  use fileio,               only: read_line, read_data, skip_lines
  use global
  use material_header,      only: Material
  use output,               only: write_message
  use string,               only: split_string, str_to_int, str_to_real, &
                                  lower_case, int_to_str

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

    integer        :: i                ! index in materials array
    integer        :: j                ! index over nuclides in material
    integer        :: index            ! index in xsdatas array
    integer        :: index_nuclides   ! index in nuclides
    integer        :: index_sab        ! index in sab_tables
    character(10)  :: key              ! name of isotope, e.g. 92235.03c
    type(Material),     pointer :: mat => null()
    type(Nuclide),      pointer :: nuc => null()
    type(SAB_Table),    pointer :: sab => null()
    type(DictionaryCI), pointer :: temp_dict => null()

    ! First we need to determine how many continuous-energy tables and how many
    ! S(a,b) thermal scattering tables there are -- this loop doesn't actually
    ! read the data, it simply counts the number of nuclides and S(a,b) tables
    index_nuclides = 0
    index_sab = 0
    do i = 1, n_materials
       ! Get pointer to material
       mat => materials(i)

       ! First go through all the nuclide and check if they exist on other
       ! materials -- if not, then increment the count for the number of
       ! nuclides and add to dictionary
       do j = 1, mat % n_nuclides
          key = mat % names(j)
          call lower_case(key)
          if (.not. dict_has_key(nuclide_dict, key)) then
             index_nuclides = index_nuclides + 1
             call dict_add_key(nuclide_dict, key, index_nuclides)
             mat % nuclide(j) = index_nuclides
          else
             mat % nuclide(j) = dict_get_key(nuclide_dict, key)
          end if
       end do

       ! Check if S(a,b) exists on other materials
       if (mat % has_sab_table) then
          key = mat % sab_name
          call lower_case(key)
          if (.not. dict_has_key(sab_dict, key)) then
             index_sab = index_sab + 1
             call dict_add_key(sab_dict, key, index_sab)
             mat % sab_table = index_sab
          else
             mat % sab_table = dict_get_key(sab_dict, key)
          end if
       else
          mat % sab_table = 0
       end if
    end do

    n_nuclides_total = index_nuclides
    n_sab_tables     = index_sab

    ! allocate arrays for ACE table storage
    allocate(nuclides(n_nuclides_total))
    allocate(sab_tables(n_sab_tables))

    ! allocate array for microscopic cross section cache
    allocate(micro_xs(n_nuclides_total))

    call dict_create(temp_dict)

    ! Now that the nuclides and sab_tables arrays have been allocated, we can
    ! repeat the same loop through each table in each material and read the
    ! cross sections
    index_nuclides = 0
    index_sab = 0
    do i = 1, n_materials
       mat => materials(i)
       do j = 1, mat % n_nuclides
          ! Get index in xsdatas array for this nuclide
          index = mat % xsdata(j)

          ! Get name of nuclide
          key = mat % names(j)
          call lower_case(key)
          
          if (.not. dict_has_key(temp_dict, key)) then
             index_nuclides = index_nuclides + 1
             call read_ACE_continuous(index_nuclides, index)
             call dict_add_key(temp_dict, key, index_nuclides)
          end if
       end do

       ! Read S(a,b) table if this material has one and it hasn't been read
       ! already
       if (mat % has_sab_table) then
          ! Get name of S(a,b) table
          key = mat % sab_name

          if (.not. dict_has_key(temp_dict, key)) then
             ! Find the entry in xsdatas for this table
             if (dict_has_key(xsdata_dict, key)) then
                index = dict_get_key(xsdata_dict, key)
             else
                message = "Cannot find cross-section " // trim(key) // &
                     " in specified xsdata file."
                call fatal_error()
             end if

             ! Read the table and add entry to dictionary
             index_sab = index_sab + 1
             call read_ACE_thermal(index_sab, index)
             call dict_add_key(temp_dict, key, index_sab)
          end if

          ! In order to know which nuclide the S(a,b) table applies to, we need
          ! to search through the list of nuclides for one which has a matching
          ! zaid
          sab => sab_tables(mat % sab_table)

          do j = 1, mat % n_nuclides
             nuc => nuclides(mat % nuclide(j))
             if (nuc % zaid == sab % zaid) then
                mat % sab_nuclide = j
             end if
          end do

          ! Check to make sure S(a,b) table matched a nuclide
          if (mat % sab_nuclide == 0) then
             message = "S(a,b) table " // trim(mat % sab_name) // " did not match " &
                  // "any nuclide on material " // trim(int_to_str(mat % uid))
             call fatal_error()
          end if
       end if
    end do

    ! delete dictionary
    call dict_delete(temp_dict)
       
  end subroutine read_xs

!===============================================================================
! READ_ACE_CONTINUOUS reads in a single ACE continuous-energy neutron cross
! section table. This routine reads the header data and then calls appropriate
! subroutines to parse the actual data.
!===============================================================================

  subroutine read_ACE_continuous(index_table, index)

    integer, intent(in) :: index_table ! index in nuclides array
    integer, intent(in) :: index       ! index in xsdatas array

    integer                 :: in = 7           ! unit to read from
    integer                 :: ioError          ! error status for file access
    integer                 :: words_per_line   ! number of words per line (data)
    integer                 :: lines            ! number of lines (data
    integer                 :: n                ! number of data values
    real(8)                 :: kT               ! ACE table temperature
    logical                 :: file_exists      ! does ACE library exist?
    logical                 :: found_xs         ! did we find table in library?
    character(7)            :: readable         ! is ACE library readable?
    character(MAX_LINE_LEN) :: line             ! single line to read
    character(MAX_WORD_LEN) :: words(MAX_WORDS) ! words on a line
    character(MAX_WORD_LEN) :: filename         ! name of ACE library file
    character(10)           :: tablename        ! name of cross section table
    type(Nuclide), pointer :: nuc => null()

    ! Check to make sure index in nuclides array and xsdata arrays are valid
    if (index_table > size(nuclides)) then
       message = "Index of table to read is greater than length of nuclides."
       call fatal_error()
    elseif (index > size(xsdatas)) then
       message = "Index of xsdata entry is greater than length of xsdatas."
       call fatal_error()
    end if

    filename = xsdatas(index)%path
    tablename = xsdatas(index)%id

    nuc => nuclides(index_table)

    ! Check if input file exists and is readable
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
    message = "Loading ACE cross section table: " // tablename
    call write_message(6)

    ! open file
    open(file=filename, unit=in, status='old', & 
         & action='read', iostat=ioError)
    if (ioError /= 0) then
       message = "Error while opening file: " // filename
       call fatal_error()
    end if

    found_xs = .false.
    do while (.not. found_xs)
       call read_line(in, line, ioError)
       if (ioError < 0) then
          message = "Could not find ACE table " // tablename // "."
          call fatal_error()
       end if
       call split_string(line, words, n)
       if (trim(words(1)) == trim(tablename)) then
          found_xs = .true.
          nuc % name = words(1)
          nuc % awr  = str_to_real(words(2))
          kT = str_to_real(words(3))
          nuc % temp = kT / K_BOLTZMANN
       end if
       
       ! Skip 5 lines
       call skip_lines(in, 5, ioError)

       ! Read NXS data
       lines = 2
       words_per_line = 8
       call read_data(in, NXS, 16, lines, words_per_line)

       ! Set ZAID of nuclide
       nuc % zaid = NXS(2)

       ! Read JXS data
       lines = 4
       call read_data(in, JXS, 32, lines, words_per_line)

       ! Calculate how many data points and lines in the XSS array
       n = NXS(1)
       lines = (n + 3)/4

       if (found_xs) then
          ! allocate storage for XSS array
          allocate(XSS(n))
          
          ! Read XSS
          words_per_line = 4
          call read_data(in, XSS, n, lines, words_per_line)
       else
          call skip_lines(in, lines, ioError)
       end if

    end do

    call read_esz(nuc)
    call read_nu_data(nuc)
    call read_reactions(nuc)
    call read_angular_dist(nuc)
    call read_energy_dist(nuc)
    call read_unr_res(nuc)

    ! Free memory from XSS array
    if(allocated(XSS)) deallocate(XSS)
    if(associated(nuc)) nullify(nuc)

    close(unit=in)

  end subroutine read_ACE_continuous

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
    allocate(nuc % absorption(NE))
    allocate(nuc % heating(NE))

    ! initialize cross sections
    nuc % total      = ZERO
    nuc % elastic    = ZERO
    nuc % fission    = ZERO
    nuc % absorption = ZERO
    nuc % heating    = ZERO

    ! Read data from XSS -- only the energy grid, elastic scattering and heating
    ! cross section values are actually read from here. The total and absorption
    ! cross sections are reconstructed from the partial reaction data.

    XSS_index = 1
    nuc % energy = get_real(NE)

    ! Skip total and absorption
    XSS_index = XSS_index + 2*NE

    ! Continue reading elastic scattering and heating
    nuc % elastic = get_real(NE)
    nuc % heating = get_real(NE)
    
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
    integer :: LNW    ! location of next energy distribution if multiple
    integer :: LAW    ! secondary energy distribution law   
    integer :: IDAT   ! location of first energy distribution for given MT
    integer :: loc    ! locator
    integer :: length             ! length of data to allocate
    integer :: length_interp_data ! length of interpolation data

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
       LNU = XSS(KNU)
       if (LNU == 1) then
          ! Polynomial data
          nuc % nu_p_type = NU_POLYNOMIAL

          ! allocate determine how many coefficients for polynomial
          NC = XSS(KNU+1)
          length = NC + 1
       elseif (LNU == 2) then
          ! Tabular data
          nuc % nu_p_type = NU_TABULAR

          ! determine number of interpolation regions and number of energies
          NR = XSS(KNU+1)
          NE = XSS(KNU+2+2*NR)
          length = 2 + 2*NR + 2*NE
       end if

       ! allocate space for nu data storage
       allocate(nuc % nu_p_data(length))

       ! read data
       XSS_index = KNU + 1
       nuc % nu_p_data = get_real(length)

       ! Now read total nu data
       KNU = JXS2 + abs(XSS(JXS2)) + 1
       LNU = XSS(KNU)
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
          LOCC = XSS(LED + i - 1)

          LNW  = XSS(LDIS + LOCC - 1)
          LAW  = XSS(LDIS + LOCC)
          IDAT = XSS(LDIS + LOCC + 1)
          NR   = XSS(LDIS + LOCC + 2)
          nuc % nu_d_edist(i) % law = LAW
          nuc % nu_d_edist(i) % n_interp = NR

          ! allocate space for ENDF interpolation parameters
          if (NR > 0) then
             allocate(nuc % nu_d_edist(i) % nbt(NR))
             allocate(nuc % nu_d_edist(i) % int(NR))
          end if

          ! read ENDF interpolation parameters
          XSS_index = LDIS + LOCC + 3
          if (NR > 0) then
             nuc % nu_d_edist(i) % nbt = get_real(NR)
             nuc % nu_d_edist(i) % int = get_real(NR)
          end if

          ! allocate space for law validity data
          NE = XSS(LDIS + LOCC + 3 + 2*NR)
          allocate(nuc % nu_d_edist(i) % energy(NE))
          allocate(nuc % nu_d_edist(i) % pvalid(NE))

          length_interp_data = 5 + 2*(NR + NE)

          ! read law validity data
          XSS_index = LDIS + LOCC + 4 + 2*NR
          nuc % nu_d_edist(i) % energy = get_real(NE)
          nuc % nu_d_edist(i) % pvalid = get_real(NE)

          ! Set index to beginning of IDAT array
          loc = LDIS + IDAT - 2

          ! determine length of energy distribution
          length = length_energy_dist(loc, LAW, LOCC, length_interp_data)

          ! allocate secondary energy distribution array
          allocate(nuc % nu_d_edist(i) % data(length))

          ! read secondary energy distribution
          XSS_index = loc + 1
          nuc % nu_d_edist(i) % data = get_real(length)
       end do

       ! =======================================================================
       ! DELAYED NEUTRON PRECUSOR YIELDS AND CONSTANTS

       ! determine length of all precursor constants/yields/interp data
       length = 0
       loc = JXS(25)
       do i = 1, NPCR
          NR = XSS(loc + length + 1)
          NE = XSS(loc + length + 2 + 2*NR)
          length = length + 3 + 2*NR + 2*NE
       end do

       ! allocate space for precusor data
       allocate(nuc % nu_d_precursor_data(length))

       ! read delayed neutron precursor data
       XSS_index = loc
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

    ! Store elastic scattering cross-section on reaction one
    rxn => nuc % reactions(1)
    rxn % MT      = 2
    rxn % Q_value = ZERO
    rxn % TY      = 1
    rxn % IE      = 1
    allocate(rxn % sigma(nuc % n_grid))
    rxn % sigma = nuc % elastic
    
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

       ! read MT number, Q-value, and neutrons produced
       rxn % MT      = XSS(LMT + i - 1)
       rxn % Q_value = XSS(JXS4 + i - 1)
       rxn % TY      = XSS(JXS5 + i - 1)

       ! read starting energy index
       LOCA = XSS(LXS + i - 1)
       IE = XSS(JXS7 + LOCA - 1)
       rxn % IE = IE

       ! read number of energies cross section values
       NE = XSS(JXS7 + LOCA)
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

          ! Keep track of this reaction for easy searching later
          i_fission = i_fission + 1
          nuc % index_fission(i_fission) = i + 1
          nuc % n_fission = nuc % n_fission + 1
       end if

       ! set defaults
       rxn % has_angle_dist  = .false.
       rxn % has_energy_dist = .false.
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
       LOCB = XSS(JXS8 + i - 1)
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
       NE = XSS(JXS9 + LOCB - 1)
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
             NP = XSS(JXS9 + abs(LC))
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
    integer :: LDIS  ! location of all energy distributions
    integer :: LOCC  ! location of energy distributions for given MT
    integer :: LNW   ! location of next energy distribution if multiple
    integer :: LAW   ! secondary energy distribution law   
    integer :: NR    ! number of interpolation regions
    integer :: NE    ! number of incoming energies
    integer :: IDAT  ! location of first energy distribution for given MT
    integer :: loc   ! locator
    integer :: length             ! length of data to allocate
    integer :: length_interp_data ! length of interpolation data
    integer :: i                  ! loop index
    type(Reaction), pointer :: rxn => null()

    LED  = JXS(10)
    LDIS = JXS(11)

    ! Loop over all reactions 
    do i = 1, NXS(5)
       rxn => nuc % reactions(i+1) ! skip over elastic scattering
       rxn % has_energy_dist = .true.

       ! find location of energy distribution data
       LOCC = XSS(LED + i - 1)
       
       LNW  = XSS(LDIS + LOCC - 1)
       LAW  = XSS(LDIS + LOCC)
       IDAT = XSS(LDIS + LOCC + 1)
       NR   = XSS(LDIS + LOCC + 2)
       rxn % edist % law = LAW
       rxn % edist % n_interp = NR

       ! allocate space for ENDF interpolation parameters
       if (NR > 0) then
          allocate(rxn % edist % nbt(NR))
          allocate(rxn % edist % int(NR))
       end if
       
       ! read ENDF interpolation parameters
       XSS_index = LDIS + LOCC + 3
       if (NR > 0) then
          rxn % edist % nbt = get_real(NR)
          rxn % edist % int = get_real(NR)
       end if

       ! allocate space for law validity data
       NE = XSS(LDIS + LOCC + 3 + 2*NR)
       allocate(rxn % edist % energy(NE))
       allocate(rxn % edist % pvalid(NE))
       
       length_interp_data = 5 + 2*(NR + NE)
       
       ! read law validity data
       XSS_index = LDIS + LOCC + 4 + 2*NR
       rxn % edist % energy = get_real(NE)
       rxn % edist % pvalid = get_real(NE)

       ! Set index to beginning of IDAT array
       loc = LDIS + IDAT - 2

       ! determine length of energy distribution
       length = length_energy_dist(loc, LAW, LOCC, length_interp_data)

       ! allocate secondary energy distribution array
       allocate(rxn % edist % data(length))

       ! read secondary energy distribution
       XSS_index = loc + 1
       rxn % edist % data = get_real(length)
    end do

  end subroutine read_energy_dist

!===============================================================================
! LENGTH_ENERGY_DIST determines how many values are contained in an LDAT energy
! distribution array based on the secondary energy law and location in XSS
!===============================================================================

  function length_energy_dist(loc, law, LOCC, lid) result(length)

    integer, intent(in) :: loc    ! location in XSS array
    integer, intent(in) :: law    ! energy distribution law
    integer, intent(in) :: LOCC   ! location of energy distribution
    integer, intent(in) :: lid    ! length of interpolation data
    integer             :: length ! length of energy distribution (LDAT)

    integer :: i,j,k
    integer :: NR    ! number of interpolation regions
    integer :: NE    ! number of incoming energies
    integer :: NP    ! number of points in outgoing energy distribution
    integer :: NMU   ! number of points in outgoing cosine distribution
    integer :: NRa   ! number of interpolation regions for Watt 'a'
    integer :: NEa   ! number of energies for Watt 'a'
    integer :: NRb   ! number of interpolation regions for Watt 'b'
    integer :: NEb   ! number of energies for Watt 'b'

    ! initialize length
    length = 0

    select case (law)
    case (1)
       ! Tabular equiprobable energy bins
       NR = XSS(loc + 1)
       NE = XSS(loc + 2 + 2*NR)
       NP = XSS(loc + 3 + 2*NR + NE)
       length = 3 + 2*NR + NE + 3*NP*NE

    case (2)
       ! Discrete photon energy
       length = 2

    case (3)
       ! Level scattering
       length = 2

    case (4)
       ! Continuous tabular distribution
       NR = XSS(loc + 1)
       NE = XSS(loc + 2 + 2*NR)
       length = length + 2 + 2*NR + 2*NE
       do i = 1,NE
          ! determine length
          NP = XSS(loc + length + 2)
          length = length + 2 + 3*NP

          ! adjust location for this block
          j = loc + 2 + 2*NR + NE + i
          XSS(j) = XSS(j) - LOCC - lid
       end do

    case (5)
       ! General evaporation spectrum
       NR = XSS(loc + 1)
       NE = XSS(loc + 2 + 2*NR)
       NP = XSS(loc + 3 + 2*NR + 2*NE)
       length = 3 + 2*NR + 2*NE + NP

    case (7)
       ! Maxwell fission spectrum
       NR = XSS(loc + 1)
       NE = XSS(loc + 2 + 2*NR)
       length = 3 + 2*NR + 2*NE

    case (9)
       ! Evaporation spectrum
       NR = XSS(loc + 1)
       NE = XSS(loc + 2 + 2*NR)
       length = 3 + 2*NR + 2*NE

    case (11)
       ! Watt spectrum
       NRa = XSS(loc + 1)
       NEa = XSS(loc + 2 + 2*NRa)
       NRb = XSS(loc + 3 + 2*(NRa+NEa))
       NEb = XSS(loc + 4 + 2*(NRa+NEa+NRb))
       length = 5 + 2*(NRa + NEa + NRb + NEb)

    case (44)
       ! Kalbach-Mann correlated scattering
       NR = XSS(loc + 1)
       NE = XSS(loc + 2 + 2*NR)
       length = length + 2 + 2*NR + 2*NE
       do i = 1,NE
          NP = XSS(loc + length + 2)
          length = length + 2 + 5*NP

          ! adjust location for this block
          j = loc + 2 + 2*NR + NE + i
          XSS(j) = XSS(j) - LOCC - lid
       end do

    case (61)
       ! Correlated energy and angle distribution
       NR = XSS(loc + 1)
       NE = XSS(loc + 2 + 2*NR)
       length = length + 2 + 2*NR + 2*NE
       do i = 1,NE
          ! outgoing energy distribution
          NP = XSS(loc + length + 2)

          ! adjust locators for angular distribution
          do j = 1, NP
             k = loc + length + 2 + 3*NP + j
             if (XSS(k) /= 0) XSS(k) = XSS(k) - LOCC - lid
          end do

          length = length + 2 + 4*NP
          do j = 1, NP
             ! outgoing angle distribution -- NMU here is actually
             ! referred to as NP in the MCNP documentation
             NMU = XSS(loc + length + 2)
             length = length + 2 + 3*NMU
          end do

          ! adjust locators for energy distribution
          j = loc + 2 + 2*NR + NE + i
          XSS(j) = XSS(j) - LOCC - lid
       end do

    case (66)
       ! N-body phase space distribution
       length = 2

    case (67)
       ! Laboratory energy-angle law
       NR = XSS(loc + 1)
       NE = XSS(loc + 2 + 2*NR)
       NMU = XSS(loc + 4 + 2*NR + 2*NE)
       length = 4 + 2*(NR + NE + NMU)

    end select

  end function length_energy_dist

!===============================================================================
! READ_UNR_RES reads in unresolved resonance probability tables if present.
!===============================================================================

  subroutine read_unr_res(nuc)

    type(Nuclide), pointer :: nuc

    integer :: JXS23 ! location of URR data
    integer :: loc   ! locator
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
       allocate(nuc % urr_data % params(6))
       loc = JXS23
    else
       nuc % urr_present = .false.
       return
    end if

    ! read parameters
    nuc % urr_data % params(1) = XSS(loc)     ! # of incident energies
    nuc % urr_data % params(2) = XSS(loc + 1) ! # of probabilities
    nuc % urr_data % params(3) = XSS(loc + 2) ! interpolation parameter
    nuc % urr_data % params(4) = XSS(loc + 3) ! inelastic competition flag
    nuc % urr_data % params(5) = XSS(loc + 4) ! other absorption flag
    nuc % urr_data % params(6) = XSS(loc + 5) ! factors flag

    ! allocate incident energies and probability tables
    N = nuc % urr_data % params(1)
    M = nuc % urr_data % params(2)
    allocate(nuc % urr_data % energy(N))
    allocate(nuc % urr_data % prob(N,6,M))

    ! read incident energies
    XSS_index = loc + 6
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

  end subroutine read_unr_res

!===============================================================================
! READ_ACE_THERMAL reads in a single ACE thermal scattering S(a,b) table. This
! routine in particular reads the header data and then calls read_thermal_data
! to parse the actual data blocks.
!===============================================================================

  subroutine read_ACE_thermal(index_table, index)

    integer, intent(in) :: index_table ! index in sab_tables array
    integer, intent(in) :: index       ! index in xsdatas array

    integer                 :: in = 7           ! unit to read from
    integer                 :: ioError          ! error status for file access
    integer                 :: words_per_line   ! number of words per line (data)
    integer                 :: lines            ! number of lines (data
    integer                 :: n                ! number of data values
    real(8)                 :: kT               ! ACE table temperature
    logical                 :: file_exists      ! does ACE library exist?
    logical                 :: found_xs         ! did we find table in library?
    character(7)            :: readable         ! is ACE library readable?
    character(MAX_LINE_LEN) :: line             ! single line to read
    character(MAX_WORD_LEN) :: words(MAX_WORDS) ! words on a line
    character(MAX_WORD_LEN) :: filename         ! name of ACE library file
    character(10)           :: tablename        ! name of cross section table
    type(SAB_Table), pointer :: table => null()

    filename = xsdatas(index)%path
    tablename = xsdatas(index)%id

    table => sab_tables(index_table)

    ! Check if input file exists and is readable
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
    message = "Loading ACE cross section table: " // tablename
    call write_message(6)

    ! open file
    open(file=filename, unit=in, status='old', & 
         & action='read', iostat=ioError)
    if (ioError /= 0) then
       message = "Error while opening file: " // filename
       call fatal_error()
    end if

    found_xs = .false.
    do while (.not. found_xs)
       call read_line(in, line, ioError)
       if (ioError < 0) then
          message = "Could not find ACE table " // tablename // "."
          call fatal_error()
       end if
       call split_string(line, words, n)
       if (trim(words(1)) == trim(tablename)) then
          found_xs = .true.
          table%name = words(1)
          table%awr = str_to_real(words(2))
          kT = str_to_real(words(3))
          table%temp = kT / K_BOLTZMANN

          ! Skip 1 line
          call skip_lines(in, 1, ioError)

          ! Read corresponding ZAID
          call read_line(in, line, ioError)
          call split_string(line, words, n)
          table % zaid = str_to_int(words(1))

          ! Skip remaining 3 lines
          call skip_lines(in, 3, ioError)
       else
          ! Skip 5 lines
          call skip_lines(in, 5, ioError)
       end if

       ! Read NXS data
       lines = 2
       words_per_line = 8
       call read_data(in, NXS, 16, lines, words_per_line)

       ! Read JXS data
       lines = 4
       call read_data(in, JXS, 32, lines, words_per_line)

       ! Calculate how many data points and lines in the XSS array
       n = NXS(1)
       lines = (n + 3)/4

       if (found_xs) then
          ! allocate storage for XSS array
          allocate(XSS(n))
          
          ! Read XSS
          words_per_line = 4
          call read_data(in, XSS, n, lines, words_per_line)
       else
          call skip_lines(in, lines, ioError)
       end if

    end do

    call read_thermal_data(table)

    ! Free memory from XSS array
    if(allocated(XSS)) deallocate(XSS)
    if(associated(table)) nullify(table)

    close(unit=in)

  end subroutine read_ACE_thermal

!===============================================================================
! READ_THERMAL_DATA reads elastic and inelastic cross sections and corresponding
! secondary energy/angle distributions derived from experimental S(a,b)
! data. Namely, in MCNP language, this routine reads the ITIE, ITCE, ITXE, and
! ITCA blocks.
!===============================================================================

  subroutine read_thermal_data(table)

    type(SAB_Table), pointer :: table

    integer :: i      ! index for incoming energies
    integer :: j      ! index for outgoing energies
    integer :: k      ! index for outoging angles
    integer :: loc    ! location in XSS array
    integer :: NE_in  ! number of incoming energies
    integer :: NE_out ! number of outgoing energies
    integer :: NMU    ! number of outgoing angles
    integer :: JXS4   ! location of elastic energy table

    ! read secondary energy mode for inelastic scattering
    table % secondary_mode = NXS(7)

    ! read number of inelastic energies and allocate arrays
    NE_in = XSS(JXS(1))
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
    NE_out = NXS(4)
    NMU = NXS(3) + 1
    table % n_inelastic_e_out = NE_out
    table % n_inelastic_mu = NMU
    allocate(table % inelastic_e_out(NE_out, NE_in))
    allocate(table % inelastic_mu(NMU, NE_out, NE_in))

    ! read outgoing energy/angle distribution for inelastic scattering
    loc = JXS(3) - 1
    do i = 1, NE_in
       do j = 1, NE_out
          ! read outgoing energy
          table % inelastic_e_out(j,i) = XSS(loc + 1)

          ! read outgoing angles for this outgoing energy
          do k = 1, NMU
             table % inelastic_mu(k,j,i) = XSS(loc + 1 + k)
          end do

          ! advance pointer
          loc = loc + 1 + NMU
       end do
    end do

    ! read number of elastic energies and allocate arrays
    JXS4 = JXS(4)
    if (JXS4 /= 0) then
       NE_in = XSS(JXS4)
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
       loc = JXS(6) - 1
       do i = 1, NE_in
          do j = 1, NMU
             table % elastic_mu(j,i) = XSS(loc + j)
          end do
          loc = loc + NMU
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

!===============================================================================
! GET_MACRO_XS
!===============================================================================

  function get_macro_xs(p, mat, MT) result(xs)

    type(Particle), pointer    :: p
    type(Material), pointer    :: mat
    integer,        intent(in) :: MT
    real(8)                    :: xs

    integer :: i, j
    integer :: n_nuclides
    integer :: IE
    real(8) :: density_i
    real(8) :: sigma_i
    real(8) :: f
    type(Nuclide),  pointer :: nuc => null()
    type(Reaction), pointer :: rxn => null()

    ! initialize xs
    xs = ZERO

    ! loop over all nuclides in material
    n_nuclides = mat % n_nuclides
    do i = 1, n_nuclides
       nuc => nuclides(mat % nuclide(i))

       ! determine nuclide atom density
       density_i = mat % atom_density(i)

       ! search nuclide energy grid
       IE = nuc%grid_index(p % IE)
       f = (p%E - nuc%energy(IE))/(nuc%energy(IE+1) - nuc%energy(IE))
       
       ! handle special case of total cross section
       if (MT == 1) then
          xs = xs + mat % density * (ONE-f) * nuc%total(IE) + & 
               & f * (nuc%total(IE+1))
          cycle
       end if

       ! loop over reactions in isotope
       do j = 1, nuc % n_reaction
          rxn => nuc % reactions(i)

          ! check for matching MT
          if (MT /= rxn % MT) cycle

          ! if energy is below threshold for this reaction, skip it
          if (IE < rxn % IE) cycle

          ! add to cumulative probability
          sigma_i = (ONE-f) * rxn%sigma(IE-rxn%IE+1) + & 
               & f * (rxn%sigma(IE-rxn%IE+2))
       end do

       ! calculate nuclide macroscopic cross-section
       xs = xs + density_i * sigma_i
    end do

  end function get_macro_xs

!===============================================================================
! READ_XSDATA reads the data in a SERPENT xsdata file and builds a dictionary to
! find cross-section information later on.
!===============================================================================

  subroutine read_xsdata(path)

    character(*), intent(in) :: path

    type(xsData), pointer    :: iso => null()
    character(MAX_LINE_LEN)  :: line
    character(MAX_WORD_LEN)  :: words(MAX_WORDS)
    character(MAX_WORD_LEN)  :: filename
    integer                  :: n
    integer                  :: in = 7
    logical                  :: file_exists
    character(7)             :: readable
    integer                  :: count
    integer                  :: index
    integer                  :: ioError

    message = "Reading cross-section summary file..."
    call write_message(5)

    ! Construct filename
    filename = trim(path)

    ! Check if xsdata exists and is readable
    inquire(FILE=filename, EXIST=file_exists, READ=readable)
    if (.not. file_exists) then
       message = "Cross section summary '" // trim(filename) // "' does not exist!"
       call fatal_error()
    elseif (readable(1:3) == 'NO') then
       message = "Cross section summary '" // trim(filename) // "' is not " &
            & // "readable! Change file permissions with chmod command."
       call fatal_error()
    end if

    ! open xsdata file
    open(FILE=filename, UNIT=in, STATUS='old', &
         & ACTION='read', IOSTAT=ioError)
    if (ioError /= 0) then
       message = "Error while opening file: " // filename
       call fatal_error()
    end if

    ! determine how many lines
    count = 0
    do
       read(UNIT=in, FMT='(A)', IOSTAT=ioError) line
       if (ioError < 0) then
          ! reached end of file
          exit
       elseif (ioError > 0) then
          message = "Unknown error while reading file: " // filename
          close(UNIT=in)
          call fatal_error()
       end if
       count = count + 1
    end do
    allocate(xsdatas(count))

    ! read actual lines
    index = 0
    rewind(in)
    do
       read(UNIT=in, FMT='(A)', IOSTAT=ioError) line
       if (ioError < 0) exit
       index = index + 1
       call split_string(line, words, n)
       if (n == 0) cycle ! skip blank line

       ! Check to make sure there are enough arguments
       if (n < 9) then
          message = "Not enough arguments on xsdata line: " // line
          close(UNIT=in)
          call fatal_error()
       end if

       iso => xsdatas(index)

       ! store data
       iso%alias = words(1)
       iso%id = words(2)
       iso%type = str_to_int(words(3))
       iso%zaid = str_to_int(words(4))
       iso%isomeric = str_to_int(words(5))
       iso%awr = str_to_real(words(6))
       iso%temp = str_to_real(words(7))
       iso%binary = str_to_int(words(8))
       iso%path = words(9)

       ! create dictionary entry
       call dict_add_key(xsdata_dict, iso%alias, index)
    end do

    close(UNIT=in)

  end subroutine read_xsdata

end module cross_section
