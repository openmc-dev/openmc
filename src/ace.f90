module ace

  use global
  use error,           only: fatal_error
  use output,          only: message
  use string,          only: lower_case
  use fileio,          only: read_line, read_data, skip_lines
  use string,          only: split_string, str_to_real
  use data_structures, only: dict_create, dict_add_key, dict_has_key, &
       &                     dict_get_key, dict_delete
  use endf,            only: reaction_name

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
! xs_continuous and xs_thermal arrays
!===============================================================================

  subroutine read_xs()

    integer        :: i                ! index in materials array
    integer        :: j                ! index over isotopes in material
    integer        :: index            ! index in xsdatas array
    integer        :: n                ! length of key
    integer        :: index_continuous ! index in xs_continuous
    integer        :: index_thermal    ! index in xs_thermal
    character(10)  :: key              ! name of isotope, e.g. 92235.03c
    character(max_line_len) :: msg              ! output/error message
    type(Material),      pointer :: mat => null()
    type(xsData),        pointer :: iso => null()
    type(AceContinuous), pointer :: ace_cont => null()
    type(AceThermal),    pointer :: ace_thermal => null()
    type(DictionaryCI),  pointer :: temp_dict => null()

    call dict_create(ace_dict)

    ! determine how many continuous-energy tables and how many S(a,b) thermal
    ! scattering tables there are
    index_continuous = 0
    index_thermal = 0
    do i = 1, n_materials
       mat => materials(i)
       do j = 1, mat%n_isotopes
          index = mat%isotopes(j)
          key = xsdatas(index)%id
          n = len_trim(key)
          call lower_case(key)
          select case (key(n:n))
          case ('c')
             if (.not. dict_has_key(ace_dict, key)) then
                index_continuous = index_continuous + 1
                call dict_add_key(ace_dict, key, index_continuous)
                mat%table(j) = index_continuous
             else
                mat%table(j) = dict_get_key(ace_dict, key)
             end if
          case ('t')
             ! Need same logic as for continuous tables
             n_thermal = n_thermal + 1
          case default
             msg = "Unknown cross section table type: " // key
             call fatal_error(msg)
          end select
       end do
    end do

    n_continuous = index_continuous
    n_thermal = index_thermal

    ! allocate arrays for ACE table storage
    allocate(xs_continuous(n_continuous))
    allocate(xs_thermal(n_thermal))

    ! loop over all nuclides in xsdata
    call dict_create(temp_dict)

    index_continuous = 0
    index_thermal = 0
    do i = 1, n_materials
       mat => materials(i)
       do j = 1, mat%n_isotopes
          index = mat%isotopes(j)
          key = xsdatas(index)%id
          n = len_trim(key)
          call lower_case(key)
          select case (key(n:n))
          case ('c')
             if (.not. dict_has_key(temp_dict, key)) then
                index_continuous = index_continuous + 1
                call read_ACE_continuous(index_continuous, index)
             end if
          case ('t')
             n_thermal = n_thermal + 1
          end select
       end do
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

    integer, intent(in) :: index_table ! index in xs_continuous array
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
    character(max_line_len) :: msg              ! output/error message
    character(max_line_len) :: line             ! single line to read
    character(max_word_len) :: words(max_words) ! words on a line
    character(max_word_len) :: filename         ! name of ACE library file
    character(10)           :: tablename        ! name of cross section table
    type(AceContinuous), pointer :: table => null()

    filename = xsdatas(index)%path
    tablename = xsdatas(index)%id

    table => xs_continuous(index_table)

    ! Check if input file exists and is readable
    inquire(FILE=filename, EXIST=file_exists, READ=readable)
    if (.not. file_exists) then
       msg = "ACE library '" // trim(filename) // "' does not exist!"
       call fatal_error(msg)
    elseif (readable(1:3) == 'NO') then
       msg = "ACE library '" // trim(filename) // "' is not readable! &
            &Change file permissions with chmod command."
       call fatal_error(msg)
    end if

    ! display message
    msg = "Loading ACE cross section table: " // tablename
    call message(msg, 6)

    ! open file
    open(file=filename, unit=in, status='old', & 
         & action='read', iostat=ioError)
    if (ioError /= 0) then
       msg = "Error while opening file: " // filename
       call fatal_error(msg)
    end if

    found_xs = .false.
    do while (.not. found_xs)
       call read_line(in, line, ioError)
       if (ioError < 0) then
          msg = "Could not find ACE table " // tablename // "."
          call fatal_error(msg)
       end if
       call split_string(line, words, n)
       if (trim(words(1)) == trim(tablename)) then
          found_xs = .true.
          table%name = words(1)
          table%awr = str_to_real(words(2))
          kT = str_to_real(words(3))
          table%temp = kT / K_BOLTZMANN
       end if
       
       ! Skip 5 lines
       call skip_lines(in, 5, ioError)

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

    call read_esz(table)
    call read_nu_data(table)
    call read_reactions(table)
    call read_angular_dist(table)
    call read_energy_dist(table)
    call read_unr_res(table)

    ! Free memory from XSS array
    if(allocated(XSS)) deallocate(XSS)
    if(associated(table)) nullify(table)

    close(unit=in)

  end subroutine read_ACE_continuous

!===============================================================================
! READ_ESZ - reads through the ESZ block. This block contains the energy grid,
! total xs, absorption xs, elastic scattering xs, and heating numbers.
!===============================================================================

  subroutine read_esz(table)

    type(AceContinuous), pointer :: table

    integer :: NE ! number of energy points for total and elastic cross sections

    ! determine number of energy points
    NE = NXS(3)
    table%n_grid = NE

    ! allocate storage for arrays
    allocate(table%energy(NE))
    allocate(table%sigma_t(NE))
    allocate(table%sigma_a(NE))
    allocate(table%sigma_el(NE))
    allocate(table%heating(NE))

    ! read data from XSS -- right now the total, absorption and elastic
    ! scattering are read in to these special arrays, but in reality, it should
    ! be necessary to only store elastic scattering and possibly total
    ! cross-section for total material xs generation.
    XSS_index = 1
    table%energy = get_real(NE)
    table%sigma_t = get_real(NE)
    table%sigma_a = get_real(NE)
    table%sigma_el = get_real(NE)
    table%heating = get_real(NE)
    
  end subroutine read_esz

!===============================================================================
! READ_NU_DATA reads data given on the number of neutrons emitted from fission
! as a function of the incoming energy of a neutron. This data may be broken
! down into prompt and delayed neutrons emitted as well.
!===============================================================================

  subroutine read_nu_data(table)

    type(AceContinuous), pointer :: table

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
       table % nu_t_type = NU_NONE
       table % nu_p_type = NU_NONE

    elseif (XSS(JXS2) > 0) then
       ! =======================================================================
       ! PROMPT OR TOTAL NU DATA
       KNU = JXS2
       LNU = int(XSS(KNU))
       if (LNU == 1) then
          ! Polynomial data
          table % nu_t_type = NU_POLYNOMIAL
          table % nu_p_type = NU_NONE

          ! allocate determine how many coefficients for polynomial
          NC = int(XSS(KNU+1))
          length = NC + 1
       elseif (LNU == 2) then
          ! Tabular data
          table % nu_t_type = NU_TABULAR
          table % nu_p_type = NU_NONE

          ! determine number of interpolation regions and number of energies
          NR = int(XSS(KNU+1))
          NE = int(XSS(KNU+2+2*NR))
          length = 2 + 2*NR + 2*NE
       end if

       ! allocate space for nu data storage
       allocate(table % nu_t_data(length))

       ! read data -- for polynomial, this is the number of coefficients and the
       ! coefficients themselves, and for tabular, this is interpolation data
       ! and tabular E/nu
       XSS_index = KNU + 1
       table % nu_t_data = get_real(length)

    elseif (XSS(JXS2) < 0) then
       ! =======================================================================
       ! PROMPT AND TOTAL NU DATA -- read prompt data first
       KNU = JXS2 + 1
       LNU = XSS(KNU)
       if (LNU == 1) then
          ! Polynomial data
          table % nu_p_type = NU_POLYNOMIAL

          ! allocate determine how many coefficients for polynomial
          NC = XSS(KNU+1)
          length = NC + 1
       elseif (LNU == 2) then
          ! Tabular data
          table % nu_p_type = NU_TABULAR

          ! determine number of interpolation regions and number of energies
          NR = XSS(KNU+1)
          NE = XSS(KNU+2+2*NR)
          length = 2 + 2*NR + 2*NE
       end if

       ! allocate space for nu data storage
       allocate(table % nu_p_data(length))

       ! read data
       XSS_index = KNU + 1
       table % nu_p_data = get_real(length)

       ! Now read total nu data
       KNU = JXS2 + abs(XSS(JXS2)) + 1
       LNU = XSS(KNU)
       if (LNU == 1) then
          ! Polynomial data
          table % nu_t_type = NU_POLYNOMIAL

          ! allocate determine how many coefficients for polynomial
          NC = int(XSS(KNU+1))
          length = NC + 1
       elseif (LNU == 2) then
          ! Tabular data
          table % nu_t_type = NU_TABULAR

          ! determine number of interpolation regions and number of energies
          NR = int(XSS(KNU+1))
          NE = int(XSS(KNU+2+2*NR))
          length = 2 + 2*NR + 2*NE
       end if

       ! allocate space for nu data storage
       allocate(table % nu_t_data(length))

       ! read data
       XSS_index = KNU + 1
       table % nu_t_data = get_real(length)
    end if

    if (JXS24 > 0) then
       ! =======================================================================
       ! DELAYED NU DATA

       table % nu_d_type = NU_TABULAR
       KNU = JXS24

       ! determine size of tabular delayed nu data
       NR = int(XSS(KNU+1))
       NE = int(XSS(KNU+2+2*NR))
       length = 2 + 2*NR + 2*NE

       ! allocate space for delayed nu data
       allocate(table % nu_d_data(length))
       
       ! read delayed nu data
       XSS_index = KNU + 1
       table % nu_d_data = get_real(length)
 
       ! =======================================================================
       ! DELAYED NEUTRON ENERGY DISTRIBUTION

       ! Allocate space for secondary energy distribution
       NPCR = NXS(8)
       table % n_precursor = NPCR
       allocate(table % nu_d_edist(NPCR))

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
          table % nu_d_edist(i) % law = LAW

          ! allocate space for ENDF interpolation parameters
          if (NR > 0) then
             allocate(table % nu_d_edist(i) % nbt(NR))
             allocate(table % nu_d_edist(i) % int(NR))
          end if

          ! read ENDF interpolation parameters
          XSS_index = LDIS + LOCC + 3
          table % nu_d_edist(i) % nbt = get_real(NR)
          table % nu_d_edist(i) % int = get_real(NR)

          ! allocate space for law validity data
          NE = XSS(LDIS + LOCC + 3 + 2*NR)
          allocate(table % nu_d_edist(i) % energy(NE))
          allocate(table % nu_d_edist(i) % pvalid(NE))

          length_interp_data = 5 + 2*(NR + NE)

          ! read law validity data
          XSS_index = LDIS + LOCC + 4 + 2*NR
          table % nu_d_edist(i) % energy = get_real(NE)
          table % nu_d_edist(i) % pvalid = get_real(NE)

          ! Set index to beginning of IDAT array
          loc = LDIS + IDAT - 2

          ! determine length of energy distribution
          length = length_energy_dist(loc, LAW, LOCC, length_interp_data)

          ! allocate secondary energy distribution array
          allocate(table % nu_d_edist(i) % data(length))

          ! read secondary energy distribution
          XSS_index = loc + 1
          table % nu_d_edist(i) % data = get_real(length)
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
       allocate(table % nu_d_precursor_data(length))

       ! read delayed neutron precursor data
       XSS_index = loc
       table % nu_d_precursor_data = get_real(length)

    else
       table % nu_d_type = NU_NONE
    end if

  end subroutine read_nu_data

!===============================================================================
! READ_REACTIONS - Get the list of reaction MTs for this cross-section
! table. The MT values are somewhat arbitrary. Also read in Q-values, neutron
! multiplicities, and cross-sections.
!===============================================================================

  subroutine read_reactions(table)

    type(AceContinuous), pointer :: table

    integer :: LMT   ! index of MT list in XSS
    integer :: NMT   ! Number of reactions
    integer :: JXS4  ! index of Q values in XSS
    integer :: JXS5  ! index of neutron multiplicities in XSS
    integer :: JXS7  ! index of reactions cross-sections in XSS
    integer :: LXS   ! location of cross-section locators
    integer :: LOCA  ! location of cross-section for given MT
    integer :: NE    ! number of energies for reaction
    type(AceReaction), pointer :: rxn => null()
    
    LMT  = JXS(3)
    JXS4 = JXS(4)
    JXS5 = JXS(5)
    LXS  = JXS(6)
    JXS7 = JXS(7)
    NMT  = NXS(4)

    ! allocate array of reactions. Add one since we need to include an elastic
    ! scattering channel
    table%n_reaction = NMT + 1
    allocate(table%reactions(NMT+1))

    ! Store elastic scattering cross-section on reaction one
    rxn => table%reactions(1)
    rxn%MT      = 2
    rxn%Q_value = ZERO
    rxn%TY      = 1
    rxn%IE      = 1
    allocate(rxn%sigma(table%n_grid))
    rxn%sigma = table%sigma_el

    do i = 1, NMT
       rxn => table%reactions(i+1)

       ! read MT number, Q-value, and neutrons produced
       rxn % MT      = XSS(LMT+i-1)
       rxn % Q_value = XSS(JXS4+i-1)
       rxn % TY      = XSS(JXS5+i-1)

       ! read cross section values
       LOCA = XSS(LXS+i-1)
       rxn % IE = XSS(JXS7 + LOCA - 1)
       NE = XSS(JXS7 + LOCA)
       allocate(rxn%sigma(NE))
       XSS_index = JXS7 + LOCA + 1
       rxn % sigma = get_real(NE)

       ! set defaults
       rxn % has_angle_dist  = .false.
       rxn % has_energy_dist = .false.
    end do

  end subroutine read_reactions

!===============================================================================
! READ_ANGULAR_DIST parses the angular distribution for each reaction with
! secondary neutrons
!===============================================================================

  subroutine read_angular_dist(table)

    type(AceContinuous), pointer :: table

    integer :: JXS8   ! location of angular distribution locators
    integer :: JXS9   ! location of angular distributions
    integer :: LOCB   ! location of angular distribution for given MT
    integer :: NE     ! number of incoming energies
    integer :: NP     ! number of points for cosine distribution
    integer :: LC     ! locator
    integer :: i      ! index in reactions array
    integer :: j      ! index over incoming energies
    integer :: length ! length of data array to allocate
    type(AceReaction), pointer :: rxn => null()

    JXS8 = JXS(8)
    JXS9 = JXS(9)

    ! loop over all reactions with secondary neutrons -- NXS(5) does not include
    ! elastic scattering
    do i = 1, NXS(5) + 1
       rxn => table%reactions(i)

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

       ! read incoming energy grid and location of tables
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

  subroutine read_energy_dist(table)

    type(AceContinuous), pointer :: table

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
    integer :: i, j, k, l         ! indices
    type(AceReaction), pointer :: rxn => null()

    LED  = JXS(10)
    LDIS = JXS(11)

    ! Loop over all reactions 
    do i = 1, NXS(5)
       rxn => table % reactions(i+1) ! skip over elastic scattering
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
       rxn % edist % nbt = get_real(NR)
       rxn % edist % int = get_real(NR)

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

  subroutine read_unr_res(table)

    type(AceContinuous), pointer :: table

    integer :: JXS23 ! location of URR data
    integer :: loc   ! locator
    integer :: N     ! # of incident energies
    integer :: M     ! # of probabilities
    integer :: i     ! index over incoming energies
    integer :: j     ! index over values
    integer :: k     ! index over cumulative probabilities

    ! determine locator for URR data
    JXS23 = JXS(23)

    ! check if URR data is present
    if (JXS23 /= 0) then
       table % urr_present = .true.
       allocate(table % urr_params(6))
       loc = JXS23
    else
       table % urr_present = .false.
       return
    end if

    ! read parameters
    table % urr_params(1) = XSS(loc)     ! # of incident energies
    table % urr_params(2) = XSS(loc + 1) ! # of probabilities
    table % urr_params(3) = XSS(loc + 2) ! interpolation parameter
    table % urr_params(4) = XSS(loc + 3) ! inelastic competition flag
    table % urr_params(5) = XSS(loc + 4) ! other absorption flag
    table % urr_params(6) = XSS(loc + 5) ! factors flag

    ! allocate incident energies and probability tables
    N = table % urr_params(1)
    M = table % urr_params(2)
    allocate(table % urr_energy(N))
    allocate(table % urr_prob(N,6,M))

    ! read incident energies
    XSS_index = loc + 6
    table % urr_energy = get_real(N)

    ! read probability tables
    do i = 1, N
       do j = 1, 6
          table % urr_prob(i,j,1:M) = get_real(M)
       end do
    end do

  end subroutine read_unr_res

!===============================================================================
! READ_ACE_THERMAL reads in a single ACE thermal scattering S(a,b) table. This
! routine in particular reads the header data and then calls read_thermal_data
! to parse the actual data blocks.
!===============================================================================

  subroutine read_ACE_thermal(index_table, index)

    integer, intent(in) :: index_table ! index in xs_thermal array
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
    character(max_line_len) :: msg              ! output/error message
    character(max_line_len) :: line             ! single line to read
    character(max_word_len) :: words(max_words) ! words on a line
    character(max_word_len) :: filename         ! name of ACE library file
    character(10)           :: tablename        ! name of cross section table
    type(AceThermal), pointer :: table => null()

    filename = xsdatas(index)%path
    tablename = xsdatas(index)%id

    table => xs_thermal(index_table)

    ! Check if input file exists and is readable
    inquire(FILE=filename, EXIST=file_exists, READ=readable)
    if (.not. file_exists) then
       msg = "ACE library '" // trim(filename) // "' does not exist!"
       call fatal_error(msg)
    elseif (readable(1:3) == 'NO') then
       msg = "ACE library '" // trim(filename) // "' is not readable! &
            &Change file permissions with chmod command."
       call fatal_error(msg)
    end if

    ! display message
    msg = "Loading ACE cross section table: " // tablename
    call message(msg, 6)

    ! open file
    open(file=filename, unit=in, status='old', & 
         & action='read', iostat=ioError)
    if (ioError /= 0) then
       msg = "Error while opening file: " // filename
       call fatal_error(msg)
    end if

    found_xs = .false.
    do while (.not. found_xs)
       call read_line(in, line, ioError)
       if (ioError < 0) then
          msg = "Could not find ACE table " // tablename // "."
          call fatal_error(msg)
       end if
       call split_string(line, words, n)
       if (trim(words(1)) == trim(tablename)) then
          found_xs = .true.
          table%name = words(1)
          table%awr = str_to_real(words(2))
          kT = str_to_real(words(3))
          table%temp = kT / K_BOLTZMANN
       end if
       
       ! Skip 5 lines
       call skip_lines(in, 5, ioError)

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

    type(AceThermal), pointer :: table

    integer :: i      ! index for incoming energies
    integer :: j      ! index for outgoing energies
    integer :: k      ! index for outoging angles
    integer :: loc    ! location in XSS array
    integer :: NE_in  ! number of incoming energies
    integer :: NE_out ! number of outgoing energies
    integer :: NMU    ! number of outgoing angles
    integer :: JXS4   ! location of elastic energy table

    ! read number of inelastic energies and allocate arrays
    NE_in = XSS(JXS(1))
    table % n_inelastic_e_in = NE_in
    allocate(table % inelastic_e_in(NE_in))
    allocate(table % inelastic_sigma(NE_in))

    ! read inelastic energies and cross-sections
    XSS_index = JXS(1) + 1
    table % inelastic_e_in = get_real(NE_in)
    table % inelastic_sigma = get_real(NE_in)

    ! allocate space for outgoing energy/angle for inelastic
    ! scattering
    NE_out = NXS(4)
    NMU = NXS(3) + 1
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

       ! determine whether sigma=P or sigma = P/E
       table % n_elastic_type = NXS(5)
    else
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

    integer :: i
    integer :: n_isotopes
    integer :: IE
    real(8) :: density
    real(8) :: density_i
    real(8) :: sigma_i
    real(8) :: f
    type(AceContinuous), pointer :: table => null()
    type(AceReaction),   pointer :: rxn => null()

    ! initialize xs
    xs = ZERO

    ! find material atom density
    density = cMaterial%atom_density

    ! loop over all isotopes in material
    n_isotopes = cMaterial % n_isotopes
    do i = 1, n_isotopes
       table => xs_continuous(cMaterial % table(i))

       ! determine nuclide atom density
       density_i = cMaterial%atom_percent(i) * density

       ! search nuclide energy grid
       IE = table%grid_index(p % IE)
       f = (p%E - table%energy(IE))/(table%energy(IE+1) - table%energy(IE))
       
       ! handle special case of total cross section
       if (MT == 1) then
          xs = xs + density * (ONE-f) * table%sigma_t(IE) + & 
               & f * (table%sigma_t(IE+1))
          cycle
       end if

       ! loop over reactions in isotope
       do j = 1, table % n_reaction
          rxn => table % reactions(i)

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

end module ace
