module ace

  use global
  use output, only: error, message
  use string, only: lower_case
  use fileio, only: read_line, read_data, skip_lines
  use string, only: split_string, str_to_real
  use data_structures, only: dict_create, dict_add_key, dict_has_key, &
       &                     dict_get_key, dict_delete
  use endf,   only: reaction_name

  integer :: NXS(16)
  integer :: JXS(32)
  real(8), allocatable :: XSS(:)
  integer :: XSS_index

  private :: NXS
  private :: JXS
  private :: XSS

contains

!=====================================================================
! READ_XS reads all the cross sections for the problem and stores them
! in xs_continuous and xs_thermal arrays
!=====================================================================

  subroutine read_xs()

    type(Material),      pointer :: mat => null()
    type(xsData),        pointer :: iso => null()
    type(AceContinuous), pointer :: ace_cont => null()
    type(AceThermal),    pointer :: ace_thermal => null()
    integer :: i, j
    integer :: index
    character(10) :: key
    character(250) :: msg
    integer :: n
    integer :: index_continuous
    integer :: index_thermal
    type(DictionaryCI),   pointer :: temp_dict

    call dict_create(ace_dict)

    ! determine how many continuous-energy tables and how many S(a,b)
    ! thermal scattering tables there are
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
             n_thermal = n_thermal + 1
          case default
             msg = "Unknown cross section table type: " // key
             call error(msg)
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

!=====================================================================
! READ_ACE_CONTINUOUS reads in a single ACE continuous-energy cross
! section table
!=====================================================================

  subroutine read_ACE_continuous(index_table, index)

    integer, intent(in) :: index_table
    integer, intent(in) :: index

    type(AceContinuous), pointer :: table => null()
    integer :: i
    integer :: in = 7
    integer :: ioError
    integer :: words_per_line
    integer :: lines
    integer :: n
    logical :: file_exists
    logical :: found_xs
    character(7) :: readable
    character(250) :: msg, line
    character(32) :: words(max_words)
    character(100) :: filename
    character(10) :: tablename
    real(8) :: kT

    filename = xsdatas(index)%path
    tablename = xsdatas(index)%id

    table => xs_continuous(index_table)

    ! Check if input file exists and is readable
    inquire(FILE=filename, EXIST=file_exists, READ=readable)
    if (.not. file_exists) then
       msg = "ACE library '" // trim(filename) // "' does not exist!"
       call error(msg)
    elseif (readable(1:3) == 'NO') then
       msg = "ACE library '" // trim(filename) // "' is not readable! &
            &Change file permissions with chmod command."
       call error(msg)
    end if

    ! display message
    msg = "Loading ACE cross section table: " // tablename
    call message(msg, 6)

    ! open file
    open(file=filename, unit=in, status='old', & 
         & action='read', iostat=ioError)
    if (ioError /= 0) then
       msg = "Error while opening file: " // filename
       call error(msg)
    end if

    found_xs = .false.
    do while (.not. found_xs)
       call read_line(in, line, ioError)
       if (ioError < 0) then
          msg = "Could not find ACE table " // tablename // "."
          call error(msg)
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

    call parse_ESZ(table)
    call parse_MTR(table)
    call read_angular_dist(table)
    call read_energy_dist(table)

    ! Free memory from XSS array
    if(allocated(XSS)) deallocate(XSS)
    if(associated(table)) nullify(table)

    close(unit=in)

  end subroutine read_ACE_continuous

!=====================================================================
! PARSE_ESZ - reads through the ESZ block. This block contains the
! energy grid, total xs, absorption xs, elastic scattering xs, and
! heating numbers.
!=====================================================================

  subroutine parse_ESZ(table)

    type(AceContinuous), pointer, intent(inout) :: table

    integer :: NE

    ! determine number of energy points
    NE = NXS(3)
    table%n_grid = NE

    ! allocate storage for arrays
    allocate(table%energy(NE))
    allocate(table%sigma_t(NE))
    allocate(table%sigma_a(NE))
    allocate(table%sigma_el(NE))
    allocate(table%heating(NE))

    ! read data from XSS -- right now the total, absorption and
    ! elastic scattering are read in to these special arrays, but in
    ! reality, it should be necessary to only store elastic scattering
    ! and possible total cross-section for total material xs
    ! generation.
    XSS_index = 1
    table%energy = get_real(NE)
    table%sigma_t = get_real(NE)
    table%sigma_a = get_real(NE)
    table%sigma_el = get_real(NE)
    table%heating = get_real(NE)
    
  end subroutine parse_ESZ

!=====================================================================
! PARSE_MTR - Get the list of reaction MTs for this cross-section
! table. The MT values are somewhat arbitrary. Also read in Q-values,
! neutron multiplicities, and cross-sections.
! =====================================================================

  subroutine parse_MTR(table)

    type(AceContinuous), pointer :: table

    type(AceReaction), pointer :: rxn
    integer :: LMT   ! index of MT list in XSS
    integer :: NMT   ! Number of reactions
    integer :: JXS4  ! index of Q values in XSS
    integer :: JXS5  ! index of neutron multiplicities in XSS
    integer :: JXS7  ! index of reactions cross-sections in XSS
    integer :: LXS   ! 
    integer :: LOCA  ! 
    integer :: NE    ! number of energies for reaction
    
    LMT = JXS(3)
    JXS4 = JXS(4)
    JXS5 = JXS(5)
    LXS  = JXS(6)
    JXS7 = JXS(7)
    NMT = NXS(4)

    ! allocate array of reactions. Add one since we need to include an
    ! elastic scattering channel
    table%n_reaction = NMT + 1
    allocate(table%reactions(NMT+1))

    ! Store elastic scattering cross-section on reaction one
    rxn => table%reactions(1)
    rxn%MT      = 2
    rxn%Q_value = 0.0_8
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

  end subroutine parse_MTR

!=====================================================================
! READ_ANGULAR_DIST parses the angular distribution for each reaction
!=====================================================================

  subroutine read_angular_dist(table)

    type(AceContinuous), pointer :: table

    type(AceReaction), pointer :: rxn
    integer :: JXS8
    integer :: JXS9
    integer :: NMT
    integer :: LOCB
    integer :: NE
    integer :: NP
    integer :: LC
    integer :: i, j
    integer :: size

    JXS8 = JXS(8)
    JXS9 = JXS(9)

    ! loop over all reactions with secondary neutrons -- NXS(5) does
    ! not include elastic scattering
    do i = 1, NXS(5) + 1
       rxn => table%reactions(i)

       ! find location of angular distribution
       LOCB = XSS(JXS8 + i - 1)
       if (LOCB == -1) then
          ! Angular distribution data are specified through LAWi = 44
          ! in the DLW block
          cycle
       elseif (LOCB == 0) then
          ! No angular distribution data are given for this reaction,
          ! isotropic scattering is asssumed (in CM if TY < 0 and in
          ! LAB if TY > 0)
          cycle
       end if
       rxn % has_angle_dist = .true.

       ! allocate space for incoming energies and locations
       NE = XSS(JXS9 + LOCB - 1)
       rxn % adist_n_energy = NE
       allocate(rxn % adist_energy(NE))
       allocate(rxn % adist_location(NE))

       ! read incoming energy grid and location of tables
       XSS_index = JXS9 + LOCB
       rxn % adist_energy   = get_real(NE)
       rxn % adist_location = get_int(NE)

       ! determine dize of data block
       size = 0
       do j = 1, NE
          LC = rxn % adist_location(j)
          if (LC > 0) then
             ! 32 equiprobable bins
             size = size + 32
          elseif (LC < 0) then
             ! tabular distribution
             NP = XSS(JXS9 + abs(LC))
             size = size + 2 + 3*NP
          end if
       end do
       
       ! allocate angular distribution data and read
       allocate(rxn % adist_data(size))

       ! read angular distribution -- currently this does not actually
       ! parse the angular distribution tables for each incoming
       ! energy, that must be done on-the-fly
       LC = rxn % adist_location(1)
       XSS_index = JXS9 + abs(LC) - 1
       rxn % adist_data = get_real(size)

!       print *, trim(table%name), rxn%MT, size

    end do

  end subroutine read_angular_dist

!=====================================================================
! READ_ENERGY_DIST
!=====================================================================

  subroutine read_energy_dist(table)

    type(AceContinuous), pointer :: table

    type(AceReaction), pointer :: rxn
    integer :: LOCC
    integer :: LED   ! location of LDLW block
    integer :: LDIS  ! location of DLW block
    integer :: LNW   ! location of next law
    integer :: LAW
    integer :: NR, NE
    integer :: NMU
    integer :: NP
    integer :: NP2
    integer :: NRa, NEa, NRb, NEb
    integer :: IDAT
    integer :: start, size, size_interp_data
    integer :: i, j, k

    LED = JXS(10)
    LDIS = JXS(11)

    ! Loop over all reactions 
    do i = 1, NXS(5)
       rxn => table % reactions(i+1) ! skip over elastic scattering

       ! find location of energy distribution data
       LOCC = XSS(LED + i - 1)
       
       LNW  = XSS(LDIS + LOCC - 1)
       LAW  = XSS(LDIS + LOCC)
       IDAT = XSS(LDIS + LOCC + 1)
       NR   = XSS(LDIS + LOCC + 2)

       ! allocate space for ENDF interpolation parameters
       if (NR > 0) then
          allocate(rxn % edist_nbt(NR))
          allocate(rxn % edist_int(NR))
       end if
       
       ! read ENDF interpolation parameters
       XSS_index = LDIS + LOCC + 3
       rxn % edist_nbt = get_real(NR)
       rxn % edist_int = get_real(NR)

       ! allocate space for law validity data
       NE = XSS(LDIS + LOCC + 3 + 2*NR)
       allocate(rxn % edist_energy(NE))
       allocate(rxn % edist_pvalid(NE))
       
       size_interp_data = 5 + 2*(NR + NE)
       
       ! read law validity data
       XSS_index = LDIS + LOCC + 4 + 2*NR
       rxn % edist_energy = get_real(NE)
       rxn % edist_pvalid = get_real(NE)

       ! Set index to beginning of IDAT array
       start = LDIS + IDAT - 2

       ! Determine size of LDAT array based on which secondary energy
       ! law it is
       size = 0
       select case (LAW)
       case (1)
          ! Tabular equiprobable energy bins
          NR = XSS(start + 1)
          NE = XSS(start + 2 + 2*NR)
          NET = XSS(start + 3 + 2*NR + NE)
          size = 3 + 2*NR + NE + 3*NET*NE

       case (2)
          ! Discrete photon energy
          size = 2
          
       case (3)
          ! Level scattering
          size = 2
         
       case (4)
          ! Continuous tabular distribution
          NR = XSS(start + 1)
          NE = XSS(start + 2 + 2*NR)
          size = size + 2 + 2*NR + 2*NE
          do j = 1,NE
             NP = XSS(start + size + 2)
             size = size + 2 + 3*NP
          end do

       case (5)
          ! General evaporation spectrum
          NR = XSS(start + 1)
          NE = XSS(start + 2 + 2*NR)
          NET = XSS(start + 3 + 2*NR + 2*NE)
          size = 3 + 2*NR + 2*NE + NET
          
       case (7)
          ! Maxwell fission spectrum
          NR = XSS(start + 1)
          NE = XSS(start + 2 + 2*NR)
          size = 3 + 2*NR + 2*NE

       case (9)
          ! Evaporation spectrum
          NR = XSS(start + 1)
          NE = XSS(start + 2 + 2*NR)
          size = 3 + 2*NR + 2*NE

       case (11)
          ! Watt spectrum
          NRa = XSS(start + 1)
          NEa = XSS(start + 2 + 2*NRa)
          NRb = XSS(start + 3 + 2*(NRa+NEa))
          NEb = XSS(start + 4 + 2*(NRa+NEa+NRb))
          size = 5 + 2*(NRa + NEa + NRb + NEb)

       case (44)
          ! Kalbach-Mann correlated scattering
          NR = XSS(start + 1)
          NE = XSS(start + 2 + 2*NR)
          size = size + 2 + 2*NR + 2*NE
          do j = 1,NE
             NP = XSS(start + size + 2)
             size = size + 2 + 5*NP
          end do

       case (61)
          ! Correlated energy and angle distribution
          NR = XSS(start + 1)
          NE = XSS(start + 2 + 2*NR)
          size = size + 2 + 2*NR + 2*NE
          do j = 1,NE
             ! outgoing energy distribution
             NP = XSS(start + size + 2)
             size = size + 2 + 4*NP
             do k = 1, NP
                ! outgoing angle distribution
                NP2 = XSS(start + size + 2)
                size = size + 2 + 3*NP2
             end do
          end do

       case (66)
          ! N-body phase space distribution
          size = 2

       case (67)
          ! Laboratory energy-angle law
          NR = XSS(start + 1)
          NE = XSS(start + 2 + 2*NR)
          NMU = XSS(start + 4 + 2*NR + 2*NE)
          size = 4 + 2*(NR + NE + NMU)
          
       end select

       ! allocate secondary energy distribution array
       allocate(rxn % edist_data(size))

       ! read secondary energy distribution
       XSS_index = start + 1
       rxn % edist_data = get_real(size)

    end do

  end subroutine read_energy_dist

!=====================================================================
! GET_INT returns an array of integers read from the current position
! in the XSS array
!=====================================================================

    function get_int(n_values) result(array)

      integer, intent(in) :: n_values
      integer :: array(n_values)

      array = int(XSS(XSS_index:XSS_index + n_values - 1))
      XSS_index = XSS_index + n_values
      
    end function get_int

!=====================================================================
! GET_REAL returns an array of real(8)s read from the current position
! in the XSS array
!=====================================================================

    function get_real(n_values) result(array)

      integer, intent(in) :: n_values
      real(8) :: array(n_values)

      array = XSS(XSS_index:XSS_index + n_values - 1)
      XSS_index = XSS_index + n_values

    end function get_real

!=====================================================================
! READ_ACE_THERMAL reads in a single ACE S(a,b) thermal scattering
! cross section table
!=====================================================================

  subroutine read_ACE_thermal(filename)

    character(*), intent(in) :: filename

  end subroutine read_ACE_thermal

end module ace
