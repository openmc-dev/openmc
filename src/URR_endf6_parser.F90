module URR_endf6_parser

  use constants
  use error,             only: fatal_error, warning
  use global
  use output,            only: write_message

  implicit none
  private
  public :: read_endf6

  integer :: in = 11 ! input unit
  character(80) :: filename ! ENDF-6 filename

contains

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_ENDF6 reads in an ENDF-6 format file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_endf6(filename_tmp, i)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    character(80) :: filename_tmp ! temporary filename
    character(80) :: rec          ! ENDF file record
    character(7)  :: readable     ! is ENDF-6 file readable?
    logical :: file_exists ! does ENDF-6 file exist?
    logical :: MF1_read    ! has MF=1 been read?
    logical :: MF2_read    ! has MF=2 been read?
    logical :: MF3_read    ! has MF=3 been read?
    integer :: i  ! isotope index
    integer :: MF ! MF file number
    integer :: MT ! MT type number
    integer :: NS ! record number

    filename = filename_tmp

    inquire(file = trim(path_endf)//trim(filename), &
      & exist = file_exists, read = readable)
    if (.not. file_exists) then
      call fatal_error('ENDF-6 file '//trim(filename)//' does not exist.')
    else if (readable(1:3) == 'NO') then
      call fatal_error('ENDF-6 file '//trim(filename)// &
        & ' is not readable.  Change file permissions with chmod command.')
    end if

    ! display message
    call write_message("Loading ENDF-6 file: "//trim(filename), 6)

    open(unit = in, &
      file = trim(path_endf)//trim(filename))

    tope => isotopes(i)

    MF1_read = .false.
    MF2_read = .false.
    MF3_read = .false.

    do
      read(in, 10) rec
10    format(A80)
      read(rec(67:70), '(I4)') tope % MAT
      read(rec(71:72), '(I2)') MF
      read(rec(73:75), '(I3)') MT
      read(rec(76:80), '(I5)') NS
      if (MF == 1 .and. MT == 451 .and. (.not. MF1_read)) then
        call read_MF1(i, rec)
        MF1_read = .true.
      else if (MF == 2 .and. MT == 151 .and. (.not. MF2_read)) then
        call read_MF2(i, rec)
        MF2_read = .true.
      else if (MF == 3 .and. MT == 1 .and. (.not. MF3_read)) then
        call read_MF3(i, rec)
        MF3_read = .true.
      else if (tope % MAT == -1) then
        exit
      end if
    end do

    close(in)

  end subroutine read_endf6

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_MF1 reads in an ENDF-6 format MF 1 file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_MF1(i, rec)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    character(80) :: rec ! ENDF-6 file record
    integer :: i ! isotope index
    real(8) :: real_ZAI

    tope => isotopes(i)

    read(rec(1:11),  '(E11.0)') real_ZAI
    tope % ZAI = int(real_ZAI)
    read(rec(12:22), '(E11.0)') tope % AWR
    read(rec(23:33), '(I11)')   tope % LRP

! TODO:    ! check ZAID agreement
    call check_zaid(tope % ZAI, tope % ZAI)

! TODO:    ! check mass ratios
    call check_mass(tope % AWR, tope % AWR)

    ! check that resonance parameters are given
    call check_parameters(tope % LRP)

  end subroutine read_MF1

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_MF2 reads in an ENDF-6 format MF 2 file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_MF2(i, rec)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    character(80) :: rec ! ENDF-6 file record
    integer :: i    ! index in global isotopes array
    integer :: i_ER ! resonance energy range index
    integer :: NIS  ! number of isotopes in material
    real(8) :: ZA
    real(8) :: A
    real(8) :: ABN

    tope => isotopes(i)

    read(rec(1:11),  '(E11.0)') ZA
    read(rec(12:22), '(E11.0)') A
    read(rec(45:55), '(I11)')   NIS

    ! check ZAID agreement
    call check_zaid(int(ZA), tope % ZAI)

    ! check that mass is consistent
    call check_mass(A, tope % AWR)

    ! check that this is a single-isotope ENDF-6 file
    call check_single_isotope(NIS)

    ! read MF=2, record=2
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') ZA
    read(rec(12:22), '(E11.0)') ABN
    read(rec(34:44), '(I11)')   tope % LFW
    read(rec(45:55), '(I11)')   tope % NER

    ! allocate energy range variables
    call tope % alloc_energy_range()

    ! check ZAID agreement
    call check_zaid(int(ZA), tope % ZAI)

    ! check abundance is unity
    call check_abundance(ABN)

    ! check URR average fission widths treatment
    call check_fission_widths(tope % LFW)

    ! loop over energy ranges
    do i_ER = 1, tope % NER

      ! read first record for this energy range
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') tope % EL(i_ER)
      read(rec(12:22), '(E11.0)') tope % EH(i_ER)
      read(rec(24:33),   '(I10)') tope % LRU(i_ER)
      read(rec(35:44),   '(I10)') tope % LRF(i_ER)
      read(rec(45:55),   '(I11)') tope % NRO(i_ER)
      read(rec(56:66),   '(I11)') tope % NAPS(i_ER)

      ! check number of resonance energy ranges and their bounding energies
      call check_energy_ranges(tope % NER, tope % EL(i_ER), tope % EH(i_ER))

      ! check channel, scattering radius energy dependence flags
      call check_radius_flags(tope % NRO(i_ER), tope % NAPS(i_ER))

      ! read energy range and formalism-dependent resonance subsection data
      call read_resonance_subsection(i, i_ER, tope % LRU(i_ER), tope %LRF(i_ER))

    end do

  end subroutine read_MF2

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_MF3 reads in an ENDF-6 format MF 3 file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_MF3(i, rec)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    integer :: i      ! index in global isotopes array
    integer :: i_rec  ! record index
    integer :: i_e    ! index in energy grid
    integer :: MF     ! ENDF-6 MF flag
    integer :: MT     ! ENDF-6 MT flag
    integer :: NP     ! number of energy-xs pairs
    integer :: NR     ! number of interpolation regions
    integer :: NBT    ! number of entries separating interp. ranges N and N+1
    integer :: INTERP ! File 3 reaction xs interpolation flag
    character(80) :: rec ! ENDF-6 file record
    real(8) :: ZA ! ZAID
    real(8) :: A  ! mass in neutron masses
    real(8) :: QI ! reaction Q-value
    logical :: read_MT ! read this MT reaction?

    tope => isotopes(i)

10  format(A80)

    do
      read(in, 10) rec
      read(rec(71:72), '(I2)') MF
      read(rec(73:75), '(I3)') MT
      if (MF == 3 .and. MT == 2) exit
      if (MT > 2) call fatal_error('Reached end of MF3 w/o reading an&
        & elastic cross section in '//trim(filename))
    end do

    read(rec(1:11),  '(E11.0)') ZA
    read(rec(12:22), '(E11.0)') A

    ! check ZAID agreement
    call check_zaid(int(ZA), tope % ZAI)

    ! check that mass is consistent
    call check_mass(A, tope % AWR)

    ! read MF=3, record=2
    read(in, 10) rec
    read(rec(45:55), '(I11)') NR
    read(rec(56:66), '(I11)') NP

    ! check number of interpolation regions
    call check_interp_regions(NR)

    ! read MF=3, record=3
    read(in, 10) rec
    read(rec( 1:11), '(I11)') NBT
    read(rec(12:22), '(I11)') tope % MF3_INT

    ! check number of energy-xs pairs in this interpolation region
    call check_n_pairs(NBT, NP)

    allocate(tope % MF3_n_e(NP))
    allocate(tope % MF3_n(NP))

    i_e = 1
    do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
      read(in, 10) rec
      read(rec( 1:11), '(E11.0)') tope % MF3_n_e(i_e)
      read(rec(12:22), '(E11.0)') tope % MF3_n(i_e)
      if (i_e == NP) exit
      i_e = i_e + 1
      read(rec(23:33), '(E11.0)') tope % MF3_n_e(i_e)
      read(rec(34:44), '(E11.0)') tope % MF3_n(i_e)
      if (i_e == NP) exit
      i_e = i_e + 1
      read(rec(45:55), '(E11.0)') tope % MF3_n_e(i_e)
      read(rec(56:66), '(E11.0)') tope % MF3_n(i_e)
      if (i_e == NP) exit
      i_e = i_e + 1
    end do

    read_MT = .false.
    do
      read(in, 10) rec
      read(rec(71:72), '(I2)') MF
      read(rec(73:75), '(I3)') MT
      if (MF == 3 .and. MT == 18) then
        read_MT = .true.
        exit
      end if
      if (MT > 18) then
        backspace(in)
        call warning('Reached end of MF3 w/o reading a&
          & fission cross section in '//trim(filename))
        exit
      end if
    end do

    if (read_MT) then
      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') A

      ! check ZAID agreement
      call check_zaid(int(ZA), tope % ZAI)

      ! check that mass is consistent
      call check_mass(A, tope % AWR)

      ! read MF=3, record=2
      read(in, 10) rec
      read(rec(45:55), '(I11)') NR
      read(rec(56:66), '(I11)') NP

      ! check number of interpolation regions
      call check_interp_regions(NR)

      ! read MF=3, record=3
      read(in, 10) rec
      read(rec( 1:11), '(I11)') NBT
      read(rec(12:22), '(I11)') INTERP

      ! check number of energy-xs pairs in this interpolation region
      call check_n_pairs(NBT, NP)

      ! check interpolation scheme is same as isotope interpolation scheme
      call check_interp_scheme(INTERP, tope % MF3_INT)

      allocate(tope % MF3_f_e(NP))
      allocate(tope % MF3_f(NP))

      i_e = 1
      do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') tope % MF3_f_e(i_e)
        read(rec(12:22), '(E11.0)') tope % MF3_f(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(23:33), '(E11.0)') tope % MF3_f_e(i_e)
        read(rec(34:44), '(E11.0)') tope % MF3_f(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(45:55), '(E11.0)') tope % MF3_f_e(i_e)
        read(rec(56:66), '(E11.0)') tope % MF3_f(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
      end do
    end if

    read_MT = .false.
    do
      read(in, 10) rec
      read(rec(71:72), '(I2)') MF
      read(rec(73:75), '(I3)') MT
      if (MF == 3 .and. MT == 51) then
        read_MT = .true.
        exit
      end if
      if (MT > 51) then
        backspace(in)
        call warning('Reached end of MF3 w/o reading an&
          & (n,n1) cross section in '//trim(filename))
        exit
      end if
    end do

    if (read_MT) then
      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') A

      ! check ZAID agreement
      call check_zaid(int(ZA), tope % ZAI)

      ! check that mass is consistent
      call check_mass(A, tope % AWR)

      ! read MF=3, record=2
      read(in, 10) rec
      read(rec(12:22), '(E11.0)') QI
      tope % E_ex1 = (A + ONE) / A * (-QI)
      read(rec(45:55), '(I11)') NR
      read(rec(56:66), '(I11)') NP

      ! check number of interpolation regions
      call check_interp_regions(NR)

      ! read MF=3, record=3
      read(in, 10) rec
      read(rec( 1:11), '(I11)') NBT
      read(rec(12:22), '(I11)') INTERP

      ! check number of energy-xs pairs in this interpolation region
      call check_n_pairs(NBT, NP)

      ! check interpolation scheme is same as isotope interpolation scheme
      call check_interp_scheme(INTERP, tope % MF3_INT)

      allocate(tope % MF3_x_e(NP))
      allocate(tope % MF3_x(NP))

      i_e = 1
      do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') tope % MF3_x_e(i_e)
        read(rec(12:22), '(E11.0)') tope % MF3_x(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(23:33), '(E11.0)') tope % MF3_x_e(i_e)
        read(rec(34:44), '(E11.0)') tope % MF3_x(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(45:55), '(E11.0)') tope % MF3_x_e(i_e)
        read(rec(56:66), '(E11.0)') tope % MF3_x(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
      end do
    end if

    read_MT = .false.
    do
      read(in, 10) rec
      read(rec(71:72), '(I2)') MF
      read(rec(73:75), '(I3)') MT
      if (MF == 3 .and. MT == 52) then
        read_MT = .true.
        exit
      end if
      if (MT > 52) then
        backspace(in)
        call warning('Reached end of MF3 w/o reading an&
          & (n,n2) Q-value in '//trim(filename))
        exit
      end if
    end do

    if (read_MT) then
      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') A

      ! check ZAID agreement
      call check_zaid(int(ZA), tope % ZAI)

      ! check that mass is consistent
      call check_mass(A, tope % AWR)

      ! read MF=3, record=2
      read(in, 10) rec
      read(rec(12:22), '(E11.0)') QI
      tope % E_ex2 = (A + ONE) / A * (-QI)
    end if

    read_MT = .false.
    do
      read(in, 10) rec
      read(rec(71:72), '(I2)') MF
      read(rec(73:75), '(I3)') MT
      if (MF == 3 .and. MT == 102) then
        read_MT = .true.
        exit
      end if
      if (MT > 102 .or. MF == 4) then
        backspace(in)
        call warning('Reached end of MF3 w/o reading a&
          & capture cross section in '//trim(filename))
        exit
      end if
    end do

    if (read_MT) then
      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') A

      ! check ZAID agreement
      call check_zaid(int(ZA), tope % ZAI)

      ! check that mass is consistent
      call check_mass(A, tope % AWR)

      ! read MF=3, record=2
      read(in, 10) rec
      read(rec(45:55), '(I11)') NR
      read(rec(56:66), '(I11)') NP

      ! check number of interpolation regions
      call check_interp_regions(NR)

      ! read MF=3, record=3
      read(in, 10) rec
      read(rec( 1:11), '(I11)') NBT
      read(rec(12:22), '(I11)') INTERP

      ! check number of energy-xs pairs in this interpolation region
      call check_n_pairs(NBT, NP)

      ! check interpolation scheme is same as isotope interpolation scheme
      call check_interp_scheme(INTERP, tope % MF3_INT)

      allocate(tope % MF3_g_e(NP))
      allocate(tope % MF3_g(NP))

      i_e = 1
      do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') tope % MF3_g_e(i_e)
        read(rec(12:22), '(E11.0)') tope % MF3_g(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(23:33), '(E11.0)') tope % MF3_g_e(i_e)
        read(rec(34:44), '(E11.0)') tope % MF3_g(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(45:55), '(E11.0)') tope % MF3_g_e(i_e)
        read(rec(56:66), '(E11.0)') tope % MF3_g(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
      end do
    end if

  end subroutine read_MF3

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_INTERP_REGIONS checks that there is an allowable number of interpolation
! regions
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_interp_regions(NR)

    integer :: NR ! number of interpolation regions

    if (NR /= 1) then
      if (master) call warning('More than 1 File 3 interpolation region&
        & in '//trim(filename))
    end if

  end subroutine check_interp_regions

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_N_PAIRS checks that there is an allowable number of energy-xs pairs
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_n_pairs(NBT, NP)

    integer :: NBT ! number of pairs between this and the next interp. region
    integer :: NP  ! number of energy-xs pairs

    if (NBT /= NP) then
      if (master) call warning('Different NBT and NP values in File 3 for '&
        &//trim(filename))
    end if

  end subroutine check_n_pairs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_INTERP_SCHEME checks that there is an allowable interpolation scheme
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_interp_scheme(INTERP, MF3_INT)

    integer :: INTERP     ! interpolation scheme for current region
    integer :: MF3_INT ! overall isotope interpolation scheme

    if (INTERP /= MF3_INT) call warning('Different interpolation schemes for &
      &different File 3 regions in '//trim(filename)//'.  Using region &
      &one interpolation scheme for all regions.')

  end subroutine check_interp_scheme

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_RESONANCE_SUBSECTION reads in the energy range and formalism-dependent
! resonance subsection data
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_resonance_subsection(i, i_ER, LRU, LRF)

    integer :: i    ! index in global isotopes array
    integer :: i_ER ! energy range index
    integer :: LRU  ! resolved (1) or unresolved (2) parameters
    integer :: LRF  ! ENDF-6 resonance formalism flag

    ! select energy range
    select case(LRU)

    ! only the scattering radius is given (LRF=0; NLS=0; LFW=0)
    case(0)

      call fatal_error('Scattering radius only (LRU = 0) not supported in &
        &'//trim(filename))

    ! resolved parameters
    case(1)

      ! select formalism
      select case(LRF)

      case(SLBW)
        call read_slbw_parameters(i, i_ER)
      case(MLBW)
        call read_mlbw_parameters(i, i_ER)
      case(REICH_MOORE)
        call read_rm_parameters(i, i_ER)
      case(ADLER_ADLER)
        call fatal_error('Adler-Adler (LRF=4) formalism not supported in '&
          & //trim(filename))
      case(R_MATRIX)
        call fatal_error('General R-Matrix (LRF=5) formalism not allowed in '&
          & //trim(filename))
      case(R_FUNCTION)
        call fatal_error('Hybrid R-Function (LRF=6) formalism not allowed in '&
          & //trim(filename))
      case(R_MATRIX_LIM)
        call fatal_error('R-Matrix Limited (LRF=7) formalism not supported in '&
          & //trim(filename))

      ! default case
      case default
        call fatal_error('LRF must be an integer between 1 and 7, inclusive in '&
          & //trim(filename))
      end select

    ! unresolved parameters
    case(2)

      isotopes(i) % i_urr = i_ER

      ! select energy-dependence of parameters
      select case(LRF)

      ! only fission widths are energy-dependent
      case(1)
        call read_urr_slbw_parameters_lrf1(i, i_ER)

      ! all parameters are energy-dependent
      case(2)
        call read_urr_slbw_parameters_lrf2(i, i_ER)

      ! default case
      case default
        call fatal_error('LRF values other than 2 are not supported in '&
          & //trim(filename))
      end select

    ! default case
    case default
      call fatal_error('LRU values other than 0, 1 or 2 are not allowed in '&
        & //trim(filename))
    end select

  end subroutine read_resonance_subsection

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_SLBW_PARAMETERS reads in Single-level Breit-Wigner resolved resonance
! region parameters
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_slbw_parameters(i, i_ER)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    character(80) :: rec ! ENDF-6 file record
    integer :: i    ! index in global isotopes array
    integer :: i_ER ! energy range index
    integer :: i_l  ! orbital quantum number index
    integer :: L    ! orbital quantum number
    integer :: NRS  ! number of resonances for this orbital quantum number
    integer :: i_R  ! resonance index
    integer :: LRX  ! competitive width flag
    real(8) :: A    ! isotope/neutron mass ratio
    real(8) :: QX   ! Q-value to be added to COM energy

    tope => isotopes(i)

    ! read first line of energy range subsection
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') tope % SPI(i_ER)
    read(rec(12:22), '(E11.0)') tope % AP(i_ER)
    call tope % channel_radius(i_ER)
    read(rec(45:55), '(I11)')   tope % NLS(i_ER)
    if (tope % NLS(i_ER) > 3) &
         call fatal_error('SLBW parameters given for a resonance higher than&
         & d-wave')

    ! allocate SLBW resonance vectors for each l
    allocate(tope % bw_resonances(tope % NLS(i_ER)))

    ! loop over orbital quantum numbers
    do i_l = 0, tope % NLS(i_ER) - 1
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') A
      read(rec(12:22), '(E11.0)') QX
      read(rec(23:33),   '(I11)') L
      read(rec(34:44),   '(I11)') LRX
      read(rec(56:66),   '(I11)') NRS

      ! allocate SLBW resonances
      call tope % bw_resonances(i_l + 1) % alloc(NRS)

      ! check mass ratios
      call check_mass(A, tope % AWR)

      ! check that Q-value is 0.0
      call check_q_value(QX)

      ! check orbital quantum number
      call check_l_number(L, i_L)

      ! loop over resonances
      do i_R = 1, NRS

        read(in, 10) rec
        read(rec(1:11),  '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % E_lam
        read(rec(12:22), '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % AJ
        read(rec(23:33), '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % GT
        read(rec(34:44), '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % GN
        read(rec(45:55), '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % GG
        read(rec(56:66), '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % GF

        if (LRX == 0) then
          tope % bw_resonances(i_l + 1) % res(i_R) % GX = ZERO
        else if (LRX == 1) then
          tope % bw_resonances(i_l + 1) % res(i_R) % GX&
               = tope % bw_resonances(i_l + 1) % res(i_R) % GT&
               - tope % bw_resonances(i_l + 1) % res(i_R) % GN&
               - tope % bw_resonances(i_l + 1) % res(i_R) % GG&
               - tope % bw_resonances(i_l + 1) % res(i_R) % GF
        else
          call fatal_error('LRX must be 0 or 1 in '//trim(filename))
        end if

        ! check sign of total angular momentum, J
        call check_j_sign(tope % bw_resonances(i_l + 1) % res(i_R) % AJ)
      end do

      if (i_ER < tope % NER - 1) then
        call tope % bw_resonances(i_l + 1) % clear()
      end if

    end do

    if (i_ER < tope % NER - 1) deallocate(tope % bw_resonances)

  end subroutine read_slbw_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_MLBW_PARAMETERS reads in Multi-level Breit-Wigner resolved resonance
! region parameters
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_mlbw_parameters(i, i_ER)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    character(80) :: rec ! ENDF-6 file record
    integer :: i    ! index in global isotopes array
    integer :: i_ER ! energy range index
    integer :: i_l  ! orbital quantum number index
    integer :: L    ! orbital quantum number
    integer :: NRS  ! number of resonances for this orbital quantum number
    integer :: i_R  ! resonance index
    integer :: LRX  ! competitive width flag
    real(8) :: A    ! isotope/neutron mass ratio
    real(8) :: QX   ! Q-value to be added to COM energy

    tope => isotopes(i)

    ! read first line of energy range subsection
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') tope % SPI(i_ER)
    read(rec(12:22), '(E11.0)') tope % AP(i_ER)
    call tope % channel_radius(i_ER)
    read(rec(45:55), '(I11)')   tope % NLS(i_ER)
    if (tope % NLS(i_ER) > 3) &
         call fatal_error('MLBW parameters given for a resonance higher than&
         & d-wave')

    ! allocate MLBW resonance vectors for each l
    allocate(tope % bw_resonances(tope % NLS(i_ER)))

    ! loop over orbital quantum numbers
    do i_l = 0, tope % NLS(i_ER) - 1
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') A
      read(rec(12:22), '(E11.0)') QX
      read(rec(23:33),   '(I11)') L
      read(rec(34:44),   '(I11)') LRX
      read(rec(56:66),   '(I11)') NRS

      ! allocate MLBW resonances
      call tope % bw_resonances(i_l + 1) % alloc(NRS)

      ! check mass ratios
      call check_mass(A, tope % AWR)

      ! check that Q-value is 0.0
      call check_q_value(QX)

      ! check orbital quantum number
      call check_l_number(L, i_L)

      ! loop over resonances
      do i_R = 1, NRS

        read(in, 10) rec
        read(rec(1:11),  '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % E_lam
        read(rec(12:22), '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % AJ
        read(rec(23:33), '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % GT
        read(rec(34:44), '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % GN
        read(rec(45:55), '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % GG
        read(rec(56:66), '(E11.0)')&
             tope % bw_resonances(i_l + 1) % res(i_R) % GF

        if (LRX == 0) then
          tope % bw_resonances(i_l + 1) % res(i_R) % GX = ZERO
        else if (LRX == 1) then
          tope % bw_resonances(i_l + 1) % res(i_R) % GX&
               = tope % bw_resonances(i_l + 1) % res(i_R) % GT&
               - tope % bw_resonances(i_l + 1) % res(i_R) % GN&
               - tope % bw_resonances(i_l + 1) % res(i_R) % GG&
               - tope % bw_resonances(i_l + 1) % res(i_R) % GF
        else
          call fatal_error('LRX must be 0 or 1 in '//trim(filename))
        end if

        ! check sign of total angular momentum, J
        call check_j_sign(tope % bw_resonances(i_l + 1) % res(i_R) % AJ)
      end do

      if (i_ER < tope % NER - 1) then
        call tope % bw_resonances(i_l + 1) % clear()
      end if
    end do

    if (i_ER < tope % NER - 1) deallocate(tope % bw_resonances)

  end subroutine read_mlbw_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_RM_PARAMETERS reads in Reich-Moore resolved resonance region
! parameters
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_rm_parameters(i, i_ER)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    character(80) :: rec ! ENDF-6 file record
    integer :: i    ! index in global isotopes array
    integer :: i_ER ! energy range index
    integer :: i_l  ! orbital quantum number index
    integer :: L    ! orbital quantum number
    integer :: NRS  ! number of resonances for this orbital quantum number
    integer :: i_R  ! resonance index
    real(8) :: A    ! isotope/neutron mass ratio
    real(8) :: APL  ! l-dependent AP value

    tope => isotopes(i)

    ! read first line of energy range subsection
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') tope % SPI(i_ER)
    read(rec(12:22), '(E11.0)') tope % AP(i_ER)
    call tope % channel_radius(i_ER)
    read(rec(45:55), '(I11)')   tope % NLS(i_ER)
    if (tope % NLS(i_ER) > 3) &
         call warning('R-M parameters given for a resonance higher than&
         & d-wave')

    ! allocate Reich-Moore resonance vectors for each l
    allocate(tope % rm_resonances(tope % NLS(i_ER)))

    ! loop over orbital quantum numbers
    do i_l = 0, tope % NLS(i_ER) - 1
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') A
      read(rec(12:22), '(E11.0)') APL
      read(rec(23:33),   '(I11)') L
      read(rec(56:66),   '(I11)') NRS

      ! allocate Reich-Moore resonances
      call tope % rm_resonances(i_l + 1) % alloc(NRS)

      ! check mass ratios
      call check_mass(A, tope % AWR)

      ! check scattering radii
      call check_scattering_radius(APL, tope % AP(i_ER))

      ! check orbital quantum number
      call check_l_number(L, i_L)

      ! loop over resonances
      do i_R = 1, NRS

        read(in, 10) rec
        read(rec(1:11),  '(E11.0)')&
             tope % rm_resonances(i_l + 1) % res(i_R) % E_lam
        read(rec(12:22), '(E11.0)')&
             tope % rm_resonances(i_l + 1) % res(i_R) % AJ
        read(rec(23:33), '(E11.0)')&
             tope % rm_resonances(i_l + 1) % res(i_R) % GN
        read(rec(34:44), '(E11.0)')&
             tope % rm_resonances(i_l + 1) % res(i_R) % GG
        read(rec(45:55), '(E11.0)')&
             tope % rm_resonances(i_l + 1) % res(i_R) % GFA
        read(rec(56:66), '(E11.0)')&
             tope % rm_resonances(i_l + 1) % res(i_R) % GFB

        ! check sign of total angular momentum, J
        call check_j_sign(tope % rm_resonances(i_l + 1) % res(i_R) % AJ)
      end do

      if (i_ER < tope % NER - 1) then
        call tope % rm_resonances(i_l + 1) % clear()
      end if

    end do

    if (i_ER < tope % NER - 1) deallocate(tope % rm_resonances)

  end subroutine read_rm_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_URR_SLBW_PARAMETERS_LRF1 reads in unresolved resonance region parameters
! for the LRF = 1 option (only fission widths are energy-dependent)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_urr_slbw_parameters_lrf1(i, i_ER)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    character(80) :: rec ! ENDF-6 file record
    integer :: i    ! index in global isotopes array
    integer :: i_ER ! energy range index
    integer :: i_l  ! orbital quantum number index
    integer :: L    ! orbital quantum number
    integer :: i_J  ! total angular momentum quantum number index
    integer :: i_ES ! tabulated fission width energy grid index
    real(8) :: A    ! isotope/neutron mass ratio

    tope => isotopes(i)

    ! this is forced by the ENDF-6 format when LRF = 1 and LFW = 1
    tope % INT = LINEAR_LINEAR

10  format(A80)
    select case(tope % LFW)
    case(0)
      ! read first line of energy range subsection
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') tope % SPI(i_ER)
      read(rec(12:22), '(E11.0)') tope % AP(i_ER)
! TODO: don't overwrite the resolved value
      call tope % channel_radius(i_ER)
      read(rec(23:33), '(I11)')   tope % LSSF
      read(rec(45:55), '(I11)')   tope % NLS(i_ER)
      if (tope % NLS(i_ER) > 3) &
           call fatal_error('URR parameters given for a spin sequence higher than&
           & d-wave')
      
      ! allocate number of total angular momenta values for each l
      allocate(tope % NJS(tope % NLS(i_ER)))

      ! allocate total angular momenta
      allocate(tope % AJ(tope % NLS(i_ER)))

      ! allocate degrees of freedom for partial widths
      allocate(tope % DOFX(tope % NLS(i_ER)))
      allocate(tope % DOFN(tope % NLS(i_ER)))
      allocate(tope % DOFG(tope % NLS(i_ER)))
      allocate(tope % DOFF(tope % NLS(i_ER)))

      ! allocate partial widths
      tope % NE = 2
      allocate(tope % ES(tope % NE))
      tope % ES(1) = tope % EL(i_ER)
      tope % ES(2) = tope % EH(i_ER)
      allocate(tope % D_mean  (tope % NLS(i_ER)))
      allocate(tope % GN0_mean(tope % NLS(i_ER)))
      allocate(tope % GG_mean (tope % NLS(i_ER)))
      allocate(tope % GF_mean (tope % NLS(i_ER)))
      allocate(tope % GX_mean (tope % NLS(i_ER)))

      if (tope%LSSF == 1 .and. (.not. write_avg_urr_xs)) call read_avg_urr_xs(i)

      ! loop over orbital quantum numbers
      do i_l = 0, tope % NLS(i_ER) - 1
        read(in, 10) rec
        read(rec(1:11),  '(E11.0)') A
        read(rec(23:33),   '(I11)') L
        read(rec(56:66),   '(I11)') tope % NJS(i_l + 1)

        ! check mass ratios
        call check_mass(A, tope % AWR)

        ! check orbital quantum number
        call check_l_number(L, i_L)

        ! allocate total angular momenta
        allocate(tope % AJ(i_l + 1) % dim1(tope % NJS(i_l + 1)))

        ! allocate degrees of freedom for partial widths
        allocate(tope % DOFX(i_l + 1) % dim1(tope % NJS(i_l + 1)))
        allocate(tope % DOFN(i_l + 1) % dim1(tope % NJS(i_l + 1)))
        allocate(tope % DOFG(i_l + 1) % dim1(tope % NJS(i_l + 1)))
        allocate(tope % DOFF(i_l + 1) % dim1(tope % NJS(i_l + 1)))

        ! allocate partial widths
        allocate(tope % D_mean  (i_l + 1) % dim2(tope % NJS(i_l + 1)))
        allocate(tope % GN0_mean(i_l + 1) % dim2(tope % NJS(i_l + 1)))
        allocate(tope % GG_mean (i_l + 1) % dim2(tope % NJS(i_l + 1)))
        allocate(tope % GF_mean (i_l + 1) % dim2(tope % NJS(i_l + 1)))
        allocate(tope % GX_mean (i_l + 1) % dim2(tope % NJS(i_l + 1)))

        ! loop over total angular momenta
        do i_J = 1, tope % NJS(i_l + 1)
          allocate(tope % D_mean  (i_l + 1) % dim2(i_J) % dim1(tope % NE))
          allocate(tope % GN0_mean(i_l + 1) % dim2(i_J) % dim1(tope % NE))
          allocate(tope % GG_mean (i_l + 1) % dim2(i_J) % dim1(tope % NE))
          allocate(tope % GF_mean (i_l + 1) % dim2(i_J) % dim1(tope % NE))
          allocate(tope % GX_mean (i_l + 1) % dim2(i_J) % dim1(2))
          read(in, 10) rec
          read(rec(1:11), '(E11.0)') tope%D_mean(i_l+1)   % dim2(i_J)%dim1(1)
          read(rec(12:22),'(E11.0)') tope%AJ(i_l + 1)     % dim1(i_J)
          read(rec(23:33),'(E11.0)') tope%DOFN(i_l + 1)   % dim1(i_J)
          read(rec(34:44),'(E11.0)') tope%GN0_mean(i_l+1) % dim2(i_J)%dim1(1)
          read(rec(45:55),'(E11.0)') tope%GG_mean (i_l+1) % dim2(i_J)%dim1(1)
          tope % DOFX(i_l + 1) % dim1(i_J) = ZERO
          tope % DOFG(i_l + 1) % dim1(i_J) = ZERO
          tope % DOFF(i_l + 1) % dim1(i_J) = ZERO
          tope % D_mean(i_l + 1) % dim2(i_J) % dim1(2)&
               = tope % D_mean(i_l + 1) % dim2(i_J) % dim1(1)
          tope % GN0_mean(i_l + 1) % dim2(i_J) % dim1(2)&
               = tope % GN0_mean(i_l + 1) % dim2(i_J) % dim1(1)
          tope % GG_mean(i_l + 1) % dim2(i_J) % dim1(2)&
               = tope % GG_mean(i_l + 1) % dim2(i_J) % dim1(1)
          tope % GF_mean(i_l + 1) % dim2(i_J) % dim1(:) = ZERO
          tope % GX_mean(i_l + 1) % dim2(i_J) % dim1(:) = ZERO

        end do
      end do

    case(1)
      ! read first line of energy range subsection
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') tope % SPI(i_ER)
      read(rec(12:22), '(E11.0)') tope % AP(i_ER)
! TODO: don't overwrite the resolved value
      call tope % channel_radius(i_ER)
      read(rec(23:33), '(I11)')   tope % LSSF
      read(rec(45:55), '(I11)')   tope % NE
      read(rec(56:66), '(I11)')   tope % NLS(i_ER)
      if (tope % NLS(i_ER) > 3) &
           call fatal_error('URR parameters given for a spin sequence higher than&
           & d-wave')

      ! allocate number of total angular momenta values for each l
      allocate(tope % NJS(tope % NLS(i_ER)))

      ! allocate total angular momenta
      allocate(tope % AJ(tope % NLS(i_ER)))

      ! allocate degrees of freedom for partial widths
      allocate(tope % DOFX(tope % NLS(i_ER)))
      allocate(tope % DOFN(tope % NLS(i_ER)))
      allocate(tope % DOFG(tope % NLS(i_ER)))
      allocate(tope % DOFF(tope % NLS(i_ER)))

      ! allocate fission width energies
      allocate(tope % ES(tope % NE))
      allocate(tope % D_mean  (tope % NLS(i_ER)))
      allocate(tope % GN0_mean(tope % NLS(i_ER)))
      allocate(tope % GG_mean (tope % NLS(i_ER)))
      allocate(tope % GF_mean (tope % NLS(i_ER)))
      allocate(tope % GX_mean (tope % NLS(i_ER)))
      
      if (tope % LSSF == 1 .and. (.not. write_avg_urr_xs)) call read_avg_urr_xs(i)

      i_ES = 1
      do
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') tope % ES(i_ES)
        if (i_ES == tope % NE) exit
        i_ES = i_ES + 1
        read(rec(12:22), '(E11.0)') tope % ES(i_ES)
        if (i_ES == tope % NE) exit
        i_ES = i_ES + 1
        read(rec(23:33), '(E11.0)') tope % ES(i_ES)
        if (i_ES == tope % NE) exit
        i_ES = i_ES + 1
        read(rec(34:44), '(E11.0)') tope % ES(i_ES)
        if (i_ES == tope % NE) exit
        i_ES = i_ES + 1
        read(rec(45:55), '(E11.0)') tope % ES(i_ES)
        if (i_ES == tope % NE) exit
        i_ES = i_ES + 1
        read(rec(56:66), '(E11.0)') tope % ES(i_ES)
        if (i_ES == tope % NE) exit
        i_ES = i_ES + 1
      end do

      ! loop over orbital quantum numbers
      do i_l = 0, tope % NLS(i_ER) - 1
        read(in, 10) rec
        read(rec(1:11),  '(E11.0)') A
        read(rec(23:33),   '(I11)') L
        read(rec(45:55),   '(I11)') tope % NJS(i_l + 1)

        ! check mass ratios
        call check_mass(A, tope % AWR)

        ! check orbital quantum number
        call check_l_number(L, i_l)

        ! allocate total angular momenta
        allocate(tope % AJ(i_l + 1) % dim1(tope % NJS(i_l + 1)))
        
        ! allocate degress of freedom for partial widths
        allocate(tope % DOFX(i_l + 1) % dim1(tope % NJS(i_l + 1)))
        allocate(tope % DOFN(i_l + 1) % dim1(tope % NJS(i_l + 1)))
        allocate(tope % DOFG(i_l + 1) % dim1(tope % NJS(i_l + 1)))
        allocate(tope % DOFF(i_l + 1) % dim1(tope % NJS(i_l + 1)))

        allocate(tope % D_mean  (i_l + 1) % dim2(tope % NJS(i_l + 1)))
        allocate(tope % GN0_mean(i_l + 1) % dim2(tope % NJS(i_l + 1)))
        allocate(tope % GG_mean (i_l + 1) % dim2(tope % NJS(i_l + 1)))
        allocate(tope % GF_mean (i_l + 1) % dim2(tope % NJS(i_l + 1)))
        allocate(tope % GX_mean (i_l + 1) % dim2(tope % NJS(i_l + 1)))

        ! loop over total angular momenta
        do i_J = 1, tope % NJS(i_l + 1)
          allocate(tope % D_mean  (i_l + 1) % dim2(i_J) % dim1(tope % NE))
          allocate(tope % GN0_mean(i_l + 1) % dim2(i_J) % dim1(tope % NE))
          allocate(tope % GG_mean (i_l + 1) % dim2(i_J) % dim1(tope % NE))
          allocate(tope % GF_mean (i_l + 1) % dim2(i_J) % dim1(tope % NE))
          allocate(tope % GX_mean (i_l + 1) % dim2(i_J) % dim1(tope % NE))

          read(in, 10) rec
          read(rec(23:33), '(I11)') L
          call check_l_number(L, i_l)
          read(rec(34:44), '(E11.0)') tope % DOFF(i_l + 1) % dim1(i_J)
          read(in, 10) rec
          read(rec( 1:11), '(E11.0)') tope % D_mean(i_l + 1) % dim2(i_J) % dim1(1)
          tope % D_mean(i_l + 1) % dim2(i_J) % dim1(:)&
               = tope % D_mean(i_l + 1) % dim2(i_J) % dim1(1)
          read(rec(12:22), '(E11.0)') tope % AJ(i_l + 1) % dim1(i_J)
          read(rec(23:33), '(E11.0)') tope % DOFN(i_l + 1) % dim1(i_J)
          read(rec(34:44), '(E11.0)') tope % GN0_mean(i_l+1) % dim2(i_J) % dim1(1)
          tope % GN0_mean(i_l + 1) % dim2(i_J) % dim1(:)&
               = tope % GN0_mean(i_l + 1) % dim2(i_J) % dim1(1)
          read(rec(45:55), '(E11.0)') tope % GG_mean (i_l+1) % dim2(i_J) % dim1(1)
          tope % GG_mean(i_l + 1) % dim2(i_J) % dim1(:)&
               = tope % GG_mean(i_l + 1) % dim2(i_J) % dim1(1)
          tope % GX_mean(i_l + 1) % dim2(i_J) % dim1(:) = ZERO
          tope % DOFG(i_l + 1) % dim1(i_J) = ZERO
          tope % DOFX(i_l + 1) % dim1(i_J) = ZERO

          ! loop over energies for which data are tabulated
          i_ES = 1
          do
            read(in, 10) rec
            read(rec( 1:11), '(E11.0)') tope % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == tope % NE) exit
            i_ES = i_ES + 1
            read(rec(12:22), '(E11.0)') tope % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == tope % NE) exit
            i_ES = i_ES + 1
            read(rec(23:33), '(E11.0)') tope % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == tope % NE) exit
            i_ES = i_ES + 1
            read(rec(34:44), '(E11.0)') tope % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == tope % NE) exit
            i_ES = i_ES + 1
            read(rec(45:55), '(E11.0)') tope % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == tope % NE) exit
            i_ES = i_ES + 1
            read(rec(56:66), '(E11.0)') tope % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == tope % NE) exit
            i_ES = i_ES + 1
          end do
        end do
      end do

    case default
      call fatal_error('LFW must be 0 or 1.')

    end select

  end subroutine read_urr_slbw_parameters_lrf1

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_URR_SLBW_PARAMETERS_LRF2 reads in unresolved resonance region parameters
! for the LRF = 2 option (allow energy-dependence for all parameters)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_urr_slbw_parameters_lrf2(i, i_ER)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    character(80) :: rec ! ENDF-6 file record
    integer :: i    ! index in global isotopes array
    integer :: i_ER ! energy range index
    integer :: i_l  ! orbital quantum number index
    integer :: L    ! orbital quantum number
    integer :: i_J  ! total angular momentum quantum number index
    integer :: i_E  ! energy region index
    real(8) :: A    ! isotope/neutron mass ratio

    tope => isotopes(i)

    ! read first line of energy range subsection
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') tope % SPI(i_ER)
    read(rec(12:22), '(E11.0)') tope % AP(i_ER)
! TODO: don't overwrite the resolved value
    call tope % channel_radius(i_ER)
    read(rec(23:33), '(I11)')   tope % LSSF
    read(rec(45:55), '(I11)')   tope % NLS(i_ER)
    if (tope % NLS(i_ER) > 3) &
         call fatal_error('URR parameters given for a spin sequence higher than&
         & d-wave')

    ! allocate number of total angular momenta values for each l
    allocate(tope % NJS(tope % NLS(i_ER)))

    ! allocate total angular momentua
    allocate(tope % AJ(tope % NLS(i_ER)))

    ! allocate degrees of freedom for partial widths
    allocate(tope % DOFX(tope % NLS(i_ER)))
    allocate(tope % DOFN(tope % NLS(i_ER)))
    allocate(tope % DOFG(tope % NLS(i_ER)))
    allocate(tope % DOFF(tope % NLS(i_ER)))

    ! loop over orbital quantum numbers
    do i_l = 0, tope % NLS(i_ER) - 1
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') A
      read(rec(23:33),   '(I11)') L
      read(rec(45:55),   '(I11)') tope % NJS(i_l + 1)

      ! check mass ratios
      call check_mass(A, tope % AWR)

      ! check orbital quantum number
      call check_l_number(L, i_L)

      ! allocate total angular momenta
      allocate(tope % AJ(i_l + 1) % dim1(tope % NJS(i_l + 1)))
      
      ! allocate degress of freedom for partial widths
      allocate(tope % DOFX(i_l + 1) % dim1(tope % NJS(i_l + 1)))
      allocate(tope % DOFN(i_l + 1) % dim1(tope % NJS(i_l + 1)))
      allocate(tope % DOFG(i_l + 1) % dim1(tope % NJS(i_l + 1)))
      allocate(tope % DOFF(i_l + 1) % dim1(tope % NJS(i_l + 1)))

      ! loop over total angular momenta
      do i_J = 1, tope % NJS(i_l + 1)

        read(in, 10) rec
        read(rec(1:11), '(E11.0)') tope % AJ(i_l + 1) % dim1(i_J)
        read(rec(23:33),  '(I11)') tope % INT
        read(rec(56:66),  '(I11)') tope % NE

        ! allocate energies, mean level spacings, partial widths, 
        ! and avgeraged URR cross section values
        if (.not. (allocated(tope % ES))) then
          allocate(tope % ES(tope % NE))
          allocate(tope % D_mean  (tope % NLS(i_ER)))
          allocate(tope % GN0_mean(tope % NLS(i_ER)))
          allocate(tope % GG_mean (tope % NLS(i_ER)))
          allocate(tope % GF_mean (tope % NLS(i_ER)))
          allocate(tope % GX_mean (tope % NLS(i_ER)))
          if ((tope % LSSF) == 1 .and. (.not. write_avg_urr_xs)) then
            call read_avg_urr_xs(i)
          end if
        end if
        if (.not. (allocated(tope % D_mean(i_l + 1) % dim2))) then
          allocate(tope % D_mean  (i_l + 1) % dim2(tope % NJS(i_l + 1)))
          allocate(tope % GN0_mean(i_l + 1) % dim2(tope % NJS(i_l + 1)))
          allocate(tope % GG_mean (i_l + 1) % dim2(tope % NJS(i_l + 1)))
          allocate(tope % GF_mean (i_l + 1) % dim2(tope % NJS(i_l + 1)))
          allocate(tope % GX_mean (i_l + 1) % dim2(tope % NJS(i_l + 1)))
        end if
        allocate(tope % D_mean  (i_l + 1) % dim2(i_J) % dim1(tope % NE))
        allocate(tope % GN0_mean(i_l + 1) % dim2(i_J) % dim1(tope % NE))
        allocate(tope % GG_mean (i_l + 1) % dim2(i_J) % dim1(tope % NE))
        allocate(tope % GF_mean (i_l + 1) % dim2(i_J) % dim1(tope % NE))
        allocate(tope % GX_mean (i_l + 1) % dim2(i_J) % dim1(tope % NE))

        ! read in degrees of freedom
        read(in, 10) rec
        read(rec(23:33), '(E11.0)') tope % DOFX(i_l + 1) % dim1(i_J)
        read(rec(34:44), '(E11.0)') tope % DOFN(i_l + 1) % dim1(i_J)
        read(rec(45:55), '(E11.0)') tope % DOFG(i_l + 1) % dim1(i_J)
        read(rec(56:66), '(E11.0)') tope % DOFF(i_l + 1) % dim1(i_J)

        ! loop over energies for which data are tabulated
        do i_E = 1, tope % NE
          read(in, 10) rec
          read(rec(1:11), '(E11.0)') tope % ES(i_E)
          read(rec(12:22),'(E11.0)') tope%D_mean  (i_l+1)% dim2(i_J) % dim1(i_E)
          read(rec(23:33),'(E11.0)') tope%GX_mean (i_l+1)% dim2(i_J) % dim1(i_E)
          read(rec(34:44),'(E11.0)') tope%GN0_mean(i_l+1)% dim2(i_J) % dim1(i_E)
          read(rec(45:55),'(E11.0)') tope%GG_mean (i_l+1)% dim2(i_J) % dim1(i_E)
          read(rec(56:66),'(E11.0)') tope%GF_mean (i_l+1)% dim2(i_J) % dim1(i_E)
        end do
      end do
    end do

  end subroutine read_urr_slbw_parameters_lrf2

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_ZAID checks that the ZAID given in the ENDF-6 file is the same as that
! given in the processed ACE file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_zaid(zaid_val, zaid_ref)

    integer :: zaid_val
    integer :: zaid_ref

    if (zaid_val /= zaid_ref .and. master)&
         call fatal_error(trim(adjustl(filename))//&
         ' and the corresponding ACE file give conflicting ZAID values')

  end subroutine check_zaid

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_MASS checks that the AWR given in the ENDF-6 file is the same value as
! that given in the processed ACE file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_mass(awr_val, awr_ref)

    real(8) :: awr_val
    real(8) :: awr_ref

    if (awr_val /= awr_ref .and. master)&
         call write_message(trim(adjustl(filename))//&
         ' and the corresponding ACE file give conflicting AWR values', 6)

  end subroutine check_mass

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_Q_VALUE checks that the Q-value is 0.0
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_q_value(QX)

    real(8) :: QX

    if (QX /= ZERO .and. master)&
         call fatal_error('Q-value is not equal to 0.0 in '//trim(filename))

 end subroutine check_q_value

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_PARAMETERS checks that resonance parameters are given in MF=2
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_parameters(lrp_val)

    integer :: lrp_val

    if (master) then
      select case(lrp_val)
      case(-1)
        call fatal_error('LRP of -1 not supported in '//trim(filename))
      case(0)
       call fatal_error('LRP of 0 not supported in '//trim(filename))
      case(1)
        continue
      case(2)
        call fatal_error('LRP of 2 not supported in '//trim(filename))
      case default
        call fatal_error('LRP must be -1, 0, 1, or 2 in '//trim(filename))
      end select
    end if

  end subroutine check_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_SINGLE_ISOTOPE checks that the given ENDF-6 file contains data for only
! a single isotope
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_single_isotope(nis_val)

    integer :: nis_val

    if (nis_val /= 1 .and. master)&
         call fatal_error(trim(filename)//&
         ' contains data for more than 1 isotope')

  end subroutine check_single_isotope

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_ABUNDANCE checks that the abundance of the isotope in the ENDF-6 file is
! unity
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_abundance(abn_val)

    real(8) :: abn_val

    if (abn_val /= ONE .and. master)&
         call fatal_error('Abundance of isotope given in '//&
         trim(filename)//' is not unity')

  end subroutine check_abundance

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_FISSION_WIDTHS checks that the treatment of average fission widths in
! the URR is supported
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_fission_widths(LFW)

    integer :: LFW

    select case(LFW)
    ! average URR fission widths are not given
    case(0)
      continue
    ! average URR fission widths are given
    case(1)
      continue
    case default
      if (master) call fatal_error('LFW must be 0 or 1 in '//trim(filename))
    end select

  end subroutine check_fission_widths

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_ENERGY_RANGES makes sure the upper energy is greater than the lower and
! that the number of resonance energy ranges is allowable
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_energy_ranges(n_ranges, e_low, e_high)

    integer :: n_ranges
    real(8) :: e_low
    real(8) :: e_high

    if (n_ranges /= 2) then
      if (master)&
           call warning('> 2 resonance energy ranges (i.e. NER > 2); see '&
           //trim(filename))
    end if

    if (e_high <= e_low) then
      if (master) call fatal_error('Upper resonance energy range bound is not &
           &greater than the lower in '//trim(filename))
    end if

  end subroutine check_energy_ranges

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_RADIUS_FLAGS checks the channel and scattering radius flags
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_radius_flags(nro_val, naps_val)

    integer :: nro_val
    integer :: naps_val

    select case(nro_val)
    case(0)
      select case(naps_val)
      case(0)
        continue
      case(1)
        continue
      case default
        if (master) call fatal_error('ENDF-6 NAPS flag must be 0 or 1 when NRO&
             & is 0 in '//trim(filename))
      end select
    case(1)
      if (master) call fatal_error('Energy-dependent scattering radius in '&
           //trim(filename)//' not yet supported')
      select case(naps_val)
      case(0)
        continue
      case(1)
        continue
      case(2)
        continue
      case default
        if (master) call fatal_error('ENDF-6 NAPS flag must be 0, 1, or 2 when&
             & NRO is 1 in '//trim(filename))
      end select
    case default
      if (master) call fatal_error('ENDF-6 NRO flag must be 0 or 1 in '&
           //trim(filename))
    end select

  end subroutine check_radius_flags

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_SCATTERING_RADIUS checks that the scattering radius does not change
! within a resonance energy range
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_scattering_radius(ap_val, ap_ref)

    real(8) :: ap_val
    real(8) :: ap_ref

    if (ap_val /= ap_ref) then
      if (master)&
           call warning('AP value changes within a resonance energy range in '&
           //trim(filename))
    end if

  end subroutine check_scattering_radius

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_L_NUMBER checks for the expected ordering of orbital quantum numbers
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_l_number(l_val, l_ref)

    integer :: l_val
    integer :: l_ref

    if (l_val /= l_ref) then
      if (master) call fatal_error('Unexpected ordering of orbital quantum &
           &numbers in '//trim(filename))
    end if

  end subroutine check_l_number

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_J_SIGN checks that the signs of the total angular momenta are all
! positive
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_j_sign(j_val)

    real(8) :: j_val

    if (j_val < ZERO) then
      if (master) call warning('Negative total angular momentum &
           &in '//trim(filename))
    end if

  end subroutine check_j_sign

end module URR_endf6_parser
