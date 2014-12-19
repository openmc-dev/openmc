module endf_reader

  use ace_header,        only: Nuclide
  use avg_urr_xs_values, only: set_avg_urr_xs
  use error,             only: fatal_error, warning
  use global
  use output,            only: write_message

  implicit none

  integer :: in = 11 ! input unit
  character(80) :: filename ! ENDF-6 filename

contains

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_ENDF6 reads in an ENDF-6 format file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_endf6(filename_tmp, i_n)

    type(Nuclide), pointer :: nuc => null()
    character(80) :: filename_tmp ! temporary filename
    character(80) :: rec          ! ENDF file record
    character(7)  :: readable     ! is ENDF-6 file readable?
    logical       :: file_exists ! does ENDF-6 file exist?
    logical       :: MF1_read    ! has MF=1 been read?
    logical       :: MF2_read    ! has MF=2 been read?
    logical       :: MF3_read    ! has MF=3 been read?
    integer       :: i_n ! index in global nuclides array
    integer       :: MF  ! MF file number
    integer       :: MT  ! MT type number
    integer       :: NS  ! record number

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

    nuc => nuclides(i_n)

    MF1_read = .false.
    MF2_read = .false.
    MF3_read = .false.

    do
      read(in, 10) rec
10    format(A80)
      read(rec(67:70), '(I4)') nuc % MAT
      read(rec(71:72), '(I2)') MF
      read(rec(73:75), '(I3)') MT
      read(rec(76:80), '(I5)') NS
      if (MF == 1 .and. MT == 451 .and. (.not. MF1_read)) then
        call read_MF1(i_n, rec)
        MF1_read = .true.
      else if (MF == 2 .and. MT == 151 .and. (.not. MF2_read)) then
        call read_MF2(i_n, rec)
        MF2_read = .true.
      else if (MF == 3 .and. MT == 1 .and. (.not. MF3_read)) then
        call read_MF3(i_n, rec)
        MF3_read = .true.
      else if (nuc % MAT == -1) then
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

  subroutine read_MF1(i_n, rec)

    type(Nuclide), pointer :: nuc => null()
    character(80) :: rec                  ! ENDF-6 file record
    integer :: i_n ! index in global nuclides array
    real(8) :: ZA
    real(8) :: AWR

    nuc => nuclides(i_n)

    read(rec(1:11),  '(E11.0)') ZA
    read(rec(12:22), '(E11.0)') AWR
    read(rec(23:33), '(I11)')   nuc % LRP

    ! check ZAID agreement
    call check_zaid(int(ZA), nuc % zaid)

    ! check mass ratios
    call check_mass(AWR, nuc % awr)

    ! check that resonance parameters are given
    call check_parameters(nuc % LRP)

  end subroutine read_MF1

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_MF2 reads in an ENDF-6 format MF 2 file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_MF2(i_n, rec)

    type(Nuclide), pointer :: nuc => null()
    character(80) :: rec ! ENDF-6 file record
    integer :: i_n  ! index in global nuclides array
    integer :: i_ER ! resonance energy range index
    real(8) :: ZA
    real(8) :: ZAI
    integer :: NIS
    real(8) :: AWR
    real(8) :: ABN

    nuc => nuclides(i_n)

    read(rec(1:11),  '(E11.0)') ZA
    read(rec(12:22), '(E11.0)') AWR
    read(rec(45:55), '(I11)')   NIS

    ! check ZAID agreement
    call check_zaid(int(ZA), nuc % zaid)

    ! check that mass is consistent
    call check_mass(AWR, nuc % awr)

    ! check that this is a single-nuclide ENDF-6 file
    call check_single_nuclide(NIS)

    ! read MF=2, record=2
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') ZAI
    read(rec(12:22), '(E11.0)') ABN
    read(rec(45:55), '(I11)')   nuc % NER

    ! allocate energy range variables
    call nuc % alloc_energy_range()

    ! check ZAID agreement
    call check_zaid(int(ZAI), nuc % zaid)

    ! check abundance is unity
    call check_abundance(ABN)

    ! loop over energy ranges
    do i_ER = 1, nuc % NER

      ! read first record for this energy range
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') nuc % EL(i_ER)
      read(rec(12:22), '(E11.0)') nuc % EH(i_ER)
      read(rec(24:33),   '(I10)') nuc % LRU(i_ER)
      read(rec(35:44),   '(I10)') nuc % LRF(i_ER)
      read(rec(45:55),   '(I11)') nuc % NRO(i_ER)
      read(rec(56:66),   '(I11)') nuc % NAPS(i_ER)

      ! check number of resonance energy ranges and their bounding energies
      call check_energy_ranges(nuc % NER, nuc % EL(i_ER), nuc % EH(i_ER))

      ! check channel, scattering radius energy dependence flags
      call check_radius_flags(nuc % NRO(i_ER), nuc % NAPS(i_ER))

      ! read energy range and formalism-dependent resonance subsection data
      call read_resonance_subsection(i_n, i_ER, &
        & nuc % LRU(i_ER), nuc % LRF(i_ER))

    end do

  end subroutine read_MF2

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_MF3 reads in an ENDF-6 format MF 3 file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_MF3(i_n, rec)

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer :: i_n   ! index in global nuclides array
    integer :: i_rec ! record index
    integer :: i_e   ! index in energy grid
    integer :: MF    ! ENDF-6 MF flag
    integer :: MT    ! ENDF-6 MT flag
    integer :: NP     ! number of energy-xs pairs
    integer :: NR     ! number of interpolation regions
    integer :: NBT    ! number of entries separating interp. ranges N and N+1
    integer :: INTERP ! interp. scheme flag for an interpolation region
    character(80) :: rec ! ENDF-6 file record
    real(8) :: ZA
    real(8) :: AWR

    nuc => nuclides(i_n)

10  format(A80)

    if (nuc % LSSF == 0) then

      do
        read(in, 10) rec
        read(rec(71:72), '(I2)') MF
        read(rec(73:75), '(I3)') MT
        if (MF == 3 .and. MT == 2) exit
        if (MF == 4) call fatal_error('Reached end of MF3 w/o reading an&
          & elastic, inelastic, fission, or capture cross section in '&
          &//trim(filename))
      end do

      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') AWR

      ! check ZAID agreement
      call check_zaid(int(ZA), nuc % zaid)

      ! check that mass is consistent
      call check_mass(AWR, nuc % awr)

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

      ! check interpolation scheme is same as nuclide interpolation scheme
      call check_interp_scheme(INTERP, nuc % INT)

      allocate(nuc % avg_urr_n_e(NP))
      allocate(nuc % avg_urr_n(NP))

      i_e = 1
      do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') nuc % avg_urr_n_e(i_e)
        read(rec(12:22), '(E11.0)') nuc % avg_urr_n(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(23:33), '(E11.0)') nuc % avg_urr_n_e(i_e)
        read(rec(34:44), '(E11.0)') nuc % avg_urr_n(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(45:55), '(E11.0)') nuc % avg_urr_n_e(i_e)
        read(rec(56:66), '(E11.0)') nuc % avg_urr_n(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
      end do

      do
        read(in, 10) rec
        read(rec(71:72), '(I2)') MF
        read(rec(73:75), '(I3)') MT
        if (MF == 3 .and. MT == 4) exit
        if (MF == 4) call fatal_error('Reached end of MF3 w/o reading an&
          & elastic, inelastic, fission, or capture cross section in '&
          &//trim(filename))
      end do

      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') AWR

      ! check ZAID agreement
      call check_zaid(int(ZA), nuc % zaid)

      ! check that mass is consistent
      call check_mass(AWR, nuc % awr)

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

      ! check interpolation scheme is same as nuclide interpolation scheme
      call check_interp_scheme(INTERP, nuc % INT)

      allocate(nuc % avg_urr_x_e(NP))
      allocate(nuc % avg_urr_x(NP))

      i_e = 1
      do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') nuc % avg_urr_x_e(i_e)
        read(rec(12:22), '(E11.0)') nuc % avg_urr_x(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(23:33), '(E11.0)') nuc % avg_urr_x_e(i_e)
        read(rec(34:44), '(E11.0)') nuc % avg_urr_x(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(45:55), '(E11.0)') nuc % avg_urr_x_e(i_e)
        read(rec(56:66), '(E11.0)') nuc % avg_urr_x(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
      end do

      do
        read(in, 10) rec
        read(rec(71:72), '(I2)') MF
        read(rec(73:75), '(I3)') MT
        if (MF == 3 .and. MT == 18) exit
        if (MF == 4) call fatal_error('Reached end of MF3 w/o reading an&
          & elastic, inelastic, fission, or capture cross section in '&
          &//trim(filename))
      end do

      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') AWR

      ! check ZAID agreement
      call check_zaid(int(ZA), nuc % zaid)

      ! check that mass is consistent
      call check_mass(AWR, nuc % awr)

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

      ! check interpolation scheme is same as nuclide interpolation scheme
      call check_interp_scheme(INTERP, nuc % INT)

      allocate(nuc % avg_urr_f_e(NP))
      allocate(nuc % avg_urr_f(NP))

      i_e = 1
      do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') nuc % avg_urr_f_e(i_e)
        read(rec(12:22), '(E11.0)') nuc % avg_urr_f(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(23:33), '(E11.0)') nuc % avg_urr_f_e(i_e)
        read(rec(34:44), '(E11.0)') nuc % avg_urr_f(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(45:55), '(E11.0)') nuc % avg_urr_f_e(i_e)
        read(rec(56:66), '(E11.0)') nuc % avg_urr_f(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
      end do

      do
        read(in, 10) rec
        read(rec(71:72), '(I2)') MF
        read(rec(73:75), '(I3)') MT
        if (MF == 3 .and. MT == 102) exit
        if (MF == 4) call fatal_error('Reached end of MF3 w/o reading an&
          & elastic, inelastic, fission, or capture cross section in '&
          &//trim(filename))
      end do

      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') AWR

      ! check ZAID agreement
      call check_zaid(int(ZA), nuc % zaid)

      ! check that mass is consistent
      call check_mass(AWR, nuc % awr)

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

      ! check interpolation scheme is same as nuclide interpolation scheme
      call check_interp_scheme(INTERP, nuc % INT)

      allocate(nuc % avg_urr_g_e(NP))
      allocate(nuc % avg_urr_g(NP))

      i_e = 1
      do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') nuc % avg_urr_g_e(i_e)
        read(rec(12:22), '(E11.0)') nuc % avg_urr_g(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(23:33), '(E11.0)') nuc % avg_urr_g_e(i_e)
        read(rec(34:44), '(E11.0)') nuc % avg_urr_g(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(45:55), '(E11.0)') nuc % avg_urr_g_e(i_e)
        read(rec(56:66), '(E11.0)') nuc % avg_urr_g(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
      end do

    else
      continue
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

    if (NR /= 1) call fatal_error('More than 1 MF3 interpolation region not&
      & supported in '//trim(filename))

  end subroutine check_interp_regions

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_N_PAIRS checks that there is an allowable number of energy-xs pairs
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_n_pairs(NBT, NP)

    integer :: NBT ! number of pairs between this and the next interp. region
    integer :: NP  ! number of energy-xs pairs

    if (NBT /= NP) call fatal_error('Different NBT and NP values in MF3 for '&
      &//trim(filename))

  end subroutine check_n_pairs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_INTERP_SCHEME checks that there is an allowable interpolation scheme
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_interp_scheme(INTERP, nuc_INT)

    integer :: INTERP     ! interpolation scheme for current region
    integer :: nuc_INT ! overall nuclide interpolation scheme

    if (INTERP /= nuc_INT) call fatal_error('Different interpolation schemes for&
      & an MF3 region and the MF2 parameters in '//trim(filename))

  end subroutine check_interp_scheme

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_RESONANCE_SUBSECTION reads in the energy range and formalism-dependent
! resonance subsection data
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_resonance_subsection(i_n, i_ER, LRU, LRF)

    integer :: i_n ! index in global nuclides array
    integer :: i_ER ! energy range index
    integer :: LRU
    integer :: LRF

    ! select energy range
    select case(LRU)

    ! resolved parameters
    case(1)

      ! select formalism
      select case(LRF)

      ! SLBW
      case(1)
        call fatal_error('SLBW (LRF=1) formalism not supported in '//trim(filename))

      ! MLBW
      case(2)
        call fatal_error('MLBW (LRF=2) formalism not supported in '//trim(filename))

      ! Reich-Moore
      case(3)
        call read_rm_parameters(i_n, i_ER)

      ! Adler-Adler
      case(4)
        call fatal_error('Adler-Adler (LRF=4) formalism not supported in '&
          & //trim(filename))

      ! General R-Matrix
      case(5)
        call fatal_error('General R-Matrix (LRF=5) formalism not allowed in '&
          & //trim(filename))

      ! Hybrid R-Function
      case(6)
        call fatal_error('Hybrid R-Function (LRF=6) formalism not allowed in '&
          & //trim(filename))

      ! R-Matrix Limited
      case(7)
        call fatal_error('R-Matrix Limited (LRF=7) formalism not supported in '&
          & //trim(filename))

      ! default case
      case default
        call fatal_error('LRF must be an integer between 1 and 7, inclusive in '&
          & //trim(filename))
      end select

    ! unresolved parameters
    case(2)

      nuclides(i_n) % i_urr = i_ER

      ! select energy-dependence of parameters
      select case(LRF)

      ! only fission widths are energy-dependent
      case(1)
        call fatal_error('LRF values other than 2 are not supported in '&
          & //trim(filename))

      ! all parameters are energy-dependent
      case(2)
        call read_urr_slbw_parameters(i_n, i_ER)

      ! default case
      case default
        call fatal_error('LRF values other than 2 are not supported in '&
          & //trim(filename))
      end select

    ! default case
    case default
      call fatal_error('LRU values other than 1 or 2 are not allowed in '&
        & //trim(filename))
    end select

  end subroutine read_resonance_subsection

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_RM_PARAMETERS reads in Reich-Moore resolved resonance region
! parameters
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_rm_parameters(i_n, i_ER)

    type(Nuclide), pointer :: nuc => null()
    character(80) :: rec ! ENDF-6 file record
    integer :: i_n  ! index in global nuclides array
    integer :: i_ER ! energy range index
    integer :: i_l  ! orbital quantum number index
    integer :: L    ! orbital quantum number
    integer :: NRS  ! number of resonances for this orbital quantum number
    integer :: i_R  ! resonance index
    real(8) :: AWRI ! isotope/neutron mass ratio
    real(8) :: APL  ! l-dependent AP value

    nuc => nuclides(i_n)

    ! read first line of energy range subsection
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') nuc % SPI(i_ER)
    read(rec(12:22), '(E11.0)') nuc % AP(i_ER)
    call nuc % pre_process(i_ER)
    read(rec(45:55), '(I11)')   nuc % NLS(i_ER)

    ! allocate Reich-Moore resonance vectors for each l
    allocate(nuc % rm_resonances(nuc % NLS(i_ER)))

    ! loop over orbital quantum numbers
    do i_l = 0, nuc % NLS(i_ER) - 1
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') AWRI
      read(rec(12:22), '(E11.0)') APL
      read(rec(23:33),   '(I11)') L
      read(rec(56:66),   '(I11)') NRS

      ! allocate Reich-Moore resonances
      call nuc % rm_resonances(i_l + 1) % alloc_rm_resonances(NRS)

      ! check mass ratios
      call check_mass(AWRI, nuc % awr)

      ! check scattering radii
      call check_scattering_radius(APL, nuc % AP(i_ER))

      ! check orbital quantum number
      call check_l_number(L, i_L)

      ! loop over resonances
      do i_R = 1, NRS

        read(in, 10) rec
        read(rec(1:11),  '(E11.0)') nuc % rm_resonances(i_l + 1) % E_lam(i_R)
        read(rec(12:22), '(E11.0)') nuc % rm_resonances(i_l + 1) % AJ(i_R)
        read(rec(23:33), '(E11.0)') nuc % rm_resonances(i_l + 1) % GN(i_R)
        read(rec(34:44), '(E11.0)') nuc % rm_resonances(i_l + 1) % GG(i_R)
        read(rec(45:55), '(E11.0)') nuc % rm_resonances(i_l + 1) % GFA(i_R)
        read(rec(56:66), '(E11.0)') nuc % rm_resonances(i_l + 1) % GFB(i_R)

        ! check sign of total angular momentum, J
        call check_j_sign(nuc % rm_resonances(i_l + 1) % AJ(i_R))
      end do

      if (i_ER < nuc % NER - 1) then
        call nuc % rm_resonances(i_l + 1) % dealloc_rm_resonances()
      end if

    end do

    if (i_ER < nuc % NER - 1) deallocate(nuc % rm_resonances)

  end subroutine read_rm_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_URR_SLBW_PARAMETERS reads in unresolved resonance region parameters
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_urr_slbw_parameters(i_n, i_ER)

    type(Nuclide), pointer, save :: nuc => null()
    character(80) :: rec ! ENDF-6 file record
    integer :: i_n  ! index in global nuclides array
    integer :: i_ER ! energy range index
    integer :: i_l  ! orbital quantum number index
    integer :: L    ! orbital quantum number
    integer :: i_J  ! total angular momentum quantum number index
    integer :: i_E  ! tabulated energy grid index
    real(8) :: AWRI ! isotope/neutron mass ratio
!$omp threadprivate(nuc)

    nuc => nuclides(i_n)

    ! read first line of energy range subsection
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') nuc % SPI(i_ER)
    read(rec(12:22), '(E11.0)') nuc % AP(i_ER)
! TODO: don't overwrite the resolved value
    call nuc % pre_process(i_ER)
    read(rec(23:33), '(I11)')   nuc % LSSF
    read(rec(45:55), '(I11)')   nuc % NLS(i_ER)

    ! allocate number of total angular momenta values for each l
    allocate(nuc % NJS(nuc % NLS(i_ER)))

    ! allocate total angular momentua
    allocate(nuc % AJ(nuc % NLS(i_ER)))

    ! allocate degrees of freedom for partial widths
    allocate(nuc % DOFX(nuc % NLS(i_ER)))
    allocate(nuc % DOFN(nuc % NLS(i_ER)))
    allocate(nuc % DOFG(nuc % NLS(i_ER)))
    allocate(nuc % DOFF(nuc % NLS(i_ER)))

    ! loop over orbital quantum numbers
    do i_l = 0, nuc % NLS(i_ER) - 1
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') AWRI
      read(rec(23:33),   '(I11)') L
      read(rec(45:55),   '(I11)') nuc % NJS(i_l + 1)

      ! check mass ratios
      call check_mass(AWRI, nuc % awr)

      ! check orbital quantum number
      call check_l_number(L, i_L)

      ! allocate total angular momenta
      allocate(nuc % AJ(i_l + 1) % data(nuc % NJS(i_l + 1)))
      
      ! allocate degress of freedom for partial widths
      allocate(nuc % DOFX(i_l + 1) % data(nuc % NJS(i_l + 1)))
      allocate(nuc % DOFN(i_l + 1) % data(nuc % NJS(i_l + 1)))
      allocate(nuc % DOFG(i_l + 1) % data(nuc % NJS(i_l + 1)))
      allocate(nuc % DOFF(i_l + 1) % data(nuc % NJS(i_l + 1)))

      ! loop over total angular momenta
      do i_J = 1, nuc % NJS(i_l + 1)

        read(in, 10) rec
        read(rec(1:11), '(E11.0)') nuc % AJ(i_l + 1) % data(i_J)
        read(rec(23:33),  '(I11)') nuc % INT
        read(rec(56:66),  '(I11)') nuc % NE

        ! allocate energies, mean level spacings, partial widths, 
        ! and avgeraged URR cross section values
        if (.not. (allocated(nuc % ES))) then
          allocate(nuc % ES(nuc % NE))
          allocate(nuc % D_mean  (nuc % NLS(i_ER)))
          allocate(nuc % GN0_mean(nuc % NLS(i_ER)))
          allocate(nuc % GG_mean (nuc % NLS(i_ER)))
          allocate(nuc % GF_mean (nuc % NLS(i_ER)))
          allocate(nuc % GX_mean (nuc % NLS(i_ER)))
          if (nuc % LSSF == 1) then
            call nuc % alloc_avg_urr(nuc % NE)
            call set_avg_urr_xs(i_n)
          end if
        end if
        if (.not. (allocated(nuc % D_mean(i_l + 1) % data))) then
          allocate(nuc % D_mean  (i_l + 1) % data(nuc % NJS(i_l + 1)))
          allocate(nuc % GN0_mean(i_l + 1) % data(nuc % NJS(i_l + 1)))
          allocate(nuc % GG_mean (i_l + 1) % data(nuc % NJS(i_l + 1)))
          allocate(nuc % GF_mean (i_l + 1) % data(nuc % NJS(i_l + 1)))
          allocate(nuc % GX_mean (i_l + 1) % data(nuc % NJS(i_l + 1)))
        end if
        allocate(nuc % D_mean  (i_l + 1) % data(i_J) % data(nuc % NE))
        allocate(nuc % GN0_mean(i_l + 1) % data(i_J) % data(nuc % NE))
        allocate(nuc % GG_mean (i_l + 1) % data(i_J) % data(nuc % NE))
        allocate(nuc % GF_mean (i_l + 1) % data(i_J) % data(nuc % NE))
        allocate(nuc % GX_mean (i_l + 1) % data(i_J) % data(nuc % NE))

        ! read in degrees of freedom
        read(in, 10) rec
        read(rec(23:33), '(E11.0)') nuc % DOFX(i_l + 1) % data(i_J)
        read(rec(34:44), '(E11.0)') nuc % DOFN(i_l + 1) % data(i_J)
        read(rec(45:55), '(E11.0)') nuc % DOFG(i_l + 1) % data(i_J)
        read(rec(56:66), '(E11.0)') nuc % DOFF(i_l + 1) % data(i_J)

        ! loop over energies for which data are tabulated
        do i_E = 1, nuc % NE
          read(in, 10) rec
          read(rec(1:11),  '(E11.0)') nuc % ES(i_E)
          read(rec(12:22), '(E11.0)') nuc % D_mean  (i_l + 1) % data(i_J) % data(i_E)
          read(rec(23:33), '(E11.0)') nuc % GX_mean (i_l + 1) % data(i_J) % data(i_E)
          read(rec(34:44), '(E11.0)') nuc % GN0_mean(i_l + 1) % data(i_J) % data(i_E)
          read(rec(45:55), '(E11.0)') nuc % GG_mean (i_l + 1) % data(i_J) % data(i_E)
          read(rec(56:66), '(E11.0)') nuc % GF_mean (i_l + 1) % data(i_J) % data(i_E)
        end do
      end do
    end do

  end subroutine read_urr_slbw_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_ZAID checks that the ZAID given in the ENDF-6 file is the same as that
! given in the processed ACE file
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_zaid(zaid_val, zaid_ref)

    integer :: zaid_val
    integer :: zaid_ref

    if (zaid_val /= zaid_ref) then
      call fatal_error(trim(adjustl(filename))//' and the corresponding ACE file give&
        & conflicting ZAID values')
    end if

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

    if (awr_val /= awr_ref) then
      call write_message(trim(adjustl(filename))//' and the corresponding ACE file give&
        & conflicting AWR values', 6)
   end if

  end subroutine check_mass

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_PARAMETERS checks that resonance parameters are given in MF=2
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_parameters(lrp_val)

    integer :: lrp_val

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

  end subroutine check_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_SINGLE_NUCLIDE checks that the given ENDF-6 file contains data for only
! a single nuclide
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_single_nuclide(nis_val)

    integer :: nis_val

    if (nis_val /= 1) then
      call fatal_error(trim(filename)//' contains data for more than 1 nuclide')
    end if

  end subroutine check_single_nuclide

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHECK_ABUNDANCE checks that the abundance of the nuclide in the ENDF-6 file is
! unity
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine check_abundance(abn_val)

    real(8) :: abn_val

    if (abn_val /= ONE) then
      call fatal_error('Abundance of nuclide given in '//trim(filename)//' is not&
        & unity')
    end if

  end subroutine check_abundance

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
      if (master) then
        call warning('More than 2 resonance energy ranges (i.e. NER > 2); see &
          &'//trim(filename))
      end if
    end if

    if (e_high <= e_low) then
      call fatal_error('Upper resonance energy range bound is not greater &
        &than the lower in '//trim(filename))
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
        call fatal_error('ENDF-6 NAPS flag must be 0 or 1 when NRO is 0 in '&
          & //trim(filename))
      end select
    case(1)
      if (master) then
        call warning('ENDF-6 NRO flag value 1 ignored in energy range for ' &
          & //trim(filename))
      end if
      select case(naps_val)
      case(0)
        continue
      case(1)
        continue
      case(2)
        continue
      case default
        call fatal_error('ENDF-6 NAPS flag must be 0, 1, or 2 when NRO is 1 in '&
          & //trim(filename))
      end select
    case default
      call fatal_error('ENDF-6 NRO flag must be 0 or 1 in '//trim(filename))
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
      if (master) call warning('AP value changes within a resonance energy range in '&
        & //trim(filename))
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
      call fatal_error('Unexpected ordering of orbital quantum numbers in '&
        & //trim(filename))
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
      call fatal_error('Negative total angular momentum not supported in '&
        //trim(filename))
    end if

  end subroutine check_j_sign

end module endf_reader
