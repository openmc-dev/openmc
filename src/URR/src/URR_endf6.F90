module URR_endf6

  use URR_constants, only: SLBW,&
                           MLBW,&
                           REICH_MOORE,&
                           ADLER_ADLER,&
                           R_MATRIX,&
                           R_FUNCTION,&
                           R_MATRIX_LIM,&
                           ZERO,&
                           ONE,&
                           LINEAR_LINEAR
  use URR_error,     only: exit_status,&
                           log_message,&
                           EXIT_FAILURE,&
                           INFO,&
                           WARNING
  use URR_io,        only: read_avg_xs
  use URR_isotope,   only: isotopes
  use URR_settings,  only: write_avg_xs,&
                           path_endf_files

  implicit none
  private
  public :: read_endf6

  integer :: in = 11 ! input unit
  character(len=:), allocatable :: filename ! ENDF-6 filename


!> ENDF-6 flags
  type ENDF6Flags

     
     
  end type ENDF6Flags
  
contains


!> Read in an ENDF-6 format file
  subroutine read_endf6(filename_tmp, i)

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

    filename = trim(adjustl(filename_tmp))

    inquire(file = trim(path_endf_files)//filename, &
      & exist = file_exists, read = readable)
    if (.not. file_exists) then
      call exit_status(EXIT_FAILURE,&
           'ENDF-6 file '//trim(path_endf_files)//filename//' does not exist.')
      return
    else if (readable(1:3) == 'NO') then
      call exit_status(EXIT_FAILURE,&
           'ENDF-6 file '//trim(path_endf_files)//filename//&
           &' is not readable.  Change file permissions with chmod command.')
      return
    end if

    open(unit = in, file = trim(path_endf_files)//filename)

    MF1_read = .false.
    MF2_read = .false.
    MF3_read = .false.

    do
      read(in, 10) rec
10    format(A80)
      read(rec(67:70), '(I4)') isotopes(i) % MAT
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
      else if (isotopes(i) % MAT == -1) then
        exit
      end if
    end do

    close(in)
    deallocate(filename)

  end subroutine read_endf6


!> Read in an ENDF-6 format MF 1 file
  subroutine read_MF1(i, rec)

    character(80) :: rec ! ENDF-6 file record
    integer :: i ! isotope index
    real(8) :: real_ZAI

    read(rec(1:11),  '(E11.0)') real_ZAI
    isotopes(i) % ZAI = int(real_ZAI)
    read(rec(12:22), '(E11.0)') isotopes(i) % AWR
    read(rec(23:33), '(I11)')   isotopes(i) % LRP

! TODO:    ! check ZAID agreement
    call check_zaid(isotopes(i) % ZAI, isotopes(i) % ZAI)

! TODO:    ! check mass ratios
    call check_mass(isotopes(i) % AWR, isotopes(i) % AWR)

    ! check that resonance parameters are given
    call check_parameters(isotopes(i) % LRP)

  end subroutine read_MF1


!> Read in an ENDF-6 format MF 2 file
  subroutine read_MF2(i, rec)

    character(80) :: rec ! ENDF-6 file record
    integer :: i    ! index in global isotopes array
    integer :: i_ER ! resonance energy range index
    integer :: NIS  ! number of isotopes in material
    real(8) :: ZA
    real(8) :: A
    real(8) :: ABN

    read(rec(1:11),  '(E11.0)') ZA
    read(rec(12:22), '(E11.0)') A
    read(rec(45:55), '(I11)')   NIS

    ! check ZAID agreement
    call check_zaid(int(ZA), isotopes(i) % ZAI)

    ! check that mass is consistent
    call check_mass(A, isotopes(i) % AWR)

    ! check that this is a single-isotope ENDF-6 file
    call check_single_isotope(NIS)

    ! read MF=2, record=2
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') ZA
    read(rec(12:22), '(E11.0)') ABN
    read(rec(34:44), '(I11)')   isotopes(i) % LFW
    read(rec(45:55), '(I11)')   isotopes(i) % NER

    ! allocate energy range variables
    call isotopes(i) % alloc_energy_range()

    ! check ZAID agreement
    call check_zaid(int(ZA), isotopes(i) % ZAI)

    ! check abundance is unity
    call check_abundance(ABN)

    ! check URR average fission widths treatment
    call check_fission_widths(isotopes(i) % LFW)

    ! loop over energy ranges
    do i_ER = 1, isotopes(i) % NER

      ! read first record for this energy range
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') isotopes(i) % EL(i_ER)
      read(rec(12:22), '(E11.0)') isotopes(i) % EH(i_ER)
      read(rec(24:33),   '(I10)') isotopes(i) % LRU(i_ER)
      read(rec(35:44),   '(I10)') isotopes(i) % LRF(i_ER)
      read(rec(45:55),   '(I11)') isotopes(i) % NRO(i_ER)
      read(rec(56:66),   '(I11)') isotopes(i) % NAPS(i_ER)

      ! check number of resonance energy ranges and their bounding energies
      call check_energy_ranges(isotopes(i) % NER, isotopes(i) % EL(i_ER), isotopes(i) % EH(i_ER))

      ! check channel, scattering radius energy dependence flags
      call check_radius_flags(isotopes(i) % NRO(i_ER), isotopes(i) % NAPS(i_ER))

      ! read energy range and formalism-dependent resonance subsection data
      call read_resonance_subsection(i, i_ER, isotopes(i) % LRU(i_ER), isotopes(i) %LRF(i_ER))

    end do

  end subroutine read_MF2


!> Read in an ENDF-6 format MF 3 file
  subroutine read_MF3(i, rec)

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

10  format(A80)

    do
      read(in, 10) rec
      read(rec(71:72), '(I2)') MF
      read(rec(73:75), '(I3)') MT
      if (MF == 3 .and. MT == 2) exit
      if (MT > 2) then
        call exit_status(EXIT_FAILURE, 'Reached end of MF3 w/o reading an&
             & elastic cross section in '//trim(filename))
        return
      end if
    end do

    read(rec(1:11),  '(E11.0)') ZA
    read(rec(12:22), '(E11.0)') A

    ! check ZAID agreement
    call check_zaid(int(ZA), isotopes(i) % ZAI)

    ! check that mass is consistent
    call check_mass(A, isotopes(i) % AWR)

    ! read MF=3, record=2
    read(in, 10) rec
    read(rec(45:55), '(I11)') NR
    read(rec(56:66), '(I11)') NP

    ! check number of interpolation regions
    call check_interp_regions(NR)

    ! read MF=3, record=3
    read(in, 10) rec
    read(rec( 1:11), '(I11)') NBT
    read(rec(12:22), '(I11)') isotopes(i) % MF3_INT

    ! check number of energy-xs pairs in this interpolation region
    call check_n_pairs(NBT, NP)

    allocate(isotopes(i) % MF3_n_e(NP))
    allocate(isotopes(i) % MF3_n(NP))

    i_e = 1
    do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
      read(in, 10) rec
      read(rec( 1:11), '(E11.0)') isotopes(i) % MF3_n_e(i_e)
      read(rec(12:22), '(E11.0)') isotopes(i) % MF3_n(i_e)
      if (i_e == NP) exit
      i_e = i_e + 1
      read(rec(23:33), '(E11.0)') isotopes(i) % MF3_n_e(i_e)
      read(rec(34:44), '(E11.0)') isotopes(i) % MF3_n(i_e)
      if (i_e == NP) exit
      i_e = i_e + 1
      read(rec(45:55), '(E11.0)') isotopes(i) % MF3_n_e(i_e)
      read(rec(56:66), '(E11.0)') isotopes(i) % MF3_n(i_e)
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
        call log_message(INFO,&
             'Reached end of MF3 w/o reading a fission cross section in '//trim(filename))
        exit
      end if
    end do

    if (read_MT) then
      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') A

      ! check ZAID agreement
      call check_zaid(int(ZA), isotopes(i) % ZAI)

      ! check that mass is consistent
      call check_mass(A, isotopes(i) % AWR)

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
      call check_interp_scheme(INTERP, isotopes(i) % MF3_INT)

      allocate(isotopes(i) % MF3_f_e(NP))
      allocate(isotopes(i) % MF3_f(NP))

      i_e = 1
      do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') isotopes(i) % MF3_f_e(i_e)
        read(rec(12:22), '(E11.0)') isotopes(i) % MF3_f(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(23:33), '(E11.0)') isotopes(i) % MF3_f_e(i_e)
        read(rec(34:44), '(E11.0)') isotopes(i) % MF3_f(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(45:55), '(E11.0)') isotopes(i) % MF3_f_e(i_e)
        read(rec(56:66), '(E11.0)') isotopes(i) % MF3_f(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
      end do
    else
      isotopes(i) % MF3_f_e = [isotopes(i) % EL(isotopes(i) % i_urr),&
                               isotopes(i) % EH(isotopes(i) % i_urr)]
      isotopes(i) % MF3_f = [ZERO, ZERO]
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
        call log_message(INFO,&
             'Reached end of MF3 w/o reading an (n,n1) cross section in '//trim(filename))
        exit
      end if
    end do

    if (read_MT) then
      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') A

      ! check ZAID agreement
      call check_zaid(int(ZA), isotopes(i) % ZAI)

      ! check that mass is consistent
      call check_mass(A, isotopes(i) % AWR)

      ! read MF=3, record=2
      read(in, 10) rec
      read(rec(12:22), '(E11.0)') QI
      isotopes(i) % E_ex1 = (A + ONE) / A * (-QI)
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
      call check_interp_scheme(INTERP, isotopes(i) % MF3_INT)

      allocate(isotopes(i) % MF3_x_e(NP))
      allocate(isotopes(i) % MF3_x(NP))

      i_e = 1
      do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') isotopes(i) % MF3_x_e(i_e)
        read(rec(12:22), '(E11.0)') isotopes(i) % MF3_x(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(23:33), '(E11.0)') isotopes(i) % MF3_x_e(i_e)
        read(rec(34:44), '(E11.0)') isotopes(i) % MF3_x(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(45:55), '(E11.0)') isotopes(i) % MF3_x_e(i_e)
        read(rec(56:66), '(E11.0)') isotopes(i) % MF3_x(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
      end do
    else
      isotopes(i) % MF3_x_e = [isotopes(i) % EL(isotopes(i) % i_urr),&
                               isotopes(i) % EH(isotopes(i) % i_urr)]
      isotopes(i) % MF3_x = [ZERO, ZERO]
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
        call log_message(INFO,&
             'Reached end of MF3 w/o reading an (n,n2) Q-value in '//trim(filename))
        exit
      end if
    end do

    if (read_MT) then
      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') A

      ! check ZAID agreement
      call check_zaid(int(ZA), isotopes(i) % ZAI)

      ! check that mass is consistent
      call check_mass(A, isotopes(i) % AWR)

      ! read MF=3, record=2
      read(in, 10) rec
      read(rec(12:22), '(E11.0)') QI
      isotopes(i) % E_ex2 = (A + ONE) / A * (-QI)
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
        call log_message(INFO,&
             'Reached end of MF3 w/o reading a capture cross section in '//trim(filename))
        exit
      end if
    end do

    if (read_MT) then
      read(rec(1:11),  '(E11.0)') ZA
      read(rec(12:22), '(E11.0)') A

      ! check ZAID agreement
      call check_zaid(int(ZA), isotopes(i) % ZAI)

      ! check that mass is consistent
      call check_mass(A, isotopes(i) % AWR)

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
      call check_interp_scheme(INTERP, isotopes(i) % MF3_INT)

      allocate(isotopes(i) % MF3_g_e(NP))
      allocate(isotopes(i) % MF3_g(NP))

      i_e = 1
      do i_rec = 1, int(NP / 3) + int(ceiling(dble(mod(NP, 3)) / 3))
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') isotopes(i) % MF3_g_e(i_e)
        read(rec(12:22), '(E11.0)') isotopes(i) % MF3_g(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(23:33), '(E11.0)') isotopes(i) % MF3_g_e(i_e)
        read(rec(34:44), '(E11.0)') isotopes(i) % MF3_g(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
        read(rec(45:55), '(E11.0)') isotopes(i) % MF3_g_e(i_e)
        read(rec(56:66), '(E11.0)') isotopes(i) % MF3_g(i_e)
        if (i_e == NP) exit
        i_e = i_e + 1
      end do
    end if

  end subroutine read_MF3


!> Check that there is an allowable number of interpolation regions
  subroutine check_interp_regions(NR)

    integer :: NR ! number of interpolation regions

    if (NR /= 1) then
      call log_message(INFO,&
           'More than 1 File 3 interpolation region in '//trim(filename))
    end if

  end subroutine check_interp_regions


!> Check that there is an allowable number of energy-xs pairs
  subroutine check_n_pairs(NBT, NP)

    integer :: NBT ! number of pairs between this and the next interp. region
    integer :: NP  ! number of energy-xs pairs

    if (NBT /= NP) then
      call log_message(WARNING,&
           'Different NBT and NP values in File 3 for '//trim(filename))
    end if

  end subroutine check_n_pairs


!> Check that there is an allowable interpolation scheme
  subroutine check_interp_scheme(INTERP, MF3_INT)

    integer :: INTERP     ! interpolation scheme for current region
    integer :: MF3_INT ! overall isotope interpolation scheme

    if (INTERP /= MF3_INT) call log_message(WARNING,&
         'Different interpolation schemes for different File 3 regions in '&
         //trim(filename)//'.  Using region 1 scheme for all regions.')

  end subroutine check_interp_scheme


!> Read in the energy range and formalism-dependent resonance subsection data
  subroutine read_resonance_subsection(i, i_ER, LRU, LRF)

    integer :: i    ! index in global isotopes array
    integer :: i_ER ! energy range index
    integer :: LRU  ! resolved (1) or unresolved (2) parameters
    integer :: LRF  ! ENDF-6 resonance formalism flag

    ! select energy range
    select case(LRU)

    ! only the scattering radius is given (LRF=0; NLS=0; LFW=0)
    case(0)

      call exit_status(EXIT_FAILURE,&
           'Scattering radius only (LRU = 0) not supported in '//trim(filename))
      return

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
        call exit_status(EXIT_FAILURE,&
             'Adler-Adler (LRF=4) formalism not supported in '//trim(filename))
        return

      case(R_MATRIX)
        call exit_status(EXIT_FAILURE,&
             'General R-Matrix (LRF=5) formalism not allowed in '//trim(filename))
        return

      case(R_FUNCTION)
        call exit_status(EXIT_FAILURE,&
             'Hybrid R-Function (LRF=6) formalism not allowed in '//trim(filename))
        return

      case(R_MATRIX_LIM)
        call exit_status(EXIT_FAILURE,&
             'R-Matrix Limited (LRF=7) formalism not supported in '//trim(filename))
        return

      ! default case
      case default
        call exit_status(EXIT_FAILURE,&
             'LRF must be an integer between 1 and 7, inclusive in '//trim(filename))
        return

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

      case default
         call exit_status(EXIT_FAILURE,&
              'LRF values other than 1 and 2 are not supported in '//trim(filename))
         return

       end select

    ! default case
    case default
      call exit_status(EXIT_FAILURE,&
           'LRU values other than 0, 1 or 2 are not allowed in '//trim(filename))
      return

    end select

  end subroutine read_resonance_subsection


!> Read in Single-level Breit-Wigner resolved resonance region parameters
  subroutine read_slbw_parameters(i, i_ER)

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

    ! read first line of energy range subsection
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') isotopes(i) % SPI(i_ER)
    read(rec(12:22), '(E11.0)') isotopes(i) % AP(i_ER)
    call isotopes(i) % channel_radius(i_ER)
    read(rec(45:55), '(I11)')   isotopes(i) % NLS(i_ER)
    if (isotopes(i) % NLS(i_ER) > 3) then
      call exit_status(EXIT_FAILURE,&
           'SLBW parameters given for a resonance higher than d-wave')
      return
    end if

    ! allocate SLBW resonance vectors for each l
    allocate(isotopes(i) % bw_resonances(isotopes(i) % NLS(i_ER)))

    ! loop over orbital quantum numbers
    do i_l = 0, isotopes(i) % NLS(i_ER) - 1
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') A
      read(rec(12:22), '(E11.0)') QX
      read(rec(23:33),   '(I11)') L
      read(rec(34:44),   '(I11)') LRX
      read(rec(56:66),   '(I11)') NRS

      ! allocate SLBW resonances
      call isotopes(i) % bw_resonances(i_l + 1) % alloc(NRS)

      ! check mass ratios
      call check_mass(A, isotopes(i) % AWR)

      ! check that Q-value is 0.0
      call check_q_value(QX)

      ! check orbital quantum number
      call check_l_number(L, i_L)

      ! loop over resonances
      do i_R = 1, NRS

        read(in, 10) rec
        read(rec(1:11),  '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % E_lam
        read(rec(12:22), '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % AJ
        read(rec(23:33), '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GT
        read(rec(34:44), '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GN
        read(rec(45:55), '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GG
        read(rec(56:66), '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GF

        if (LRX == 0) then
          isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GX = ZERO
        else if (LRX == 1) then
          isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GX&
               = isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GT&
               - isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GN&
               - isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GG&
               - isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GF
        else
          call exit_status(EXIT_FAILURE, 'LRX must be 0 or 1 in '//trim(filename))
          return
        end if

        ! check sign of total angular momentum, J
        call check_j_sign(isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % AJ)
      end do

      if (i_ER < isotopes(i) % NER - 1) then
        call isotopes(i) % bw_resonances(i_l + 1) % dealloc()
      end if

    end do

    if (i_ER < isotopes(i) % NER - 1) deallocate(isotopes(i) % bw_resonances)

  end subroutine read_slbw_parameters


!> Read in Multi-level Breit-Wigner resolved resonance region parameters
  subroutine read_mlbw_parameters(i, i_ER)

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

    ! read first line of energy range subsection
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') isotopes(i) % SPI(i_ER)
    read(rec(12:22), '(E11.0)') isotopes(i) % AP(i_ER)
    call isotopes(i) % channel_radius(i_ER)
    read(rec(45:55), '(I11)')   isotopes(i) % NLS(i_ER)
    if (isotopes(i) % NLS(i_ER) > 3) then
      call exit_status(EXIT_FAILURE,&
           'MLBW parameters given for a resonance higher than d-wave')
      return
    end if

    ! allocate MLBW resonance vectors for each l
    allocate(isotopes(i) % bw_resonances(isotopes(i) % NLS(i_ER)))

    ! loop over orbital quantum numbers
    do i_l = 0, isotopes(i) % NLS(i_ER) - 1
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') A
      read(rec(12:22), '(E11.0)') QX
      read(rec(23:33),   '(I11)') L
      read(rec(34:44),   '(I11)') LRX
      read(rec(56:66),   '(I11)') NRS

      ! allocate MLBW resonances
      call isotopes(i) % bw_resonances(i_l + 1) % alloc(NRS)

      ! check mass ratios
      call check_mass(A, isotopes(i) % AWR)

      ! check that Q-value is 0.0
      call check_q_value(QX)

      ! check orbital quantum number
      call check_l_number(L, i_L)

      ! loop over resonances
      do i_R = 1, NRS

        read(in, 10) rec
        read(rec(1:11),  '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % E_lam
        read(rec(12:22), '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % AJ
        read(rec(23:33), '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GT
        read(rec(34:44), '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GN
        read(rec(45:55), '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GG
        read(rec(56:66), '(E11.0)')&
             isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GF

        if (LRX == 0) then
          isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GX = ZERO
        else if (LRX == 1) then
          isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GX&
               = isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GT&
               - isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GN&
               - isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GG&
               - isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % GF
        else
          call exit_status(EXIT_FAILURE, 'LRX must be 0 or 1 in '//trim(filename))
          return
        end if

        ! check sign of total angular momentum, J
        call check_j_sign(isotopes(i) % bw_resonances(i_l + 1) % res(i_R) % AJ)
      end do

      if (i_ER < isotopes(i) % NER - 1) then
        call isotopes(i) % bw_resonances(i_l + 1) % dealloc()
      end if
    end do

    if (i_ER < isotopes(i) % NER - 1) deallocate(isotopes(i) % bw_resonances)

  end subroutine read_mlbw_parameters


!> Read in Reich-Moore resolved resonance region parameters
  subroutine read_rm_parameters(i, i_ER)

    character(80) :: rec ! ENDF-6 file record
    integer :: i    ! index in global isotopes array
    integer :: i_ER ! energy range index
    integer :: i_l  ! orbital quantum number index
    integer :: L    ! orbital quantum number
    integer :: NRS  ! number of resonances for this orbital quantum number
    integer :: i_R  ! resonance index
    real(8) :: A    ! isotope/neutron mass ratio
    real(8) :: APL  ! l-dependent AP value

    ! read first line of energy range subsection
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') isotopes(i) % SPI(i_ER)
    read(rec(12:22), '(E11.0)') isotopes(i) % AP(i_ER)
    call isotopes(i) % channel_radius(i_ER)
    read(rec(45:55), '(I11)')   isotopes(i) % NLS(i_ER)
    if (isotopes(i) % NLS(i_ER) > 3) &
         call log_message(WARNING,&
         'R-M parameters given for a resonance higher than d-wave')

    ! allocate Reich-Moore resonance vectors for each l
    allocate(isotopes(i) % rm_resonances(isotopes(i) % NLS(i_ER)))

    ! loop over orbital quantum numbers
    do i_l = 0, isotopes(i) % NLS(i_ER) - 1
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') A
      read(rec(12:22), '(E11.0)') APL
      read(rec(23:33),   '(I11)') L
      read(rec(56:66),   '(I11)') NRS

      ! allocate Reich-Moore resonances
      call isotopes(i) % rm_resonances(i_l + 1) % alloc(NRS)

      ! check mass ratios
      call check_mass(A, isotopes(i) % AWR)

      ! check scattering radii
      call check_scattering_radius(APL, isotopes(i) % AP(i_ER))

      ! check orbital quantum number
      call check_l_number(L, i_L)

      ! loop over resonances
      do i_R = 1, NRS

        read(in, 10) rec
        read(rec(1:11),  '(E11.0)')&
             isotopes(i) % rm_resonances(i_l + 1) % res(i_R) % E_lam
        read(rec(12:22), '(E11.0)')&
             isotopes(i) % rm_resonances(i_l + 1) % res(i_R) % AJ
        read(rec(23:33), '(E11.0)')&
             isotopes(i) % rm_resonances(i_l + 1) % res(i_R) % GN
        read(rec(34:44), '(E11.0)')&
             isotopes(i) % rm_resonances(i_l + 1) % res(i_R) % GG
        read(rec(45:55), '(E11.0)')&
             isotopes(i) % rm_resonances(i_l + 1) % res(i_R) % GFA
        read(rec(56:66), '(E11.0)')&
             isotopes(i) % rm_resonances(i_l + 1) % res(i_R) % GFB

        ! check sign of total angular momentum, J
        call check_j_sign(isotopes(i) % rm_resonances(i_l + 1) % res(i_R) % AJ)
      end do

      if (i_ER < isotopes(i) % NER - 1) then
        call isotopes(i) % rm_resonances(i_l + 1) % dealloc()
      end if

    end do

    if (i_ER < isotopes(i) % NER - 1) deallocate(isotopes(i) % rm_resonances)

  end subroutine read_rm_parameters


!> Read in unresolved resonance region parameters for the LRF = 1 option
!! (only fission widths are energy-dependent)
  subroutine read_urr_slbw_parameters_lrf1(i, i_ER)

    character(80) :: rec ! ENDF-6 file record
    integer :: i    ! index in global isotopes array
    integer :: i_ER ! energy range index
    integer :: i_l  ! orbital quantum number index
    integer :: L    ! orbital quantum number
    integer :: i_J  ! total angular momentum quantum number index
    integer :: i_ES ! tabulated fission width energy grid index
    real(8) :: A    ! isotope/neutron mass ratio

    ! this is forced by the ENDF-6 format when LRF = 1 and LFW = 1
    isotopes(i) % INT = LINEAR_LINEAR

10  format(A80)
    select case(isotopes(i) % LFW)
    case(0)
      ! read first line of energy range subsection
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') isotopes(i) % SPI(i_ER)
      read(rec(12:22), '(E11.0)') isotopes(i) % AP(i_ER)
! TODO: don't overwrite the resolved value
      call isotopes(i) % channel_radius(i_ER)
      read(rec(23:33), '(I11)')   isotopes(i) % LSSF
      read(rec(45:55), '(I11)')   isotopes(i) % NLS(i_ER)
      if (isotopes(i) % NLS(i_ER) > 3) then
        call exit_status(EXIT_FAILURE,&
             'URR parameters given for a spin sequence higher than d-wave')
        return
      end if
      
      ! allocate number of total angular momenta values for each l
      allocate(isotopes(i) % NJS(isotopes(i) % NLS(i_ER)))

      ! allocate total angular momenta
      allocate(isotopes(i) % AJ(isotopes(i) % NLS(i_ER)))

      ! allocate degrees of freedom for partial widths
      allocate(isotopes(i) % DOFX(isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % DOFN(isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % DOFG(isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % DOFF(isotopes(i) % NLS(i_ER)))

      ! allocate partial widths
      isotopes(i) % NE = 2
      allocate(isotopes(i) % ES(isotopes(i) % NE))
      isotopes(i) % ES(1) = isotopes(i) % EL(i_ER)
      isotopes(i) % ES(2) = isotopes(i) % EH(i_ER)
      allocate(isotopes(i) % D_mean  (isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % GN0_mean(isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % GG_mean (isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % GF_mean (isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % GX_mean (isotopes(i) % NLS(i_ER)))

      if (isotopes(i)%LSSF == 1 .and. (.not. write_avg_xs)) call read_avg_xs(i)

      ! loop over orbital quantum numbers
      do i_l = 0, isotopes(i) % NLS(i_ER) - 1
        read(in, 10) rec
        read(rec(1:11),  '(E11.0)') A
        read(rec(23:33),   '(I11)') L
        read(rec(56:66),   '(I11)') isotopes(i) % NJS(i_l + 1)

        ! check mass ratios
        call check_mass(A, isotopes(i) % AWR)

        ! check orbital quantum number
        call check_l_number(L, i_L)

        ! allocate total angular momenta
        allocate(isotopes(i) % AJ(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))

        ! allocate degrees of freedom for partial widths
        allocate(isotopes(i) % DOFX(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % DOFN(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % DOFG(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % DOFF(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))

        ! allocate partial widths
        allocate(isotopes(i) % D_mean  (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % GN0_mean(i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % GG_mean (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % GF_mean (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % GX_mean (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))

        ! loop over total angular momenta
        do i_J = 1, isotopes(i) % NJS(i_l + 1)
          allocate(isotopes(i) % D_mean  (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
          allocate(isotopes(i) % GN0_mean(i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
          allocate(isotopes(i) % GG_mean (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
          allocate(isotopes(i) % GF_mean (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
          allocate(isotopes(i) % GX_mean (i_l + 1) % dim2(i_J) % dim1(2))
          read(in, 10) rec
          read(rec(1:11), '(E11.0)') isotopes(i)%D_mean(i_l+1)   % dim2(i_J)%dim1(1)
          read(rec(12:22),'(E11.0)') isotopes(i)%AJ(i_l + 1)     % dim1(i_J)
          read(rec(23:33),'(E11.0)') isotopes(i)%DOFN(i_l + 1)   % dim1(i_J)
          read(rec(34:44),'(E11.0)') isotopes(i)%GN0_mean(i_l+1) % dim2(i_J)%dim1(1)
          read(rec(45:55),'(E11.0)') isotopes(i)%GG_mean (i_l+1) % dim2(i_J)%dim1(1)
          isotopes(i) % DOFX(i_l + 1) % dim1(i_J) = ZERO
          isotopes(i) % DOFG(i_l + 1) % dim1(i_J) = ZERO
          isotopes(i) % DOFF(i_l + 1) % dim1(i_J) = ZERO
          isotopes(i) % D_mean(i_l + 1) % dim2(i_J) % dim1(2)&
               = isotopes(i) % D_mean(i_l + 1) % dim2(i_J) % dim1(1)
          isotopes(i) % GN0_mean(i_l + 1) % dim2(i_J) % dim1(2)&
               = isotopes(i) % GN0_mean(i_l + 1) % dim2(i_J) % dim1(1)
          isotopes(i) % GG_mean(i_l + 1) % dim2(i_J) % dim1(2)&
               = isotopes(i) % GG_mean(i_l + 1) % dim2(i_J) % dim1(1)
          isotopes(i) % GF_mean(i_l + 1) % dim2(i_J) % dim1(:) = ZERO
          isotopes(i) % GX_mean(i_l + 1) % dim2(i_J) % dim1(:) = ZERO

        end do
      end do

    case(1)
      ! read first line of energy range subsection
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') isotopes(i) % SPI(i_ER)
      read(rec(12:22), '(E11.0)') isotopes(i) % AP(i_ER)
! TODO: don't overwrite the resolved value
      call isotopes(i) % channel_radius(i_ER)
      read(rec(23:33), '(I11)')   isotopes(i) % LSSF
      read(rec(45:55), '(I11)')   isotopes(i) % NE
      read(rec(56:66), '(I11)')   isotopes(i) % NLS(i_ER)
      if (isotopes(i) % NLS(i_ER) > 3) then
        call exit_status(EXIT_FAILURE,&
             'URR parameters given for a spin sequence higher than d-wave')
        return
      end if

      ! allocate number of total angular momenta values for each l
      allocate(isotopes(i) % NJS(isotopes(i) % NLS(i_ER)))

      ! allocate total angular momenta
      allocate(isotopes(i) % AJ(isotopes(i) % NLS(i_ER)))

      ! allocate degrees of freedom for partial widths
      allocate(isotopes(i) % DOFX(isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % DOFN(isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % DOFG(isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % DOFF(isotopes(i) % NLS(i_ER)))

      ! allocate fission width energies
      allocate(isotopes(i) % ES(isotopes(i) % NE))
      allocate(isotopes(i) % D_mean  (isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % GN0_mean(isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % GG_mean (isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % GF_mean (isotopes(i) % NLS(i_ER)))
      allocate(isotopes(i) % GX_mean (isotopes(i) % NLS(i_ER)))
      
      if (isotopes(i) % LSSF == 1 .and. (.not. write_avg_xs)) call read_avg_xs(i)

      i_ES = 1
      do
        read(in, 10) rec
        read(rec( 1:11), '(E11.0)') isotopes(i) % ES(i_ES)
        if (i_ES == isotopes(i) % NE) exit
        i_ES = i_ES + 1
        read(rec(12:22), '(E11.0)') isotopes(i) % ES(i_ES)
        if (i_ES == isotopes(i) % NE) exit
        i_ES = i_ES + 1
        read(rec(23:33), '(E11.0)') isotopes(i) % ES(i_ES)
        if (i_ES == isotopes(i) % NE) exit
        i_ES = i_ES + 1
        read(rec(34:44), '(E11.0)') isotopes(i) % ES(i_ES)
        if (i_ES == isotopes(i) % NE) exit
        i_ES = i_ES + 1
        read(rec(45:55), '(E11.0)') isotopes(i) % ES(i_ES)
        if (i_ES == isotopes(i) % NE) exit
        i_ES = i_ES + 1
        read(rec(56:66), '(E11.0)') isotopes(i) % ES(i_ES)
        if (i_ES == isotopes(i) % NE) exit
        i_ES = i_ES + 1
      end do

      ! loop over orbital quantum numbers
      do i_l = 0, isotopes(i) % NLS(i_ER) - 1
        read(in, 10) rec
        read(rec(1:11),  '(E11.0)') A
        read(rec(23:33),   '(I11)') L
        read(rec(45:55),   '(I11)') isotopes(i) % NJS(i_l + 1)

        ! check mass ratios
        call check_mass(A, isotopes(i) % AWR)

        ! check orbital quantum number
        call check_l_number(L, i_l)

        ! allocate total angular momenta
        allocate(isotopes(i) % AJ(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
        
        ! allocate degress of freedom for partial widths
        allocate(isotopes(i) % DOFX(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % DOFN(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % DOFG(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % DOFF(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))

        allocate(isotopes(i) % D_mean  (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % GN0_mean(i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % GG_mean (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % GF_mean (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
        allocate(isotopes(i) % GX_mean (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))

        ! loop over total angular momenta
        do i_J = 1, isotopes(i) % NJS(i_l + 1)
          allocate(isotopes(i) % D_mean  (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
          allocate(isotopes(i) % GN0_mean(i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
          allocate(isotopes(i) % GG_mean (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
          allocate(isotopes(i) % GF_mean (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
          allocate(isotopes(i) % GX_mean (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))

          read(in, 10) rec
          read(rec(23:33), '(I11)') L
          call check_l_number(L, i_l)
          read(rec(34:44), '(E11.0)') isotopes(i) % DOFF(i_l + 1) % dim1(i_J)
          read(in, 10) rec
          read(rec( 1:11), '(E11.0)') isotopes(i) % D_mean(i_l + 1) % dim2(i_J) % dim1(1)
          isotopes(i) % D_mean(i_l + 1) % dim2(i_J) % dim1(:)&
               = isotopes(i) % D_mean(i_l + 1) % dim2(i_J) % dim1(1)
          read(rec(12:22), '(E11.0)') isotopes(i) % AJ(i_l + 1) % dim1(i_J)
          read(rec(23:33), '(E11.0)') isotopes(i) % DOFN(i_l + 1) % dim1(i_J)
          read(rec(34:44), '(E11.0)') isotopes(i) % GN0_mean(i_l+1) % dim2(i_J) % dim1(1)
          isotopes(i) % GN0_mean(i_l + 1) % dim2(i_J) % dim1(:)&
               = isotopes(i) % GN0_mean(i_l + 1) % dim2(i_J) % dim1(1)
          read(rec(45:55), '(E11.0)') isotopes(i) % GG_mean (i_l+1) % dim2(i_J) % dim1(1)
          isotopes(i) % GG_mean(i_l + 1) % dim2(i_J) % dim1(:)&
               = isotopes(i) % GG_mean(i_l + 1) % dim2(i_J) % dim1(1)
          isotopes(i) % GX_mean(i_l + 1) % dim2(i_J) % dim1(:) = ZERO
          isotopes(i) % DOFG(i_l + 1) % dim1(i_J) = ZERO
          isotopes(i) % DOFX(i_l + 1) % dim1(i_J) = ZERO

          ! loop over energies for which data are tabulated
          i_ES = 1
          do
            read(in, 10) rec
            read(rec( 1:11), '(E11.0)') isotopes(i) % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == isotopes(i) % NE) exit
            i_ES = i_ES + 1
            read(rec(12:22), '(E11.0)') isotopes(i) % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == isotopes(i) % NE) exit
            i_ES = i_ES + 1
            read(rec(23:33), '(E11.0)') isotopes(i) % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == isotopes(i) % NE) exit
            i_ES = i_ES + 1
            read(rec(34:44), '(E11.0)') isotopes(i) % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == isotopes(i) % NE) exit
            i_ES = i_ES + 1
            read(rec(45:55), '(E11.0)') isotopes(i) % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == isotopes(i) % NE) exit
            i_ES = i_ES + 1
            read(rec(56:66), '(E11.0)') isotopes(i) % GF_mean(i_l+1)%dim2(i_J)%dim1(i_ES)
            if (i_ES == isotopes(i) % NE) exit
            i_ES = i_ES + 1
          end do
        end do
      end do

    case default
      call exit_status(EXIT_FAILURE, 'LFW must be 0 or 1.')
      return

    end select

  end subroutine read_urr_slbw_parameters_lrf1


!> Read in unresolved resonance region parameters for the LRF = 2 option
!! (allow energy-dependence for all parameters)
  subroutine read_urr_slbw_parameters_lrf2(i, i_ER)

    character(80) :: rec ! ENDF-6 file record
    integer :: i    ! index in global isotopes array
    integer :: i_ER ! energy range index
    integer :: i_l  ! orbital quantum number index
    integer :: L    ! orbital quantum number
    integer :: i_J  ! total angular momentum quantum number index
    integer :: i_E  ! energy region index
    real(8) :: A    ! isotope/neutron mass ratio

    ! read first line of energy range subsection
    read(in, 10) rec
10  format(A80)
    read(rec(1:11),  '(E11.0)') isotopes(i) % SPI(i_ER)
    read(rec(12:22), '(E11.0)') isotopes(i) % AP(i_ER)
! TODO: don't overwrite the resolved value
    call isotopes(i) % channel_radius(i_ER)
    read(rec(23:33), '(I11)')   isotopes(i) % LSSF
    read(rec(45:55), '(I11)')   isotopes(i) % NLS(i_ER)
    if (isotopes(i) % NLS(i_ER) > 3) then
      call exit_status(EXIT_FAILURE,&
           'URR parameters given for a spin sequence higher than d-wave')
      return
    end if

    ! allocate number of total angular momenta values for each l
    allocate(isotopes(i) % NJS(isotopes(i) % NLS(i_ER)))

    ! allocate total angular momentua
    allocate(isotopes(i) % AJ(isotopes(i) % NLS(i_ER)))

    ! allocate degrees of freedom for partial widths
    allocate(isotopes(i) % DOFX(isotopes(i) % NLS(i_ER)))
    allocate(isotopes(i) % DOFN(isotopes(i) % NLS(i_ER)))
    allocate(isotopes(i) % DOFG(isotopes(i) % NLS(i_ER)))
    allocate(isotopes(i) % DOFF(isotopes(i) % NLS(i_ER)))

    ! loop over orbital quantum numbers
    do i_l = 0, isotopes(i) % NLS(i_ER) - 1
      read(in, 10) rec
      read(rec(1:11),  '(E11.0)') A
      read(rec(23:33),   '(I11)') L
      read(rec(45:55),   '(I11)') isotopes(i) % NJS(i_l + 1)

      ! check mass ratios
      call check_mass(A, isotopes(i) % AWR)

      ! check orbital quantum number
      call check_l_number(L, i_L)

      ! allocate total angular momenta
      allocate(isotopes(i) % AJ(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
      
      ! allocate degress of freedom for partial widths
      allocate(isotopes(i) % DOFX(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
      allocate(isotopes(i) % DOFN(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
      allocate(isotopes(i) % DOFG(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))
      allocate(isotopes(i) % DOFF(i_l + 1) % dim1(isotopes(i) % NJS(i_l + 1)))

      ! loop over total angular momenta
      do i_J = 1, isotopes(i) % NJS(i_l + 1)

        read(in, 10) rec
        read(rec(1:11), '(E11.0)') isotopes(i) % AJ(i_l + 1) % dim1(i_J)
        read(rec(23:33),  '(I11)') isotopes(i) % INT
        read(rec(56:66),  '(I11)') isotopes(i) % NE

        ! allocate energies, mean level spacings, partial widths, 
        ! and avgeraged URR cross section values
        if (.not. (allocated(isotopes(i) % ES))) then
          allocate(isotopes(i) % ES(isotopes(i) % NE))
          allocate(isotopes(i) % D_mean  (isotopes(i) % NLS(i_ER)))
          allocate(isotopes(i) % GN0_mean(isotopes(i) % NLS(i_ER)))
          allocate(isotopes(i) % GG_mean (isotopes(i) % NLS(i_ER)))
          allocate(isotopes(i) % GF_mean (isotopes(i) % NLS(i_ER)))
          allocate(isotopes(i) % GX_mean (isotopes(i) % NLS(i_ER)))
          if ((isotopes(i) % LSSF) == 1 .and. (.not. write_avg_xs)) then
            call read_avg_xs(i)
          end if
        end if
        if (.not. (allocated(isotopes(i) % D_mean(i_l + 1) % dim2))) then
          allocate(isotopes(i) % D_mean  (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
          allocate(isotopes(i) % GN0_mean(i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
          allocate(isotopes(i) % GG_mean (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
          allocate(isotopes(i) % GF_mean (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
          allocate(isotopes(i) % GX_mean (i_l + 1) % dim2(isotopes(i) % NJS(i_l + 1)))
        end if
        allocate(isotopes(i) % D_mean  (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
        allocate(isotopes(i) % GN0_mean(i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
        allocate(isotopes(i) % GG_mean (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
        allocate(isotopes(i) % GF_mean (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))
        allocate(isotopes(i) % GX_mean (i_l + 1) % dim2(i_J) % dim1(isotopes(i) % NE))

        ! read in degrees of freedom
        read(in, 10) rec
        read(rec(23:33), '(E11.0)') isotopes(i) % DOFX(i_l + 1) % dim1(i_J)
        read(rec(34:44), '(E11.0)') isotopes(i) % DOFN(i_l + 1) % dim1(i_J)
        read(rec(45:55), '(E11.0)') isotopes(i) % DOFG(i_l + 1) % dim1(i_J)
        read(rec(56:66), '(E11.0)') isotopes(i) % DOFF(i_l + 1) % dim1(i_J)

        ! loop over energies for which data are tabulated
        do i_E = 1, isotopes(i) % NE
          read(in, 10) rec
          read(rec(1:11), '(E11.0)') isotopes(i) % ES(i_E)
          read(rec(12:22),'(E11.0)') isotopes(i)%D_mean  (i_l+1)% dim2(i_J) % dim1(i_E)
          read(rec(23:33),'(E11.0)') isotopes(i)%GX_mean (i_l+1)% dim2(i_J) % dim1(i_E)
          read(rec(34:44),'(E11.0)') isotopes(i)%GN0_mean(i_l+1)% dim2(i_J) % dim1(i_E)
          read(rec(45:55),'(E11.0)') isotopes(i)%GG_mean (i_l+1)% dim2(i_J) % dim1(i_E)
          read(rec(56:66),'(E11.0)') isotopes(i)%GF_mean (i_l+1)% dim2(i_J) % dim1(i_E)
        end do
      end do
    end do

  end subroutine read_urr_slbw_parameters_lrf2


!> Check that the ZAID given in the ENDF-6 file is the same as that
!! given in the processed ACE file
  subroutine check_zaid(zaid_val, zaid_ref)

    integer :: zaid_val
    integer :: zaid_ref

    if (zaid_val /= zaid_ref) then
      call exit_status(EXIT_FAILURE, trim(adjustl(filename))//&
           ' and the corresponding ACE file give conflicting ZAID values')
      return
    end if

  end subroutine check_zaid


!> Check that the AWR values are consistent within an ENDF-6 file
  subroutine check_mass(awr_val, awr_ref)

    real(8) :: awr_val
    real(8) :: awr_ref

    if (awr_val /= awr_ref) then
      call log_message(WARNING,&
           trim(adjustl(filename))//' contains inconsistent AWR values')
    end if

  end subroutine check_mass


!> Check that the Q-value is 0.0
  subroutine check_q_value(QX)

    real(8) :: QX

    if (QX /= ZERO) then
      call exit_status(EXIT_FAILURE,&
           'Q-value is not equal to 0.0 in '//trim(filename))
      return
    end if

 end subroutine check_q_value


!> Check that resonance parameters are given in MF=2
  subroutine check_parameters(lrp_val)

    integer :: lrp_val

    select case(lrp_val)
    case(-1)
      call exit_status(EXIT_FAILURE, 'LRP of -1 not supported in '//trim(filename))
      return
    case(0)
      call exit_status(EXIT_FAILURE, 'LRP of 0 not supported in '//trim(filename))
      return
    case(1)
      continue
    case(2)
      call exit_status(EXIT_FAILURE, 'LRP of 2 not supported in '//trim(filename))
      return
    case default
      call exit_status(EXIT_FAILURE, 'LRP must be -1, 0, 1, or 2 in '//trim(filename))
      return
    end select

  end subroutine check_parameters


!> Check that the given ENDF-6 file contains data for only a single isotope
  subroutine check_single_isotope(nis_val)

    integer :: nis_val

    if (nis_val /= 1) then
      call exit_status(EXIT_FAILURE,&
           trim(filename)//' contains data for more than 1 isotope')
      return
    end if

  end subroutine check_single_isotope


!> Check that the abundance of the isotope in the ENDF-6 file is unity
  subroutine check_abundance(abn_val)

    real(8) :: abn_val

    if (abn_val /= ONE) then
      call exit_status(EXIT_FAILURE,&
           'Abundance of isotope given in '//trim(filename)//' is not unity')
      return
    end if

  end subroutine check_abundance


!> Check that the treatment of average fission widths in the URR is supported
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
      call exit_status(EXIT_FAILURE, 'LFW must be 0 or 1 in '//trim(filename))
      return
    end select

  end subroutine check_fission_widths


!> Check that the upper energy is greater than the lower and
!! that the number of resonance energy ranges is allowable
  subroutine check_energy_ranges(num_ranges, e_low, e_high)

    integer :: num_ranges
    real(8) :: e_low
    real(8) :: e_high

    if (num_ranges > 2) then
      call log_message(WARNING,&
           'NER > 2: more than 2 resonance energy regions in '//trim(filename))
    end if

    if (e_high <= e_low) then
      call exit_status(EXIT_FAILURE,&
           'Upper resonance energy range bound is not greater than the lower in '//trim(filename))
      return
    end if

  end subroutine check_energy_ranges


!> Check the channel and scattering radius flags
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
        call exit_status(EXIT_FAILURE,&
             'ENDF-6 NAPS flag must be 0 or 1 when NRO is 0 in '//trim(filename))
        return
      end select
    case(1)
      call exit_status(EXIT_FAILURE,&
           'Energy-dependent scattering radius in '//trim(filename)//' not yet supported')
      return

      select case(naps_val)
      case(0)
        continue
      case(1)
        continue
      case(2)
        continue
      case default
        call exit_status(EXIT_FAILURE,&
             'ENDF-6 NAPS flag must be 0, 1, or 2 when NRO is 1 in '//trim(filename))
        return
      end select
    case default
      call exit_status(EXIT_FAILURE,&
           'ENDF-6 NRO flag must be 0 or 1 in '//trim(filename))
      return
    end select

  end subroutine check_radius_flags


!> Checks that the scattering radius is constant within a resonance energy range
  subroutine check_scattering_radius(ap_val, ap_ref)

    real(8) :: ap_val
    real(8) :: ap_ref

    if (ap_val /= ap_ref) then
      call log_message(WARNING,&
           'AP value changes within a resonance energy range in '//trim(filename))
    end if

  end subroutine check_scattering_radius


!> Checks for the expected ordering of orbital quantum numbers
  subroutine check_l_number(l_val, l_ref)

    integer :: l_val
    integer :: l_ref

    if (l_val /= l_ref) then
      call exit_status(EXIT_FAILURE,&
           'Unexpected ordering of orbital quantum numbers in '//trim(filename))
      return
    end if

  end subroutine check_l_number


!> Checks that the signs of the total angular momenta are all positive
  subroutine check_j_sign(j_val)

    real(8) :: j_val

    if (j_val < ZERO) then
      call exit_status(EXIT_FAILURE,&
           'Negative total angular momentum in '//trim(filename))
      return
    end if

  end subroutine check_j_sign


end module URR_endf6
