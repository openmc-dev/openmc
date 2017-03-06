module URR_io

  use URR_constants, only: ZERO,&
                           ONE
  use URR_error,     only: exit_status,&
                           log_message,&
                           EXIT_FAILURE,&
                           INFO
  use URR_isotope,   only: isotopes
  use URR_resonance, only: BreitWignerResonance
  use URR_settings,  only: endf_filenames,&
                           path_avg_xs,&
                           path_endf_files,&
                           path_prob_tables

  implicit none
  private
  public :: read_avg_xs,&
            read_prob_tables,&
            write_MF2

contains


!> Reads in pre-computed average URR cross sections needed for LSSF=1 treatment
  subroutine read_avg_xs(i_iso)

    integer, intent(in) :: i_iso ! isotope index

    logical :: file_exists ! does avg URR xs file exist?
    character(7)   :: readable ! is avg URR xs file readable?
    character(7)   :: zaid_str ! ZA number as a string
    character(255) :: rec      ! file record
    integer :: E_depend      ! parameter energy dependence flag
    integer :: avg_unit = 12 ! input unit
    integer :: ir            ! record index
    integer :: ZAI           ! isotope ZA number
    real(8) :: tol ! tolerance used in generating avg xs values

    if (isotopes(i_iso) % metastable) then
      write(zaid_str, '(I7)') isotopes(i_iso) % ZAI
      zaid_str = trim(adjustl(zaid_str)) // 'm'
    else
      write(zaid_str, '(I7)') isotopes(i_iso) % ZAI
    end if

    ! check that file exists and is readable
    inquire(file = trim(path_avg_xs)//trim(adjustl(zaid_str))//&
         &'-avg-urr-xs.dat', exist = file_exists, read = readable)
    if (.not. file_exists) then
      call exit_status(EXIT_FAILURE, 'Averaged URR cross sections file '&
           //trim(adjustl(zaid_str))//'-avg-urr-xs.dat does not exist.')
    else if (readable(1:3) == 'NO') then
      call exit_status(EXIT_FAILURE, 'Averaged URR cross sections file '&
           //trim(adjustl(zaid_str))//'-avg-urr-xs.dat&
           & is not readable.  Change file permissions with chmod command.')
    end if

    ! open file with average xs values
    call log_message(INFO, 'Loading average URR cross sections for ZA '//zaid_str)
    open(unit=avg_unit,&
         file=trim(path_avg_xs)//trim(adjustl(zaid_str))//'-avg-urr-xs.dat')

10  format(A255)
    ! ENDF-6 filepath
    read(avg_unit, 10) rec

    ! resonance formalism
    read(avg_unit, 10) rec

    ! number of contributing s-wave resonances
    read(avg_unit, 10) rec

    ! number of contributing p-wave resonances
    read(avg_unit, 10) rec

    ! number of contributing d-wave resonances
    read(avg_unit, 10) rec

    ! number of contributing f-wave resonances
    read(avg_unit, 10) rec

    ! model competitive resonance structure
    read(avg_unit, 10) rec

    ! parameter energy dependence
    read(avg_unit, 10) rec

    ! Faddeeva function evaluation
    read(avg_unit, 10) rec

    ! target tolerance (1sigma/mean) for averaged partial cross sections
    read(avg_unit, 10) rec

    ! number of energies
    read(avg_unit, 10) rec

    read(rec(10:80), '(i71)') isotopes(i_iso) % nE_tabs

    isotopes(i_iso) % num_avg_xs_grid = isotopes(i_iso) % nE_tabs

    ! column labels
    read(avg_unit, 10) rec

    ! allocate average xs grids
    allocate(isotopes(i_iso) % E_avg_xs(isotopes(i_iso) % num_avg_xs_grid))
    allocate(isotopes(i_iso) % avg_xs(isotopes(i_iso) % num_avg_xs_grid))

20  format(6ES24.16)
    ! read in average xs values

    do ir = 1, isotopes(i_iso) % num_avg_xs_grid
      read(avg_unit, 20)&
           isotopes(i_iso) % E_avg_xs(ir),&
           isotopes(i_iso) % avg_xs(ir) % t,&
           isotopes(i_iso) % avg_xs(ir) % n,&
           isotopes(i_iso) % avg_xs(ir) % g,&
           isotopes(i_iso) % avg_xs(ir) % f,&
           isotopes(i_iso) % avg_xs(ir) % x
    end do

    ! read max 1sigma/mean values
    read(avg_unit, *) rec

    close(avg_unit)

  end subroutine read_avg_xs


!> Reads in pre-generated URR probability tables
  subroutine read_prob_tables(i_iso)

    integer, intent(in) :: i_iso ! isotope index

    logical :: file_exists ! does probability table file exist?
    integer :: i_E           ! energy index
    integer :: i_T           ! temperature index
    integer :: i_b           ! probability band index
    integer :: tab_unit = 99 ! tables output file unit
    character(7) :: readable ! is probability table file readable?
    character(7) :: zaid_str ! ZA number as a string
    character(len=:), allocatable :: rec ! file record

    if (isotopes(i_iso) % metastable) then
      write(zaid_str, '(I7)') isotopes(i_iso) % ZAI
      zaid_str = trim(adjustl(zaid_str)) // 'm'
    else
      write(zaid_str, '(I7)') isotopes(i_iso) % ZAI
    end if

    ! check that file exists and is readable
    inquire(file = trim(path_prob_tables)//trim(adjustl(zaid_str))//&
         &'-urr-tables.dat', exist = file_exists, read = readable)
    if (.not. file_exists) then
      call exit_status(EXIT_FAILURE, 'Probability table file '//trim(adjustl(zaid_str))//&
           &'-urr-tables.dat does not exist.')
    else if (readable(1:3) == 'NO') then
      call exit_status(EXIT_FAILURE, 'Probability table file '//trim(adjustl(zaid_str))//&
           &'-urr-tables.dat is not readable.  Change file permissions with &
           &chmod command.')
    end if

    ! open probability table file
    call log_message(INFO, 'Loading probability tables for ZA '//zaid_str)
    open(unit = tab_unit, file =&
         trim(path_prob_tables)//trim(adjustl(zaid_str))//'-urr-tables.dat')

10  format(A255)
    ! ENDF-6 filepath
    read(tab_unit, 10) rec
    
    ! resonance formalism
    read(tab_unit, 10) rec

    ! number of contributing s-wave resonances
    read(tab_unit, 10) rec

    ! number of contributing p-wave resonances
    read(tab_unit, 10) rec

    ! number of contributing d-wave resonances
    read(tab_unit, 10) rec

    ! number of contributing f-wave resonances
    read(tab_unit, 10) rec

    ! model competitive resonance structure
    read(tab_unit, 10) rec

    ! parameter energy dependence
    read(tab_unit, 10) rec

    ! Faddeeva function evaluation
    read(tab_unit, 10) rec

    ! target tolerance (1sigma/mean) for averaged partial cross sections
    read(tab_unit, 10) rec

    ! number of energies
    read(tab_unit, 10) rec
    read(rec(10:80), '(i71)') isotopes(i_iso) % nE_tabs

    ! number of temperatures
    read(tab_unit, 10) rec
    read(rec(14:80), '(i67)') isotopes(i_iso) % nT_tabs

    ! number of probability-xs bands
    read(tab_unit, 10) rec
    read(rec(7:80), '(i74)') isotopes(i_iso) % n_bands

    ! allocate probability tables
    allocate(isotopes(i_iso) % E_tabs(isotopes(i_iso) % nE_tabs))
    allocate(isotopes(i_iso) % T_tabs(isotopes(i_iso) % nT_tabs))
    allocate(isotopes(i_iso) % prob_tables(isotopes(i_iso) % nE_tabs, isotopes(i_iso) % nT_tabs))
    do i_E = 1, isotopes(i_iso) % nE_tabs
      do i_T = 1, isotopes(i_iso) % nT_tabs
        allocate(isotopes(i_iso) % prob_tables(i_E, i_T) % t(isotopes(i_iso) % n_bands))
        allocate(isotopes(i_iso) % prob_tables(i_E, i_T) % n(isotopes(i_iso) % n_bands))
        allocate(isotopes(i_iso) % prob_tables(i_E, i_T) % g(isotopes(i_iso) % n_bands))
        allocate(isotopes(i_iso) % prob_tables(i_E, i_T) % f(isotopes(i_iso) % n_bands))
        allocate(isotopes(i_iso) % prob_tables(i_E, i_T) % x(isotopes(i_iso) % n_bands))
      end do
    end do

    ! read probability tables
    do i_E = 1, isotopes(i_iso) % nE_tabs
      do i_T = 1, isotopes(i_iso) % nT_tabs
        ! read energy        
        read(tab_unit, 10) rec
        read(rec(25:48), '(es24.16)') isotopes(i_iso) % E_tabs(i_E)

        ! read temperature
        read(tab_unit, 10) rec
        read(rec(25:48), '(es24.16)') isotopes(i_iso) % T_tabs(i_T)

        ! read column labels
        read(tab_unit, 10) rec

        ! read cross sections
        do i_b = 1, isotopes(i_iso) % n_bands
          read(tab_unit, 10) rec
          read(rec(49:192), '(6es24.16)')&
               isotopes(i_iso) % prob_tables(i_E, i_T) % t(i_b) % cnt_mean,&
               isotopes(i_iso) % prob_tables(i_E, i_T) % t(i_b) % xs_mean,&
               isotopes(i_iso) % prob_tables(i_E, i_T) % n(i_b) % xs_mean,&
               isotopes(i_iso) % prob_tables(i_E, i_T) % g(i_b) % xs_mean,&
               isotopes(i_iso) % prob_tables(i_E, i_T) % f(i_b) % xs_mean,&
               isotopes(i_iso) % prob_tables(i_E, i_T) % x(i_b) % xs_mean
        end do

        ! read number of batches
        read(tab_unit, 10) rec

        ! read mean cross sections
        read(tab_unit, 10) rec

        ! read 1sigma values
        read(tab_unit, 10) rec

        ! read 1sigma/mean values
        read(tab_unit, 10) rec

      end do
    end do

    close(tab_unit)

  end subroutine read_prob_tables
  

!> Writes out resonance parameters and other data for a realization of
!! unresolved resonances
  subroutine write_MF2(i_iso)

    integer, intent(in) :: i_iso ! isotope index

    integer :: res_unit = 12 ! resonance realization file unit
    integer :: i_l           ! l loop index
    integer :: i_J           ! J loop index
    integer :: i_res_l       ! resonance index for a given l, any J
    integer :: i_res_lJ      ! resonance index for a given l, J
    integer :: NRS           ! number of resonances for a given l
    integer :: num_sorted_resonances ! number of l-wave resonances that have been sorted by increasing energyg
    type(BreitWignerResonance), allocatable :: l_wave_bw_resonances(:) ! all Breit-Wigner resonances for a given l, sorted by energy
    character(7) :: zaid_str ! ZA number as a string
    character(80) :: rec     ! ENDF-6 file record

    if (isotopes(i_iso) % metastable) then
      write(zaid_str, '(I7)') isotopes(i_iso) % ZAI
      zaid_str = trim(adjustl(zaid_str)) // 'm'
    else
      write(zaid_str, '(I7)') isotopes(i_iso) % ZAI
    end if

    open(unit = res_unit, file = trim(adjustl(zaid_str))//'-urr-realization.dat')
    write(res_unit, '("ENDF-6 File:")', advance='no')
    write(res_unit, *) trim(adjustl(path_endf_files))//trim(adjustl(endf_filenames(i_iso)))

10  format(A80)

    ! HEAD
    write(rec, 10) ''
    write(rec(1:11),  '(I11)') isotopes(i_iso) % ZAI
    write(rec(12:22), '(es11.5e1)') isotopes(i_iso) % AWR
    write(rec(23:33), '(I11)') 0
    write(rec(34:44), '(I11)') 0
    write(rec(45:55), '(I11)') 1 ! NIS = 1 --> 1 isotope present in this ENDF-6 file
    write(rec(56:66), '(I11)') 0
    write(res_unit, 10) rec

    ! CONT
    write(rec, 10) ''
    write(rec(1:11),  '(i11)') isotopes(i_iso) % ZAI
    write(rec(12:22), '(es11.5e1)') ONE
    write(rec(23:33), '(I11)')   0
    write(rec(34:44), '(I11)')   0 ! LFW = 0 since URR is being replaced with a realization
    write(rec(45:55), '(I11)')   isotopes(i_iso) % NER ! same NER, URR is just being replaced by a realization
    write(rec(56:66), '(I11)')   0
    write(res_unit, 10) rec

    ! CONT
    write(rec, 10) ''
    write(rec(1:11),  '(es11.5e1)') isotopes(i_iso) % EL(isotopes(i_iso) % i_urr)
    write(rec(12:22), '(es11.5e1)') isotopes(i_iso) % EH(isotopes(i_iso) % i_urr)
    write(rec(24:33),   '(I10)') 1 ! isotopes(i_iso) % LRU(isotopes(i_iso) % i_urr): 2 --> 1 because the URR parameters are converted to a resolved realization
    write(rec(35:44),   '(I10)') 1 ! isotopes(i_iso) % LRF(isotopes(i_iso) % i_urr): new meaning of LRF in the RRR --> 1 indicates SLBW parameters
    write(rec(45:55),   '(I11)') isotopes(i_iso) % NRO(isotopes(i_iso) % i_urr)
    write(rec(56:66),   '(I11)') isotopes(i_iso) % NAPS(isotopes(i_iso) % i_urr)
    write(res_unit, 10) rec

    ! CONT (handles SLBW, LRF = 1 and MLBW, LRF = 2)
    write(rec, 10) ''
    write(rec(1:11),  '(es11.5e1)') isotopes(i_iso) % SPI(isotopes(i_iso) % i_urr)
    write(rec(12:22), '(es11.5e1)') isotopes(i_iso) % AP(isotopes(i_iso) % i_urr)
    write(rec(23:33), '(I11)')   0
    write(rec(34:44), '(I11)')   0
    write(rec(45:55), '(I11)')   isotopes(i_iso) % NLS(isotopes(i_iso) % i_urr)
    write(rec(56:66), '(I11)')   0
    write(res_unit, 10) rec

    ! LIST
    do i_l = 1, isotopes(i_iso) % NLS(isotopes(i_iso) % i_urr)
      write(rec, 10) ''
      write(rec(1:11),  '(es11.5e1)') isotopes(i_iso) % AWR
      write(rec(12:22), '(es11.5e1)') ZERO ! QX = 0.0 because no competitive width is given
      write(rec(23:33),   '(I11)') i_l - 1
      write(rec(34:44),   '(I11)') 0 ! LRX = 0 --> no competitive width given
      NRS = 0
      do i_J = 1, isotopes(i_iso) % NJS(i_l)
        NRS = NRS + size(isotopes(i_iso) % urr_resonances(i_l, 1) % J(i_J) % res(:))
      end do
      write(rec(45:55),   '(I11)') 6 * NRS
      write(rec(56:66),   '(I11)') NRS
      write(res_unit, 10) rec
      allocate(l_wave_bw_resonances(NRS))
      num_sorted_resonances = size(isotopes(i_iso) % urr_resonances(i_l, 1) % J(1) % res(:))
      l_wave_bw_resonances(1:num_sorted_resonances) =&
           isotopes(i_iso) % urr_resonances(i_l, 1) % J(1) % res(:)
      do i_J = 2, isotopes(i_iso) % NJS(i_l)
        i_res_l = 1
        do i_res_lJ = 1, size(isotopes(i_iso) % urr_resonances(i_l, 1) % J(i_J) % res(:))
          do
            if (isotopes(i_iso) % urr_resonances(i_l, 1) % J(i_J) % res(i_res_lJ) % E_lam&
                 < l_wave_bw_resonances(i_res_l) % E_lam) exit
            i_res_l = i_res_l + 1
            if (i_res_l > num_sorted_resonances) exit
          end do
          if (i_res_l == NRS) then
            l_wave_bw_resonances(i_res_l) = isotopes(i_iso) % urr_resonances(i_l, 1) % J(i_J) % res(i_res_lJ)
          else
            l_wave_bw_resonances(i_res_l+1:NRS) = l_wave_bw_resonances(i_res_l:NRS-1)
            l_wave_bw_resonances(i_res_l) = isotopes(i_iso) % urr_resonances(i_l, 1) % J(i_J) % res(i_res_lJ)
          end if
          num_sorted_resonances = num_sorted_resonances + 1
        end do
      end do
      do i_res_l = 1, NRS
        write(rec, 10) ''
        write(rec(1:66), '(6es11.5e1)')&
             l_wave_bw_resonances(i_res_l) % E_lam,&
             l_wave_bw_resonances(i_res_l) % AJ,&
             l_wave_bw_resonances(i_res_l) % GT,&
             l_wave_bw_resonances(i_res_l) % GN,&
             l_wave_bw_resonances(i_res_l) % GG,&
             l_wave_bw_resonances(i_res_l) % GF
        write(res_unit, 10) rec
      end do
      deallocate(l_wave_bw_resonances)
    end do

    ! SEND
    write(rec, 10) ''
    write(rec(1:11),  '(es11.5e1)') ZERO
    write(rec(12:22), '(es11.5e1)') ZERO
    write(rec(23:33), '(I11)')   0
    write(rec(34:44), '(I11)')   0
    write(rec(45:55), '(I11)')   0
    write(rec(56:66), '(I11)')   0
    write(res_unit, 10) rec

    close(res_unit)

  end subroutine write_MF2


end module URR_io
