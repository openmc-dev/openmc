module URR_io

  use URR_error,     only: exit_status,&
                           log_message,&
                           EXIT_FAILURE,&
                           INFO
  use URR_isotope,   only: isotopes
  use URR_settings,  only: path_avg_xs,&
                           path_prob_tables

  implicit none
  private
  public :: read_avg_xs,&
       read_prob_tables

contains


!> Reads in pre-computed average URR cross sections needed for LSSF=1 treatment
  subroutine read_avg_xs(i)

    logical :: file_exists ! does avg URR xs file exist?
    character(7) :: readable ! is avg URR xs file readable?
    character(6) :: zaid_str ! ZAID number as a string
    character(255) :: rec    ! file record
    integer :: E_depend    ! parameter energy dependence flag
    integer :: i           ! isotope index
    integer :: i_unit = 12 ! input unit
    integer :: ir          ! record index
    integer :: ZAI         ! isotope ZAID number
    real(8) :: tol ! tolerance used in generating avg xs values

    write(zaid_str, '(I6)') isotopes(i) % ZAI

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
    call log_message(INFO, 'Loading average URR cross sections for ZAI '//zaid_str)
    open(unit=i_unit,&
         file=trim(path_avg_xs)//trim(adjustl(zaid_str))//'-avg-urr-xs.dat')

10  format(A255)
    ! ENDF-6 filepath
    read(i_unit, 10) rec
    
    ! resonance formalism
    read(i_unit, 10) rec

    ! number of contributing s-wave resonances
    read(i_unit, 10) rec

    ! number of contributing p-wave resonances
    read(i_unit, 10) rec

    ! number of contributing d-wave resonances
    read(i_unit, 10) rec

    ! number of contributing f-wave resonances
    read(i_unit, 10) rec

    ! model competitive resonance structure
    read(i_unit, 10) rec

    ! parameter energy dependence
    read(i_unit, 10) rec

    ! Faddeeva function evaluation
    read(i_unit, 10) rec

    ! target tolerance (1sigma/mean) for averaged partial cross sections
    read(i_unit, 10) rec

    ! number of energies
    read(i_unit, 10) rec
    read(rec(10:255), '(i246)') isotopes(i) % nE_tabs
    isotopes(i) % num_avg_xs_grid = isotopes(i) % nE_tabs

    ! column labels
    read(i_unit, 10) rec

    ! allocate average xs grids
    allocate(isotopes(i) % E_avg_xs(isotopes(i) % num_avg_xs_grid))
    allocate(isotopes(i) % avg_xs(isotopes(i) % num_avg_xs_grid))

20  format(6ES24.16)
    ! read in average xs values
    do ir = 1, isotopes(i) % num_avg_xs_grid
      read(i_unit, 20)&
           isotopes(i) % E_avg_xs(ir),&
           isotopes(i) % avg_xs(ir) % t,&
           isotopes(i) % avg_xs(ir) % n,&
           isotopes(i) % avg_xs(ir) % g,&
           isotopes(i) % avg_xs(ir) % f,&
           isotopes(i) % avg_xs(ir) % x
    end do

    ! read max 1sigma/mean values
    read(i_unit, *) rec

    close(i_unit)

  end subroutine read_avg_xs


!> Reads in pre-generated URR probability tables
  subroutine read_prob_tables(i)

    logical :: file_exists ! does probability table file exist?
    integer :: i             ! isotope index
    integer :: i_E           ! energy index
    integer :: i_T           ! temperature index
    integer :: i_b           ! probability band index
    integer :: tab_unit = 99 ! tables output file unit
    character(7) :: readable ! is probability table file readable?
    character(6) :: zaid_str ! ZAID number as a string
    character(:), allocatable :: rec ! file record

    write(zaid_str, '(i6)') isotopes(i) % ZAI

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
    call log_message(INFO, 'Loading probability tables for ZAID ='//zaid_str)
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
    read(rec(10:80), '(i71)') isotopes(i) % nE_tabs

    ! number of temperatures
    read(tab_unit, 10) rec
    read(rec(14:80), '(i67)') isotopes(i) % nT_tabs

    ! number of probability-xs bands
    read(tab_unit, 10) rec
    read(rec(7:80), '(i74)') isotopes(i) % n_bands

    ! allocate probability tables
    allocate(isotopes(i) % E_tabs(isotopes(i) % nE_tabs))
    allocate(isotopes(i) % T_tabs(isotopes(i) % nT_tabs))
    allocate(isotopes(i) % prob_tables(isotopes(i) % nE_tabs, isotopes(i) % nT_tabs))
    do i_E = 1, isotopes(i) % nE_tabs
      do i_T = 1, isotopes(i) % nT_tabs
        allocate(isotopes(i) % prob_tables(i_E, i_T) % t(isotopes(i) % n_bands))
        allocate(isotopes(i) % prob_tables(i_E, i_T) % n(isotopes(i) % n_bands))
        allocate(isotopes(i) % prob_tables(i_E, i_T) % g(isotopes(i) % n_bands))
        allocate(isotopes(i) % prob_tables(i_E, i_T) % f(isotopes(i) % n_bands))
        allocate(isotopes(i) % prob_tables(i_E, i_T) % x(isotopes(i) % n_bands))
      end do
    end do

    ! read probability tables
    do i_E = 1, isotopes(i) % nE_tabs
      do i_T = 1, isotopes(i) % nT_tabs
        ! read energy        
        read(tab_unit, 10) rec
        read(rec(25:48), '(es24.16)') isotopes(i) % E_tabs(i_E)

        ! read temperature
        read(tab_unit, 10) rec
        read(rec(25:48), '(es24.16)') isotopes(i) % T_tabs(i_T)

        ! read column labels
        read(tab_unit, 10) rec

        ! read cross sections
        do i_b = 1, isotopes(i) % n_bands
          read(tab_unit, 10) rec
          read(rec(49:192), '(6es24.16)')&
               isotopes(i) % prob_tables(i_E, i_T) % t(i_b) % cnt_mean,&
               isotopes(i) % prob_tables(i_E, i_T) % t(i_b) % xs_mean,&
               isotopes(i) % prob_tables(i_E, i_T) % n(i_b) % xs_mean,&
               isotopes(i) % prob_tables(i_E, i_T) % g(i_b) % xs_mean,&
               isotopes(i) % prob_tables(i_E, i_T) % f(i_b) % xs_mean,&
               isotopes(i) % prob_tables(i_E, i_T) % x(i_b) % xs_mean
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

  end subroutine read_prob_tables
  

end module URR_io
