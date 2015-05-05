module avg_urr_xs_values

  use error,      only: fatal_error
  use global
  use unresolved, only: Isotope, &
                        isotopes, &
                        path_avg_urr_xs, &
                        represent_params

  implicit none

contains

  subroutine read_avg_urr_xs(i)

    type(Isotope), pointer :: tope => null() ! isotope object pointer
    logical :: file_exists ! does avg URR xs file exist?
    character(7)  :: readable ! is avg URR xs file readable?
    character(6)  :: zaid_str ! ENDF-6 MAT number as a string
    character(80) :: rec      ! file record
    integer :: E_depend ! parameter energy dependence flag
    integer :: i        ! isotope index
    integer :: in = 12  ! input unit
    integer :: ir       ! record index
    integer :: ZAI      ! isotope ZAID number
    real(8) :: tol ! tolerance used in generating avg xs values
!$omp threadprivate(tope)

    tope => isotopes(i)

    write(zaid_str, '(I6)') tope % ZAI
    inquire(file = trim(path_avg_urr_xs)&
      //trim(adjustl(zaid_str))//'-avg-urr-xs.dat',&
      exist = file_exists, read = readable)
    if (.not. file_exists) then
      call fatal_error('Averaged URR cross sections file '&
        //trim(adjustl(zaid_str))//'-avg-urr-xs.dat'//' does not exist.')
    else if (readable(1:3) == 'NO') then
      call fatal_error('Averaged URR cross sections file '&
        //trim(adjustl(zaid_str))//'-avg-urr-xs.dat'&
        //' is not readable.  Change file permissions with chmod command.')
    end if

    open(unit = in, file = &
      trim(path_avg_urr_xs)//trim(adjustl(zaid_str))//'-avg-urr-xs.dat')
    read(in, 10) rec
10  format(A80)
    read(rec(1:6), '(I6)') ZAI
    read(in, 10) rec
    read(rec(1:1), '(I1)') E_depend
    read(in, 10) rec
    read(rec(1:13), '(ES13.6)') tol
    read(in, 10) rec
    if (ZAI /= tope % ZAI)&
      call fatal_error('ZAID number disagreement')
    if (master)&
      write(*, '(A56,ES13.6)') trim(adjustl(zaid_str))//'-avg-urr-xs.dat'&
      //' was generated with a tolerance of ',tol
20  format(ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6)
    allocate(tope % Eavg(50))
    allocate(tope % avg_urr_t(50))
    allocate(tope % avg_urr_n(50))
    allocate(tope % avg_urr_g(50))
    allocate(tope % avg_urr_f(50))
    allocate(tope % avg_urr_x(50))
    do ir = 1, 50
      read(in, 20)&
        tope % Eavg(ir),&
        tope % avg_urr_t(ir),&
        tope % avg_urr_n(ir),&
        tope % avg_urr_g(ir),&
        tope % avg_urr_f(ir),&
        tope % avg_urr_x(ir)
    end do

    tope % nEavg = size(tope % Eavg)

    close(in)

  end subroutine read_avg_urr_xs

end module avg_urr_xs_values
