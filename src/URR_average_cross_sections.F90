module URR_average_cross_sections

  use error,      only: fatal_error
  use global

  implicit none
  private
  public :: read_avg_urr_xs

contains

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! READ_AVG_URR_XS reads in pre-computed average URR cross sections needed for
! LSSF = 1 treatment
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_avg_urr_xs(i)

    type(Isotope), pointer :: tope ! isotope object pointer
    logical :: file_exists ! does avg URR xs file exist?
    character(7)  :: readable ! is avg URR xs file readable?
    character(6)  :: zaid_str ! ZAID number as a string
    character(80) :: rec      ! file record
    integer :: E_depend ! parameter energy dependence flag
    integer :: i        ! isotope index
    integer :: in = 12  ! input unit
    integer :: ir       ! record index
    integer :: ZAI      ! isotope ZAID number
    real(8) :: tol ! tolerance used in generating avg xs values

    tope => isotopes(i)
    write(zaid_str, '(I6)') tope % ZAI

    ! check that file exists and is readable
    inquire(file = trim(path_avg_urr_xs)//trim(adjustl(zaid_str))//&
         &'-avg-urr-xs.dat', exist = file_exists, read = readable)
    if (.not. file_exists) then
      call fatal_error('Averaged URR cross sections file '&
           //trim(adjustl(zaid_str))//'-avg-urr-xs.dat does not exist.')
    else if (readable(1:3) == 'NO') then
      call fatal_error('Averaged URR cross sections file '&
           //trim(adjustl(zaid_str))//'-avg-urr-xs.dat&
           & is not readable.  Change file permissions with chmod command.')
    end if

    ! open file with average xs values
    if (master)&
         write(*,*) 'Loading averaged URR cross sections for ZAID ='//zaid_str
    open(unit=in,&
         file=trim(path_avg_urr_xs)//trim(adjustl(zaid_str))//'-avg-urr-xs.dat')
    
10  format(A80)
    ! SHA1
    read(in, 10) rec

    ! ZAID
    read(in, 10) rec
    read(rec(1:6), '(I6)') ZAI

    ! parameter energy dependence
    read(in, 10) rec
    read(rec(1:1), '(I1)') E_depend

    ! tolerance on average partial xs max relative uncertainty (1sigma/mean)
    read(in, 10) rec
    read(rec(1:13), '(ES13.6)') tol
    if (ZAI /= tope % ZAI)&
         call fatal_error('ZAID number disagreement')
    if (master)&
        write(*, '(A55,ES13.6)') adjustl(trim(adjustl(zaid_str)))//&
             &'-avg-urr-xs.dat was generated with a tolerance of',tol

    ! read column labels
    read(in, 10) rec

    ! allocate average xs grids
    tope % nEavg = N_AVG_URR_GRID
    allocate(tope % Eavg(N_AVG_URR_GRID))
    allocate(tope % avg_urr_t(N_AVG_URR_GRID))
    allocate(tope % avg_urr_n(N_AVG_URR_GRID))
    allocate(tope % avg_urr_g(N_AVG_URR_GRID))
    allocate(tope % avg_urr_f(N_AVG_URR_GRID))
    allocate(tope % avg_urr_x(N_AVG_URR_GRID))


20  format(6ES13.6)
    ! read in average xs values
    do ir = 1, N_AVG_URR_GRID
      read(in, 20)&
           tope % Eavg(ir),&
           tope % avg_urr_t(ir),&
           tope % avg_urr_n(ir),&
           tope % avg_urr_g(ir),&
           tope % avg_urr_f(ir),&
           tope % avg_urr_x(ir)
    end do

    close(in)

  end subroutine read_avg_urr_xs

end module URR_average_cross_sections
