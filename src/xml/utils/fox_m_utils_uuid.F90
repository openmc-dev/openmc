module fox_m_utils_uuid

  !This generates UUIDs according to RFC 4122

  ! Only types 1 (time-based) and 4 (pseudo-RNG-based) are implemented.

#ifndef DUMMYLIB
  use fox_m_utils_mtprng, only : mtprng_state, mtprng_init, mtprng_rand64
#endif

  implicit none
  private
  
#ifndef DUMMYLIB
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: i8b = selected_int_kind(18)
  
  character, parameter :: hexdigits(0:15) = &
    (/'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'/)
  
  type(mtprng_state), save :: rng_state
  logical, save :: initialized = .false.
  integer, save :: values_save ! must be default for date_and_time
  integer(kind=i4b), save :: hires_count = 0

! clock-seq holds a random number constant for the lifetime of the program
! using this module. That's the best we can do per S 4.1.5
  integer, save :: clock_seq = 0

#endif

  public :: generate_uuid
  
contains
  
  function generate_uuid(version) result(uuid)
    integer, intent(in), optional :: version
    character(len=36) :: uuid

#ifdef DUMMYLIB
    uuid = ""
#else
    integer(kind=i8b) :: timestamp, node
    integer(kind=i4b) :: clock_sequence

    integer(kind=i4b) :: time_low, time_mid, time_hi_and_version
    integer(kind=i4b) :: clk_seq_hi_res, clk_seq_low

    integer :: values(8) ! must be default for date_and_time
    integer(kind=i4b) :: variant, v


    if (.not.initialized) then
      ! Use the current date and time to init mtprng
      ! but this gives limited varaibility, so mix 
      ! the result up.  Can we do better? In any
      ! case, this gets passed through a quick 
      ! generator inside mtprng_init.
      call date_and_time(values=values)
      values(7) = values(7)*1000+values(5)*100+values(3)*10+values(1)
      values(8) = values(2)*1000+values(4)*100+values(6)*10+values(8)
      call mtprng_init(int(values(7)*10000+values(8), i4b), rng_state)
      clock_seq = int(mtprng_rand64(rng_state), i4b)
      initialized = .true.
    endif

    variant = 1

    if (present(version)) then
      v = version
    else
      v = 4
    endif

    select case (v)
    case (0)
      ! Nil UUID  - S 4.1.7
      uuid = repeat('0',8)//'-'//repeat('0',4)//'-'//repeat('0',4)// &
        '-'//repeat('0',4)//'-'//repeat('0',12)
      return
    case(1)
      call date_and_time(values=values)
      ! In case of too-frequent requests, we will replace time_low
      ! with the count below ...
      if (all(values==values_save)) then
        hires_count = hires_count + 1
      else
        hires_count = 0
      endif
    case(2-3)
      !Unimplemented
      uuid = ''
      return
    case(4)
      continue
    case(5)
      !Unimplemented
      uuid = ''
      return
    case default
      !Unspecified
      uuid = ''
      return
    end select

!4.1.4 Timestamp

    select case(v)
    case(1)
      timestamp = get_utc_since_1582(values)
    case(4)
      timestamp = ior(mtprng_rand64(rng_state), ishft(mtprng_rand64(rng_state), 28))
    end select

!4.1.5 Clock Sequence
    ! 14 bits
    select case(v)
    case(1)
      clock_sequence = clock_seq
    case(4)
      clock_sequence = int(mtprng_rand64(rng_state), i4b)
    end select

!4.1.6 Node
    ! 48 bits
    select case(v)
    case(1)
      node = ior(mtprng_rand64(rng_state), ishft(mtprng_rand64(rng_state), 16))
      ! No MAC address accessible - see section 4.5 !FIXME
    case(4)
      node = ior(mtprng_rand64(rng_state), ishft(mtprng_rand64(rng_state), 16))
    end select

    time_low = ibits(timestamp, 0, 32)
    time_mid = ibits(timestamp, 32, 16)
    if (hires_count==0) then
      time_hi_and_version = ior(int(ibits(timestamp, 48, 12), i4b), ishft(v, 12))
    else
      time_hi_and_version = ior(hires_count, ishft(v, 12))
    endif

    clk_seq_low = ibits(clock_sequence, 0, 8)
    clk_seq_hi_res = ior(ibits(clock_sequence, 8, 6), ishft(variant, 6))

    uuid = int32ToHexOctets(time_low, 4)//"-"// &
      int32ToHexOctets(time_mid, 2)//"-"// &
      int32ToHexOctets(time_hi_and_version, 2)//"-"// & 
      int32ToHexOctets(clk_seq_hi_res, 1)// &
      int32ToHexOctets(clk_seq_low, 1)//"-"// &
      int64ToHexOctets(node, 6)

  contains

    function int32ToHexOctets(b, n) result(s)
      integer(i4b), intent(in) :: b
      integer, intent(in) :: n ! number of octets to print
      character(len=2*n) :: s
      
      integer :: i
      
      do i = 0, 2*n-1
        s(2*n-i:2*n-i) = hexdigits(ibits(b, i*4, 4))
      enddo
      
    end function int32ToHexOctets
    function int64ToHexOctets(b, n) result(s)
      integer(i8b), intent(in) :: b
      integer, intent(in) :: n ! number of octets to print
      character(len=2*n) :: s
      
      integer :: i
      
      do i = 0, 2*n-1
        s(2*n-i:2*n-i) = hexdigits(ibits(b, i*4, 4))
      enddo
      
    end function int64ToHexOctets

#endif 
  end function generate_uuid

#ifndef DUMMYLIB
  function get_utc_since_1582(values) result(ns)
    ! This subroutine is a little broken. It only works
    ! for times after 1/1/2006 and takes no account
    ! of any future leapseconds. It ought to serve regardless.

    ! It returns the number of 100-ns intervals since 1582-10-15-00-00-00

    integer, dimension(8), intent(in) :: values
    integer(kind=i8b) :: ns

    integer :: days
    integer :: years

    integer, parameter :: days_in_normal_year(12) = &
      (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    ns = 23_i8b * 1000_i8b * 1000_i8b * 10_i8b ! 23 leap seconds until 24:00:00 31/12/2005

    ! A count of the 100-nanosecond intervals since the
    ! beginning of the day.
    ns = ns &
      ! milliseconds
      + int(values(8), i8b)             * 10_i8b * 1000_i8b &
      ! seconds
      + int(values(7), i8b)             * 10_i8b * 1000_i8b * 1000_i8b &
      ! minutes (with timezone adjustment)
      + int(values(6) + values(4), i8b) * 10_i8b * 1000_i8b * 1000_i8b * 60_i8b &
      ! hours
      + int(values(5), i8b)             * 10_i8b * 1000_i8b * 1000_i8b * 60_i8b * 60_i8b

    ! Number of days this year:
    days = sum(days_in_normal_year(:values(2)-1))
    days = days + values(3) - 1 !add days in current month
    if (values(2)>2 .and. isLeapYear(values(1))) then
      days = days + 1
    endif
    !That's all the time since the turn of this year

    days = days + 78 ! From the start of 15th Oct to the end of 31st Dec in 1582
    !That's the additional time before the turn of the year 1583

    days = days + 102  ! 102 leap years from 1584 to 2000 inclusive
    ! That's all the intercalataed days until 2000
    
    years = values(1) - 2000 - 1 ! years since 2000 - not including this year

    days = days + years/4 - years/100 + years/400 !Add extra leap days to this total:
    ! That's all out intercalated days - remaining years are all 365 days long.

    years = years + 418 ! Add the years from 1583-2000 inclusive back on.

    ! Multiply by number of time units in one day & add to today's total.
    ns = ns + 864000000000_i8b * (int(days,i8b) + 365_i8b * int(years,i8b))

  contains
    function isLeapYear(y) result(p)
      integer, intent(in) :: y
      logical :: p
      p = (mod(y,4)==0 .and. .not.mod(y,100)==0 .or. mod(y,400)==0)
    end function isLeapYear

  end function get_utc_since_1582

#endif
end module fox_m_utils_uuid
