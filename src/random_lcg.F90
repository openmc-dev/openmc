module random_lcg

  use, intrinsic :: ISO_C_BINDING

  use constants

  implicit none

  private
  save

  ! Starting seed
  integer(C_INT64_T), public, bind(C) :: seed = 1_8

  ! LCG parameters
  integer(C_INT64_T), parameter :: prn_mult = 2806196910506780709_8 ! multiplication factor, g
  integer(C_INT64_T), parameter :: prn_add = 1_8                    ! additive factor, c
  integer,            parameter :: prn_bits = 63                    ! number of bits, M
  integer(C_INT64_T), parameter :: prn_mod = ibset(0_8, prn_bits)   ! 2^M
  integer(C_INT64_T), parameter :: prn_mask = not(prn_mod)          ! 2^M - 1
  integer(C_INT64_T), parameter :: prn_stride = 152917_8            ! stride between particles
  real(C_DOUBLE),     parameter :: prn_norm = 2._8**(-prn_bits)     ! 2^(-M)

  ! Current PRNG state
  integer(C_INT64_T) :: prn_seed(N_STREAMS) ! current seed
  integer            :: stream     ! current RNG stream
!$omp threadprivate(prn_seed, stream)

  public :: prn
  public :: future_prn
  public :: set_particle_seed
  public :: advance_prn_seed
  public :: prn_set_stream
  public :: openmc_set_seed

contains

!===============================================================================
! PRN generates a pseudo-random number using a linear congruential generator
!===============================================================================

  function prn() result(pseudo_rn)

    real(C_DOUBLE) :: pseudo_rn

    ! This algorithm uses bit-masking to find the next integer(C_INT64_T) value
    ! to be used to calculate the random number

    prn_seed(stream) = iand(prn_mult*prn_seed(stream) + prn_add, prn_mask)

    ! Once the integer is calculated, we just need to divide by 2**m,
    ! represented here as multiplying by a pre-calculated factor

    pseudo_rn = prn_seed(stream) * prn_norm

  end function prn

!===============================================================================
! FUTURE_PRN generates a pseudo-random number which is 'n' times ahead from the
! current seed.
!===============================================================================

  function future_prn(n) result(pseudo_rn)

    integer(C_INT64_T), intent(in) :: n        ! number of prns to skip

    real(C_DOUBLE) :: pseudo_rn

    pseudo_rn  = future_seed(n, prn_seed(stream)) * prn_norm

  end function future_prn

!===============================================================================
! SET_PARTICLE_SEED sets the seed to a unique value based on the ID of the
! particle
!===============================================================================

  subroutine set_particle_seed(id)

    integer(C_INT64_T), intent(in) :: id

    integer :: i

    do i = 1, N_STREAMS
      prn_seed(i) = future_seed(id*prn_stride, seed + i - 1)
    end do

  end subroutine set_particle_seed

!===============================================================================
! ADVANCE_PRN_SEED advances the random number seed 'n' times from the current
! seed.
!===============================================================================

  subroutine advance_prn_seed(n)

    integer(C_INT64_T), intent(in) :: n ! number of seeds to skip

    prn_seed(stream) = future_seed(n, prn_seed(stream))

  end subroutine advance_prn_seed

!===============================================================================
! FUTURE_SEED advances the random number seed 'skip' times. This is usually
! used to skip a fixed number of random numbers (the stride) so that a given
! particle always has the same starting seed regardless of how many processors
! are used
!===============================================================================

  function future_seed(n, seed) result(new_seed)

    integer(C_INT64_T), intent(in) :: n        ! number of seeds to skip
    integer(C_INT64_T), intent(in) :: seed     ! original seed
    integer(C_INT64_T)             :: new_seed ! new seed

    integer(C_INT64_T) :: nskip ! positive number of seeds to skip
    integer(C_INT64_T) :: g     ! original multiplicative constant
    integer(C_INT64_T) :: c     ! original additive constnat
    integer(C_INT64_T) :: g_new ! new effective multiplicative constant
    integer(C_INT64_T) :: c_new ! new effective additive constant

    ! In cases where we want to skip backwards, we add the period of the random
    ! number generator until the number of PRNs to skip is positive since
    ! skipping ahead that much is the same as skipping backwards by the original
    ! amount

    nskip = n
    do while (nskip < 0_8)
      nskip = nskip + prn_mod
    end do

    ! Make sure nskip is less than 2^M
    nskip = iand(nskip, prn_mask)

    ! The algorithm here to determine the parameters used to skip ahead is
    ! described in F. Brown, "Random Number Generation with Arbitrary Stride,"
    ! Trans. Am. Nucl. Soc. (Nov. 1994). This algorithm is able to skip ahead in
    ! O(log2(N)) operations instead of O(N). Basically, it computes parameters G
    ! and C which can then be used to find x_N = G*x_0 + C mod 2^M.

    ! Initialize constants
    g     = prn_mult
    c     = prn_add
    g_new = 1
    c_new = 0
    BIT_LOOP: do while (nskip > 0_8)
      ! Check if least significant bit is 1
      if (btest(nskip,0)) then
        g_new = iand(g_new*g, prn_mask)
        c_new = iand(c_new*g + c, prn_mask)
      endif
      c = iand((g+1)*c, prn_mask)
      g = iand(g*g, prn_mask)

      ! Move bits right, dropping least significant bit
      nskip = ishft(nskip, -1)
    end do BIT_LOOP

    ! With G and C, we can now find the new seed
    new_seed = iand(g_new*seed + c_new, prn_mask)

  end function future_seed

!===============================================================================
! PRN_SET_STREAM changes the random number stream. If random numbers are needed
! in routines not used directly for tracking (e.g. physics), this allows the
! numbers to be generated without affecting reproducibility of the physics.
!===============================================================================

  subroutine prn_set_stream(i)

    integer, intent(in) :: i

    stream = i

  end subroutine prn_set_stream

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_set_seed(new_seed) result(err) bind(C)
    ! Saves the starting seed and sets up the PRNG thread state
    integer(C_INT64_T), value, intent(in) :: new_seed
    integer(C_INT) :: err

    integer :: i

    err = 0
    seed = new_seed
!$omp parallel
    do i = 1, N_STREAMS
      prn_seed(i) = seed + i - 1
    end do
    stream = STREAM_TRACKING
!$omp end parallel

  end function openmc_set_seed

end module random_lcg
