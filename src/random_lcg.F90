module random_lcg

  implicit none

  private
  save

  integer(8) :: prn_seed0  ! original seed
  integer(8) :: prn_seed   ! current seed
  integer(8) :: prn_mult   ! multiplication factor, g
  integer(8) :: prn_add    ! additive factor, c
  integer    :: prn_bits   ! number of bits, M
  integer(8) :: prn_mod    ! 2^M
  integer(8) :: prn_mask   ! 2^M - 1
  integer(8) :: prn_stride ! stride between particles
  real(8)    :: prn_norm   ! 2^(-M)

  public :: prn
  public :: initialize_prng
  public :: set_particle_seed
  public :: prn_skip

contains

!===============================================================================
! PRN generates a pseudo-random number using a linear congruential generator
!===============================================================================

  function prn() result(pseudo_rn)

    real(8) :: pseudo_rn

    ! This algorithm uses bit-masking to find the next integer(8) value to be
    ! used to calculate the random number

    prn_seed = iand(prn_mult*prn_seed + prn_add, prn_mask)

    ! Once the integer is calculated, we just need to divide by 2**m,
    ! represented here as multiplying by a pre-calculated factor

    pseudo_rn = prn_seed * prn_norm

  end function prn

!===============================================================================
! INITIALIZE_PRNG sets up the random number generator, determining the seed and
! values for g, c, and m.
!===============================================================================

  subroutine initialize_prng()

    use global, only: seed

    prn_seed0  = seed
    prn_seed   = prn_seed
    prn_mult   = 2806196910506780709_8
    prn_add    = 1_8
    prn_bits   = 63
    prn_mod    = ibset(0_8, prn_bits)  ! clever way of calculating 2**bits
    prn_mask   = prn_mod - 1_8
    prn_stride = 152917_8
    prn_norm   = 2._8**(-prn_bits)

  end subroutine initialize_prng

!===============================================================================
! SET_PARTICLE_SEED sets the seed to a unique value based on the ID of the
! particle
!===============================================================================

  subroutine set_particle_seed(id)

    integer(8), intent(in) :: id

    prn_seed = prn_skip_ahead(id*prn_stride, prn_seed0)

  end subroutine set_particle_seed
    
!===============================================================================
! PRN_SKIP advances the random number seed 'n' times from the current seed
!===============================================================================

  subroutine prn_skip(n)

    integer(8), intent(in) :: n ! number of seeds to skip

    prn_seed = prn_skip_ahead(n, prn_seed)

  end subroutine prn_skip

!===============================================================================
! PRN_SKIP_AHEAD advances the random number seed 'skip' times. This is usually
! used to skip a fixed number of random numbers (the stride) so that a given
! particle always has the same starting seed regardless of how many processors
! are used
!===============================================================================

  function prn_skip_ahead(n, seed) result(new_seed)

    integer(8), intent(in) :: n        ! number of seeds to skip
    integer(8), intent(in) :: seed     ! original seed
    integer(8)             :: new_seed ! new seed

    integer(8) :: nskip
    integer(8) :: gen
    integer(8) :: g
    integer(8) :: inc
    integer(8) :: c
    integer(8) :: gp

    ! In cases where we want to skip backwards, we add the period of the random
    ! number generator until the number of PRNs to skip is positive since
    ! skipping ahead that much is the same as skipping backwards by the original
    ! amount

    nskip = n
    do while (nskip < 0_8)
       nskip = nskip + prn_mod
    enddo

    ! The algorithm here to determine the parameters used to skip ahead is
    ! described in F. Brown, "Random Number Generation with Arbitrary Stride,"
    ! Trans. Am. Nucl. Soc. (Nov. 1994). This algorithm is able to skip ahead in
    ! O(log2(N)) operations instead of O(N). Basically, it computes parameters G
    ! and C which can then be used to find x_N = G*x_0 + C mod 2^M.

    nskip = iand(nskip, prn_mask)
    gen   = 1
    g     = prn_mult
    inc   = 0
    c     = prn_add
    do while (nskip > 0_8)
       if (btest(nskip,0)) then
          gen = iand(gen*g, prn_mask)
          inc = iand(inc*g, prn_mask)
          inc = iand(inc+c, prn_mask)
       endif
       gp    = iand(g+1,  prn_mask)
       g     = iand(g*g,  prn_mask)
       c     = iand(gp*c, prn_mask)
       nskip = ishft(nskip, -1)
    enddo

    ! With G and C, we can now find the new seed
    new_seed = iand(gen*seed + inc, prn_mask)

  end function prn_skip_ahead

end module random_lcg
