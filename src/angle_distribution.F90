module angle_distribution

  use constants, only: ZERO, ONE
  use distribution_univariate, only: DistributionContainer
  use random_lcg, only: prn
  use search, only: binary_search

  implicit none
  private

!===============================================================================
! ANGLEDISTRIBUTION represents an angular distribution that is to be used in an
! uncorrelated angle-energy distribution. This occurs whenever the angle
! distrbution is given in File 4 in an ENDF file. The distribution of angles
! depends on the incoming energy of the neutron, so this type stores a
! distribution for each of a set of incoming energies.
!===============================================================================

  type, public :: AngleDistribution
    real(8), allocatable :: energy(:)
    type(DistributionContainer), allocatable :: distribution(:)
  contains
    procedure :: sample => angle_sample
  end type AngleDistribution

contains

  function angle_sample(this, E) result(mu)
    class(AngleDistribution), intent(in) :: this
    real(8), intent(in) :: E  ! incoming energy
    real(8)             :: mu ! sampled cosine of scattering angle

    integer :: i  ! index on incoming energy grid
    integer :: n  ! number of incoming energies
    real(8) :: r  ! interpolation factor on incoming energy grid

    ! Determine number of incoming energies
    n = size(this%energy)

    ! Find energy bin and calculate interpolation factor -- if the energy is
    ! outside the range of the tabulated energies, choose the first or last bins
    if (E < this%energy(1)) then
      i = 1
      r = ZERO
    elseif (E > this%energy(n)) then
      i = n - 1
      r = ONE
    else
      i = binary_search(this%energy, n, E)
      r = (E - this%energy(i))/(this%energy(i+1) - this%energy(i))
    end if

    ! Sample between the ith and (i+1)th bin
    if (r > prn()) i = i + 1

    ! Sample i-th distribution
    mu = this%distribution(i)%obj%sample()

    ! Make sure mu is in range [-1,1]
    if (abs(mu) > ONE) mu = sign(ONE, mu)
  end function angle_sample

end module angle_distribution
