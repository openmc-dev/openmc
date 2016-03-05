module secondary_header

  use constants, only: ZERO
  use endf_header, only: Tab1
  use interpolation, only: interpolate_tab1
  use random_lcg, only: prn

!===============================================================================
! ANGLEENERGY (abstract) defines a correlated or uncorrelated angle-energy
! distribution that is a function of incoming energy. Each derived type must
! implement a sample() subroutine that returns an outgoing energy and scattering
! cosine given an incoming energy.
!===============================================================================

  type, abstract :: AngleEnergy
  contains
    procedure(angleenergy_sample_), deferred :: sample
  end type AngleEnergy

  abstract interface
    subroutine angleenergy_sample_(this, E_in, E_out, mu)
      import AngleEnergy
      class(AngleEnergy), intent(in) :: this
      real(8), intent(in) :: E_in
      real(8), intent(out) :: E_out
      real(8), intent(out) :: mu
    end subroutine angleenergy_sample_
  end interface

  type :: AngleEnergyContainer
    class(AngleEnergy), allocatable :: obj
  end type AngleEnergyContainer

!===============================================================================
! SECONDARYDISTRIBUTION stores multiple angle-energy distributions, each of
! which has a given probability of occurring for a given incoming energy. In
! general, most secondary distributions only have one angle-energy distribution,
! but for some cases (e.g., (n,2n) in certain nuclides) multiple distinct
! distributions exist.
!===============================================================================

  type :: SecondaryDistribution
    type(Tab1), allocatable :: applicability(:)
    type(AngleEnergyContainer), allocatable :: distribution(:)
  contains
    procedure :: sample => secondary_sample
  end type SecondaryDistribution

contains

  subroutine secondary_sample(this, E_in, E_out, mu)
    class(SecondaryDistribution), intent(in) :: this
    real(8), intent(in)  :: E_in  ! incoming energy
    real(8), intent(out) :: E_out ! sampled outgoing energy
    real(8), intent(out) :: mu    ! sampled scattering cosine

    integer :: i       ! loop counter
    integer :: n       ! number of angle-energy distributions
    real(8) :: prob    ! cumulative probability
    real(8) :: c       ! sampled cumulative probability

    n = size(this%applicability)
    if (n > 1) then
      prob = ZERO
      c = prn()
      do i = 1, n
        ! Determine probability that i-th energy distribution is sampled
        prob = prob + interpolate_tab1(this%applicability(i), E_in)

        ! If i-th distribution is sampled, sample energy from the distribution
        if (c <= prob) then
          call this%distribution(i)%obj%sample(E_in, E_out, mu)
          exit
        end if
      end do
    else
      ! If only one distribution is present, go ahead and sample it
      call this%distribution(1)%obj%sample(E_in, E_out, mu)
    end if

  end subroutine secondary_sample

end module secondary_header
