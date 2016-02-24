module product_header

  use angleenergy_header, only: AngleEnergyContainer
  use constants, only: ZERO
  use endf_header, only: Tab1
  use interpolation, only: interpolate_tab1
  use random_lcg, only: prn

!===============================================================================
! REACTIONPRODUCT stores a data for a reaction product including its yield and
! angle-energy distributions, each of which has a given probability of occurring
! for a given incoming energy. In general, most products only have one
! angle-energy distribution, but for some cases (e.g., (n,2n) in certain
! nuclides) multiple distinct distributions exist.
!===============================================================================

  type :: ReactionProduct
    integer :: yield                    ! Number of secondary particles released
    logical :: yield_with_E = .false.   ! Flag to indicate E-dependent yield
    type(Tab1), pointer :: yield_E => null()  ! Energy-dependent neutron yield
    type(Tab1), allocatable :: applicability(:)
    type(AngleEnergyContainer), allocatable :: distribution(:)
  contains
    procedure :: sample => reactionproduct_sample
  end type ReactionProduct

contains

  subroutine reactionproduct_sample(this, E_in, E_out, mu)
    class(ReactionProduct), intent(in) :: this
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

  end subroutine reactionproduct_sample

end module product_header
