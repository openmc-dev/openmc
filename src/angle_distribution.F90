module angle_distribution

  use hdf5, only: HID_T, HSIZE_T

  use algorithm, only: binary_search
  use constants, only: ZERO, ONE, HISTOGRAM, LINEAR_LINEAR
  use distribution_univariate, only: DistributionContainer, Tabular
  use hdf5_interface, only: read_attribute, get_shape, read_dataset, &
       open_dataset, close_dataset
  use random_lcg, only: prn

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
    procedure :: from_hdf5 => angle_from_hdf5
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

  subroutine angle_from_hdf5(this, group_id)
    class(AngleDistribution), intent(inout) :: this
    integer(HID_T),           intent(in)    :: group_id

    integer :: i, j
    integer :: n
    integer :: n_energy
    integer(HID_T) :: dset_id
    integer(HSIZE_T) :: dims(1), dims2(2)
    integer, allocatable :: offsets(:)
    integer, allocatable :: interp(:)
    real(8), allocatable :: temp(:,:)

    ! Get incoming energies
    dset_id = open_dataset(group_id, 'energy')
    call get_shape(dset_id, dims)
    n_energy = int(dims(1), 4)
    allocate(this % energy(n_energy))
    allocate(this % distribution(n_energy))
    call read_dataset(this % energy, dset_id)
    call close_dataset(dset_id)

    ! Get outgoing energy distribution data
    dset_id = open_dataset(group_id, 'mu')
    call read_attribute(offsets, dset_id, 'offsets')
    call read_attribute(interp, dset_id, 'interpolation')
    call get_shape(dset_id, dims2)
    allocate(temp(dims2(1), dims2(2)))
    call read_dataset(temp, dset_id)
    call close_dataset(dset_id)

    do i = 1, n_energy
      ! Determine number of outgoing energies
      j = offsets(i)
      if (i < n_energy) then
        n = offsets(i+1) - j
      else
        n = size(temp, 1) - j
      end if

      ! Create and initialize tabular distribution
      allocate(Tabular :: this % distribution(i) % obj)
      select type (mudist => this % distribution(i) % obj)
      type is (Tabular)
        mudist % interpolation = interp(i)
        allocate(mudist % x(n), mudist % p(n), mudist % c(n))
        mudist % x(:) = temp(j+1:j+n, 1)
        mudist % p(:) = temp(j+1:j+n, 2)

        ! To get answers that match ACE data, for now we still use the tabulated
        ! CDF values that were passed through to the HDF5 library. At a later
        ! time, we can remove the CDF values from the HDF5 library and
        ! reconstruct them using the PDF
        if (.true.) then
          mudist % c(:) = temp(j+1:j+n, 3)
        else
          call mudist % initialize(temp(j+1:j+n, 1), temp(j+1:j+n, 2), interp(i))
        end if
      end select

      j = j + n
    end do
  end subroutine angle_from_hdf5

end module angle_distribution
