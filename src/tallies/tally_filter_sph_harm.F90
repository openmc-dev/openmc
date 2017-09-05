module tally_filter_sph_harm

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T

  use constants
  use error
  use hdf5_interface
  use math,                only: calc_pn, calc_rn
  use particle_header,     only: Particle
  use string,              only: to_str, to_lower
  use tally_filter_header
  use xml_interface

  implicit none
  private

  integer, parameter :: COSINE_SCATTER = 1
  integer, parameter :: COSINE_PARTICLE = 2

!===============================================================================
! LEGENDREFILTER gives Legendre moments of the change in scattering angle
!===============================================================================

  type, public, extends(TallyFilter) :: SphericalHarmonicsFilter
    integer :: order
    integer :: cosine
  contains
    procedure :: from_xml
    procedure :: get_all_bins
    procedure :: to_statepoint
    procedure :: text_label
  end type SphericalHarmonicsFilter

contains

!===============================================================================
! SphericalHarmonicsFilter methods
!===============================================================================

  subroutine from_xml(this, node)
    class(SphericalHarmonicsFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    character(MAX_WORD_LEN) :: temp_str

    ! Get specified order
    call get_node_value(node, "order", this % order)
    this % n_bins = (this % order + 1)**2

    ! Determine how cosine term is to be treated
    if (check_for_node(node, "cosine")) then
      call get_node_value(node, "cosine", temp_str)
      select case (to_lower(temp_str))
      case ('scatter')
        this % cosine = COSINE_SCATTER
      case ('particle')
        this % cosine = COSINE_PARTICLE
      end select
    else
      this % cosine = COSINE_PARTICLE
    end if
  end subroutine from_xml

  subroutine get_all_bins(this, p, estimator, match)
    class(SphericalHarmonicsFilter), intent(in)  :: this
    type(Particle),      intent(in)  :: p
    integer,             intent(in)  :: estimator
    type(TallyFilterMatch),   intent(inout) :: match

    integer :: i, j, n
    integer :: num_nm
    real(8) :: wgt
    real(8) :: rn(2*this % order + 1)

    ! TODO: Use recursive formula to calculate higher orders
    j = 0
    do n = 0, this % order
      ! Determine cosine term for scatter expansion if necessary
      if (this % cosine == COSINE_SCATTER) then
        wgt = calc_pn(n, p % mu)
      else
        wgt = ONE
      end if

      ! Calculate n-th order spherical harmonics for (u,v,w)
      num_nm = 2*n + 1
      rn(1:num_nm) = calc_rn(n, p % last_uvw)

      ! Append matching (bin,weight) for each moment
      do i = 1, num_nm
        j = j + 1
        call match % bins % push_back(j)
        call match % weights % push_back(wgt * rn(i))
      end do
    end do
  end subroutine get_all_bins

  subroutine to_statepoint(this, filter_group)
    class(SphericalHarmonicsFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "sphericalharmonics")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "order", this % order)
    if (this % cosine == COSINE_SCATTER) then
      call write_dataset(filter_group, "cosine", "scatter")
    else
      call write_dataset(filter_group, "cosine", "particle")
    end if
  end subroutine to_statepoint

  function text_label(this, bin) result(label)
    class(SphericalHarmonicsFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    integer :: n, m

    do n = 0, this % order
      if (bin <= (n + 1)**2) then
        m = (bin - n**2 - 1) - n
        label = "Spherical harmonic expansion, Y" // trim(to_str(n)) // &
             "," // trim(to_str(m))
        exit
      end if
    end do
  end function text_label

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

end module tally_filter_sph_harm
