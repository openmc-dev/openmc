module mgxs_header

  use constants,       only: MAX_FILE_LEN, ZERO, ONE, TWO, PI
  use error,           only: fatal_error
  use list_header,     only: ListInt
  use material_header, only: material
  use math,            only: calc_pn, calc_rn, expand_harmonic, &
                             evaluate_legendre, find_angle
  use nuclide_header,  only: NuclideMGContainer, MaterialMacroXS
  use random_lcg,      only: prn
  use scattdata_header
  use string
  use xml_interface

!===============================================================================
! MGXS contains the base mgxs data for a nuclide/material
!===============================================================================

  type, abstract :: Mgxs
    character(12) :: name    ! name of dataset, e.g. 92235.03c
    integer       :: zaid    ! Z and A identifier, e.g. 92235
    real(8)       :: awr     ! Atomic Weight Ratio
    integer       :: listing ! index in xs_listings
    real(8)       :: kT      ! temperature in MeV (k*T)

    ! Fission information
    logical :: fissionable   ! mgxs object is fissionable?
    integer :: scatt_type    ! either legendre, histogram, or tabular.

  contains
    procedure(mgxs_print_),    deferred :: print    ! Writes object info
    procedure(mgxs_init_xml_), deferred :: init_xml ! Initialize the data
    procedure(mgxs_get_xs_),   deferred :: get_xs   ! Get the requested xs
    procedure(mgxs_combine_),  deferred :: combine  ! initializes object
    ! Sample the outgoing energy from a fission event
    procedure(mgxs_sample_fission_), deferred :: sample_fission_energy
    ! Sample the outgoing energy and angle from a scatter event
    procedure(mgxs_sample_scatter_), deferred :: sample_scatter
    ! Calculate the material specific MGXS data from the nuclides
    procedure(mgxs_calculate_xs_), deferred   :: calculate_xs
  end type Mgxs

  abstract interface
    subroutine mgxs_print_(this, unit)
      import Mgxs
      class(Mgxs),intent(in)     :: this
      integer, optional, intent(in) :: unit
    end subroutine mgxs_print_

    subroutine mgxs_init_xml_(this, node_xsdata, groups, get_kfiss, get_fiss, &
                              max_order)
      import Mgxs, Node
      class(Mgxs), intent(inout)      :: this        ! Working Object
      type(Node), pointer, intent(in) :: node_xsdata ! Data from MGXS xml
      integer, intent(in)             :: groups      ! Number of Energy groups
      logical, intent(in)             :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)             :: get_fiss    ! Should we get fiss data?
      integer, intent(in)             :: max_order ! Maximum requested order
    end subroutine mgxs_init_xml_

    function mgxs_get_xs_(this, xstype, gin, gout, uvw, mu, iazi, ipol) &
         result(xs)
      import Mgxs
      class(Mgxs), intent(in)       :: this
      character(*), intent(in)      :: xstype ! Cross Section Type
      integer, intent(in)           :: gin    ! Incoming Energy group
      integer, optional, intent(in) :: gout   ! Outgoing Group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      integer, optional, intent(in) :: iazi  ! Azimuthal Index
      integer, optional, intent(in) :: ipol  ! Polar Index
      real(8)                       :: xs     ! Resultant xs
    end function mgxs_get_xs_

    pure function mgxs_calc_f_(this, gin, gout, mu, uvw, iazi, ipol) result(f)
      import Mgxs
      class(Mgxs), intent(in)       :: this
      integer, intent(in)           :: gin   ! Incoming Energy Group
      integer, intent(in)           :: gout  ! Outgoing Energy Group
      real(8), intent(in)           :: mu    ! Angle of interest
      real(8), intent(in), optional :: uvw(3) ! Direction vector
      integer, intent(in), optional :: iazi ! Incoming Energy Group
      integer, intent(in), optional :: ipol ! Outgoing Energy Group
      real(8)                       :: f     ! Return value of f(mu)

    end function mgxs_calc_f_

    subroutine mgxs_combine_(this, mat, nuclides, groups, get_kfiss, get_fiss, &
                             max_order, scatt_type)
      import Mgxs, Material, NuclideMGContainer, MAX_LINE_LEN
      class(Mgxs),           intent(inout) :: this ! The Mgxs to initialize
      type(Material), pointer,  intent(in) :: mat  ! base material
      type(NuclideMGContainer), intent(in) :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                  :: groups ! Number of E groups
      logical, intent(in)                  :: get_kfiss ! Should we get kfiss data?
      logical, intent(in)                  :: get_fiss ! Should we get fiss data?
      integer, intent(in)                  :: max_order ! Maximum requested order
      integer, intent(in)                  :: scatt_type ! Legendre or Tabular Scatt?
    end subroutine mgxs_combine_

    function mgxs_sample_fission_(this, gin, uvw) result(gout)
      import Mgxs
      class(Mgxs), intent(in) :: this   ! Data to work with
      integer, intent(in)     :: gin    ! Incoming energy group
      real(8), intent(in)     :: uvw(3) ! Particle Direction
      integer                 :: gout   ! Sampled outgoing group

    end function mgxs_sample_fission_

    subroutine mgxs_sample_scatter_(this, uvw, gin, gout, mu, wgt)
      import Mgxs
      class(Mgxs),    intent(in)    :: this
      real(8),        intent(in)    :: uvw(3) ! Incoming neutron direction
      integer,        intent(in)    :: gin    ! Incoming neutron group
      integer,        intent(out)   :: gout   ! Sampled outgoin group
      real(8),        intent(out)   :: mu     ! Sampled change in angle
      real(8),        intent(inout) :: wgt    ! Particle weight
    end subroutine mgxs_sample_scatter_

    subroutine mgxs_calculate_xs_(this, gin, uvw, xs)
      import Mgxs, MaterialMacroXS
      class(Mgxs),           intent(in)    :: this
      integer,               intent(in)    :: gin         ! Incoming neutron group
      real(8),               intent(in)    :: uvw(3)      ! Incoming neutron direction
      type(MaterialMacroXS), intent(inout) :: xs
    end subroutine mgxs_calculate_xs_
  end interface

end module mgxs_header