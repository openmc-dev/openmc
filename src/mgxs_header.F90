module mgxs_header

  use constants,       only: MAX_FILE_LEN, ZERO, ONE, TWO, PI
  use error,           only: fatal_error
  use, intrinsic :: ISO_FORTRAN_ENV,      only: OUTPUT_UNIT
  use list_header,     only: ListInt
  use material_header, only: material
  use math,            only: calc_pn, calc_rn, expand_harmonic, &
                             evaluate_legendre, find_angle
  use nuclide_header,  only: MaterialMacroXS
  use random_lcg,      only: prn
  use scattdata_header
  use string
  use xml_interface

!===============================================================================
! MGXS contains the base mgxs data for a nuclide/material
!===============================================================================

  type, abstract :: Mgxs
    character(len=104) :: name    ! name of dataset, e.g. 92235.03c
    integer            :: zaid    ! Z and A identifier, e.g. 92235
    real(8)            :: awr     ! Atomic Weight Ratio
    integer            :: listing ! index in xs_listings
    real(8)            :: kT      ! temperature in MeV (k*T)

    ! Fission information
    logical :: fissionable   ! mgxs object is fissionable?
    integer :: scatt_type    ! either legendre, histogram, or tabular.

  contains
    procedure(mgxs_init_file_), deferred :: init_file ! Initialize the data
    procedure(mgxs_print_),    deferred  :: print     ! Writes object info
    procedure(mgxs_get_xs_),   deferred  :: get_xs    ! Get the requested xs
    procedure(mgxs_combine_),  deferred  :: combine   ! initializes object
    ! Sample the outgoing energy from a fission event
    procedure(mgxs_sample_fission_), deferred :: sample_fission_energy
    ! Sample the outgoing energy and angle from a scatter event
    procedure(mgxs_sample_scatter_), deferred :: sample_scatter
    ! Calculate the material specific MGXS data from the nuclides
    procedure(mgxs_calculate_xs_), deferred   :: calculate_xs
  end type Mgxs

!===============================================================================
! MGXSCONTAINER pointer array for storing Nuclides
!===============================================================================

  type MgxsContainer
    class(Mgxs), pointer :: obj
  end type MgxsContainer

!===============================================================================
! Interfaces for MGXS
!===============================================================================

  abstract interface
    subroutine mgxs_init_file_(this,node_xsdata,groups,get_kfiss,get_fiss, &
                               max_order,i_listing)
      import Mgxs, Node
      class(Mgxs), intent(inout)      :: this        ! Working Object
      type(Node), pointer, intent(in) :: node_xsdata ! Data from MGXS xml
      integer, intent(in)             :: groups      ! Number of Energy groups
      logical, intent(in)             :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)             :: get_fiss    ! Should we get fiss data?
      integer, intent(in)             :: max_order   ! Maximum requested order
      integer, intent(in)             :: i_listing   ! Index of listings array
    end subroutine mgxs_init_file_

    subroutine mgxs_print_(this, unit)
      import Mgxs
      class(Mgxs),intent(in)     :: this
      integer, optional, intent(in) :: unit
    end subroutine mgxs_print_

    pure function mgxs_get_xs_(this,xstype,gin,gout,uvw,mu) result(xs)
      import Mgxs
      class(Mgxs), intent(in)       :: this
      character(*), intent(in)      :: xstype ! Cross Section Type
      integer, intent(in)           :: gin    ! Incoming Energy group
      integer, optional, intent(in) :: gout   ! Outgoing Group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      real(8)                       :: xs     ! Resultant xs
    end function mgxs_get_xs_

    pure function mgxs_calc_f_(this,gin,gout,mu,uvw,iazi,ipol) result(f)
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

    subroutine mgxs_combine_(this,mat,nuclides,groups,max_order,scatt_type, &
                             i_listing)
      import Mgxs, Material, MgxsContainer
      class(Mgxs),           intent(inout) :: this ! The Mgxs to initialize
      type(Material), pointer,  intent(in) :: mat  ! base material
      type(MgxsContainer), intent(in)      :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                  :: groups     ! Number of E groups
      integer, intent(in)                  :: max_order  ! Maximum requested order
      integer, intent(in)                  :: scatt_type ! Legendre or Tabular Scatt?
      integer, intent(in)                  :: i_listing  ! Index in listings
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
      integer,               intent(in)    :: gin    ! Incoming neutron group
      real(8),               intent(in)    :: uvw(3) ! Incoming neutron direction
      type(MaterialMacroXS), intent(inout) :: xs     ! Resultant Mgxs Data
    end subroutine mgxs_calculate_xs_
  end interface

!===============================================================================
! MGXSISO contains the base MGXS data specifically for
! isotropically weighted MGXS
!===============================================================================

  type, extends(Mgxs) :: MgxsIso

    ! Microscopic cross sections
    real(8), allocatable :: total(:)        ! total cross section
    real(8), allocatable :: absorption(:)   ! absorption cross section
    class(ScattData), allocatable :: scatter ! scattering information
    real(8), allocatable :: nu_fission(:)   ! fission matrix (Gout x Gin)
    real(8), allocatable :: k_fission(:)    ! kappa-fission
    real(8), allocatable :: fission(:)      ! neutron production
    real(8), allocatable :: chi(:,:)        ! Fission Spectra

  contains
    procedure :: init_file   => mgxsiso_init_file ! Initialize Nuclidic MGXS Data
    procedure :: print       => mgxsiso_print   ! Writes nuclide info
    procedure :: get_xs      => mgxsiso_get_xs  ! Gets Size of Data w/in Object
    procedure :: combine     => mgxsiso_combine ! inits object
    procedure :: sample_fission_energy => mgxsiso_sample_fission_energy
    procedure :: sample_scatter => mgxsiso_sample_scatter
    procedure :: calculate_xs => mgxsiso_calculate_xs
  end type MgxsIso

!===============================================================================
! MGXSANGLE contains the base MGXS data specifically for
! angular flux weighted MGXS
!===============================================================================

  type, extends(Mgxs) :: MgxsAngle

    ! Microscopic cross sections
    real(8), allocatable :: total(:,:,:)        ! total cross section
    real(8), allocatable :: absorption(:,:,:)   ! absorption cross section
    type(ScattDataContainer), allocatable :: scatter(:,:) ! scattering information
    real(8), allocatable :: nu_fission(:,:,:)   ! fission matrix (Gout x Gin)
    real(8), allocatable :: k_fission(:,:,:)    ! kappa-fission
    real(8), allocatable :: fission(:,:,:)      ! neutron production
    real(8), allocatable :: chi(:,:,:,:)        ! Fission Spectra
    ! In all cases, right-most indices are theta, phi
    integer              :: n_pol         ! Number of polar angles
    integer              :: n_azi         ! Number of azimuthal angles
    real(8), allocatable :: polar(:)     ! polar angles
    real(8), allocatable :: azimuthal(:) ! azimuthal angles

  contains
    procedure :: init_file   => mgxsang_init_file ! Initialize Nuclidic MGXS Data
    procedure :: print       => mgxsang_print   ! Writes nuclide info
    procedure :: get_xs      => mgxsang_get_xs  ! Gets Size of Data w/in Object
    procedure :: combine     => mgxsang_combine ! inits object
    procedure :: sample_fission_energy => mgxsang_sample_fission_energy
    procedure :: sample_scatter => mgxsang_sample_scatter
    procedure :: calculate_xs => mgxsang_calculate_xs
  end type MgxsAngle

  contains

!===============================================================================
! MGXS*_INIT reads in the data from the XML file.  At the point of entry
! the file would have been opened and metadata read.  This routine begins with
! the xsdata object node itself.
!===============================================================================

    subroutine mgxs_init_file(this,node_xsdata,i_listing)
      class(Mgxs), intent(inout)      :: this        ! Working Object
      type(Node), pointer, intent(in) :: node_xsdata ! Data from MGXS xml
      integer, intent(in)             :: i_listing   ! Index in listings array

      character(MAX_LINE_LEN) :: temp_str

      ! Load the nuclide metadata
      call get_node_value(node_xsdata, "name", this % name)
      this % name = to_lower(this % name)
      if (check_for_node(node_xsdata, "kT")) then
        call get_node_value(node_xsdata, "kT", this % kT)
      else
        this % kT = ZERO
      end if
      if (check_for_node(node_xsdata, "zaid")) then
        call get_node_value(node_xsdata, "zaid", this % zaid)
      else
        this % zaid = 0
      end if
      if (check_for_node(node_xsdata, "scatt_type")) then
        call get_node_value(node_xsdata, "scatt_type", temp_str)
        temp_str = trim(to_lower(temp_str))
        if (temp_str == 'legendre') then
          this % scatt_type = ANGLE_LEGENDRE
        else if (temp_str == 'histogram') then
          this % scatt_type = ANGLE_HISTOGRAM
        else if (temp_str == 'tabular') then
          this % scatt_type = ANGLE_TABULAR
        else
          call fatal_error("Invalid Scatt Type Option!")
        end if
      else
        this % scatt_type = ANGLE_LEGENDRE
      end if

      if (check_for_node(node_xsdata, "fissionable")) then
        call get_node_value(node_xsdata, "fissionable", temp_str)
        temp_str = to_lower(temp_str)
        if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') then
          this % fissionable = .true.
        else
          this % fissionable = .false.
        end if
      else
        call fatal_error("Fissionable element must be set!")
      end if

      ! Keep track of what listing is associated with this nuclide
      this % listing = i_listing

    end subroutine mgxs_init_file

    subroutine mgxsiso_init_file(this,node_xsdata,groups,get_kfiss,get_fiss, &
                                 max_order,i_listing)
      class(MgxsIso), intent(inout)   :: this        ! Working Object
      type(Node), pointer, intent(in) :: node_xsdata ! Data from MGXS xml
      integer, intent(in)             :: groups      ! Number of Energy groups
      logical, intent(in)             :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)             :: get_fiss    ! Need fiss data?
      integer, intent(in)             :: max_order   ! Maximum requested order
      integer, intent(in)             :: i_listing   ! Index in listings array

      type(Node), pointer     :: node_legendre_mu
      character(MAX_LINE_LEN) :: temp_str
      logical                 :: enable_leg_mu
      real(8), allocatable    :: temp_arr(:), temp_2d(:,:)
      real(8), allocatable    :: temp_mult(:,:)
      real(8), allocatable    :: scatt_coeffs(:,:,:)
      real(8), allocatable    :: input_scatt(:,:,:)
      real(8), allocatable    :: temp_scatt(:,:,:)
      real(8)                 :: dmu, mu, norm
      integer                 :: order, order_dim, gin, gout, l, arr_len
      integer                 :: legendre_mu_points, imu

      ! Call generic data gathering routine (will populate the metadata)
      call mgxs_init_file(this,node_xsdata,i_listing)

      ! Load the more specific data
      allocate(this % nu_fission(groups))
      allocate(this % chi(groups,groups))
      if (this % fissionable) then
        if (check_for_node(node_xsdata,"chi")) then
          ! Chi was provided, that means they are giving chi and nu-fission
          ! vectors
          ! Get chi
          allocate(temp_arr(1 * groups))
          call get_node_array(node_xsdata,"chi",temp_arr)
          do gin = 1, groups
            do gout = 1, groups
              this % chi(gout,gin) = temp_arr(gout)
            end do
            ! Normalize chi so its CDF goes to 1
            this % chi(:,gin) = this % chi(:,gin) / sum(this % chi(:,gin))
          end do
          deallocate(temp_arr)

          ! Get nu_fission (as a vector)
          if (check_for_node(node_xsdata,"nu_fission")) then
            call get_node_array(node_xsdata,"nu_fission",this % nu_fission)
          else
            call fatal_error("If fissionable, must provide nu_fission!")
          end if

        else
          ! chi isnt provided but is within nu_fission, existing as a matrix
          ! So, get nu_fission (as a matrix)
          if (check_for_node(node_xsdata,"nu_fission")) then
            allocate(temp_arr(groups*groups))
            call get_node_array(node_xsdata,"nu_fission",temp_arr)
            allocate(temp_2d(groups,groups))
            temp_2d = reshape(temp_arr,(/groups,groups/))
            deallocate(temp_arr)
          else
            call fatal_error("If fissionable, must provide nu_fission!")
          end if

          ! Set the vector nu-fission from the matrix nu-fission
          do gin = 1, groups
            this % nu_fission(gin) = sum(temp_2d(:,gin))
          end do

          ! Now pull out information needed for chi
          this % chi = temp_2d
          ! Normalize chi so its CDF goes to 1
          do gin = 1, groups
            this % chi(:,gin) = this % chi(:,gin) / sum(this % chi(:,gin))
          end do
          deallocate(temp_2d)
        end if
        ! If we have a need* for the fission and kappa-fission x/s, get them
        ! (*Need is defined as will be using it to tally)
        if (get_fiss) then
          allocate(this % fission(groups))
          if (check_for_node(node_xsdata,"fission")) then
            call get_node_array(node_xsdata,"fission",this % fission)
          else
            call fatal_error("Fission data missing, required due to fission&
                             & tallies in tallies.xml file!")
          end if
        end if
        if (get_kfiss) then
          allocate(this % k_fission(groups))
          if (check_for_node(node_xsdata,"kappa_fission")) then
            call get_node_array(node_xsdata,"kappa_fission",this % k_fission)
          else
            call fatal_error("kappa_fission data missing, required due to &
                             &kappa-fission tallies in tallies.xml file!")
          end if
        end if
      else
        this % nu_fission = ZERO
        this % chi = ZERO
      end if

      allocate(this % absorption(groups))
      if (check_for_node(node_xsdata,"absorption")) then
        call get_node_array(node_xsdata,"absorption",this % absorption)
      else
        call fatal_error("Must provide absorption!")
      end if

      ! Get multiplication data if present
      allocate(temp_mult(groups, groups))
      if (check_for_node(node_xsdata,"multiplicity")) then
        arr_len = get_arraysize_double(node_xsdata,"multiplicity")
        if (arr_len == groups * groups) then
          allocate(temp_arr(arr_len))
          call get_node_array(node_xsdata,"multiplicity",temp_arr)
          temp_mult = reshape(temp_arr, (/groups, groups/))
          deallocate(temp_arr)
        else
          call fatal_error("Multiplicity length not same as number of groups&
                           & squared!")
        end if
      else
        temp_mult = ONE
      end if

      ! Get scattering treatment information
      ! Tabular_legendre tells us if we are to treat the provided
      ! Legendre polynomials as tabular data (if enable is true) or leaving
      ! them as Legendres (if enable is false, or the default)

      ! Set the default (leave as Legendre polynomials)
      enable_leg_mu = .false.
      if (check_for_node(node_xsdata,"tabular_legendre")) then
        call get_node_ptr(node_xsdata,"tabular_legendre",node_legendre_mu)
        if (check_for_node(node_legendre_mu, "enable")) then
          call get_node_value(node_legendre_mu,"enable",temp_str)
          temp_str = trim(to_lower(temp_str))
          if (temp_str == 'true' .or. temp_str == '1') then
            enable_leg_mu = .true.
          elseif (temp_str == 'false' .or. temp_str == '0') then
            enable_leg_mu = .false.
          else
            call fatal_error("Unrecognized tabular_legendre/enable: " // temp_str)
          end if
        end if
        ! Ok, so if we need to convert to a tabular form, get the user provided
        ! number of points
        if (enable_leg_mu) then
          if (check_for_node(node_legendre_mu,"num_points")) then
            call get_node_value(node_legendre_mu,"num_points", &
                 legendre_mu_points)
            if (legendre_mu_points <= 0) &
                 call fatal_error("num_points element must be positive&
                                  & and non-zero!")
          else
            ! Set the default number of points (0.0625 spacing)
            legendre_mu_points = 33
          end if
        end if
      end if

      ! Get the library's value for the order
      if (check_for_node(node_xsdata,"order")) then
        call get_node_value(node_xsdata,"order",order)
      else
        call fatal_error("Order Must Be Provided!")
      end if

      ! Before retrieving the data, store the dimensionality of the data in
      ! order_dim.  For Legendre data, we usually refer to it as Pn where
      ! n is the order.  However Pn has n+1 sets of points (since you need to
      ! the count the P0 moment).  Adjust for that.  Histogram and Tabular
      ! formats dont need this adjustment.
      if (this % scatt_type == ANGLE_LEGENDRE) then
        order_dim = order + 1
      else
        order_dim = order
      end if

      ! The input is gathered in the more user-friendly facing format of
      ! Gout x Gin x Order.  We will get it in that format in input_scatt,
      ! but then need to convert it to a more useful ordering for processing
      ! (Order x Gout x Gin).
      allocate(input_scatt(groups, groups, order_dim))
      if (check_for_node(node_xsdata,"scatter")) then
        allocate(temp_arr(groups * groups * order_dim))
        call get_node_array(node_xsdata,"scatter",temp_arr)
        input_scatt = reshape(temp_arr,(/groups,groups,order_dim/))
        deallocate(temp_arr)

        ! Compare the number of orders given with the maximum order of the
        ! problem.  Strip off the supefluous orders if needed.
        if (this % scatt_type == ANGLE_LEGENDRE) then
          order = min(order_dim - 1, max_order)
          order_dim = order + 1
        end if
        allocate(temp_scatt(groups,groups,order_dim))
        temp_scatt(:,:,:) = input_scatt(:,:,1:order_dim)

        ! Take input format (groups, groups, order) and convert to
        ! the more useful format needed for scattdata: (order, groups, groups)
        ! However, if scatt_type was ANGLE_LEGENDRE (i.e., the data was
        ! provided as Legendre coefficients), and the user requested that
        ! these legendres be converted to tabular form (note this is also
        ! the default behavior), convert that now.
        if (this % scatt_type == ANGLE_LEGENDRE .and. enable_leg_mu) then
          ! Convert input parameters to what we need for the rest.
          this % scatt_type = ANGLE_TABULAR
          order_dim = legendre_mu_points
          order = order_dim
          dmu = TWO / real(order - 1,8)

          allocate(scatt_coeffs(order_dim,groups,groups))
          do gin = 1, groups
            do gout = 1, groups
              norm = ZERO
              do imu = 1, order_dim
                if (imu == 1) then
                  mu = -ONE
                else if (imu == order_dim) then
                  mu = ONE
                else
                  mu = -ONE + real(imu - 1,8) * dmu
                end if
                scatt_coeffs(imu,gout,gin) = &
                     evaluate_legendre(temp_scatt(gout,gin,:),mu)
                ! Ensure positivity of distribution
                if (scatt_coeffs(imu,gout,gin) < ZERO) &
                     scatt_coeffs(imu,gout,gin) = ZERO
                ! And accrue the integral
                if (imu > 1) then
                  norm = norm + HALF * dmu * (scatt_coeffs(imu-1,gout,gin) + &
                                              scatt_coeffs(imu,gout,gin))
                end if
              end do
              ! Now that we have the integral, lets ensure that the distribution
              ! is normalized such that it preserves the original scattering xs
              if (norm > ZERO) then
                scatt_coeffs(:,gout,gin) = scatt_coeffs(:,gout,gin) * &
                     temp_scatt(gout,gin,1) / norm
              end if
            end do
          end do
        else
          ! Sticking with current representation, carry forward but change
          ! the array ordering
          allocate(scatt_coeffs(order_dim,groups,groups))
          do gin = 1, groups
            do gout = 1, groups
              do l = 1, order_dim
                scatt_coeffs(l,gout,gin) = temp_scatt(gout,gin,l)
              end do
            end do
          end do
        end if
        deallocate(temp_scatt)
      else
        call fatal_error("Must provide scatter!")
      end if

      ! Allocate and initialize our ScattData Object.
      if (this % scatt_type == ANGLE_HISTOGRAM) then
        allocate(ScattDataHistogram :: this % scatter)
      else if (this % scatt_type == ANGLE_TABULAR) then
        allocate(ScattDataTabular :: this % scatter)
      else if (this % scatt_type == ANGLE_LEGENDRE) then
        allocate(ScattDataLegendre :: this % scatter)
      end if

      ! Initialize the ScattData Object
      call this % scatter % init(temp_mult, scatt_coeffs)

      ! Check sigA to ensure it is not 0 since it is
      ! often divided by in the tally routines
      ! (This may happen with Helium data)
      do gin = 1, groups
        if (this % absorption(gin) == ZERO) this % absorption(gin) = 1E-10_8
      end do

      ! Get, or infer, total xs data.
      allocate(this % total(groups))
      if (check_for_node(node_xsdata,"total")) then
        call get_node_array(node_xsdata,"total",this % total)
      else
        this % total = this % absorption + this % scatter % scattxs
      end if

      ! Deallocate temporaries for the next material
      deallocate(input_scatt,scatt_coeffs,temp_mult)

      ! Finally, check sigT to ensure it is not 0 since it is
      ! often divided by in the tally routines
      do gin = 1, groups
        if (this % total(gin) == ZERO) this % total(gin) = 1E-10_8
      end do


    end subroutine mgxsiso_init_file

    subroutine mgxsang_init_file(this,node_xsdata,groups,get_kfiss,get_fiss, &
                                 max_order,i_listing)
      class(MgxsAngle), intent(inout) :: this        ! Working Object
      type(Node), pointer, intent(in) :: node_xsdata ! Data from MGXS xml
      integer, intent(in)             :: groups      ! Number of Energy groups
      logical, intent(in)             :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)             :: get_fiss    ! Should we get fiss data?
      integer, intent(in)             :: max_order   ! Maximum requested order
      integer, intent(in)             :: i_listing   ! Index in listings array

      type(Node), pointer     :: node_legendre_mu
      character(MAX_LINE_LEN) :: temp_str
      logical                 :: enable_leg_mu
      real(8), allocatable    :: temp_arr(:), temp_4d(:,:,:,:)
      real(8), allocatable    :: temp_mult(:,:,:,:)
      real(8), allocatable    :: scatt_coeffs(:,:,:,:,:)
      real(8), allocatable    :: input_scatt(:,:,:,:,:)
      real(8), allocatable    :: temp_scatt(:,:,:,:,:)
      real(8)                 :: dmu, mu, norm, dangle
      integer                 :: order, order_dim, gin, gout, l, arr_len
      integer                 :: legendre_mu_points, imu, ipol, iazi

      ! Call generic data gathering routine (will populate the metadata)
      call mgxs_init_file(this,node_xsdata,i_listing)

      if (check_for_node(node_xsdata, "num_polar")) then
        call get_node_value(node_xsdata, "num_polar", this % n_pol)
      else
        call fatal_error("num_polar Must Be Provided!")
      end if

      if (check_for_node(node_xsdata, "num_azimuthal")) then
        call get_node_value(node_xsdata, "num_azimuthal", this % n_azi)
      else
        call fatal_error("num_azimuthal Must Be Provided!")
      end if

      ! Load angle data, if present (else equally spaced)
      allocate(this % polar(this % n_pol))
      allocate(this % azimuthal(this % n_azi))
      if (check_for_node(node_xsdata, "polar")) then
        call fatal_error("User-Specified polar angle bins not yet supported!")
        ! When this feature is supported, this line will be activated
        call get_node_array(node_xsdata, "polar", this % polar)
      else
        dangle = PI / real(this % n_pol,8)
        do ipol = 1, this % n_pol
          this % polar(ipol) = (real(ipol,8) - HALF) * dangle
        end do
      end if
      if (check_for_node(node_xsdata, "azimuthal")) then
        call fatal_error("User-Specified azimuthal angle bins not yet supported!")
        ! When this feature is supported, this line will be activated
        call get_node_array(node_xsdata, "azimuthal", this % azimuthal)
      else
        dangle = TWO * PI / real(this % n_azi,8)
        do iazi = 1, this % n_azi
          this % azimuthal(iazi) = -PI + (real(iazi,8) - HALF) * dangle
        end do
      end if

      ! Load the more specific data
      allocate(this % nu_fission(groups,this % n_azi,this % n_pol))
      allocate(this % chi(groups,groups,this % n_azi,this % n_pol))
      if (this % fissionable) then
        if (check_for_node(node_xsdata,"chi")) then
          ! Chi was provided, that means they are giving chi and nu-fission
          ! vectors
          ! Get chi
          allocate(temp_arr(1 * groups * this % n_azi * this % n_pol))
          call get_node_array(node_xsdata,"chi",temp_arr)
          ! Initialize counter for temp_arr
          l = 0
          gin = 1
          do ipol = 1, this % n_pol
            do iazi = 1, this % n_azi
              do gout = 1, groups
                l = l + 1
                this % chi(gout,gin,iazi,ipol) = temp_arr(l)
              end do
              ! Normalize chi so its CDF goes to 1
              this % chi(:,gin,iazi,ipol) = this % chi(:,gin,iazi,ipol) / &
                   sum(this % chi(:,gin,iazi,ipol))
            end do
          end do

          ! Now set all the other gin values
          do ipol = 1, this % n_pol
            do iazi = 1, this % n_azi
              do gin = 2, groups
                this % chi(:,gin,iazi,ipol) = this % chi(:,1,iazi,ipol)
              end do
            end do
          end do
          deallocate(temp_arr)

          ! Get nu_fission (as a vector)
          if (check_for_node(node_xsdata,"nu_fission")) then
            allocate(temp_arr(groups * this % n_azi * this % n_pol))
            call get_node_array(node_xsdata,"nu_fission",temp_arr)
            this % nu_fission = reshape(temp_arr,(/groups,this % n_azi,this % n_pol/))
            deallocate(temp_arr)
          else
            call fatal_error("If fissionable, must provide nu_fission!")
          end if

        else
          ! chi isnt provided but is within nu_fission, existing as a matrix
          ! So, get nu_fission (as a matrix)
          if (check_for_node(node_xsdata,"nu_fission")) then
            allocate(temp_arr(groups * groups * this % n_azi * this % n_pol))
            call get_node_array(node_xsdata,"nu_fission",temp_arr)
            allocate(temp_4d(groups,groups,this % n_azi,this % n_pol))
            temp_4d = reshape(temp_arr,(/groups,groups,this % n_azi,this % n_pol/))
            deallocate(temp_arr)
          else
            call fatal_error("If fissionable, must provide nu_fission!")
          end if

          ! Set the vector nu-fission from the matrix nu-fission
          do ipol = 1, this % n_pol
            do iazi = 1, this % n_azi
              do gin = 1, groups
                this % nu_fission(gin,iazi,ipol) = sum(temp_4d(:,gin,iazi,ipol))
              end do
            end do
          end do

          ! Now pull out information needed for chi
          this % chi = temp_4d
          ! Normalize chi so its CDF goes to 1
          do ipol = 1, this % n_pol
            do iazi = 1, this % n_azi
              do gin = 1, groups
                this % chi(:,gin,iazi,ipol) = this % chi(:,gin,iazi,ipol) / &
                     sum(this % chi(:,gin,iazi,ipol))
              end do
            end do
          end do
          deallocate(temp_4d)
        end if

        ! If we have a need* for the fission and kappa-fission x/s, get them
        ! (*Need is defined as will be using it to tally)
        if (get_fiss) then
          if (check_for_node(node_xsdata,"fission")) then
            allocate(temp_arr(groups * this % n_azi * this % n_pol))
            call get_node_array(node_xsdata,"fission",temp_arr)
            allocate(this % fission(groups,this % n_azi,this % n_pol))
            this % fission = reshape(temp_arr,(/groups,this % n_azi,this % n_pol/))
            deallocate(temp_arr)
          else
            call fatal_error("Fission data missing, required due to fission&
                             & tallies in tallies.xml file!")
          end if
        end if
        if (get_kfiss) then
          if (check_for_node(node_xsdata,"kappa_fission")) then
            allocate(temp_arr(groups * this % n_azi * this % n_pol))
            call get_node_array(node_xsdata,"kappa_fission",temp_arr)
            allocate(this % k_fission(groups,this % n_azi,this % n_pol))
            this % k_fission = reshape(temp_arr,(/groups, this % n_azi,this % n_pol/))
            deallocate(temp_arr)
          else
            call fatal_error("kappa_fission data missing, required due to &
                             &kappa-fission tallies in tallies.xml file!")
          end if
        end if
      else
        this % nu_fission = ZERO
        this % chi = ZERO
      end if

      if (check_for_node(node_xsdata,"absorption")) then
        allocate(temp_arr(groups * this % n_azi * this % n_pol))
        call get_node_array(node_xsdata,"absorption",temp_arr)
        allocate(this % absorption(groups,this % n_azi,this % n_pol))
        this % absorption = reshape(temp_arr,(/groups,this % n_azi,this % n_pol/))
        deallocate(temp_arr)
      else
        call fatal_error("Must provide absorption!")
      end if

      ! Get multiplication data if present
      allocate(temp_mult(groups,groups,this % n_azi,this % n_pol))
      if (check_for_node(node_xsdata,"multiplicity")) then
        arr_len = get_arraysize_double(node_xsdata,"multiplicity")
        if (arr_len == groups * groups * this % n_azi * this % n_pol) then
          allocate(temp_arr(arr_len))
          call get_node_array(node_xsdata,"multiplicity",temp_arr)
          temp_mult = reshape(temp_arr,(/groups,groups,this % n_azi,this % n_pol/))
          deallocate(temp_arr)
        else
          call fatal_error("Multiplicity length not same as number of groups&
                           & squared!")
        end if
      else
        temp_mult = ONE
      end if

      ! Get scattering treatment information
      ! Tabular_legendre tells us if we are to treat the provided
      ! Legendre polynomials as tabular data (if enable is true) or leaving
      ! them as Legendres (if enable is false, or the default)

      ! Set the default (leave as Legendre polynomials)
      enable_leg_mu = .false.
      if (check_for_node(node_xsdata,"tabular_legendre")) then
        call get_node_ptr(node_xsdata,"tabular_legendre",node_legendre_mu)
        if (check_for_node(node_legendre_mu, "enable")) then
          call get_node_value(node_legendre_mu,"enable",temp_str)
          temp_str = trim(to_lower(temp_str))
          if (temp_str == 'true' .or. temp_str == '1') then
            enable_leg_mu = .true.
          elseif (temp_str == 'false' .or. temp_str == '0') then
            enable_leg_mu = .false.
          else
            call fatal_error("Unrecognized tabular_legendre/enable: " // temp_str)
          end if
        end if
        ! Ok, so if we need to convert to a tabular form, get the user provided
        ! number of points
        if (enable_leg_mu) then
          if (check_for_node(node_legendre_mu,"num_points")) then
            call get_node_value(node_legendre_mu,"num_points", &
                 legendre_mu_points)
            if (legendre_mu_points <= 0) &
                 call fatal_error("num_points element must be positive&
                                  & and non-zero!")
          else
            ! Set the default number of points (0.0625 spacing)
            legendre_mu_points = 33
          end if
        end if
      end if

      ! Get the library's value for the order
      if (check_for_node(node_xsdata,"order")) then
        call get_node_value(node_xsdata,"order",order)
      else
        call fatal_error("Order Must Be Provided!")
      end if

      ! Before retrieving the data, store the dimensionality of the data in
      ! order_dim.  For Legendre data, we usually refer to it as Pn where
      ! n is the order.  However Pn has n+1 sets of points (since you need to
      ! the count the P0 moment).  Adjust for that.  Histogram and Tabular
      ! formats dont need this adjustment.
      if (this % scatt_type == ANGLE_LEGENDRE) then
        order_dim = order + 1
      else
        order_dim = order
      end if

      ! The input is gathered in the more user-friendly facing format of
      ! Gout x Gin x Order x Azi x Pol.  We will get it in that format in
      ! input_scatt, but then need to convert it to a more useful ordering
      ! for processing (Order x Gout x Gin x Azi x Pol).
      allocate(input_scatt(groups,groups,order_dim,this % n_azi,this % n_pol))
      if (check_for_node(node_xsdata,"scatter")) then
        allocate(temp_arr(groups * groups * order_dim * this % n_azi * &
                          this % n_pol))
        call get_node_array(node_xsdata,"scatter",temp_arr)
        input_scatt = reshape(temp_arr,(/groups,groups,order_dim,this % n_azi, &
                                        this % n_pol/))
        deallocate(temp_arr)

        ! Compare the number of orders given with the maximum order of the
        ! problem.  Strip off the supefluous orders if needed.
        if (this % scatt_type == ANGLE_LEGENDRE) then
          order = min(order_dim - 1, max_order)
          order_dim = order + 1
        end if

        allocate(temp_scatt(groups,groups,order_dim,this % n_azi,this % n_pol))
        temp_scatt(:,:,:,:,:) = input_scatt(:,:,1:order_dim,:,:)

        ! Take input format (groups, groups, order) and convert to
        ! the more useful format needed for scattdata: (order, groups, groups)
        ! However, if scatt_type was ANGLE_LEGENDRE (i.e., the data was
        ! provided as Legendre coefficients), and the user requested that
        ! these legendres be converted to tabular form (note this is also
        ! the default behavior), convert that now.
        if (this % scatt_type == ANGLE_LEGENDRE .and. enable_leg_mu) then

          ! Convert input parameters to what we need for the rest.
          this % scatt_type = ANGLE_TABULAR
          order_dim = legendre_mu_points
          order = order_dim
          dmu = TWO / real(order - 1,8)

          allocate(scatt_coeffs(order_dim,groups,groups,this % n_azi,this % n_pol))
          do ipol = 1, this % n_pol
            do iazi = 1, this % n_azi
              do gin = 1, groups
                do gout = 1, groups
                  norm = ZERO
                  do imu = 1, order_dim
                    if (imu == 1) then
                      mu = -ONE
                    else if (imu == order_dim) then
                      mu = ONE
                    else
                      mu = -ONE + real(imu - 1,8) * dmu
                    end if
                    scatt_coeffs(imu,gout,gin,iazi,ipol) = &
                         evaluate_legendre(temp_scatt(gout,gin,:,iazi,ipol),mu)
                    ! Ensure positivity of distribution
                    if (scatt_coeffs(imu,gout,gin,iazi,ipol) < ZERO) &
                         scatt_coeffs(imu,gout,gin,iazi,ipol) = ZERO
                    ! And accrue the integral
                    if (imu > 1) then
                      norm = norm + HALF * dmu * &
                           (scatt_coeffs(imu-1,gout,gin,iazi,ipol) + &
                            scatt_coeffs(imu,gout,gin,iazi,ipol))
                    end if
                  end do
                  ! Now that we have the integral, lets ensure that the distribution
                  ! is normalized such that it preserves the original scattering xs
                  if (norm > ZERO) then
                    scatt_coeffs(:,gout,gin,iazi,ipol) = &
                         scatt_coeffs(:,gout,gin,iazi,ipol) * &
                         temp_scatt(gout,gin,1,iazi,ipol) / norm
                  end if
                end do
              end do
            end do
          end do
        else
          ! Sticking with current representation, carry forward but change
          ! the array ordering
          allocate(scatt_coeffs(order_dim,groups,groups,this % n_azi,this % n_pol))
          do ipol = 1, this % n_pol
            do iazi = 1, this % n_azi
              do gin = 1, groups
                do gout = 1, groups
                  do l = 1, order_dim
                    scatt_coeffs(l,gout,gin,iazi,ipol) = &
                         temp_scatt(gout,gin,l,iazi,ipol)
                  end do
                end do
              end do
            end do
          end do
        end if
        deallocate(temp_scatt)
      else
        call fatal_error("Must provide scatter!")
      end if

      allocate(this % scatter(this % n_azi, this % n_pol))
      do ipol = 1, this % n_pol
        do iazi = 1, this % n_azi
          ! Allocate and initialize our ScattData Object.
          if (this % scatt_type == ANGLE_HISTOGRAM) then
            allocate(ScattDataHistogram :: this % scatter(iazi,ipol) % obj)
          else if (this % scatt_type == ANGLE_TABULAR) then
            allocate(ScattDataTabular :: this % scatter(iazi,ipol) % obj)
          else if (this % scatt_type == ANGLE_LEGENDRE) then
            allocate(ScattDataLegendre :: this % scatter(iazi,ipol) % obj)
          end if

          ! Initialize the ScattData Object
          call this % scatter(iazi,ipol) % obj % init(&
               temp_mult(:,:,iazi,ipol), scatt_coeffs(:,:,:,iazi,ipol))
        end do
      end do
      ! Deallocate temporaries for the next material
      deallocate(input_scatt,scatt_coeffs,temp_mult)

      allocate(this % total(groups,this % n_azi,this % n_pol))
      if (check_for_node(node_xsdata,"total")) then
        allocate(temp_arr(groups * this % n_azi * this % n_pol))
        call get_node_array(node_xsdata,"total",temp_arr)
        this % total = reshape(temp_arr,(/groups,this % n_azi,this % n_pol/))
        deallocate(temp_arr)
      else
        do ipol = 1, this % n_pol
          do iazi = 1, this % n_azi
            this % total(:,iazi,ipol) = this % absorption(:,iazi,ipol) + &
                 this % scatter(iazi,ipol) % obj % scattxs(:)
          end do
        end do
      end if

    end subroutine mgxsang_init_file

!===============================================================================
! MGXS*_PRINT displays information about a continuous-energy neutron
! cross_section table and its reactions and secondary angle/energy distributions
!===============================================================================

    subroutine mgxs_print(this, unit_)
      class(Mgxs), intent(in) :: this
      integer, intent(in)           :: unit_

      character(MAX_LINE_LEN) :: temp_str

      ! Basic nuclide information
      write(unit_,*) 'MGXS Entry: ' // trim(this % name)
      if (this % zaid > 0) then
        write(unit_,*) '  ZAID = ' // trim(to_str(this % zaid))
      else if (this % zaid < 0) then
        write(unit_,*) '  Material id = ' // trim(to_str(-this % zaid))
      end if
      if (this % awr > ZERO) then
        write(unit_,*) '  AWR = ' // trim(to_str(this % awr))
      end if
      if (this % kT > ZERO) then
        write(unit_,*) '  kT = ' // trim(to_str(this % kT))
      end if
      if (this % scatt_type == ANGLE_LEGENDRE) then
        temp_str = "Legendre"
        write(unit_,*) '  Scattering Type = ' // trim(temp_str)
        select type(this)
        type is (MgxsIso)
          temp_str = to_str(size(this % scatter % dist(1) % data,dim=1) - 1)
        end select
        write(unit_,*) '    Scattering Order = ' // trim(temp_str)
      else if (this % scatt_type == ANGLE_HISTOGRAM) then
        temp_str = "Histogram"
        write(unit_,*) '  Scattering Type = ' // trim(temp_str)
        select type(this)
        type is (MgxsIso)
          temp_str = to_str(size(this % scatter % dist(1) % data,dim=1))
        end select
        write(unit_,*) '    Num. Distribution Bins = ' // trim(temp_str)
      else if (this % scatt_type == ANGLE_TABULAR) then
        temp_str = "Tabular"
        write(unit_,*) '  Scattering Type = ' // trim(temp_str)
        select type(this)
        type is (MgxsIso)
          temp_str = to_str(size(this % scatter % dist(1) % data,dim=1))
        end select
        write(unit_,*) '    Num. Distribution Points = ' // trim(temp_str)
      end if
      write(unit_,*) '  Fissionable = ', this % fissionable

    end subroutine mgxs_print

    subroutine mgxsiso_print(this, unit)

      class(MgxsIso), intent(in) :: this
      integer, optional, intent(in)  :: unit

      integer :: unit_             ! unit to write to
      integer :: size_total, size_scattmat, size_mgxs
      integer :: gin

      ! set default unit for writing information
      if (present(unit)) then
        unit_ = unit
      else
        unit_ = OUTPUT_UNIT
      end if

      ! Write Basic Nuclide Information
      call mgxs_print(this, unit_)

      ! Determine size of mgxs and scattering matrices
      size_scattmat = 0
      do gin = 1, size(this % scatter % energy)
        size_scattmat = size_scattmat + &
             2 * size(this % scatter % energy(gin) % data) + &
             size(this % scatter % dist(gin) % data)
      end do
      size_scattmat = size_scattmat + size(this % scatter % scattxs)
      size_scattmat = size_scattmat * 8

      size_mgxs = size(this % total) + size(this % absorption) + &
           size(this % nu_fission) + size(this % k_fission) + &
           size(this % fission) + size(this % chi)
      size_mgxs = size_mgxs * 8

      ! Calculate total memory
      size_total = size_scattmat + size_mgxs

      ! Write memory used
      write(unit_,*) '  Memory Requirements'
      write(unit_,*) '    Cross sections = ' // trim(to_str(size_mgxs)) // ' bytes'
      write(unit_,*) '    Scattering Matrices = ' // &
           trim(to_str(size_scattmat)) // ' bytes'
      write(unit_,*) '    Total = ' // trim(to_str(size_total)) // ' bytes'

      ! Blank line at end of nuclide
      write(unit_,*)

    end subroutine mgxsiso_print

    subroutine mgxsang_print(this, unit)

      class(MgxsAngle), intent(in)  :: this
      integer, optional, intent(in) :: unit

      integer :: unit_             ! unit to write to
      integer :: size_total, size_scattmat, size_mgxs
      integer :: ipol, iazi, gin

      ! set default unit for writing information
      if (present(unit)) then
        unit_ = unit
      else
        unit_ = OUTPUT_UNIT
      end if

      ! Write Basic Nuclide Information
      call mgxs_print(this, unit_)

      write(unit_,*) '  # of Polar Angles = ' // trim(to_str(this % n_pol))
      write(unit_,*) '  # of Azimuthal Angles = ' // trim(to_str(this % n_azi))

      ! Determine size of mgxs and scattering matrices
      size_scattmat = 0
      do ipol = 1, this % n_pol
        do iazi = 1, this % n_azi
          do gin = 1, size(this % scatter(iazi,ipol) % obj % energy)
            size_scattmat = size_scattmat + &
                 2 * size(this % scatter(iazi,ipol) % obj % energy(gin) % data) + &
                 size(this % scatter(iazi,ipol) % obj % dist(gin) % data)
          end do
          size_scattmat = size_scattmat + &
               size(this % scatter(iazi,ipol) % obj % scattxs)
        end do
      end do
      size_scattmat = size_scattmat * 8

      size_mgxs = size(this % total) + size(this % absorption) + &
           size(this % nu_fission) + size(this % k_fission) + &
           size(this % fission) + size(this % chi)
      size_mgxs = size_mgxs * 8

      ! Calculate total memory
      size_total = size_scattmat + size_mgxs

      ! Write memory used
      write(unit_,*) '  Memory Requirements'
      write(unit_,*) '    Cross sections = ' // trim(to_str(size_mgxs)) // ' bytes'
      write(unit_,*) '    Scattering Matrices = ' // &
           trim(to_str(size_scattmat)) // ' bytes'
      write(unit_,*) '    Total = ' // trim(to_str(size_total)) // ' bytes'

      ! Blank line at end of nuclide
      write(unit_,*)
    end subroutine mgxsang_print

!===============================================================================
! MGXS*_GET_XS returns the requested data cross section data
!===============================================================================

    pure function mgxsiso_get_xs(this, xstype, gin, gout, uvw, mu) result(xs)
      class(MgxsIso), intent(in)    :: this   ! The Mgxs to initialize
      character(*) , intent(in)     :: xstype ! Type of xs requested
      integer, intent(in)           :: gin    ! Incoming Energy group
      integer, optional, intent(in) :: gout   ! Outgoing Energy group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      real(8)                       :: xs     ! Requested x/s

      select case(xstype)
      case('total')
        xs = this % total(gin)
      case('absorption')
        xs = this % absorption(gin)
      case('fission')
        if (allocated(this % fission)) then
          xs = this % fission(gin)
        else
          xs = ZERO
        end if
      case('kappa_fission')
        if (allocated(this % k_fission)) then
          xs = this % k_fission(gin)
        else
          xs = ZERO
        end if
      case('nu_fission')
        xs = this % nu_fission(gin)
      case('chi')
        if (present(gout)) then
          xs = this % chi(gout,gin)
        else
          ! Not sure youd want a 1 or a 0, but here you go!
          xs = sum(this % chi(:,gin))
        end if
      case('scatter')
        if (present(gout)) then
          if (gout < this % scatter % gmin(gin) .or. &
               gout > this % scatter % gmax(gin)) then
            xs = ZERO
          else
            xs = this % scatter % scattxs(gin) * &
                 this % scatter % energy(gin) % data(gout)
          end if
        else
          xs = this % scatter % scattxs(gin)
        end if
      case('scatter/mult')
        if (present(gout)) then
          if (gout < this % scatter % gmin(gin) .or. &
               gout > this % scatter % gmax(gin)) then
            xs = ZERO
          else
            xs = this % scatter % scattxs(gin) * &
                 this % scatter % energy(gin) % data(gout) / &
                 this % scatter % mult(gin) % data(gout)
          end if
        else
          xs = this % scatter % scattxs(gin) / &
               (dot_product(this % scatter % mult(gin) % data, &
                this % scatter % energy(gin) % data))
        end if
      case('scatter*f_mu/mult','scatter*f_mu')
        if (present(gout)) then
          if (gout < this % scatter % gmin(gin) .or. &
               gout > this % scatter % gmax(gin)) then
            xs = ZERO
          else
            xs = this % scatter % scattxs(gin) * &
                 this % scatter % energy(gin) % data(gout) * &
                 this % scatter % calc_f(gin, gout, mu)
            if (xstype == 'scatter*f_mu/mult') then
              xs = xs / this % scatter % mult(gin) % data(gout)
            end if
          end if
        else
          xs = ZERO
          ! TODO (Not likely needed)
          ! (asking for f_mu without asking for a group or mu would mean the
          ! user of this code wants the complete 1-outgoing group distribution
          ! which Im not sure what they would do with that.
        end if
      case default
        xs = ZERO
      end select

    end function mgxsiso_get_xs

    pure function mgxsang_get_xs(this, xstype, gin, gout, uvw, mu) result(xs)
      class(MgxsAngle), intent(in)  :: this   ! The Mgxs to initialize
      character(*) , intent(in)     :: xstype ! Type of xs requested
      integer, intent(in)           :: gin    ! Incoming Energy group
      integer, optional, intent(in) :: gout   ! Outgoing Energy group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      real(8)                       :: xs     ! Requested x/s

      integer :: iazi, ipol

      if (present(uvw)) then
        call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)
        select case(xstype)
        case('total')
          xs = this % total(gin,iazi,ipol)
        case('absorption')
          xs = this % absorption(gin,iazi,ipol)
        case('fission')
          if (allocated(this % fission)) then
            xs = this % fission(gin,iazi,ipol)
          else
            xs = ZERO
          end if
        case('kappa_fission')
          if (allocated(this % k_fission)) then
            xs = this % k_fission(gin,iazi,ipol)
          else
            xs = ZERO
          end if
        case('nu_fission')
          xs = this % nu_fission(gin,iazi,ipol)
        case('chi')
          if (present(gout)) then
            xs = this % chi(gout,gin,iazi,ipol)
          else
            ! Not sure youd want a 1 or a 0, but here you go!
            xs = sum(this % chi(:,gin,iazi,ipol))
          end if
        case('scatter')
          if (present(gout)) then
            if (gout < this % scatter(iazi,ipol) % obj % gmin(gin) .or. &
                 gout > this % scatter(iazi,ipol) % obj % gmax(gin)) then
              xs = ZERO
            else
              xs = this % scatter(iazi,ipol) % obj % scattxs(gin) * &
                   this % scatter(iazi,ipol) % obj % energy(gin) % data(gout)
            end if
          else
            xs = this % scatter(iazi,ipol) % obj % scattxs(gin)
          end if
        case('scatter/mult')
          if (present(gout)) then
            if (gout < this % scatter(iazi,ipol) % obj % gmin(gin) .or. &
                 gout > this % scatter(iazi,ipol) % obj % gmax(gin)) then
              xs = ZERO
            else
              xs = this % scatter(iazi,ipol) % obj % scattxs(gin) * &
                   this % scatter(iazi,ipol) % obj % energy(gin) % data(gout) / &
                   this % scatter(iazi,ipol) % obj % mult(gin) % data(gout)
            end if
          else
            xs = this % scatter(iazi,ipol) % obj % scattxs(gin) / &
                 (dot_product(this % scatter(iazi,ipol) % obj % mult(gin) % data, &
                  this % scatter(iazi,ipol) % obj % energy(gin) % data))
          end if
        case('scatter*f_mu/mult','scatter*f_mu')
          if (present(gout)) then
            if (gout < this % scatter(iazi,ipol) % obj % gmin(gin) .or. &
                 gout > this % scatter(iazi,ipol) % obj % gmax(gin)) then
              xs = ZERO
            else
              xs = this % scatter(iazi,ipol) % obj % scattxs(gin) * &
                   this % scatter(iazi,ipol) % obj % energy(gin) % data(gout)
              xs = xs * this % scatter(iazi,ipol) % obj % calc_f(gin, gout, mu)
              if (xstype == 'scatter*f_mu/mult') then
                xs = xs / this % scatter(iazi,ipol) % obj % mult(gin) % data(gout)
              end if
            end if
          else
            xs = ZERO
            ! TODO (Not likely needed)
            ! (asking for f_mu without asking for a group or mu would mean the
            ! user of this code wants the complete 1-outgoing group distribution
            ! which Im not sure what they would do with that.
          end if
        case default
          xs = ZERO
        end select
      else
        xs = ZERO
      end if

    end function mgxsang_get_xs

!===============================================================================
! MACROXS*_COMBINE Builds a macroscopic Mgxs object from microscopic Mgxs
! objects
!===============================================================================

    subroutine mgxs_combine(this,mat,scatt_type,i_listing)
      class(Mgxs), intent(inout)          :: this ! The Mgxs to initialize
      type(Material), pointer, intent(in) :: mat  ! base material
      integer, intent(in)                 :: scatt_type ! How is data presented
      integer, intent(in)                 :: i_listing  ! Index in listings

      ! Fill in meta-data from material information
      if (mat % name == "") then
        this % name      = trim(to_str(mat % id))
      else
        this % name      = mat % name
      end if
      this % zaid        = -mat % id
      this % listing     = i_listing
      this % fissionable = mat % fissionable
      this % scatt_type  = scatt_type

      ! The following info we should initialize, but we dont need it nor
      ! does it have guaranteed meaning.
      this % awr = -ONE
      this % kT  = -ONE

    end subroutine mgxs_combine

    subroutine mgxsiso_combine(this,mat,nuclides,groups,max_order,scatt_type, &
                               i_listing)
      class(MgxsIso), intent(inout)       :: this ! The Mgxs to initialize
      type(Material), pointer, intent(in) :: mat  ! base material
      type(MgxsContainer), intent(in)     :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                 :: groups     ! Number of E groups
      integer, intent(in)                 :: max_order  ! Maximum requested order
      integer, intent(in)                 :: scatt_type ! How is data presented
      integer, intent(in)                 :: i_listing  ! Index in listings

      integer :: i             ! loop index over nuclides
      integer :: gin, gout     ! group indices
      real(8) :: atom_density  ! atom density of a nuclide
      real(8) :: norm
      integer :: mat_max_order, order, order_dim, nuc_order_dim
      real(8), allocatable :: temp_mult(:,:)
      real(8), allocatable :: scatt_coeffs(:,:,:)

      ! Set the meta-data
      call mgxs_combine(this,mat,scatt_type,i_listing)

      ! Determine the scattering type of our data and ensure all scattering orders
      ! are the same.
      select type(nuc => nuclides(mat % nuclide(1)) % obj)
      type is (MgxsIso)
        order = size(nuc % scatter % dist(1) % data, dim=1)
      end select
      ! If we have tabular only data, then make sure all datasets have same size
      if (scatt_type == ANGLE_HISTOGRAM) then
        ! Check all scattering data to ensure it is the same size
        ! order = size(nuclides(mat % nuclide(1)) % obj % scatter % data,dim=1)
        do i = 2, mat % n_nuclides
          select type(nuc => nuclides(mat % nuclide(i)) % obj)
          type is (MgxsIso)
            if (order /= size(nuc % scatter % dist(1) % data,dim=1)) &
                 call fatal_error("All Histogram Scattering Entries Must Be&
                                  & Same Length!")
          end select
        end do
        ! Ok, got our order, store the dimensionality
        order_dim = order

        ! Set our Scatter Object Type
        allocate(ScattDataHistogram :: this % scatter)

      else if (scatt_type == ANGLE_TABULAR) then
        ! Check all scattering data to ensure it is the same size
        do i = 2, mat % n_nuclides
          select type(nuc => nuclides(mat % nuclide(i)) % obj)
          type is (MgxsIso)
            if (order /= size(nuc % scatter % dist(1) % data,dim=1)) &
                 call fatal_error("All Tabular Scattering Entries Must Be&
                                  & Same Length!")
          end select
        end do
        ! Ok, got our order, store the dimensionality
        order_dim = order

        ! Set our Scatter Object Type
        allocate(ScattDataTabular :: this % scatter)

      else if (scatt_type == ANGLE_LEGENDRE) then
        ! Need to determine the maximum scattering order of all data in this material
        mat_max_order = 0
        do i = 1, mat % n_nuclides
          select type(nuc => nuclides(mat % nuclide(i)) % obj)
          type is (MgxsIso)
            if (size(nuc % scatter % dist(1) % data,dim=1) > mat_max_order) &
                 mat_max_order = size(nuc % scatter % dist(1) % data,dim=1)
          end select
        end do

        ! Now need to compare this material maximum scattering order with
        ! the problem wide max scatt order and use whichever is lower
        order = min(mat_max_order, max_order)
        ! Ok, got our order, store the dimensionality
        order_dim = order + 1

        ! Set our Scatter Object Type
        allocate(ScattDataLegendre :: this % scatter)
      end if

      ! Allocate and initialize data needed for macro_xs(i_mat) object
      allocate(this % total(groups))
      this % total = ZERO
      allocate(this % absorption(groups))
      this % absorption = ZERO
      allocate(this % fission(groups))
      this % fission = ZERO
      allocate(this % k_fission(groups))
      this % k_fission = ZERO
      allocate(this % nu_fission(groups))
      this % nu_fission = ZERO
      allocate(this % chi(groups,groups))
      this % chi = ZERO
      allocate(temp_mult(groups,groups))
      temp_mult = ZERO
      allocate(scatt_coeffs(order_dim,groups,groups))
      scatt_coeffs = ZERO

      ! Add contribution from each nuclide in material
      do i = 1, mat % n_nuclides
        ! Copy atom density of nuclide in material
        atom_density = mat % atom_density(i)

        ! Perform our operations which depend upon the type
        select type(nuc => nuclides(mat % nuclide(i)) % obj)
        type is (MgxsIso)
          ! Add contributions to total, absorption, and fission data (if necessary)
          this % total = this % total + atom_density * nuc % total
          this % absorption = this % absorption + &
               atom_density * nuc % absorption
          if (nuc % fissionable) then
            this % chi = this % chi + atom_density * nuc % chi
            this % nu_fission = this % nu_fission + atom_density * &
                 nuc % nu_fission
            if (allocated(nuc % fission)) then
              this % fission = this % fission + atom_density * nuc % fission
            end if
            if (allocated(nuc % k_fission)) then
              this % k_fission = this % k_fission + atom_density * nuc % k_fission
            end if
          end if

          ! Get the multiplication matrix
          do gin = 1, groups
            do gout = nuc % scatter % gmin(gin), nuc % scatter % gmax(gin)
              temp_mult(gout,gin) = temp_mult(gout,gin) + atom_density * &
                   nuc % scatter % mult(gin) % data(gout)
            end do
          end do

          ! Get the complete scattering matrix
          nuc_order_dim = size(nuc % scatter % dist(1) % data,dim=1)
          scatt_coeffs(1:min(nuc_order_dim, order_dim),:,:) = &
               scatt_coeffs(1:min(nuc_order_dim, order_dim),:,:) + &
               atom_density * &
               nuc % scatter % get_matrix(min(nuc_order_dim,order_dim))

        type is (MgxsAngle)
          call fatal_error("Invalid Passing of MgxsAngle to MgxsIso Object")
        end select
      end do

      ! Initialize the ScattData Object
      call this % scatter % init(temp_mult,scatt_coeffs)

      ! Now normalize chi
      if (mat % fissionable) then
        do gin = 1, groups
          norm =  sum(this % chi(:,gin))
          if (norm > ZERO) then
            this % chi(:,gin) = this % chi(:,gin) / norm
          end if
        end do
      end if

      ! Deallocate temporaries
      deallocate(scatt_coeffs, temp_mult)

    end subroutine mgxsiso_combine

    subroutine mgxsang_combine(this,mat,nuclides,groups,max_order,scatt_type,&
                               i_listing)
      class(MgxsAngle), intent(inout)     :: this ! The Mgxs to initialize
      type(Material), pointer, intent(in) :: mat  ! base material
      type(MgxsContainer), intent(in)     :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                 :: groups     ! Number of E groups
      integer, intent(in)                 :: max_order  ! Maximum requested order
      integer, intent(in)                 :: scatt_type ! Legendre or Tabular Scatt?
      integer, intent(in)                 :: i_listing  ! Index in listings

      integer :: i             ! loop index over nuclides
      integer :: gin, gout     ! group indices
      real(8) :: atom_density  ! atom density of a nuclide
      integer :: ipol, iazi, n_pol, n_azi
      real(8) :: norm
      integer :: mat_max_order, order, order_dim, nuc_order_dim
      real(8), allocatable :: temp_mult(:,:,:,:)
      real(8), allocatable :: scatt_coeffs(:,:,:,:,:)

      ! Set the meta-data
      call mgxs_combine(this,mat,scatt_type,i_listing)

      ! Get the number of each polar and azi angles and make sure all the
      ! NuclideAngle types have the same number of these angles
      n_pol = -1
      n_azi = -1
      do i = 1, mat % n_nuclides
        select type(nuc => nuclides(mat % nuclide(i)) % obj)
        type is (MgxsAngle)
          if (n_pol == -1) then
            n_pol = nuc % n_pol
            n_azi = nuc % n_azi
            allocate(this % polar(n_pol))
            this % polar = nuc % polar
            allocate(this % azimuthal(n_azi))
            this % azimuthal = nuc % azimuthal
          else
            if ((n_pol /= nuc % n_pol) .or. (n_azi /= nuc % n_azi)) then
              call fatal_error("All Angular Data Must Be Same Length!")
            end if
          end if
        end select
      end do

      ! Determine the scattering type of our data and ensure all scattering orders
      ! are the same.
      select type(nuc => nuclides(mat % nuclide(1)) % obj)
      type is (MgxsAngle)
        order = size(nuc % scatter(1,1) % obj % dist(1) % data, dim=1)
      end select
      ! If we have tabular only data, then make sure all datasets have same size
      if (scatt_type == ANGLE_HISTOGRAM) then
        ! Check all scattering data to ensure it is the same size
        ! order = size(nuclides(mat % nuclide(1)) % obj % scatter % data,dim=1)
        do i = 2, mat % n_nuclides
          select type(nuc => nuclides(mat % nuclide(i)) % obj)
          type is (MgxsAngle)
            if (order /= size(nuc % scatter(1,1) % obj % dist(1) % data,dim=1)) &
                 call fatal_error("All Histogram Scattering Entries Must Be&
                                  & Same Length!")
          end select
        end do
        ! Ok, got our order, store the dimensionality
        order_dim = order

        ! Set our Scatter Object Type
        allocate(this % scatter(n_azi, n_pol))
        do ipol = 1, n_pol
          do iazi = 1, n_azi
            allocate(ScattDataHistogram :: this % scatter(iazi, ipol) % obj)
          end do
        end do

      else if (scatt_type == ANGLE_TABULAR) then
        ! Check all scattering data to ensure it is the same size
        do i = 2, mat % n_nuclides
          select type(nuc => nuclides(mat % nuclide(i)) % obj)
          type is (MgxsAngle)
            if (order /= size(nuc % scatter(1,1) % obj % dist(1) % data,dim=1)) &
                 call fatal_error("All Tabular Scattering Entries Must Be&
                                  & Same Length!")
          end select
        end do
        ! Ok, got our order, store the dimensionality
        order_dim = order

        ! Set our Scatter Object Type
        allocate(this % scatter(n_azi, n_pol))
        do ipol = 1, n_pol
          do iazi = 1, n_azi
            allocate(ScattDataTabular :: this % scatter(iazi, ipol) % obj)
          end do
        end do

      else if (scatt_type == ANGLE_LEGENDRE) then
        ! Need to determine the maximum scattering order of all data in this material
        mat_max_order = 0
        do i = 1, mat % n_nuclides
          select type(nuc => nuclides(mat % nuclide(i)) % obj)
          type is (MgxsAngle)
            if (size(nuc % scatter(1,1) % obj % dist(1) % data,dim=1) > mat_max_order) &
                 mat_max_order = size(nuc % scatter(1,1) % obj% dist(1) % data,dim=1)
          end select
        end do

        ! Now need to compare this material maximum scattering order with
        ! the problem wide max scatt order and use whichever is lower
        order = min(mat_max_order, max_order)
        ! Ok, got our order, store the dimensionality
        order_dim = order + 1

        ! Set our Scatter Object Type
        allocate(this % scatter(n_azi, n_pol))
        do ipol = 1, n_pol
          do iazi = 1, n_azi
            allocate(ScattDataLegendre :: this % scatter(iazi, ipol) % obj)
          end do
        end do
      end if

      ! Allocate and initialize data within macro_xs(i_mat) object
      allocate(this % total(groups,n_azi,n_pol))
      this % total = ZERO
      allocate(this % absorption(groups,n_azi,n_pol))
      this % absorption = ZERO
      allocate(this % fission(groups,n_azi,n_pol))
      this % fission = ZERO
      allocate(this % k_fission(groups,n_azi,n_pol))
      this % k_fission = ZERO
      allocate(this % nu_fission(groups,n_azi,n_pol))
      this % nu_fission = ZERO
      allocate(this % chi(groups,groups,n_azi,n_pol))
      this % chi = ZERO
      allocate(temp_mult(groups,groups,n_azi,n_pol))
      temp_mult = ZERO
      allocate(scatt_coeffs(order_dim,groups,groups,n_azi,n_pol))
      scatt_coeffs = ZERO

      ! Add contribution from each nuclide in material
      do i = 1, mat % n_nuclides
        ! Copy atom density of nuclide in material
        atom_density = mat % atom_density(i)

        ! Perform our operations which depend upon the type
        select type(nuc => nuclides(mat % nuclide(i)) % obj)
        type is (MgxsIso)
          call fatal_error("Invalid Passing of MgxsIso to MgxsAngle Object")
        type is (MgxsAngle)
          ! Add contributions to total, absorption, and fission data (if necessary)
          this % total = this % total + atom_density * nuc % total
          this % absorption = this % absorption + &
               atom_density * nuc % absorption
          if (nuc % fissionable) then
            this % chi = this % chi + atom_density * nuc % chi
            this % nu_fission = this % nu_fission + atom_density * &
                 nuc % nu_fission
            if (allocated(nuc % fission)) then
              this % fission = this % fission + atom_density * nuc % fission
            end if
            if (allocated(nuc % k_fission)) then
              this % k_fission = this % k_fission + atom_density * nuc % k_fission
            end if
          end if

          ! Get the multiplication matrix
          do ipol = 1, n_pol
            do iazi = 1, n_azi
              do gin = 1, groups
                do gout = nuc % scatter(iazi,ipol) % obj % gmin(gin), &
                     nuc % scatter(iazi,ipol) % obj % gmax(gin)
                  temp_mult(gout,gin,iazi,ipol) = temp_mult(gout,gin,iazi,ipol) + &
                       atom_density * &
                       nuc % scatter(iazi,ipol) % obj % mult(gin) % data(gout)
                end do
              end do
            end do
          end do

          ! Get the complete scattering matrix
          nuc_order_dim = size(nuc % scatter(1,1) % obj % dist(1) % data,dim=1)
          do ipol = 1, n_pol
            do iazi = 1, n_azi
              scatt_coeffs(1:min(nuc_order_dim,order_dim),:,:,iazi,ipol) = &
                   scatt_coeffs(1:min(nuc_order_dim, order_dim),:,:,iazi,ipol) + &
                   atom_density * &
                   nuc % scatter(iazi,ipol) % obj % get_matrix(&
                   min(nuc_order_dim,order_dim))
            end do
          end do
        end select
      end do

      ! Initialize the ScattData Object
      do ipol = 1, n_pol
        do iazi = 1, n_azi
          call this % scatter(iazi,ipol) % obj % init( &
               temp_mult(:,:,iazi,ipol), scatt_coeffs(:,:,:,iazi,ipol))
        end do
      end do

      ! Now normalize chi
      if (mat % fissionable) then
        do ipol = 1, n_pol
          do iazi = 1, n_azi
            do gin = 1, groups
              norm =  sum(this % chi(:,gin,iazi,ipol))
              if (norm > ZERO) then
                this % chi(:,gin,iazi,ipol) = this % chi(:,gin,iazi,ipol) / norm
              end if
            end do
          end do
        end do
      end if

      ! Deallocate temporaries for the next material
      deallocate(scatt_coeffs, temp_mult)

    end subroutine mgxsang_combine

!===============================================================================
! MGXS*_SAMPLE_FISSION_ENERGY samples the outgoing energy from a fission event
!===============================================================================

    function mgxsiso_sample_fission_energy(this, gin, uvw) result(gout)
      class(MgxsIso), intent(in) :: this   ! Data to work with
      integer, intent(in)           :: gin    ! Incoming energy group
      real(8), intent(in)           :: uvw(3) ! Particle Direction
      integer                       :: gout   ! Sampled outgoing group
      real(8) :: xi               ! Our random number
      real(8) :: prob             ! Running probability

      xi = prn()
      gout = 1
      prob = this % chi(gout,gin)

      do while (prob < xi)
        gout = gout + 1
        prob = prob + this % chi(gout,gin)
      end do

    end function mgxsiso_sample_fission_energy

    function mgxsang_sample_fission_energy(this, gin, uvw) result(gout)
      class(MgxsAngle), intent(in) :: this  ! Data to work with
      integer, intent(in)             :: gin    ! Incoming energy group
      real(8), intent(in)             :: uvw(3) ! Particle Direction
      integer                         :: gout   ! Sampled outgoing group
      real(8) :: xi               ! Our random number
      real(8) :: prob             ! Running probability
      integer :: iazi, ipol

      call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)

      xi = prn()
      gout = 1
      prob = this % chi(gout,gin,iazi,ipol)

      do while (prob < xi)
        gout = gout + 1
        prob = prob + this % chi(gout,gin,iazi,ipol)
      end do

    end function mgxsang_sample_fission_energy

!===============================================================================
! MGXS*_SAMPLE_SCATTER Selects outgoing energy and angle after a scatter event
!===============================================================================

    subroutine mgxsiso_sample_scatter(this, uvw, gin, gout, mu, wgt)
      class(MgxsIso), intent(in)    :: this
      real(8),           intent(in)    :: uvw(3) ! Incoming neutron direction
      integer,           intent(in)    :: gin    ! Incoming neutron group
      integer,           intent(out)   :: gout   ! Sampled outgoin group
      real(8),           intent(out)   :: mu     ! Sampled change in angle
      real(8),           intent(inout) :: wgt    ! Particle weight

      call this % scatter % sample(gin, gout, mu, wgt)

    end subroutine mgxsiso_sample_scatter

    subroutine mgxsang_sample_scatter(this, uvw, gin, gout, mu, wgt)
      class(MgxsAngle), intent(in)    :: this
      real(8),             intent(in)    :: uvw(3) ! Incoming neutron direction
      integer,             intent(in)    :: gin    ! Incoming neutron group
      integer,             intent(out)   :: gout   ! Sampled outgoin group
      real(8),             intent(out)   :: mu     ! Sampled change in angle
      real(8),             intent(inout) :: wgt    ! Particle weight

      integer :: iazi, ipol ! Angular indices

      call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)
      call this % scatter(iazi,ipol) % obj % sample(gin,gout,mu,wgt)

    end subroutine mgxsang_sample_scatter

!===============================================================================
! MGXS*_CALCULATE_XS determines the multi-group cross sections
! for the material the particle is currently traveling through.
!===============================================================================

    subroutine mgxsiso_calculate_xs(this, gin, uvw, xs)
      class(MgxsIso),        intent(in)    :: this
      integer,               intent(in)    :: gin    ! Incoming neutron group
      real(8),               intent(in)    :: uvw(3) ! Incoming neutron direction
      type(MaterialMacroXS), intent(inout) :: xs     ! Resultant Mgxs Data

      xs % total         = this % total(gin)
      xs % elastic       = this % scatter % scattxs(gin)
      xs % absorption    = this % absorption(gin)
      xs % fission       = this % fission(gin)
      xs % nu_fission    = this % nu_fission(gin)

    end subroutine mgxsiso_calculate_xs

    subroutine mgxsang_calculate_xs(this, gin, uvw, xs)
      class(MgxsAngle),      intent(in)    :: this
      integer,               intent(in)    :: gin    ! Incoming neutron group
      real(8),               intent(in)    :: uvw(3) ! Incoming neutron direction
      type(MaterialMacroXS), intent(inout) :: xs     ! Resultant Mgxs Data

      integer :: iazi, ipol

      call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)
      xs % total         = this % total(gin,iazi,ipol)
      xs % elastic       = this % scatter(iazi,ipol) % obj % scattxs(gin)
      xs % absorption    = this % absorption(gin,iazi,ipol)
      xs % fission       = this % fission(gin,iazi,ipol)
      xs % nu_fission    = this % nu_fission(gin,iazi,ipol)

    end subroutine mgxsang_calculate_xs

!!!TODO:
! Move find_angle from math to here after we fully implement this and are ready
! to delete macroxs_header and relevant portions from nuclide_header.
end module mgxs_header