module nuclide_header

  use, intrinsic :: ISO_FORTRAN_ENV

  use ace_header
  use constants
  use endf,        only: reaction_name
  use error,       only: fatal_error
  use list_header, only: ListInt
  use math,        only: evaluate_legendre, find_angle
  use scattdata_header
  use string
  use xml_interface

  implicit none

!===============================================================================
! Nuclide contains the base nuclidic data for a nuclide, which does not depend
! upon how the nuclear data is represented (i.e., CE, or any variant of MG).
! The extended types, NuclideCE and NuclideMG deal with the rest
!===============================================================================

  type, abstract :: Nuclide
    character(12) :: name    ! name of nuclide, e.g. 92235.03c
    integer       :: zaid    ! Z and A identifier, e.g. 92235
    real(8)       :: awr     ! Atomic Weight Ratio
    integer       :: listing ! index in xs_listings
    real(8)       :: kT      ! temperature in MeV (k*T)

    ! Linked list of indices in nuclides array of instances of this same nuclide
    type(VectorInt) :: nuc_list

    ! Fission information
    logical :: fissionable         ! nuclide is fissionable?

  contains
    procedure(nuclide_print_), deferred :: print ! Writes nuclide info
  end type Nuclide

  abstract interface
    subroutine nuclide_print_(this, unit)
      import Nuclide
      class(Nuclide),intent(in)     :: this
      integer, optional, intent(in) :: unit
    end subroutine nuclide_print_
  end interface

  type, extends(Nuclide) :: NuclideCE
    ! Energy grid information
    integer :: n_grid                     ! # of nuclide grid points
    integer, allocatable :: grid_index(:) ! log grid mapping indices
    real(8), allocatable :: energy(:)     ! energy values corresponding to xs

    ! Microscopic cross sections
    real(8), allocatable :: total(:)      ! total cross section
    real(8), allocatable :: elastic(:)    ! elastic scattering
    real(8), allocatable :: fission(:)    ! fission
    real(8), allocatable :: nu_fission(:) ! neutron production
    real(8), allocatable :: absorption(:) ! absorption (MT > 100)
    real(8), allocatable :: heating(:)    ! heating

    ! Resonance scattering info
    logical              :: resonant = .false. ! resonant scatterer?
    character(10)        :: name_0K = '' ! name of 0K nuclide, e.g. 92235.00c
    character(16)        :: scheme ! target velocity sampling scheme
    integer              :: n_grid_0K ! number of 0K energy grid points
    real(8), allocatable :: energy_0K(:)  ! energy grid for 0K xs
    real(8), allocatable :: elastic_0K(:) ! Microscopic elastic cross section
    real(8), allocatable :: xs_cdf(:) ! CDF of v_rel times cross section
    real(8)              :: E_min ! lower cutoff energy for res scattering
    real(8)              :: E_max ! upper cutoff energy for res scattering

    ! Fission information
    logical :: has_partial_fission ! nuclide has partial fission reactions?
    integer :: n_fission           ! # of fission reactions
    integer, allocatable :: index_fission(:) ! indices in reactions

    ! Total fission neutron emission
    integer :: nu_t_type
    real(8), allocatable :: nu_t_data(:)

    ! Prompt fission neutron emission
    integer :: nu_p_type
    real(8), allocatable :: nu_p_data(:)

    ! Delayed fission neutron emission
    integer :: nu_d_type
    integer :: n_precursor ! # of delayed neutron precursors
    real(8), allocatable :: nu_d_data(:)
    real(8), allocatable :: nu_d_precursor_data(:)
    type(AngleEnergyContainer), allocatable :: nu_d_edist(:)

    ! Unresolved resonance data
    logical                :: urr_present
    integer                :: urr_inelastic
    type(UrrData), pointer :: urr_data => null()

    ! Reactions
    integer :: n_reaction ! # of reactions
    type(Reaction), allocatable :: reactions(:)
    type(DictIntInt) :: reaction_index ! map MT values to index in reactions
                                       ! array; used at tally-time

  contains
    procedure :: clear => nuclidece_clear
    procedure :: print => nuclidece_print
  end type NuclideCE

  type, abstract, extends(Nuclide) :: NuclideMG
    integer :: scatt_type ! either legendre, histogram, or tabular.
  contains
    procedure(nuclidemg_init_),   deferred :: init   ! Initialize the data
    procedure(nuclidemg_get_xs_), deferred :: get_xs ! Get the requested xs
  end type NuclideMG

  abstract interface

    subroutine nuclidemg_init_(this, node_xsdata, groups, get_kfiss, get_fiss, &
                               max_order)
      import NuclideMG, Node
      class(NuclideMG), intent(inout) :: this        ! Working Object
      type(Node), pointer, intent(in) :: node_xsdata ! Data from MGXS xml
      integer, intent(in)             :: groups      ! Number of Energy groups
      logical, intent(in)             :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)             :: get_fiss    ! Should we get fiss data?
      integer, intent(in)             :: max_order ! Maximum requested order
    end subroutine nuclidemg_init_

    function nuclidemg_get_xs_(this, g, xstype, gout, uvw, mu, iazi, ipol) &
         result(xs)
      import NuclideMG
      class(NuclideMG), intent(in) :: this
      integer, intent(in)           :: g      ! Incoming Energy group
      character(*), intent(in)      :: xstype ! Cross Section Type
      integer, optional, intent(in) :: gout   ! Outgoing Group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      integer, optional, intent(in) :: iazi  ! Azimuthal Index
      integer, optional, intent(in) :: ipol  ! Polar Index
      real(8)                       :: xs     ! Resultant xs
    end function nuclidemg_get_xs_

    pure function nuclidemg_calc_f_(this, gin, gout, mu, uvw, iazi, ipol) result(f)
      import NuclideMG
      class(NuclideMG), intent(in) :: this
      integer, intent(in)           :: gin   ! Incoming Energy Group
      integer, intent(in)           :: gout  ! Outgoing Energy Group
      real(8), intent(in)           :: mu    ! Angle of interest
      real(8), intent(in), optional :: uvw(3) ! Direction vector
      integer, intent(in), optional :: iazi ! Incoming Energy Group
      integer, intent(in), optional :: ipol ! Outgoing Energy Group
      real(8)                       :: f     ! Return value of f(mu)

    end function nuclidemg_calc_f_
  end interface

!===============================================================================
! NuclideIso contains the base MGXS data for a nuclide specifically for
! isotropically weighted MGXS
!===============================================================================

  type, extends(NuclideMG) :: NuclideIso

    ! Microscopic cross sections
    real(8), allocatable :: total(:)        ! total cross section
    real(8), allocatable :: absorption(:)   ! absorption cross section
    class(ScattData), allocatable :: scatter ! scattering information
    real(8), allocatable :: nu_fission(:,:) ! fission matrix (Gout x Gin)
    real(8), allocatable :: k_fission(:)    ! kappa-fission
    real(8), allocatable :: fission(:)      ! neutron production
    real(8), allocatable :: chi(:)          ! Fission Spectra

  contains
    procedure :: init   => nuclideiso_init   ! Initialize Nuclidic MGXS Data
    procedure :: print  => nuclideiso_print  ! Writes nuclide info
    procedure :: get_xs => nuclideiso_get_xs ! Gets Size of Data w/in Object
  end type NuclideIso

!===============================================================================
! NuclideAngle contains the base MGXS data for a nuclide specifically for
! explicit angle-dependent weighted MGXS
!===============================================================================

  type, extends(NuclideMG) :: NuclideAngle

    ! Microscopic cross sections. Dimensions are: (n_pol, n_azi, Nl, Ng, Ng)
    real(8), allocatable :: total(:,:,:)        ! total cross section
    real(8), allocatable :: absorption(:,:,:)   ! absorption cross section
    type(ScattDataContainer), allocatable :: scatter(:,:) ! scattering information
    real(8), allocatable :: nu_fission(:,:,:,:) ! fission matrix (Gout x Gin)
    real(8), allocatable :: k_fission(:,:,:)    ! kappa-fission
    real(8), allocatable :: fission(:,:,:)      ! neutron production
    real(8), allocatable :: chi(:,:,:)          ! Fission Spectra
    real(8), allocatable :: mult(:,:,:,:)       ! Scatter multiplicity (Gout x Gin)

    ! In all cases, right-most indices are theta, phi
    integer              :: n_pol         ! Number of polar angles
    integer              :: n_azi         ! Number of azimuthal angles
    real(8), allocatable :: polar(:)     ! polar angles
    real(8), allocatable :: azimuthal(:) ! azimuthal angles

  contains
    procedure :: init   => nuclideangle_init   ! Initialize Nuclidic MGXS Data
    procedure :: print  => nuclideangle_print  ! Gets Size of Data w/in Object
    procedure :: get_xs => nuclideangle_get_xs ! Gets Size of Data w/in Object
  end type NuclideAngle

!===============================================================================
! NUCLIDEMGCONTAINER pointer array for storing Nuclides
!===============================================================================

  type NuclideMGContainer
    class(NuclideMG), pointer :: obj
  end type NuclideMGContainer

!===============================================================================
! NUCLIDE0K temporarily contains all 0K cross section data and other parameters
! needed to treat resonance scattering before transferring them to NuclideCE
!===============================================================================

  type Nuclide0K
    character(10) :: nuclide             ! name of nuclide, e.g. U-238
    character(16) :: scheme = 'ares'     ! target velocity sampling scheme
    character(10) :: name                ! name of nuclide, e.g. 92235.03c
    character(10) :: name_0K             ! name of 0K nuclide, e.g. 92235.00c
    real(8)       :: E_min = 0.01e-6_8   ! lower cutoff energy for res scattering
    real(8)       :: E_max = 1000.0e-6_8 ! upper cutoff energy for res scattering
  end type Nuclide0K

!===============================================================================
! NUCLIDEMICROXS contains cached microscopic cross sections for a
! particular nuclide at the current energy
!===============================================================================

  type NuclideMicroXS
    integer :: index_grid      ! index on nuclide energy grid
    integer :: index_temp      ! temperature index for nuclide
    real(8) :: last_E = ZERO   ! last evaluated energy
    real(8) :: interp_factor   ! interpolation factor on nuc. energy grid
    real(8) :: total           ! microscropic total xs
    real(8) :: elastic         ! microscopic elastic scattering xs
    real(8) :: absorption      ! microscopic absorption xs
    real(8) :: fission         ! microscopic fission xs
    real(8) :: nu_fission      ! microscopic production xs

    ! Information for S(a,b) use
    integer :: index_sab          ! index in sab_tables (zero means no table)
    integer :: last_index_sab = 0 ! index in sab_tables last used by this nuclide
    real(8) :: elastic_sab        ! microscopic elastic scattering on S(a,b) table

    ! Information for URR probability table use
    logical :: use_ptable  ! in URR range with probability tables?
    real(8) :: last_prn
  end type NuclideMicroXS

!===============================================================================
! MATERIALMACROXS contains cached macroscopic cross sections for the material a
! particle is traveling through
!===============================================================================

  type MaterialMacroXS
    real(8) :: total         ! macroscopic total xs
    real(8) :: elastic       ! macroscopic elastic scattering xs
    real(8) :: absorption    ! macroscopic absorption xs
    real(8) :: fission       ! macroscopic fission xs
    real(8) :: nu_fission    ! macroscopic production xs
  end type MaterialMacroXS

!===============================================================================
! XSLISTING contains data read from a CE or MG cross_sections.xml file
! (or equivalent)
!===============================================================================

  type XsListing
    character(12) :: name       ! table name, e.g. 92235.70c
    character(12) :: alias      ! table alias, e.g. U-235.70c
    integer       :: type       ! type of table (cont-E neutron, S(A,b), etc)
    integer       :: zaid       ! ZAID identifier = 1000*Z + A
    integer       :: filetype   ! ASCII or BINARY
    integer       :: location   ! location of table within library
    integer       :: recl       ! record length for library
    integer       :: entries    ! number of entries per record
    real(8)       :: awr        ! atomic weight ratio (# of neutron masses)
    real(8)       :: kT         ! Boltzmann constant * temperature (MeV)
    logical       :: metastable ! is this nuclide metastable?
    character(MAX_FILE_LEN) :: path ! path to library containing table
  end type XsListing

  contains

!===============================================================================
! NUCLIDE_*_INIT reads in the data from the XML file, as already accessed
!===============================================================================

    subroutine nuclidemg_init(this, node_xsdata)
      class(NuclideMG), intent(inout) :: this        ! Working Object
      type(Node), pointer, intent(in) :: node_xsdata ! Data from MGXS xml

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
        this % zaid = -1
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

    end subroutine nuclidemg_init

    subroutine nuclideiso_init(this, node_xsdata, groups, get_kfiss, get_fiss, &
                               max_order)
      class(NuclideIso), intent(inout) :: this        ! Working Object
      type(Node), pointer, intent(in)  :: node_xsdata ! Data from MGXS xml
      integer, intent(in)              :: groups      ! Number of Energy groups
      logical, intent(in)              :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)              :: get_fiss    ! Need fiss data?
      integer, intent(in)              :: max_order   ! Maximum requested order

      type(Node), pointer     :: node_legendre_mu
      character(MAX_LINE_LEN) :: temp_str
      logical                 :: enable_leg_mu
      real(8), allocatable    :: temp_arr(:)
      real(8), allocatable    :: temp_mult(:,:)
      real(8), allocatable    :: scatt_coeffs(:,:,:)
      real(8), allocatable    :: input_scatt(:,:,:)
      real(8), allocatable    :: temp_scatt(:,:,:)
      real(8)                 :: dmu, mu, norm
      integer                 :: order, order_dim, gin, gout, l, arr_len
      integer                 :: legendre_mu_points, imu

      ! Call generic data gathering routine (will populate the metadata)
      call nuclidemg_init(this, node_xsdata)

      ! Load the more specific data
      if (this % fissionable) then

        if (check_for_node(node_xsdata,"chi")) then
          ! Get chi
          allocate(this % chi(groups))
          call get_node_array(node_xsdata,"chi",this % chi)

          ! Get nu_fission (as a vector)
          if (check_for_node(node_xsdata,"nu_fission")) then
            allocate(temp_arr(groups * 1))
            call get_node_array(node_xsdata,"nu_fission",temp_arr)
            allocate(this % nu_fission(groups,1))
            this % nu_fission = reshape(temp_arr,(/groups,1/))
            deallocate(temp_arr)
          else
            call fatal_error("If fissionable, must provide nu_fission!")
          end if

        else
          ! Get nu_fission (as a matrix)
          if (check_for_node(node_xsdata,"nu_fission")) then

            allocate(temp_arr(groups*groups))
            call get_node_array(node_xsdata,"nu_fission",temp_arr)
            allocate(this % nu_fission(groups, groups))
            this % nu_fission = reshape(temp_arr,(/groups,groups/))
            deallocate(temp_arr)
          else
            call fatal_error("If fissionable, must provide nu_fission!")
          end if
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

      ! Get, or infer, total xs data.
      allocate(this % total(groups))
      if (check_for_node(node_xsdata,"total")) then
        call get_node_array(node_xsdata,"total",this % total)
      else
        this % total = this % absorption + this % scatter % scattxs
      end if

      ! Deallocate temporaries for the next material
      deallocate(input_scatt,scatt_coeffs,temp_mult)

    end subroutine nuclideiso_init

    subroutine nuclideangle_init(this, node_xsdata, groups, get_kfiss, get_fiss, &
                                 max_order)
      class(NuclideAngle), intent(inout) :: this        ! Working Object
      type(Node), pointer, intent(in)    :: node_xsdata ! Data from MGXS xml
      integer, intent(in)                :: groups      ! Number of Energy groups
      logical, intent(in)                :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)                :: get_fiss    ! Should we get fiss data?
      integer, intent(in)                :: max_order   ! Maximum requested order

      type(Node), pointer     :: node_legendre_mu
      character(MAX_LINE_LEN) :: temp_str
      logical                 :: enable_leg_mu
      real(8), allocatable    :: temp_arr(:)
      real(8), allocatable    :: temp_mult(:,:,:,:)
      real(8), allocatable    :: scatt_coeffs(:,:,:,:,:)
      real(8), allocatable    :: input_scatt(:,:,:,:,:)
      real(8), allocatable    :: temp_scatt(:,:,:,:,:)
      real(8)                 :: dmu, mu, norm, dangle
      integer                 :: order, order_dim, gin, gout, l, arr_len
      integer                 :: legendre_mu_points, imu, ipol, iazi

      ! Call generic data gathering routine (will populate the metadata)
      call nuclidemg_init(this, node_xsdata)

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
      if (this % fissionable) then

        if (check_for_node(node_xsdata,"chi")) then
          ! Get chi
          allocate(temp_arr(groups * this % n_azi * this % n_pol))
          call get_node_array(node_xsdata,"chi",temp_arr)
          allocate(this % chi(groups,this % n_azi,this % n_pol))
          this % chi = reshape(temp_arr,(/groups,this % n_azi,this % n_pol/))
          deallocate(temp_arr)

          ! Get nu_fission (as a vector)
          if (check_for_node(node_xsdata,"nu_fission")) then
            allocate(temp_arr(groups * this % n_azi * this % n_pol))
            call get_node_array(node_xsdata,"nu_fission", temp_arr)
            allocate(this % nu_fission(groups,1,this % n_azi,this % n_pol))
            this % nu_fission = reshape(temp_arr, (/groups,1,this % n_azi, &
                                                    this % n_pol/))
            deallocate(temp_arr)
          else
            call fatal_error("If fissionable, must provide nu_fission!")
          end if

        else
          ! Get nu_fission (as a matrix)
          if (check_for_node(node_xsdata,"nu_fission")) then

            allocate(temp_arr(groups * this % n_azi * this % n_pol))
            call get_node_array(node_xsdata,"nu_fission",temp_arr)
            allocate(this % nu_fission(groups,groups,this % n_azi,this % n_pol))
            this % nu_fission = reshape(temp_arr,(/groups,groups, &
                                                    this % n_azi,this % n_pol/))
            deallocate(temp_arr)
          else
            call fatal_error("If fissionable, must provide nu_fission!")
          end if
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

    end subroutine nuclideangle_init

!===============================================================================
! NUCLIDECE_CLEAR resets and deallocates data in Nuclide, NuclideIso
! or NuclideAngle
!===============================================================================

    subroutine nuclidece_clear(this)

      class(NuclideCE), intent(inout) :: this ! The Nuclide object to clear

      integer :: i ! Loop counter

      if (associated(this % urr_data)) deallocate(this % urr_data)

      if (allocated(this % reactions)) then
        do i = 1, size(this % reactions)
          call this % reactions(i) % clear()
        end do
      end if

      call this % reaction_index % clear()

    end subroutine nuclidece_clear

!===============================================================================
! NUCLIDE*_PRINT displays information about a continuous-energy neutron
! cross_section table and its reactions and secondary angle/energy distributions
!===============================================================================

    subroutine nuclidece_print(this, unit)
      class(NuclideCE), intent(in) :: this
      integer, intent(in), optional :: unit

      integer :: i                 ! loop index over nuclides
      integer :: unit_             ! unit to write to
      integer :: size_xs           ! memory used for cross-sections (bytes)
      integer :: size_urr          ! memory used for probability tables (bytes)
      type(UrrData),  pointer :: urr

      ! set default unit for writing information
      if (present(unit)) then
        unit_ = unit
      else
        unit_ = OUTPUT_UNIT
      end if

      ! Initialize totals
      size_urr = 0
      size_xs = 0

      ! Basic nuclide information
      write(unit_,*) 'Nuclide ' // trim(this % name)
      write(unit_,*) '  zaid = ' // trim(to_str(this % zaid))
      write(unit_,*) '  awr = ' // trim(to_str(this % awr))
      write(unit_,*) '  kT = ' // trim(to_str(this % kT))
      write(unit_,*) '  # of grid points = ' // trim(to_str(this % n_grid))
      write(unit_,*) '  Fissionable = ', this % fissionable
      write(unit_,*) '  # of fission reactions = ' // trim(to_str(this % n_fission))
      write(unit_,*) '  # of reactions = ' // trim(to_str(this % n_reaction))

      ! Information on each reaction
      write(unit_,*) '  Reaction     Q-value  COM    IE'
      do i = 1, this % n_reaction
        associate (rxn => this % reactions(i))
          write(unit_,'(3X,A11,1X,F8.3,3X,L1,3X,I6)') &
               reaction_name(rxn % MT), rxn % Q_value, rxn % scatter_in_cm, &
               rxn % threshold

          ! Accumulate data size
          size_xs = size_xs + (this % n_grid - rxn%threshold + 1) * 8
        end associate
      end do

      ! Add memory required for summary reactions (total, absorption, fission,
      ! nu-fission)
      size_xs = 8 * this % n_grid * 4

      ! Write information about URR probability tables
      size_urr = 0
      if (this % urr_present) then
        urr => this % urr_data
        write(unit_,*) '  Unresolved resonance probability table:'
        write(unit_,*) '    # of energies = ' // trim(to_str(urr % n_energy))
        write(unit_,*) '    # of probabilities = ' // trim(to_str(urr % n_prob))
        write(unit_,*) '    Interpolation =  ' // trim(to_str(urr % interp))
        write(unit_,*) '    Inelastic flag = ' // trim(to_str(urr % inelastic_flag))
        write(unit_,*) '    Absorption flag = ' // trim(to_str(urr % absorption_flag))
        write(unit_,*) '    Multiply by smooth? ', urr % multiply_smooth
        write(unit_,*) '    Min energy = ', trim(to_str(urr % energy(1)))
        write(unit_,*) '    Max energy = ', trim(to_str(urr % energy(urr % n_energy)))

        ! Calculate memory used by probability tables and add to total
        size_urr = urr % n_energy * (urr % n_prob * 6 + 1) * 8
      end if

      ! Write memory used
      write(unit_,*) '  Memory Requirements'
      write(unit_,*) '    Cross sections = ' // trim(to_str(size_xs)) // ' bytes'
      write(unit_,*) '    Probability Tables = ' // &
           trim(to_str(size_urr)) // ' bytes'

      ! Blank line at end of nuclide
      write(unit_,*)
    end subroutine nuclidece_print

    subroutine nuclidemg_print(this, unit_)
      class(NuclideMG), intent(in) :: this
      integer, intent(in)           :: unit_

      character(MAX_LINE_LEN) :: temp_str

      ! Basic nuclide information
      write(unit_,*) 'Nuclide ' // trim(this % name)
      if (this % zaid > 0) then
        ! Dont print if data was macroscopic and thus zaid & AWR would be nonsense
        write(unit_,*) '  zaid = ' // trim(to_str(this % zaid))
        write(unit_,*) '  awr = ' // trim(to_str(this % awr))
      end if
      write(unit_,*) '  kT = ' // trim(to_str(this % kT))
      if (this % scatt_type == ANGLE_LEGENDRE) then
        temp_str = "Legendre"
        write(unit_,*) '  Scattering Type = ' // trim(temp_str)
        select type(this)
        type is (NuclideIso)
          temp_str = to_str(size(this % scatter % dist(1) % data,dim=1) - 1)
        end select
        write(unit_,*) '  Scattering Order = ' // trim(temp_str)
      else if (this % scatt_type == ANGLE_HISTOGRAM) then
        temp_str = "Histogram"
        write(unit_,*) '  Scattering Type = ' // trim(temp_str)
        select type(this)
        type is (NuclideIso)
          temp_str = to_str(size(this % scatter % dist(1) % data,dim=1))
        end select
        write(unit_,*) '  Num. Distribution Bins = ' // trim(temp_str)
      else if (this % scatt_type == ANGLE_TABULAR) then
        temp_str = "Tabular"
        write(unit_,*) '  Scattering Type = ' // trim(temp_str)
        select type(this)
        type is (NuclideIso)
          temp_str = to_str(size(this % scatter % dist(1) % data,dim=1))
        end select
        write(unit_,*) '  Num. Distribution Points = ' // trim(temp_str)
      end if
      write(unit_,*) '  Fissionable = ', this % fissionable

    end subroutine nuclidemg_print

    subroutine nuclideiso_print(this, unit)

      class(NuclideIso), intent(in) :: this
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
      call nuclidemg_print(this, unit_)

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

    end subroutine nuclideiso_print

    subroutine nuclideangle_print(this, unit)

      class(NuclideAngle), intent(in) :: this
      integer, optional, intent(in)    :: unit

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
      call nuclidemg_print(this, unit_)
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

      size_scattmat = (size(this % scatter) + size(this % mult)) * 8
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


    end subroutine nuclideangle_print

!===============================================================================
! NUCLIDE*_GET_XS Returns the requested data type
!===============================================================================

    function nuclideiso_get_xs(this, g, xstype, gout, uvw, mu, iazi, ipol) &
         result(xs)
      class(NuclideIso), intent(in) :: this
      integer, intent(in)            :: g      ! Incoming Energy group
      character(*), intent(in)       :: xstype ! Cross Section Type
      integer, optional, intent(in)  :: gout   ! Outgoing Group
      real(8), optional, intent(in)  :: uvw(3) ! Requested Angle
      real(8), optional, intent(in)  :: mu     ! Change in angle
      integer, optional, intent(in)  :: iazi  ! Azimuthal Index
      integer, optional, intent(in)  :: ipol  ! Polar Index
      real(8)                        :: xs     ! Resultant xs

      xs = ZERO

      if ((xstype == 'nu_fission' .or. xstype == 'fission' .or. xstype =='chi' &
           .or. xstype =='k_fission') .and. (.not. this % fissionable)) then
        return
      end if

      if (present(gout)) then
        select case(xstype)
        case('mult')
          xs = this % scatter % mult(g) % data(gout)
        case('nu_fission')
          xs = this % nu_fission(gout,g)
        case('f_mu', 'f_mu/mult')
          xs = this % scatter % calc_f(g, gout, mu)
          if (xstype == 'f_mu/mult') then
            xs = xs / this % scatter % mult(g) % data(gout)
          end if
        end select
      else
        select case(xstype)
        case('total')
          xs = this % total(g)
        case('absorption')
          xs = this % absorption(g)
        case('fission')
          xs = this % fission(g)
        case('k_fission')
          if (allocated(this % k_fission)) then
            xs = this % k_fission(g)
          end if
        case('chi')
          xs = this % chi(g)
        case('scatter')
          xs = this % scatter % scattxs(g)
        end select
      end if
    end function nuclideiso_get_xs

    function nuclideangle_get_xs(this, g, xstype, gout, uvw, mu, iazi, ipol) &
         result(xs)
      class(NuclideAngle), intent(in) :: this
      integer, intent(in)              :: g      ! Incoming Energy group
      character(*), intent(in)         :: xstype ! Cross Section Type
      integer, optional, intent(in)    :: gout   ! Outgoing Group
      real(8), optional, intent(in)    :: mu     ! Change in angle
      real(8), optional, intent(in)    :: uvw(3) ! Requested Angle
      integer, optional, intent(in)    :: iazi  ! Azimuthal Index
      integer, optional, intent(in)    :: ipol  ! Polar Index
      real(8)                          :: xs     ! Resultant xs

      integer :: iazi_, ipol_

      xs = ZERO

      if ((xstype == 'nu_fission' .or. xstype == 'fission' .or. xstype =='chi' &
           .or. xstype =='k_fission') .and. (.not. this % fissionable)) then
        return
      end if

      if (present(iazi) .and. present(ipol)) then
        iazi_ = iazi
        ipol_ = ipol
      else
        call find_angle(this % polar, this % azimuthal, uvw, iazi_, ipol_)
      end if

      if (present(gout)) then
        select case(xstype)
        case('mult')
          xs = this % scatter(iazi_,ipol_) % obj % mult(g) % data(gout)
        case('nu_fission')
          xs = this % nu_fission(gout,g,iazi_,ipol_)
        case('chi')
          xs = this % chi(gout,iazi_,ipol_)
        case('f_mu', 'f_mu/mult')
          xs = this % scatter(iazi_,ipol_) % obj % calc_f(g,gout,mu)
          if (xstype == 'f_mu/mult') then
            xs = xs / this % scatter(iazi_,ipol_) % obj % mult(g) % data(gout)
          end if
        end select
      else
        select case(xstype)
        case('total')
          xs = this % total(g,iazi_,ipol_)
        case('absorption')
          xs = this % absorption(g,iazi_,ipol_)
        case('fission')
          xs = this % fission(g,iazi_,ipol_)
        case('k_fission')
          if (allocated(this % k_fission)) then
            xs = this % k_fission(g,iazi_,ipol_)
          end if
        case('chi')
          xs = this % chi(g,iazi_,ipol_)
        case('scatter')
          xs = this % scatter(iazi_,ipol_) % obj % scattxs(g)
        end select
      end if

    end function nuclideangle_get_xs

end module nuclide_header
