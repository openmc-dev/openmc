module mgxs_header

  use, intrinsic :: ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T, HSIZE_T, SIZE_T

  use algorithm,       only: find, sort
  use constants,       only: MAX_WORD_LEN, ZERO, ONE, TWO, PI
  use error,           only: fatal_error
  use hdf5_interface
  use material_header, only: material
  use math,            only: evaluate_legendre
  use nuclide_header,  only: MaterialMacroXS
  use random_lcg,      only: prn
  use scattdata_header
  use string
  use stl_vector,      only: VectorInt, VectorReal

!===============================================================================
! XS* contains the temperature-dependent cross section data for an MGXS
!===============================================================================

  type :: XsDataIso
    ! Microscopic cross sections
    real(8), allocatable :: total(:)         ! total cross section
    real(8), allocatable :: absorption(:)    ! absorption cross section
    class(ScattData), allocatable :: scatter ! scattering info
    real(8), allocatable :: delayed_nu_fission(:,:) ! Delayed fission matrix (Dg x Gin)
    real(8), allocatable :: prompt_nu_fission(:)    ! Prompt fission vector (Gin)
    real(8), allocatable :: kappa_fission(:)        ! Kappa fission
    real(8), allocatable :: fission(:)              ! Neutron production
    real(8), allocatable :: decay_rate(:)           ! Delayed neutron precursor decay rate
    real(8), allocatable :: inverse_velocity(:)     ! Inverse neutron velocity
    real(8), allocatable :: chi_delayed(:, :, :)    ! Delayed fission spectra
    real(8), allocatable :: chi_prompt(:, :)        ! Prompt fission spectra
  end type XsDataIso

  type :: XsDataAngle
    ! Microscopic cross sections
    ! In all cases, right-most indices are theta, phi
    real(8), allocatable :: total(:, :, :)        ! total cross section
    real(8), allocatable :: absorption(:, :, :)   ! absorption cross section
    type(ScattDataContainer), allocatable :: scatter(:, :) ! scattering info
    real(8), allocatable :: delayed_nu_fission(:, :, :, :) ! Delayed fission matrix (Gout x Gin)
    real(8), allocatable :: prompt_nu_fission(:, :, :)     ! Prompt fission matrix (Gout x Gin)
    real(8), allocatable :: kappa_fission(:, :, :)         ! Kappa fission
    real(8), allocatable :: fission(:, :, :)               ! Neutron production
    real(8), allocatable :: decay_rate(:, :, :)            ! Delayed neutron precursor decay rate
    real(8), allocatable :: inverse_velocity(:, :, :)      ! Inverse neutron velocity
    real(8), allocatable :: chi_delayed(:, :, :, :, :)     ! Delayed fission spectra
    real(8), allocatable :: chi_prompt(:, :, :, :)         ! Prompt fission spectra
  end type XsDataAngle

!===============================================================================
! MGXS contains the base mgxs data for a nuclide/material
!===============================================================================

  type, abstract :: Mgxs
    character(len=MAX_WORD_LEN) :: name   ! name of dataset, e.g. UO2
    real(8)                     :: awr    ! Atomic Weight Ratio
    real(8), allocatable        :: kTs(:) ! temperature in eV (k*T)

    ! Fission information
    logical :: fissionable  ! mgxs object is fissionable?
    integer :: scatter_format ! either legendre, histogram, or tabular.
    integer :: num_delayed_groups ! Num delayed groups

    ! Caching information
    integer :: index_temp ! temperature index for nuclide

  contains
    procedure(mgxs_from_hdf5_), deferred :: from_hdf5 ! Load the data
    procedure(mgxs_combine_),  deferred  :: combine   ! initializes object
    procedure(mgxs_get_xs_), deferred    :: get_xs  ! Get the requested xs

    ! Sample the outgoing energy from a fission event
    procedure(mgxs_sample_fission_), deferred :: sample_fission_energy

    ! Sample the outgoing energy and angle from a scatter event
    procedure(mgxs_sample_scatter_), deferred :: sample_scatter

    ! Calculate the material specific MGXS data from the nuclides
    procedure(mgxs_calculate_xs_), deferred   :: calculate_xs

    ! Find the temperature
    procedure :: find_temperature => mgxs_find_temperature
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
    subroutine mgxs_from_hdf5_(this, xs_id, energy_groups, delayed_groups, &
         temperature, method, tolerance, max_order, legendre_to_tabular, &
         legendre_to_tabular_points)
      import Mgxs, HID_T, VectorReal
      class(Mgxs), intent(inout)   :: this        ! Working Object
      integer(HID_T), intent(in)   :: xs_id       ! Library data
      integer, intent(in)          :: energy_groups  ! Number of energy groups
      integer, intent(in)          :: delayed_groups ! Number of delayed groups
      type(VectorReal), intent(in) :: temperature ! list of desired temperatures
      integer, intent(inout)       :: method      ! Type of temperature access
      real(8), intent(in)          :: tolerance   ! Tolerance on method
      integer, intent(in)          :: max_order   ! Maximum requested order
      logical, intent(in)          :: legendre_to_tabular ! Convert Legendres to Tabular?
      integer, intent(in)          :: legendre_to_tabular_points ! Number of points to use
                                                                 ! in that  conversion
    end subroutine mgxs_from_hdf5_

    subroutine mgxs_combine_(this, temps, mat, nuclides, energy_groups, &
         delayed_groups, max_order, tolerance, method)
      import Mgxs, Material, MgxsContainer, VectorReal
      class(Mgxs), intent(inout)          :: this ! The Mgxs to initialize
      type(VectorReal), intent(in)        :: temps ! Temperatures to obtain
      type(Material), pointer, intent(in) :: mat  ! base material
      type(MgxsContainer), intent(in)     :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                 :: energy_groups ! Number of energy groups
      integer, intent(in)                 :: delayed_groups ! Number of delayed groups
      integer, intent(in)                 :: max_order  ! Maximum requested order
      real(8), intent(in)                 :: tolerance  ! Tolerance on method
      integer, intent(in)                 :: method     ! Type of temperature access
    end subroutine mgxs_combine_

    pure function mgxs_get_xs_(this, xstype, gin, gout, uvw, mu, dg) result(xs_val)
      import Mgxs
      class(Mgxs), intent(in)       :: this
      character(*), intent(in)      :: xstype ! Cross Section Type
      integer, intent(in)           :: gin    ! Incoming Energy group
      integer, optional, intent(in) :: gout   ! Outgoing Group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      integer, optional, intent(in) :: dg     ! Delayed group
      real(8)                       :: xs_val ! Resultant xs
    end function mgxs_get_xs_

    subroutine mgxs_sample_fission_(this, gin, uvw, dg, gout)
      import Mgxs
      class(Mgxs), intent(in) :: this
      integer, intent(in)     :: gin    ! Incoming energy group
      real(8), intent(in)     :: uvw(3) ! Particle Direction
      integer, intent(out)    :: dg     ! Delayed group
      integer, intent(out)    :: gout   ! Sampled outgoing group

    end subroutine mgxs_sample_fission_

    subroutine mgxs_sample_scatter_(this, uvw, gin, gout, mu, wgt)
      import Mgxs
      class(Mgxs), intent(in)       :: this
      real(8),        intent(in)    :: uvw(3) ! Incoming neutron direction
      integer,        intent(in)    :: gin    ! Incoming neutron group
      integer,        intent(out)   :: gout   ! Sampled outgoin group
      real(8),        intent(out)   :: mu     ! Sampled change in angle
      real(8),        intent(inout) :: wgt    ! Particle weight
    end subroutine mgxs_sample_scatter_

    subroutine mgxs_calculate_xs_(this, gin, sqrtkT, uvw, xs)
      import Mgxs, MaterialMacroXS
      class(Mgxs),           intent(inout) :: this
      integer,               intent(in)    :: gin    ! Incoming neutron group
      real(8),               intent(in)    :: sqrtkT ! Material temperature
      real(8),               intent(in)    :: uvw(3) ! Incoming neutron direction
      type(MaterialMacroXS), intent(inout) :: xs     ! Resultant Mgxs Data
    end subroutine mgxs_calculate_xs_
  end interface

!===============================================================================
! MGXSISO contains the base MGXS data specifically for
! isotropically weighted MGXS
!===============================================================================

  type, extends(Mgxs) :: MgxsIso
    type(XsDataIso), allocatable :: xs(:) ! One for every temperature
  contains
    procedure :: from_hdf5   => mgxsiso_from_hdf5 ! Initialize Nuclidic MGXS Data
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
    type(XsDataAngle), allocatable :: xs(:)    ! One for every temperature
    integer                        :: n_pol    ! Number of polar angles
    integer                        :: n_azi    ! Number of azimuthal angles
    real(8), allocatable           :: polar(:) ! polar angles
    real(8), allocatable           :: azimuthal(:) ! azimuthal angles

  contains
    procedure :: from_hdf5   => mgxsang_from_hdf5 ! Initialize Nuclidic MGXS Data
    procedure :: get_xs      => mgxsang_get_xs  ! Gets Size of Data w/in Object
    procedure :: combine     => mgxsang_combine ! inits object
    procedure :: sample_fission_energy => mgxsang_sample_fission_energy
    procedure :: sample_scatter => mgxsang_sample_scatter
    procedure :: calculate_xs => mgxsang_calculate_xs
  end type MgxsAngle

  ! Cross section arrays
  type(MgxsContainer), allocatable, target :: nuclides_MG(:)

  ! Cross section caches
  type(MgxsContainer), target, allocatable :: macro_xs(:)

  ! Number of energy groups
  integer :: num_energy_groups

  ! Number of delayed groups
  integer :: num_delayed_groups

  ! Energy group structure with decreasing energy
  real(8), allocatable :: energy_bins(:)

  ! Midpoint of the energy group structure
  real(8), allocatable :: energy_bin_avg(:)

  ! Energy group structure with increasing energy
  real(8), allocatable :: rev_energy_bins(:)

contains

!===============================================================================
! MGXS*_FROM_HDF5 reads in the data from the HDF5 Library. At the point of entry
! the file would have been opened and metadata read.
!===============================================================================

    subroutine mgxs_from_hdf5(this, xs_id, temperature, method, tolerance, &
                              temps_to_read, order_dim)
      class(Mgxs), intent(inout)     :: this          ! Working Object
      integer(HID_T), intent(in)     :: xs_id         ! Group in H5 file
      type(VectorReal), intent(in)   :: temperature   ! list of desired temperatures
      integer, intent(inout)         :: method        ! Type of temperature access
      real(8), intent(in)            :: tolerance     ! Tolerance on method
      type(VectorInt), intent(out)   :: temps_to_read ! Temperatures to read
      integer, intent(out)           :: order_dim     ! Scattering data order size

      integer(SIZE_T) :: name_len
      integer(HID_T) :: kT_group
      character(MAX_WORD_LEN), allocatable :: dset_names(:)
      real(8), allocatable :: temps_available(:) ! temperatures available
      real(8) :: temp_desired
      real(8) :: temp_actual
      character(MAX_WORD_LEN) :: temp_str
      real(8) :: dangle
      integer :: ipol, iazi

      ! Get name of dataset from group
      name_len = len(this % name)
      this % name = get_name(xs_id, name_len)

      ! Get rid of leading '/'
      this % name = trim(this % name(2:))

      if (attribute_exists(xs_id, "atomic_weight_ratio")) then
        call read_attribute(this % awr, xs_id, "atomic_weight_ratio")
      else
        this % awr = -ONE
      end if

      ! Determine temperatures available
      kT_group = open_group(xs_id, 'kTs')
      call get_datasets(kT_group, dset_names)
      allocate(temps_available(size(dset_names)))
      do i = 1, size(dset_names)
        ! Read temperature value
        call read_dataset(temps_available(i), kT_group, trim(dset_names(i)))
        ! Convert eV to Kelvin
        temps_available(i) = temps_available(i) / K_BOLTZMANN
      end do
      call sort(temps_available)

      ! If only one temperature is available, revert to nearest temperature
      if (size(temps_available) == 1 .and. &
           method == TEMPERATURE_INTERPOLATION) then
        call warning("Cross sections for " // trim(this % name) // " are only &
             &available at one temperature. Reverting to nearest temperature &
             &method.")
        method = TEMPERATURE_NEAREST
      end if

      select case (method)
      case (TEMPERATURE_NEAREST)
        ! Determine actual temperatures to read
        TEMP_LOOP: do i = 1, temperature % size()
          temp_desired = temperature % data(i)
          i_closest = minloc(abs(temps_available - temp_desired), dim=1)
          temp_actual = temps_available(i_closest)
          if (abs(temp_actual - temp_desired) < tolerance) then
            if (find(temps_to_read, nint(temp_actual)) == -1) then
              call temps_to_read % push_back(nint(temp_actual))
            end if
          else
            call fatal_error("MGXS library does not contain cross sections &
                 &for " // trim(this % name) // " at or near " // &
                 trim(to_str(nint(temp_desired))) // " K.")
          end if
        end do TEMP_LOOP

      case (TEMPERATURE_INTERPOLATION)
        ! If temperature interpolation or multipole is selected, get a list of
        ! bounding temperatures for each actual temperature present in the model
        TEMPS_LOOP: do i = 1, temperature % size()
          temp_desired = temperature % data(i)

          do j = 1, size(temps_available) - 1
            if (temps_available(j) <= temp_desired .and. &
                 temp_desired < temps_available(j + 1)) then
              if (find(temps_to_read, nint(temps_available(j))) == -1) then
                call temps_to_read % push_back(nint(temps_available(j)))
              end if
              if (find(temps_to_read, nint(temps_available(j + 1))) == -1) then
                call temps_to_read % push_back(nint(temps_available(j + 1)))
              end if
              cycle TEMPS_LOOP
            end if
          end do

          call fatal_error("MGXS library does not contain cross sections &
               &for " // trim(this % name) // " at temperatures that bound " // &
               trim(to_str(nint(temp_desired))) // " K.")
        end do TEMPS_LOOP
      end select

      ! Sort temperatures to read
      call sort(temps_to_read)

      ! Get temperatures
      n_temperature = temps_to_read % size()
      allocate(this % kTs(n_temperature))
      do i = 1, n_temperature
        ! Get temperature as a string
        temp_str = trim(to_str(temps_to_read % data(i))) // "K"

        ! Read exact temperature value
        call read_dataset(this % kTs(i), kT_group, trim(temp_str))
      end do
      call close_group(kT_group)

      ! Allocate the XS object for the number of temperatures
      select type(this)
      type is (MgxsIso)
        allocate(this % xs(n_temperature))
      type is (MgxsAngle)
        allocate(this % xs(n_temperature))
      end select

      ! Load the remaining metadata
      if (attribute_exists(xs_id, "scatter_format")) then
        call read_attribute(temp_str, xs_id, "scatter_format")
        temp_str = trim(temp_str)
        if (to_lower(temp_str) == 'legendre') then
          this % scatter_format = ANGLE_LEGENDRE
        else if (to_lower(temp_str) == 'histogram') then
          this % scatter_format = ANGLE_HISTOGRAM
        else if (to_lower(temp_str) == 'tabular') then
          this % scatter_format = ANGLE_TABULAR
        else
          call fatal_error("Invalid scatter_format option!")
        end if
      else
        this % scatter_format = ANGLE_LEGENDRE
      end if
      if (attribute_exists(xs_id, "scatter_shape")) then
        call read_attribute(temp_str, xs_id, "scatter_shape")
        temp_str = trim(temp_str)
        if (to_lower(temp_str) /= "[g][g'][order]") then
          call fatal_error("Invalid scatter_shape option!")
        end if
      end if
      if (attribute_exists(xs_id, "fissionable")) then
        call read_attribute(this % fissionable, xs_id, "fissionable")
      else
        call fatal_error("Fissionable element must be set!")
      end if

      ! Get the library's value for the order
      if (attribute_exists(xs_id, "order")) then
        call read_attribute(order_dim, xs_id, "order")
      else
        call fatal_error("Order must be provided!")
      end if

      ! Store the dimensionality of the data in order_dim.
      ! For Legendre data, we usually refer to it as Pn where n is the order.
      ! However Pn has n+1 sets of points (since you need to
      ! the count the P0 moment).  Adjust for that.  Histogram and Tabular
      ! formats dont need this adjustment.
      if (this % scatter_format == ANGLE_LEGENDRE) then
        order_dim = order_dim + 1
      else
        order_dim = order_dim
      end if

      ! Get angular meta-data and allocate as needed based off of the
      ! information therein
      select type(this)
      type is (MgxsAngle)
        if (attribute_exists(xs_id, "num_polar")) then
          call read_attribute(this % n_pol, xs_id, "num_polar")
        else
          call fatal_error("num_polar must be provided!")
        end if

        if (attribute_exists(xs_id, "num_azimuthal")) then
          call read_attribute(this % n_azi, xs_id, "num_azimuthal")
        else
          call fatal_error("num_azimuthal must be provided!")
        end if

        ! Set angle data to use equally-spaced bins
        allocate(this % polar(this % n_pol))
        dangle = PI / real(this % n_pol, 8)
        do ipol = 1, this % n_pol
          this % polar(ipol) = (real(ipol, 8) - HALF) * dangle
        end do
        allocate(this % azimuthal(this % n_azi))
        dangle = TWO * PI / real(this % n_azi, 8)
        do iazi = 1, this % n_azi
          this % azimuthal(iazi) = -PI + (real(iazi, 8) - HALF) * dangle
        end do
      end select

    end subroutine mgxs_from_hdf5

    subroutine mgxsiso_from_hdf5(this, xs_id, energy_groups, delayed_groups, &
         temperature, method, tolerance, max_order, &
         legendre_to_tabular, legendre_to_tabular_points)
      class(MgxsIso), intent(inout) :: this        ! Working Object
      integer(HID_T), intent(in)    :: xs_id       ! Group in H5 file
      integer, intent(in)           :: energy_groups  ! Number of energy groups
      integer, intent(in)           :: delayed_groups ! Number of delayed groups
      type(VectorReal), intent(in)  :: temperature ! list of desired temperatures
      integer, intent(inout)        :: method      ! Type of temperature access
      real(8), intent(in)           :: tolerance   ! Tolerance on method
      integer, intent(in)           :: max_order   ! Maximum requested order
      logical, intent(in)           :: legendre_to_tabular ! Convert Legendres to Tabular?
      integer, intent(in)           :: legendre_to_tabular_points ! Number of points to use
                                                                  ! in that  conversion

      character(MAX_LINE_LEN)     :: temp_str
      integer(HID_T)              :: xsdata, xsdata_grp, scatt_grp
      integer                     :: ndims
      integer(HSIZE_T)            :: dims(2)
      real(8), allocatable        :: temp_arr(:), temp_2d(:, :)
      real(8), allocatable        :: temp_beta(:, :), temp_3d(:, :, :)
      real(8)                     :: dmu, mu, norm, chi_sum
      integer                     :: order, order_dim, gin, gout, l, imu, length
      type(VectorInt)             :: temps_to_read
      integer                     :: t, dg, order_data
      type(Jagged2D), allocatable :: input_scatt(:), scatt_coeffs(:)
      type(Jagged1D), allocatable :: temp_mult(:)
      integer, allocatable        :: gmin(:), gmax(:)

      ! Call generic data gathering routine (will populate the metadata)
      call mgxs_from_hdf5(this, xs_id, temperature, method, tolerance, &
                          temps_to_read, order_data)

      ! Set the number of delayed groups
      this % num_delayed_groups = delayed_groups

      ! Load the more specific data
      do t = 1, temps_to_read % size()
        associate(xs => this % xs(t))

          ! Get temperature as a string
          temp_str = trim(to_str(temps_to_read % data(t))) // "K"
          xsdata_grp = open_group(xs_id, trim(temp_str))

          ! Allocate data for all the cross sections
          allocate(xs % total(energy_groups))
          allocate(xs % absorption(energy_groups))
          allocate(xs % delayed_nu_fission(delayed_groups, energy_groups))
          allocate(xs % prompt_nu_fission(energy_groups))
          allocate(xs % fission(energy_groups))
          allocate(xs % kappa_fission(energy_groups))
          allocate(xs % decay_rate(delayed_groups))
          allocate(xs % inverse_velocity(energy_groups))
          allocate(xs % chi_delayed(delayed_groups, energy_groups, &
               energy_groups))
          allocate(xs % chi_prompt(energy_groups, energy_groups))

          ! Set all fissionable terms to zero
          xs % delayed_nu_fission = ZERO
          xs % prompt_nu_fission  = ZERO
          xs % fission            = ZERO
          xs % kappa_fission      = ZERO
          xs % chi_delayed        = ZERO
          xs % chi_prompt         = ZERO
          xs % decay_rate         = ZERO
          xs % inverse_velocity   = ZERO

          if (this % fissionable) then

            ! Allocate temporary array for beta
            allocate(temp_beta(delayed_groups, energy_groups))

            ! Set beta
            if (object_exists(xsdata_grp, "beta")) then

              ! Get the dimensions of the beta dataset
              xsdata = open_dataset(xsdata_grp, "beta")
              call get_ndims(xsdata, ndims)

              ! Beta is input as (delayed_groups)
              if (ndims == 1) then

                ! Allocate temporary array for beta
                allocate(temp_arr(delayed_groups))

                ! Read beta
                call read_dataset(temp_arr, xsdata_grp, "beta")

                do dg = 1, delayed_groups
                  do gin = 1, energy_groups
                    temp_beta(dg, gin) = temp_arr(dg)
                  end do
                end do

                ! Deallocate temporary beta array
                deallocate(temp_arr)

                ! Beta is input as (delayed_groups, energy_groups)
              else if (ndims == 2) then

                ! Allocate temporary array for beta
                allocate(temp_arr(delayed_groups * energy_groups))

                ! Read beta
                call read_dataset(temp_arr, xsdata_grp, "beta")

                ! Reshape array and set to dedicated beta array
                temp_beta = reshape(temp_arr, (/delayed_groups, energy_groups/))

                ! Deallocate temporary beta array
                deallocate(temp_arr)

              else
                call fatal_error("beta must be provided as a 1D or 2D array")
              end if
            else
              temp_beta = ZERO
            end if

            ! If chi provided, set chi-prompt and chi-delayed
            if (object_exists(xsdata_grp, "chi")) then

              ! Allocate temporary array for chi
              allocate(temp_arr(energy_groups))

              ! Read chi
              call read_dataset(temp_arr, xsdata_grp, "chi")

              do gin = 1, energy_groups
                do gout = 1, energy_groups
                  xs % chi_prompt(gout, gin) = temp_arr(gout)
                end do

                ! Normalize chi-prompt so its CDF goes to 1
                chi_sum =sum(xs % chi_prompt(:, gin))
                if (chi_sum == ZERO) then
                  call fatal_error("Encountered chi for a group that sums to &
                       &zero")
                else
                  xs % chi_prompt(:, gin) = xs % chi_prompt(:, gin) / chi_sum
                end if
              end do

              ! Set chi-delayed to chi-prompt
              do dg = 1, delayed_groups
                xs % chi_delayed(dg, :, :) = xs % chi_prompt(:, :)
              end do

              ! Deallocate temporary chi array
              deallocate(temp_arr)
            end if

            ! If nu-fission provided, set prompt-nu_-ission and
            ! delayed-nu-fission. If nu fission is a matrix, set chi-prompt and
            ! chi-delayed.
            if (object_exists(xsdata_grp, "nu-fission")) then

              ! Get the dimensions of the nu-fission dataset
              xsdata = open_dataset(xsdata_grp, "nu-fission")
              call get_ndims(xsdata, ndims)

              ! If nu-fission is a vector
              if (ndims == 1) then

                ! Get nu-fission
                call read_dataset(xs % prompt_nu_fission, xsdata_grp, &
                     "nu-fission")

                ! Set delayed-nu-fission and correct prompt-nu-fission with
                ! beta
                do gin = 1, energy_groups
                  do dg = 1, delayed_groups

                    ! Set delayed-nu-fission using delayed neutron fraction
                    xs % delayed_nu_fission(dg, gin) = temp_beta(dg, gin) * &
                         xs % prompt_nu_fission(gin)
                  end do

                  ! Correct prompt-nu-fission using delayed neutron fraction
                  if (delayed_groups > 0) then
                    xs % prompt_nu_fission(gin) = (1 - sum(temp_beta(:, gin))) &
                         * xs % prompt_nu_fission(gin)
                  end if
                end do

                ! If nu-fission is a matrix, set prompt-nu-fission,
                ! delayed-nu-fission, chi-prompt, and chi-delayed.
              else if (ndims == 2) then

                ! chi is embedded in nu-fission -> extract chi
                allocate(temp_arr(energy_groups * energy_groups))
                call read_dataset(temp_arr, xsdata_grp, "nu-fission")
                allocate(temp_2d(energy_groups, energy_groups))
                temp_2d = reshape(temp_arr, (/energy_groups, energy_groups/))

                ! Deallocate temporary 1D array for nu-fission matrix
                deallocate(temp_arr)

                ! Set the vector nu-fission from the matrix nu-fission
                do gin = 1, energy_groups
                  xs % prompt_nu_fission(gin) = sum(temp_2d(:, gin))
                end do

                ! Set delayed-nu-fission and correct prompt-nu-fission with
                ! beta
                do gin = 1, energy_groups
                  do dg = 1, delayed_groups

                    ! Set delayed-nu-fission using delayed neutron fraction
                    xs % delayed_nu_fission(dg, gin) = temp_beta(dg, gin) * &
                         xs % prompt_nu_fission(gin)
                  end do

                  ! Correct prompt-nu-fission using delayed neutron fraction
                  if (delayed_groups > 0) then
                    xs % prompt_nu_fission(gin) = (1 - sum(temp_beta(:, gin))) &
                         * xs % prompt_nu_fission(gin)
                  end if
                end do

                ! Now pull out information needed for chi
                xs % chi_prompt(:, :) = temp_2d

                ! Deallocate temporary 2D array for nu-fission matrix
                deallocate(temp_2d)

                ! Normalize chi so its CDF goes to 1
                do gin = 1, energy_groups
                  chi_sum = sum(xs % chi_prompt(:, gin))
                  if (chi_sum == ZERO) then
                    call fatal_error("Encountered chi for a group that sums to &
                         &zero")
                  else
                    xs % chi_prompt(:, gin) = xs % chi_prompt(:, gin) / chi_sum
                  end if
                end do

                ! Set chi-delayed to chi-prompt
                do dg = 1, delayed_groups
                  xs % chi_delayed(dg, :, :) = xs % chi_prompt(:, :)
                end do
              else
                call fatal_error("nu-fission must be provided as a 1D or 2D &
                     &array")
              end if
            end if

            ! If chi-prompt provided, set chi-prompt
            if (object_exists(xsdata_grp, "chi-prompt")) then

              ! Allocate temporary array for chi-prompt
              allocate(temp_arr(energy_groups))

              ! Get array with chi-prompt
              call read_dataset(temp_arr, xsdata_grp, "chi-prompt")

              do gin = 1, energy_groups
                do gout = 1, energy_groups
                  xs % chi_prompt(gout, gin) = temp_arr(gout)
                end do

                ! Normalize chi so its CDF goes to 1
                chi_sum = sum(xs % chi_prompt(:, gin))
                if (chi_sum == ZERO) then
                  call fatal_error("Encountered chi prompt for a group that &
                       &sums to zero")
                else
                  xs % chi_prompt(:, gin) = xs % chi_prompt(:, gin) / chi_sum
                end if
              end do

              ! Deallocate temporary array for chi-prompt
              deallocate(temp_arr)
            end if

            ! If chi-delayed provided, set chi-delayed
            if (object_exists(xsdata_grp, "chi-delayed")) then

              ! Get the dimensions of the chi-delayed dataset
              xsdata = open_dataset(xsdata_grp, "chi-delayed")
              call get_ndims(xsdata, ndims)

              ! If chi-delayed is a vector
              if (ndims == 1) then

                ! Allocate temporary array for chi-delayed
                allocate(temp_arr(energy_groups))

                ! Get chi-delayed
                call read_dataset(temp_arr, xsdata_grp, "chi-delayed")

                do dg = 1, delayed_groups
                  do gin = 1, energy_groups
                    do gout = 1, energy_groups
                      xs % chi_delayed(dg, gout, gin) = temp_arr(gout)
                    end do

                    ! Normalize chi so its CDF goes to 1
                    chi_sum = sum(xs % chi_delayed(dg, :, gin))
                    if (chi_sum == ZERO) then
                      call fatal_error("Encountered chi delayed for a group &
                           &that sums to zero")
                    else
                      xs % chi_delayed(dg, :, gin) = &
                           xs % chi_delayed(dg, :, gin) / chi_sum
                    end if
                  end do
                end do

                ! Deallocate temporary array for chi-delayed
                deallocate(temp_arr)

              else if (ndims == 2) then

                ! Allocate temporary array for chi-delayed
                allocate(temp_arr(delayed_groups * energy_groups))

                ! Get chi-delayed
                call read_dataset(temp_arr, xsdata_grp, "chi-delayed")
                allocate(temp_2d(delayed_groups, energy_groups))
                temp_2d = reshape(temp_arr, (/delayed_groups, energy_groups/))

                do dg = 1, delayed_groups
                  do gin = 1, energy_groups
                    do gout = 1, energy_groups
                      xs % chi_delayed(dg, gout, gin) = temp_2d(dg, gout)
                    end do

                    ! Normalize chi so its CDF goes to 1
                    chi_sum = sum(xs % chi_delayed(dg, :, gin))
                    if (chi_sum == ZERO) then
                      call fatal_error("Encountered chi delayed for a group &
                           &that sums to zero")
                    else
                      xs % chi_delayed(dg, :, gin) = &
                           xs % chi_delayed(dg, :, gin) / chi_sum
                    end if
                  end do
                end do

                ! Deallocate temporary arrays for chi-delayed
                deallocate(temp_arr)
                deallocate(temp_2d)

              else
                call fatal_error("chi-delayed must be provided as a 1D or 2D &
                     &array")
              end if
            end if

            ! If prompt-nu-fission present, set prompt-nu-fission
            if (object_exists(xsdata_grp, "prompt-nu-fission")) then

              ! Get the dimensions of the prompt-nu-fission dataset
              xsdata = open_dataset(xsdata_grp, "prompt-nu-fission")
              call get_ndims(xsdata, ndims)

              ! If prompt-nu-fission is a vector
              if (ndims == 1) then

                ! Set prompt_nu_fission
                call read_dataset(xs % prompt_nu_fission, xsdata_grp, &
                     "prompt-nu-fission")

                ! If prompt-nu-fission is a matrix, set prompt_nu_fission and
                ! chi_prompt.
              else if (ndims == 2) then

                ! chi_prompt is embedded in prompt_nu_fission -> extract
                ! chi_prompt
                allocate(temp_arr(energy_groups * energy_groups))
                call read_dataset(temp_arr, xsdata_grp, "prompt-nu-fission")
                allocate(temp_2d(energy_groups, energy_groups))
                temp_2d = reshape(temp_arr, (/energy_groups, energy_groups/))

                ! Deallocate temporary 1D array for prompt_nu_fission matrix
                deallocate(temp_arr)

                ! Set the vector prompt-nu-fission from the matrix
                ! prompt-nu-fission
                do gin = 1, energy_groups
                  xs % prompt_nu_fission(gin) = sum(temp_2d(:, gin))
                end do

                ! Now pull out information needed for chi
                xs % chi_prompt(:, :) = temp_2d

                ! Deallocate temporary 2D array for nu_fission matrix
                deallocate(temp_2d)

                ! Normalize chi so its CDF goes to 1
                do gin = 1, energy_groups
                  chi_sum = sum(xs % chi_prompt(:, gin))
                  if (chi_sum == ZERO) then
                    call fatal_error("Encountered chi prompt for a group &
                         &that sums to zero")
                  else
                    xs % chi_prompt(:, gin) = xs % chi_prompt(:, gin) / chi_sum
                  end if
                end do
              else
                call fatal_error("prompt-nu-fission must be provided as a 1D &
                     &or 2D array")
              end if
            end if

            ! If delayed-nu-fission provided, set delayed-nu-fission. If
            ! delayed-nu-fission is a matrix, set chi-delayed.
            if (object_exists(xsdata_grp, "delayed-nu-fission")) then

              ! Get the dimensions of the delayed-nu-fission dataset
              xsdata = open_dataset(xsdata_grp, "delayed-nu-fission")
              call get_ndims(xsdata, ndims)

              ! If delayed-nu-fission is a vector
              if (ndims == 1) then

                ! If beta is zeros, raise error
                if (temp_beta(1,1) == ZERO) then
                  call fatal_error("cannot set delayed-nu-fission with a 1D &
                       &array if beta not provided")
                end if

                ! Allocate temporary array for delayed-nu-fission
                allocate(temp_arr(energy_groups))

                ! Get delayed-nu-fission
                call read_dataset(temp_arr, xsdata_grp, "delayed-nu-fission")

                do gin = 1, energy_groups
                  do dg = 1, delayed_groups

                    ! Set delayed-nu-fission using delayed neutron fraction
                    xs % delayed_nu_fission(dg, gin) = temp_beta(dg, gin) * &
                         temp_arr(gin)
                  end do
                end do

                ! Deallocate temporary delayed-nu-fission array
                deallocate(temp_arr)

                ! If delayed-nu-fission is a (delayed_group, energy_group)
                ! matrix, set delayed-nu-fission separately for each delayed
                ! group.
              else if (ndims == 2) then

                ! Get the shape of delayed-nu-fission
                call get_shape(xsdata, dims)

                ! Issue error if 1st dimension not correct
                if (dims(1) /= delayed_groups) then
                  call fatal_error("The delayed-nu-fission matrix was input &
                       &with a 1st dimension not equal to the number of &
                       &delayed groups.")
                end if

                ! Issue error if 2nd dimension not correct
                if (dims(2) /= energy_groups) then
                  call fatal_error("The delayed-nu-fission matrix was input &
                       &with a 2nd dimension not equal to the number of &
                       &energy groups.")
                end if

                ! Issue warning if delayed_groups == energy_groups
                if (delayed_groups == energy_groups) then
                  call warning("delayed-nu-fission was input as a dimension &
                       &2 matrix with the same number of delayed groups and &
                       &groups. It is important to know that OpenMC assumes &
                       &the dimensions in the matrix are (delayed_groups, &
                       &energy_groups). Currently, delayed-nu-fission cannot &
                       &be set as a group by group matrix.")
                end if

                ! Get delayed-nu-fission
                allocate(temp_arr(delayed_groups * energy_groups))
                call read_dataset(temp_arr, xsdata_grp, "delayed-nu-fission")
                xs % delayed_nu_fission = reshape(temp_arr, (/delayed_groups, &
                     energy_groups/))

                ! Deallocate temporary array for delayed-nu-fission matrix
                deallocate(temp_arr)

                ! If delayed nu-fission is a 3D matrix, set delayed_nu_fission
                ! and chi_delayed.
              else if (ndims == 3) then

                ! chi_delayed is embedded in delayed_nu_fission -> extract
                ! chi_delayed
                allocate(temp_arr(delayed_groups * energy_groups * &
                     energy_groups))
                call read_dataset(temp_arr, xsdata_grp, "delayed-nu-fission")
                allocate(temp_3d(delayed_groups, energy_groups, energy_groups))
                temp_3d = reshape(temp_arr, (/delayed_groups, energy_groups, &
                     energy_groups/))

                ! Deallocate temporary 1D array for delayed_nu_fission matrix
                deallocate(temp_arr)

                ! Set the 2D delayed-nu-fission matrix and 3D chi_dealyed matrix
                ! from the 3D delayed-nu-fission matrix
                do dg = 1, delayed_groups
                  do gin = 1, energy_groups
                    xs % delayed_nu_fission(dg, gin) = sum(temp_3d(dg, :, gin))
                    do gout = 1, energy_groups
                      xs % chi_delayed(dg, gout, gin) = temp_3d(dg, gout, gin)
                    end do
                  end do
                end do

                ! Normalize chi_delayed so its CDF goes to 1
                do dg = 1, delayed_groups
                  do gin = 1, energy_groups
                    chi_sum = sum(xs % chi_delayed(dg, :, gin))
                    if (chi_sum == ZERO) then
                      call fatal_error("Encountered chi delayed for a group &
                           &that sums to zero")
                    else
                      xs % chi_delayed(dg, :, gin) = &
                           xs % chi_delayed(dg, :, gin) / chi_sum
                    end if
                  end do
                end do

                ! Deallocate temporary 3D matrix for delayed_nu_fission
                deallocate(temp_3d)
              else
                call fatal_error("delayed-nu-fission must be provided as a &
                     &1D, 2D, or 3D array")
              end if
            end if

            ! Deallocate temporary beta array
            deallocate(temp_beta)

            ! chi-prompt, chi-delayed, prompt-nu-fission, and delayed-nu-fission
            ! have been set; Now we will check for the rest of the XS that are
            ! unique to fissionable isotopes

            ! Get fission xs
            if (object_exists(xsdata_grp, "fission")) then
              call read_dataset(xs % fission, xsdata_grp, "fission")
            end if

            ! Get kappa-fission xs
            if (object_exists(xsdata_grp, "kappa-fission")) then
              call read_dataset(xs % kappa_fission, xsdata_grp, "kappa-fission")
            end if

            ! Get decay rate xs
            if (object_exists(xsdata_grp, "decay rate")) then
              call read_dataset(xs % decay_rate, xsdata_grp, "decay rate")
            end if
          end if

          ! All the XS unique to fissionable isotopes have been set; Now set all
          ! the generation XS

          if (object_exists(xsdata_grp, "absorption")) then
            call read_dataset(xs % absorption, xsdata_grp, "absorption")
          else
            call fatal_error("Must provide absorption!")
          end if

          ! Get inverse velocity
          if (object_exists(xsdata_grp, "inverse-velocity")) then
            call read_dataset(xs % inverse_velocity, xsdata_grp, &
                 "inverse-velocity")
          end if

          ! Get scattering data
          if (.not. object_exists(xsdata_grp, "scatter_data")) &
               call fatal_error("Must provide 'scatter_data'")

          scatt_grp = open_group(xsdata_grp, 'scatter_data')

          ! First get the outgoing group boundary indices
          if (object_exists(scatt_grp, "g_min")) then
            allocate(gmin(energy_groups))
            call read_dataset(gmin, scatt_grp, "g_min")
          else
            call fatal_error("'g_min' for the scatter_data must be provided")
          end if

          if (object_exists(scatt_grp, "g_max")) then
            allocate(gmax(energy_groups))
            call read_dataset(gmax, scatt_grp, "g_max")
          else
            call fatal_error("'g_max' for the scatter_data must be provided")
          end if

          ! Now use this information to find the length of a container array
          ! to hold the flattened data
          length = 0

          do gin = 1, energy_groups
            length = length + order_data * (gmax(gin) - gmin(gin) + 1)
          end do

          ! Allocate flattened array
          allocate(temp_arr(length))

          if (.not. object_exists(scatt_grp, 'scatter_matrix')) &
               call fatal_error("'scatter_matrix' must be provided")
          call read_dataset(temp_arr, scatt_grp, "scatter_matrix")

          ! Compare the number of orders given with the maximum order of the
          ! problem.  Strip off the supefluous orders if needed.
          if (this % scatter_format == ANGLE_LEGENDRE) then
            order = min(order_data - 1, max_order)
            order_dim = order + 1
          else
            order_dim = order_data
          end if

          ! Convert temp_arr to a jagged array ((gin) % data(l, gout)) for
          ! passing to ScattData
          allocate(input_scatt(energy_groups))

          index = 1
          do gin = 1, energy_groups
            allocate(input_scatt(gin) % data(order_dim, gmin(gin):gmax(gin)))
            do gout = gmin(gin), gmax(gin)
              do l = 1, order_dim
                input_scatt(gin) % data(l, gout) = temp_arr(index)
                index = index + 1
              end do
              ! Adjust index for the orders we didnt take
              index = index + (order_data - order_dim)
            end do
          end do

          deallocate(temp_arr)

          ! Finally convert the legendre to tabular if needed
          allocate(scatt_coeffs(energy_groups))

          if (this % scatter_format == ANGLE_LEGENDRE .and. &
               legendre_to_tabular) then

            this % scatter_format = ANGLE_TABULAR
            order_dim = legendre_to_tabular_points
            order = order_dim
            dmu = TWO / real(order - 1, 8)

            do gin = 1, energy_groups
              allocate(scatt_coeffs(gin) % data(order_dim, gmin(gin):gmax(gin)))
              do gout = gmin(gin), gmax(gin)

                norm = ZERO

                do imu = 1, order_dim

                  if (imu == 1) then
                    mu = -ONE
                  else if (imu == order_dim) then
                    mu = ONE
                  else
                    mu = -ONE + real(imu - 1, 8) * dmu
                  end if

                  scatt_coeffs(gin) % data(imu, gout) = &
                       evaluate_legendre(input_scatt(gin) % data(:, gout), mu)

                  ! Ensure positivity of distribution
                  if (scatt_coeffs(gin) % data(imu, gout) < ZERO) &
                       scatt_coeffs(gin) % data(imu, gout) = ZERO

                  ! And accrue the integral
                  if (imu > 1) then
                    norm = norm + HALF * dmu * &
                         (scatt_coeffs(gin) % data(imu - 1, gout) + &
                          scatt_coeffs(gin) % data(imu, gout))
                  end if
                end do ! mu

                ! Now that we have the integral, lets ensure that the
                ! distribution is normalized such that it preserves the original
                ! scattering xs
                if (norm > ZERO) then
                  scatt_coeffs(gin) % data(:, gout) = &
                       scatt_coeffs(gin) % data(:, gout) * &
                       input_scatt(gin) % data(1, gout) / norm
                end if
              end do ! gout
            end do ! gin
          else

            ! Sticking with current representation
            do gin = 1, energy_groups
              allocate(scatt_coeffs(gin) % data(order_dim, gmin(gin):gmax(gin)))
              scatt_coeffs(gin) % data(:, :) = &
                   input_scatt(gin) % data(1:order_dim, :)
            end do
          end if

          deallocate(input_scatt)

          ! Now get the multiplication matrix
          if (object_exists(scatt_grp, 'multiplicity_matrix')) then

            ! Now use this information to find the length of a container array
            ! to hold the flattened data
            length = 0

            do gin = 1, energy_groups
              length = length + (gmax(gin) - gmin(gin) + 1)
            end do

            ! Allocate flattened array
            allocate(temp_arr(length))
            call read_dataset(temp_arr, scatt_grp, "multiplicity_matrix")

            ! Convert temp_arr to a jagged array ((gin) % data(gout)) for
            ! passing to ScattData
            allocate(temp_mult(energy_groups))

            index = 1
            do gin = 1, energy_groups

              allocate(temp_mult(gin) % data(gmin(gin):gmax(gin)))

              do gout = gmin(gin), gmax(gin)
                temp_mult(gin) % data(gout) = temp_arr(index)
                index = index + 1
              end do
            end do
            deallocate(temp_arr)
          else

            ! Default to multiplicities of 1.0
            allocate(temp_mult(energy_groups))

            do gin = 1, energy_groups
              allocate(temp_mult(gin) % data(gmin(gin):gmax(gin)))
              temp_mult(gin) % data = ONE
            end do
          end if

          ! Allocate and initialize our ScattData Object.
          if (this % scatter_format == ANGLE_HISTOGRAM) then
            allocate(ScattDataHistogram :: xs % scatter)
          else if (this % scatter_format == ANGLE_TABULAR) then
            allocate(ScattDataTabular :: xs % scatter)
          else if (this % scatter_format == ANGLE_LEGENDRE) then
            allocate(ScattDataLegendre :: xs % scatter)
          end if

          ! Initialize the ScattData Object
          call xs % scatter % init(gmin, gmax, temp_mult, scatt_coeffs)

          ! Check sigA to ensure it is not 0 since it is
          ! often divided by in the tally routines
          ! (This may happen with Helium data)
          do gin = 1, energy_groups
            if (xs % absorption(gin) == ZERO) xs % absorption(gin) = 1E-10_8
          end do

          ! Get, or infer, total xs data.
          if (object_exists(xsdata_grp, "total")) then
            call read_dataset(xs % total, xsdata_grp, "total")
          else
            xs % total(:) = xs % absorption(:) + xs % scatter % scattxs(:)
          end if

          ! Check sigT to ensure it is not 0 since it is
          ! often divided by in the tally routines
          do gin = 1, energy_groups
            if (xs % total(gin) == ZERO) xs % total(gin) = 1E-10_8
          end do

          ! Close the groups we have opened and deallocate
          call close_group(xsdata_grp)
          call close_group(scatt_grp)
          deallocate(scatt_coeffs, temp_mult)
        end associate ! xs
      end do ! Temperatures

    end subroutine mgxsiso_from_hdf5

    subroutine mgxsang_from_hdf5(this, xs_id, energy_groups, delayed_groups, &
         temperature, method, tolerance, max_order, legendre_to_tabular, &
         legendre_to_tabular_points)
      class(MgxsAngle), intent(inout) :: this        ! Working Object
      integer(HID_T), intent(in)      :: xs_id       ! Group in H5 file
      integer, intent(in)             :: energy_groups  ! Number of energy groups
      integer, intent(in)             :: delayed_groups ! Number of energy groups
      type(VectorReal), intent(in)    :: temperature ! list of desired temperatures
      integer, intent(inout)          :: method      ! Type of temperature access
      real(8), intent(in)             :: tolerance   ! Tolerance on method
      integer, intent(in)             :: max_order   ! Maximum requested order
      logical, intent(in)             :: legendre_to_tabular ! Convert Legendres to Tabular?
      integer, intent(in)             :: legendre_to_tabular_points ! Number of points to use
                                                                    ! in that  conversion

      character(MAX_LINE_LEN)     :: temp_str
      integer(HID_T)              :: xsdata, xsdata_grp, scatt_grp
      integer                     :: ndims
      integer(HSIZE_T)            :: dims(4)
      integer, allocatable        :: int_arr(:)
      real(8), allocatable        :: temp_1d(:), temp_3d(:, :, :)
      real(8), allocatable        :: temp_4d(:, :, :, :), temp_5d(:, :, :, :, :)
      real(8), allocatable        :: temp_beta(:, :, :, :)
      real(8)                     :: dmu, mu, norm, chi_sum
      integer                     :: order, order_dim, gin, gout, l, imu, dg
      type(VectorInt)             :: temps_to_read
      integer                     :: t, length, ipol, iazi, order_data
      type(Jagged2D), allocatable :: input_scatt(:, :, :), scatt_coeffs(:, :, :)
      type(Jagged1D), allocatable :: temp_mult(:, :, :)
      integer, allocatable        :: gmin(:, :, :), gmax(:, :, :)

      ! Call generic data gathering routine (will populate the metadata)
      call mgxs_from_hdf5(this, xs_id, temperature, method, tolerance, &
                          temps_to_read, order_data)

      ! Set the number of delayed groups
      this % num_delayed_groups = delayed_groups

      ! Load the more specific data
      do t = 1, temps_to_read % size()
        associate(xs => this % xs(t))

          ! Get temperature as a string
          temp_str = trim(to_str(temps_to_read % data(t))) // "K"
          xsdata_grp = open_group(xs_id, trim(temp_str))

          ! Load the more specific data
          allocate(xs % prompt_nu_fission(energy_groups, this % n_azi, &
               this % n_pol))
          allocate(xs % delayed_nu_fission(delayed_groups, energy_groups, &
               this % n_azi, this % n_pol))
          allocate(xs % chi_prompt(energy_groups, energy_groups, this % n_azi, &
               this % n_pol))
          allocate(xs % chi_delayed(delayed_groups, energy_groups, &
               energy_groups, this % n_azi, this % n_pol))
          allocate(xs % total(energy_groups, this % n_azi, this % n_pol))
          allocate(xs % absorption(energy_groups, this % n_azi, this % n_pol))
          allocate(xs % fission(energy_groups, this % n_azi, this % n_pol))
          allocate(xs % kappa_fission(energy_groups, this % n_azi, &
               this % n_pol))
          allocate(xs % decay_rate(delayed_groups, this % n_azi, this % n_pol))
          allocate(xs % inverse_velocity(energy_groups, this % n_azi, &
               this % n_pol))

          ! Set all fissionable terms to zero
          xs % delayed_nu_fission = ZERO
          xs % prompt_nu_fission  = ZERO
          xs % fission            = ZERO
          xs % kappa_fission      = ZERO
          xs % chi_delayed        = ZERO
          xs % chi_prompt         = ZERO
          xs % decay_rate         = ZERO
          xs % inverse_velocity   = ZERO

          if (this % fissionable) then

            ! Allocate temporary array for beta
            allocate(temp_beta(delayed_groups, energy_groups, this % n_azi, &
                 this % n_pol))

            ! Set beta
            if (object_exists(xsdata_grp, "beta")) then

              ! Get the dimensions of the beta dataset
              xsdata = open_dataset(xsdata_grp, "beta")
              call get_ndims(xsdata, ndims)

              ! Beta is input as (delayed_groups, n_azi, n_pol)
              if (ndims == 3) then

                ! Allocate temporary arrays for beta
                allocate(temp_1d(delayed_groups * this % n_azi * this % n_pol))
                allocate(temp_3d(delayed_groups, this % n_azi, this % n_pol))

                ! Read beta
                call read_dataset(temp_1d, xsdata_grp, "beta")
                temp_3d = reshape(temp_1d, (/delayed_groups, this % n_azi, &
                     this % n_pol/))

                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do dg = 1, delayed_groups
                      do gin = 1, energy_groups
                        temp_beta(dg, gin, iazi, ipol) = temp_3d(dg, iazi, ipol)
                      end do
                    end do
                  end do
                end do

                ! Deallocate temporary beta arrays
                deallocate(temp_1d)
                deallocate(temp_3d)

                ! Beta is input as (delayed_groups, energy_groups, n_azi, n_pol)
              else if (ndims == 4) then

                ! Allocate temporary array for beta
                allocate(temp_1d(delayed_groups * energy_groups * this % n_azi &
                     * this % n_pol))

                ! Read beta
                call read_dataset(temp_1d, xsdata_grp, "beta")

                ! Reshape array and set to dedicated beta array
                temp_beta = reshape(temp_1d, (/delayed_groups, &
                     energy_groups, this % n_azi, this % n_pol/))

                ! Deallocate temporary beta array
                deallocate(temp_1d)

              else
                call fatal_error("beta must be provided as a 3D or 4D array")
              end if
            else
              temp_beta = ZERO
            end if

            ! If chi provided, set chi-prompt and chi-delayed
            if (object_exists(xsdata_grp, "chi")) then

              ! Allocate temporary array for chi
              allocate(temp_1d(energy_groups * this % n_azi * this % n_pol))
              allocate(temp_3d(energy_groups, this % n_azi, this % n_pol))

              ! Read chi
              call read_dataset(temp_1d, xsdata_grp, "chi")
              temp_3d = reshape(temp_1d, (/energy_groups, this % n_azi, &
                   this % n_pol/))

              do ipol = 1, this % n_pol
                do iazi = 1, this % n_azi
                  do gin = 1, energy_groups
                    do gout = 1, energy_groups
                      xs % chi_prompt(gout, gin, iazi, ipol) = &
                           temp_3d(gout, iazi, ipol)
                    end do

                    ! Normalize chi-prompt so its CDF goes to 1
                    chi_sum = sum(xs % chi_prompt(:, gin, iazi, ipol))
                    if (chi_sum == ZERO) then
                      call fatal_error("Encountered chi for a group that sums&
                           & to zero")
                    else
                      xs % chi_prompt(:, gin, iazi, ipol) = &
                           xs % chi_prompt(:, gin, iazi, ipol) / chi_sum
                    end if
                  end do
                end do
              end do

              ! Set chi-delayed to chi-prompt
              do ipol = 1, this % n_pol
                do iazi = 1, this % n_azi
                  do dg = 1, delayed_groups
                    xs % chi_delayed(dg, :, :, iazi, ipol) = &
                         xs % chi_prompt(:, :, iazi, ipol)
                  end do
                end do
              end do

              ! Deallocate temporary chi arrays
              deallocate(temp_1d)
              deallocate(temp_3d)
            end if

            ! If nu-fission provided, set prompt-nu_-ission and
            ! delayed-nu-fission. If nu fission is a matrix, set chi-prompt and
            ! chi-delayed.
            if (object_exists(xsdata_grp, "nu-fission")) then

              ! Get the dimensions of the nu-fission dataset
              xsdata = open_dataset(xsdata_grp, "nu-fission")
              call get_ndims(xsdata, ndims)

              ! If nu-fission is a 3D array
              if (ndims == 3) then

                ! Get nu-fission
                call read_dataset(xs % prompt_nu_fission, xsdata_grp, &
                     "nu-fission")

                ! Set delayed-nu-fission and correct prompt-nu-fission with
                ! beta
                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do gin = 1, energy_groups
                      do dg = 1, delayed_groups

                        ! Set delayed-nu-fission using delayed neutron fraction
                        xs % delayed_nu_fission(dg, gin, iazi, ipol) = &
                             temp_beta(dg, gin, iazi, ipol) * &
                             xs % prompt_nu_fission(gin, iazi, ipol)
                      end do

                      ! Correct prompt-nu-fission using delayed neutron fraction
                      if (delayed_groups > 0) then
                        xs % prompt_nu_fission(gin, iazi, ipol) = &
                             (1 - sum(temp_beta(:, gin, iazi, ipol))) * &
                             xs % prompt_nu_fission(gin, iazi, ipol)
                      end if
                    end do
                  end do
                end do

                ! If nu-fission is a matrix, set prompt-nu-fission,
                ! delayed-nu-fission, chi-prompt, and chi-delayed.
              else if (ndims == 4) then

                ! chi is embedded in nu-fission -> extract chi
                allocate(temp_1d(energy_groups * energy_groups * &
                     this % n_azi * this % n_pol))
                call read_dataset(temp_1d, xsdata_grp, "nu-fission")
                allocate(temp_4d(energy_groups, energy_groups, this % n_azi, &
                     this % n_pol))
                temp_4d = reshape(temp_1d, (/energy_groups, energy_groups, &
                     this % n_azi, this % n_pol /))

                ! Deallocate temporary 1D array for nu-fission matrix
                deallocate(temp_1d)

                ! Set the vector nu-fission from the matrix nu-fission
                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do gin = 1, energy_groups
                      xs % prompt_nu_fission(gin, iazi, ipol) = &
                           sum(temp_4d(:, gin, iazi, ipol))
                    end do
                  end do
                end do

                ! Set delayed-nu-fission and correct prompt-nu-fission with
                ! beta
                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do gin = 1, energy_groups
                      do dg = 1, delayed_groups

                        ! Set delayed-nu-fission using delayed neutron fraction
                        xs % delayed_nu_fission(dg, gin, iazi, ipol) = &
                             temp_beta(dg, gin, iazi, ipol) * &
                             xs % prompt_nu_fission(gin, iazi, ipol)
                      end do

                      ! Correct prompt-nu-fission using delayed neutron fraction
                      if (delayed_groups > 0) then
                        xs % prompt_nu_fission(gin, iazi, ipol) = &
                             (1 - sum(temp_beta(:, gin, iazi, ipol))) * &
                             xs % prompt_nu_fission(gin, iazi, ipol)
                      end if
                    end do
                  end do
                end do

                ! Now pull out information needed for chi
                xs % chi_prompt(:, :, :, :) = temp_4d

                ! Deallocate temporary 4D array for nu-fission matrix
                deallocate(temp_4d)

                ! Normalize chi so its CDF goes to 1
                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do gin = 1, energy_groups
                      chi_sum = sum(xs % chi_prompt(:, gin, iazi, ipol))
                      if (chi_sum == ZERO) then
                        call fatal_error("Encountered chi for a group that &
                             &sums to zero")
                      else
                        xs % chi_prompt(:, gin, iazi, ipol) = &
                             xs % chi_prompt(:, gin, iazi, ipol) / chi_sum
                      end if
                    end do

                    ! Set chi-delayed to chi-prompt
                    do dg = 1, delayed_groups
                      xs % chi_delayed(dg, :, :, iazi, ipol) = &
                           xs % chi_prompt(:, :, iazi, ipol)
                    end do
                  end do
                end do
              else
                call fatal_error("nu-fission must be provided as a 3D or &
                     &4D array")
              end if
            end if

            ! If chi-prompt provided, set chi-prompt
            if (object_exists(xsdata_grp, "chi-prompt")) then

              ! Allocate temporary array for chi-prompt
              allocate(temp_1d(energy_groups * this % n_azi * this % n_pol))
              allocate(temp_3d(energy_groups, this % n_azi, this % n_pol))

              ! Get array with chi-prompt
              call read_dataset(temp_1d, xsdata_grp, "chi-prompt")
              temp_3d = reshape(temp_1d, (/energy_groups, this % n_azi, &
                   this % n_pol/))

              do ipol = 1, this % n_pol
                do iazi = 1, this % n_azi
                  do gin = 1, energy_groups
                    do gout = 1, energy_groups
                      xs % chi_prompt(gout, gin, iazi, ipol) = &
                           temp_3d(gout, iazi, ipol)
                    end do

                    ! Normalize chi so its CDF goes to 1
                    chi_sum = sum(xs % chi_prompt(:, gin, iazi, ipol))
                    if (chi_sum == ZERO) then
                      call fatal_error("Encountered chi prompt for a group that&
                           & sums to zero")
                    else
                      xs % chi_prompt(:, gin, iazi, ipol) = &
                           xs % chi_prompt(:, gin, iazi, ipol) / chi_sum
                    end if
                  end do
                end do
              end do

              ! Deallocate temporary arrays for chi-prompt
              deallocate(temp_1d)
              deallocate(temp_3d)
            end if

            ! If chi-delayed provided, set chi-delayed
            if (object_exists(xsdata_grp, "chi-delayed")) then

              ! Get the dimensions of the chi-delayed dataset
              xsdata = open_dataset(xsdata_grp, "chi-delayed")
              call get_ndims(xsdata, ndims)

              ! chi-delayed is input as (energy_groups, n_azi, n_pol)
              if (ndims == 3) then

                ! Allocate temporary array for chi-prompt
                allocate(temp_1d(energy_groups * this % n_azi * this % n_pol))
                allocate(temp_3d(energy_groups, this % n_azi, this % n_pol))

                ! Get array with chi-prompt
                call read_dataset(temp_1d, xsdata_grp, "chi-delayed")
                temp_3d = reshape(temp_1d, (/energy_groups, this % n_azi, &
                     this % n_pol/))

                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do dg = 1, delayed_groups
                      do gin = 1, energy_groups
                        do gout = 1, energy_groups
                          xs % chi_delayed(dg, gout, gin, iazi, ipol) = &
                               temp_3d(gout, iazi, ipol)
                        end do

                        ! Normalize chi so its CDF goes to 1
                        chi_sum = sum(xs % chi_delayed(dg, :, gin, iazi, ipol))
                        if (chi_sum == ZERO) then
                          call fatal_error("Encountered chi delayed for a group&
                               & that sums to zero")
                        else
                          xs % chi_delayed(dg, :, gin, iazi, ipol) = &
                               xs % chi_delayed(dg, :, gin, iazi, ipol) / &
                               chi_sum
                        end if
                      end do
                    end do
                  end do
                end do

                ! Deallocate temporary arrays for chi-delayed
                deallocate(temp_1d)
                deallocate(temp_3d)

                ! chi-delayed is input as (delayed_groups, energy_groups, n_azi,
                ! n_pol)
              else if (ndims == 4) then

                ! Allocate temporary array for chi-delayed
                allocate(temp_1d(delayed_groups * energy_groups * this % n_azi &
                     * this % n_pol))
                allocate(temp_4d(delayed_groups, energy_groups, this % n_azi, &
                     this % n_pol))

                ! Get chi-delayed
                call read_dataset(temp_1d, xsdata_grp, "chi-delayed")
                temp_4d = reshape(temp_1d, (/delayed_groups, energy_groups, &
                     this % n_azi, this % n_pol/))

                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do dg = 1, delayed_groups
                      do gin = 1, energy_groups
                        do gout = 1, energy_groups
                          xs % chi_delayed(dg, gout, gin, iazi, ipol) = &
                               temp_4d(dg, gout, iazi, ipol)
                        end do

                        ! Normalize chi so its CDF goes to 1
                        chi_sum = sum(xs % chi_delayed(dg, :, gin, iazi, ipol))
                        if (chi_sum == ZERO) then
                          call fatal_error("Encountered chi delayed for a group&
                               & that sums to zero")
                        else
                          xs % chi_delayed(dg, :, gin, iazi, ipol) = &
                               xs % chi_delayed(dg, :, gin, iazi, ipol) / &
                               chi_sum
                        end if
                      end do
                    end do
                  end do
                end do

                ! Deallocate temporary arrays for chi-delayed
                deallocate(temp_1d)
                deallocate(temp_4d)

              else
                call fatal_error("chi-delayed must be provided as a 3D or 4D &
                     &array")
              end if
            end if

            ! If prompt-nu-fission present, set prompt-nu-fission
            if (object_exists(xsdata_grp, "prompt-nu-fission")) then

              ! Get the dimensions of the prompt-nu-fission dataset
              xsdata = open_dataset(xsdata_grp, "prompt-nu-fission")
              call get_ndims(xsdata, ndims)

              ! If prompt-nu-fission is a vector for each azi and pol
              if (ndims == 3) then

                ! Set prompt_nu_fission
                call read_dataset(xs % prompt_nu_fission, xsdata_grp, &
                     "prompt-nu-fission")

                ! If prompt-nu-fission is a matrix for each azi and pol,
                ! set prompt_nu_fission and chi_prompt.
              else if (ndims == 4) then

                ! chi_prompt is embedded in prompt_nu_fission -> extract
                ! chi_prompt
                allocate(temp_1d(energy_groups * energy_groups &
                     * this % n_azi * this % n_pol))
                allocate(temp_4d(energy_groups, energy_groups, this % n_azi, &
                     this % n_pol))
                call read_dataset(temp_1d, xsdata_grp, "prompt-nu-fission")
                temp_4d = reshape(temp_1d, (/energy_groups, energy_groups, &
                     this % n_azi, this % n_pol/))

                ! Deallocate temporary 1D array for prompt_nu_fission matrix
                deallocate(temp_1d)

                ! Set the vector prompt-nu-fission from the matrix
                ! prompt-nu-fission
                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do gin = 1, energy_groups
                      xs % prompt_nu_fission(gin, iazi, ipol) = &
                           sum(temp_4d(:, gin, iazi, ipol))
                    end do
                  end do
                end do

                ! Now pull out information needed for chi
                xs % chi_prompt(:, :, :, :) = temp_4d

                ! Deallocate temporary 4D array for nu_fission matrix
                deallocate(temp_4d)

                ! Normalize chi so its CDF goes to 1
                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do gin = 1, energy_groups
                      chi_sum = sum(xs % chi_prompt(:, gin, iazi, ipol))
                      if (chi_sum == ZERO) then
                        call fatal_error("Encountered chi prompt for a group &
                             &that sums to zero")
                      else
                        xs % chi_prompt(:, gin, iazi, ipol) = &
                             xs % chi_prompt(:, gin, iazi, ipol) / chi_sum
                      end if
                    end do
                  end do
                end do
              else
                call fatal_error("prompt-nu-fission must be provided as a 3D &
                     &or 4D array")
              end if
            end if

            ! If delayed-nu-fission provided, set delayed-nu-fission. If
            ! delayed-nu-fission is a matrix, set chi-delayed.
            if (object_exists(xsdata_grp, "delayed-nu-fission")) then

              ! Get the dimensions of the delayed-nu-fission dataset
              xsdata = open_dataset(xsdata_grp, "delayed-nu-fission")
              call get_ndims(xsdata, ndims)

              ! delayed-nu-fission is input as (energy_groups, n_azi, n_pol)
              if (ndims == 3) then

                ! If beta is zeros, raise error
                if (temp_beta(1,1,1,1) == ZERO) then
                  call fatal_error("cannot set delayed-nu-fission with a 3D &
                       &array if beta not provided")
                end if

                ! Allocate temporary arrays for delayed-nu-fission
                allocate(temp_1d(energy_groups * this % n_azi * this % n_pol))
                allocate(temp_3d(energy_groups, this % n_azi, this % n_pol))

                ! Get delayed-nu-fission
                call read_dataset(temp_1d, xsdata_grp, "delayed-nu-fission")
                temp_3d = reshape(temp_1d, (/energy_groups, this % n_azi, &
                     this % n_pol/))

                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do gin = 1, energy_groups
                      do dg = 1, delayed_groups

                        ! Set delayed-nu-fission using delayed neutron fraction
                        xs % delayed_nu_fission(dg, gin, iazi, ipol) = &
                             temp_beta(dg, gin, iazi, ipol) * &
                             temp_3d(gin, iazi, ipol)
                      end do
                    end do
                  end do
                end do

                ! Deallocate temporary delayed-nu-fission arrays
                deallocate(temp_1d)
                deallocate(temp_3d)

                ! If delayed-nu-fission is a (delayed_group, energy_group,
                ! n_azi, n_pol) matrix, set delayed-nu-fission separately for
                ! each delayed group.
              else if (ndims == 4) then

                ! Get the shape of delayed-nu-fission
                call get_shape(xsdata, dims)

                ! Issue error if 1st dimension not correct
                if (dims(1) /= delayed_groups) then
                  call fatal_error("The delayed-nu-fission matrix was input &
                       &with a 1st dimension not equal to the number of &
                       &delayed groups.")
                end if

                ! Issue error if 2nd dimension not correct
                if (dims(2) /= energy_groups) then
                  call fatal_error("The delayed-nu-fission matrix was input &
                       &with a 2nd dimension not equal to the number of &
                       &energy groups.")
                end if

                ! Issue warning if delayed_groups == energy_groups
                if (delayed_groups == energy_groups) then
                  call warning("delayed-nu-fission was input as a dimension &
                       &4 matrix with the same number of delayed groups and &
                       &groups. It is important to know that OpenMC assumes &
                       &the dimensions in the matrix are (delayed_groups, &
                       &energy_groups, n_azi, n_pol). Currently, &
                       &delayed-nu-fission cannot be set as a group by group &
                       &matrix.")
                end if

                ! Get delayed-nu-fission
                allocate(temp_1d(delayed_groups * energy_groups * this % n_azi &
                     * this % n_pol))
                call read_dataset(temp_1d, xsdata_grp, "delayed-nu-fission")
                xs % delayed_nu_fission = reshape(temp_1d, (/delayed_groups, &
                     energy_groups, this % n_azi, this % n_pol /))

                ! Deallocate temporary array for delayed-nu-fission matrix
                deallocate(temp_1d)

                ! If delayed nu-fission is a 5D matrix, set delayed_nu_fission
                ! and chi_delayed.
              else if (ndims == 5) then

                ! chi_delayed is embedded in delayed_nu_fission -> extract
                ! chi_delayed
                allocate(temp_1d(delayed_groups * energy_groups * &
                     energy_groups * this % n_azi * this % n_pol))
                allocate(temp_5d(delayed_groups, energy_groups, energy_groups, &
                     this % n_azi, this % n_pol))
                call read_dataset(temp_1d, xsdata_grp, "delayed-nu-fission")
                temp_5d = reshape(temp_1d, (/delayed_groups, energy_groups, &
                     energy_groups, this % n_azi, this % n_pol/))

                ! Deallocate temporary 1D array for delayed_nu_fission matrix
                deallocate(temp_1d)

                ! Set the 4D delayed-nu-fission matrix and 5D chi_delayed matrix
                ! from the 5D delayed-nu-fission matrix
                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do dg = 1, delayed_groups
                      do gin = 1, energy_groups
                        xs % delayed_nu_fission(dg, gin, iazi, ipol) = &
                             sum(temp_5d(dg, :, gin, iazi, ipol))
                        do gout = 1, energy_groups
                          xs % chi_delayed(dg, gout, gin, iazi, ipol) = &
                               temp_5d(dg, gout, gin, iazi, ipol)
                        end do
                      end do
                    end do
                  end do
                end do

                ! Normalize chi_delayed so its CDF goes to 1
                do ipol = 1, this % n_pol
                  do iazi = 1, this % n_azi
                    do dg = 1, delayed_groups
                      do gin = 1, energy_groups
                        chi_sum = sum(xs % chi_delayed(dg, :, gin, iazi, ipol))
                        if (chi_sum == ZERO) then
                          call fatal_error("Encountered chi delayed for a group&
                               & that sums to zero")
                        else
                          xs % chi_delayed(dg, :, gin, iazi, ipol) = &
                               xs % chi_delayed(dg, :, gin, iazi, ipol) / &
                               chi_sum
                        end if
                      end do
                    end do
                  end do
                end do

                ! Deallocate temporary 5D matrix for delayed_nu_fission
                deallocate(temp_5d)
              else
                call fatal_error("delayed-nu-fission must be provided as a &
                     &3D, 4D, or 5D array")
              end if
            end if

            ! Deallocate temporary beta array
            deallocate(temp_beta)

            ! chi-prompt, chi-delayed, prompt-nu-fission, and delayed-nu-fission
            ! have been set; Now we will check for the rest of the XS that are
            ! unique to fissionable isotopes

            ! Set fission xs
            if (object_exists(xsdata_grp, "fission")) then

              ! Allocate temporary array for fission
              allocate(temp_1d(energy_groups * this % n_azi * this % n_pol))

              ! Get fission array
              call read_dataset(temp_1d, xsdata_grp, "fission")
              xs % fission(:, :, :) = reshape(temp_1d, (/energy_groups, &
                   this % n_azi, this % n_pol/))

              ! Deallocate temporary array for fission
              deallocate(temp_1d)
            end if

            ! Set kappa-fission xs
            if (object_exists(xsdata_grp, "kappa-fission")) then

              ! Allocate temporary array for kappa-fission
              allocate(temp_1d(energy_groups * this % n_azi * this % n_pol))

              ! Get kappa-fission array
              call read_dataset(temp_1d, xsdata_grp, "kappa-fission")
              xs % kappa_fission(:, :, :) = reshape(temp_1d, (/energy_groups, &
                   this % n_azi, this % n_pol/))

              ! Deallocate temporary array for kappa-fission
              deallocate(temp_1d)
            end if

            ! Set decay rate
            if (object_exists(xsdata_grp, "decay rate")) then

              ! Allocate temporary array for decay rate
              allocate(temp_1d(this % n_azi * this % n_pol * delayed_groups))

              ! Get decay rate array
              call read_dataset(temp_1d, xsdata_grp, "decay rate")
              xs % decay_rate(:, :, :) = reshape(temp_1d, (/delayed_groups, &
                   this % n_azi, this % n_pol/))

              ! Deallocate temporary array for decay rate
              deallocate(temp_1d)
            end if
          end if

          ! All the XS unique to fissionable isotopes have been set; Now set all
          ! the generation XS

          if (object_exists(xsdata_grp, "absorption")) then

            ! Allocate temporary array for absorption xs
            allocate(temp_1d(energy_groups * this % n_azi * this % n_pol))

            ! Read in absorption xs
            call read_dataset(temp_1d, xsdata_grp, "absorption")

            xs % absorption = reshape(temp_1d, (/energy_groups, this % n_azi, &
                 this % n_pol/))

            ! Deallocate temporary array for absorption xs
            deallocate(temp_1d)
          else
            call fatal_error("Must provide absorption!")
          end if

          if (object_exists(xsdata_grp, "inverse-velocity")) then

            ! Allocate temporary array for inverse velocity
            allocate(temp_1d(energy_groups * this % n_azi * this % n_pol))

            ! Read in inverse velocity
            call read_dataset(temp_1d, xsdata_grp, "inverse-velocity")

            xs % inverse_velocity = reshape(temp_1d, (/energy_groups, &
                 this % n_azi, this % n_pol/))

            ! Deallocate temporary array for inverse velocity
            deallocate(temp_1d)
          end if

          ! Get scattering data
          if (.not. object_exists(xsdata_grp, "scatter_data")) &
               call fatal_error("Must provide 'scatter_data'")

          scatt_grp = open_group(xsdata_grp, 'scatter_data')

          ! First get the outgoing group boundary indices
          if (object_exists(scatt_grp, "g_min")) then

            allocate(int_arr(energy_groups * this % n_azi * this % n_pol))

            call read_dataset(int_arr, scatt_grp, "g_min")
            allocate(gmin(energy_groups, this % n_azi, this % n_pol))
            gmin = reshape(int_arr, (/energy_groups, this % n_azi, &
                 this % n_pol/))

            deallocate(int_arr)
          else
            call fatal_error("'g_min' for the scatter_data must be provided")
          end if

          if (object_exists(scatt_grp, "g_max")) then

            allocate(int_arr(energy_groups * this % n_azi * this % n_pol))

            call read_dataset(int_arr, scatt_grp, "g_max")
            allocate(gmax(energy_groups, this % n_azi, this % n_pol))
            gmax = reshape(int_arr, (/energy_groups, this % n_azi, &
                 this % n_pol/))

            deallocate(int_arr)
          else
            call fatal_error("'g_max' for the scatter_data must be provided")
          end if

          ! Now use this information to find the length of a container array
          ! to hold the flattened data
          length = 0
          do ipol = 1, this % n_pol
            do iazi = 1, this % n_azi
              do gin = 1, energy_groups
                length = length + order_data * (gmax(gin, iazi, ipol) - &
                     gmin(gin, iazi, ipol) + 1)
              end do
            end do
          end do

          ! Allocate flattened array
          allocate(temp_1d(length))

          if (.not. object_exists(scatt_grp, 'scatter_matrix')) &
               call fatal_error("'scatter_matrix' must be provided")
          call read_dataset(temp_1d, scatt_grp, "scatter_matrix")

          ! Compare the number of orders given with the maximum order of the
          ! problem.  Strip off the superfluous orders if needed.
          if (this % scatter_format == ANGLE_LEGENDRE) then
            order = min(order_data - 1, max_order)
            order_dim = order + 1
          else
            order_dim = order_data
          end if

          ! Convert temp_1d to a jagged array ((gin) % data(l, gout)) for
          ! passing to ScattData
          allocate(input_scatt(energy_groups, this % n_azi, this % n_pol))

          index = 1
          do ipol = 1, this % n_pol
            do iazi = 1, this % n_azi
              do gin = 1, energy_groups
                allocate(input_scatt(gin, iazi, ipol) % data(order_dim, &
                     gmin(gin, iazi, ipol):gmax(gin, iazi, ipol)))
                do gout = gmin(gin, iazi, ipol), gmax(gin, iazi, ipol)
                  do l = 1, order_dim
                    input_scatt(gin, iazi, ipol) % data(l, gout) = &
                         temp_1d(index)
                    index = index + 1
                  end do ! gout
                  ! Adjust index for the orders we didnt take
                  index = index + (order_data - order_dim)
                end do ! order
              end do ! gin
            end do ! iazi
          end do ! ipol

          deallocate(temp_1d)

          ! Finally convert the legendre to tabular if needed
          allocate(scatt_coeffs(energy_groups, this % n_azi, this % n_pol))

          if (this % scatter_format == ANGLE_LEGENDRE .and. &
               legendre_to_tabular) then

            this % scatter_format = ANGLE_TABULAR
            order_dim = legendre_to_tabular_points
            order = order_dim
            dmu = TWO / real(order - 1, 8)

            do ipol = 1, this % n_pol
              do iazi = 1, this % n_azi
                do gin = 1, energy_groups
                  allocate(scatt_coeffs(gin, iazi, ipol) % data(&
                       order_dim, &
                       gmin(gin, iazi, ipol):gmax(gin, iazi, ipol)))

                  do gout = gmin(gin, iazi, ipol), gmax(gin, iazi, ipol)
                    norm = ZERO
                    do imu = 1, order_dim
                      if (imu == 1) then
                        mu = -ONE
                      else if (imu == order_dim) then
                        mu = ONE
                      else
                        mu = -ONE + real(imu - 1, 8) * dmu
                      end if

                      scatt_coeffs(gin, iazi, ipol) % data(imu, gout) = &
                           evaluate_legendre(&
                           input_scatt(gin, iazi, ipol) % data(:, gout), mu)

                      ! Ensure positivity of distribution
                      if (scatt_coeffs(gin, iazi, ipol) % data(imu, gout) < ZERO) &
                           scatt_coeffs(gin, iazi, ipol) % data(imu, gout) = ZERO

                      ! And accrue the integral
                      if (imu > 1) then
                        norm = norm + HALF * dmu * &
                             (scatt_coeffs(gin, iazi, ipol) % data(imu - 1, gout) + &
                              scatt_coeffs(gin, iazi, ipol) % data(imu, gout))
                      end if
                    end do ! mu

                    ! Now that we have the integral, lets ensure that the distribution
                    ! is normalized such that it preserves the original scattering xs
                    if (norm > ZERO) then
                      scatt_coeffs(gin, iazi, ipol) % data(:, gout) = &
                           scatt_coeffs(gin, iazi, ipol) % data(:, gout) * &
                           input_scatt(gin, iazi, ipol) % data(1, gout) / &
                           norm
                    end if
                  end do ! gout
                end do ! gin
              end do ! iazi
            end do ! ipol
          else
            ! Sticking with current representation, carry forward
            do ipol = 1, this % n_pol
              do iazi = 1, this % n_azi
                do gin = 1, energy_groups
                  allocate(scatt_coeffs(gin, iazi, ipol) % data(order_dim, &
                       gmin(gin, iazi, ipol):gmax(gin, iazi, ipol)))
                  scatt_coeffs(gin, iazi, ipol) % data(:, :) = &
                       input_scatt(gin, iazi, ipol) % data(1:order_dim, :)
                end do
              end do
            end do
          end if

          deallocate(input_scatt)

          ! Now get the multiplication matrix
          if (object_exists(scatt_grp, 'multiplicity_matrix')) then

            ! Now use this information to find the length of a container array
            ! to hold the flattened data
            length = 0

            do ipol = 1, this % n_pol
              do iazi = 1, this % n_azi
                do gin = 1, energy_groups
                  length = length + (gmax(gin, iazi, ipol) - gmin(gin, iazi, ipol) + 1)
                end do
              end do
            end do

            ! Allocate flattened array
            allocate(temp_1d(length))
            call read_dataset(temp_1d, scatt_grp, "multiplicity_matrix")

            ! Convert temp_1d to a jagged array ((gin) % data(gout)) for passing
            ! to ScattData
            allocate(temp_mult(energy_groups, this % n_azi, this % n_pol))

            index = 1
            do ipol = 1, this % n_pol
              do iazi = 1, this % n_azi
                do gin = 1, energy_groups
                  allocate(temp_mult(gin, iazi, ipol) % data( &
                       gmin(gin, iazi, ipol):gmax(gin, iazi, ipol)))
                  do gout = gmin(gin, iazi, ipol), gmax(gin, iazi, ipol)
                    temp_mult(gin, iazi, ipol) % data(gout) = temp_1d(index)
                    index = index + 1
                  end do
                end do
              end do
            end do
            deallocate(temp_1d)
          else

            allocate(temp_mult(energy_groups, this % n_azi, this % n_pol))

            ! Default to multiplicities of 1.0
            do ipol = 1, this % n_pol
              do iazi = 1, this % n_azi
                do gin = 1, energy_groups
                  allocate(temp_mult(gin, iazi, ipol) % data( &
                       gmin(gin, iazi, ipol):gmax(gin, iazi, ipol)))
                  temp_mult(gin, iazi, ipol) % data = ONE
                end do
              end do
            end do
          end if

          ! Allocate and initialize our ScattData Object.
          allocate(xs % scatter(this % n_azi, this % n_pol))

          do ipol = 1, this % n_pol
            do iazi = 1,  this % n_azi

              ! Allocate and initialize our ScattData Object.
              if (this % scatter_format == ANGLE_HISTOGRAM) then
                allocate(ScattDataHistogram :: xs % scatter(iazi, ipol) % obj)
              else if (this % scatter_format == ANGLE_TABULAR) then
                allocate(ScattDataTabular :: xs % scatter(iazi, ipol) % obj)
              else if (this % scatter_format == ANGLE_LEGENDRE) then
                allocate(ScattDataLegendre :: xs % scatter(iazi, ipol) % obj)
              end if

              ! Initialize the ScattData Object
              call xs % scatter(iazi, ipol) % obj % init(gmin(:, iazi, ipol), &
                   gmax(:, iazi, ipol), temp_mult(:, iazi, ipol), &
                   scatt_coeffs(:, iazi, ipol))
            end do
          end do

          ! Check sigA to ensure it is not 0 since it is
          ! often divided by in the tally routines
          ! (This may happen with Helium data)
          do ipol = 1, this % n_pol
            do iazi = 1,  this % n_azi
              do gin = 1, energy_groups
                if (xs % absorption(gin, iazi, ipol) == ZERO) then
                  xs % absorption(gin, iazi, ipol) = 1E-10_8
                end if
              end do
            end do
          end do

          if (object_exists(xsdata_grp, "total")) then

            allocate(temp_1d(energy_groups * this % n_azi * this % n_pol))
            call read_dataset(temp_1d, xsdata_grp, "total")
            xs % total = reshape(temp_1d, (/energy_groups, this % n_azi, &
                 this % n_pol/))

            deallocate(temp_1d)
          else
            do ipol = 1, this % n_pol
              do iazi = 1, this % n_azi
                xs % total(:, iazi, ipol) = xs % absorption(:, iazi, ipol) + &
                     xs % scatter(iazi, ipol) % obj % scattxs(:)
              end do
            end do
          end if

          ! Check sigT to ensure it is not 0 since it is often divided by in
          ! the tally routines
          do ipol = 1, this % n_pol
            do iazi = 1,  this % n_azi
              do gin = 1, energy_groups
                if (xs % total(gin, iazi, ipol) == ZERO) then
                  xs % total(gin, iazi, ipol) = 1E-10_8
                end if
              end do
            end do
          end do

          ! Close the groups we have opened and deallocate
          call close_group(xsdata_grp)
          call close_group(scatt_grp)
          deallocate(scatt_coeffs, temp_mult)

        end associate ! xs
      end do ! Temperatures
    end subroutine mgxsang_from_hdf5

!===============================================================================
! MGXS*_COMBINE Builds a macroscopic Mgxs object from microscopic Mgxs objects
!===============================================================================

    subroutine mgxs_combine(this, temps, mat, nuclides, max_order, &
                            scatter_format, order_dim)
      class(Mgxs), intent(inout)          :: this ! The Mgxs to initialize
      type(VectorReal), intent(in)        :: temps ! Temperatures to obtain
      type(Material), pointer, intent(in) :: mat  ! base material
      type(MgxsContainer), intent(in)     :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                 :: max_order  ! Maximum requested order
      integer, intent(out)                :: scatter_format ! Type of scatter
      integer, intent(out)                :: order_dim    ! Scattering data order size

      integer :: t, mat_max_order, order

      ! Fill in meta-data from material information
      if (mat % name == "") then
        this % name = trim(to_str(mat % id))
      else
        this % name = trim(mat % name)
      end if

      ! Set whether this material is fissionable
      this % fissionable = mat % fissionable

      ! The following info we should initialize, but we dont need it nor
      ! does it have guaranteed meaning.
      this % awr = -ONE

      allocate(this % kTs(temps % size()))

      do t = 1, temps % size()
        this % kTs(t) = temps % data(t)
      end do

      ! Allocate the XS object for the number of temperatures
      select type(this)
      type is (MgxsIso)
        allocate(this % xs(temps % size()))
      type is (MgxsAngle)
        allocate(this % xs(temps % size()))
      end select

      ! Determine the scattering type of our data and ensure all scattering orders
      ! are the same.
      scatter_format = nuclides(mat % nuclide(1)) % obj % scatter_format

      select type(nuc => nuclides(mat % nuclide(1)) % obj)
      type is (MgxsIso)
        order = size(nuc % xs(1) % scatter % dist(1) % data, dim=1)
      type is (MgxsAngle)
        order = size(nuc % xs(1) % scatter(1, 1) % obj % dist(1) % data, dim=1)
      end select

      ! If we have tabular only data, then make sure all datasets have same size
      if (scatter_format == ANGLE_HISTOGRAM) then
        ! Check all scattering data to ensure it is the same size
        do i = 2, mat % n_nuclides
          select type(nuc => nuclides(mat % nuclide(i)) % obj)
          type is (MgxsIso)
            if (order /= size(nuc % xs(1) % scatter % dist(1) % data, dim=1)) &
                 call fatal_error("All histogram scattering entries must be&
                                  & same length!")
          type is (MgxsAngle)
            if (order /= size(nuc % xs(1) % scatter(1, 1) % obj % dist(1) % data, dim=1)) &
                 call fatal_error("All histogram scattering entries must be&
                                  & same length!")
          end select
        end do

        ! Ok, got our order, store the dimensionality
        order_dim = order

      else if (scatter_format == ANGLE_TABULAR) then
        ! Check all scattering data to ensure it is the same size
        do i = 2, mat % n_nuclides
          select type(nuc => nuclides(mat % nuclide(i)) % obj)
          type is (MgxsIso)
            if (order /= size(nuc % xs(1) % scatter % dist(1) % data, dim=1)) &
                 call fatal_error("All tabular scattering entries must be&
                                  & same length!")
          type is (MgxsAngle)
            if (order /= size(nuc % xs(1) % scatter(1, 1) % obj % dist(1) % data, dim=1)) &
                 call fatal_error("All tabular scattering entries must be&
                                  & same length!")
          end select
        end do

        ! Ok, got our order, store the dimensionality
        order_dim = order

      else if (scatter_format == ANGLE_LEGENDRE) then

        ! Need to determine the maximum scattering order of all data in this material
        mat_max_order = 0

        do i = 1, mat % n_nuclides
          select type(nuc => nuclides(mat % nuclide(i)) % obj)
          type is (MgxsIso)
            if (size(nuc % xs(1) % scatter % dist(1) % data, &
                     dim=1) > mat_max_order) &
                 mat_max_order = size(nuc % xs(1) % scatter % dist(1) % data, &
                                      dim=1)
          type is (MgxsAngle)
            if (size(nuc % xs(1) % scatter(1, 1) % obj % dist(1) % data, &
                     dim=1) > mat_max_order) &
                 mat_max_order = &
                     size(nuc % xs(1) % scatter(1, 1) % obj % dist(1) % data, &
                          dim=1)
          end select
        end do

        ! Now need to compare this material maximum scattering order with
        ! the problem wide max scatt order and use whichever is lower
        order = min(mat_max_order, max_order + 1)

        ! Ok, got our order, store the dimensionality
        order_dim = order
      end if

    end subroutine mgxs_combine

    subroutine mgxsiso_combine(this, temps, mat, nuclides, energy_groups, &
         delayed_groups, max_order, tolerance, method)
      class(MgxsIso), intent(inout)       :: this  ! The Mgxs to initialize
      type(VectorReal), intent(in)        :: temps ! Temperatures to obtain [MeV]
      type(Material), pointer, intent(in) :: mat   ! base material
      type(MgxsContainer), intent(in)     :: nuclides(:)    ! List of nuclides to harvest from
      integer, intent(in)                 :: energy_groups  ! Number of energy groups
      integer, intent(in)                 :: delayed_groups ! Number of delayed groups
      integer, intent(in)                 :: max_order  ! Maximum requested order
      real(8), intent(in)                 :: tolerance  ! Tolerance on method
      integer, intent(in)                 :: method     ! Type of temperature access

      integer :: i            ! loop index over nuclides
      integer :: t            ! Index in to temps
      integer :: gin, gout    ! group indices
      integer :: dg           ! delayed group index
      real(8) :: atom_density ! atom density of a nuclide
      real(8) :: norm, nuscatt
      integer :: order_dim, nuc_order_dim
      real(8), allocatable :: temp_mult(:, :), mult_num(:, :), mult_denom(:, :)
      real(8), allocatable :: scatt_coeffs(:, :, :)
      integer :: nuc_t
      integer, allocatable :: nuc_ts(:)
      real(8) :: temp_actual, temp_desired, interp
      integer :: scatter_format
      type(Jagged2D), allocatable :: nuc_matrix(:)
      integer, allocatable :: gmin(:), gmax(:)
      type(Jagged2D), allocatable :: jagged_scatt(:)
      type(Jagged1D), allocatable :: jagged_mult(:)

      ! Set the meta-data
      call mgxs_combine(this, temps, mat, nuclides, max_order, scatter_format, &
                        order_dim)

      ! Set the number of delayed groups
      this % num_delayed_groups = delayed_groups

      ! Create the Xs Data for each temperature
      TEMP_LOOP: do t = 1, temps % size()

        ! Allocate and initialize the data needed for macro_xs(i_mat) object
        allocate(this % xs(t) % total(energy_groups))
        this % xs(t) % total(:) = ZERO

        allocate(this % xs(t) % absorption(energy_groups))
        this % xs(t) % absorption(:) = ZERO

        allocate(this % xs(t) % fission(energy_groups))
        this % xs(t) % fission(:) = ZERO

        allocate(this % xs(t) % kappa_fission(energy_groups))
        this % xs(t) % kappa_fission(:) = ZERO

        allocate(this % xs(t) % prompt_nu_fission(energy_groups))
        this % xs(t) % prompt_nu_fission(:) = ZERO

        allocate(this % xs(t) % delayed_nu_fission(delayed_groups, &
             energy_groups))
        this % xs(t) % delayed_nu_fission(:, :) = ZERO

        allocate(this % xs(t) % chi_prompt(energy_groups, energy_groups))
        this % xs(t) % chi_prompt(:, :) = ZERO

        allocate(this % xs(t) % chi_delayed(delayed_groups, energy_groups, &
             energy_groups))
        this % xs(t) % chi_delayed(:, :, :) = ZERO

        allocate(this % xs(t) % inverse_velocity(energy_groups))
        this % xs(t) % inverse_velocity(:) = ZERO

        allocate(this % xs(t) % decay_rate(delayed_groups))
        this % xs(t) % decay_rate(:) = ZERO

        allocate(temp_mult(energy_groups, energy_groups))
        temp_mult(:, :) = ZERO

        allocate(mult_num(energy_groups, energy_groups))
        mult_num(:, :) = ZERO

        allocate(mult_denom(energy_groups, energy_groups))
        mult_denom(:, :) = ZERO

        allocate(scatt_coeffs(order_dim, energy_groups, energy_groups))
        scatt_coeffs(:, :, :) = ZERO

        this % scatter_format = scatter_format

        if (scatter_format == ANGLE_LEGENDRE) then
          allocate(ScattDataLegendre :: this % xs(t) % scatter)
        else if (scatter_format == ANGLE_TABULAR) then
          allocate(ScattDataTabular :: this % xs(t) % scatter)
        else if (scatter_format == ANGLE_HISTOGRAM) then
          allocate(ScattDataHistogram :: this % xs(t) % scatter)
        end if

        ! Add contribution from each nuclide in material
        NUC_LOOP: do i = 1, mat % n_nuclides
          associate(nuc => nuclides(mat % nuclide(i)) % obj)

            ! Copy atom density of nuclide in material
            atom_density = mat % atom_density(i)

            select case (method)
            case (TEMPERATURE_NEAREST)

              ! Determine actual temperatures to read
              temp_desired = temps % data(i)
              allocate(nuc_ts(1))

              nuc_ts(1) = minloc(abs(nuc % kTs - temp_desired), dim=1)
              temp_actual = nuc % kTs(nuc_ts(1))

              if (abs(temp_actual - temp_desired) >= K_BOLTZMANN * tolerance) then
                call fatal_error("MGXS library does not contain cross sections &
                     &for " // trim(this % name) // " at or near " // &
                     trim(to_str(nint(temp_desired / K_BOLTZMANN))) // " K.")
              end if

            case (TEMPERATURE_INTERPOLATION)

              ! If temperature interpolation or multipole is selected, get a
              ! list of bounding temperatures for each actual temperature
              ! present in the model
              temp_desired = temps % data(i)
              allocate(nuc_ts(2))

              do j = 1, size(nuc % kTs) - 1
                if (nuc % kTs(j) <= temp_desired .and. &
                     temp_desired < nuc % kTs(j + 1)) then
                  nuc_ts(1) = j
                  nuc_ts(2) = j + 1
                end if
              end do

              call fatal_error("Nuclear data library does not contain cross sections &
                   &for " // trim(this % name) // " at temperatures that bound " // &
                   trim(to_str(nint(temp_desired / K_BOLTZMANN))) // " K.")
            end select

            select type(nuc)
            type is (MgxsIso)
              do j = 1, size(nuc_ts)

                nuc_t = nuc_ts(j)

                if (size(nuc_ts) == 1) then
                  interp = ONE
                else if (j == 1) then
                  interp = (ONE - (temp_desired - nuc % kTs(nuc_ts(1))) / &
                       (nuc % kTs(nuc_ts(2)) - nuc % kTs(nuc_ts(1))))
                else
                  interp = ONE - interp
                end if

                ! Perform our operations which depend upon the type
                ! Add contributions to total, absorption, and fission data (if necessary)
                this % xs(t) % total = this % xs(t) % total + &
                     atom_density * nuc % xs(nuc_t) % total * interp

                this % xs(t) % absorption = this % xs(t) % absorption + &
                     atom_density * nuc % xs(nuc_t) % absorption * interp

                this % xs(t) % decay_rate = this % xs(t) % decay_rate + &
                     atom_density * nuc % xs(nuc_t) % decay_rate * interp

                this % xs(t) % inverse_velocity = &
                     this % xs(t) % inverse_velocity + &
                     atom_density * nuc % xs(nuc_t) % inverse_velocity * interp

                if (nuc % fissionable) then

                  this % xs(t) % chi_prompt = this % xs(t) % chi_prompt + &
                       atom_density * nuc % xs(nuc_t) % chi_prompt * interp

                  this % xs(t) % chi_delayed = this % xs(t) % chi_delayed + &
                       atom_density * nuc % xs(nuc_t) % chi_delayed * interp

                  this % xs(t) % prompt_nu_fission = this % xs(t) % &
                       prompt_nu_fission + atom_density * nuc % xs(nuc_t) % &
                       prompt_nu_fission * interp

                  this % xs(t) % delayed_nu_fission = this % xs(t) % &
                       delayed_nu_fission + atom_density * nuc % xs(nuc_t) % &
                       delayed_nu_fission * interp

                  this % xs(t) % fission = this % xs(t) % fission + &
                       atom_density * nuc % xs(nuc_t) % fission * interp

                  this % xs(t) % kappa_fission = this % xs(t) % kappa_fission +&
                       atom_density * nuc % xs(nuc_t) % kappa_fission * interp
                end if

                ! We will next gather the multiplicity and scattering matrices.
                ! To avoid multiple re-allocations as we resize the storage
                ! matrix (and/or to avoidlots of duplicate code), we will use a
                ! dense matrix for this storage, with a reduction to the sparse
                ! format at the end.

                ! Get the multiplicity_matrix
                ! To combine from nuclidic data we need to use the final relationship
                ! mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) /
                !              sum_i(N_i*(nuscatt_{i,g,g'} / mult_{i,g,g'}))
                ! Developed as follows:
                ! mult_{gg'} = nuScatt{g,g'} / Scatt{g,g'}
                ! mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) / sum(N_i*scatt_{i,g,g'})
                ! mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) /
                !              sum_i(N_i*(nuscatt_{i,g,g'} / mult_{i,g,g'}))
                ! nuscatt_{i,g,g'} can be reconstructed from scatter % energy and
                ! scatter % scattxs
                do gin = 1, energy_groups
                  do gout = nuc % xs(nuc_t) % scatter % gmin(gin), nuc % xs(nuc_t) % scatter % gmax(gin)

                    nuscatt = nuc % xs(nuc_t) % scatter % scattxs(gin) * &
                         nuc % xs(nuc_t) % scatter % energy(gin) % data(gout)

                    mult_num(gout, gin) = mult_num(gout, gin) + atom_density * &
                         nuscatt * interp

                    if (nuc % xs(nuc_t) % scatter % mult(gin) % data(gout) > ZERO) then
                      mult_denom(gout, gin) = mult_denom(gout,gin) + atom_density * &
                           nuscatt / nuc % xs(nuc_t) % scatter % mult(gin) % data(gout) * &
                           interp
                    else
                      ! Avoid division by zero
                      mult_denom(gout, gin) = mult_denom(gout,gin) + atom_density * &
                           interp
                    end if
                  end do
                end do

                ! Get the complete scattering matrix
                nuc_order_dim = size(nuc % xs(nuc_t) % scatter % dist(1) % data, dim=1)
                nuc_order_dim = min(nuc_order_dim, order_dim)

                call nuc % xs(nuc_t) % scatter % get_matrix(nuc_order_dim, &
                     nuc_matrix)

                do gin = 1, energy_groups
                  do gout = nuc % xs(nuc_t) % scatter % gmin(gin), &
                       nuc % xs(nuc_t) % scatter % gmax(gin)
                    scatt_coeffs(1:nuc_order_dim, gout, gin) = &
                         scatt_coeffs(1: nuc_order_dim, gout, gin) + &
                         atom_density * interp * &
                         nuc_matrix(gin) % data(1:nuc_order_dim, gout)
                  end do
                end do
              end do

            type is (MgxsAngle)
              call fatal_error("Invalid passing of MgxsAngle to MgxsIso object")
            end select

            ! Obtain temp_mult
            do gin = 1, energy_groups
              do gout = 1, energy_groups
                if (mult_denom(gout, gin) > ZERO) then
                  temp_mult(gout, gin) = mult_num(gout, gin) / mult_denom(gout, gin)
                else
                  temp_mult(gout, gin) = ONE
                end if
              end do
            end do

            ! Now create our jagged data from the dense data
            call jagged_from_dense_2D(scatt_coeffs, jagged_scatt, gmin, gmax)
            call jagged_from_dense_1D(temp_mult, jagged_mult)

            ! Initialize the ScattData Object
            call this % xs(t) % scatter % init(gmin, gmax, jagged_mult, &
                 jagged_scatt)

            ! Now normalize chi
            if (mat % fissionable) then
              do gin = 1, energy_groups
                norm =  sum(this % xs(t) % chi_prompt(:, gin))
                if (norm > ZERO) then
                  this % xs(t) % chi_prompt(:, gin) = &
                       this % xs(t) % chi_prompt(:, gin) / norm
                end if
              end do

              do dg = 1, delayed_groups
                do gin = 1, energy_groups
                  norm =  sum(this % xs(t) % chi_delayed(dg, :, gin))
                  if (norm > ZERO) then
                    this % xs(t) % chi_delayed(dg, :, gin) = &
                         this % xs(t) % chi_delayed(dg, :, gin) / norm
                  end if
                end do
              end do
            end if

            ! Deallocate temporaries
            deallocate(jagged_mult, jagged_scatt, gmin, gmax, scatt_coeffs, &
                       temp_mult, mult_num, mult_denom)
          end associate ! nuc
        end do NUC_LOOP
      end do TEMP_LOOP

    end subroutine mgxsiso_combine

    subroutine mgxsang_combine(this, temps, mat, nuclides, energy_groups, &
         delayed_groups, max_order, tolerance, method)
      class(MgxsAngle), intent(inout)     :: this ! The Mgxs to initialize
      type(VectorReal), intent(in)        :: temps ! Temperatures to obtain
      type(Material), pointer, intent(in) :: mat  ! base material
      type(MgxsContainer), intent(in)     :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                 :: energy_groups  ! Number of energy groups
      integer, intent(in)                 :: delayed_groups ! Number of delayed groups
      integer, intent(in)                 :: max_order  ! Maximum requested order
      real(8), intent(in)                 :: tolerance  ! Tolerance on method
      integer, intent(in)                 :: method     ! Type of temperature access

      integer :: i             ! loop index over nuclides
      integer :: t             ! temperature loop index
      integer :: gin, gout     ! group indices
      integer :: dg            ! delayed group index
      real(8) :: atom_density  ! atom density of a nuclide
      integer :: ipol, iazi, n_pol, n_azi
      real(8) :: norm, nuscatt
      integer :: order_dim, nuc_order_dim
      real(8), allocatable :: temp_mult(:, :, :, :), mult_num(:, :, :, :)
      real(8), allocatable :: mult_denom(:, :, :, :), scatt_coeffs(:, :, :, :, :)
      integer :: nuc_t
      integer, allocatable :: nuc_ts(:)
      real(8) :: temp_actual, temp_desired, interp
      integer :: scatter_format
      type(Jagged2D), allocatable :: nuc_matrix(:)
      integer, allocatable :: gmin(:), gmax(:)
      type(Jagged2D), allocatable :: jagged_scatt(:)
      type(Jagged1D), allocatable :: jagged_mult(:)

      ! Set the meta-data
      call mgxs_combine(this, temps, mat, nuclides, max_order, scatter_format, &
           order_dim)

      ! Set the number of delayed groups
      this % num_delayed_groups = delayed_groups

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
            this % polar(:) = nuc % polar(:)

            allocate(this % azimuthal(n_azi))
            this % azimuthal(:) = nuc % azimuthal(:)
          else
            if ((n_pol /= nuc % n_pol) .or. (n_azi /= nuc % n_azi)) then
              call fatal_error("All angular data must be same length!")
            end if
          end if
        end select
      end do

      ! Create the Xs Data for each temperature
      TEMP_LOOP: do t = 1, temps % size()

        ! Allocate and initialize the data needed for macro_xs(i_mat) object
        allocate(this % xs(t) % total(energy_groups, n_azi, n_pol))
        this % xs(t) % total = ZERO

        allocate(this % xs(t) % absorption(energy_groups, n_azi, n_pol))
        this % xs(t) % absorption = ZERO

        allocate(this % xs(t) % fission(energy_groups, n_azi, n_pol))
        this % xs(t) % fission = ZERO

        allocate(this % xs(t) % decay_rate(delayed_groups, n_azi, n_pol))
        this % xs(t) % decay_rate = ZERO

        allocate(this % xs(t) % inverse_velocity(energy_groups, n_azi, n_pol))
        this % xs(t) % inverse_velocity = ZERO

        allocate(this % xs(t) % kappa_fission(energy_groups, n_azi, n_pol))
        this % xs(t) % kappa_fission = ZERO

        allocate(this % xs(t) % prompt_nu_fission(energy_groups, n_azi, n_pol))
        this % xs(t) % prompt_nu_fission = ZERO

        allocate(this % xs(t) % delayed_nu_fission(delayed_groups, &
             energy_groups, n_azi, n_pol))
        this % xs(t) % delayed_nu_fission = ZERO

        allocate(this % xs(t) % chi_prompt(energy_groups, energy_groups, &
             n_azi, n_pol))
        this % xs(t) % chi_prompt = ZERO

        allocate(this % xs(t) % chi_delayed(delayed_groups, energy_groups, &
             energy_groups, n_azi, n_pol))
        this % xs(t) % chi_delayed = ZERO

        allocate(temp_mult(energy_groups, energy_groups, n_azi, n_pol))
        temp_mult = ZERO

        allocate(mult_num(energy_groups, energy_groups, n_azi, n_pol))
        mult_num = ZERO

        allocate(mult_denom(energy_groups, energy_groups, n_azi, n_pol))
        mult_denom = ZERO

        allocate(scatt_coeffs(order_dim, energy_groups, energy_groups, n_azi, n_pol))
        scatt_coeffs = ZERO

        allocate(this % xs(t) % scatter(n_azi, n_pol))

        do ipol = 1, n_pol
          do iazi = 1, n_azi
            if (scatter_format == ANGLE_LEGENDRE) then
              allocate(ScattDataLegendre :: &
                       this % xs(t) % scatter(iazi, ipol) % obj)
            else if (scatter_format == ANGLE_TABULAR) then
              allocate(ScattDataTabular :: &
                       this % xs(t) % scatter(iazi, ipol) % obj)
            else if (scatter_format == ANGLE_HISTOGRAM) then
              allocate(ScattDataHistogram :: &
                       this % xs(t) % scatter(iazi, ipol) % obj)
            end if
          end do
        end do

        ! Add contribution from each nuclide in material
        NUC_LOOP: do i = 1, mat % n_nuclides
          associate(nuc => nuclides(mat % nuclide(i)) % obj)

          select case (method)
          case (TEMPERATURE_NEAREST)

            ! Determine actual temperatures to read
            temp_desired = temps % data(i)
            allocate(nuc_ts(1))

            nuc_ts(1) = minloc(abs(nuc % kTs - temp_desired), dim=1)
            temp_actual = nuc % kTs(nuc_ts(1))

            if (abs(temp_actual - temp_desired) >= K_BOLTZMANN * tolerance) then
              call fatal_error("MGXS library does not contain cross sections &
                   &for " // trim(this % name) // " at or near " // &
                   trim(to_str(nint(temp_desired / K_BOLTZMANN))) // " K.")
            end if

          case (TEMPERATURE_INTERPOLATION)

            ! If temperature interpolation or multipole is selected, get a
            ! list of bounding temperatures for each actual temperature
            ! present in the model
            temp_desired = temps % data(i)
            allocate(nuc_ts(2))

            do j = 1, size(nuc % kTs) - 1
              if (nuc % kTs(j) <= temp_desired .and. &
                   temp_desired < nuc % kTs(j + 1)) then
                nuc_ts(1) = j
                nuc_ts(2) = j + 1
              end if
            end do

            call fatal_error("Nuclear data library does not contain cross sections &
                 &for " // trim(this % name) // " at temperatures that bound " // &
                 trim(to_str(nint(temp_desired / K_BOLTZMANN))) // " K.")
          end select

             ! Copy atom density of nuclide in material
            atom_density = mat % atom_density(i)

            select type(nuc)
            type is (MgxsAngle)
              do j = 1, size(nuc_ts)

                nuc_t = nuc_ts(j)

                if (size(nuc_ts) == 1) then
                  interp = ONE
                else if (j == 1) then
                  interp = (ONE - (temp_desired - nuc % kTs(nuc_ts(1))) / &
                       (nuc % kTs(nuc_ts(2)) - nuc % kTs(nuc_ts(1))))
                else
                  interp = ONE - interp
                end if

                ! Perform our operations which depend upon the type
                ! Add contributions to total, absorption, and fission data
                ! (if necessary)
                this % xs(t) % total = this % xs(t) % total + &
                     atom_density * nuc % xs(nuc_t) % total * interp

                this % xs(t) % absorption = this % xs(t) % absorption + &
                     atom_density * nuc % xs(nuc_t) % absorption * interp

                this % xs(t) % decay_rate = this % xs(t) % decay_rate + &
                     atom_density * nuc % xs(nuc_t) % decay_rate * interp

                this % xs(t) % inverse_velocity = &
                     this % xs(t) % inverse_velocity + &
                     atom_density * nuc % xs(nuc_t) % inverse_velocity * interp

                if (nuc % fissionable) then

                  this % xs(t) % chi_prompt = this % xs(t) % chi_prompt + &
                       atom_density * nuc % xs(nuc_t) % chi_prompt * interp

                  this % xs(t) % chi_delayed = this % xs(t) % chi_delayed + &
                       atom_density * nuc % xs(nuc_t) % chi_delayed * interp

                  this % xs(t) % prompt_nu_fission = &
                       this % xs(t) % prompt_nu_fission + atom_density * &
                       nuc % xs(nuc_t) % prompt_nu_fission * interp

                  this % xs(t) % delayed_nu_fission = &
                       this % xs(t) % delayed_nu_fission + atom_density * &
                       nuc % xs(nuc_t) % delayed_nu_fission * interp

                  this % xs(t) % fission = this % xs(t) % fission + &
                       atom_density * nuc % xs(nuc_t) % fission * interp

                  this % xs(t) % kappa_fission = this % xs(t) % kappa_fission &
                       + atom_density * nuc % xs(nuc_t) % kappa_fission * interp

                end if

                ! We will next gather the multiplicity and scattering matrices.
                ! To avoid multiple re-allocations as we resize the storage
                ! matrix (and/or to avoidlots of duplicate code), we will use a
                ! dense matrix for this storage, with a reduction to the sparse
                ! format at the end.

                ! Get the multiplicity_matrix
                ! To combine from nuclidic data we need to use the final relationship
                ! mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) /
                !              sum_i(N_i*(nuscatt_{i,g,g'} / mult_{i,g,g'}))
                ! Developed as follows:
                ! mult_{gg'} = nuScatt{g,g'} / Scatt{g,g'}
                ! mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) / sum(N_i*scatt_{i,g,g'})
                ! mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) /
                !              sum_i(N_i*(nuscatt_{i,g,g'} / mult_{i,g,g'}))
                ! nuscatt_{i,g,g'} can be reconstructed from scatter % energy and
                ! scatter % scattxs
                do ipol = 1, n_pol
                  do iazi = 1, n_azi
                    do gin = 1, energy_groups
                      do gout = nuc % xs(nuc_t) % scatter(iazi, ipol) %obj % gmin(gin), &
                             nuc % xs(nuc_t) % scatter(iazi, ipol) % obj % gmax(gin)

                        nuscatt = nuc % xs(nuc_t) % scatter(iazi, ipol) % obj % scattxs(gin) * &
                             nuc % xs(nuc_t) % scatter(iazi, ipol) % obj % energy(gin) % data(gout)

                        mult_num(gout, gin, iazi, ipol) = &
                             mult_num(gout, gin, iazi, ipol) + atom_density * &
                             nuscatt * interp

                        if (nuc % xs(nuc_t) % scatter(iazi, ipol) % obj % mult(gin) % data(gout) > ZERO) then
                          mult_denom(gout, gin, iazi, ipol) = &
                               mult_denom(gout, gin, iazi, ipol) + atom_density * &
                               interp * nuscatt / &
                               nuc % xs(nuc_t) % scatter(iazi, ipol) % obj % mult(gin) % data(gout)
                        else
                          ! Avoid division by zero
                          mult_denom(gout, gin, iazi, ipol) = &
                               mult_denom(gout, gin, iazi, ipol) + atom_density * &
                               interp
                        end if
                      end do
                    end do

                    ! Get the complete scattering matrix
                    nuc_order_dim = &
                         size(nuc % xs(nuc_t) % scatter(iazi, ipol) % obj % &
                         dist(1) % data, dim=1)
                    nuc_order_dim = min(nuc_order_dim, order_dim)

                    call nuc % xs(nuc_t) % scatter(iazi, ipol) % obj % &
                         get_matrix(nuc_order_dim, nuc_matrix)

                    do gin = 1, energy_groups
                      do gout = &
                           nuc % xs(nuc_t) % scatter(iazi, ipol) % obj % gmin(gin), &
                           nuc % xs(nuc_t) % scatter(iazi, ipol) % obj % gmax(gin)
                        scatt_coeffs(1:nuc_order_dim, gout, gin, iazi, ipol) = &
                             scatt_coeffs(1: nuc_order_dim, gout, gin, iazi, ipol) + &
                             atom_density * interp * &
                             nuc_matrix(gin) % data(1:nuc_order_dim, gout)
                      end do ! gout
                    end do ! gin
                  end do ! iazi
                end do ! ipol
              end do
            type is (MgxsIso)
              call fatal_error("Invalid passing of MgxsIso to MgxsAngle object")
            end select

            ! Obtain temp_mult, create jaged arrays and initialize the
            ! ScattData object.
            do ipol = 1, n_pol
              do iazi = 1, n_azi
                ! Obtain temp_mult
                do gin = 1, energy_groups
                  do gout = 1, energy_groups
                    if (mult_denom(gout, gin, iazi, ipol) > ZERO) then
                      temp_mult(gout, gin, iazi, ipol) = &
                           mult_num(gout, gin, iazi, ipol) / &
                           mult_denom(gout, gin, iazi, ipol)
                    else
                      temp_mult(gout, gin, iazi, ipol) = ONE
                    end if
                  end do
                end do

                ! Now create our jagged data from the dense data
                call jagged_from_dense_2D(scatt_coeffs(:, :, :, iazi, ipol), &
                                          jagged_scatt, gmin, gmax)
                call jagged_from_dense_1D(temp_mult(:, :, iazi, ipol), &
                                          jagged_mult)

                ! Initialize the ScattData Object
                call this % xs(t) % scatter(iazi, ipol) % obj % init(gmin, &
                     gmax, jagged_mult, jagged_scatt)
                deallocate(jagged_scatt, jagged_mult, gmin, gmax)
              end do
            end do

            ! Now normalize chi
            if (mat % fissionable) then
              do ipol = 1, n_pol
                do iazi = 1, n_azi
                  do gin = 1, energy_groups
                    norm =  sum(this % xs(t) % chi_prompt(:, gin, iazi, ipol))
                    if (norm > ZERO) then
                      this % xs(t) % chi_prompt(:, gin, iazi, ipol) = &
                           this % xs(t) % chi_prompt(:, gin, iazi, ipol) / norm
                    end if
                  end do
                end do
              end do

              do dg = 1, delayed_groups
                do ipol = 1, n_pol
                  do iazi = 1, n_azi
                    do gin = 1, energy_groups
                      norm =  sum(this % xs(t) % chi_delayed(dg, :, gin, iazi, ipol))
                      if (norm > ZERO) then
                        this % xs(t) % chi_delayed(dg, :, gin, iazi, ipol) = &
                             this % xs(t) % chi_delayed(dg, :, gin, iazi, ipol)&
                             / norm
                      end if
                    end do
                  end do
                end do
              end do

            end if

            ! Deallocate temporaries
            deallocate(scatt_coeffs, temp_mult, mult_num, mult_denom)
          end associate ! nuc
        end do NUC_LOOP

      end do TEMP_LOOP

    end subroutine mgxsang_combine

!===============================================================================
! MGXS*_GET_XS returns the requested data cross section data
!===============================================================================

    pure function mgxsiso_get_xs(this, xstype, gin, gout, uvw, mu, dg) result(xs)
      class(MgxsIso), intent(in)    :: this   ! The Xs to get data from
      character(*) , intent(in)     :: xstype ! Type of xs requested
      integer, intent(in)           :: gin    ! Incoming Energy group
      integer, optional, intent(in) :: gout   ! Outgoing Energy group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      integer, optional, intent(in) :: dg     ! Delayed group
      real(8)                       :: xs ! Requested x/s
      integer                       :: t ! temperature index

      t = this % index_temp

      select case(xstype)

      case('total')
        xs = this % xs(t) % total(gin)

      case('absorption')
        xs = this % xs(t) % absorption(gin)

      case('fission')
        xs = this % xs(t) % fission(gin)

      case('kappa-fission')
        xs = this % xs(t) % kappa_fission(gin)

      case('inverse-velocity')
        xs = this % xs(t) % inverse_velocity(gin)

      case('decay rate')
        if (present(dg)) then
          xs = this % xs(t) % decay_rate(dg)
        else
          xs = this % xs(t) % decay_rate(1)
        end if

      case('prompt-nu-fission')
        xs = this % xs(t) % prompt_nu_fission(gin)

      case('delayed-nu-fission')
        if (present(dg)) then
          xs = this % xs(t) % delayed_nu_fission(dg, gin)
        else
          xs = sum(this % xs(t) % delayed_nu_fission(:, gin))
        end if

      case('nu-fission')
        xs = this % xs(t) % prompt_nu_fission(gin) + &
             sum(this % xs(t) % delayed_nu_fission(:, gin))

      case('chi-prompt')
        if (present(gout)) then
          xs = this % xs(t) % chi_prompt(gout,gin)
        else
          ! Not sure youd want a 1 or a 0, but here you go!
          xs = sum(this % xs(t) % chi_prompt(:, gin))
        end if

      case('chi-delayed')
        if (present(gout)) then
          if (present(dg)) then
            xs = this % xs(t) % chi_delayed(dg, gout, gin)
          else
            xs = this % xs(t) % chi_delayed(1, gout, gin)
          end if
        else
          if (present(dg)) then
            xs = sum(this % xs(t) % chi_delayed(dg, :, gin))
          else
            xs = sum(this % xs(t) % chi_delayed(dg, :, gin))
          end if
        end if

      case('scatter')
        if (present(gout)) then
          if (gout < this % xs(t) % scatter % gmin(gin) .or. &
               gout > this % xs(t) % scatter % gmax(gin)) then
            xs = ZERO
          else
            xs = this % xs(t) % scatter % scattxs(gin) * &
                 this % xs(t) % scatter % energy(gin) % data(gout)
          end if
        else
          xs = this % xs(t) % scatter % scattxs(gin)
        end if

      case('scatter/mult')
        if (present(gout)) then
          if (gout < this % xs(t) % scatter % gmin(gin) .or. &
               gout > this % xs(t) % scatter % gmax(gin)) then
            xs = ZERO
          else
            xs = this % xs(t) % scatter % scattxs(gin) * &
                 this % xs(t) % scatter % energy(gin) % data(gout) / &
                 this % xs(t) % scatter % mult(gin) % data(gout)
          end if
        else
          xs = this % xs(t) % scatter % scattxs(gin) / &
               (dot_product(this % xs(t) % scatter % mult(gin) % data, &
                this % xs(t) % scatter % energy(gin) % data))
        end if

      case('scatter*f_mu/mult','scatter*f_mu')
        if (present(gout)) then
          if (gout < this % xs(t) % scatter % gmin(gin) .or. &
               gout > this % xs(t) % scatter % gmax(gin)) then
            xs = ZERO
          else
            xs = this % xs(t) % scatter % scattxs(gin) * &
                 this % xs(t) % scatter % energy(gin) % data(gout) * &
                 this % xs(t) % scatter % calc_f(gin, gout, mu)
            if (xstype == 'scatter*f_mu/mult') then
              xs = xs / this % xs(t) % scatter % mult(gin) % data(gout)
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

    pure function mgxsang_get_xs(this, xstype, gin, gout, uvw, mu, dg) result(xs)
      class(MgxsAngle), intent(in)  :: this   ! The Mgxs to initialize
      character(*) , intent(in)     :: xstype ! Type of xs requested
      integer, intent(in)           :: gin    ! Incoming Energy group
      integer, optional, intent(in) :: gout   ! Outgoing Energy group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      integer, optional, intent(in) :: dg     ! Delayed group
      real(8)                       :: xs ! Requested x/s

      integer :: iazi, ipol, t

      t = this % index_temp

      if (present(uvw)) then

        call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)

        select case(xstype)

        case('total')
          xs = this % xs(t) % total(gin, iazi, ipol)

        case('absorption')
          xs = this % xs(t) % absorption(gin, iazi, ipol)

        case('fission')
          xs = this % xs(t) % fission(gin, iazi, ipol)

        case('kappa-fission')
          xs = this % xs(t) % kappa_fission(gin, iazi, ipol)

        case('prompt-nu-fission')
          xs = this % xs(t) % prompt_nu_fission(gin, iazi, ipol)

        case('delayed-nu-fission')
          if (present(dg)) then
            xs = this % xs(t) % delayed_nu_fission(dg, gin, iazi, ipol)
          else
            xs = sum(this % xs(t) % delayed_nu_fission(:, gin, iazi, ipol))
          end if

        case('nu-fission')
          xs = this % xs(t) % prompt_nu_fission(gin, iazi, ipol) + &
               sum(this % xs(t) % delayed_nu_fission(:, gin, iazi, ipol))

        case('chi-prompt')
          if (present(gout)) then
            xs = this % xs(t) % chi_prompt(gout, gin, iazi, ipol)
          else
            ! Not sure you would want a 1 or a 0, but here you go!
            xs = sum(this % xs(t) % chi_prompt(:, gin, iazi, ipol))
          end if

        case('chi-delayed')
          if (present(gout)) then
            if (present(dg)) then
              xs = this % xs(t) % chi_delayed(dg, gout, gin, iazi, ipol)
            else
              xs = this % xs(t) % chi_delayed(1, gout, gin, iazi, ipol)
            end if
          else
            if (present(dg)) then
              xs = sum(this % xs(t) % chi_delayed(dg, :, gin, iazi, ipol))
            else
              xs = sum(this % xs(t) % chi_delayed(1, :, gin, iazi, ipol))
            end if
          end if

        case('decay rate')
          if (present(dg)) then
            xs = this % xs(t) % decay_rate(iazi, ipol, dg)
          else
            xs = this % xs(t) % decay_rate(iazi, ipol, 1)
          end if

        case('inverse-velocity')
          xs = this % xs(t) % inverse_velocity(gin, iazi, ipol)

        case('scatter')
          if (present(gout)) then
            if (gout < this % xs(t) % scatter(iazi, ipol) % obj % gmin(gin) .or. &
                 gout > this % xs(t) % scatter(iazi, ipol) % obj % gmax(gin)) then
              xs = ZERO
            else
              xs = this % xs(t) % scatter(iazi, ipol) % obj % scattxs(gin) * &
                   this % xs(t) % scatter(iazi, ipol) % obj % energy(gin) % data(gout)
            end if
          else
            xs = this % xs(t) % scatter(iazi, ipol) % obj % scattxs(gin)
          end if

        case('scatter/mult')
          if (present(gout)) then
            if (gout < this % xs(t) % scatter(iazi, ipol) % obj % gmin(gin) .or. &
                 gout > this % xs(t) % scatter(iazi, ipol) % obj % gmax(gin)) then
              xs = ZERO
            else
              xs = this % xs(t) % scatter(iazi, ipol) % obj % scattxs(gin) * &
                   this % xs(t) % scatter(iazi, ipol) % obj % energy(gin) % data(gout) / &
                   this % xs(t) % scatter(iazi, ipol) % obj % mult(gin) % data(gout)
            end if
          else
            xs = this % xs(t) % scatter(iazi, ipol) % obj % scattxs(gin) / &
                 (dot_product(this % xs(t) % scatter(iazi, ipol) % obj % mult(gin) % data, &
                  this % xs(t) % scatter(iazi, ipol) % obj % energy(gin) % data))
          end if

        case('scatter*f_mu/mult','scatter*f_mu')
          if (present(gout)) then
            if (gout < this % xs(t) % scatter(iazi, ipol) % obj % gmin(gin) .or. &
                 gout > this % xs(t) % scatter(iazi, ipol) % obj % gmax(gin)) then
              xs = ZERO
            else
              xs = this % xs(t) % scatter(iazi, ipol) % obj % scattxs(gin) * &
                   this % xs(t) % scatter(iazi, ipol) % obj % energy(gin) % data(gout)
              xs = xs * this % xs(t) % scatter(iazi, ipol) % obj % calc_f(gin, gout, mu)
              if (xstype == 'scatter*f_mu/mult') then
                xs = xs / &
                     this % xs(t) % scatter(iazi, ipol) % obj % mult(gin) % data(gout)
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
! MGXS*_SAMPLE_FISSION_ENERGY samples the outgoing energy from a fission event
!===============================================================================

    subroutine mgxsiso_sample_fission_energy(this, gin, uvw, dg, gout)

      class(MgxsIso), intent(in)    :: this   ! Data to work with
      integer, intent(in)           :: gin    ! Incoming energy group
      real(8), intent(in)           :: uvw(3) ! Particle Direction
      integer, intent(out)          :: dg     ! Delayed group
      integer, intent(out)          :: gout   ! Sampled outgoing group
      real(8) :: xi_pd            ! Our random number for prompt/delayed
      real(8) :: xi_gout          ! Our random number for gout
      real(8) :: prob_gout        ! Running probability for gout

      ! Get nu and nu_prompt
      real(8) :: prob_prompt

      prob_prompt = this % get_xs('prompt-nu-fission', gin) / &
           this % get_xs('nu-fission', gin)

      ! Sample random numbers
      xi_pd = prn()
      xi_gout = prn()

      ! Neutron is born prompt
      if (xi_pd <= prob_prompt) then

        ! set the delayed group for the particle born from fission to 0
        dg = 0

        gout = 1
        prob_gout = this % get_xs('chi-prompt', gin, gout)

        do while (prob_gout < xi_gout)
          gout = gout + 1
          prob_gout = prob_gout + this % get_xs('chi-prompt', gin, gout)
        end do

        ! Neutron is born delayed
      else

        ! Get the delayed group
        dg = 0

        do while (xi_pd >= prob_prompt)
          dg = dg + 1
          prob_prompt = prob_prompt + &
               this % get_xs('delayed-nu-fission', gin, dg=dg) &
               / this % get_xs('nu-fission', gin)
        end do

        ! Adjust dg in case of round off error
        dg = min(dg, this % num_delayed_groups)

        ! Get the outgoing group
        gout = 1
        prob_gout = this % get_xs('chi-delayed', gin, gout, dg=dg)

        do while (prob_gout < xi_gout)
          gout = gout + 1
          prob_gout = prob_gout + this % get_xs('chi-delayed', gin, gout, dg=dg)
        end do
      end if

    end subroutine mgxsiso_sample_fission_energy

    subroutine mgxsang_sample_fission_energy(this, gin, uvw, dg, gout)
      class(MgxsAngle), intent(in) :: this  ! Data to work with
      integer, intent(in)          :: gin    ! Incoming energy group
      real(8), intent(in)          :: uvw(3) ! Direction vector
      integer, intent(out)         :: dg     ! Delayed group
      integer, intent(out)         :: gout   ! Sampled outgoing group
      real(8) :: xi_pd            ! Our random number for prompt/delayed
      real(8) :: xi_gout          ! Our random number for gout
      real(8) :: prob_gout        ! Running probability for gout
      real(8) :: prob_prompt

      ! Get nu and nu_prompt
      prob_prompt = this % get_xs('prompt-nu-fission', gin, uvw=uvw) / &
           this % get_xs('nu-fission', gin, uvw=uvw)

      ! Sample random numbers
      xi_pd = prn()
      xi_gout = prn()

      ! Neutron is born prompt
      if (xi_pd <= prob_prompt) then

        ! set the delayed group for the particle born from fission to 0
        dg = 0

        gout = 1
        prob_gout = this % get_xs('chi-prompt', gin, gout, uvw=uvw)

        do while (prob_gout < xi_gout)
          gout = gout + 1
          prob_gout = prob_gout + &
               this % get_xs('chi-prompt', gin, gout, uvw=uvw)
        end do

        ! Neutron is born delayed
      else

        ! Get the delayed group
        dg = 0

        do while (xi_pd >= prob_prompt)
          dg = dg + 1
          prob_prompt = prob_prompt + &
               this % get_xs('delayed-nu-fission', gin, uvw=uvw, dg=dg) / &
               this % get_xs('nu-fission', gin, uvw=uvw)
        end do

        ! Adjust dg in case of round off error
        dg = min(dg, this % num_delayed_groups)

        ! Get the outgoing group
        gout = 1
        prob_gout = this % get_xs('chi-delayed', gin, gout, uvw=uvw, dg=dg)

        do while (prob_gout < xi_gout)
          gout = gout + 1
          prob_gout = prob_gout + &
               this % get_xs('chi-delayed', gin, gout, uvw=uvw, dg=dg)
        end do
      end if

    end subroutine mgxsang_sample_fission_energy

!===============================================================================
! MGXS*_SAMPLE_SCATTER Selects outgoing energy and angle after a scatter event
!===============================================================================

    subroutine mgxsiso_sample_scatter(this, uvw, gin, gout, mu, wgt)
      class(MgxsIso), intent(in)    :: this
      real(8),        intent(in)    :: uvw(3) ! Incoming neutron direction
      integer,        intent(in)    :: gin    ! Incoming neutron group
      integer,        intent(out)   :: gout   ! Sampled outgoin group
      real(8),        intent(out)   :: mu     ! Sampled change in angle
      real(8),        intent(inout) :: wgt    ! Particle weight

      call this % xs(this % index_temp) % scatter % sample(gin, gout, mu, wgt)

    end subroutine mgxsiso_sample_scatter

    subroutine mgxsang_sample_scatter(this, uvw, gin, gout, mu, wgt)
      class(MgxsAngle), intent(in)    :: this
      real(8),          intent(in)    :: uvw(3) ! Incoming neutron direction
      integer,          intent(in)    :: gin    ! Incoming neutron group
      integer,          intent(out)   :: gout   ! Sampled outgoin group
      real(8),          intent(out)   :: mu     ! Sampled change in angle
      real(8),          intent(inout) :: wgt    ! Particle weight

      integer :: iazi, ipol ! Angular indices

      call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)
      call this % xs(this % index_temp) % scatter(iazi, ipol) % obj % sample( &
           gin, gout, mu, wgt)

    end subroutine mgxsang_sample_scatter

!===============================================================================
! MGXS*_CALCULATE_XS determines the multi-group cross sections
! for the material the particle is currently traveling through.
!===============================================================================

    subroutine mgxsiso_calculate_xs(this, gin, sqrtkT, uvw, xs)
      class(MgxsIso),        intent(inout) :: this
      integer,               intent(in)    :: gin    ! Incoming neutron group
      real(8),               intent(in)    :: sqrtkT ! Material temperature
      real(8),               intent(in)    :: uvw(3) ! Incoming neutron direction
      type(MaterialMacroXS), intent(inout) :: xs     ! Resultant Mgxs Data

      ! Update the temperature index
      call this % find_temperature(sqrtkT)

      xs % total         = this % xs(this % index_temp) % total(gin)
      xs % absorption    = this % xs(this % index_temp) % absorption(gin)
      xs % nu_fission    = &
           this % xs(this % index_temp) % prompt_nu_fission(gin) + &
           sum(this % xs(this % index_temp) % delayed_nu_fission(:, gin))

    end subroutine mgxsiso_calculate_xs

    subroutine mgxsang_calculate_xs(this, gin, sqrtkT, uvw, xs)
      class(MgxsAngle),      intent(inout) :: this
      integer,               intent(in)    :: gin    ! Incoming neutron group
      real(8),               intent(in)    :: sqrtkT ! Material temperature
      real(8),               intent(in)    :: uvw(3) ! Incoming neutron direction
      type(MaterialMacroXS), intent(inout) :: xs     ! Resultant Mgxs Data

      integer :: iazi, ipol

      ! Update the temperature and angle indices
      call this % find_temperature(sqrtkT)
      call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)

      xs % total         = this % xs(this % index_temp) % &
           total(gin, iazi, ipol)
      xs % absorption    = this % xs(this % index_temp) % &
           absorption(gin, iazi, ipol)
      xs % nu_fission    = this % xs(this % index_temp) % &
           prompt_nu_fission(gin, iazi, ipol) + &
           sum(this % xs(this % index_temp) % &
           delayed_nu_fission(:, gin, iazi, ipol))

    end subroutine mgxsang_calculate_xs

!===============================================================================
! MGXS_FIND_TEMPERATURE sets the temperature index for the given
! sqrt(temperature), (with temperature in units of eV)
!===============================================================================

    subroutine mgxs_find_temperature(this, sqrtkT)
      class(Mgxs), intent(inout) :: this
      real(8), intent(in)        :: sqrtkT    ! Temperature (in units of eV)

      this % index_temp = minloc(abs(this % kTs - (sqrtkT * sqrtkT)), dim=1)

    end subroutine mgxs_find_temperature

!===============================================================================
! FIND_ANGLE finds the closest angle on the data grid and returns that index
!===============================================================================

    pure subroutine find_angle(polar, azimuthal, uvw, i_azi, i_pol)
      real(8), intent(in) :: polar(:)     ! Polar angles [0,pi]
      real(8), intent(in) :: azimuthal(:) ! Azi. angles [-pi,pi]
      real(8), intent(in) :: uvw(3)       ! Direction of motion
      integer, intent(inout) :: i_pol     ! Closest polar bin
      integer, intent(inout) :: i_azi     ! Closest azi bin

      real(8) :: my_pol, my_azi, dangle

      ! Convert uvw to polar and azi

      my_pol = acos(uvw(3))
      my_azi = atan2(uvw(2), uvw(1))

      ! Search for equi-binned angles
      dangle = PI / real(size(polar),8)
      i_pol  = floor(my_pol / dangle + ONE)
      dangle = TWO * PI / real(size(azimuthal),8)
      i_azi  = floor((my_azi + PI) / dangle + ONE)

    end subroutine find_angle

!===============================================================================
! FREE_MEMORY_MGXS deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_mgxs()
    if (allocated(nuclides_MG)) deallocate(nuclides_MG)
    if (allocated(macro_xs)) deallocate(macro_xs)
    if (allocated(energy_bins)) deallocate(energy_bins)
    if (allocated(energy_bin_avg)) deallocate(energy_bin_avg)
  end subroutine free_memory_mgxs

end module mgxs_header
