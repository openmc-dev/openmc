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
    real(8), allocatable :: nu_fission(:)    ! fission production
    real(8), allocatable :: k_fission(:)     ! kappa-fission
    real(8), allocatable :: fission(:)       ! fission
    real(8), allocatable :: chi(:, :)        ! Fission Spectra
    real(8), allocatable :: inv_vel(:)       ! Inverse velocities
  end type XsDataIso

  type :: XsDataAngle
    ! Microscopic cross sections
    ! In all cases, right-most indices are theta, phi
    real(8), allocatable :: total(:, :, :)        ! total cross section
    real(8), allocatable :: absorption(:, :, :)   ! absorption cross section
    type(ScattDataContainer), allocatable :: scatter(:, :) ! scattering info
    real(8), allocatable :: nu_fission(:, :, :)   ! fission production
    real(8), allocatable :: k_fission(:, :, :)    ! kappa-fission
    real(8), allocatable :: fission(:, :, :)      ! fission
    real(8), allocatable :: chi(:, :, :, :)       ! Fission Spectra
    real(8), allocatable :: inv_vel(:, :, :)      ! Inverse velocities
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
    subroutine mgxs_from_hdf5_(this, xs_id, groups, temperature, method, &
                               tolerance, get_kfiss, get_fiss, max_order, &
                               legendre_to_tabular, legendre_to_tabular_points)
      import Mgxs, HID_T, VectorReal
      class(Mgxs), intent(inout)   :: this        ! Working Object
      integer(HID_T), intent(in)   :: xs_id       ! Library data
      integer, intent(in)          :: groups      ! Number of Energy groups
      type(VectorReal), intent(in) :: temperature ! list of desired temperatures
      integer, intent(inout)       :: method      ! Type of temperature access
      real(8), intent(in)          :: tolerance   ! Tolerance on method
      logical, intent(in)          :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)          :: get_fiss    ! Should we get fiss data?
      integer, intent(in)          :: max_order   ! Maximum requested order
      logical, intent(in)          :: legendre_to_tabular ! Convert Legendres to Tabular?
      integer, intent(in)          :: legendre_to_tabular_points ! Number of points to use
                                                                 ! in that  conversion
    end subroutine mgxs_from_hdf5_

    subroutine mgxs_combine_(this, temps, mat, nuclides, groups, max_order, &
                             tolerance, method)
      import Mgxs, Material, MgxsContainer, VectorReal
      class(Mgxs), intent(inout)          :: this ! The Mgxs to initialize
      type(VectorReal), intent(in)        :: temps ! Temperatures to obtain
      type(Material), pointer, intent(in) :: mat  ! base material
      type(MgxsContainer), intent(in)     :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                 :: groups     ! Number of E groups
      integer, intent(in)                 :: max_order  ! Maximum requested order
      real(8), intent(in)                 :: tolerance  ! Tolerance on method
      integer, intent(in)                 :: method     ! Type of temperature access
    end subroutine mgxs_combine_

    pure function mgxs_get_xs_(this, xstype, gin, gout, uvw, mu) result(xs_val)
      import Mgxs
      class(Mgxs), intent(in)       :: this
      character(*), intent(in)      :: xstype ! Cross Section Type
      integer, intent(in)           :: gin    ! Incoming Energy group
      integer, optional, intent(in) :: gout   ! Outgoing Group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      real(8)                       :: xs_val ! Resultant xs
    end function mgxs_get_xs_

    function mgxs_sample_fission_(this, gin, uvw) result(gout)
      import Mgxs
      class(Mgxs), intent(in) :: this
      integer, intent(in)     :: gin    ! Incoming energy group
      real(8), intent(in)     :: uvw(3) ! Particle Direction
      integer                 :: gout   ! Sampled outgoing group

    end function mgxs_sample_fission_

    subroutine mgxs_sample_scatter_(this, uvw, gin, gout, mu, wgt)
      import Mgxs
      class(Mgxs), intent(in)       :: this
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
        if (to_lower(temp_str) /= "[order][g][g']") then
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

    subroutine mgxsiso_from_hdf5(this, xs_id, groups, temperature, method, &
                                 tolerance, get_kfiss, get_fiss, max_order, &
                                 legendre_to_tabular, legendre_to_tabular_points)
      class(MgxsIso), intent(inout) :: this        ! Working Object
      integer(HID_T), intent(in)    :: xs_id       ! Group in H5 file
      integer, intent(in)           :: groups      ! Number of Energy groups
      type(VectorReal), intent(in)  :: temperature ! list of desired temperatures
      integer, intent(inout)        :: method      ! Type of temperature access
      real(8), intent(in)           :: tolerance   ! Tolerance on method
      logical, intent(in)           :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)           :: get_fiss    ! Should we get fiss data?
      integer, intent(in)           :: max_order   ! Maximum requested order
      logical, intent(in)           :: legendre_to_tabular ! Convert Legendres to Tabular?
      integer, intent(in)           :: legendre_to_tabular_points ! Number of points to use
                                                                  ! in that  conversion

      character(MAX_LINE_LEN)     :: temp_str
      integer(HID_T)              :: xsdata_grp, scatt_grp
      real(8), allocatable        :: temp_arr(:), temp_2d(:, :)
      real(8)                     :: dmu, mu, norm
      integer                     :: order, order_dim, gin, gout, l, imu, length
      type(VectorInt)             :: temps_to_read
      integer                     :: t
      type(Jagged2D), allocatable :: input_scatt(:), scatt_coeffs(:)
      type(Jagged1D), allocatable :: temp_mult(:)
      integer, allocatable        :: gmin(:), gmax(:)

      ! Call generic data gathering routine (will populate the metadata)
      call mgxs_from_hdf5(this, xs_id, temperature, method, tolerance, &
                          temps_to_read, order_dim)

      ! Load the more specific data
      do t = 1, temps_to_read % size()
        associate(xs => this % xs(t))
          ! Get temperature as a string
          temp_str = trim(to_str(temps_to_read % data(t))) // "K"
          xsdata_grp = open_group(xs_id, trim(temp_str))
          allocate(xs % nu_fission(groups))
          allocate(xs % chi(groups, groups))
          if (this % fissionable) then
            if (object_exists(xsdata_grp, "chi")) then
              ! Chi was provided, that means we need chi and nu-fission vectors
              ! Get chi
              allocate(temp_arr(groups))
              call read_dataset(temp_arr, xsdata_grp, "chi")
              do gin = 1, groups
                do gout = 1, groups
                  xs % chi(gout, gin) = temp_arr(gout)
                end do
                ! Normalize chi so its CDF goes to 1
                xs % chi(:, gin) = xs % chi(:, gin) / sum(xs % chi(:, gin))
              end do
              deallocate(temp_arr)

              ! Get nu_fission (as a vector)
              if (object_exists(xsdata_grp, "nu-fission")) then
                call read_dataset(xs % nu_fission, xsdata_grp, "nu-fission")
              else
                call fatal_error("If fissionable, must provide nu-fission!")
              end if

            else
              ! chi isnt provided but is within nu_fission, existing as a matrix
              ! So, get nu_fission (as a matrix)
              if (object_exists(xsdata_grp, "nu-fission")) then
                allocate(temp_2d(groups, groups))
                call read_dataset(temp_2d, xsdata_grp, "nu-fission")
              else
                call fatal_error("If fissionable, must provide nu-fission!")
              end if

              ! Set the vector nu-fission from the matrix nu-fission
              do gin = 1, groups
                xs % nu_fission(gin) = sum(temp_2d(:, gin))
              end do

              ! Now pull out information needed for chi
              xs % chi(:, :) = temp_2d
              ! Normalize chi so its CDF goes to 1
              do gin = 1, groups
                xs % chi(:, gin) = xs % chi(:, gin) / sum(xs % chi(:, gin))
              end do
              deallocate(temp_2d)
            end if

            ! If we have a need* for the fission and kappa-fission x/s, get them
            ! (*Need is defined as will be using it to tally)
            if (get_fiss) then
              if (object_exists(xsdata_grp, "fission")) then
                allocate(xs % fission(groups))
                call read_dataset(xs % fission, xsdata_grp, "fission")
              else
                call fatal_error("Fission data missing, required due to fission&
                                 & tallies in tallies.xml file!")
              end if
            end if
            if (get_kfiss) then
              if (object_exists(xsdata_grp, "kappa-fission")) then
                allocate(xs % k_fission(groups))
                call read_dataset(xs % k_fission, xsdata_grp, "kappa-fission")
              else
                call fatal_error("kappa-fission data missing, required due to &
                                 &kappa-fission tallies in tallies.xml file!")
              end if
            end if
          else
            xs % nu_fission = ZERO
            xs % chi = ZERO
          end if

          if (object_exists(xsdata_grp, "absorption")) then
            allocate(xs % absorption(groups))
            call read_dataset(xs % absorption, xsdata_grp, "absorption")
          else
            call fatal_error("Must provide absorption!")
          end if

          ! Get scattering data
          if (.not. object_exists(xsdata_grp, "scatter_data")) &
               call fatal_error("Must provide 'scatter_data'")
          scatt_grp = open_group(xsdata_grp, 'scatter_data')
          ! First get the outgoing group boundary indices
          if (object_exists(scatt_grp, "g_min")) then
            allocate(gmin(groups))
            call read_dataset(gmin, scatt_grp, "g_min")
          else
            call fatal_error("'g_min' for the scatter_data must be provided")
          end if
          if (object_exists(scatt_grp, "g_max")) then
            allocate(gmax(groups))
            call read_dataset(gmax, scatt_grp, "g_max")
          else
            call fatal_error("'g_max' for the scatter_data must be provided")
          end if

          ! Now use this information to find the length of a container array
          ! to hold the flattened data
          length = 0
          do gin = 1, groups
            length = length + order_dim * (gmax(gin) - gmin(gin) + 1)
          end do
          ! Allocate flattened array
          allocate(temp_arr(length))
          if (.not. object_exists(scatt_grp, 'scatter_matrix')) &
               call fatal_error("'scatter_matrix' must be provided")
          call read_dataset(temp_arr, scatt_grp, "scatter_matrix")

          ! Compare the number of orders given with the maximum order of the
          ! problem.  Strip off the supefluous orders if needed.
          if (this % scatter_format == ANGLE_LEGENDRE) then
            order = min(order_dim - 1, max_order)
            order_dim = order + 1
          end if

          ! Convert temp_arr to a jagged array ((gin) % data(l, gout)) for passing
          ! to ScattData
          allocate(input_scatt(groups))
          index = 1
          do gin = 1, groups
            allocate(input_scatt(gin) % data(order_dim, gmin(gin):gmax(gin)))
            do gout = gmin(gin), gmax(gin)
              do l = 1, order_dim
                input_scatt(gin) % data(l, gout) = temp_arr(index)
                index = index + 1
              end do
            end do
          end do
          deallocate(temp_arr)

          ! Finally convert the legendre to tabular if needed
          allocate(scatt_coeffs(groups))
          if (this % scatter_format == ANGLE_LEGENDRE .and. &
               legendre_to_tabular) then
            this % scatter_format = ANGLE_TABULAR
            order_dim = legendre_to_tabular_points
            order = order_dim
            dmu = TWO / real(order - 1, 8)
            do gin = 1, groups
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
                ! Now that we have the integral, lets ensure that the distribution
                ! is normalized such that it preserves the original scattering xs
                if (norm > ZERO) then
                  scatt_coeffs(gin) % data(:, gout) = &
                       scatt_coeffs(gin) % data(:, gout) * &
                       input_scatt(gin) % data(1, gout) / norm
                end if
              end do ! gout
            end do ! gin
          else
            ! Sticking with current representation
            do gin = 1, groups
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
            do gin = 1, groups
              length = length + (gmax(gin) - gmin(gin) + 1)
            end do
            ! Allocate flattened array
            allocate(temp_arr(length))
            call read_dataset(temp_arr, scatt_grp, "multiplicity_matrix")

            ! Convert temp_arr to a jagged array ((gin) % data(gout)) for passing
            ! to ScattData
            allocate(temp_mult(groups))
            index = 1
            do gin = 1, groups
              allocate(temp_mult(gin) % data(gmin(gin):gmax(gin)))
              do gout = gmin(gin), gmax(gin)
                temp_mult(gin) % data(gout) = temp_arr(index)
                index = index + 1
              end do
            end do
            deallocate(temp_arr)
          else
            ! Default to multiplicities of 1.0
            allocate(temp_mult(groups))
            do gin = 1, groups
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
          do gin = 1, groups
            if (xs % absorption(gin) == ZERO) xs % absorption(gin) = 1E-10_8
          end do

          ! Get, or infer, total xs data.
          allocate(xs % total(groups))
          if (object_exists(xsdata_grp, "total")) then
            call read_dataset(xs % total, xsdata_grp, "total")
          else
            xs % total(:) = xs % absorption(:) + xs % scatter % scattxs(:)
          end if

          ! Check sigT to ensure it is not 0 since it is
          ! often divided by in the tally routines
          do gin = 1, groups
            if (xs % total(gin) == ZERO) xs % total(gin) = 1E-10_8
          end do

          ! Get kinetics data
          if (object_exists(xsdata_grp, "inverse_velocities")) then
            allocate(xs % inv_vel(groups))
            call read_dataset(xs % inv_vel, xsdata_grp, "inverse_velocities")
          end if

          ! Close the groups we have opened and deallocate
          call close_group(xsdata_grp)
          call close_group(scatt_grp)
          deallocate(scatt_coeffs, temp_mult)
        end associate ! xs
      end do ! Temperatures
    end subroutine mgxsiso_from_hdf5

    subroutine mgxsang_from_hdf5(this, xs_id, groups, temperature, method, &
                                 tolerance, get_kfiss, get_fiss, max_order, &
                                 legendre_to_tabular, legendre_to_tabular_points)
      class(MgxsAngle), intent(inout) :: this        ! Working Object
      integer(HID_T), intent(in)      :: xs_id       ! Group in H5 file
      integer, intent(in)             :: groups      ! Number of Energy groups
      type(VectorReal), intent(in)    :: temperature ! list of desired temperatures
      integer, intent(inout)          :: method      ! Type of temperature access
      real(8), intent(in)             :: tolerance   ! Tolerance on method
      logical, intent(in)             :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)             :: get_fiss    ! Should we get fiss data?
      integer, intent(in)             :: max_order   ! Maximum requested order
      logical, intent(in)             :: legendre_to_tabular ! Convert Legendres to Tabular?
      integer, intent(in)             :: legendre_to_tabular_points ! Number of points to use
                                                                    ! in that  conversion

      character(MAX_LINE_LEN)     :: temp_str
      integer(HID_T)              :: xsdata_grp, scatt_grp
      integer, allocatable        :: int_arr(:)
      real(8), allocatable        :: temp_arr(:), temp_4d(:, :, :, :)
      real(8)                     :: dmu, mu, norm
      integer                     :: order, order_dim, gin, gout, l, imu
      type(VectorInt)             :: temps_to_read
      integer                     :: t, length, ipol, iazi
      type(Jagged2D), allocatable :: input_scatt(:, :, :), scatt_coeffs(:, :, :)
      type(Jagged1D), allocatable :: temp_mult(:, :, :)
      integer, allocatable        :: gmin(:, :, :), gmax(:, :, :)

      ! Call generic data gathering routine (will populate the metadata)
      call mgxs_from_hdf5(this, xs_id, temperature, method, tolerance, &
                          temps_to_read, order_dim)

      ! Load the more specific data
      do t = 1, temps_to_read % size()
        associate(xs => this % xs(t))
          ! Get temperature as a string
          temp_str = trim(to_str(temps_to_read % data(t))) // "K"
          xsdata_grp = open_group(xs_id, trim(temp_str))
          allocate(xs % nu_fission(groups, this % n_azi, this % n_pol))
          allocate(xs % chi(groups, groups, this % n_azi, this % n_pol))
          if (this % fissionable) then
            if (object_exists(xsdata_grp, "chi")) then
              ! Chi was provided, that means we need chi and nu-fission vectors
              ! Get chi
              allocate(temp_arr(groups * this % n_azi * this % n_pol))
              call read_dataset(temp_arr, xsdata_grp, "chi")
              ! Initialize counter for temp_arr
              l = 0
              gin = 1
              do ipol = 1, this % n_pol
                do iazi = 1, this % n_azi
                  do gout = 1, groups
                    l = l + 1
                    xs % chi(gout, gin, iazi, ipol) = temp_arr(l)
                  end do
                  ! Normalize chi so its CDF goes to 1
                  xs % chi(:, gin, iazi, ipol) = &
                       xs % chi(:, gin, iazi, ipol) / &
                       sum(xs % chi(:, gin, iazi, ipol))
                end do
              end do

              ! Now set all the other gin values
              do ipol = 1, this % n_pol
                do iazi = 1, this % n_azi
                  do gin = 2, groups
                    xs % chi(:, gin, iazi, ipol) = &
                         xs % chi(:, 1, iazi, ipol)
                  end do
                end do
              end do
              deallocate(temp_arr)

              ! Get nu_fission (as a vector)
              if (object_exists(xsdata_grp, "nu-fission")) then
                allocate(temp_arr(groups * this % n_azi * this % n_pol))
                call read_dataset(temp_arr, xsdata_grp, "nu-fission")
                xs % nu_fission = reshape(temp_arr, (/groups, this % n_azi, &
                                                      this % n_pol/))
                deallocate(temp_arr)
              else
                call fatal_error("If fissionable, must provide nu-fission!")
              end if

            else
              ! chi isnt provided but is within nu_fission, existing as a matrix
              ! So, get nu_fission (as a matrix)
              if (object_exists(xsdata_grp, "nu-fission")) then
                allocate(temp_arr(groups * groups * this % n_azi * this % n_pol))
                call read_dataset(temp_arr, xsdata_grp, "nu-fission")
                allocate(temp_4d(groups, groups, this % n_azi, this % n_pol))
                temp_4d(:, :, :, :) = reshape(temp_arr, (/groups, groups, &
                     this % n_azi, this % n_pol/))
                deallocate(temp_arr)
              else
                call fatal_error("If fissionable, must provide nu-fission!")
              end if

              ! Set the vector nu-fission from the matrix nu-fission
              do ipol = 1, this % n_pol
                do iazi = 1, this % n_azi
                  do gin = 1, groups
                    xs % nu_fission(gin, iazi, ipol) = &
                         sum(temp_4d(:, gin, iazi, ipol))
                  end do
                end do
              end do

              ! Now pull out information needed for chi
              xs % chi = temp_4d
              ! Normalize chi so its CDF goes to 1
              do ipol = 1, this % n_pol
                do iazi = 1, this % n_azi
                  do gin = 1, groups
                    xs % chi(:, gin, iazi, ipol) = &
                         xs % chi(:, gin, iazi, ipol) / &
                         sum(xs % chi(:, gin, iazi, ipol))
                  end do
                end do
              end do
              deallocate(temp_4d)
            end if

            ! If we have a need* for the fission and kappa-fission x/s, get them
            ! (*Need is defined as will be using it to tally)
            if (get_fiss) then
              if (object_exists(xsdata_grp, "fission")) then
                allocate(temp_arr(groups * this % n_azi * this % n_pol))
                call read_dataset(temp_arr, xsdata_grp, "fission")
                allocate(xs % fission(groups, this % n_azi, this % n_pol))
                xs % fission = reshape(temp_arr, (/groups, this % n_azi, &
                                                   this % n_pol/))
                deallocate(temp_arr)
              else
                call fatal_error("Fission data missing, required due to fission&
                                 & tallies in tallies.xml file!")
              end if
            end if
            if (get_kfiss) then
              if (object_exists(xsdata_grp, "kappa-fission")) then
                allocate(temp_arr(groups * this % n_azi * this % n_pol))
                call read_dataset(temp_arr, xsdata_grp, "kappa-fission")
                allocate(xs % k_fission(groups, this % n_azi, this % n_pol))
                xs % k_fission = reshape(temp_arr, (/groups, this % n_azi, &
                                                   this % n_pol/))
                deallocate(temp_arr)
              else
                call fatal_error("kappa_fission data missing, required due to &
                                 &kappa-fission tallies in tallies.xml file!")
              end if
            end if
          else
            xs % nu_fission = ZERO
            xs % chi = ZERO
          end if

          if (object_exists(xsdata_grp, "absorption")) then
            allocate(temp_arr(groups * this % n_azi * this % n_pol))
            call read_dataset(temp_arr, xsdata_grp, "absorption")
            allocate(xs % absorption(groups, this % n_azi, this % n_pol))
            xs % absorption = reshape(temp_arr, (/groups, this % n_azi, &
                                                  this % n_pol/))
            deallocate(temp_arr)
          else
            call fatal_error("Must provide absorption!")
          end if

          ! Get scattering data
          if (.not. object_exists(xsdata_grp, "scatter_data")) &
               call fatal_error("Must provide 'scatter_data'")
          scatt_grp = open_group(xsdata_grp, 'scatter_data')
          ! First get the outgoing group boundary indices
          if (object_exists(scatt_grp, "g_min")) then
            allocate(int_arr(groups * this % n_azi * this % n_pol))
            call read_dataset(int_arr, scatt_grp, "g_min")
            allocate(gmin(groups, this % n_azi, this % n_pol))
            gmin = reshape(int_arr, (/groups, this % n_azi, this % n_pol/))
            deallocate(int_arr)
          else
            call fatal_error("'g_min' for the scatter_data must be provided")
          end if
          if (object_exists(scatt_grp, "g_max")) then
            allocate(int_arr(groups * this % n_azi * this % n_pol))
            call read_dataset(int_arr, scatt_grp, "g_max")
            allocate(gmax(groups, this % n_azi, this % n_pol))
            gmax = reshape(int_arr, (/groups, this % n_azi, this % n_pol/))
            deallocate(int_arr)
          else
            call fatal_error("'g_max' for the scatter_data must be provided")
          end if

          ! Now use this information to find the length of a container array
          ! to hold the flattened data
          length = 0
          do ipol = 1, this % n_pol
            do iazi = 1, this % n_azi
              do gin = 1, groups
                length = length + order_dim * (gmax(gin, iazi, ipol) - &
                                               gmin(gin, iazi, ipol) + 1)
              end do
            end do
          end do
          ! Allocate flattened array
          allocate(temp_arr(length))
          if (.not. object_exists(scatt_grp, 'scatter_matrix')) &
               call fatal_error("'scatter_matrix' must be provided")
          call read_dataset(temp_arr, scatt_grp, "scatter_matrix")

          ! Compare the number of orders given with the maximum order of the
          ! problem.  Strip off the superfluous orders if needed.
          if (this % scatter_format == ANGLE_LEGENDRE) then
            order = min(order_dim - 1, max_order)
            order_dim = order + 1
          end if

          ! Convert temp_arr to a jagged array ((gin) % data(l, gout)) for passing
          ! to ScattData
          allocate(input_scatt(groups, this % n_azi, this % n_pol))
          index = 1
          do ipol = 1, this % n_pol
            do iazi = 1, this % n_azi
              do gin = 1, groups
                allocate(input_scatt(gin, iazi, ipol) % data(order_dim, &
                     gmin(gin, iazi, ipol):gmax(gin, iazi, ipol)))
                do gout = gmin(gin, iazi, ipol), gmax(gin, iazi, ipol)
                  do l = 1, order_dim
                    input_scatt(gin, iazi, ipol) % data(l, gout) = &
                         temp_arr(index)
                    index = index + 1
                  end do ! gout
                end do ! order
              end do ! gin
            end do ! iazi
          end do ! ipol
          deallocate(temp_arr)

          ! Finally convert the legendre to tabular if needed
          allocate(scatt_coeffs(groups, this % n_azi, this % n_pol))
          if (this % scatter_format == ANGLE_LEGENDRE .and. &
               legendre_to_tabular) then
            this % scatter_format = ANGLE_TABULAR
            order_dim = legendre_to_tabular_points
            order = order_dim
            dmu = TWO / real(order - 1, 8)
            do ipol = 1, this % n_pol
              do iazi = 1, this % n_azi
                do gin = 1, groups
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
                do gin = 1, groups
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
                do gin = 1, groups
                  length = length + (gmax(gin, iazi, ipol) - gmin(gin, iazi, ipol) + 1)
                end do
              end do
            end do
            ! Allocate flattened array
            allocate(temp_arr(length))
            call read_dataset(temp_arr, scatt_grp, "multiplicity_matrix")
            ! Convert temp_arr to a jagged array ((gin) % data(gout)) for passing
            ! to ScattData
            allocate(temp_mult(groups, this % n_azi, this % n_pol))
            index = 1
            do ipol = 1, this % n_pol
              do iazi = 1, this % n_azi
                do gin = 1, groups
                  allocate(temp_mult(gin, iazi, ipol) % data( &
                       gmin(gin, iazi, ipol):gmax(gin, iazi, ipol)))
                  do gout = gmin(gin, iazi, ipol), gmax(gin, iazi, ipol)
                    temp_mult(gin, iazi, ipol) % data(gout) = temp_arr(index)
                    index = index + 1
                  end do
                end do
              end do
            end do
            deallocate(temp_arr)
          else
            allocate(temp_mult(groups, this % n_azi, this % n_pol))
            ! Default to multiplicities of 1.0
            do ipol = 1, this % n_pol
              do iazi = 1, this % n_azi
                do gin = 1, groups
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
              do gin = 1, groups
                if (xs % absorption(gin, iazi, ipol) == ZERO) then
                  xs % absorption(gin, iazi, ipol) = 1E-10_8
                end if
              end do
            end do
          end do

          allocate(xs % total(groups, this % n_azi, this % n_pol))
          if (object_exists(xsdata_grp, "total")) then
            allocate(temp_arr(groups * this % n_azi * this % n_pol))
            call read_dataset(temp_arr, xsdata_grp, "total")
            xs % total = reshape(temp_arr, (/groups, this % n_azi, &
                                             this % n_pol/))
            deallocate(temp_arr)
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
              do gin = 1, groups
                if (xs % total(gin, iazi, ipol) == ZERO) then
                  xs % total(gin, iazi, ipol) = 1E-10_8
                end if
              end do
            end do
          end do

          ! Get kinetics data
          if (object_exists(xsdata_grp, "inverse_velocities")) then
            allocate(xs % inv_vel(groups, this % n_azi, this % n_pol))
            allocate(temp_arr(groups * this % n_azi * this % n_pol))
            call read_dataset(temp_arr, xsdata_grp, "inverse_velocities")
            xs % inv_vel = reshape(temp_arr, (/groups, this % n_azi, &
                                               this % n_pol/))
            deallocate(temp_arr)
          end if

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
        order = min(mat_max_order, max_order)
        ! Ok, got our order, store the dimensionality
        order_dim = order
      end if

    end subroutine mgxs_combine

    subroutine mgxsiso_combine(this, temps, mat, nuclides, groups, max_order, &
                               tolerance, method)
      class(MgxsIso), intent(inout)       :: this ! The Mgxs to initialize
      type(VectorReal), intent(in)        :: temps ! Temperatures to obtain [eV]
      type(Material), pointer, intent(in) :: mat  ! base material
      type(MgxsContainer), intent(in)     :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                 :: groups     ! Number of E groups
      integer, intent(in)                 :: max_order  ! Maximum requested order
      real(8), intent(in)                 :: tolerance  ! Tolerance on method
      integer, intent(in)                 :: method     ! Type of temperature access

      integer :: i            ! loop index over nuclides
      integer :: t            ! Index in to temps
      integer :: gin, gout    ! group indices
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

      ! Create the Xs Data for each temperature
      TEMP_LOOP: do t = 1, temps % size()
        ! Allocate and initialize the data needed for macro_xs(i_mat) object
        allocate(this % xs(t) % total(groups))
        this % xs(t) % total(:) = ZERO
        allocate(this % xs(t) % absorption(groups))
        this % xs(t) % absorption(:) = ZERO
        allocate(this % xs(t) % fission(groups))
        this % xs(t) % fission(:) = ZERO
        allocate(this % xs(t) % k_fission(groups))
        this % xs(t) % k_fission(:) = ZERO
        allocate(this % xs(t) % nu_fission(groups))
        this % xs(t) % nu_fission(:) = ZERO
        allocate(this % xs(t) % chi(groups,groups))
        this % xs(t) % chi(:, :) = ZERO
        allocate(this % xs(t) % inv_vel(groups))
        this % xs(t) % inv_vel(:) = ZERO
        allocate(temp_mult(groups,groups))
        temp_mult(:, :) = ZERO
        allocate(mult_num(groups,groups))
        mult_num(:, :) = ZERO
        allocate(mult_denom(groups,groups))
        mult_denom(:, :) = ZERO
        allocate(scatt_coeffs(order_dim,groups,groups))
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
                if (nuc % fissionable) then
                  this % xs(t) % chi = this % xs(t) % chi + &
                       atom_density * nuc % xs(nuc_t) % chi * interp
                  this % xs(t) % nu_fission = this % xs(t) % nu_fission + &
                       atom_density * nuc % xs(nuc_t) % nu_fission * interp
                  if (allocated(nuc % xs(nuc_t) % fission)) then
                    this % xs(t) % fission = this % xs(t) % fission + &
                         atom_density * nuc % xs(nuc_t) % fission * interp
                  end if
                  if (allocated(nuc % xs(nuc_t) % k_fission)) then
                    this % xs(t) % k_fission = this % xs(t) % k_fission + &
                         atom_density * nuc % xs(nuc_t) % k_fission * interp
                  end if
                end if
                if (allocated(nuc % xs(nuc_t) % inv_vel)) then
                  this % xs(t) % inv_vel = this % xs(t) % inv_vel + &
                       atom_density * nuc % xs(nuc_t) % inv_vel * interp
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
                do gin = 1, groups
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
                do gin = 1, groups
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
            do gin = 1, groups
              do gout = 1, groups
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
              do gin = 1, groups
                norm =  sum(this % xs(t) % chi(:, gin))
                if (norm > ZERO) then
                  this % xs(t) % chi(:, gin) = this % xs(t) % chi(:, gin) / norm
                end if
              end do
            end if

            ! Deallocate temporaries
            deallocate(jagged_mult, jagged_scatt, gmin, gmax, scatt_coeffs, &
                       temp_mult, mult_num, mult_denom)
          end associate ! nuc
        end do NUC_LOOP
      end do TEMP_LOOP

    end subroutine mgxsiso_combine

    subroutine mgxsang_combine(this, temps, mat, nuclides, groups, max_order, &
                               tolerance, method)
      class(MgxsAngle), intent(inout)     :: this ! The Mgxs to initialize
      type(VectorReal), intent(in)        :: temps ! Temperatures to obtain
      type(Material), pointer, intent(in) :: mat  ! base material
      type(MgxsContainer), intent(in)     :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                 :: groups     ! Number of E groups
      integer, intent(in)                 :: max_order  ! Maximum requested order
      real(8), intent(in)                 :: tolerance  ! Tolerance on method
      integer, intent(in)                 :: method     ! Type of temperature access

      integer :: i             ! loop index over nuclides
      integer :: t             ! temperature loop index
      integer :: gin, gout     ! group indices
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
        allocate(this % xs(t) % total(groups, n_azi, n_pol))
        this % xs(t) % total = ZERO
        allocate(this % xs(t) % absorption(groups, n_azi, n_pol))
        this % xs(t) % absorption = ZERO
        allocate(this % xs(t) % fission(groups, n_azi, n_pol))
        this % xs(t) % fission = ZERO
        allocate(this % xs(t) % k_fission(groups, n_azi, n_pol))
        this % xs(t) % k_fission = ZERO
        allocate(this % xs(t) % nu_fission(groups, n_azi, n_pol))
        this % xs(t) % nu_fission = ZERO
        allocate(this % xs(t) % chi(groups, groups, n_azi, n_pol))
        this % xs(t) % chi = ZERO
        allocate(this % xs(t) % inv_vel(groups, n_azi, n_pol))
        this % xs(t) % inv_vel = ZERO
        allocate(temp_mult(groups, groups, n_azi, n_pol))
        temp_mult = ZERO
        allocate(mult_num(groups, groups, n_azi, n_pol))
        mult_num = ZERO
        allocate(mult_denom(groups, groups, n_azi, n_pol))
        mult_denom = ZERO
        allocate(scatt_coeffs(order_dim, groups, groups, n_azi, n_pol))
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
                if (allocated(nuc % xs(nuc_t) % inv_vel)) then
                  this % xs(t) % inv_vel = this % xs(t) % inv_vel + &
                       atom_density * nuc % xs(nuc_t) % inv_vel * interp
                end if
                if (nuc % fissionable) then
                  this % xs(t) % chi = this % xs(t) % chi + &
                       atom_density * nuc % xs(nuc_t) % chi * interp
                  this % xs(t) % nu_fission = this % xs(t) % nu_fission + &
                       atom_density * nuc % xs(nuc_t) % nu_fission * interp
                  if (allocated(nuc % xs(nuc_t) % fission)) then
                    this % xs(t) % fission = this % xs(t) % fission + &
                         atom_density * nuc % xs(nuc_t) % fission * interp
                  end if
                  if (allocated(nuc % xs(nuc_t) % k_fission)) then
                    this % xs(t) % k_fission = this % xs(t) % k_fission + &
                         atom_density * nuc % xs(nuc_t) % k_fission * interp
                  end if
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
                    do gin = 1, groups
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
                    do gin = 1, groups
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
                do gin = 1, groups
                  do gout = 1, groups
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
                                          jagged_scatt)
                call jagged_from_dense_1D(temp_mult(:, :, iazi, ipol), &
                                          jagged_mult, gmin, gmax)

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
                  do gin = 1, groups
                    norm =  sum(this % xs(t) % chi(:, gin, iazi, ipol))
                    if (norm > ZERO) then
                      this % xs(t) % chi(:, gin, iazi, ipol) = &
                           this % xs(t) % chi(:, gin, iazi, ipol) / norm
                    end if
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

    pure function mgxsiso_get_xs(this, xstype, gin, gout, uvw, mu) result(xs)
      class(MgxsIso), intent(in)    :: this   ! The Xs to get data from
      character(*) , intent(in)     :: xstype ! Type of xs requested
      integer, intent(in)           :: gin    ! Incoming Energy group
      integer, optional, intent(in) :: gout   ! Outgoing Energy group
      real(8), optional, intent(in) :: uvw(3) ! Requested Angle
      real(8), optional, intent(in) :: mu     ! Change in angle
      real(8)                       :: xs ! Requested x/s
      integer                       :: t ! temperature index

      t = this % index_temp

      select case(xstype)
      case('total')
        xs = this % xs(t) % total(gin)
      case('absorption')
        xs = this % xs(t) % absorption(gin)
      case('fission')
        if (allocated(this % xs(t) % fission)) then
          xs = this % xs(t) % fission(gin)
        else
          xs = ZERO
        end if
      case('kappa_fission')
        if (allocated(this % xs(t) % k_fission)) then
          xs = this % xs(t) % k_fission(gin)
        else
          xs = ZERO
        end if
      case('nu_fission')
        xs = this % xs(t) % nu_fission(gin)
      case('chi')
        if (present(gout)) then
          xs = this % xs(t) % chi(gout,gin)
        else
          ! Not sure youd want a 1 or a 0, but here you go!
          xs = sum(this % xs(t) % chi(:, gin))
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
      case('inv_vel')
        xs = this % xs(t) % inv_vel(gin)
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
          if (allocated(this % xs(t) % fission)) then
            xs = this % xs(t) % fission(gin, iazi, ipol)
          else
            xs = ZERO
          end if
        case('kappa_fission')
          if (allocated(this % xs(t) % k_fission)) then
            xs = this % xs(t) % k_fission(gin, iazi, ipol)
          else
            xs = ZERO
          end if
        case('nu_fission')
          xs = this % xs(t) % nu_fission(gin, iazi, ipol)
        case('chi')
          if (present(gout)) then
            xs = this % xs(t) % chi(gout, gin, iazi, ipol)
          else
            ! Not sure you would want a 1 or a 0, but here you go!
            xs = sum(this % xs(t) % chi(:, gin, iazi, ipol))
          end if
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
        case('inv_vel')
          xs = this % xs(t) % inv_vel(gin, iazi, ipol)
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

    function mgxsiso_sample_fission_energy(this, gin, uvw) result(gout)
      class(MgxsIso), intent(in) :: this   ! Data to work with
      integer, intent(in)        :: gin    ! Incoming energy group
      real(8), intent(in)        :: uvw(3) ! Particle Direction
      integer                    :: gout   ! Sampled outgoing group
      real(8) :: xi   ! Our random number
      real(8) :: prob ! Running probability

      xi = prn()
      gout = 1
      prob = this % xs(this % index_temp) % chi(gout,gin)

      do while (prob < xi)
        gout = gout + 1
        prob = prob + this % xs(this % index_temp) % chi(gout,gin)
      end do

    end function mgxsiso_sample_fission_energy

    function mgxsang_sample_fission_energy(this, gin, uvw) result(gout)
      class(MgxsAngle), intent(in) :: this  ! Data to work with
      integer, intent(in)          :: gin    ! Incoming energy group
      real(8), intent(in)          :: uvw(3) ! Particle Direction
      integer                      :: gout   ! Sampled outgoing group
      real(8) :: xi         ! Our random number
      real(8) :: prob       ! Running probability
      integer :: iazi, ipol ! Angle indices

      call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)

      xi = prn()
      gout = 1
      prob = this % xs(this % index_temp) % chi(gout, gin, iazi, ipol)

      do while (prob < xi)
        gout = gout + 1
        prob = prob + this % xs(this % index_temp) % chi(gout, gin, iazi, ipol)
      end do

    end function mgxsang_sample_fission_energy

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

    subroutine mgxsiso_calculate_xs(this, gin, uvw, xs)
      class(MgxsIso),        intent(in)    :: this
      integer,               intent(in)    :: gin    ! Incoming neutron group
      real(8),               intent(in)    :: uvw(3) ! Incoming neutron direction
      type(MaterialMacroXS), intent(inout) :: xs     ! Resultant Mgxs Data

      xs % total         = this % xs(this % index_temp) % total(gin)
      xs % elastic       = this % xs(this % index_temp) % scatter % scattxs(gin)
      xs % absorption    = this % xs(this % index_temp) % absorption(gin)
      xs % fission       = this % xs(this % index_temp) % fission(gin)
      xs % nu_fission    = this % xs(this % index_temp) % nu_fission(gin)

    end subroutine mgxsiso_calculate_xs

    subroutine mgxsang_calculate_xs(this, gin, uvw, xs)
      class(MgxsAngle),      intent(in)    :: this
      integer,               intent(in)    :: gin    ! Incoming neutron group
      real(8),               intent(in)    :: uvw(3) ! Incoming neutron direction
      type(MaterialMacroXS), intent(inout) :: xs     ! Resultant Mgxs Data

      integer :: iazi, ipol

      call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)
      xs % total         = this % xs(this % index_temp) % &
           total(gin, iazi, ipol)
      xs % elastic       = this % xs(this % index_temp) % &
           scatter(iazi, ipol) % obj % scattxs(gin)
      xs % absorption    = this % xs(this % index_temp) % &
           absorption(gin, iazi, ipol)
      xs % fission       = this % xs(this % index_temp) % &
           fission(gin, iazi, ipol)
      xs % nu_fission    = this % xs(this % index_temp) % &
           nu_fission(gin, iazi, ipol)

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

end module mgxs_header
