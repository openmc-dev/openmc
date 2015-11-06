module scattdata_header

  use math
  use constants

  implicit none

!===============================================================================
! SCATTDATA contains all the data to describe the scattering energy and
! angular distribution
!===============================================================================

  type, abstract :: ScattData_Base
    ! p0 matrix on its own for sampling energy
    real(8), allocatable :: energy(:,:) ! (Gout x Gin)
    real(8), allocatable :: mult(:,:)   ! (Gout x Gin)
    real(8), allocatable :: data(:,:,:) ! (Order/Nmu x Gout x Gin)

    ! Type-Bound procedures
    contains
      procedure(init_), deferred, pass   :: init   ! Initializes ScattData
      procedure(calc_f_), deferred, pass :: calc_f ! Calculates f, given mu
      procedure(clear_), deferred, pass  :: clear  ! Deallocates ScattData
  end type ScattData_Base

  abstract interface
    subroutine init_(this, order, energy, mult, coeffs)
      import ScattData_Base
      class(ScattData_Base), intent(inout) :: this          ! Object to work on
      integer, intent(in)                  :: order         ! Data Order
      real(8), intent(in)                  :: energy(:,:) ! Energy Transfer Matrix
      real(8), intent(in)                  :: mult(:,:)     ! Scatter Prod'n Matrix
      real(8), intent(in)                  :: coeffs(:,:,:) ! Coefficients to use
    end subroutine init_

    pure function calc_f_(this, gin, gout, mu) result(f)
      import ScattData_Base
      class(ScattData_Base), intent(in) :: this ! The ScattData to evaluate
      integer, intent(in)               :: gin  ! Incoming Energy Group
      integer, intent(in)               :: gout ! Outgoing Energy Group
      real(8), intent(in)               :: mu   ! Angle of interest
      real(8)                           :: f    ! Return value of f(mu)

    end function calc_f_

    subroutine clear_(this)
      import ScattData_Base
      class(ScattData_Base), intent(inout) :: this ! The ScattData to clear
    end subroutine clear_
  end interface

  type, extends(ScattData_Base) :: ScattData_Legendre
    contains
      procedure, pass :: init   => scattdata_legendre_init
      procedure, pass :: calc_f => scattdata_legendre_calc_f
      procedure, pass :: clear  => scattdata_legendre_clear
  end type ScattData_Legendre

  type, extends(ScattData_Base) :: ScattData_Histogram
    real(8), allocatable :: mu(:) ! Mu bins
    real(8)              :: dmu   ! Mu spacing
    contains
      procedure, pass :: init   => scattdata_histogram_init
      procedure, pass :: calc_f => scattdata_histogram_calc_f
      procedure, pass :: clear  => scattdata_histogram_clear
  end type ScattData_Histogram

  type, extends(ScattData_Base) :: ScattData_Tabular
    real(8), allocatable :: mu(:)      ! Mu bins
    real(8)              :: dmu        ! Mu spacing
    real(8), allocatable :: fmu(:,:,:) ! PDF of f(mu)
    contains
      procedure, pass :: init   => scattdata_tabular_init
      procedure, pass :: calc_f => scattdata_tabular_calc_f
      procedure, pass :: clear  => scattdata_tabular_clear
  end type ScattData_Tabular

!===============================================================================
! SCATTDATACONTAINER allocatable array for storing ScattData Objects (for angle)
!===============================================================================

  type ScattDataContainer
    class(ScattData_Base), allocatable :: obj
  end type ScattDataContainer

contains

!===============================================================================
! SCATTDATA_INIT builds the scattdata object
!===============================================================================

    subroutine scattdata_base_init(this, order, energy, mult)
      class(ScattData_Base), intent(inout) :: this        ! Object to work on
      integer, intent(in)                  :: order       ! Data Order
      real(8), intent(in)                  :: energy(:,:) ! Energy Transfer Matrix
      real(8), intent(in)                  :: mult(:,:)   ! Scatter Prod'n Matrix

      integer :: groups

      groups = size(energy, dim=1)

      allocate(this % energy(groups, groups))
      this % energy = energy
      allocate(this % mult(groups, groups))
      this % mult = mult
      allocate(this % data(order, groups, groups))
      this % data = ZERO

    end subroutine scattdata_base_init

    subroutine scattdata_legendre_init(this, order, energy, mult, coeffs)
      class(ScattData_Legendre), intent(inout) :: this   ! Object to work on
      integer, intent(in)                      :: order  ! Data Order
      real(8), intent(in)              :: energy(:,:) ! Energy Transfer Matrix
      real(8), intent(in)              :: mult(:,:)   ! Scatter Prod'n Matrix
      real(8), intent(in)              :: coeffs(:,:,:) ! Coefficients to use

      call scattdata_base_init(this, order, energy, mult)

      this % data = coeffs

    end subroutine scattdata_legendre_init

    subroutine scattdata_histogram_init(this, order, energy, mult, coeffs)
      class(ScattData_Histogram), intent(inout) :: this   ! Object to work on
      integer, intent(in)                       :: order  ! Data Order
      real(8), intent(in)              :: energy(:,:) ! Energy Transfer Matrix
      real(8), intent(in)              :: mult(:,:)   ! Scatter Prod'n Matrix
      real(8), intent(in)              :: coeffs(:,:,:) ! Coefficients to use

      integer :: imu, gin, gout, groups
      real(8) :: norm

      groups = size(energy,dim=1)

      call scattdata_base_init(this, order, energy, mult)

      allocate(this % mu(order))
      this % dmu = TWO / (real(order,8))
      this % mu(1) = -ONE
      do imu = 2, order
        this % mu(imu) = -ONE + (imu - 1) * this % dmu
      end do

      ! Best to integrate this histogram so we can avoid rejection sampling
      do gin = 1, groups
        do gout = 1, groups
          if (energy(gout,gin) > ZERO) then
            ! Integrate the histogram
            this % data(1,gout,gin) = this % dmu * coeffs(1,gout,gin)
            do imu = 2, order
              this % data(imu,gout,gin) = this % dmu * coeffs(imu,gout,gin) + &
                   this % data(imu-1,gout,gin)
            end do
            ! Now make sure integral norms to zero
            norm = this % data(order,gout,gin)
            if (norm > ZERO) then
              this % data(:,gout,gin) = this % data(:,gout,gin) / norm
            end if
          end if
        end do
      end do

    end subroutine scattdata_histogram_init

    subroutine scattdata_tabular_init(this, order, energy, mult, coeffs)
      class(ScattData_Tabular), intent(inout) :: this   ! Object to work on
      integer, intent(in)                     :: order  ! Data Order
      real(8), intent(in)              :: energy(:,:) ! Energy Transfer Matrix
      real(8), intent(in)              :: mult(:,:)   ! Scatter Prod'n Matrix
      real(8), intent(in)              :: coeffs(:,:,:) ! Coefficients to use

      integer :: imu, gin, gout, groups
      real(8) :: norm
      logical :: legendre_flag
      integer :: this_order

      if (order < 0) then
        legendre_flag = .true.
        this_order = -1 * order
      else
        legendre_flag = .false.
        this_order = order
      end if

      groups = size(energy,dim=1)

      call scattdata_base_init(this, this_order, energy, mult)

      allocate(this % mu(this_order))
      this % dmu = TWO / (real(this_order) - 1)
      do imu = 1, this_order - 1
        this % mu(imu) = -ONE + real(imu - 1) * this % dmu
      end do
      this % mu(this_order) = ONE

      ! Best to integrate this histogram so we can avoid rejection sampling
      allocate(this % fmu(this_order,groups,groups))
      do gin = 1, groups
        do gout = 1, groups
          if (energy(gout,gin) > ZERO) then
            if (legendre_flag) then
              ! Coeffs are legendre coeffs.  Need to build f(mu) then integrate
              ! and store the integral in this % data
              ! Ensure the coeffs are normalized
              norm = ONE / coeffs(1,gout,gin)
              do imu = 1, this_order
                this % fmu(imu,gout,gin) = evaluate_legendre(norm * coeffs(:,gout,gin), this % mu(imu))
                ! Force positivity
                if (this % fmu(imu,gout,gin) < ZERO) then
                  this % fmu(imu,gout,gin) = ZERO
                end if
              end do
            else
              ! Coeffs contain f(mu), put in f(mu) to save duplicate.
              this % fmu(:,gout,gin) = this % data(:,gout,gin)
            end if

            ! Re-normalize fmu for numerical integration issues and in case
            ! the negative fix-up introduced un-normalized data
            norm = ZERO
            do imu = 2, this_order
              norm = norm + HALF * this % dmu * (this % fmu(imu-1,gout,gin) + this % fmu(imu,gout,gin))
            end do
            if (norm > ZERO) then
              this % fmu(:,gout,gin) = this % fmu(:,gout,gin) / norm
            end if

            ! Now create CDF from fmu with trapezoidal rule
            this % data(1,gout,gin) = ZERO
            do imu = 2, this_order - 1
              this % data(imu,gout,gin) = this % data(imu-1,gout,gin) + &
                   HALF * this % dmu * (this % fmu(imu-1,gout,gin) + this % fmu(imu,gout,gin))
            end do
            this % data(this_order,gout,gin) = ONE
          end if
        end do
      end do

    end subroutine scattdata_tabular_init

!===============================================================================
! SCATTDATA_CLEAR resets and deallocates data in ScattData.
!===============================================================================

    subroutine scattdata_base_clear(this)
      class(ScattData_Base), intent(inout) :: this

      if (allocated(this % energy)) then
        deallocate(this % energy)
      end if

      if (allocated(this % mult)) then
        deallocate(this % mult)
      end if

      if (allocated(this % data)) then
        deallocate(this % data)
      end if

    end subroutine scattdata_base_clear

    subroutine scattdata_legendre_clear(this)
      class(ScattData_Legendre), intent(inout) :: this

      call scattdata_base_clear(this)

    end subroutine scattdata_legendre_clear

    subroutine scattdata_histogram_clear(this)
      class(ScattData_Histogram), intent(inout) :: this

      call scattdata_base_clear(this)

      if (allocated(this % mu)) then
        deallocate(this % mu)
      end if

    end subroutine scattdata_histogram_clear

    subroutine scattdata_tabular_clear(this)
      class(ScattData_Tabular), intent(inout) :: this

      call scattdata_base_clear(this)

      if (allocated(this % mu)) then
        deallocate(this % mu)
      end if

    end subroutine scattdata_tabular_clear

!===============================================================================
! SCATTDATA_*_CALC_F Calculates the value of f given mu (and gin,gout pair)
!===============================================================================

    pure function scattdata_legendre_calc_f(this, gin, gout, mu) result(f)
      class(ScattData_Legendre), intent(in) :: this ! The ScattData to evaluate
      integer, intent(in)                   :: gin  ! Incoming Energy Group
      integer, intent(in)                   :: gout ! Outgoing Energy Group
      real(8), intent(in)                   :: mu   ! Angle of interest
      real(8)                               :: f    ! Return value of f(mu)

      ! Plug mu in to the legendre expansion and go from there
      f = evaluate_legendre(this % data(:, gout, gin), mu)

    end function scattdata_legendre_calc_f

    pure function scattdata_histogram_calc_f(this, gin, gout, mu) result(f)
      class(ScattData_Histogram), intent(in) :: this ! The ScattData to evaluate
      integer, intent(in)                    :: gin  ! Incoming Energy Group
      integer, intent(in)                    :: gout ! Outgoing Energy Group
      real(8), intent(in)                    :: mu   ! Angle of interest
      real(8)                                :: f    ! Return value of f(mu)

      integer :: imu

      ! Find mu bin
      imu = floor((mu + ONE)/ this % dmu + ONE)
      ! Adjust so interpolation works on the last bin if necessary
      if (imu == size(this % data, dim=1)) then
        imu = imu - 1
      end if

      ! Use histogram interpolation to find f(mu)
      f = this % data(imu, gout, gin)

    end function scattdata_histogram_calc_f

    pure function scattdata_tabular_calc_f(this, gin, gout, mu) result(f)
      class(ScattData_Tabular), intent(in) :: this ! The ScattData to evaluate
      integer, intent(in)                  :: gin  ! Incoming Energy Group
      integer, intent(in)                  :: gout ! Outgoing Energy Group
      real(8), intent(in)                  :: mu   ! Angle of interest
      real(8)                              :: f    ! Return value of f(mu)

      integer :: imu
      real(8) :: r

      ! Find mu bin
      imu = floor((mu + ONE)/ this % dmu + ONE)
      ! Adjust so interpolation works on the last bin if necessary
      if (imu == size(this % data, dim=1)) then
        imu = imu - 1
      end if

      ! ! Now interpolate to find f(mu)
      r = (mu - this % mu(imu)) / (this % mu(imu + 1) - this % mu(imu))
      f = (ONE - r) * this % data(imu, gout, gin) + &
           r * this % data(imu + 1, gout, gin)

    end function scattdata_tabular_calc_f



end module scattdata_header