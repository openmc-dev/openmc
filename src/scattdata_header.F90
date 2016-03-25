module scattdata_header

  use constants
  use error,      only: fatal_error
  use math
  use random_lcg, only: prn
  use search,     only: binary_search

  implicit none

!===============================================================================
! SCATTDATA contains all the data to describe the scattering energy and
! angular distribution
!===============================================================================

  type, abstract :: ScattData
    ! p0 matrix on its own for sampling energy
    real(8), allocatable :: energy(:,:) ! (Gout x Gin)
    real(8), allocatable :: mult(:,:)   ! (Gout x Gin)
    real(8), allocatable :: data(:,:,:) ! (Order/Nmu x Gout x Gin)

  contains
    procedure(scattdata_init_), deferred   :: init   ! Initializes ScattData
    procedure(scattdata_calc_f_), deferred :: calc_f ! Calculates f, given mu
    procedure(scattdata_sample_), deferred :: sample ! sample the scatter event
  end type ScattData

  abstract interface
    subroutine scattdata_init_(this, order, energy, mult, coeffs)
      import ScattData
      class(ScattData), intent(inout) :: this          ! Object to work on
      integer, intent(in)             :: order         ! Data Order
      real(8), intent(in)             :: energy(:,:)   ! Energy Transfer Matrix
      real(8), intent(in)             :: mult(:,:)     ! Scatter Prod'n Matrix
      real(8), intent(in)             :: coeffs(:,:,:) ! Coefficients to use
    end subroutine scattdata_init_

    pure function scattdata_calc_f_(this, gin, gout, mu) result(f)
      import ScattData
      class(ScattData), intent(in) :: this ! The ScattData to evaluate
      integer, intent(in)          :: gin  ! Incoming Energy Group
      integer, intent(in)          :: gout ! Outgoing Energy Group
      real(8), intent(in)          :: mu   ! Angle of interest
      real(8)                      :: f    ! Return value of f(mu)

    end function scattdata_calc_f_

    subroutine scattdata_sample_(this, gin, gout, mu, wgt)
      import ScattData
      class(ScattData), intent(in)    :: this ! Scattering Object to Use
      integer,          intent(in)    :: gin  ! Incoming neutron group
      integer,          intent(out)   :: gout ! Sampled outgoin group
      real(8),          intent(out)   :: mu   ! Sampled change in angle
      real(8),          intent(inout) :: wgt  ! Particle weight
    end subroutine scattdata_sample_
  end interface

  type, extends(ScattData) :: ScattDataLegendre
    ! Maximal value for rejection sampling from rectangle
    real(8), allocatable  :: max_val(:,:)
  contains
    procedure :: init   => scattdatalegendre_init
    procedure :: calc_f => scattdatalegendre_calc_f
    procedure :: sample => scattdatalegendre_sample
  end type ScattDataLegendre

  type, extends(ScattData) :: ScattDataHistogram
    real(8), allocatable   :: mu(:) ! Mu bins
    real(8)                :: dmu   ! Mu spacing
  contains
    procedure :: init   => scattdatahistogram_init
    procedure :: calc_f => scattdatahistogram_calc_f
    procedure :: sample => scattdatahistogram_sample
  end type ScattDataHistogram

  type, extends(ScattData) :: ScattDataTabular
    real(8), allocatable   :: mu(:)      ! Mu bins
    real(8)                :: dmu        ! Mu spacing
    real(8), allocatable   :: fmu(:,:,:) ! PDF of f(mu)
  contains
    procedure :: init   => scattdatatabular_init
    procedure :: calc_f => scattdatatabular_calc_f
    procedure :: sample => scattdatatabular_sample
  end type ScattDataTabular

!===============================================================================
! SCATTDATACONTAINER allocatable array for storing ScattData Objects (for angle)
!===============================================================================

  type ScattDataContainer
    class(ScattData), allocatable :: obj
  end type ScattDataContainer

contains

!===============================================================================
! SCATTDATA_INIT builds the scattdata object
!===============================================================================

    subroutine scattdata_init(this, order, energy, mult)
      class(ScattData), intent(inout) :: this        ! Object to work on
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

    end subroutine scattdata_init

    subroutine scattdatalegendre_init(this, order, energy, mult, coeffs)
      class(ScattDataLegendre), intent(inout) :: this   ! Object to work on
      integer, intent(in)                     :: order  ! Data Order
      real(8), intent(in)              :: energy(:,:)   ! Energy Transfer Matrix
      real(8), intent(in)              :: mult(:,:)     ! Scatter Prod'n Matrix
      real(8), intent(in)              :: coeffs(:,:,:) ! Coefficients to use

      real(8) :: dmu, mu, f
      integer :: imu, Nmu, gout, gin, groups

      call scattdata_init(this, order, energy, mult)

      this % data = coeffs

      groups = size(this % energy,dim=1)

      allocate(this % max_val(groups, groups))
      this % max_val = ZERO
      ! Step through the polynomial with fixed number of points to identify
      ! the maximal value.
      Nmu = 1001
      dmu = TWO / real(Nmu,8)
      do imu = 1, Nmu
        ! Update mu. Do first and last seperate to avoid float errors
        if (imu == 1) then
          mu = -ONE
        else if (imu == Nmu) then
          mu = ONE
        end if
        mu = -ONE + real(imu - 1,8) * dmu
        do gin = 1, groups
          do gout = 1, groups
            ! Calculate probability
            f = this % calc_f(gin,gout,mu)
            ! If this is a new max, store it.
            if (f > this % max_val(gout,gin)) this % max_val(gout,gin) = f
          end do
        end do
      end do

      ! Finally, since we may not have caught the exact max, add 10% margin
      this % max_val = this % max_val * 1.1_8

    end subroutine scattdatalegendre_init

    subroutine scattdatahistogram_init(this, order, energy, mult, coeffs)
      class(ScattDataHistogram), intent(inout) :: this   ! Object to work on
      integer, intent(in)                       :: order  ! Data Order
      real(8), intent(in)              :: energy(:,:) ! Energy Transfer Matrix
      real(8), intent(in)              :: mult(:,:)   ! Scatter Prod'n Matrix
      real(8), intent(in)              :: coeffs(:,:,:) ! Coefficients to use

      integer :: imu, gin, gout, groups
      real(8) :: norm

      groups = size(energy,dim=1)

      call scattdata_init(this, order, energy, mult)

      allocate(this % mu(order))
      this % dmu = TWO / real(order,8)
      this % mu(1) = -ONE
      do imu = 2, order
        this % mu(imu) = -ONE + real(imu - 1,8) * this % dmu
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

    end subroutine scattdatahistogram_init

    subroutine scattdatatabular_init(this, order, energy, mult, coeffs)
      class(ScattDataTabular), intent(inout) :: this   ! Object to work on
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

      call scattdata_init(this, this_order, energy, mult)

      allocate(this % mu(this_order))
      this % dmu = TWO / real(this_order - 1)
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

    end subroutine scattdatatabular_init

!===============================================================================
! SCATTDATA_*_CALC_F Calculates the value of f given mu (and gin,gout pair)
!===============================================================================

    pure function scattdatalegendre_calc_f(this, gin, gout, mu) result(f)
      class(ScattDataLegendre), intent(in) :: this ! The ScattData to evaluate
      integer, intent(in)                   :: gin  ! Incoming Energy Group
      integer, intent(in)                   :: gout ! Outgoing Energy Group
      real(8), intent(in)                   :: mu   ! Angle of interest
      real(8)                               :: f    ! Return value of f(mu)

      ! Plug mu in to the legendre expansion and go from there
      f = evaluate_legendre(this % data(:, gout, gin), mu)

    end function scattdatalegendre_calc_f

    pure function scattdatahistogram_calc_f(this, gin, gout, mu) result(f)
      class(ScattDataHistogram), intent(in) :: this ! The ScattData to evaluate
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

    end function scattdatahistogram_calc_f

    pure function scattdatatabular_calc_f(this, gin, gout, mu) result(f)
      class(ScattDataTabular), intent(in) :: this ! The ScattData to evaluate
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

    end function scattdatatabular_calc_f

!===============================================================================
! SCATTDATA*_SCATTER Samples the outgoing energy and change in angle.
!===============================================================================

  subroutine scattdatalegendre_sample(this, gin, gout, mu, wgt)
    class(ScattDataLegendre), intent(in)    :: this ! Scattering object to use
    integer,                  intent(in)    :: gin  ! Incoming neutron group
    integer,                  intent(out)   :: gout ! Sampled outgoin group
    real(8),                  intent(out)   :: mu   ! Sampled change in angle
    real(8),                  intent(inout) :: wgt  ! Particle weight

    real(8) :: xi     ! Our random number
    real(8) :: prob   ! Running probability
    real(8) :: u, f, M
    integer :: samples

    xi = prn()
    prob = ZERO
    gout = 0

    do while (prob < xi)
      gout = gout + 1
      prob = prob + this % energy(gout,gin)
    end do

    ! Now we can sample mu using the legendre representation of the thisering
    ! kernel in data(1:this % order)

    ! Do with rejection sampling
    ! Set maximal value
    M = this % max_val(gout,gin)
    samples = 0
    do
      mu = TWO * prn() - ONE
      f = this % calc_f(gin,gout,mu)
      if (f > ZERO) then
        u = prn() * M
        if (u <= f) then
          exit
        end if
      end if
      samples = samples + 1
      if (samples > MAX_SAMPLE) then
        call fatal_error("Maximum number of Legendre expansion samples reached!")
      end if
    end do

    wgt = wgt * this % mult(gout,gin)

  end subroutine scattdatalegendre_sample

  subroutine scattdatahistogram_sample(this, gin, gout, mu, wgt)
    class(ScattDataHistogram), intent(in)    :: this ! Scattering object to use
    integer,                   intent(in)    :: gin  ! Incoming neutron group
    integer,                   intent(out)   :: gout ! Sampled outgoin group
    real(8),                   intent(out)   :: mu   ! Sampled change in angle
    real(8),                   intent(inout) :: wgt  ! Particle weight

    real(8) :: xi     ! Our random number
    real(8) :: prob   ! Running probability
    integer :: imu

    xi = prn()
    prob = ZERO
    gout = 0

    do while (prob < xi)
      gout = gout + 1
      prob = prob + this % energy(gout,gin)
    end do

    xi = prn()
    if (xi < this % data(1,gout,gin)) then
      imu = 1
    else
      imu = binary_search(this % data(:,gout,gin), &
                          size(this % data(:,gout,gin)), xi)
    end if

    ! Randomly select a mu in this bin.
    mu = prn() * this % dmu + this % mu(imu)

    wgt = wgt * this % mult(gout,gin)

  end subroutine scattdatahistogram_sample

  subroutine scattdatatabular_sample(this, gin, gout, mu, wgt)
    class(ScattDataTabular), intent(in)    :: this ! Scattering object to use
    integer,                  intent(in)    :: gin  ! Incoming neutron group
    integer,                  intent(out)   :: gout ! Sampled outgoin group
    real(8),                  intent(out)   :: mu   ! Sampled change in angle
    real(8),                  intent(inout) :: wgt  ! Particle weight

    real(8) :: xi     ! Our random number
    real(8) :: prob   ! Running probability
    real(8) :: mu0, frac, mu1
    real(8) :: c_k, c_k1, p0, p1
    integer :: k, NP

    xi = prn()
    prob = ZERO
    gout = 0

    do while (prob < xi)
      gout = gout + 1
      prob = prob + this % energy(gout,gin)
    end do

    ! determine outgoing cosine bin
    NP = size(this % data(:,gout,gin))
    xi = prn()

    c_k = this % data(1,gout,gin)
    do k = 1, NP - 1
      c_k1 = this % data(k+1,gout,gin)
      if (xi < c_k1) exit
      c_k = c_k1
    end do

    ! check to make sure k is <= NP - 1
    k = min(k, NP - 1)

    p0  = this % fmu(k,gout,gin)
    mu0 = this % mu(k)
    ! Linear-linear interpolation to find mu value w/in bin.
    p1  = this % fmu(k+1,gout,gin)
    mu1 = this % mu(k+1)

    frac = (p1 - p0)/(mu1 - mu0)

    if (frac == ZERO) then
      mu = mu0 + (xi - c_k)/p0
    else
      mu = mu0 + (sqrt(max(ZERO, p0*p0 + TWO*frac*(xi - c_k))) - p0)/frac
    end if

    if (mu <= -ONE) then
      mu = -ONE
    else if (mu >= ONE) then
      mu = ONE
    end if

    wgt = wgt * this % mult(gout,gin)

  end subroutine scattdatatabular_sample

end module scattdata_header
