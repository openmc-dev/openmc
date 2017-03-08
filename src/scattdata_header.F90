module scattdata_header

  use algorithm,  only: binary_search
  use constants
  use error,      only: fatal_error
  use math
  use random_lcg, only: prn

  implicit none


!===============================================================================
! JAGGED1D and JAGGED2D is a type which allows for jagged 1-D or 2-D array.
!===============================================================================

  type :: Jagged2D
    real(8), allocatable :: data(:, :)
  end type Jagged2D

  type :: Jagged1D
    real(8), allocatable :: data(:)
  end type Jagged1D

!===============================================================================
! SCATTDATA contains all the data to describe the scattering energy and
! angular distribution
!===============================================================================

  type, abstract :: ScattData
    ! The data attribute of the energy, mult, and dist arrays
    ! are not necessarily 1-indexed as they instead will be allocated
    ! from a minimum outgoing group to an outgoing minimum group.
    ! Normalized p0 matrix on its own for sampling energy
    type(Jagged1D), allocatable :: energy(:) ! (Gin % data(Gout))
    ! Nu-scatter multiplication (i.e. nu-scatt/scatt)
    type(Jagged1D), allocatable :: mult(:)   ! (Gin % data(Gout))
    ! Angular distribution
    type(Jagged2D), allocatable :: dist(:)   ! (Gin % data(Order/Nmu, Gout)
    integer, allocatable :: gmin(:)     ! Minimum outgoing group
    integer, allocatable :: gmax(:)     ! Maximum outgoing group
    real(8), allocatable :: scattxs(:)  ! Isotropic Sigma_{s,g_{in}}

  contains
    procedure(scattdata_init_),   deferred :: init   ! Initializes ScattData
    procedure(scattdata_calc_f_), deferred :: calc_f ! Calculates f, given mu
    procedure(scattdata_sample_), deferred :: sample ! sample the scatter event
    procedure :: get_matrix => scattdata_get_matrix  ! Rebuild scattering matrix
  end type ScattData

  abstract interface
    subroutine scattdata_init_(this, gmin, gmax, mult, coeffs)
      import ScattData, Jagged1D, Jagged2D
      class(ScattData), intent(inout) :: this      ! Object to work with
      integer,          intent(in)    :: gmin(:)   ! Min Gout
      integer,          intent(in)    :: gmax(:)   ! Max Gout
      type(Jagged1D),   intent(in)    :: mult(:)   ! Scatter Prod'n Matrix
      type(Jagged2D),   intent(in)    :: coeffs(:) ! Coefficients to use
    end subroutine scattdata_init_

    pure function scattdata_calc_f_(this, gin, gout, mu) result(f)
      import ScattData
      class(ScattData), intent(in) :: this ! Scattering Object to work with
      integer, intent(in)          :: gin  ! Incoming Energy Group
      integer, intent(in)          :: gout ! Outgoing Energy Group
      real(8), intent(in)          :: mu   ! Angle of interest
      real(8)                      :: f    ! Return value of f(mu)

    end function scattdata_calc_f_

    subroutine scattdata_sample_(this, gin, gout, mu, wgt)
      import ScattData
      class(ScattData), intent(in)    :: this ! Scattering Object to work with
      integer,          intent(in)    :: gin  ! Incoming neutron group
      integer,          intent(out)   :: gout ! Sampled outgoin group
      real(8),          intent(out)   :: mu   ! Sampled change in angle
      real(8),          intent(inout) :: wgt  ! Particle weight
    end subroutine scattdata_sample_
  end interface

  type, extends(ScattData) :: ScattDataLegendre
    ! Maximal value for rejection sampling from rectangle
    type(Jagged1D), allocatable :: max_val(:) ! (Gin % data(Gout))
  contains
    procedure :: init => scattdatalegendre_init
    procedure :: calc_f          => scattdatalegendre_calc_f
    procedure :: sample          => scattdatalegendre_sample
  end type ScattDataLegendre

  type, extends(ScattData)      :: ScattDataHistogram
    real(8), allocatable        :: mu(:) ! Mu bins
    real(8)                     :: dmu   ! Mu spacing
    ! Histogram of f(mu) (dist has CDF)
    type(Jagged2D), allocatable :: fmu(:) ! (Gin % data(Order/Nmu x Gout)
  contains
    procedure :: init => scattdatahistogram_init
    procedure :: calc_f          => scattdatahistogram_calc_f
    procedure :: sample          => scattdatahistogram_sample
    procedure :: get_matrix      => scattdatahistogram_get_matrix
  end type ScattDataHistogram

  type, extends(ScattData)      :: ScattDataTabular
    real(8), allocatable        :: mu(:)       ! Mu bins
    real(8)                     :: dmu         ! Mu spacing
    ! PDF of f(mu) (dist has CDF)
    type(Jagged2D), allocatable :: fmu(:) ! (Gin % data(Order/Nmu x Gout)
  contains
    procedure :: init => scattdatatabular_init
    procedure :: calc_f          => scattdatatabular_calc_f
    procedure :: sample          => scattdatatabular_sample
    procedure :: get_matrix      => scattdatatabular_get_matrix
  end type ScattDataTabular

!===============================================================================
! SCATTDATACONTAINER allocatable array for storing ScattData Objects (for angle)
!===============================================================================

  type ScattDataContainer
    class(ScattData), allocatable :: obj
  end type ScattDataContainer

contains

!===============================================================================
! SCATTDATA*_INIT builds the scattdata object
!===============================================================================

    subroutine scattdata_init(this, order, gmin, gmax, energy, mult)
      class(ScattData), intent(inout) :: this      ! Object to work on
      integer,          intent(in)    :: order     ! Data Order
      integer,          intent(in)    :: gmin(:)   ! Min Gout
      integer,          intent(in)    :: gmax(:)   ! Max Gout
      type(Jagged1D),   intent(inout) :: energy(:) ! Energy Transfer Matrix
      type(Jagged1D),   intent(in)    :: mult(:)   ! Scatter Prod'n Matrix

      integer :: groups, gin
      real(8) :: norm

      groups = size(energy, dim=1)

      allocate(this % gmin(groups))
      allocate(this % gmax(groups))
      allocate(this % energy(groups))
      allocate(this % mult(groups))
      allocate(this % dist(groups))

      this % gmin = gmin
      this % gmax = gmax

      ! Set the outgoing energy PDF values
      do gin = 1, groups
        ! Make sure energy is normalized (i.e., CDF is 1)
        norm = sum(energy(gin) % data(:))
        if (norm /= ZERO) energy(gin) % data(:) = energy(gin) % data(:) / norm
        ! Set the values
        allocate(this % energy(gin) % data(gmin(gin):gmax(gin)))
        this % energy(gin) % data(:) = energy(gin) % data(:)
        allocate(this % mult(gin) % data(gmin(gin):gmax(gin)))
        this % mult(gin) % data(gmin(gin):gmax(gin)) = &
             mult(gin) % data(gmin(gin):gmax(gin))
        allocate(this % dist(gin) % data(order, gmin(gin):gmax(gin)))
        this % dist(gin) % data = ZERO
      end do
    end subroutine scattdata_init

    subroutine scattdatalegendre_init(this, gmin, gmax, mult, coeffs)
      class(ScattDataLegendre), intent(inout) :: this      ! Object to work on
      integer,                  intent(in)    :: gmin(:)   ! Min Gout
      integer,                  intent(in)    :: gmax(:)   ! Max Gout
      type(Jagged1D),           intent(in)    :: mult(:)   ! Scatter Prod'n Matrix
      type(Jagged2D),           intent(in)    :: coeffs(:) ! Coefficients to use

      real(8) :: dmu, mu, f, norm
      integer :: imu, Nmu, gout, gin, groups, order
      type(Jagged1D), allocatable :: energy(:)
      type(Jagged2D), allocatable :: matrix(:)

      groups = size(coeffs)
      order = size(coeffs(1) % data, dim=1)

      ! make a copy of coeffs that we can use to extract data and normalize
      allocate(matrix(groups))
      do gin = 1, groups
        allocate(matrix(gin) % data(order, gmin(gin):gmax(gin)))
        matrix(gin) % data = coeffs(gin) % data
      end do

      ! Get scattxs value
      allocate(this % scattxs(groups))
      ! Get this by summing the un-normalized P0 coefficient in matrix
      ! over all outgoing groups
      do gin = 1, groups
        this % scattxs(gin) = sum(matrix(gin) % data(1, :), dim=1)
      end do

      allocate(energy(groups))
      ! Build energy transfer probability matrix from data in matrix
      ! while also normalizing matrix itself (making CDF of f(mu=1)=1)
      do gin = 1, groups
        allocate(energy(gin) % data(gmin(gin):gmax(gin)))
        energy(gin) % data = ZERO
        do gout = gmin(gin), gmax(gin)
          norm = matrix(gin) % data(1, gout)
          energy(gin) % data(gout) = norm
          if (norm /= ZERO) then
            matrix(gin) % data(:, gout) = matrix(gin) % data(:, gout) / norm
          end if
        end do
      end do

      call scattdata_init(this, order, gmin, gmax, energy, mult)

      allocate(this % max_val(groups))
      ! Set dist values from matrix and initialize max_val
      do gin = 1, groups
        do gout = gmin(gin), gmax(gin)
          this % dist(gin) % data(:, gout) = matrix(gin) % data(:, gout)
        end do
        allocate(this % max_val(gin) % data(gmin(gin):gmax(gin)))
        this % max_val(gin) % data(:) = ZERO
      end do

      ! Step through the polynomial with fixed number of points to identify
      ! the maximal value.
      Nmu = 1001
      dmu = TWO / real(Nmu - 1, 8)
      do gin = 1, groups
        do gout = gmin(gin), gmax(gin)
          do imu = 1, Nmu
            ! Update mu. Do first and last seperate to avoid float errors
            if (imu == 1) then
              mu = -ONE
            else if (imu == Nmu) then
              mu = ONE
            else
              mu = -ONE + real(imu - 1, 8) * dmu
            end if
            ! Calculate probability
            f = this % calc_f(gin,gout,mu)
            ! If this is a new max, store it.
            if (f > this % max_val(gin) % data(gout)) &
                 this % max_val(gin) % data(gout) = f
          end do
          ! Finally, since we may not have caught the exact max, add 10% margin
          this % max_val(gin) % data(gout) = &
               this % max_val(gin) % data(gout) * 1.1_8
        end do
      end do
    end subroutine scattdatalegendre_init

    subroutine scattdatahistogram_init(this, gmin, gmax, mult, coeffs)
      class(ScattDataHistogram), intent(inout) :: this      ! Object to work on
      integer,                   intent(in)    :: gmin(:)   ! Min Gout
      integer,                   intent(in)    :: gmax(:)   ! Max Gout
      type(Jagged1D),            intent(in)    :: mult(:)   ! Scatter Prod'n Matrix
      type(Jagged2D),            intent(in)    :: coeffs(:) ! Coefficients to use

      integer :: imu, gin, gout, groups, order
      real(8) :: norm
      type(Jagged1D), allocatable :: energy(:)
      type(Jagged2D), allocatable :: matrix(:)

      groups = size(coeffs)
      order = size(coeffs(1) % data, dim=1)

      ! make a copy of coeffs that we can use to extract data and normalize
      allocate(matrix(groups))
      do gin = 1, groups
        allocate(matrix(gin) % data(order, gmin(gin):gmax(gin)))
        matrix(gin) % data(:, :) = coeffs(gin) % data(:, :)
      end do

      ! Get scattxs value
      allocate(this % scattxs(groups))
      ! Get this by summing the un-normalized angular distribution in matrix
      ! over all outgoing groups
      do gin = 1, groups
        this % scattxs(gin) = sum(matrix(gin) % data(:, :))
      end do

      allocate(energy(groups))
      ! Build energy transfer probability matrix from data in matrix
      ! while also normalizing matrix itself (making CDF of f(mu=1)=1)
      do gin = 1, groups
        allocate(energy(gin) % data(gmin(gin):gmax(gin)))
        do gout = gmin(gin), gmax(gin)
          norm = sum(matrix(gin) % data(:, gout))
          energy(gin) % data(gout) = norm
          if (norm /= ZERO) then
            matrix(gin) % data(:, gout) = matrix(gin) % data(:, gout) / norm
          end if
        end do
      end do

      call scattdata_init(this, order, gmin, gmax, energy, mult)

      allocate(this % mu(order))
      this % dmu = TWO / real(order, 8)
      this % mu(1) = -ONE
      do imu = 2, order
        this % mu(imu) = -ONE + real(imu - 1, 8) * this % dmu
      end do

      ! Integrate this histogram so we can avoid rejection sampling while
      ! also saving the original histogram in fmu
      allocate(this % fmu(groups))
      do gin = 1, groups
        allocate(this % fmu(gin) % data(order, gmin(gin):gmax(gin)))
        do gout = gmin(gin), gmax(gin)
          ! Store the histogram
          this % fmu(gin) % data(:, gout) = matrix(gin) % data(:, gout)
          ! Integrate the histogram
          this % dist(gin) % data(1, gout) = &
               this % dmu * matrix(gin) % data(1, gout)
          do imu = 2, order
            this % dist(gin) % data(imu, gout) = &
                 this % dmu * matrix(gin) % data(imu, gout) + &
                 this % dist(gin) % data(imu - 1, gout)
          end do

          ! Normalize the integral to unity
          norm = this % dist(gin) % data(order, gout)
          if (norm > ZERO) then
            this % fmu(gin) % data(:, gout) = &
                 this % fmu(gin) % data(:, gout) / norm
            this % dist(gin) % data(:, gout) = &
                 this % dist(gin) % data(:, gout) / norm
          end if
        end do
      end do

    end subroutine scattdatahistogram_init

    subroutine scattdatatabular_init(this, gmin, gmax, mult, coeffs)
      class(ScattDataTabular), intent(inout) :: this      ! Object to work on
      integer,                 intent(in)    :: gmin(:)   ! Min Gout
      integer,                 intent(in)    :: gmax(:)   ! Max Gout
      type(Jagged1D),          intent(in)    :: mult(:)   ! Scatter Prod'n Matrix
      type(Jagged2D),          intent(in)    :: coeffs(:) ! Coefficients to use

      integer :: imu, gin, gout, groups, order
      real(8) :: norm
      type(Jagged1D), allocatable :: energy(:)
      type(Jagged2D), allocatable :: matrix(:)

      groups = size(coeffs)
      order = size(coeffs(1) % data, dim=1)

      ! make a copy of coeffs that we can use to extract data and normalize
      allocate(matrix(groups))
      do gin = 1, groups
        allocate(matrix(gin) % data(order, gmin(gin):gmax(gin)))
        matrix(gin) % data = coeffs(gin) % data
      end do

      ! Build the angular distribution mu values
      allocate(this % mu(order))
      this % dmu = TWO / real(order - 1, 8)
      this % mu(1) = -ONE
      do imu = 2, order - 1
        this % mu(imu) = -ONE + real(imu - 1, 8) * this % dmu
      end do
      this % mu(order) = ONE

      ! Get scattxs
      allocate(this % scattxs(groups))
      ! Get this by integrating the scattering distribution over all mu points
      ! and then combining over all outgoing groups
      ! over all outgoing groups
      do gin = 1, groups
        norm = ZERO
        do gout = gmin(gin), gmax(gin)
          do imu = 2, order
            norm = norm + HALF * this % dmu * &
                 (matrix(gin) % data(imu - 1, gout) + &
                  matrix(gin) % data(imu, gout))
          end do
        end do
        this % scattxs(gin) = norm
      end do

      allocate(energy(groups))
      ! Build energy transfer probability matrix from data in matrix
      do gin = 1, groups
        allocate(energy(gin) % data(gmin(gin):gmax(gin)))
        do gout = gmin(gin), gmax(gin)
          norm = ZERO
          do imu = 2, order
            norm = norm + HALF * this % dmu * &
                 (matrix(gin) % data(imu - 1, gout) + &
                  matrix(gin) % data(imu, gout))
          end do
          energy(gin) % data(gout) = norm
        end do
      end do
      call scattdata_init(this, order, gmin, gmax, energy, mult)

      ! Calculate f(mu) and integrate it so we can avoid rejection sampling
      allocate(this % fmu(groups))
      do gin = 1, groups
        allocate(this % fmu(gin) % data(order, gmin(gin):gmax(gin)))
        do gout = gmin(gin), gmax(gin)
          ! Coeffs contain f(mu), put in f(mu) as that is where the
          ! PDF lives
          this % fmu(gin) % data(:, gout) = matrix(gin) % data(:, gout)

          ! Force positivity
          do imu = 1, order
            if (this % fmu(gin) % data(imu, gout) < ZERO) then
              this % fmu(gin) % data(imu, gout) = ZERO
            end if
          end do

          ! Re-normalize fmu for numerical integration issues and in case
          ! the negative fix-up introduced un-normalized data while
          ! accruing the CDF
          norm = ZERO
          do imu = 2, order
            norm = norm + HALF * this % dmu * &
                 (this % fmu(gin) % data(imu - 1, gout) + &
                  this % fmu(gin) % data(imu, gout))
            this % dist(gin) % data(imu, gout) = norm
          end do
          if (norm > ZERO) then
            this % fmu(gin) % data(:, gout) = &
                 this % fmu(gin) % data(:, gout) / norm
            this % dist(gin) % data(:, gout) = &
                 this % dist(gin) % data(:, gout) / norm
          end if
        end do
      end do
    end subroutine scattdatatabular_init

!===============================================================================
! SCATTDATA_*_CALC_F Calculates the value of f given mu (and gin,gout pair)
!===============================================================================

    pure function scattdatalegendre_calc_f(this, gin, gout, mu) result(f)
      class(ScattDataLegendre), intent(in) :: this ! The ScattData to evaluate
      integer, intent(in)                  :: gin  ! Incoming Energy Group
      integer, intent(in)                  :: gout ! Outgoing Energy Group
      real(8), intent(in)                  :: mu   ! Angle of interest
      real(8)                              :: f    ! Return value of f(mu)

      ! Plug mu in to the legendre expansion and go from there
      if (gout < this % gmin(gin) .or. gout > this % gmax(gin)) then
        f = ZERO
      else
        f = evaluate_legendre(this % dist(gin) % data(:, gout), mu)
      end if

    end function scattdatalegendre_calc_f

    pure function scattdatahistogram_calc_f(this, gin, gout, mu) result(f)
      class(ScattDataHistogram), intent(in) :: this ! The ScattData to evaluate
      integer, intent(in)                   :: gin  ! Incoming Energy Group
      integer, intent(in)                   :: gout ! Outgoing Energy Group
      real(8), intent(in)                   :: mu   ! Angle of interest
      real(8)                               :: f    ! Return value of f(mu)

      integer :: imu

      if (gout < this % gmin(gin) .or. gout > this % gmax(gin)) then
        f = ZERO
      else
        ! Find mu bin
        if (mu == ONE) then
          imu = size(this % fmu(gin) % data, dim=1)
        else
          imu = floor((mu + ONE) / this % dmu + ONE)
        end if

        f = this % fmu(gin) % data(imu, gout)
      end if

    end function scattdatahistogram_calc_f

    pure function scattdatatabular_calc_f(this, gin, gout, mu) result(f)
      class(ScattDataTabular), intent(in) :: this ! The ScattData to evaluate
      integer, intent(in)                 :: gin  ! Incoming Energy Group
      integer, intent(in)                 :: gout ! Outgoing Energy Group
      real(8), intent(in)                 :: mu   ! Angle of interest
      real(8)                             :: f    ! Return value of f(mu)

      integer :: imu
      real(8) :: r

      if (gout < this % gmin(gin) .or. gout > this % gmax(gin)) then
        f = ZERO
      else
        ! Find mu bin
        if (mu == ONE) then
          imu = size(this % fmu(gin) % data, dim=1) - 1
        else
          imu = floor((mu + ONE) / this % dmu + ONE)
        end if

        ! Now interpolate to find f(mu)
        r = (mu - this % mu(imu)) / (this % mu(imu + 1) - this % mu(imu))
        f = (ONE - r) * this % fmu(gin) % data(imu, gout) + &
             r * this % fmu(gin) % data(imu + 1, gout)
      end if

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
    gout = this % gmin(gin)
    prob = this % energy(gin) % data(gout)

    do while ((prob < xi) .and. (gout < this % gmax(gin)))
      gout = gout + 1
      prob = prob + this % energy(gin) % data(gout)
    end do

    ! Now we can sample mu using the legendre representation of the scattering
    ! kernel in data(1:this % order)

    ! Do with rejection sampling from a rectangular bounding box
    ! Set maximal value
    M = this % max_val(gin) % data(gout)
    samples = 0
    do
      mu = TWO * prn() - ONE
      f = this % calc_f(gin, gout, mu)
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

    wgt = wgt * this % mult(gin) % data(gout)

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
    gout = this % gmin(gin)
    prob = this % energy(gin) % data(gout)

    do while ((prob < xi) .and. (gout < this % gmax(gin)))
      gout = gout + 1
      prob = prob + this % energy(gin) % data(gout)
    end do

    xi = prn()
    if (xi < this % dist(gin) % data(1, gout)) then
      imu = 1
    else
      imu = binary_search(this % dist(gin) % data(:, gout), &
                          size(this % dist(gin) % data(:, gout)), xi) + 1
    end if

    ! Randomly select a mu in this bin.
    mu = prn() * this % dmu + this % mu(imu)

    wgt = wgt * this % mult(gin) % data(gout)

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
    gout = this % gmin(gin)
    prob = this % energy(gin) % data(gout)

    do while ((prob < xi) .and. (gout < this % gmax(gin)))
      gout = gout + 1
      prob = prob + this % energy(gin) % data(gout)
    end do

    ! determine outgoing cosine bin
    NP = size(this % dist(gin) % data(:, gout))
    xi = prn()

    c_k = this % dist(gin) % data(1, gout)
    do k = 1, NP - 1
      c_k1 = this % dist(gin) % data(k + 1, gout)
      if (xi < c_k1) exit
      c_k = c_k1
    end do

    ! check to make sure k is <= NP - 1
    k = min(k, NP - 1)

    p0  = this % fmu(gin) % data(k, gout)
    mu0 = this % mu(k)
    ! Linear-linear interpolation to find mu value w/in bin.
    p1  = this % fmu(gin) % data(k + 1, gout)
    mu1 = this % mu(k + 1)

    if (p0 == p1) then
      mu = mu0 + (xi - c_k) / p0
    else
      frac = (p1 - p0) / (mu1 - mu0)
      mu = mu0 + &
           (sqrt(max(ZERO, p0 * p0 + TWO * frac * (xi - c_k))) - p0) / frac
    end if

    if (mu <= -ONE) then
      mu = -ONE
    else if (mu >= ONE) then
      mu = ONE
    end if

    wgt = wgt * this % mult(gin) % data(gout)

  end subroutine scattdatatabular_sample

!===============================================================================
! SCATTDATA*_GET_MATRIX Reproduces the original scattering matrix (densely)
! using ScattData's information of fmu/dist, energy, and scattxs
!===============================================================================

    subroutine scattdata_get_matrix(this, req_order, matrix)
      class(ScattData), intent(in) :: this      ! Scattering Object to work with
      integer, intent(in)          :: req_order ! Requested order of matrix
      type(Jagged2D), allocatable, intent(inout) :: matrix(:) ! Resultant matrix just built

      integer :: order, groups, gin, gout

      groups = size(this % energy)
      ! Set gin and gout for getting the order
      order = min(req_order, size(this % dist(1) % data, dim=1))

      if (allocated(matrix)) deallocate(matrix)
      allocate(matrix(groups))
      ! Initialize to 0; this way the zero entries in the dense matrix dont
      ! need to be explicitly set, requiring a significant increase in the
      ! lines of code.
      do gin = 1, groups
        allocate(matrix(gin) % data(order, groups))
        do gout = this % gmin(gin), this % gmax(gin)
          matrix(gin) % data(:, gout) = this % scattxs(gin) * &
               this % energy(gin) % data(gout) * &
               this % dist(gin) % data(1:order, gout)
        end do
      end do
    end subroutine scattdata_get_matrix

    subroutine scattdatahistogram_get_matrix(this, req_order, matrix)
      class(ScattDataHistogram), intent(in) :: this      ! Scattering Object to work with
      integer, intent(in)                   :: req_order ! Requested order of matrix
      type(Jagged2D), allocatable, intent(inout) :: matrix(:) ! Resultant matrix just built

      integer :: order, groups, gin, gout

      groups = size(this % energy)
      order = min(req_order, size(this % dist(1) % data, dim=1))

      if (allocated(matrix)) deallocate(matrix)
      allocate(matrix(groups))
      ! Initialize to 0; this way the zero entries in the dense matrix dont
      ! need to be explicitly set, requiring a significant increase in the
      ! lines of code.
      do gin = 1, groups
        allocate(matrix(gin) % data(order, groups))
        do gout = this % gmin(gin), this % gmax(gin)
          matrix(gin) % data(:, gout) = this % scattxs(gin) * &
               this % energy(gin) % data(gout) * &
               this % fmu(gin) % data(1:order, gout)
        end do
      end do
    end subroutine scattdatahistogram_get_matrix

    subroutine scattdatatabular_get_matrix(this, req_order, matrix)
      class(ScattDataTabular), intent(in) :: this      ! Scattering Object to work with
      integer, intent(in)                 :: req_order ! Requested order of matrix
      type(Jagged2D), allocatable, intent(inout) :: matrix(:) ! Resultant matrix just built

      integer :: order, groups, gin, gout

      groups = size(this % energy)
      order = min(req_order, size(this % dist(1) % data, dim=1))

      if (allocated(matrix)) deallocate(matrix)
      allocate(matrix(groups))
      ! Initialize to 0; this way the zero entries in the dense matrix dont
      ! need to be explicitly set, requiring a significant increase in the
      ! lines of code.
      do gin = 1, groups
        allocate(matrix(gin) % data(order, groups))
        do gout = this % gmin(gin), this % gmax(gin)
          matrix(gin) % data(:, gout) = this % scattxs(gin) * &
               this % energy(gin) % data(gout) * &
               this % fmu(gin) % data(1:order, gout)
        end do
      end do
    end subroutine scattdatatabular_get_matrix

!===============================================================================
! JAGGED_FROM_DENSE_*D Creates a jagged array from a sparse dense matrix.
! The user can supply a key which indicates the values to remove, but the
! default is ZERO
!===============================================================================

    subroutine jagged_from_dense_1D(dense, jagged, lo_bounds_, hi_bounds_, key_)
      real(8), intent(in)                        :: dense(:, :)
      type(Jagged1D), allocatable, intent(inout) :: jagged(:)
      real(8), intent(in), optional              :: key_
      integer, intent(inout), allocatable, optional :: lo_bounds_(:)
      integer, intent(inout), allocatable, optional :: hi_bounds_(:)

      real(8) :: key
      integer :: i, jmin, jmax
      integer, allocatable :: lo_bounds(:), hi_bounds(:)

      if (present(key_)) then
        key = key_
      else
        key = ZERO
      end if

      allocate(lo_bounds(size(dense, dim=2)))
      allocate(hi_bounds(size(dense, dim=2)))

      if (allocated(jagged)) deallocate(jagged)
      allocate(jagged(size(dense, dim=2)))
      do i = 1, size(dense, dim=2)
        ! Find the min and max j values
        do jmin = 1, size(dense, dim=1)
          if (dense(jmin, i) /= key) exit
        end do
        do jmax = size(dense, dim=1), 1, -1
          if (dense(jmax, i) /= key) exit
        end do
        ! Treat the case of all values matching the key
        if (jmin > jmax) then
          jmin = i
          jmax = i
        end if

        ! Now store the jagged row
        allocate(jagged(i) % data(jmin:jmax))
        jagged(i) % data(jmin:jmax) = dense(jmin:jmax, i)

        lo_bounds(i) = jmin
        hi_bounds(i) = jmax
      end do

      if (present(lo_bounds_)) then
        if (allocated(lo_bounds_)) deallocate(lo_bounds_)
        allocate(lo_bounds_(size(dense, dim=2)))
        lo_bounds_ = lo_bounds
      end if
      if (present(hi_bounds_)) then
        if (allocated(hi_bounds_)) deallocate(hi_bounds_)
        allocate(hi_bounds_(size(dense, dim=2)))
        hi_bounds_ = hi_bounds
      end if

    end subroutine jagged_from_dense_1D

    subroutine jagged_from_dense_2D(dense, jagged, lo_bounds_, hi_bounds_, key_)
      real(8), intent(in)                        :: dense(:, :, :)
      type(Jagged2D), allocatable, intent(inout) :: jagged(:)
      real(8), intent(in), optional              :: key_
      integer, intent(inout), allocatable, optional :: lo_bounds_(:)
      integer, intent(inout), allocatable, optional :: hi_bounds_(:)

      real(8) :: key
      integer :: i, jmin, jmax
      integer, allocatable :: lo_bounds(:), hi_bounds(:)

      if (present(key_)) then
        key = key_
      else
        key = ZERO
      end if

      allocate(lo_bounds(size(dense, dim=3)))
      allocate(hi_bounds(size(dense, dim=3)))

      if (allocated(jagged)) deallocate(jagged)
      allocate(jagged(size(dense, dim=3)))
      do i = 1, size(dense, dim=3)
        ! Find the min and max j values
        do jmin = 1, size(dense, dim=2)
          if (any(dense(:, jmin, i) /= key)) exit
        end do
        do jmax = size(dense, dim=2), 1, -1
          if (any(dense(:, jmax, i) /= key)) exit
        end do
        ! Treat the case of all values matching the key
        if (jmin > jmax) then
          jmin = i
          jmax = i
        end if

        ! Now store the jagged row
        allocate(jagged(i) % data(size(dense, dim=1), jmin:jmax))
        jagged(i) % data(:, jmin:jmax) = dense(:, jmin:jmax, i)

        lo_bounds(i) = jmin
        hi_bounds(i) = jmax
      end do

      if (present(lo_bounds_)) then
        if (allocated(lo_bounds_)) deallocate(lo_bounds_)
        allocate(lo_bounds_(size(dense, dim=3)))
        lo_bounds_ = lo_bounds
      end if
      if (present(hi_bounds_)) then
        if (allocated(hi_bounds_)) deallocate(hi_bounds_)
        allocate(hi_bounds_(size(dense, dim=3)))
        hi_bounds_ = hi_bounds
      end if

    end subroutine jagged_from_dense_2D

end module scattdata_header
