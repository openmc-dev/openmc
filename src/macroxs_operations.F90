module macroxs_operations

  use constants
  use macroxs_header,  only: MacroXS, MacroXSIso, MacroXSAngle, &
                             expand_harmonic
  use material_header, only: Material
  use math
  use nuclide_header,  only: find_angle, MaterialMacroXS, NuclideMicroXS, &
                             NuclideMG, NuclideMGContainer
  use random_lcg,      only: prn
  use scattdata_header
  use search

  implicit none

contains

!===============================================================================
! UPDATE_XS stores the xs to work with
!===============================================================================

  subroutine calculate_mgxs(this, gin, uvw, xs)
    class(MacroXS),        intent(in)    :: this
    integer,               intent(in)    :: gin         ! Incoming neutron group
    real(8),               intent(in)    :: uvw(3)      ! Incoming neutron direction
    type(MaterialMacroXS), intent(inout) :: xs

    integer :: iazi, ipol

    select type(this)
    type is (MacroXSIso)
      xs % total         = this % total(gin)
      xs % elastic       = this % scattxs(gin)
      xs % absorption    = this % absorption(gin)
      xs % nu_fission    = this % nu_fission(gin)

    type is (MacroXSAngle)
      call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)
      xs % total         = this % total(gin, iazi, ipol)
      xs % elastic       = this % scattxs(gin, iazi, ipol)
      xs % absorption    = this % absorption(gin, iazi, ipol)
      xs % nu_fission    = this % nu_fission(gin, iazi, ipol)
    end select

  end subroutine calculate_mgxs


!===============================================================================
! SAMPLE_FISSION_ENERGY acts as a templating code for macroxs_*_sample_fission_energy
!===============================================================================

  function sample_fission_energy(this, gin, uvw) result(gout)
    class(MacroXS),      intent(in) :: this   ! Data to work with
    integer, intent(in)             :: gin    ! Incoming energy group
    real(8), intent(in)             :: uvw(3) ! Particle Direction
    integer                         :: gout   ! Sampled outgoing group

    select type(this)
    type is (MacroXSIso)
      gout = macroxsiso_sample_fission_energy(this, gin, uvw)
    type is (MacroXSAngle)
      gout = macroxsangle_sample_fission_energy(this, gin, uvw)
    end select

  end function sample_fission_energy

!===============================================================================
! MACROXS_*_SAMPLE_FISSION_ENERGY samples the outgoing energy and mu from a scatter event.
! Implemented as % scatter.
!===============================================================================

  function macroxsiso_sample_fission_energy(this, gin, uvw) result(gout)
    class(MacroXSIso), intent(in) :: this   ! Data to work with
    integer, intent(in)           :: gin    ! Incoming energy group
    real(8), intent(in)           :: uvw(3) ! Particle Direction
    integer                       :: gout   ! Sampled outgoing group
    real(8) :: xi               ! Our random number
    real(8) :: prob             ! Running probability

    xi = prn()
    prob = ZERO
    gout = 0

    do while (prob < xi)
      gout = gout + 1
      prob = prob + this % chi(gout,gin)
    end do

  end function macroxsiso_sample_fission_energy

  function macroxsangle_sample_fission_energy(this, gin, uvw) result(gout)
    class(MacroXSAngle), intent(in) :: this  ! Data to work with
    integer, intent(in)             :: gin    ! Incoming energy group
    real(8), intent(in)             :: uvw(3) ! Particle Direction
    integer                         :: gout   ! Sampled outgoing group
    real(8) :: xi               ! Our random number
    real(8) :: prob             ! Running probability
    integer :: iazi, ipol

    call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)

    xi = prn()
    prob = ZERO
    gout = 0

    do while (prob < xi)
      gout = gout + 1
      prob = prob + this % chi(gout,gin,iazi,ipol)
    end do

  end function macroxsangle_sample_fission_energy

!===============================================================================
! SAMPLE_SCATTER acts as a templating code for macroxs_*_sample_scatter
!===============================================================================

  subroutine sample_scatter(this, uvw, gin, gout, mu, wgt)
    class(MacroXS), intent(in)        :: this
    real(8),             intent(in)   :: uvw(3) ! Incoming neutron direction
    integer,             intent(in)   :: gin    ! Incoming neutron group
    integer,            intent(out)   :: gout   ! Sampled outgoin group
    real(8),            intent(out)   :: mu     ! Sampled change in angle
    real(8),            intent(inout) :: wgt    ! Particle weight

    integer :: iazi, ipol ! Angular indices

    select type(this)
    type is (MacroXSIso)
      call macroxs_sample_scatter(this % scatter, gin, gout, mu, wgt)
    type is (MacroXSAngle)
      call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)
      call macroxs_sample_scatter(this % scatter(iazi,ipol) % obj,gin,gout,mu,wgt)
    end select

  end subroutine sample_scatter

!===============================================================================
! MACROXS_SAMPLE_SCATTER performs the work with ScattData to sample outgoing
! energy and change in angle.
!===============================================================================

  subroutine macroxs_sample_scatter(scatt, gin, gout, mu, wgt)
    class(ScattData), intent(in)    :: scatt ! Scattering Object to Use
    integer,          intent(in)    :: gin   ! Incoming neutron group
    integer,          intent(out)   :: gout  ! Sampled outgoin group
    real(8),          intent(out)   :: mu    ! Sampled change in angle
    real(8),          intent(inout) :: wgt   ! Particle weight

    real(8) :: xi     ! Our random number
    real(8) :: prob   ! Running probability
    integer :: imu
    real(8) :: u, f, M
    real(8) :: mu0, frac, mu1
    real(8) :: c_k, c_k1, p0, p1
    integer :: k, NP, samples

    xi = prn()
    prob = ZERO
    gout = 0

    do while (prob < xi)
      gout = gout + 1
      prob = prob + scatt % energy(gout,gin)
    end do

    select type (scatt)
    type is (ScattDataHistogram)
      xi = prn()
      if (xi < scatt % data(1,gout,gin)) then
        imu = 1
      else
        imu = binary_search(scatt % data(:,gout,gin), &
                            size(scatt % data(:,gout,gin)), xi)
      end if

      ! Randomly select a mu in this bin.
      mu = prn() * scatt % dmu + scatt % mu(imu)

    type is (ScattDataTabular)
      ! determine outgoing cosine bin
      NP = size(scatt % data(:,gout,gin))
      xi = prn()

      c_k = scatt % data(1,gout,gin)
      do k = 1, NP - 1
        c_k1 = scatt % data(k+1,gout,gin)
        if (xi < c_k1) exit
        c_k = c_k1
      end do

      ! check to make sure k is <= NP - 1
      k = min(k, NP - 1)

      p0  = scatt % fmu(k,gout,gin)
      mu0 = scatt % mu(k)
      ! Linear-linear interpolation to find mu value w/in bin.
      p1  = scatt % fmu(k+1,gout,gin)
      mu1 = scatt % mu(k+1)

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

    type is (ScattDataLegendre)
      ! Now we can sample mu using the legendre representation of the scattering
      ! kernel in data(1:this % order)

      ! Do with rejection sampling
      ! Set upper bound (instead of searching for max - though this is inefficient)
      M = 4.0_8
      samples = 0
      do
        mu = TWO * prn() - ONE
        f = scatt % calc_f(gin,gout,mu)
        if (f > ZERO) then
          u = prn() * M
          if (u <= f) then
            exit
          end if
        end if
        samples = samples + 1
        if (samples > MAX_SAMPLE) then
          ! Exit with an isotropic event.
          exit
        end if
      end do
    end select

    wgt = wgt * scatt % mult(gout,gin)

  end subroutine macroxs_sample_scatter

end module macroxs_operations