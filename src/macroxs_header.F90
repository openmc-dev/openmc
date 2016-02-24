module macroxs_header

  use constants,       only: MAX_FILE_LEN, ZERO, ONE, TWO, PI
  use list_header,     only: ListInt
  use material_header, only: material
  use math,            only: calc_pn, calc_rn, expand_harmonic, find_angle
  use nuclide_header
  use random_lcg,      only: prn
  use scattdata_header

  implicit none

!===============================================================================
! MACROXS_* contains cached macroscopic cross sections for the material a
! particle is traveling through
!===============================================================================

  type, abstract :: MacroXS
    ! Data Order
    integer :: order

  contains
    procedure(macroxs_init_),   deferred :: init   ! initializes object
    procedure(macroxs_get_xs_), deferred :: get_xs ! Return xs
    ! Sample the outgoing energy from a fission event
    procedure(macroxs_sample_fission_), deferred :: sample_fission_energy
    ! Sample the outgoing energy and angle from a scatter event
    procedure(macroxs_sample_scatter_), deferred :: sample_scatter
    ! Calculate the material specific MGXS data from the nuclides
    procedure(macroxs_calculate_xs_), deferred   :: calculate_xs
  end type MacroXS

  abstract interface
    subroutine macroxs_init_(this, mat, nuclides, groups, get_kfiss, get_fiss, &
                             max_order, scatt_type, legendre_mu_points, &
                             error_code, error_text)
      import MacroXS, Material, NuclideMGContainer, MAX_LINE_LEN
      class(MacroXS), intent(inout)        :: this ! The MacroXS to initialize
      type(Material), pointer, intent(in)  :: mat  ! base material
      type(NuclideMGContainer), intent(in) :: nuclides(:) ! List of nuclides to harvest from
      integer, intent(in)                  :: groups ! Number of E groups
      logical, intent(in)                  :: get_kfiss ! Should we get kfiss data?
      logical, intent(in)                  :: get_fiss ! Should we get fiss data?
      integer, intent(in)                  :: max_order ! Maximum requested order
      integer, intent(in)                  :: scatt_type ! Legendre or Tabular Scatt?
      integer, intent(in)                  :: legendre_mu_points ! Treat as Leg or Tabular?
      integer, intent(inout)               :: error_code  ! Code signifying error
      character(MAX_LINE_LEN), intent(inout) :: error_text ! Error message to print
    end subroutine macroxs_init_

    function macroxs_get_xs_(this, g, xstype, gout, uvw) result(xs)
      import MacroXS
      class(MacroXS), intent(in)      :: this   ! The MacroXS to initialize
      integer, intent(in)             :: g      ! Incoming Energy group
      character(*) , intent(in)       :: xstype ! Cross Section Type
      integer, optional, intent(in)   :: gout   ! Outgoing Energy group
      real(8), optional, intent(in)   :: uvw(3) ! Requested Angle
      real(8)                         :: xs     ! Resultant xs
    end function macroxs_get_xs_

    function macroxs_sample_fission_(this, gin, uvw) result(gout)
      import MacroXS
      class(MacroXS), intent(in) :: this   ! Data to work with
      integer, intent(in)        :: gin    ! Incoming energy group
      real(8), intent(in)        :: uvw(3) ! Particle Direction
      integer                    :: gout   ! Sampled outgoing group

    end function macroxs_sample_fission_

    subroutine macroxs_sample_scatter_(this, uvw, gin, gout, mu, wgt)
      import MacroXS
      class(MacroXS), intent(in)    :: this
      real(8),        intent(in)    :: uvw(3) ! Incoming neutron direction
      integer,        intent(in)    :: gin    ! Incoming neutron group
      integer,        intent(out)   :: gout   ! Sampled outgoin group
      real(8),        intent(out)   :: mu     ! Sampled change in angle
      real(8),        intent(inout) :: wgt    ! Particle weight
    end subroutine macroxs_sample_scatter_

    subroutine macroxs_calculate_xs_(this, gin, uvw, xs)
      import MacroXS, MaterialMacroXS
      class(MacroXS),        intent(in)    :: this
      integer,               intent(in)    :: gin         ! Incoming neutron group
      real(8),               intent(in)    :: uvw(3)      ! Incoming neutron direction
      type(MaterialMacroXS), intent(inout) :: xs
    end subroutine macroxs_calculate_xs_
  end interface

  type, extends(MacroXS) :: MacroXSIso
    ! Microscopic cross sections
    real(8), allocatable :: total(:)      ! total cross section
    real(8), allocatable :: absorption(:) ! absorption cross section
    class(ScattData), allocatable :: scatter ! scattering information
    real(8), allocatable :: nu_fission(:) ! nu-fission
    real(8), allocatable :: k_fission(:)  ! kappa-fission
    real(8), allocatable :: fission(:)    ! fission x/s
    real(8), allocatable :: scattxs(:)    ! scattering xs
    real(8), allocatable :: chi(:,:)      ! fission spectra

  contains
    procedure :: init        => macroxsiso_init   ! inits object
    procedure :: get_xs      => macroxsiso_get_xs ! Returns xs
    procedure :: sample_fission_energy => macroxsiso_sample_fission_energy
    procedure :: sample_scatter => macroxsiso_sample_scatter
    procedure :: calculate_xs => macroxsiso_calculate_xs
  end type MacroXSIso

  type, extends(MacroXS) :: MacroXSAngle
    ! Macroscopic cross sections
    real(8), allocatable :: total(:,:,:)      ! total cross section
    real(8), allocatable :: absorption(:,:,:) ! absorption cross section
    type(ScattDataContainer), allocatable :: scatter(:,:) ! scattering information
    real(8), allocatable :: nu_fission(:,:,:) ! nu-fission
    real(8), allocatable :: k_fission(:,:,:)  ! kappa-fission
    real(8), allocatable :: fission(:,:,:)    ! fission x/s
    real(8), allocatable :: chi(:,:,:,:)      ! fission spectra
    real(8), allocatable :: scattxs(:,:,:)    ! scattering xs
    real(8), allocatable :: polar(:)          ! polar angles
    real(8), allocatable :: azimuthal(:)      ! azimuthal angles

  contains
    procedure :: init     => macroxsangle_init   ! inits object
    procedure :: get_xs   => macroxsangle_get_xs ! Returns xs
    procedure :: sample_fission_energy => macroxsangle_sample_fission_energy
    procedure :: sample_scatter => macroxsangle_sample_scatter
    procedure :: calculate_xs => macroxsangle_calculate_xs
  end type MacroXSAngle

!===============================================================================
! MACROXSCONTAINER pointer array for storing MacroXS objects.
!===============================================================================

  type MacroXSContainer
    class(MacroXS), allocatable :: obj
  end type MacroXSContainer

contains

!===============================================================================
! MACROXS*_INIT sets the MacroXS Data
!===============================================================================

  subroutine macroxsiso_init(this, mat, nuclides, groups, get_kfiss, get_fiss, &
       max_order, scatt_type, legendre_mu_points, error_code, error_text)
    class(MacroXSIso), intent(inout)     :: this ! The MacroXS to initialize
    type(Material), pointer, intent(in)  :: mat  ! base material
    type(NuclideMGContainer), intent(in) :: nuclides(:) ! List of nuclides to harvest from
    integer, intent(in)                  :: groups ! Number of E groups
    logical, intent(in)                  :: get_kfiss ! Should we get kfiss data?
    logical, intent(in)                  :: get_fiss ! Should we get fiss data?
    integer, intent(in)                  :: max_order ! Maximum requested order
    integer, intent(in)                  :: scatt_type ! How is data presented
    integer, intent(in)                  :: legendre_mu_points ! Treat as Leg or Tabular?
    integer, intent(inout)               :: error_code  ! Code signifying error
    character(MAX_LINE_LEN), intent(inout) :: error_text ! Error message to print

    integer :: i             ! loop index over nuclides
    integer :: gin, gout     ! group indices
    real(8) :: atom_density  ! atom density of a nuclide
    integer :: imu
    real(8) :: norm
    integer :: mat_max_order, order, l
    real(8), allocatable :: temp_mult(:,:)
    real(8), allocatable :: temp_energy(:,:)
    real(8), allocatable :: scatt_coeffs(:,:,:)

    ! Initialize error data
    error_code = 0
    error_text = ''

    ! If we have tabular only data, then make sure all datasets have same size
    if (scatt_type == ANGLE_HISTOGRAM) then
      ! Check all scattering data of same size
      order = nuclides(mat % nuclide(1)) % obj % order
      do i = 2, mat % n_nuclides
        if (order /= nuclides(mat % nuclide(i)) % obj % order) then
          error_code = 1
          error_text = "All Histogram Scattering Entries Must Be Same Length!"
          return
        end if
      end do
      ! Ok, got our order, store it
      this % order = order

      ! Allocate stuff for later
      allocate(scatt_coeffs(order, groups, groups))
      scatt_coeffs = ZERO
      allocate(ScattDataHistogram :: this % scatter)

    else if (scatt_type == ANGLE_TABULAR) then
      ! Check all scattering data of same size
      order = nuclides(mat % nuclide(1)) % obj % order
      do i = 2, mat % n_nuclides
        if (order /= nuclides(mat % nuclide(i)) % obj % order) then
          error_code = 1
          error_text = "All Tabular Scattering Entries Must Be Same Length!"
          return
        end if
      end do
      ! Ok, got our order, store it
      this % order = order

      ! Allocate stuff for later
      allocate(scatt_coeffs(order, groups, groups))
      scatt_coeffs = ZERO
      allocate(ScattDataTabular :: this % scatter)

    else if (scatt_type == ANGLE_LEGENDRE) then
      ! Otherwise find the maximum scattering order
      ! Need to determine the maximum scattering order of all data in this material
      mat_max_order = 0
      do i = 1, mat % n_nuclides
        if (nuclides(mat % nuclide(i)) % obj % order > mat_max_order) then
          mat_max_order = nuclides(mat % nuclide(i)) % obj % order
        end if
      end do

      ! Now need to compare this material maximum scattering order with
      ! the problem wide max scatt order and use whichever is lower
      order = min(mat_max_order, max_order)
      this % order = order + 1

      ! Now we can allocate our scatt_coeffs object accordingly
      allocate(scatt_coeffs(order + 1, groups, groups))
      scatt_coeffs = ZERO
      if (legendre_mu_points == 1) then
        allocate(ScattDataLegendre :: this % scatter)
      else
        allocate(ScattDataTabular :: this % scatter)
      end if
    end if

    ! Allocate and initialize data within macro_xs(i_mat) object
    allocate(this % total(groups))
    this % total = ZERO
    allocate(this % absorption(groups))
    this % absorption = ZERO
    if (get_fiss) then
      allocate(this % fission(groups))
      this % fission = ZERO
    end if
    if (get_kfiss) then
      allocate(this % k_fission(groups))
      this % k_fission = ZERO
    end if
    allocate(this % nu_fission(groups))
    this % nu_fission = ZERO
    allocate(this % chi(groups, groups))
    this % chi = ZERO
    allocate(temp_energy(groups, groups))
    temp_energy = ZERO
    allocate(temp_mult(groups, groups))
    temp_mult = ZERO
    allocate(this % scattxs(groups))

    ! Add contribution from each nuclide in material
    do i = 1, mat % n_nuclides
      ! Copy atom density of nuclide in material
      atom_density = mat % atom_density(i)

      ! Perform our operations which depend upon the type
      select type(nuc => nuclides(mat % nuclide(i)) % obj)
      type is (NuclideIso)

        ! Add contributions to total, absorption, and fission data (if necessary)
        this % total = this % total + atom_density * nuc % total
        this % absorption = this % absorption + &
             atom_density * nuc % absorption
        if (nuc % fissionable) then
          if (allocated(nuc % chi)) then
            do gin = 1, groups
              do gout = 1, groups
                this % chi(gout,gin) = this % chi(gout,gin) + atom_density * &
                     nuc % chi(gout) * nuc % nu_fission(gin,1)
              end do
            end do
            this % nu_fission = this % nu_fission + atom_density * &
                 nuc % nu_fission(:,1)
          else
            this % chi = this % chi + atom_density * nuc % nu_fission
            do gin = 1, groups
              this % nu_fission(gin) = this % nu_fission(gin) + atom_density * &
                   sum(nuc % nu_fission(:,gin))
            end do
          end if
          if (get_fiss) then
            this % fission = this % fission + atom_density * nuc % fission
          end if
          if (get_kfiss) then
            this % k_fission = this % k_fission + atom_density * nuc % k_fission
          end if
        end if

        ! Now time to do the scattering
        do gin = 1, groups
          do gout = 1, groups
            if (scatt_type == ANGLE_HISTOGRAM .or. scatt_type == ANGLE_TABULAR) then
              ! Transfer matrix
              temp_energy(gout,gin) = temp_energy(gout,gin) + atom_density * &
                   sum(nuc % scatter(gout,gin,:))

              ! Determine the angular distribution
              do imu = 1, order
                scatt_coeffs(imu, gout, gin) = scatt_coeffs(imu, gout, gin) + &
                     nuc % scatter(gout,gin,imu) * &
                     atom_density
              end do

            else if (scatt_type == ANGLE_LEGENDRE) then
              ! Transfer matrix
              temp_energy(gout,gin) = temp_energy(gout,gin) + atom_density * &
                   nuc % scatter(gout,gin,1)

              ! Determine the angular distribution coefficients so we can later
              ! expand do the complete distribution
              do l = 1, min(nuc % order, order) + 1
                scatt_coeffs(l, gout, gin) = scatt_coeffs(l, gout, gin) + &
                     nuc % scatter(gout,gin,l) * &
                     atom_density
              end do

            end if

            ! Multiplicity matrix
            temp_mult(gout,gin) = temp_mult(gout,gin) + atom_density * &
                 nuc % mult(gout,gin)
          end do
        end do
      type is (NuclideAngle)
        error_code = 1
        error_text = "Invalid Passing of NuclideAngle to MacroXSIso Object"
        return
      end select
    end do

    ! Store the scattering xs
    if (scatt_type == ANGLE_HISTOGRAM .or. scatt_type == ANGLE_TABULAR) then
      this % scattxs(:) = sum(sum(scatt_coeffs(:,:,:),dim=1),dim=1)
    else if (scatt_type == ANGLE_LEGENDRE) then
      this % scattxs(:) = sum(scatt_coeffs(1,:,:),dim=1)
    end if

    ! Normalize the scatt_coeffs
    do gin = 1, groups
      do gout = 1, groups
        if (scatt_type == ANGLE_HISTOGRAM .or. scatt_type == ANGLE_TABULAR) then
          norm = sum(scatt_coeffs(:,gout,gin))
        else if (scatt_type == ANGLE_LEGENDRE) then
          norm = scatt_coeffs(1,gout,gin)
        end if
        if (norm /= ZERO) then
          scatt_coeffs(:, gout, gin) = scatt_coeffs(:, gout,gin) / norm
        end if
      end do
      ! Now normalize temp_energy (outgoing scattering energy probabilities)
      norm = sum(temp_energy(:,gin))
      if (norm > ZERO) then
        temp_energy(:,gin) = temp_energy(:,gin) / norm
      end if
    end do

    if (scatt_type == ANGLE_LEGENDRE .and. legendre_mu_points /= 1) then
      call this % scatter % init(legendre_mu_points, temp_energy, temp_mult, &
                               scatt_coeffs)
    else
      call this % scatter % init(this % order, temp_energy, temp_mult, &
                               scatt_coeffs)
    end if

    ! Now normalize chi
    if (mat % fissionable) then
      do gin = 1, groups
        ! Normalize Chi
        norm =  sum(this % chi(:,gin))
        if (norm > ZERO) then
          this % chi(:,gin) = this % chi(:,gin) / norm
        end if
      end do
    end if

    ! Deallocate temporaries for the next material
    deallocate(scatt_coeffs, temp_energy, temp_mult)

  end subroutine macroxsiso_init

  subroutine macroxsangle_init(this, mat, nuclides, groups, get_kfiss, get_fiss, &
       max_order, scatt_type, legendre_mu_points, error_code, error_text)
    class(MacroXSAngle), intent(inout)   :: this ! The MacroXS to initialize
    type(Material), pointer, intent(in)  :: mat  ! base material
    type(NuclideMGContainer), intent(in) :: nuclides(:) ! List of nuclides to harvest from
    integer, intent(in)                  :: groups ! Number of E groups
    logical, intent(in)                  :: get_kfiss ! Should we get kfiss data?
    logical, intent(in)                  :: get_fiss ! Should we get fiss data?
    integer, intent(in)                  :: max_order ! Maximum requested order
    integer, intent(in)                  :: scatt_type ! Legendre or Tabular Scatt?
    integer, intent(in)                  :: legendre_mu_points ! Treat as Leg or Tabular?
    integer, intent(inout)               :: error_code  ! Code signifying error
    character(MAX_LINE_LEN), intent(inout) :: error_text ! Error message to print

    integer :: i             ! loop index over nuclides
    integer :: gin, gout     ! group indices
    real(8) :: atom_density  ! atom density of a nuclide
    integer :: ipol, iazi, npol, nazi
    integer :: imu
    real(8) :: norm
    integer :: mat_max_order, order, l
    real(8), allocatable :: temp_mult(:,:,:,:)
    real(8), allocatable :: temp_energy(:,:,:,:)
    real(8), allocatable :: scatt_coeffs(:,:,:,:,:)

    ! Initialize error data
    error_code = 0
    error_text = ''

    ! Get the number of each polar and azi angles and make sure all the
    ! NuclideAngle types have the same number of these angles
    npol = -1
    nazi = -1
    do i = 1, mat % n_nuclides
      select type(nuc => nuclides(mat % nuclide(i)) % obj)
      type is (NuclideAngle)
        if (npol == -1) then
          npol = nuc % n_pol
          nazi = nuc % n_azi
          allocate(this % polar(npol))
          this % polar = nuc % polar
          allocate(this % azimuthal(nazi))
          this % azimuthal = nuc % azimuthal
        else
          if ((npol /= nuc % n_pol) .or. (nazi /= nuc % n_azi)) then
            error_code = 1
            error_text = "All Angular Data Must Be Same Length!"
          end if
        end if
      end select
    end do

    ! If we have tabular only data, then make sure all datasets have same size
    if (scatt_type == ANGLE_HISTOGRAM) then
      ! Check all scattering data of same size
      order = nuclides(mat % nuclide(1)) % obj % order
      do i = 2, mat % n_nuclides
        if (order /= nuclides(mat % nuclide(i)) % obj % order) then
          error_code = 1
          error_text = "All Histogram Scattering Entries Must Be Same Length!"
          return
        end if
      end do
      ! Ok, got our order, store it
      this % order = order

      ! Allocate stuff for later
      allocate(scatt_coeffs(order, groups, groups, nazi, npol))
      scatt_coeffs = ZERO
      allocate(this % scatter(nazi, npol))
      do ipol = 1, npol
        do iazi = 1, nazi
          allocate(ScattDataHistogram :: this % scatter(iazi, ipol) % obj)
        end do
      end do

    else if (scatt_type == ANGLE_TABULAR) then
      ! Check all scattering data of same size
      order = nuclides(mat % nuclide(1)) % obj % order
      do i = 2, mat % n_nuclides
        if (order /= nuclides(mat % nuclide(i)) % obj % order) then
          error_code = 1
          error_text = "All Tabular Scattering Entries Must Be Same Length!"
          return
        end if
      end do
      ! Ok, got our order, store it
      this % order = order

      ! Allocate stuff for later
      allocate(scatt_coeffs(order, groups, groups, nazi, npol))
      scatt_coeffs = ZERO
      allocate(this % scatter(nazi, npol))
      do ipol = 1, npol
        do iazi = 1, nazi
          allocate(ScattDataTabular :: this % scatter(iazi, ipol) % obj)
        end do
      end do

    else if (scatt_type == ANGLE_LEGENDRE) then
      ! Otherwise find the maximum scattering order
      ! Need to determine the maximum scattering order of all data in this material
      mat_max_order = 0
      do i = 1, mat % n_nuclides
        if (nuclides(mat % nuclide(i)) % obj % order > mat_max_order) then
          mat_max_order = nuclides(mat % nuclide(i)) % obj % order
        end if
      end do

      ! Now need to compare this material maximum scattering order with
      ! the problem wide max scatt order and use whichever is lower
      order = min(mat_max_order, max_order)
      this % order = order + 1

      ! Now we can allocate our scatt_coeffs object accordingly
      allocate(scatt_coeffs(order + 1, groups, groups, nazi, npol))
      scatt_coeffs = ZERO
      allocate(this % scatter(nazi, npol))
      do ipol = 1, npol
        do iazi = 1, nazi
          if (legendre_mu_points == 1) then
            allocate(ScattDataLegendre :: this % scatter(iazi, ipol) % obj)
          else
            allocate(ScattDataTabular :: this % scatter(iazi, ipol) % obj)
          end if
        end do
      end do
    end if

    ! Allocate and initialize data within macro_xs(i_mat) object
    allocate(this % total(groups,nazi,npol))
    this % total = ZERO
    allocate(this % absorption(groups,nazi,npol))
    this % absorption = ZERO
    if (get_fiss) then
      allocate(this % fission(groups,nazi,npol))
      this % fission = ZERO
    end if
    if (get_kfiss) then
      allocate(this % k_fission(groups,nazi,npol))
      this % k_fission = ZERO
    end if
    allocate(this % nu_fission(groups,nazi,npol))
    this % nu_fission = ZERO
    allocate(this % chi(groups, groups, nazi, npol))
    this % chi = ZERO
    allocate(temp_energy(groups,groups,nazi,npol))
    temp_energy = ZERO
    allocate(temp_mult(groups,groups,nazi,npol))
    temp_mult = ZERO
    allocate(this % scattxs(groups,nazi,npol))

    ! Add contribution from each nuclide in material
    do i = 1, mat % n_nuclides
      ! Copy atom density of nuclide in material
      atom_density = mat % atom_density(i)

      ! Perform our operations which depend upon the type
      select type(nuc => nuclides(mat % nuclide(i)) % obj)
      type is (NuclideIso)
        error_code = 1
        error_text = "Invalid Passing of NuclideIso to MacroXSAngle Object"
        return
      type is (NuclideAngle)
        ! Add contributions to total, absorption, and fission data (if necessary)
        this % total = this % total + atom_density * nuc % total
        this % absorption = this % absorption + &
             atom_density * nuc % absorption
        if (nuc % fissionable) then
          if (allocated(nuc % chi)) then
            do gin = 1, groups
              do gout = 1, groups
                this % chi(gout,gin,:,:) = this % chi(gout,gin,:,:) + atom_density * &
                     nuc % chi(gout,:,:) * nuc % nu_fission(gin,1,:,:)
              end do
            end do
            this % nu_fission = this % nu_fission + atom_density * &
                 nuc % nu_fission(:,1,:,:)
          else
            this % chi = this % chi + atom_density * nuc % nu_fission
            do gin = 1, groups
              this % nu_fission(gin,:,:) = this % nu_fission(gin,:,:) + atom_density * &
                   sum(nuc % nu_fission(:,gin,:,:),dim=1)
            end do
          end if
          if (get_fiss) then
            this % fission = this % fission + atom_density * nuc % fission
          end if
          if (get_kfiss) then
            this % k_fission = this % k_fission + atom_density * nuc % k_fission
          end if
        end if

        ! Now time to do the scattering
        do gin = 1, groups
          do gout = 1, groups
            if (scatt_type == ANGLE_HISTOGRAM .or. scatt_type == ANGLE_TABULAR) then
              ! Transfer matrix
              temp_energy(gout,gin,:,:) = temp_energy(gout,gin,:,:) + atom_density * &
                   sum(nuc % scatter(gout,gin,:,:,:),dim=1)

              ! Determine the angular distribution
              do imu = 1, order
                scatt_coeffs(imu,gout,gin,:,:) = scatt_coeffs(imu,gout,gin,:,:) + &
                     nuc % scatter(gout,gin,imu,:,:) * &
                     atom_density
              end do
            else if (scatt_type == ANGLE_LEGENDRE) then
              ! Transfer matrix
              temp_energy(gout,gin,:,:) = temp_energy(gout,gin,:,:) + atom_density * &
                   nuc % scatter(gout,gin,1,:,:)

              ! Determine the angular distribution coefficients so we can later
              ! expand do the complete distribution
              do l = 1, min(nuc % order, order) + 1
                scatt_coeffs(l, gout, gin,:,:) = scatt_coeffs(l, gout, gin,:,:) + &
                     nuc % scatter(gout,gin,l,:,:) * &
                     atom_density
              end do
            end if

            ! Multiplicity matrix
            temp_mult(gout,gin,:,:) = temp_mult(gout,gin,:,:) + atom_density * &
                 nuc % mult(gout,gin,:,:)
          end do
        end do
      end select
    end do

    ! Store the scattering xs
    if (scatt_type == ANGLE_HISTOGRAM .or. scatt_type == ANGLE_TABULAR) then
      this % scattxs(:,:,:) = sum(sum(scatt_coeffs(:,:,:,:,:),dim=1),dim=1)
    else if (scatt_type == ANGLE_LEGENDRE) then
      this % scattxs(:,:,:) = sum(scatt_coeffs(1,:,:,:,:),dim=1)
    end if

    ! Normalize the scatt_coeffs
    do ipol = 1, npol
      do iazi = 1, nazi
        do gin = 1, groups
          do gout = 1, groups
            if (scatt_type == ANGLE_HISTOGRAM .or. scatt_type == ANGLE_TABULAR) then
              norm = sum(scatt_coeffs(:,gout,gin,iazi,ipol))
            else if (scatt_type == ANGLE_LEGENDRE) then
              norm = scatt_coeffs(1,gout,gin,iazi,ipol)
            end if
            if (norm /= ZERO) then
              scatt_coeffs(:,gout,gin,iazi,ipol) = &
                   scatt_coeffs(:,gout,gin,iazi,ipol) / norm
            end if
          end do
          ! Now normalize temp_energy (outgoing scattering energy probabilities)
          norm = sum(temp_energy(:,gin,iazi,ipol))
          if (norm > ZERO) then
            temp_energy(:,gin,iazi,ipol) = temp_energy(:,gin,iazi,ipol) / norm
          end if
        end do

        if (scatt_type == ANGLE_LEGENDRE .and. legendre_mu_points /= 1) then
          call this % scatter(iazi, ipol) % obj % init(legendre_mu_points, &
                 temp_energy(:,:,iazi,ipol), temp_mult(:,:,iazi,ipol), &
                 scatt_coeffs(:,:,:,iazi,ipol))
        else
          call this % scatter(iazi, ipol) % obj % init(this % order, &
                 temp_energy(:,:,iazi,ipol), temp_mult(:,:,iazi,ipol), &
                 scatt_coeffs(:,:,:,iazi,ipol))
        end if

      end do
    end do

    ! Now go through and normalize chi
    if (mat % fissionable) then
      do ipol = 1, npol
        do iazi = 1, nazi
          do gin = 1, groups
            ! Normalize Chi
            norm =  sum(this % chi(:,gin,iazi,ipol))
            if (norm > ZERO) then
              this % chi(:,gin,iazi,ipol) = this % chi(:,gin,iazi,ipol) / norm
            end if
          end do
        end do
      end do
    end if

    ! Deallocate temporaries for the next material
    deallocate(scatt_coeffs, temp_energy, temp_mult)

  end subroutine macroxsangle_init

!===============================================================================
! MACROXS_*_GET_XS returns the requested data type
!===============================================================================

  function macroxsiso_get_xs(this, g, xstype, gout, uvw) result(xs)
    class(MacroXSIso), intent(in)  :: this   ! The MacroXS to initialize
    integer, intent(in)            :: g      ! Incoming Energy group
    character(*) , intent(in)      :: xstype ! Type of xs requested
    integer, optional, intent(in)  :: gout   ! Outgoing Energy group
    real(8), optional, intent(in)  :: uvw(3) ! Requested Angle
    real(8)                        :: xs     ! Requested x/s

    select case(xstype)
    case('total')
      xs = this % total(g)
    case('absorption')
      xs = this % absorption(g)
    case('fission')
      xs = this % fission(g)
    case('k_fission')
      xs = this % k_fission(g)
    case('nu_fission')
      xs = this % nu_fission(g)
    case('scatter')
      xs = this % scattxs(g)
    case('mult')
      if (present(gout)) then
        xs = this % scatter % mult(gout,g)
      else
        xs = sum(this % scatter % mult(:,g))
      end if
    end select

  end function macroxsiso_get_xs

  function macroxsangle_get_xs(this, g, xstype, gout,uvw) result(xs)
    class(MacroXSAngle), intent(in)  :: this   ! The MacroXS to initialize
    integer, intent(in)              :: g      ! Incoming Energy group
    character(*) , intent(in)        :: xstype ! Type of xs requested
    integer, optional, intent(in)    :: gout   ! Outgoing Energy group
    real(8), optional, intent(in)    :: uvw(3) ! Requested Angle
    real(8)                          :: xs     ! Requested x/s

    integer :: iazi, ipol

    if (present(uvw)) then
      call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)
      select case(xstype)
      case('total')
        xs = this % total(g,iazi,ipol)
      case('absorption')
        xs = this % absorption(g,iazi,ipol)
      case('fission')
        xs = this % fission(g,iazi,ipol)
      case('k_fission')
        xs = this % k_fission(g,iazi,ipol)
      case('nu_fission')
        xs = this % nu_fission(g,iazi,ipol)
      case('scatter')
        xs = this % scattxs(g,iazi,ipol)
      case('mult')
        if (present(gout)) then
          xs = this % scatter(iazi,ipol) % obj % mult(gout,g)
        else
          xs = sum(this % scatter(iazi,ipol) % obj % mult(:,g))
        end if
      end select
    end if

  end function macroxsangle_get_xs

!===============================================================================
! MACROXS_*_SAMPLE_FISSION_ENERGY samples the outgoing energy from a fission
! event
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
! MACROXS*_SAMPLE_SCATTER Selects outgoing energy and angle after a scatter
! event
!===============================================================================

  subroutine macroxsiso_sample_scatter(this, uvw, gin, gout, mu, wgt)
    class(MacroXSIso), intent(in)    :: this
    real(8),           intent(in)    :: uvw(3) ! Incoming neutron direction
    integer,           intent(in)    :: gin    ! Incoming neutron group
    integer,           intent(out)   :: gout   ! Sampled outgoin group
    real(8),           intent(out)   :: mu     ! Sampled change in angle
    real(8),           intent(inout) :: wgt    ! Particle weight

    call this % scatter % sample(gin, gout, mu, wgt)

  end subroutine macroxsiso_sample_scatter

  subroutine macroxsangle_sample_scatter(this, uvw, gin, gout, mu, wgt)
    class(MacroXSAngle), intent(in)    :: this
    real(8),             intent(in)    :: uvw(3) ! Incoming neutron direction
    integer,             intent(in)    :: gin    ! Incoming neutron group
    integer,             intent(out)   :: gout   ! Sampled outgoin group
    real(8),             intent(out)   :: mu     ! Sampled change in angle
    real(8),             intent(inout) :: wgt    ! Particle weight

    integer :: iazi, ipol ! Angular indices

    call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)
    call this % scatter(iazi,ipol) % obj % sample(gin,gout,mu,wgt)

  end subroutine macroxsangle_sample_scatter

!===============================================================================
! MACROXS*_CALCULATE_XS determines the multi-group macroscopic cross sections
! for the material the particle is currently traveling through.
!===============================================================================

  subroutine macroxsiso_calculate_xs(this, gin, uvw, xs)
    class(MacroXSIso),     intent(in)    :: this
    integer,               intent(in)    :: gin    ! Incoming neutron group
    real(8),               intent(in)    :: uvw(3) ! Incoming neutron direction
    type(MaterialMacroXS), intent(inout) :: xs     ! Resultant MacroXS Data

    xs % total         = this % total(gin)
    xs % elastic       = this % scattxs(gin)
    xs % absorption    = this % absorption(gin)
    xs % nu_fission    = this % nu_fission(gin)

  end subroutine macroxsiso_calculate_xs

  subroutine macroxsangle_calculate_xs(this, gin, uvw, xs)
    class(MacroXSAngle),   intent(in)    :: this
    integer,               intent(in)    :: gin    ! Incoming neutron group
    real(8),               intent(in)    :: uvw(3) ! Incoming neutron direction
    type(MaterialMacroXS), intent(inout) :: xs     ! Resultant MacroXS Data

    integer :: iazi, ipol

    call find_angle(this % polar, this % azimuthal, uvw, iazi, ipol)
    xs % total         = this % total(gin, iazi, ipol)
    xs % elastic       = this % scattxs(gin, iazi, ipol)
    xs % absorption    = this % absorption(gin, iazi, ipol)
    xs % nu_fission    = this % nu_fission(gin, iazi, ipol)

  end subroutine macroxsangle_calculate_xs

end module macroxs_header
