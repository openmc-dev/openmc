module macroxs_header

  use constants,       only: MAX_FILE_LEN, ZERO, ONE, TWO, PI
  use list_header,     only: ListInt
  use material_header, only: material
  use math,            only: calc_pn, calc_rn, expand_harmonic
  use nuclide_header
  use scattdata_header

  implicit none

!===============================================================================
! MACROXS_* contains cached macroscopic cross sections for the material a
! particle is traveling through
!===============================================================================

  type, abstract :: MacroXS_Base
    ! Data Order
    integer :: order

    ! Type-Bound procedures
    contains
      procedure(macroxs_init_),   deferred, pass :: init     ! initializes object
      procedure(macroxs_clear_),  deferred, pass :: clear    ! Deallocates object
      procedure(macroxs_get_xs_), deferred, pass :: get_xs   ! Return xs
  end type MacroXS_Base

  abstract interface
    subroutine macroxs_init_(this, mat, nuclides, groups, get_kfiss, get_fiss, &
                             max_order, scatt_type, legendre_mu_points, &
                             error_code, error_text)

      import MacroXS_Base
      import Material
      import NuclideMGContainer
      import MAX_LINE_LEN
      class(MacroXS_Base), intent(inout)   :: this ! The MacroXS to initialize
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
      import MacroXS_Base
      class(MacroXS_Base), intent(in) :: this   ! The MacroXS to initialize
      integer, intent(in)             :: g      ! Incoming Energy group
      character(*) , intent(in)       :: xstype ! Cross Section Type
      integer, optional, intent(in)   :: gout   ! Outgoing Energy group
      real(8), optional, intent(in)   :: uvw(3) ! Requested Angle
      real(8)                         :: xs     ! Resultant xs

    end function macroxs_get_xs_

    subroutine macroxs_clear_(this)

      import MacroXS_Base
      class(MacroXS_Base), intent(inout) :: this ! The MacroXS to clear

    end subroutine macroxs_clear_

  end interface

  type, extends(MacroXS_Base) :: MacroXS_Iso
    ! Microscopic cross sections
    real(8), allocatable :: total(:)      ! total cross section
    real(8), allocatable :: absorption(:) ! absorption cross section
    class(ScattData_Base), allocatable :: scatter ! scattering information
    real(8), allocatable :: nu_fission(:) ! nu-fission
    real(8), allocatable :: k_fission(:)  ! kappa-fission
    real(8), allocatable :: fission(:)    ! fission x/s
    real(8), allocatable :: scattxs(:)    ! scattering xs
    real(8), allocatable :: chi(:,:)      ! fission spectra

    ! Type-Bound procedures
    contains
      procedure, pass :: init        => macroxs_iso_init        ! inits object
      procedure, pass :: clear       => macroxs_iso_clear       ! Deallocates object
      procedure, pass :: get_xs      => macroxs_iso_get_xs      ! Returns xs
  end type MacroXS_Iso

  type, extends(MacroXS_Base) :: MacroXS_Angle
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

    ! Type-Bound procedures
    contains
      procedure, pass :: init     => macroxs_angle_init   ! inits object
      procedure, pass :: clear    => macroxs_angle_clear  ! Deallocates object
      procedure, pass :: get_xs   => macroxs_angle_get_xs ! Returns xs
  end type MacroXS_Angle

!===============================================================================
! MACROXSCONTAINER pointer array for storing MacroXS objects.
!===============================================================================

  type MacroXSContainer
    class(MacroXS_Base), allocatable :: obj
  end type MacroXSContainer

contains

!===============================================================================
! MACROXS*_INIT sets the MacroXS Data
!===============================================================================

  subroutine macroxs_iso_init(this, mat, nuclides, groups, get_kfiss, get_fiss, &
       max_order, scatt_type, legendre_mu_points, error_code, error_text)

    class(MacroXS_Iso), intent(inout)    :: this ! The MacroXS to initialize
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
    ! class(Nuclide_Base), pointer  :: nuc ! current nuclide
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
      allocate(ScattData_Histogram :: this % scatter)

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
      allocate(ScattData_Tabular :: this % scatter)

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
        allocate(ScattData_Legendre :: this % scatter)
      else
        allocate(ScattData_Tabular :: this % scatter)
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
      type is (Nuclide_Iso)

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
      type is (Nuclide_Angle)
        error_code = 1
        error_text = "Invalid Passing of Nuclide_Angle to MacroXS_Iso Object"
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

  end subroutine macroxs_iso_init

  subroutine macroxs_angle_init(this, mat, nuclides, groups, get_kfiss, get_fiss, &
       max_order, scatt_type, legendre_mu_points, error_code, error_text)

    class(MacroXS_Angle), intent(inout)  :: this ! The MacroXS to initialize
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
    ! Nuclide_Angle types have the same number of these angles
    npol = -1
    nazi = -1
    do i = 1, mat % n_nuclides
      select type(nuc => nuclides(mat % nuclide(i)) % obj)
      type is (Nuclide_Angle)
        if (npol == -1) then
          npol = nuc % Npol
          nazi = nuc % Nazi
          allocate(this % polar(npol))
          this % polar = nuc % polar
          allocate(this % azimuthal(nazi))
          this % azimuthal = nuc % azimuthal
        else
          if ((npol /= nuc % Npol) .or. (nazi /= nuc % Nazi)) then
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
          allocate(ScattData_Histogram :: this % scatter(iazi, ipol) % obj)
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
          allocate(ScattData_Tabular :: this % scatter(iazi, ipol) % obj)
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
            allocate(ScattData_Legendre :: this % scatter(iazi, ipol) % obj)
          else
            allocate(ScattData_Tabular :: this % scatter(iazi, ipol) % obj)
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
      type is (Nuclide_Iso)
        error_code = 1
        error_text = "Invalid Passing of Nuclide_Iso to MacroXS_Angle Object"
        return
      type is (Nuclide_Angle)
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

  end subroutine macroxs_angle_init

!===============================================================================
! MACROXS*_CLEAR resets and deallocates data in MacroXS.
!===============================================================================

  subroutine macroxs_iso_clear(this)

    class(MacroXS_Iso), intent(inout) :: this ! The MacroXS to clear

    if (allocated(this % total)) then
      deallocate(this % total, this % absorption, &
                 this % nu_fission)
    end if

    if (allocated(this % fission)) then
      deallocate(this % fission)
    end if

    if (allocated(this % k_fission)) then
      deallocate(this % k_fission)
    end if

    call this % scatter % clear()

    if (allocated(this % chi)) then
      deallocate(this % chi)
    end if

  end subroutine macroxs_iso_clear

  subroutine macroxs_angle_clear(this)

    class(MacroXS_Angle), intent(inout) :: this ! The MacroXS to clear
    integer :: i, j

    if (allocated(this % total)) then
      deallocate(this % total, this % absorption, &
                 this % nu_fission)
    end if

    if (allocated(this % fission)) then
      deallocate(this % fission)
    end if

    if (allocated(this % k_fission)) then
      deallocate(this % k_fission)
    end if

    do i = 1, size(this % scatter,dim=2)
      do j = 1, size(this % scatter,dim=1)
        call this % scatter(j,i) % obj % clear()
      end do
    end do
    if (allocated(this % scatter)) then
      deallocate(this % scatter)
    end if

    if (allocated(this % chi)) then
      deallocate(this % chi)
    end if

  end subroutine macroxs_angle_clear

!===============================================================================
! MACROXS_*_GET_XS returns the requested data type
!===============================================================================

  function macroxs_iso_get_xs(this, g, xstype, gout, uvw) result(xs)
    class(MacroXS_Iso), intent(in) :: this   ! The MacroXS to initialize
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

  end function macroxs_iso_get_xs

  function macroxs_angle_get_xs(this, g, xstype, gout,uvw) result(xs)
    class(MacroXS_Angle), intent(in) :: this   ! The MacroXS to initialize
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

  end function macroxs_angle_get_xs

!===============================================================================
! THIN_GRID thins an (x,y) set while also thinning an associated y2
!===============================================================================

  subroutine thin_grid(xout, yout, yout2, tol, compression, maxerr)
    real(8), allocatable, intent(inout) :: xout(:)     ! Resultant x grid
    real(8), allocatable, intent(inout) :: yout(:)     ! Resultant y values
    real(8), allocatable, intent(inout) :: yout2(:)     ! Secondary y values
    real(8), intent(in)                 :: tol         ! Desired fractional error to maintain
    real(8), intent(out)                :: compression ! Data reduction fraction
    real(8), intent(inout)              :: maxerr      ! Maximum error due to compression

    real(8), allocatable :: xin(:)      ! Incoming x grid
    real(8), allocatable :: yin(:)      ! Incoming y values
    real(8), allocatable :: yin2(:)     ! Secondary Incoming y values
    integer :: k, klo, khi
    integer :: all_ok
    real(8) :: x1, y1, x2, y2, x, y, testval
    integer :: num_keep, remove_it
    real(8) :: initial_size
    real(8) :: error
    real(8) :: x_frac

    initial_size = real(size(xout), 8)

    allocate(xin(size(xout)))
    xin = xout
    allocate(yin(size(yout)))
    yin = yout
    allocate(yin2(size(yout2)))
    yin2 = yout2

    all_ok = size(yin)
    maxerr = 0.0_8

    ! This loop will step through each entry in dim==3 and check to see if
    ! all of the values in other 2 dims can be replaced with linear interp.
    ! If not, the value will be saved to a new array, if so, it will be
    ! skipped.

    xout = 0.0_8
    yout = 0.0_8

    ! Keep first point's data
    xout(1) = xin(1)
    yout(1) = yin(1)
    yout2(1) = yin2(1)

    ! Initialize data
    num_keep = 1
    klo = 1
    khi = 3
    k = 2
    do while (khi <= size(xin))
      remove_it = 0
      x1 = xin(klo)
      x2 = xin(khi)
      x  = xin(k)
      x_frac = 1.0_8 / (x2 - x1) * (x - x1)  ! Linear interp.

      ! Check for removal. Otherwise, it stays. This is accomplished by leaving
      ! remove_it as 0, entering the else portion of if(remove_it==all_ok)
      y1 = yin(klo)
      y2 = yin(khi)
      y  = yin(k)

      testval = y1 + (y2 - y1) * x_frac
      error = abs(testval - y)
      if (y /= 0.0_8) then
        error = error / y
      end if
      if (error <= tol) then
        remove_it = remove_it + 1
        if (error > maxerr) then
          maxerr = abs(testval - y)
        end if
      end if
      ! Now place the point in to the proper bin and advance iterators.
      if (remove_it /= 0) then
        ! Then don't put it in the new grid but advance iterators
        k = k + 1
        khi = khi + 1
      else
        ! Put it in new grid and advance iterators accordingly
        num_keep = num_keep + 1
        xout(num_keep) = xin(k)
        yout(num_keep) = yin(k)
        yout2(num_keep) = yin2(k)
        klo = k
        k = k + 1
        khi = khi + 1
      end if
    end do
    ! Save the last point's data
    num_keep = num_keep + 1
    xout(num_keep) = xin(size(xin))
    yout(num_keep) = yin(size(xin))
    yout2(num_keep) = yin2(size(xin))

    ! Finally, xout and yout were sized to match xin and yin since we knew
    ! they would be no larger than those.  Now we must resize these arrays
    ! and copy only the useful data in. Will use xin/yin for temp arrays.
    xin = xout(1:num_keep)
    yin = yout(1:num_keep)
    yin2 = yout2(1:num_keep)

    deallocate(xout)
    deallocate(yout)
    deallocate(yout2)
    allocate(xout(num_keep))
    allocate(yout(size(yin)))
    allocate(yout2(size(yin2)))

    xout = xin(1:num_keep)
    yout = yin(1:num_keep)
    yout2 = yin2(1:num_keep)

    ! Clean up
    deallocate(xin)
    deallocate(yin)
    deallocate(yin2)

    compression = (initial_size - real(size(xout),8)) / initial_size

  end subroutine thin_grid

end module macroxs_header
