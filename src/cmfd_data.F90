module cmfd_data

!==============================================================================
! CMFD_DATA -- This module processes the cmfd tally object to generate
! parameters for CMFD calculation.
!==============================================================================

  use cmfd_header,         only: allocate_cmfd, cmfd, cmfd_coremap, &
                                 cmfd_downscatter, cmfd_tallies, dhat_reset
  use constants
  use tally_filter_mesh, only: MeshFilter

  implicit none
  private
  public :: set_up_cmfd, neutron_balance

contains

!==============================================================================
! SET_UP_CMFD configures cmfd object for a CMFD eigenvalue calculation
!==============================================================================

  subroutine set_up_cmfd()

    use constants,           only: CMFD_NOACCEL

    ! Check for core map and set it up
    if ((cmfd_coremap) .and. (cmfd%mat_dim == CMFD_NOACCEL)) call set_coremap()

    ! Calculate all cross sections based on reaction rates from last batch
    call compute_xs()

    ! Compute effective downscatter cross section
    if (cmfd_downscatter) call compute_effective_downscatter()

    ! Check neutron balance
    call neutron_balance()

    ! Calculate dtilde
    call compute_dtilde()

    ! Calculate dhat
    call compute_dhat()

  end subroutine set_up_cmfd

!===============================================================================
! COMPUTE_XS takes tallies and computes macroscopic cross sections
!===============================================================================

  subroutine compute_xs()

    use constants,    only: FILTER_MESH, FILTER_ENERGYIN, FILTER_ENERGYOUT,    &
                           FILTER_SURFACE, OUT_LEFT, OUT_RIGHT, OUT_BACK,      &
                           OUT_FRONT, OUT_BOTTOM, OUT_TOP, IN_LEFT, IN_RIGHT,  &
                           IN_BACK, IN_FRONT, IN_BOTTOM, IN_TOP, CMFD_NOACCEL, &
                           ZERO, ONE, TINY_BIT
    use error,        only: fatal_error
    use mesh_header,  only: RegularMesh, meshes
    use string,       only: to_str
    use tally_filter_header, only: filters, filter_matches

    integer :: nx            ! number of mesh cells in x direction
    integer :: ny            ! number of mesh cells in y direction
    integer :: nz            ! number of mesh cells in z direction
    integer :: ng            ! number of energy groups
    integer :: i             ! iteration counter for x
    integer :: j             ! iteration counter for y
    integer :: k             ! iteration counter for z
    integer :: g             ! iteration counter for g
    integer :: h             ! iteration counter for outgoing groups
    integer :: l             ! iteration counter for tally filters
    integer :: ital          ! tally object index
    integer :: ijk(3)        ! indices for mesh cell
    integer :: score_index   ! index to pull from tally object
    integer :: i_filter_mesh ! index for mesh filter
    integer :: i_filter_ein  ! index for incoming energy filter
    integer :: i_filter_eout ! index for outgoing energy filter
    integer :: i_mesh        ! flattend index for mesh
    logical :: energy_filters! energy filters present
    real(8) :: flux          ! temp variable for flux
    type(RegularMesh), pointer :: m ! pointer for mesh object

    ! Extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! Set flux object and source distribution to all zeros
    cmfd % flux = ZERO
    cmfd % openmc_src = ZERO

    ! Associate tallies and mesh
    associate (t => cmfd_tallies(1) % obj)
      i_filter_mesh = t % filter(t % find_filter(FILTER_MESH))
    end associate

    select type(filt => filters(i_filter_mesh) % obj)
    type is (MeshFilter)
      m => meshes(filt % mesh)
    end select

    ! Set mesh widths
    cmfd % hxyz(1,:,:,:) = m % width(1) ! set x width
    cmfd % hxyz(2,:,:,:) = m % width(2) ! set y width
    cmfd % hxyz(3,:,:,:) = m % width(3) ! set z width

    cmfd % keff_bal = ZERO

    ! Begin loop around tallies
    TAL: do ital = 1, size(cmfd_tallies)

      ! Associate tallies and mesh
      associate (t => cmfd_tallies(ital) % obj)

      if (ital < 3) then
        i_filter_mesh = t % filter(t % find_filter(FILTER_MESH))
      else
        i_filter_mesh = t % filter(t % find_filter(FILTER_MESHSURFACE))
      end if

      ! Check for energy filters
      energy_filters = (t % find_filter(FILTER_ENERGYIN) > 0)

      if (energy_filters) then
        i_filter_ein  = t % filter(t % find_filter(FILTER_ENERGYIN))
        i_filter_eout = t % filter(t % find_filter(FILTER_ENERGYOUT))
      end if

      ! Begin loop around space
      ZLOOP: do k = 1,nz

        YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

            ! Check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
                cycle
              end if
            end if

            ! Loop around energy groups
            OUTGROUP: do h = 1,ng

              ! Start tally 1
              TALLY: if (ital == 1) then

                ! Reset all bins to 1
                do l = 1, size(t % filter)
                  call filter_matches(t % filter(l)) % bins % clear()
                  call filter_matches(t % filter(l)) % bins % push_back(1)
                end do

                ! Set ijk as mesh indices
                ijk = (/ i, j, k /)

                ! Get bin number for mesh indices
                filter_matches(i_filter_mesh) % bins % data(1) = &
                     m % get_bin_from_indices(ijk)

                ! Apply energy in filter
                if (energy_filters) then
                  filter_matches(i_filter_ein) % bins % data(1) = ng - h + 1
                end if

                ! Calculate score index from bins
                score_index = 1
                do l = 1, size(t % filter)
                  score_index = score_index + (filter_matches(t % filter(l)) &
                       % bins % data(1) - 1) * t % stride(l)
                end do

                ! Get flux
                flux = t % results(RESULT_SUM,1,score_index)
                cmfd % flux(h,i,j,k) = flux

                ! Detect zero flux, abort if located
                if ((flux - ZERO) < TINY_BIT) then
                  call fatal_error('Detected zero flux without coremap overlay &
                       &at: (' // to_str(i) // ',' // to_str(j) // ',' // &
                       &to_str(k) // ') in group ' // to_str(h))
                end if

                ! Get total rr and convert to total xs
                cmfd % totalxs(h,i,j,k) = t % results(RESULT_SUM,2,score_index) / flux

                ! Get p1 scatter rr and convert to p1 scatter xs
                cmfd % p1scattxs(h,i,j,k) = t % results(RESULT_SUM,3,score_index) / flux

                ! Calculate diffusion coefficient
                cmfd % diffcof(h,i,j,k) = ONE/(3.0_8*(cmfd % totalxs(h,i,j,k) - &
                     cmfd % p1scattxs(h,i,j,k)))

              else if (ital == 2) then

                ! Begin loop to get energy out tallies
                INGROUP: do g = 1, ng

                  ! Reset all bins to 1
                  do l = 1, size(t % filter)
                    call filter_matches(t % filter(l)) % bins % clear()
                    call filter_matches(t % filter(l)) % bins % push_back(1)
                  end do

                  ! Set ijk as mesh indices
                  ijk = (/ i, j, k /)

                  ! Get bin number for mesh indices
                  filter_matches(i_filter_mesh) % bins % data(1) = &
                       m % get_bin_from_indices(ijk)

                  if (energy_filters) then
                    ! Apply energy in filter
                    filter_matches(i_filter_ein) % bins % data(1) = ng - h + 1

                    ! Set energy out bin
                    filter_matches(i_filter_eout) % bins % data(1) = ng - g + 1
                  end if

                  ! Calculate score index from bins
                  score_index = 1
                  do l = 1, size(t % filter)
                    score_index = score_index + (filter_matches(t % filter(l)) &
                         % bins % data(1) - 1) * t % stride(l)
                  end do

                  ! Get scattering
                  cmfd % scattxs(h,g,i,j,k) = t % results(RESULT_SUM,1,score_index) /&
                       cmfd % flux(h,i,j,k)

                  ! Get nu-fission
                  cmfd % nfissxs(h,g,i,j,k) = t % results(RESULT_SUM,2,score_index) /&
                       cmfd % flux(h,i,j,k)

                  ! Bank source
                  cmfd % openmc_src(g,i,j,k) = cmfd % openmc_src(g,i,j,k) + &
                       t % results(RESULT_SUM,2,score_index)
                  cmfd % keff_bal = cmfd % keff_bal + &
                       t % results(RESULT_SUM,2,score_index) / t % n_realizations

                end do INGROUP

              else if (ital == 3) then

                ! Initialize and filter for energy
                do l = 1, size(t % filter)
                  call filter_matches(t % filter(l)) % bins % clear()
                  call filter_matches(t % filter(l)) % bins % push_back(1)
                end do

                ! Set the bin for this mesh cell
                i_mesh = m % get_bin_from_indices([ i, j, k ])
                filter_matches(i_filter_mesh) % bins % data(1) = 12*(i_mesh - 1) + 1

                ! Set the energy bin if needed
                if (energy_filters) then
                  filter_matches(i_filter_ein) % bins % data(1) = ng - h + 1
                end if

                score_index = 0
                do l = 1, size(t % filter)
                  score_index = score_index + (filter_matches(t % filter(l)) &
                       % bins % data(1) - 1) * t % stride(l)
                end do

                ! Left surface
                cmfd % current(1,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + OUT_LEFT)
                cmfd % current(2,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + IN_LEFT)

                ! Right surface
                cmfd % current(3,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + IN_RIGHT)
                cmfd % current(4,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + OUT_RIGHT)

                ! Back surface
                cmfd % current(5,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + OUT_BACK)
                cmfd % current(6,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + IN_BACK)

                ! Front surface
                cmfd % current(7,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + IN_FRONT)
                cmfd % current(8,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + OUT_FRONT)

                ! Left surface
                cmfd % current(9,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + OUT_BOTTOM)
                cmfd % current(10,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + IN_BOTTOM)

                ! Right surface
                cmfd % current(11,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + IN_TOP)
                cmfd % current(12,h,i,j,k) = t % results(RESULT_SUM, 1, &
                     score_index + OUT_TOP)
              end if TALLY

            end do OUTGROUP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

      end associate
    end do TAL

    ! Normalize openmc source distribution
    cmfd % openmc_src = cmfd % openmc_src/sum(cmfd % openmc_src)*cmfd%norm

    ! Nullify all pointers
    if (associated(m)) nullify(m)

  end subroutine compute_xs

!===============================================================================
! SET_COREMAP is a routine that sets the core mapping information
!===============================================================================

  subroutine set_coremap()

    use constants,  only: CMFD_NOACCEL

    integer :: counter=1 ! counter for unique fuel assemblies
    integer :: nx        ! number of mesh cells in x direction
    integer :: ny        ! number of mesh cells in y direction
    integer :: nz        ! number of mesh cells in z direction
    integer :: i         ! iteration counter for x
    integer :: j         ! iteration counter for y
    integer :: k         ! iteration counter for z

    ! Extract spatial indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)

    ! Count how many fuel assemblies exist
    cmfd % mat_dim = sum(cmfd % coremap - 1)

    ! Allocate indexmap
    if (.not. allocated(cmfd % indexmap)) &
         allocate(cmfd % indexmap(cmfd % mat_dim,3))

    ! Begin loops over spatial indices
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          ! Check for reflector
          if (cmfd % coremap(i,j,k) == 1) then

            ! reset value to CMFD no acceleration constant
            cmfd % coremap(i,j,k) = CMFD_NOACCEL

          else

            ! Must be a fuel --> give unique id number
            cmfd % coremap(i,j,k) = counter
            cmfd % indexmap(counter,1) = i
            cmfd % indexmap(counter,2) = j
            cmfd % indexmap(counter,3) = k
            counter = counter + 1

          end if

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine set_coremap

!===============================================================================
! NEUTRON_BALANCE computes the RMS neutron balance over the CMFD mesh
!===============================================================================

  subroutine neutron_balance()

    use constants,    only: ONE, ZERO, CMFD_NOACCEL, CMFD_NORES
    use simulation_header, only: keff, current_batch

    integer :: nx           ! number of mesh cells in x direction
    integer :: ny           ! number of mesh cells in y direction
    integer :: nz           ! number of mesh cells in z direction
    integer :: ng           ! number of energy groups
    integer :: i            ! iteration counter for x
    integer :: j            ! iteration counter for y
    integer :: k            ! iteration counter for z
    integer :: g            ! iteration counter for g
    integer :: h            ! iteration counter for outgoing groups
    integer :: l            ! iteration counter for leakage
    integer :: cnt          ! number of locations to count neutron balance
    real(8) :: leakage      ! leakage term in neutron balance
    real(8) :: interactions ! total number of interactions in balance
    real(8) :: scattering   ! scattering term in neutron balance
    real(8) :: fission      ! fission term in neutron balance
    real(8) :: res          ! residual of neutron balance
    real(8) :: rms          ! RMS of the residual

    ! Extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! Allocate res dataspace
    if (.not. allocated(cmfd%resnb)) allocate(cmfd%resnb(ng,nx,ny,nz))

    ! Reset rms and cnt
    rms = ZERO
    cnt = 0

    ! Begin loop around space and energy groups
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          GROUPG: do g = 1, ng

            ! Check for active mesh
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
                cmfd%resnb(g,i,j,k) = CMFD_NORES
                cycle
              end if
            end if

            ! Get leakage
            leakage = ZERO
            LEAK: do l = 1, 3
              leakage = leakage + ((cmfd % current(4*l,g,i,j,k) - &
                   cmfd % current(4*l-1,g,i,j,k))) - &
                   ((cmfd % current(4*l-2,g,i,j,k) - &
                   cmfd % current(4*l-3,g,i,j,k)))

            end do LEAK

            ! Interactions
            interactions = cmfd % totalxs(g,i,j,k) * cmfd % flux(g,i,j,k)

            ! Get scattering and fission
            scattering = ZERO
            fission = ZERO
            GROUPH: do h = 1, ng

              scattering = scattering + cmfd % scattxs(h,g,i,j,k) * &
                   cmfd % flux(h,i,j,k)

              fission = fission + cmfd % nfissxs(h,g,i,j,k) * &
                   cmfd % flux(h,i,j,k)

            end do GROUPH

            ! Compute residual
            res = leakage + interactions - scattering - (ONE/keff)*fission

            ! Normalize by flux
            res = res/cmfd%flux(g,i,j,k)

            ! Bank res in cmfd object
            cmfd%resnb(g,i,j,k) = res

            ! Take square for RMS calculation
            rms = rms + res**2
            cnt = cnt + 1

          end do GROUPG

        end do XLOOP

      end do YLOOP

    end do ZLOOP

    ! Calculate RMS and record in vector for this batch
    cmfd % balance(current_batch) = sqrt(ONE/dble(cnt)*rms)

  end subroutine neutron_balance

!===============================================================================
! COMPUTE_DTILDE precomputes the diffusion coupling coefficient
!===============================================================================

  subroutine compute_dtilde()

    use constants,  only: CMFD_NOACCEL, ZERO_FLUX, TINY_BIT

    integer :: nx           ! maximum number of cells in x direction
    integer :: ny           ! maximum number of cells in y direction
    integer :: nz           ! maximum number of cells in z direction
    integer :: ng           ! maximum number of energy groups
    integer :: nxyz(3,2)    ! single vector containing boundary locations
    integer :: i            ! iteration counter for x
    integer :: j            ! iteration counter for y
    integer :: k            ! iteration counter for z
    integer :: g            ! iteration counter for groups
    integer :: l            ! iteration counter for leakages
    integer :: xyz_idx      ! index for determining if x,y or z leakage
    integer :: dir_idx      ! index for determining - or + face of cell
    integer :: shift_idx    ! parameter to shift index by +1 or -1
    integer :: neig_idx(3)  ! spatial indices of neighbour
    integer :: bound(6)     ! vector containing indices for boudary check
    real(8) :: albedo(6)    ! albedo vector with global boundaries
    real(8) :: cell_dc      ! diffusion coef of current cell
    real(8) :: cell_hxyz(3) ! cell dimensions of current ijk cell
    real(8) :: neig_dc      ! diffusion coefficient of neighbor cell
    real(8) :: neig_hxyz(3) ! cell dimensions of neighbor cell
    real(8) :: dtilde       ! finite difference coupling parameter
    real(8) :: ref_albedo   ! albedo to reflector

    ! Get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! Create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! Get boundary condition information
    albedo = cmfd%albedo

    ! Loop over group and spatial indices
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          GROUP: do g = 1, ng

            ! Check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) cycle
            end if

            ! Get cell data
            cell_dc = cmfd%diffcof(g,i,j,k)
            cell_hxyz = cmfd%hxyz(:,i,j,k)

            ! Setup of vector to identify boundary conditions
            bound = (/i,i,j,j,k,k/)

            ! Begin loop around sides of cell for leakage
            LEAK: do l = 1, 6

              ! Define xyz and +/- indices
              xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
              dir_idx = 2 - mod(l,2) ! -=1, +=2
              shift_idx = -2*mod(l,2) + 1 ! shift neig by -1 or +1

              ! Check if at a boundary
              if (bound(l) == nxyz(xyz_idx,dir_idx)) then

                ! Compute dtilde with albedo boundary condition
                dtilde = (2*cell_dc*(1-albedo(l)))/(4*cell_dc*(1+albedo(l)) + &
                     (1-albedo(l))*cell_hxyz(xyz_idx))

                ! Check for zero flux
                if (abs(albedo(l) - ZERO_FLUX) < TINY_BIT) dtilde = 2*cell_dc / &
                     cell_hxyz(xyz_idx)

              else  ! not a boundary

                ! Compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! Get neigbor cell data
                neig_dc = cmfd%diffcof(g,neig_idx(1),neig_idx(2),neig_idx(3))
                neig_hxyz = cmfd%hxyz(:,neig_idx(1),neig_idx(2),neig_idx(3))

                ! Check for fuel-reflector interface
                if (cmfd_coremap) then

                  if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) == &
                       CMFD_NOACCEL .and. cmfd % coremap(i,j,k) /= CMFD_NOACCEL) then

                    ! Get albedo
                    ref_albedo = get_reflector_albedo(l,g,i,j,k)

                    ! Compute dtilde
                    dtilde = (2*cell_dc*(1-ref_albedo))/(4*cell_dc*(1+ &
                         ref_albedo)+(1-ref_albedo)*cell_hxyz(xyz_idx))

                  else ! Not next to a reflector or no core map

                    ! Compute dtilde to neighbor cell
                    dtilde = (2*cell_dc*neig_dc)/(neig_hxyz(xyz_idx)*cell_dc + &
                         cell_hxyz(xyz_idx)*neig_dc)

                  end if

                else ! no core map

                  ! Compute dtilde to neighbor cell
                  dtilde = (2*cell_dc*neig_dc)/(neig_hxyz(xyz_idx)*cell_dc + &
                       cell_hxyz(xyz_idx)*neig_dc)

                end if

              end if

              ! Record dtilde in cmfd object
              cmfd%dtilde(l,g,i,j,k) = dtilde

            end do LEAK

          end do GROUP

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine compute_dtilde

!===============================================================================
! COMPUTE_DHAT computes the nonlinear coupling coefficient
!===============================================================================

  subroutine compute_dhat()

    use constants, only: CMFD_NOACCEL, ZERO
    use error,     only: write_message
    use string,    only: to_str

    integer :: nx             ! maximum number of cells in x direction
    integer :: ny             ! maximum number of cells in y direction
    integer :: nz             ! maximum number of cells in z direction
    integer :: ng             ! maximum number of energy groups
    integer :: nxyz(3,2)      ! single vector containing boundary locations
    integer :: i              ! iteration counter for x
    integer :: j              ! iteration counter for y
    integer :: k              ! iteration counter for z
    integer :: g              ! iteration counter for groups
    integer :: l              ! iteration counter for leakages
    integer :: xyz_idx        ! index for determining if x,y or z leakage
    integer :: dir_idx        ! index for determining - or + face of cell
    integer :: shift_idx      ! parameter to shift index by +1 or -1
    integer :: neig_idx(3)    ! spatial indices of neighbour
    integer :: bound(6)       ! vector containing indices for boudary check
    real(8) :: cell_dtilde(6) ! cell dtilde for each face
    real(8) :: cell_flux      ! flux in current cell
    real(8) :: current(12)    ! area integrated cell current at each face
    real(8) :: net_current    ! net current on a face
    real(8) :: neig_flux      ! flux in neighbor cell
    real(8) :: dhat           ! dhat equivalence parameter

    ! Get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! Create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! Geting loop over group and spatial indices
    ZLOOP:  do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          GROUP: do g = 1,ng

            ! Check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
                cycle
              end if
            end if

            ! Get cell data
            cell_dtilde = cmfd%dtilde(:,g,i,j,k)
            cell_flux = cmfd%flux(g,i,j,k)/product(cmfd%hxyz(:,i,j,k))
            current = cmfd%current(:,g,i,j,k)

            ! Setup of vector to identify boundary conditions
            bound = (/i,i,j,j,k,k/)

            ! Begin loop around sides of cell for leakage
            LEAK: do l = 1,6

              ! Define xyz and +/- indices
              xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
              dir_idx = 2 - mod(l,2) ! -=1, +=2
              shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1

              ! Calculate net current on l face (divided by surf area)
              net_current = (current(2*l) - current(2*l-1)) / &
                   product(cmfd%hxyz(:,i,j,k)) * cmfd%hxyz(xyz_idx,i,j,k)

              ! Check if at a boundary
              if (bound(l) == nxyz(xyz_idx,dir_idx)) then

                ! Compute dhat
                dhat = (net_current - shift_idx*cell_dtilde(l)*cell_flux) / &
                     cell_flux

              else  ! not a boundary

                ! Compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! Get neigbor flux
                neig_flux = cmfd%flux(g,neig_idx(1),neig_idx(2),neig_idx(3)) / &
                     product(cmfd%hxyz(:,neig_idx(1),neig_idx(2),neig_idx(3)))

                ! Check for fuel-reflector interface
                if (cmfd_coremap) then

                  if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) == &
                       CMFD_NOACCEL .and. cmfd % coremap(i,j,k) /= CMFD_NOACCEL) then

                    ! compute dhat
                    dhat = (net_current - shift_idx*cell_dtilde(l)*cell_flux) /&
                         cell_flux

                  else ! not a fuel-reflector interface

                    ! Compute dhat
                    dhat = (net_current + shift_idx*cell_dtilde(l)* &
                         (neig_flux - cell_flux))/(neig_flux + cell_flux)

                  end if

                else ! not for fuel-reflector case

                  ! Compute dhat
                  dhat = (net_current + shift_idx*cell_dtilde(l)* &
                       (neig_flux - cell_flux))/(neig_flux + cell_flux)

                end if

              end if

              ! record dhat in cmfd object
              cmfd%dhat(l,g,i,j,k) = dhat

              ! check for dhat reset
              if (dhat_reset) then
                cmfd%dhat(l,g,i,j,k) = ZERO
              end if

            end do LEAK

          end do GROUP

        end do XLOOP

      end do YLOOP

    end do ZLOOP

    ! write that dhats are zero
    if (dhat_reset) then
      call write_message('Dhats reset to zero.', 8)
    end if

  end subroutine compute_dhat

!===============================================================================
! GET_REFLECTOR_ALBEDO is a function that calculates the albedo to the reflector
!===============================================================================

  function get_reflector_albedo(l, g, i, j, k)

    use constants,  only: ONE

    real(8) :: get_reflector_albedo ! reflector albedo
    integer, intent(in) :: i ! iteration counter for x
    integer, intent(in) :: j ! iteration counter for y
    integer, intent(in) :: k ! iteration counter for z
    integer, intent(in) :: g ! iteration counter for groups
    integer, intent(in) :: l ! iteration counter for leakages

    integer :: shift_idx   ! parameter to shift index by +1 or -1
    real(8) :: current(12) ! partial currents for all faces of mesh cell
    real(8) :: albedo      ! the albedo

    ! Get partial currents from object
    current = cmfd%current(:,g,i,j,k)

    ! Define xyz and +/- indices
    shift_idx = -2*mod(l,2) + 1          ! shift neig by -1 or +1

    ! Calculate albedo
    if ((shift_idx ==  1 .and. current(2*l  ) < 1.0e-10_8) .or. &
         (shift_idx == -1 .and. current(2*l-1) < 1.0e-10_8)) then
      albedo = ONE
    else
      albedo = (current(2*l-1)/current(2*l))**(shift_idx)
    end if

    ! Assign to function variable
    get_reflector_albedo = albedo

  end function get_reflector_albedo

!===============================================================================
! COMPUTE_EFFECTIVE_DOWNSCATTER changes downscatter rate for zero upscatter
!===============================================================================

  subroutine compute_effective_downscatter()

    use constants, only: ZERO, CMFD_NOACCEL

    integer :: nx                ! number of mesh cells in x direction
    integer :: ny                ! number of mesh cells in y direction
    integer :: nz                ! number of mesh cells in z direction
    integer :: ng                ! number of energy groups
    integer :: i                 ! iteration counter for x
    integer :: j                 ! iteration counter for y
    integer :: k                 ! iteration counter for z
    real(8) :: flux1             ! group 1 volume int flux
    real(8) :: flux2             ! group 2 volume int flux
    real(8) :: sigt1             ! group 1 total xs
    real(8) :: sigt2             ! group 2 total xs
    real(8) :: sigs11            ! scattering transfer 1 --> 1
    real(8) :: sigs21            ! scattering transfer 2 --> 1
    real(8) :: sigs12            ! scattering transfer 1 --> 2
    real(8) :: sigs22            ! scattering transfer 2 --> 2
    real(8) :: siga1             ! group 1 abs xs
    real(8) :: siga2             ! group 2 abs xs
    real(8) :: sigs12_eff        ! effective downscatter xs

    ! Extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! Return if not two groups
    if (ng /= 2) return

    ! Begin loop around space and energy groups
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          ! Check for active mesh
          if (allocated(cmfd%coremap)) then
            if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) cycle
          end if

          ! Extract cross sections and flux from object
          flux1 = cmfd % flux(1,i,j,k)
          flux2 = cmfd % flux(2,i,j,k)
          sigt1 = cmfd % totalxs(1,i,j,k)
          sigt2 = cmfd % totalxs(2,i,j,k)
          sigs11 = cmfd % scattxs(1,1,i,j,k)
          sigs21 = cmfd % scattxs(2,1,i,j,k)
          sigs12 = cmfd % scattxs(1,2,i,j,k)
          sigs22 = cmfd % scattxs(2,2,i,j,k)

          ! Compute absorption xs
          siga1 = sigt1 - sigs11 - sigs12
          siga2 = sigt2 - sigs22 - sigs21

          ! Compute effective downscatter xs
          sigs12_eff = sigs12 - sigs21*flux2/flux1

          ! Recompute total cross sections (use effective and no upscattering)
          sigt1 = siga1 + sigs11 + sigs12_eff
          sigt2 = siga2 + sigs22

          ! Record total xs
          cmfd % totalxs(1,i,j,k) = sigt1
          cmfd % totalxs(2,i,j,k) = sigt2

          ! Record effective downscatter xs
          cmfd % scattxs(1,2,i,j,k) = sigs12_eff

          ! Zero out upscatter cross section
          cmfd % scattxs(2,1,i,j,k) = ZERO

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine compute_effective_downscatter

end module cmfd_data
