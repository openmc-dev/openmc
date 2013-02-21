module cmfd_data

!==============================================================================
! CMFD_DATA -- This module processes the cmfd tally object to generate
! parameters for CMFD calculation.
!==============================================================================


  implicit none
  private
  public :: set_up_cmfd, neutron_balance

  logical :: dhat_reset = .false.

contains

!==============================================================================
! SET_UP_CMFD
!==============================================================================

  subroutine set_up_cmfd() 

    use cmfd_header,         only: allocate_cmfd
    use constants,           only: CMFD_NOACCEL
    use global,              only: cmfd, cmfd_coremap, cmfd_run_2grp

    ! initialize cmfd object
    if (.not.allocated(cmfd%flux)) call allocate_cmfd(cmfd)

    ! check for core map and set it up
    if ((cmfd_coremap) .and. (cmfd%mat_dim == CMFD_NOACCEL)) call set_coremap()

    ! calculate all cross sections based on reaction rates from last batch
    if (.not. cmfd_run_2grp) call compute_xs()

    ! check neutron balance
!   call neutron_balance(670)

    ! fix 2 grp cross sections
    if (cmfd_run_2grp) call fix_2_grp()

    ! calculate dtilde
    call compute_dtilde()

    ! calucate dhat
    call compute_dhat()

  end subroutine set_up_cmfd 

!===============================================================================
! COMPUTE_XS takes tallies and computes macroscopic cross sections
!===============================================================================

  subroutine compute_xs()

    use constants,    only: FILTER_MESH, FILTER_ENERGYIN, FILTER_ENERGYOUT,     &
                            FILTER_SURFACE, IN_RIGHT, OUT_RIGHT, IN_FRONT,      &
                            OUT_FRONT, IN_TOP, OUT_TOP, CMFD_NOACCEL, ZERO, ONE
    use error,        only: fatal_error
    use global,       only: cmfd, message, n_cmfd_tallies, cmfd_tallies, meshes
    use mesh,         only: mesh_indices_to_bin
    use mesh_header,  only: StructuredMesh
    use tally_header, only: TallyObject

    integer :: nx            ! number of mesh cells in x direction
    integer :: ny            ! number of mesh cells in y direction
    integer :: nz            ! number of mesh cells in z direction
    integer :: ng            ! number of energy groups
    integer :: i             ! iteration counter for x
    integer :: j             ! iteration counter for y
    integer :: k             ! iteration counter for z
    integer :: g             ! iteration counter for g
    integer :: h             ! iteration counter for outgoing groups
    integer :: ital          ! tally object index
    integer :: ijk(3)        ! indices for mesh cell
    integer :: score_index   ! index to pull from tally object
    integer :: i_mesh        ! index in meshes array
    integer :: i_filter_mesh ! index for mesh filter
    integer :: i_filter_ein  ! index for incoming energy filter
    integer :: i_filter_eout ! index for outgoing energy filter
    integer :: i_filter_surf ! index for surface filter
    real(8) :: flux          ! temp variable for flux
    type(TallyObject),    pointer :: t => null() ! pointer for tally object
    type(StructuredMesh), pointer :: m => null() ! pointer for mesh object

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! set flux object and source distribution to all zeros
    cmfd % flux = ZERO
    cmfd % openmc_src = ZERO

    ! associate tallies and mesh
    t => cmfd_tallies(1)
    i_mesh = t % filters(t % find_filter(FILTER_MESH)) % int_bins(1)
    m => meshes(i_mesh)

    ! set mesh widths
    cmfd % hxyz(1,:,:,:) = m % width(1) ! set x width
    cmfd % hxyz(2,:,:,:) = m % width(2) ! set y width
    cmfd % hxyz(3,:,:,:) = m % width(3) ! set z width

   ! begin loop around tallies
   TAL: do ital = 1, n_cmfd_tallies

     ! associate tallies and mesh
     t => cmfd_tallies(ital)
     i_mesh = t % filters(t % find_filter(FILTER_MESH)) % int_bins(1)
     m => meshes(i_mesh)

     i_filter_mesh = t % find_filter(FILTER_MESH)
     i_filter_ein  = t % find_filter(FILTER_ENERGYIN)
     i_filter_eout = t % find_filter(FILTER_ENERGYOUT)
     i_filter_surf = t % find_filter(FILTER_SURFACE)

     ! begin loop around space
     ZLOOP: do k = 1,nz

       YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx
 
            ! check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
                cycle
              end if
            end if

            ! loop around energy groups
            OUTGROUP: do h = 1,ng

              ! start tally 1
              TALLY: if (ital == 1) then

                ! reset all bins to 1
                t % matching_bins = 1

                ! set ijk as mesh indices
                ijk = (/ i, j, k /)

                ! get bin number for mesh indices
                t % matching_bins(i_filter_mesh) = mesh_indices_to_bin(m,ijk)

                ! apply energy in filter
                if (i_filter_ein > 0) then
                  t % matching_bins(i_filter_ein) = ng - h + 1
                end if

                ! calculate score index from bins
                score_index = sum((t % matching_bins - 1) * t%stride) + 1

                ! get flux
                flux = t % results(1,score_index) % sum
                cmfd % flux(h,i,j,k) = flux

                ! detect zero flux
                if ((flux - 0.0D0) < 1.0E-10_8) then
                  print *,h,i,j,k,flux
                  message = 'Detected zero flux without coremap overlay'
                  call fatal_error()
                end if

                ! get total rr and convert to total xs
                cmfd % totalxs(h,i,j,k) = t % results(2,score_index) % sum / flux

                ! get p1 scatter rr and convert to p1 scatter xs
                cmfd % p1scattxs(h,i,j,k) = t % results(3,score_index) % sum / flux

                ! extract diffusion coefficient tally
                cmfd % diffusion(h,i,j,k) = t % results(4,score_index) % sum / flux

                ! calculate diffusion coefficient
!               cmfd % diffcof(h,i,j,k) = ONE/(3.0_8*cmfd%totalxs(h,i,j,k))
                cmfd % diffcof(h,i,j,k) = ONE/(3.0_8*(cmfd % totalxs(h,i,j,k) - &
                     cmfd % p1scattxs(h,i,j,k)))
!               cmfd % diffcof(h,i,j,k) = cmfd % diffusion(h,i,j,k)

              else if (ital == 2) then

                ! begin loop to get energy out tallies
                INGROUP: do g = 1, ng

                  ! reset all bins to 1
                  t % matching_bins = 1

                  ! set ijk as mesh indices
                  ijk = (/ i, j, k /)

                  ! get bin number for mesh indices
                  t % matching_bins(i_filter_mesh) = mesh_indices_to_bin(m,ijk)

                  if (i_filter_ein > 0) then
                    ! apply energy in filter
                    t % matching_bins(i_filter_ein) = ng - h + 1

                    ! set energy out bin
                    t % matching_bins(i_filter_eout) = ng - g + 1
                  end if

                  ! calculate score index from bins
                  score_index = sum((t % matching_bins - 1) * t%stride) + 1

                  ! get scattering
                  cmfd % scattxs(h,g,i,j,k) = t % results(1,score_index) % sum /&
                       cmfd % flux(h,i,j,k)

                  ! get nu-fission
                  cmfd % nfissxs(h,g,i,j,k) = t % results(2,score_index) % sum /&
                       cmfd % flux(h,i,j,k)

                  ! bank source
                  cmfd % openmc_src(g,i,j,k) = cmfd % openmc_src(g,i,j,k) + &
                       t % results(2,score_index) % sum

                end do INGROUP

              else if (ital == 3) then

                ! initialize and filter for energy
                t % matching_bins = 1
                if (i_filter_ein > 0) then
                  t % matching_bins(i_filter_ein) = ng - h + 1
                end if

                ! left surface
                t % matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i-1, j, k /) + 1, .true.)
                t % matching_bins(i_filter_surf) = IN_RIGHT
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! outgoing
                cmfd % current(1,h,i,j,k) = t % results(1,score_index) % sum
                t % matching_bins(i_filter_surf) = OUT_RIGHT
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! incoming
                cmfd % current(2,h,i,j,k) = t % results(1,score_index) % sum

                ! right surface
                t % matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j, k /) + 1, .true.)
                t % matching_bins(i_filter_surf) = IN_RIGHT
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! incoming
                cmfd % current(3,h,i,j,k) = t % results(1,score_index) % sum
                t % matching_bins(i_filter_surf) = OUT_RIGHT
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! outgoing
                cmfd % current(4,h,i,j,k) = t % results(1,score_index) % sum

                ! back surface
                t % matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j-1, k /) + 1, .true.)
                t % matching_bins(i_filter_surf) = IN_FRONT
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! outgoing
                cmfd % current(5,h,i,j,k) = t % results(1,score_index) % sum
                t % matching_bins(i_filter_surf) = OUT_FRONT
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! incoming
                cmfd % current(6,h,i,j,k) = t % results(1,score_index) % sum

                ! front surface
                t % matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j, k /) + 1, .true.)
                t % matching_bins(i_filter_surf) = IN_FRONT
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! incoming
                cmfd % current(7,h,i,j,k) = t % results(1,score_index) % sum
                t % matching_bins(i_filter_surf) = OUT_FRONT
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! outgoing
                cmfd % current(8,h,i,j,k) = t % results(1,score_index) % sum

                ! bottom surface
                t % matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j, k-1 /) + 1, .true.)
                t % matching_bins(i_filter_surf) = IN_TOP
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! outgoing
                cmfd % current(9,h,i,j,k) = t % results(1,score_index) % sum
                t % matching_bins(i_filter_surf) = OUT_TOP
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! incoming
                cmfd % current(10,h,i,j,k) = t % results(1,score_index) % sum

                ! top surface
                t % matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j, k /) + 1, .true.)
                t % matching_bins(i_filter_surf) = IN_TOP
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! incoming
                cmfd % current(11,h,i,j,k) = t % results(1,score_index) % sum
                t % matching_bins(i_filter_surf) = OUT_TOP
                score_index = sum((t % matching_bins - 1) * t % stride) + 1 ! outgoing
                cmfd % current(12,h,i,j,k) = t % results(1,score_index) % sum

              end if TALLY

            end do OUTGROUP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

    end do TAL

    ! normalize openmc source distribution
    cmfd % openmc_src = cmfd % openmc_src/sum(cmfd % openmc_src)*cmfd%norm

    ! nullify all pointers
    if (associated(t)) nullify(t)
    if (associated(m)) nullify(m)

  end subroutine compute_xs

!===============================================================================
! SET_COREMAP is a routine that sets the core mapping information
!===============================================================================

  subroutine set_coremap()

    use constants,  only: CMFD_NOACCEL
    use global,     only: cmfd

    integer :: kount=1           ! counter for unique fuel assemblies
    integer :: nx                ! number of mesh cells in x direction
    integer :: ny                ! number of mesh cells in y direction
    integer :: nz                ! number of mesh cells in z direction
    integer :: i                 ! iteration counter for x
    integer :: j                 ! iteration counter for y
    integer :: k                 ! iteration counter for z

    ! extract spatial indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)

    ! count how many fuel assemblies exist
    cmfd % mat_dim = sum(cmfd % coremap - 1)

    ! allocate indexmap
    if (.not. allocated(cmfd % indexmap)) &
         allocate(cmfd % indexmap(cmfd % mat_dim,3))

    ! begin loops over spatial indices
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          ! check for reflector
          if (cmfd % coremap(i,j,k) == 1) then

            ! reset value to 99999
            cmfd % coremap(i,j,k) = CMFD_NOACCEL 

          else

            ! must be a fuel --> give unique id number
            cmfd % coremap(i,j,k) = kount
            cmfd % indexmap(kount,1) = i
            cmfd % indexmap(kount,2) = j
            cmfd % indexmap(kount,3) = k
            kount = kount + 1

          end if

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine set_coremap

!===============================================================================
! NEUTRON_BALANCE writes a file that contains n. bal. info for all cmfd mesh
!===============================================================================

  subroutine neutron_balance(uid)

    use constants,    only: ONE, ZERO, CMFD_NOACCEL, CMFD_NORES
    use global,       only: cmfd, keff

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
    integer :: uid
    real(8) :: leakage      ! leakage term in neutron balance
    real(8) :: interactions ! total number of interactions in balance
    real(8) :: scattering   ! scattering term in neutron balance
    real(8) :: fission      ! fission term in neutron balance
    real(8) :: res          ! residual of neutron balance

    ! check if keff is close to 0 (happens on first active batch)
    if (keff < 1e-8_8) return

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! allocate res dataspace
    if (.not. allocated(cmfd%resnb)) allocate(cmfd%resnb(ng,nx,ny,nz))

    ! begin loop around space and energy groups
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          GROUPG: do g = 1, ng

            ! check for active mesh
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
                cmfd%resnb(g,i,j,k) = CMFD_NORES
                cycle
              end if
            end if

            ! get leakage
            leakage = ZERO
            LEAK: do l = 1, 3

              leakage = leakage + ((cmfd % current(4*l,g,i,j,k) - &
                   cmfd % current(4*l-1,g,i,j,k))) - &
                   ((cmfd % current(4*l-2,g,i,j,k) - &
                   cmfd % current(4*l-3,g,i,j,k)))

            end do LEAK

            ! interactions
            interactions = cmfd % totalxs(g,i,j,k) * cmfd % flux(g,i,j,k)

            ! get scattering and fission
            scattering = ZERO
            fission = ZERO
            GROUPH: do h = 1, ng

              scattering = scattering + cmfd % scattxs(h,g,i,j,k) * &
                   cmfd % flux(h,i,j,k)

              fission = fission + cmfd % nfissxs(h,g,i,j,k) * &
                   cmfd % flux(h,i,j,k)

            end do GROUPH

            ! compute residual
            res = leakage + interactions - scattering - (ONE/keff)*fission

            ! normalize by flux
            res = res/cmfd%flux(g,i,j,k)

            ! bank res in cmfd object
            cmfd%resnb(g,i,j,k) = res

            ! write out info to file
            write(uid,'(A,1X,I0,1X,I0,1X,I0,1X,A,1X,I0,1X,A,1PE11.4)') &
                          'Location',i,j,k,' Group:',g,'Balance:',res
            write(uid,100) 'Leakage:',leakage
            write(uid,100) 'Interactions:',interactions
            write(uid,100) 'Scattering:',scattering
            write(uid,100) 'Fission:',fission
            write(uid,100) 'k-eff:',keff
            write(uid,100) 'k-eff balanced:',fission/(leakage+interactions-scattering)
            write(uid,100) 'Balance:',keff - fission/(leakage+interactions-scattering)

          end do GROUPG

        end do XLOOP

      end do YLOOP

    end do ZLOOP

 100 FORMAT(A,1X,1PE11.4)

  end subroutine neutron_balance

!===============================================================================
! COMPUTE_DTILDE computes the diffusion coupling coefficient
!===============================================================================

  subroutine compute_dtilde()

    use constants,  only: CMFD_NOACCEL, ZERO_FLUX, TINY_BIT
    use global,     only: cmfd, cmfd_coremap

    integer :: nx                 ! maximum number of cells in x direction
    integer :: ny                 ! maximum number of cells in y direction
    integer :: nz                 ! maximum number of cells in z direction
    integer :: ng                 ! maximum number of energy groups
    integer :: nxyz(3,2)          ! single vector containing boundary locations
    integer :: i                  ! iteration counter for x
    integer :: j                  ! iteration counter for y
    integer :: k                  ! iteration counter for z
    integer :: g                  ! iteration counter for groups
    integer :: l                  ! iteration counter for leakages
    integer :: xyz_idx            ! index for determining if x,y or z leakage
    integer :: dir_idx            ! index for determining - or + face of cell
    integer :: shift_idx          ! parameter to shift index by +1 or -1
    integer :: neig_idx(3)        ! spatial indices of neighbour
    integer :: bound(6)           ! vector containing indices for boudary check
    real(8) :: albedo(6)          ! albedo vector with global boundaries
    real(8) :: cell_dc            ! diffusion coef of current cell
    real(8) :: cell_hxyz(3)       ! cell dimensions of current ijk cell
    real(8) :: neig_dc            ! diffusion coefficient of neighbor cell
    real(8) :: neig_hxyz(3)       ! cell dimensions of neighbor cell
    real(8) :: dtilde             ! finite difference coupling parameter 
    real(8) :: ref_albedo         ! albedo to reflector

    ! get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! get boundary condition information
    albedo = cmfd%albedo

    ! geting loop over group and spatial indices
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          GROUP: do g = 1, ng

            ! check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) cycle
            end if

            ! get cell data
            cell_dc = cmfd%diffcof(g,i,j,k)
            cell_hxyz = cmfd%hxyz(:,i,j,k)

            ! setup of vector to identify boundary conditions
            bound = (/i,i,j,j,k,k/)

            ! begin loop around sides of cell for leakage
            LEAK: do l = 1, 6

              ! define xyz and +/- indices
              xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
              dir_idx = 2 - mod(l,2) ! -=1, +=2
              shift_idx = -2*mod(l,2) + 1          ! shift neig by -1 or +1

              ! check if at a boundary
              if (bound(l) == nxyz(xyz_idx,dir_idx)) then

                ! compute dtilde
                dtilde = (2*cell_dc*(1-albedo(l)))/(4*cell_dc*(1+albedo(l)) + &
                     (1-albedo(l))*cell_hxyz(xyz_idx))

                ! check for zero flux
                if (abs(albedo(l) - ZERO_FLUX) < TINY_BIT) dtilde = 2*cell_dc / &
                     cell_hxyz(xyz_idx)

              else  ! not a boundary

                ! compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! get neigbor cell data
                neig_dc = cmfd%diffcof(g,neig_idx(1),neig_idx(2),neig_idx(3))
                neig_hxyz = cmfd%hxyz(:,neig_idx(1),neig_idx(2),neig_idx(3))

                ! check for fuel-reflector interface
                if (cmfd_coremap) then

                  if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) == &
                       CMFD_NOACCEL .and. cmfd % coremap(i,j,k) /= CMFD_NOACCEL) then

                    ! get albedo
                    ref_albedo = get_reflector_albedo(l,g,i,j,k)

                    ! compute dtilde
                    dtilde = (2*cell_dc*(1-ref_albedo))/(4*cell_dc*(1+ &
                         ref_albedo)+(1-ref_albedo)*cell_hxyz(xyz_idx))

                  else ! not next to a reflector or no core map

                    ! compute dtilde
                    dtilde = (2*cell_dc*neig_dc)/(neig_hxyz(xyz_idx)*cell_dc + &
                         cell_hxyz(xyz_idx)*neig_dc)

                  end if

                else ! no core map

                  ! compute dtilde
                  dtilde = (2*cell_dc*neig_dc)/(neig_hxyz(xyz_idx)*cell_dc + &
                       cell_hxyz(xyz_idx)*neig_dc)

               end if

              end if

              ! record dtilde in cmfd object
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

    use constants,  only: CMFD_NOACCEL, ZERO
    use global,     only: cmfd, cmfd_coremap

    integer :: nx                 ! maximum number of cells in x direction
    integer :: ny                 ! maximum number of cells in y direction
    integer :: nz                 ! maximum number of cells in z direction
    integer :: ng                 ! maximum number of energy groups
    integer :: nxyz(3,2)          ! single vector containing boundary locations
    integer :: i                  ! iteration counter for x
    integer :: j                  ! iteration counter for y
    integer :: k                  ! iteration counter for z
    integer :: g                  ! iteration counter for groups
    integer :: l                  ! iteration counter for leakages
    integer :: xyz_idx            ! index for determining if x,y or z leakage
    integer :: dir_idx            ! index for determining - or + face of cell
    integer :: shift_idx          ! parameter to shift index by +1 or -1
    integer :: neig_idx(3)        ! spatial indices of neighbour
    integer :: bound(6)           ! vector containing indices for boudary check
    real(8) :: cell_dtilde(6)     ! cell dtilde for each face
    real(8) :: cell_flux          ! flux in current cell
    real(8) :: current(12)        ! area integrated cell current at each face
    real(8) :: net_current        ! net current on a face
    real(8) :: neig_flux          ! flux in neighbor cell
    real(8) :: dhat               ! dhat equivalence parameter

    ! get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! geting loop over group and spatial indices
    ZLOOP:  do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          GROUP: do g = 1,ng

            ! check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
                cycle
              end if
            end if

            ! get cell data
            cell_dtilde = cmfd%dtilde(:,g,i,j,k)
            cell_flux = cmfd%flux(g,i,j,k)/product(cmfd%hxyz(:,i,j,k))
            current = cmfd%current(:,g,i,j,k)

            ! setup of vector to identify boundary conditions
            bound = (/i,i,j,j,k,k/)

            ! begin loop around sides of cell for leakage
            LEAK: do l = 1,6

              ! define xyz and +/- indices
              xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
              dir_idx = 2 - mod(l,2) ! -=1, +=2
              shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1

              ! calculate net current on l face (divided by surf area)
              net_current = (current(2*l) - current(2*l-1)) / &
                   product(cmfd%hxyz(:,i,j,k)) * cmfd%hxyz(xyz_idx,i,j,k)

              ! check if at a boundary
              if (bound(l) == nxyz(xyz_idx,dir_idx)) then

                ! compute dhat
                dhat = (net_current - shift_idx*cell_dtilde(l)*cell_flux) / &
                     cell_flux

              else  ! not a boundary

                ! compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! get neigbor flux 
                neig_flux = cmfd%flux(g,neig_idx(1),neig_idx(2),neig_idx(3)) / &
                     product(cmfd%hxyz(:,neig_idx(1),neig_idx(2),neig_idx(3)))

                ! check for fuel-reflector interface
                if (cmfd_coremap) then

                  if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) == &
                       CMFD_NOACCEL .and. cmfd % coremap(i,j,k) /= CMFD_NOACCEL) then

                    ! compute dhat
                    dhat = (net_current - shift_idx*cell_dtilde(l)*cell_flux) /&
                         cell_flux

                  else ! not a fuel-reflector interface

                    ! compute dhat 
                    dhat = (net_current + shift_idx*cell_dtilde(l)* &
                         (neig_flux - cell_flux))/(neig_flux + cell_flux)

                  end if

                else ! not for fuel-reflector case

                  ! compute dhat 
                  dhat = (net_current + shift_idx*cell_dtilde(l)* &
                       (neig_flux - cell_flux))/(neig_flux + cell_flux)

                end if

              end if

              ! check for zero net current
!             if ((abs(current(2*l)-current(2*l-1)) < 1e-8_8).and.xyz_idx/=3) then
!               print *,'Zero net current interface',g,i,j,k
!               print *,current(2*l),current(2*l-1),net_current
!               dhat = ZERO
!             end if

              ! record dhat in cmfd object
              cmfd%dhat(l,g,i,j,k) = dhat

              ! check for dhat reset
              if (dhat_reset) cmfd%dhat(l,g,i,j,k) = ZERO

            end do LEAK

          end do GROUP

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine compute_dhat

!===============================================================================
! GET_REFLECTOR_ALBEDO is a function that calculates the albedo to the reflector
!===============================================================================

  function get_reflector_albedo(l, g, i, j, k)

    use constants,  only: ALBEDO_REJECT
    use global,     only: cmfd, cmfd_hold_weights

    real(8) :: get_reflector_albedo ! reflector albedo
    integer :: i                    ! iteration counter for x
    integer :: j                    ! iteration counter for y
    integer :: k                    ! iteration counter for z
    integer :: g                    ! iteration counter for groups
    integer :: l                    ! iteration counter for leakages

    integer :: shift_idx            ! parameter to shift index by +1 or -1
    real(8) :: current(12)          ! partial currents for all faces of mesh cell            
    real(8) :: albedo               ! the albedo

    ! get partial currents from object
    current = cmfd%current(:,g,i,j,k)

    ! define xyz and +/- indices
    shift_idx = -2*mod(l,2) + 1          ! shift neig by -1 or +1

    ! calculate albedo
    if ((shift_idx ==  1 .and. current(2*l  ) < 1.0e-10_8) .or. &
        (shift_idx == -1 .and. current(2*l-1) < 1.0e-10_8)) then
      albedo = ALBEDO_REJECT 
      cmfd_hold_weights = .true. 
    else
      albedo = (current(2*l-1)/current(2*l))**(shift_idx)
    end if

    ! assign to function variable
    get_reflector_albedo = albedo

  end function get_reflector_albedo

!===============================================================================
! FIX_2_GRP modifies xs for benchmark with stand alone (homogeneous)
!===============================================================================

  subroutine fix_2_grp()

    use global,  only: cmfd

    ! overwrite cross sections
    cmfd % totalxs(1,:,:,:) = 0.02597_8
    cmfd % totalxs(2,:,:,:) = 0.06669_8
    cmfd % scattxs(1,1,:,:,:) = 0.0_8
    cmfd % scattxs(1,2,:,:,:) = 0.01742_8
    cmfd % scattxs(2,1,:,:,:) = 0.0_8
    cmfd % scattxs(2,2,:,:,:) = 0.0_8
    cmfd % nfissxs(1,1,:,:,:) = 0.00536_8
    cmfd % nfissxs(1,2,:,:,:) = 0.0_8
    cmfd % nfissxs(2,1,:,:,:) = 0.10433_8
    cmfd % nfissxs(2,2,:,:,:) = 0.0_8
    cmfd % diffcof(1,:,:,:) = 1.4176_8
    cmfd % diffcof(2,:,:,:) = 0.37336_8
    cmfd % hxyz(1,:,:,:) = 0.5_8
    cmfd % hxyz(2,:,:,:) = 0.5_8
    cmfd % hxyz(3,:,:,:) = 0.5_8

    ! set dhat reset to true
    dhat_reset = .true.

  end subroutine fix_2_grp

!===============================================================================
! FIX_NEUTRON_BALANCE
!===============================================================================

  subroutine fix_neutron_balance()

    use constants,  only: ONE, ZERO, CMFD_NOACCEL
    use global,     only: cmfd, keff
    use, intrinsic :: ISO_FORTRAN_ENV

    integer :: nx                ! number of mesh cells in x direction
    integer :: ny                ! number of mesh cells in y direction
    integer :: nz                ! number of mesh cells in z direction
    integer :: ng                ! number of energy groups
    integer :: i                 ! iteration counter for x
    integer :: j                 ! iteration counter for y
    integer :: k                 ! iteration counter for z
    integer :: l                 ! iteration counter for surface
    real(8) :: leak1             ! leakage rate in group 1
    real(8) :: leak2             ! leakage rate in group 2
    real(8) :: flux1             ! group 1 volume int flux
    real(8) :: flux2             ! group 2 volume int flux
    real(8) :: sigt1             ! group 1 total xs
    real(8) :: sigt2             ! group 2 total xs
    real(8) :: sigs11            ! scattering transfer 1 --> 1
    real(8) :: sigs21            ! scattering transfer 2 --> 1
    real(8) :: sigs12            ! scattering transfer 1 --> 2
    real(8) :: sigs22            ! scattering transfer 2 --> 2
    real(8) :: nsigf11           ! fission transfer 1 --> 1
    real(8) :: nsigf21           ! fission transfer 2 --> 1
    real(8) :: nsigf12           ! fission transfer 1 --> 2
    real(8) :: nsigf22           ! fission transfer 2 --> 2
    real(8) :: siga1             ! group 1 abs xs
    real(8) :: siga2             ! group 2 abs xs
    real(8) :: sigs12_eff        ! effective downscatter xs

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! return if not two groups
    if (ng /= 2) return

    ! begin loop around space and energy groups
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          ! check for active mesh
          if (allocated(cmfd%coremap)) then
            if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) cycle 
          end if

          ! compute leakage in groups 1 and 2
          leak1 = ZERO 
          leak2 = ZERO
          LEAK: do l = 1, 3

            leak1 = leak1 + ((cmfd % current(4*l,1,i,j,k) - &
                 cmfd % current(4*l-1,1,i,j,k))) - &
                 ((cmfd % current(4*l-2,1,i,j,k) - &
                 cmfd % current(4*l-3,1,i,j,k)))

            leak2 = leak2 + ((cmfd % current(4*l,2,i,j,k) - &
                 cmfd % current(4*l-1,2,i,j,k))) - &
                 ((cmfd % current(4*l-2,2,i,j,k) - &
                 cmfd % current(4*l-3,2,i,j,k)))


          end do LEAK

          ! extract cross sections and flux from object
          flux1 = cmfd % flux(1,i,j,k)
          flux2 = cmfd % flux(2,i,j,k)
          sigt1 = cmfd % totalxs(1,i,j,k)
          sigt2 = cmfd % totalxs(2,i,j,k)
          sigs11 = cmfd % scattxs(1,1,i,j,k)
          sigs21 = cmfd % scattxs(2,1,i,j,k)
          sigs12 = cmfd % scattxs(1,2,i,j,k)
          sigs22 = cmfd % scattxs(2,2,i,j,k)
          nsigf11 = cmfd % nfissxs(1,1,i,j,k)
          nsigf21 = cmfd % nfissxs(2,1,i,j,k)
          nsigf12 = cmfd % nfissxs(1,2,i,j,k)
          nsigf22 = cmfd % nfissxs(2,2,i,j,k)

          ! check for no fission into group 2
          if (.not.(nsigf12 < 1e-6_8 .and. nsigf22 < 1e-6_8)) then
            write(OUTPUT_UNIT,'(A,1PE11.4,1X,1PE11.4)') 'Fission in G=2', &
                  nsigf12,nsigf22
          end if

          ! compute absorption xs
          siga1 = sigt1 - sigs11 - sigs12
          siga2 = sigt2 - sigs22 - sigs21

          ! compute effective downscatter xs
          sigs12_eff = (ONE/keff*nsigf11*flux1 - leak1 - siga1*flux1 &
               - ONE/keff*nsigf21/siga2*leak2 ) / ( flux1*(ONE &
               - ONE/keff*nsigf21/siga2))

          ! redefine flux 2
          flux2 = (sigs12_eff*flux1 - leak2)/siga2
          cmfd % flux(2,i,j,k) = flux2
 
          ! recompute total cross sections (use effective and no upscattering)
          sigt1 = siga1 + sigs11 + sigs12_eff
          sigt2 = siga2 + sigs22

          ! record total xs
          cmfd % totalxs(1,i,j,k) = sigt1
          cmfd % totalxs(2,i,j,k) = sigt2

          ! record effective downscatter xs
          cmfd % scattxs(1,2,i,j,k) = sigs12_eff

          ! zero out upscatter cross section
          cmfd % scattxs(2,1,i,j,k) = ZERO 

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine fix_neutron_balance

!===============================================================================
! FIX_NEUTRON_BALANCE
!===============================================================================

  subroutine compute_effective_downscatter()

    use constants, only: ZERO, CMFD_NOACCEL
    use global,    only: cmfd, cmfd_downscatter

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

    ! check for balance
    if (.not. cmfd_downscatter) return

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! return if not two groups
    if (ng /= 2) return

    ! begin loop around space and energy groups
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          ! check for active mesh
          if (allocated(cmfd%coremap)) then
            if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) cycle
          end if

          ! extract cross sections and flux from object
          flux1 = cmfd % flux(1,i,j,k)
          flux2 = cmfd % flux(2,i,j,k)
          sigt1 = cmfd % totalxs(1,i,j,k)
          sigt2 = cmfd % totalxs(2,i,j,k)
          sigs11 = cmfd % scattxs(1,1,i,j,k)
          sigs21 = cmfd % scattxs(2,1,i,j,k)
          sigs12 = cmfd % scattxs(1,2,i,j,k)
          sigs22 = cmfd % scattxs(2,2,i,j,k)

          ! compute absorption xs
          siga1 = sigt1 - sigs11 - sigs12
          siga2 = sigt2 - sigs22 - sigs21

          ! compute effective downscatter xs
          sigs12_eff = sigs12 - sigs21*flux2/flux1 

          ! recompute total cross sections (use effective and no upscattering)
          sigt1 = siga1 + sigs11 + sigs12_eff
          sigt2 = siga2 + sigs22
    
          ! record total xs
          cmfd % totalxs(1,i,j,k) = sigt1
          cmfd % totalxs(2,i,j,k) = sigt2

          ! record effective downscatter xs
          cmfd % scattxs(1,2,i,j,k) = sigs12_eff

          ! zero out upscatter cross section
          cmfd % scattxs(2,1,i,j,k) = ZERO 

        end do XLOOP

      end do YLOOP
    
    end do ZLOOP

  end subroutine compute_effective_downscatter 

end module cmfd_data
