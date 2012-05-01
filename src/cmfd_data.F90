module cmfd_data

  implicit none
  private
  public :: set_up_cmfd

contains

!===============================================================================
! SET_UP_CMFD
!===============================================================================

  subroutine set_up_cmfd()

    use global,       only: cmfd,cmfd_coremap,current_batch,n_inactive,        &
   &                        cmfd_balance,cmfd_downscatter
    use cmfd_header,  only: allocate_cmfd
    use cmfd_output,  only: neutron_balance

    ! initialize data
    call allocate_cmfd(cmfd)

    ! check for core map
    if ((cmfd_coremap) .and. (current_batch == n_inactive+1)) then
      call set_coremap()
    end if

    ! calculate all cross sections based on reaction rates from last batch
    call compute_xs()

    ! overwrite xs
!   call fix_1_grp()
!   call fix_2_grp()

    ! fix balance of neutrons
    if (cmfd_balance) call fix_neutron_balance()

    ! compute effective downscatter
    if (cmfd_downscatter) call compute_effective_downscatter()

    ! compute neutron balance 
    call neutron_balance()

    ! compute dtilde terms
    call compute_dtilde()

    ! compute dhat terms
    call compute_dhat()

  end subroutine set_up_cmfd

!===============================================================================
! COMPUTE_XS takes tallies and computes macroscopic cross sections
!===============================================================================

  subroutine compute_xs()

    use global
    use mesh,          only: mesh_indices_to_bin
    use mesh_header,   only: StructuredMesh
    use tally_header,  only: TallyObject, TallyScore
 
    integer :: nx                ! number of mesh cells in x direction
    integer :: ny                ! number of mesh cells in y direction
    integer :: nz                ! number of mesh cells in z direction
    integer :: ng                ! number of energy groups
    integer :: i                 ! iteration counter for x
    integer :: j                 ! iteration counter for y
    integer :: k                 ! iteration counter for z
    integer :: g                 ! iteration counter for g
    integer :: h                 ! iteration counter for outgoing groups
    integer :: ital              ! tally object index
    integer :: ijk(3)            ! indices for mesh cell
    integer :: score_index       ! index to pull from tally object
    integer :: bins(N_FILTER_TYPES) ! bins for filters

    real(8) :: flux   ! temp variable for flux

    type(TallyObject), pointer :: t    ! pointer for tally object
    type(StructuredMesh), pointer :: m ! pointer for mesh object

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! set flux object and source distribution to all zeros
    cmfd % flux = 0.0_8
    cmfd % openmc_src = 0.0_8

    ! associate tallies and mesh
    t => tallies(n_user_tallies + 1)
    m => meshes(t % mesh)

    ! set mesh widths
    cmfd % hxyz(1,:,:,:) = m % width(1) ! set x width
    cmfd % hxyz(2,:,:,:) = m % width(2) ! set y width
    cmfd % hxyz(3,:,:,:) = m % width(3) ! set z width

    ! begin loop around tallies
    TAL: do ital = n_user_tallies + 1,n_tallies 

      ! associate tallies and mesh
      t => tallies(ital)
      m => meshes(t % mesh)

      ! begin loop around space
      ZLOOP: do k = 1,nz

        YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx
 
            ! check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == 99999) then
                cycle
              end if
            end if

            ! loop around energy groups
            OUTGROUP: do h = 1,ng

              ! start tally 1
              TALLY: if (ital == n_user_tallies + 1) then

                ! reset all bins to 1
                bins = 1

                ! set ijk as mesh indices
                ijk = (/ i, j, k /)

                ! get bin number for mesh indices
                bins(FILTER_MESH) = mesh_indices_to_bin(m,ijk)

                ! apply energy in filter
                bins(FILTER_ENERGYIN) = ng - h + 1

                ! calculate score index from bins
                score_index = sum((bins - 1) * t%stride) + 1

                ! get flux
                flux = t % scores(1,score_index) % sum
                cmfd % flux(h,i,j,k) = flux

                ! detect zero flux
                if ((flux - 0.0D0) < 1.0D-10) then
                  write(*,*) 'Fatal: detected zero flux without coremap'
                  stop
                end if

                ! get total rr and convert to total xs
                cmfd % totalxs(h,i,j,k) = t % scores(2,score_index) % sum / flux

                ! get p1 scatter rr and convert to p1 scatter xs
                cmfd % p1scattxs(h,i,j,k) = t % scores(3,score_index) % sum / flux

                ! extract diffusion coefficient tally
                cmfd % diffusion(h,i,j,k) = t % scores(4,score_index) % sum / flux

                ! calculate diffusion coefficient
!               cmfd % diffcof(h,i,j,k) = 1.0_8/(3.0_8*cmfd%totalxs(h,i,j,k))
!               cmfd % diffcof(h,i,j,k) = 1.0_8/(3.0_8*(cmfd % totalxs(h,i,j,k) -&
!              &                                cmfd % p1scattxs(h,i,j,k)))
                cmfd % diffcof(h,i,j,k) = cmfd % diffusion(h,i,j,k)

              else if (ital == n_user_tallies + 2) then

                ! begin loop to get energy out tallies
                INGROUP: do g = 1,ng

                  ! reset all bins to 1
                  bins = 1

                  ! set ijk as mesh indices
                  ijk = (/ i, j, k /)

                  ! get bin number for mesh indices
                  bins(FILTER_MESH) = mesh_indices_to_bin(m,ijk)

                  ! apply energy in filter
                  bins(FILTER_ENERGYIN) = ng - h + 1

                  ! set energy out bin
                  bins(FILTER_ENERGYOUT) = ng - g + 1

                  ! calculate score index from bins
                  score_index = sum((bins - 1) * t%stride) + 1

                  ! get scattering
                  cmfd % scattxs(h,g,i,j,k) = t % scores(1,score_index) % sum /&
                 &                            cmfd % flux(h,i,j,k)

                  ! get nu-fission
                  cmfd % nfissxs(h,g,i,j,k) = t % scores(2,score_index) % sum /&
                 &                            cmfd % flux(h,i,j,k)

                  ! bank source
                  cmfd % openmc_src(g,i,j,k) = cmfd % openmc_src(g,i,j,k) +    &
                 &                             t % scores(2,score_index) % sum 

                end do INGROUP

              else

                ! initialize and filter for energy
                bins = 1
                bins(SURF_FILTER_ENERGYIN) = ng - h + 1

                ! left surface
                bins(1:3) = (/ i-1, j, k /) + 1
                bins(SURF_FILTER_SURFACE) = IN_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1 ! outgoing
                cmfd % current(1,h,i,j,k) = t % scores(1,score_index) % sum
                bins(SURF_FILTER_SURFACE) = OUT_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
                cmfd % current(2,h,i,j,k) = t % scores(1,score_index) % sum

                ! right surface
                bins(1:3) = (/ i, j, k /) + 1
                bins(SURF_FILTER_SURFACE) = IN_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
                cmfd % current(3,h,i,j,k) = t % scores(1,score_index) % sum
                bins(SURF_FILTER_SURFACE) = OUT_RIGHT
                score_index = sum((bins - 1) * t % stride) + 1 ! outgoing 
                cmfd % current(4,h,i,j,k) = t % scores(1,score_index) % sum

                ! back surface
                bins(1:3) = (/ i, j-1, k /) + 1
                bins(SURF_FILTER_SURFACE) = IN_FRONT
                score_index = sum((bins - 1) * t % stride) + 1 ! outgoing
                cmfd % current(5,h,i,j,k) = t % scores(1,score_index) % sum
                bins(SURF_FILTER_SURFACE) = OUT_FRONT
                score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
                cmfd % current(6,h,i,j,k) = t % scores(1,score_index) % sum

                ! front surface
                bins(1:3) = (/ i, j, k /) + 1
                bins(SURF_FILTER_SURFACE) = IN_FRONT
                score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
                cmfd % current(7,h,i,j,k) = t % scores(1,score_index) % sum
                bins(SURF_FILTER_SURFACE) = OUT_FRONT
                score_index = sum((bins - 1) * t % stride) + 1 ! outgoing 
                cmfd % current(8,h,i,j,k) = t % scores(1,score_index) % sum

                ! bottom surface
                bins(1:3) = (/ i, j, k-1 /) + 1
                bins(SURF_FILTER_SURFACE) = IN_TOP
                score_index = sum((bins - 1) * t % stride) + 1 ! outgoing
                cmfd % current(9,h,i,j,k) = t % scores(1,score_index) % sum
                bins(SURF_FILTER_SURFACE) = OUT_TOP
                score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
                cmfd % current(10,h,i,j,k) = t % scores(1,score_index) % sum

                ! top surface
                bins(1:3) = (/ i, j, k /) + 1
                bins(SURF_FILTER_SURFACE) = IN_TOP
                score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
                cmfd % current(11,h,i,j,k) = t % scores(1,score_index) % sum
                bins(SURF_FILTER_SURFACE) = OUT_TOP
                score_index = sum((bins - 1) * t % stride) + 1 ! outgoing 
                cmfd % current(12,h,i,j,k) = t % scores(1,score_index) % sum

              end if TALLY

            end do OUTGROUP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

    end do TAL

    ! normalize openmc source distribution
    cmfd % openmc_src = cmfd % openmc_src/sum(cmfd % openmc_src)*cmfd%norm

  end subroutine compute_xs

!===============================================================================
! COMPUTE_DTILDE computes the diffusion coupling coefficient
!===============================================================================

  subroutine compute_dtilde()

    use global, only: cmfd,cmfd_coremap

    ! local variables
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
    ZLOOP:  do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          GROUP: do g = 1,ng

            ! check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == 99999) then
                cycle
              end if
            end if

            ! get cell data
            cell_dc = cmfd%diffcof(g,i,j,k)
            cell_hxyz = cmfd%hxyz(:,i,j,k)

            ! setup of vector to identify boundary conditions
            bound = (/i,i,j,j,k,k/)

            ! begin loop around sides of cell for leakage
            LEAK: do l = 1,6

              ! define xyz and +/- indices
              xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
              dir_idx = 2 - mod(l,2) ! -=1, +=2
              shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1

              ! check if at a boundary
              if (bound(l) == nxyz(xyz_idx,dir_idx)) then

                ! compute dtilde
                dtilde = (2*cell_dc*(1-albedo(l)))/(4*cell_dc*(1+albedo(l)) +  &
               &         (1-albedo(l))*cell_hxyz(xyz_idx))

              else  ! not a boundary

                ! compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! get neigbor cell data
                neig_dc = cmfd%diffcof(g,neig_idx(1),neig_idx(2),neig_idx(3))
                neig_hxyz = cmfd%hxyz(:,neig_idx(1),neig_idx(2),neig_idx(3))

                ! check for fuel-reflector interface
                if (cmfd_coremap) then

                  if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) ==   &
                 &    99999 .and. cmfd % coremap(i,j,k) /= 99999) then

                    ! get albedo
                    ref_albedo = get_reflector_albedo(l,g,i,j,k)

                    ! compute dtilde
                    dtilde = (2*cell_dc*(1-ref_albedo))/(4*cell_dc*(1+         &
                 &         ref_albedo)+(1-ref_albedo)*cell_hxyz(xyz_idx))
      !             dtilde = 0.0_8
                  else ! not next to a reflector or no core map

                    ! compute dtilde
                    dtilde = (2*cell_dc*neig_dc)/(neig_hxyz(xyz_idx)*cell_dc + &
                   &         cell_hxyz(xyz_idx)*neig_dc)

                  end if

                else ! no core map

                  ! compute dtilde
                  dtilde = (2*cell_dc*neig_dc)/(neig_hxyz(xyz_idx)*cell_dc +   &
                 &         cell_hxyz(xyz_idx)*neig_dc)

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

    use global, only:cmfd,cmfd_coremap,dhat_reset 

    ! local variables
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
              if (cmfd%coremap(i,j,k) == 99999) then
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
              net_current = (current(2*l) - current(2*l-1)) /                  &
             &               product(cmfd%hxyz(:,i,j,k)) *                     &
             &               cmfd%hxyz(xyz_idx,i,j,k)

              ! check if at a boundary
              if (bound(l) == nxyz(xyz_idx,dir_idx)) then

                ! compute dhat
                dhat = (net_current - shift_idx*cell_dtilde(l)*cell_flux) /    &
               &        cell_flux

              else  ! not a boundary

                ! compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! get neigbor flux 
                neig_flux = cmfd%flux(g,neig_idx(1),neig_idx(2),neig_idx(3)) / &
                     product(cmfd%hxyz(:,neig_idx(1),neig_idx(2),neig_idx(3)))

                ! check for fuel-reflector interface
                if (cmfd_coremap) then

                  if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) ==   &
                 &    99999 .and. cmfd % coremap(i,j,k) /= 99999) then

                    ! compute dhat
                    dhat = (net_current - shift_idx*cell_dtilde(l)*cell_flux) /&
                   &        cell_flux

                  else ! not a fuel-reflector interface

                    ! compute dhat 
                    dhat = (net_current + shift_idx*cell_dtilde(l)*            &
                   &       (neig_flux - cell_flux))/(neig_flux + cell_flux)

                  end if

                else ! not for fuel-reflector case

                  ! compute dhat 
                  dhat = (net_current + shift_idx*cell_dtilde(l)*              &
                 &       (neig_flux - cell_flux))/(neig_flux + cell_flux)

                end if

              end if

              ! record dtilde in cmfd object
              cmfd%dhat(l,g,i,j,k) = dhat

              ! check for dhat reset
              if (dhat_reset) cmfd%dhat(l,g,i,j,k) = 0.0_8

            end do LEAK

          end do GROUP

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine compute_dhat

!===============================================================================
! SET_COREMAP is a routine that sets the core mapping information
!===============================================================================

  subroutine set_coremap()

    use global, only: cmfd

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
    if (.not. allocated(cmfd % indexmap)) allocate                             &
   &                                      (cmfd % indexmap(cmfd % mat_dim,3))

    ! begin loops over spatial indices
    ZLOOP: do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          ! check for reflector
          if (cmfd % coremap(i,j,k) == 1) then

            ! reset value to 99999
            cmfd % coremap(i,j,k) = 99999

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
! GET_REFLECTOR_ALBEDO is a function that calculates the albedo to the reflector
!===============================================================================

  function get_reflector_albedo(l,g,i,j,k)

    use global, only: cmfd

    ! function variable
    real(8) :: get_reflector_albedo ! reflector albedo

    ! local variable
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
    albedo = (current(2*l-1)/current(2*l))**(shift_idx)

    ! assign to function variable
    get_reflector_albedo = albedo

  end function get_reflector_albedo

!===============================================================================
! FIX_1_GRP modifies xs for benchmark with stand alone (homogeneous)
!===============================================================================

  subroutine fix_1_grp()

    use global, only: cmfd,dhat_reset

    ! overwrite cross sections
    cmfd % totalxs = 1.0_8
    cmfd % scattxs = 0.5_8
    cmfd % nfissxs = 0.48_8
    cmfd % diffcof = 0.3_8
    cmfd % hxyz(1,:,:,:) = 2.0_8
    cmfd % hxyz(2,:,:,:) = 0.1_8
    cmfd % hxyz(3,:,:,:) = 0.1_8

    ! set dhat reset to true
    dhat_reset = .true.

  end subroutine fix_1_grp

!===============================================================================
! FIX_2_GRP modifies xs for benchmark with stand alone (homogeneous)
!===============================================================================

  subroutine fix_2_grp()

    use global, only: cmfd,dhat_reset

    ! overwrite cross sections
    cmfd % totalxs(1,:,:,:) = 0.6451_8
    cmfd % totalxs(2,:,:,:) = 1.4050_8
    cmfd % scattxs(1,1,:,:,:) = 0.6134_8
    cmfd % scattxs(1,2,:,:,:) = 0.0219_8
    cmfd % scattxs(2,1,:,:,:) = 0.0_8
    cmfd % scattxs(2,2,:,:,:) = 1.2929_8
    cmfd % nfissxs(1,1,:,:,:) = 0.0064_8
    cmfd % nfissxs(1,2,:,:,:) = 0.0_8
    cmfd % nfissxs(2,1,:,:,:) = 0.1923_8
    cmfd % nfissxs(2,2,:,:,:) = 0.0_8
    cmfd % diffcof(1,:,:,:) = 1.3050_8
    cmfd % diffcof(2,:,:,:) = 0.5112_8
    cmfd % hxyz(1,:,:,:) = 2.0_8
    cmfd % hxyz(2,:,:,:) = 0.1_8
    cmfd % hxyz(3,:,:,:) = 0.1_8

    ! set dhat reset to true
    dhat_reset = .true.

  end subroutine fix_2_grp

!===============================================================================
! FIX_NEUTRON_BALANCE
!===============================================================================

  subroutine fix_neutron_balance()

    use constants, only: ONE,ZERO
    use global, only: cmfd,cmfd_balance,keff 

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

    ! check for balance
    if (.not. cmfd_balance) return

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! return if not two groups
    if (ng /= 2) return

    ! begin loop around space and energy groups
    ZLOOP: do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          ! check for active mesh
          if (allocated(cmfd%coremap)) then
            if (cmfd%coremap(i,j,k) == 99999) cycle 
          end if

          ! compute leakage in groups 1 and 2
          leak1 = ZERO 
          leak2 = ZERO
          LEAK: do l = 1,3

            leak1 = leak1 + ((cmfd % current(4*l,1,i,j,k) - &
             & cmfd % current(4*l-1,1,i,j,k))) - &
             & ((cmfd % current(4*l-2,1,i,j,k) - &
             & cmfd % current(4*l-3,1,i,j,k)))

            leak2 = leak2 + ((cmfd % current(4*l,2,i,j,k) - &
             & cmfd % current(4*l-1,2,i,j,k))) - &
             & ((cmfd % current(4*l-2,2,i,j,k) - &
             & cmfd % current(4*l-3,2,i,j,k)))


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
          if (.not.(nsigf12 < 1e-8_8 .and. nsigf22 < 1e-8_8)) then
            write(*,*) 'Downscatter rebalance not valid for this group struct'
            stop
          end if

          ! compute absorption xs
          siga1 = sigt1 - sigs11 - sigs12
          siga2 = sigt2 - sigs22 - sigs21

          ! compute effective downscatter xs
          sigs12_eff = ( ONE/keff*nsigf11*flux1 - leak1 - siga1*flux1          &
         &               - ONE/keff*nsigf21/siga2*leak2 ) / ( flux1*(ONE       &
         &               - ONE/keff*nsigf21/siga2) )

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

!         print *,'Leakage g=1:',leak1
!         print *,'Leakage g=2:',leak2
!         print *,'Flux g=1:',flux1
!         print *,'Flux g=2:',flux2
!         print *,'Total g=1:',sigt1
!         print *,'Total g=2:',sigt2
!         print *,'Scattering 1-->1:',sigs11
!         print *,'Scattering 2-->1:',sigs21
!         print *,'Scattering 1-->2:',sigs12
!         print *,'Scattering 2-->2:',sigs22
!         print *,'Fission 1-->1:',nsigf11
!         print *,'Fission 2-->1:',nsigf21
!         print *,'Fission 1-->2:',nsigf12
!         print *,'Fission 2-->2:',nsigf22
!         print *,'keff:',keff
!         print *,'Effective Downscatter:',sigs12_eff
!
!         print *,'Group 1 balance:',leak1+sigt1*flux1-sigs11*flux1-ONE/keff*(nsigf11*flux1+nsigf21*flux2)
!         print *,'Group 2 balance:',leak2+sigt2*flux2-sigs12_eff*flux1-sigs22*flux2
! stop
        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine fix_neutron_balance

!===============================================================================
! FIX_NEUTRON_BALANCE
!===============================================================================

  subroutine compute_effective_downscatter()

    use constants, only: ZERO
    use global, only: cmfd,cmfd_downscatter,keff

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
    ZLOOP: do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          ! check for active mesh
          if (allocated(cmfd%coremap)) then
            if (cmfd%coremap(i,j,k) == 99999) cycle
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

!==============================================================================
! EXTRACT_ACCUM_FSRC
!==============================================================================
!
!  subroutine extract_accum_fsrc(i,j,k)
!
!    use global
!    use mesh,          only: mesh_indices_to_bin
!    use mesh_header,   only: StructuredMesh
!    use tally_header,  only: TallyObject, TallyScore
!
!    integer :: nx                ! number of mesh cells in x direction
!    integer :: ny                ! number of mesh cells in y direction
!    integer :: nz                ! number of mesh cells in z direction
!    integer :: ng                ! number of energy groups
!    integer :: i                 ! iteration counter for x
!    integer :: j                 ! iteration counter for y
!    integer :: k                 ! iteration counter for z
!    integer :: g                 ! iteration counter for g
!    integer :: h                 ! iteration counter for outgoing groups
!    integer :: ital              ! tally object index
!    integer :: ijk(3)            ! indices for mesh cell
!    integer :: score_index       ! index to pull from tally object
!    integer :: bins(N_FILTER_TYPES) ! bins for filters
!
!    real(8) :: flux   ! temp variable for flux
!
!    type(TallyObject), pointer :: t    ! pointer for tally object
!    type(StructuredMesh), pointer :: m ! pointer for mesh object
!
!    ! associate tallies and mesh
!    t => tallies(1)
!    m => meshes(t % mesh)
!
!    ! reset all bins to 1
!    bins = 1
!
!    ! set ijk as mesh indices
!    ijk = (/ i, j, k /)
!
!    ! get bin number for mesh indices
!    bins(FILTER_MESH) = mesh_indices_to_bin(m,ijk)
!
!    ! calculate score index from bins
!    score_index = sum((bins - 1) * t%stride) + 1
!
!    ! bank source
!    cmfd % openmc_accum_src(i,j,k) = t % scores(1,score_index) % sum
!
!  end subroutine extract_accum_fsrc

end module cmfd_data
