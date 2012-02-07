module cmfd_data

  implicit none
  private
  public :: set_up_cmfd

contains

!==============================================================================
! SET_UP_CMFD
!===============================================================================

  subroutine set_up_cmfd()

    use global,       only: cmfd,cmfd_coremap,current_cycle,n_inactive,time_cmfd
    use cmfd_header,  only: allocate_cmfd
    use cmfd_output,  only: neutron_balance,write_cmfd_hdf5
    use timing

    ! initialize data
    call allocate_cmfd(cmfd)

    ! check for core map
    if ((cmfd_coremap) .and. (current_cycle == n_inactive+1)) then
      call set_coremap()
    end if

    ! calculate all cross sections based on reaction rates from last batch
    call compute_xs()

    ! write out the neutron balance file
    call neutron_balance()

    ! compute dtilde terms
    call compute_dtilde()

    ! compute dhat terms
    call compute_dhat()

#ifdef HDF5
    ! write out hdf5 file for cmfd object
    call write_cmfd_hdf5()
#endif

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

    ! set flux object to all zeros
    cmfd % flux = 0.0_8

    ! begin loop around space and energy groups
    ZLOOP: do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          OUTGROUP: do h = 1,ng

            ! check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == 99999) then
                cycle
              end if
            end if

            ! begin with first tally 
            t => tallies(1)
            m => meshes(t % mesh)

            ! set mesh widths
            cmfd % hxyz(1,:,:,:) = m % width(1) ! set x width
            cmfd % hxyz(2,:,:,:) = m % width(2) ! set y width
            cmfd % hxyz(3,:,:,:) = m % width(3) ! set z width

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
            flux = t % scores(score_index,1) % val
            cmfd % flux(h,i,j,k) = flux

            ! detect zero flux
            if ((flux - 0.0D0) < 1.0D-10) then
              if (.not. cmfd_coremap) then
                write(*,*) 'Fatal: detected zero flux without coremap'
                stop
              else
!               write(*,*) 'Warning: detected zero flux at:',i,j,k
                flux = 99999.0D0
                cycle
!               write(*,*) 'Core map location is:',cmfd%coremap(i,j,k)
                if (.not. cmfd%coremap(i,j,k) == 99999) then
!                 write(*,*) 'Fatal: need to check core map with zero flux'
                  stop
                end if
              end if
            end if

            ! get total rr and convert to total xs
            cmfd % totalxs(h,i,j,k) = t % scores(score_index,2) % val / flux

            ! get p1 scatter rr and convert to p1 scatter xs
            cmfd % p1scattxs(h,i,j,k) = t % scores(score_index,3) % val / flux

            ! calculate diffusion coefficient
            cmfd % diffcof(h,i,j,k) = 1.0_8/(3.0_8*(cmfd % totalxs(h,i,j,k) -  &
           &                                cmfd % p1scattxs(h,i,j,k)))
            cmfd % diffcof(h,i,j,k) = 1.0_8/(3.0_8*(cmfd % totalxs(h,i,j,k)))

            ! begin loop to get energy out tallies
            INGROUP: do g = 1,ng

              ! associate tally pointer to energy out tally object
              t => tallies(2)

              ! set energy out bin
              bins(FILTER_ENERGYOUT) = ng - g + 1

              ! calculate score index from bins
              score_index = sum((bins - 1) * t%stride) + 1

              ! get scattering
              cmfd % scattxs(h,g,i,j,k) = t % scores(score_index,1) % val / flux

              ! get nu-fission
              cmfd % nfissxs(h,g,i,j,k) = t % scores(score_index,2) % val / flux

            end do INGROUP

            ! extract surface currents 
            t => tallies(3)

            ! initialize and filter for energy
            bins = 1
            bins(SURF_FILTER_ENERGYIN) = ng - h + 1

            ! left surface
            bins(1:3) = (/ i-1, j, k /) + 1
            bins(SURF_FILTER_SURFACE) = IN_RIGHT
            score_index = sum((bins - 1) * t % stride) + 1 ! outgoing
            cmfd % current(1,h,i,j,k) = t % scores(score_index,1) % val
            bins(SURF_FILTER_SURFACE) = OUT_RIGHT
            score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
            cmfd % current(2,h,i,j,k) = t % scores(score_index,1) % val

            ! right surface
            bins(1:3) = (/ i, j, k /) + 1
            bins(SURF_FILTER_SURFACE) = IN_RIGHT
            score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
            cmfd % current(3,h,i,j,k) = t % scores(score_index,1) % val
            bins(SURF_FILTER_SURFACE) = OUT_RIGHT
            score_index = sum((bins - 1) * t % stride) + 1 ! outgoing 
            cmfd % current(4,h,i,j,k) = t % scores(score_index,1) % val

            ! back surface
            bins(1:3) = (/ i, j-1, k /) + 1
            bins(SURF_FILTER_SURFACE) = IN_FRONT
            score_index = sum((bins - 1) * t % stride) + 1 ! outgoing
            cmfd % current(5,h,i,j,k) = t % scores(score_index,1) % val
            bins(SURF_FILTER_SURFACE) = OUT_FRONT
            score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
            cmfd % current(6,h,i,j,k) = t % scores(score_index,1) % val

            ! front surface
            bins(1:3) = (/ i, j, k /) + 1
            bins(SURF_FILTER_SURFACE) = IN_FRONT
            score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
            cmfd % current(7,h,i,j,k) = t % scores(score_index,1) % val
            bins(SURF_FILTER_SURFACE) = OUT_FRONT
            score_index = sum((bins - 1) * t % stride) + 1 ! outgoing 
            cmfd % current(8,h,i,j,k) = t % scores(score_index,1) % val

            ! bottom surface
            bins(1:3) = (/ i, j, k-1 /) + 1
            bins(SURF_FILTER_SURFACE) = IN_TOP
            score_index = sum((bins - 1) * t % stride) + 1 ! outgoing
            cmfd % current(9,h,i,j,k) = t % scores(score_index,1) % val
            bins(SURF_FILTER_SURFACE) = OUT_TOP
            score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
            cmfd % current(10,h,i,j,k) = t % scores(score_index,1) % val

            ! top surface
            bins(1:3) = (/ i, j, k /) + 1
            bins(SURF_FILTER_SURFACE) = IN_TOP
            score_index = sum((bins - 1) * t % stride) + 1 ! incoming 
            cmfd % current(11,h,i,j,k) = t % scores(score_index,1) % val
            bins(SURF_FILTER_SURFACE) = OUT_TOP
            score_index = sum((bins - 1) * t % stride) + 1 ! outgoing 
            cmfd % current(12,h,i,j,k) = t % scores(score_index,1) % val

          end do OUTGROUP

        end do XLOOP

      end do YLOOP

    end do ZLOOP

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
    real(8) :: cell_totxs         ! total cross section of current ijk cell
    real(8) :: cell_dc            ! diffusion coef of current cell
    real(8) :: cell_hxyz(3)       ! cell dimensions of current ijk cell
    real(8) :: neig_totxs         ! total xs of neighbor cell
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

    use global, only:cmfd,cmfd_coremap 

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
    integer :: ng                ! number of energy groups
    integer :: i                 ! iteration counter for x
    integer :: j                 ! iteration counter for y
    integer :: k                 ! iteration counter for z

    ! extract spatial indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)

    ! count how many fuel assemblies exist
    cmfd % mat_dim = sum(cmfd % coremap - 1)

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

end module cmfd_data 
