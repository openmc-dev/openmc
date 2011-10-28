module cmfd_utils

  use global

 implicit none

contains

!===============================================================================
! READ_INPUT reads the CMFD input file and organizes it into a data structure
!===============================================================================

  subroutine read_input()

    use xml_data_cmfd_t

    ! local variables
    integer :: nx                    ! number of volumes in x-direction
    integer :: ny                    ! number of volumes in y-direction
    integer :: nz                    ! number of volumes in z-direction
    integer :: ng                    ! number of energy groups
    integer :: nnx                   ! number of sub meshes in x-direction
    integer :: nny                   ! number of sub meshes in y-direction
    integer :: nnz                   ! number of sub meshes in z-direction
    integer :: i                     ! x iteration counter
    integer :: j                     ! y iteration counter
    integer :: k                     ! z iteration counter
    integer :: g                     ! group iteration counter
    integer :: h                     ! group iteration counter
    integer :: map_idx               ! vector location for core map
    integer :: matid                 ! material id number
    integer :: hg_idx                ! energy h-->g vector location 
    integer :: dim_idx               ! vector location for dimensions
    integer :: ii                    ! inner loop x iteration
    integer :: jj                    ! inner loop y iteration
    integer :: kk                    ! inner loop z iteration
    integer :: ix                    ! shifted x index
    integer :: jy                    ! shifted y index
    integer :: kz                    ! shifted z index
    integer, allocatable :: xgrid(:) ! grid in the x direction
    integer, allocatable :: ygrid(:) ! grid in the y direction
    integer, allocatable :: zgrid(:) ! grid in the z direciton

    ! read xml input file
    call read_xml_file_cmfd_t('cmfd.xml')

    ! get mesh and group indices
    nnx = geometry%nx
    nny = geometry%ny
    nnz = geometry%nz
    ng = geometry%ng

    ! allocate grid varibles
    allocate(xgrid(nnx))
    allocate(ygrid(nny))
    allocate(zgrid(nnz))

    ! get grid variables
    xgrid = geometry%xgrid
    ygrid = geometry%ygrid
    zgrid = geometry%zgrid

    ! get total in each direction
    nx = sum(xgrid)
    ny = sum(ygrid)
    nz = sum(zgrid)

    ! allocate cmfd object
    allocate(cmfd%totalxs(ng,nx,ny,nz))
    allocate(cmfd%scattxs(ng,ng,nx,ny,nz))
    allocate(cmfd%nfissxs(ng,ng,nx,ny,nz))
    allocate(cmfd%diffcof(ng,nx,ny,nz))
    allocate(cmfd%dtilda(6,ng,nx,ny,nz))
    allocate(cmfd%hxyz(nx,ny,nz,3))
    allocate(cmfd%coremap(nx,ny,nz))

    ! record indices in object
    cmfd%indices(1) = nx
    cmfd%indices(2) = ny
    cmfd%indices(3) = nz
    cmfd%indices(4) = ng

    ! set boundary conditions
    cmfd%albedo = geometry%bc

    ! check core map dimensions
    if (size(geometry%mesh,1) /= nnx*nny*nnz) then
    
      ! write out fatal error
      print *,'FATAL ===> core map dimensions not consistent'
      STOP

    end if

    ! read in core map and xs
    ZLOOP:  do k = 1,nnz

      YLOOP: do j = 1,nny

        XLOOP: do i = 1,nnx

          GROUP: do g = 1,ng

            ! get vector idx for core map
            map_idx = get_matrix_idx(1,i,j,k,1,nnx,nny)

            ! extract material identifier
            matid = geometry%mesh(map_idx)

            ! record to core map
            ZZLOOP: do kk = 1,zgrid(k)

              YYLOOP: do jj = 1,ygrid(j)

                XXLOOP: do ii = 1,xgrid(i)

                  ! compute shifted indices
                  ix = sum(xgrid(1:i)) - xgrid(i) + ii
                  jy = sum(ygrid(1:j)) - ygrid(j) + jj
                  kz = sum(zgrid(1:k)) - zgrid(k) + kk

                  ! record in object
                  cmfd%coremap(ix,jy,kz) = matid

                  ! check to see if matid is there
                  if (matid > size(mat,1)) then

                    ! write out fatal error
                    print *, 'Fatal Error ===> MATERIAL ID',matid,' NOT SET!'
                    STOP

                  end if

                  ! check if remxs is set
                  if (associated(mat(matid)%remxs)) then

                    ! set arbitrary totalxs
                    cmfd%totalxs(g,ix,jy,kz) = 0.5

                  else 

                    ! set tot xs since it is given
                    cmfd%totalxs(g,ix,jy,kz) = mat(matid)%totalxs(g)

                  end if

                  ! set diffusion coefficient
                  cmfd%diffcof(g,ix,jy,kz) = mat(matid)%diffcoef(g)

                  ! loop around outgoing energy groups 
                  ELOOP: do h = 1,ng

                    ! get vector h-->g index
                    ! two group --> 1-->1,1-->2,2-->1,2-->2
                    hg_idx = g + ng*(h - 1)

                    ! check if remxs is set and if so adjust within group
                    if (associated(mat(matid)%remxs) .and. g == h) then

                      ! set within group scattering
                      cmfd%scattxs(h,g,ix,jy,kz) = cmfd%totalxs(g,ix,jy,kz) -  &
                     &                           mat(matid)%remxs(g)

                    else

                      ! set regular of just off-scattering
                      cmfd%scattxs(h,g,ix,jy,kz) = mat(matid)%scattxs(hg_idx)

                    end if

                    ! check if chi was entered
                    if (associated(mat(matid)%chi)) then
 
                      ! set nfissxs transfer based on chi and nfissxs vector
                      cmfd%nfissxs(h,g,ix,jy,kz) = mat(matid)%chi(g)*mat(matid)%nfissxs(h)

                    else

                      ! user entered in transfer matrix just record
                      cmfd%nfissxs(h,g,ix,jy,kz) = mat(matid)%nfissxs(hg_idx)

                    end if

                  end do ELOOP

                end do XXLOOP

              end do YYLOOP

            end do ZZLOOP

          end do GROUP 

        end do XLOOP

      end do YLOOP

    end do ZLOOP 

    ! get dimensions
    if (associated(geometry%uniform)) then

      ! record uniform dimensions
      cmfd%hxyz(:,:,:,1) = geometry%uniform(1)
      cmfd%hxyz(:,:,:,2) = geometry%uniform(2)
      cmfd%hxyz(:,:,:,3) = geometry%uniform(3)

    else if (associated(geometry%dx)) then

      ! loop through to get nonuniform dimensions
      ZLOOP2: do k = 1,nnz

        YLOOP2: do j = 1,nny

          XLOOP2: do i = 1,nnx

            ! record to core map
            ZZLOOP2: do kk = 1,zgrid(k)

              YYLOOP2: do jj = 1,ygrid(j)

                XXLOOP2: do ii = 1,xgrid(i)

                  ! compute shifted indices
                  ix = sum(xgrid(1:i)) - xgrid(i) + ii
                  jy = sum(ygrid(1:j)) - ygrid(j) + jj
                  kz = sum(zgrid(1:k)) - zgrid(k) + kk

                  ! record dimension
                  cmfd%hxyz(ix,jy,kz,1) = geometry%dx(i)/xgrid(i)
                  cmfd%hxyz(ix,jy,kz,2) = geometry%dy(j)/ygrid(j)
                  cmfd%hxyz(ix,jy,kz,3) = geometry%dz(k)/zgrid(k)

                end do XXLOOP2

              end do YYLOOP2

            end do ZZLOOP2

          end do XLOOP2

        end do YLOOP2

      end do ZLOOP2

    else

      ! user forgot dimensions
      print *,'Fatal Error ===> Dimensions not entered correctly!'
      STOP

    end if

    ! write out core map and deallocate
    write(200,*) cmfd%coremap
    deallocate(cmfd%coremap)

    ! echo input
!   print *, 'Dimensions:'
!   print *,cmfd%indices
!   print *, 'CORE MAP:'
!   print *,cmfd%coremap
!   print *, 'TOTAL XS:'
!   print *,minval(cmfd%totalxs),maxval(cmfd%totalxs)
!   print *, 'SCATTERING XS:'
!   print *,minval(cmfd%scattxs),maxval(cmfd%scattxs)
!   print *, 'Nu-FISSION XS:'
!   print *,minval(cmfd%nfissxs),maxval(cmfd%nfissxs)
!   print *, 'DIFFUSION COEFFICIENT:'
!   print *,minval(cmfd%diffcof),maxval(cmfd%diffcof)
!   print *, 'BOUNDARY CONDITIONS:'
!   print *,cmfd%albedo
!   print *, 'CORE CELL DIMENSIONS X:'
!   print *,cmfd%hxyz(:,:,:,1)
!   print *, 'CORE CELL DIMENSIONS Y:'
!   print *,cmfd%hxyz(:,:,:,2)
!   print *, 'CORE CELL DIMENSIONS Z:'
!   print *,cmfd%hxyz(:,:,:,3)

  end subroutine read_input

!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix 
!===============================================================================

  function get_matrix_idx(g,i,j,k,ng,nx,ny)

    ! arguments
    integer :: get_matrix_idx  ! the index location in matrix
    integer :: i               ! current x index
    integer :: j               ! current y index
    integer :: k               ! current z index
    integer :: g               ! current group index
    integer :: ng              ! max energy groups
    integer :: nx              ! maximum cells in x direction
    integer :: ny              ! maximum cells in y direction

    ! local variables
    integer :: nidx            ! index in matrix

    ! compute index
    nidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    ! record value to function
    get_matrix_idx = nidx

  end function get_matrix_idx

end module cmfd_utils
