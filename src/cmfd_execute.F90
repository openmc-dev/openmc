module cmfd_execute

  use cmfd_data,         only: set_up_cmfd
  use cmfd_output,       only: write_cmfd_hdf5
  use global,            only: cmfd,cmfd_only,time_cmfd,master,rank,mpi_err,   &
 &                             current_batch,n_inactive,n_batches,n_procs,     &
 &                             n_procs_cmfd,neut_feedback
  use timing,            only: timer_start,timer_stop,timer_reset


#ifdef PETSC
  use cmfd_power_solver, only: cmfd_power_execute
  use cmfd_slepc_solver, only: cmfd_slepc_execute
  use cmfd_snes_solver,  only: cmfd_snes_execute
#endif
 
  implicit none

#ifdef PETSC
# include <finclude/petsc.h90>
# include <finclude/slepcsys.h>
# include <finclude/slepceps.h>

contains

!==============================================================================
! EXECUTE_CMFD
!==============================================================================

  subroutine execute_cmfd()

    integer :: ierr  ! petsc error code

    ! set how many processors to run cmfd on
    n_procs_cmfd = min(n_procs_cmfd,n_procs) 

    ! initialize mpi communicator
    if(current_batch == n_inactive + 1) call petsc_init_mpi()

    ! filter procs
    if (rank < n_procs_cmfd) then

      ! initialize slepc/petsc (communicates to world)
      if(current_batch == n_inactive + 1) call SlepcInitialize                 &
     &                                       (PETSC_NULL_CHARACTER,ierr)

      ! only run if master process
      if (master) then

        ! begin timer
        call timer_start(time_cmfd)

        ! set up cmfd
        if(.not. cmfd_only) call set_up_cmfd()

      end if

      ! broadcast cmfd object to all procs
      call cmfd_bcast()

      ! execute snes solver
      call cmfd_snes_execute()
!     call cmfd_slepc_execute()

      ! only run if master process
      if (master) call timer_stop(time_cmfd)

      ! finalize slepc
      if (current_batch == n_batches) call SlepcFinalize(ierr)

    end if

    ! sync up procs
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

    ! compute fission source for reweighting
    call calc_fission_source()

    ! perform cmfd re-weighting
    if (neut_feedback) call cmfd_reweight()

    ! write out hdf5 file for cmfd object
    if (master) then
#     ifdef HDF5
        call write_cmfd_hdf5()
#     endif
    end if

  end subroutine execute_cmfd

!==============================================================================
! CMFD_BCAST 
!==============================================================================

  subroutine cmfd_bcast()

    use cmfd_header, only: allocate_cmfd
    use global,      only: cmfd_coremap

    integer :: nx  ! number of mesh cells in x direction
    integer :: ny  ! number of mesh cells in y direction
    integer :: nz  ! number of mesh cells in z direction
    integer :: ng  ! number of energy groups

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! initialize data
    call allocate_cmfd(cmfd)

    ! sync up procs
    call MPI_Barrier(PETSC_COMM_WORLD,mpi_err)

    ! broadcast all data
    call MPI_BCAST(cmfd%flux,ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%totalxs,ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%p1scattxs,ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%scattxs,ng*ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%nfissxs,ng*ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%diffcof,ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%dtilde,6*ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%dhat,6*ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%hxyz,3*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%current,12*ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)

    ! broadcast coremap info
    if (cmfd_coremap) then
      call MPI_BCAST(cmfd%coremap,nx*ny*nz,MPI_INT,0,PETSC_COMM_WORLD,mpi_err)
      call MPI_BCAST(cmfd%mat_dim,1,MPI_INT,0,PETSC_COMM_WORLD,mpi_err)
      if (.not. allocated(cmfd % indexmap)) allocate                           &
     &           (cmfd % indexmap(cmfd % mat_dim,3))
      call MPI_BCAST(cmfd%indexmap,cmfd%mat_dim*3,MPI_INT,0,PETSC_COMM_WORLD,mpi_err)
    end if

  end subroutine cmfd_bcast

!===============================================================================
! PETSC_INIT_MPI
!===============================================================================
 
  subroutine petsc_init_mpi()
 
    integer               :: new_comm   ! new communicator
    integer               :: orig_group ! original MPI group for MPI_COMM_WORLD
    integer               :: new_group  ! new MPI group subset of orig_group
    integer,allocatable   :: ranks(:)   ! ranks to include for petsc
    integer               :: k          ! iteration counter
 
    ! set ranks
    if (.not. allocated(ranks)) allocate(ranks(0:n_procs_cmfd - 1))
    ranks = (/(k,k=0,n_procs_cmfd - 1)/)
 
    ! get the original mpi group
    call MPI_COMM_GROUP(MPI_COMM_WORLD,orig_group,mpi_err)
 
    ! new group init
    call MPI_GROUP_INCL(orig_group,size(ranks),ranks,new_group,mpi_err)
 
    ! create new communicator
    call MPI_COMM_CREATE(MPI_COMM_WORLD,new_group,new_comm,mpi_err)
 
    ! deallocate ranks
    if (allocated(ranks)) deallocate(ranks)
 
    ! set PETSC_COMM_WORLD to this subset
    PETSC_COMM_WORLD = new_comm

  end subroutine petsc_init_mpi

!===============================================================================
! CALC_FISSION_SOURCE calculates the cmfd fission source
!===============================================================================

  subroutine calc_fission_source()

    use global, only: cmfd,cmfd_coremap,entropy_on 
    use tally,  only: tally_reset

    ! local variables
    integer :: nx ! maximum number of cells in x direction
    integer :: ny ! maximum number of cells in y direction
    integer :: nz ! maximum number of cells in z direction
    integer :: ng ! maximum number of energy groups
    integer :: i ! iteration counter for x
    integer :: j ! iteration counter for y
    integer :: k ! iteration counter for z
    integer :: g ! iteration counter for groups
    integer :: idx ! index in vector
    real(8) :: hxyz(3) ! cell dimensions of current ijk cell
    real(8) :: vol     ! volume of cell
    real(8),allocatable :: source(:,:,:,:)  ! tmp source array for entropy

    ! get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! allocate cmfd source if not already allocated and allocate buffer
    if (.not. allocated(cmfd%source)) allocate(cmfd%source(ng,nx,ny,nz))

    ! reset cmfd source to 0
    cmfd%source = 0.0_8

    ! only perform for master
    if (master) then

      ! loop around indices to map to cmfd object
      ZLOOP: do k = 1,nz

        YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

            GROUP: do g = 1,ng

              ! check for core map
              if (cmfd_coremap) then
                if (cmfd%coremap(i,j,k) == 99999) then
                  cycle
                end if
              end if

              ! get dimensions of cell
              hxyz = cmfd%hxyz(:,i,j,k)

              ! calculate volume
              vol = hxyz(1)*hxyz(2)*hxyz(3)

              ! get first index
              idx = get_matrix_idx(1,i,j,k,ng,nx,ny)

              ! compute fission source
              cmfd%source(g,i,j,k) = sum(cmfd%nfissxs(:,g,i,j,k)*cmfd%phi(idx:idx+(ng-1)))*vol 

            end do GROUP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

      ! normalize source such that it sums to 1.0
      cmfd%source = cmfd%source/sum(cmfd%source)

      ! compute entropy
      if (entropy_on) then

        ! allocate tmp array
        if (.not.allocated(source)) allocate(source(ng,nx,ny,nz))

        ! initialize the source
        source = 0.0_8

        ! compute log
        where (cmfd%source > 0.0_8)
          source = cmfd%source*log(cmfd%source)/log(2.0_8)
        end where

        ! sum that source
        cmfd%entropy = -1.0_8*sum(source)

        ! deallocate tmp array
        if (allocated(source)) deallocate(source)

      end if

      ! normalize source so average is 1.0
      cmfd%source = cmfd%source*cmfd%norm

    end if

    ! broadcast full source to all procs 
    call MPI_BCAST(cmfd%source,ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)

  end subroutine calc_fission_source 

!===============================================================================
! CMFD_REWEIGHT
!===============================================================================

  subroutine cmfd_reweight()

    use global,      only: n_particles,meshes,source_bank,work
    use mesh_header, only: StructuredMesh
    use mesh,        only: count_bank_sites,get_mesh_indices
    use search,      only: binary_search

    ! local variables
    integer :: nx ! maximum number of cells in x direction
    integer :: ny ! maximum number of cells in y direction
    integer :: nz ! maximum number of cells in z direction
    integer :: ng ! maximum number of energy groups
    integer :: i ! iteration counter
    integer :: ijk(3) ! spatial bin location
    integer :: e_bin ! energy bin of source particle
    integer :: n_groups ! number of energy groups
    logical :: outside ! any source sites outside mesh
    logical :: in_mesh ! source site is inside mesh
    type(StructuredMesh), pointer :: m ! point to mesh
    real(8), allocatable :: egrid(:)

    ! associate pointer
    m => meshes(1)

    ! get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! allocate arrays in cmfd object (can take out later extend to multigroup)
    if (.not.allocated(cmfd%sourcecounts))                                     &
   &         allocate(cmfd%sourcecounts(ng,nx,ny,nz))
    if (.not.allocated(cmfd%weightfactors))                                    &
   &         allocate(cmfd%weightfactors(ng,nx,ny,nz))

    ! allocate energy grid and reverse cmfd energy grid
    if (.not. allocated(egrid)) allocate(egrid(ng+1))
    egrid = (/(cmfd%egrid(ng-i+2),i = 1,ng+1)/)

    ! zero out weights
    cmfd%weightfactors = 0.0_8

    ! count bank sites in mesh
    call count_bank_sites(m,source_bank,cmfd%sourcecounts,egrid,sites_outside=outside)

    ! have master compute weight factors
    if (master) then
      where(cmfd%source > 0.0_8)
        cmfd%weightfactors = cmfd%source/sum(cmfd%source)*dble(n_particles) /  &
     &                      dble(cmfd%sourcecounts)
      end where
    end if

    ! broadcast weight factors to all procs
    call MPI_BCAST(cmfd%weightfactors,ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,  &
   &               mpi_err)

    ! begin loop over source bank
    do i = 1, int(work,4)

      ! determine spatial bin
      call get_mesh_indices(m,source_bank(i)%xyz,ijk,in_mesh)

      ! determine energy bin
      n_groups = size(cmfd%egrid) - 1
      if (source_bank(i) % E < cmfd%egrid(1)) then
        e_bin = 1
        write(*,*) 'Warning, Source pt below energy grid'
      elseif (source_bank(i) % E > cmfd%egrid(n_groups+1)) then
        e_bin = n_groups
        write(*,*) 'Warning, Source pt above energy grid'
      else
        e_bin = binary_search(cmfd%egrid, n_groups + 1, source_bank(i) % E)
      end if

      ! reverese energy bin (lowest grp is highest energy bin)
      e_bin = n_groups - e_bin + 1

      ! check for outside of mesh
      if (.not. in_mesh) then
        write(*,*) 'FATAL: source site found outside of mesh'
        write(*,*) 'XYZ:',source_bank(i)%xyz
        stop
      end if

      ! reweight particle
      source_bank(i)%wgt = source_bank(i)%wgt *                                &
     &                   cmfd%weightfactors(e_bin,ijk(1),ijk(2),ijk(3))

    end do

    ! deallocate
    if (allocated(egrid)) deallocate(egrid)

  end subroutine cmfd_reweight

!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix 
!===============================================================================
  
  function get_matrix_idx(g,i,j,k,ng,nx,ny) result (matidx)

    use global, only: cmfd,cmfd_coremap

    integer :: matidx         ! the index location in matrix
    integer :: i               ! current x index
    integer :: j               ! current y index
    integer :: k               ! current z index
    integer :: g               ! current group index
    integer :: nx ! maximum number of cells in x direction
    integer :: ny ! maximum number of cells in y direction
    integer :: ng ! maximum number of energy groups

    ! check if coremap is used
    if (cmfd_coremap) then

      ! get idx from core map
      matidx = ng*(cmfd % coremap(i,j,k)) - (ng - g)

    else

      ! compute index
      matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if
  
  end function get_matrix_idx

#endif

end module cmfd_execute
