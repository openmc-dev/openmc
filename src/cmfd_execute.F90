module cmfd_execute 

!==============================================================================
! CMFD_EXECUTE -- This module is the highest level cmfd module that controls the
! cross section generation, diffusion calculation, and source re-weighting
!==============================================================================

  implicit none
  private
  public :: execute_cmfd, cmfd_init_batch

# ifdef PETSC
#   include <finclude/petsc.h90>
# endif

contains

!==============================================================================
! EXECUTE_CMFD
!==============================================================================

  subroutine execute_cmfd()

# ifdef PETSC

    use cmfd_data,              only: set_up_cmfd
    use cmfd_message_passing,   only: petsc_init_mpi, cmfd_bcast
    use cmfd_power_solver,      only: cmfd_power_execute
    use cmfd_snes_solver,       only: cmfd_snes_execute
    use error,                  only: warning, fatal_error 
    use global,                 only: n_procs_cmfd, cmfd,                       &
                                      cmfd_solver_type, time_cmfd,              &
                                      cmfd_run_adjoint, cmfd_write_hdf5,        &
                                      cmfd_feedback,cmfd_hold_weights,          &
                                      cmfd_inact_flush, cmfd_keff_tol,          &
                                      cmfd_act_flush, current_batch, keff,      &
                                      n_batches, message, master, mpi_err, rank

    logical :: leave_cmfd

    ! set leave cmfd to false
    leave_cmfd = .false.

    ! stop cmfd timer
    if (master) then
      call time_cmfd % start()
    end if

    ! filter processors (lowest PETSc group)
    if (rank < n_procs_cmfd) then

      ! set up cmfd data (master only)
      if (master) call set_up_cmfd()

      ! broadcast cmfd to all petsc procs
      call cmfd_bcast()

      ! process solver options
      call process_cmfd_options()

    end if

    ! check to hold weights
    if (cmfd_hold_weights) then
      message = 'Not Modifying Weights - Albedo estimate not good, increase batch size.'
      call warning() 
      cmfd_hold_weights = .false.
      if (cmfd_feedback) call cmfd_reweight(.false.)
      leave_cmfd = .true. 
    end if
    call MPI_BCAST(leave_cmfd, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_err)
    if (leave_cmfd) return

    ! filter processors (lowest PETSc group)
    if (rank < n_procs_cmfd) then

      ! call solver
      if (trim(cmfd_solver_type) == 'power') then
        call cmfd_power_execute()
      elseif (trim(cmfd_solver_type) == 'jfnk') then
        call cmfd_snes_execute()
      else
        message = 'solver type became invalid after input processing'
        call fatal_error() 
      end if

      ! perform any last batch tasks 
      if (current_batch == n_batches) then

        ! check for adjoint run
        if (cmfd_run_adjoint) then
          if (trim(cmfd_solver_type) == 'power') then
            call cmfd_power_execute(adjoint = .true.)
          elseif (trim(cmfd_solver_type) == 'jfnk') then
            call cmfd_snes_execute(adjoint = .true.)
          end if
        end if

      end if

    end if

    ! check to hold weights
    if ((abs(cmfd%keff-keff)/keff > cmfd_keff_tol)) then
      if (current_batch >= cmfd_inact_flush(1) .or. &
           current_batch >= cmfd_act_flush - 1 ) then
        message = 'Not Modifying Weights - keff %diff > 0.005, up batch size'
        call warning() 
        if (cmfd_feedback) call cmfd_reweight(.false.)
        leave_cmfd = .true.
      end if
    end if
    call MPI_BCAST(leave_cmfd, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_err)
    if (leave_cmfd) return

    ! calculate fission source
    call calc_fission_source()

    ! calculate weight factors
    if (cmfd_feedback) call cmfd_reweight(.true.)

    ! stop cmfd timer
    if (master) then
      call time_cmfd % stop()
    end if

    ! wait here for all procs
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

# endif

  end subroutine execute_cmfd

!==============================================================================
! CMFD_INIT_BATCH
!==============================================================================

  subroutine cmfd_init_batch()

    use global,            only: cmfd_begin, cmfd_on, cmfd_tally_on,         &
                                 cmfd_inact_flush, cmfd_act_flush, cmfd_run, &
                                 current_batch, cmfd_hold_weights

    ! check to activate CMFD diffusion and possible feedback
    ! this guarantees that when cmfd begins at least one batch of tallies are
    ! accumulated
    if (cmfd_run .and. cmfd_begin == current_batch) then
      cmfd_on = .true.
      cmfd_tally_on = .true.
    end if

    ! check to flush cmfd tallies for active batches, no more inactive flush
    if (cmfd_run .and. cmfd_act_flush == current_batch) then
      call cmfd_tally_reset()
      cmfd_tally_on = .true.
      cmfd_inact_flush(2) = -1
    end if

    ! check to flush cmfd tallies during inactive batches (>= on number of
    ! flushes important as the code will flush on the first batch which we
    ! dont want to count)
!    if (cmfd_run .and. current_batch < n_inactive .and. mod(current_batch-1,cmfd_inact_flush(1))   &
!       == 0 .and. cmfd_inact_flush(2) >= 0) then
    if (cmfd_run .and. mod(current_batch,cmfd_inact_flush(1))   &
       == 0 .and. cmfd_inact_flush(2) > 0 .and. cmfd_begin < current_batch) then
        cmfd_hold_weights = .true.
        call cmfd_tally_reset()
        cmfd_inact_flush(2) = cmfd_inact_flush(2) - 1
    end if

  end subroutine cmfd_init_batch

# ifdef PETSC

!==============================================================================
! PROCESS_CMFD_OPTIONS 
!==============================================================================

  subroutine process_cmfd_options()

    use global,       only: cmfd_snes_monitor, cmfd_ksp_monitor, mpi_err

    ! check for snes monitor
    if (cmfd_snes_monitor) call PetscOptionsSetValue("-snes_monitor", &
         "stdout", mpi_err)

    ! check for ksp monitor
    if (cmfd_ksp_monitor) call PetscOptionsSetValue("-ksp_monitor", &
         "stdout", mpi_err)

    end subroutine process_cmfd_options

!===============================================================================
! CALC_FISSION_SOURCE calculates the cmfd fission source
!===============================================================================

  subroutine calc_fission_source()

    use constants,  only: CMFD_NOACCEL, ZERO, TWO
    use global,     only: cmfd, cmfd_coremap, master, mpi_err, entropy_on

    integer :: nx ! maximum number of cells in x direction
    integer :: ny ! maximum number of cells in y direction
    integer :: nz ! maximum number of cells in z direction
    integer :: ng ! maximum number of energy groups
    integer :: n  ! total size
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
    n  = ng*nx*ny*nz

    ! allocate cmfd source if not already allocated and allocate buffer
    if (.not. allocated(cmfd%cmfd_src)) allocate(cmfd%cmfd_src(ng,nx,ny,nz))

    ! reset cmfd source to 0
    cmfd%cmfd_src = ZERO

    ! only perform for master
    if (master) then

      ! loop around indices to map to cmfd object
      ZLOOP: do k = 1, nz

        YLOOP: do j = 1, ny

          XLOOP: do i = 1, nx

            GROUP: do g = 1, ng

              ! check for core map
              if (cmfd_coremap) then
                if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
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
              cmfd%cmfd_src(g,i,j,k) = sum(cmfd%nfissxs(:,g,i,j,k) * &
                   cmfd%phi(idx:idx+(ng-1)))*vol

            end do GROUP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

      ! normalize source such that it sums to 1.0
      cmfd%cmfd_src = cmfd%cmfd_src/sum(cmfd%cmfd_src)

      ! compute entropy
      if (entropy_on) then

        ! allocate tmp array
        if (.not.allocated(source)) allocate(source(ng,nx,ny,nz))

        ! initialize the source
        source = ZERO

        ! compute log
        where (cmfd%cmfd_src > ZERO)
          source = cmfd%cmfd_src*log(cmfd%cmfd_src)/log(TWO)
        end where

        ! sum that source
        cmfd%entropy = -sum(source)

        ! deallocate tmp array
        if (allocated(source)) deallocate(source)

      end if

      ! normalize source so average is 1.0
      cmfd%cmfd_src = cmfd%cmfd_src/sum(cmfd%cmfd_src)*cmfd%norm

    end if

    ! broadcast full source to all procs
    call MPI_BCAST(cmfd%cmfd_src, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)

  end subroutine calc_fission_source

!===============================================================================
! CMFD_REWEIGHT
!===============================================================================

  subroutine cmfd_reweight(new_weights)

    use constants,   only: ZERO, ONE
    use error,       only: warning, fatal_error
    use global,      only: n_particles, meshes, source_bank, work,             &
                           n_user_meshes, message, cmfd, master, mpi_err,      &
                           bank_first, bank_last
    use mesh_header, only: StructuredMesh
    use mesh,        only: count_bank_sites, get_mesh_indices
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
    integer(8) :: size_bank ! size of source bank
    logical :: outside ! any source sites outside mesh
    logical :: in_mesh ! source site is inside mesh
    logical :: new_weights ! calcualte new weights
    type(StructuredMesh), pointer :: m ! point to mesh
    real(8), allocatable :: egrid(:)

    ! associate pointer
    m => meshes(n_user_meshes + 1)

    ! get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! compute size of source bank
    size_bank = bank_last - bank_first + 1_8 

    ! allocate arrays in cmfd object (can take out later extend to multigroup)
    if (.not.allocated(cmfd%sourcecounts)) then 
      allocate(cmfd%sourcecounts(ng,nx,ny,nz))
      cmfd % sourcecounts = 0
    end if
    if (.not.allocated(cmfd%weightfactors)) then 
      allocate(cmfd%weightfactors(ng,nx,ny,nz))
      cmfd % weightfactors = ONE
    end if

    ! allocate energy grid and reverse cmfd energy grid
    if (.not. allocated(egrid)) allocate(egrid(ng+1))
    egrid = (/(cmfd%egrid(ng-i+2),i = 1,ng+1)/)

    ! compute new weight factors
    if (new_weights) then

      ! zero out weights
      cmfd%weightfactors = ZERO

      ! count bank sites in mesh
      call count_bank_sites(m, source_bank, cmfd%sourcecounts, egrid, &
           sites_outside=outside, size_bank = size_bank)

      ! check for sites outside of the mesh
      if (master .and. outside) then
        message = "Source sites outside of the CMFD mesh!"
        call fatal_error()
      end if

      ! have master compute weight factors
      if (master) then
        where(cmfd%cmfd_src > ZERO .and. cmfd%sourcecounts > ZERO)
          cmfd%weightfactors = cmfd%cmfd_src/sum(cmfd%cmfd_src)* &
                               sum(cmfd%sourcecounts) / cmfd%sourcecounts
        end where
      end if

      ! broadcast weight factors to all procs
      call MPI_BCAST(cmfd%weightfactors, ng*nx*ny*nz, MPI_REAL8, 0, &
           MPI_COMM_WORLD, mpi_err)

   end if

    ! begin loop over source bank
    do i = 1, int(size_bank, 4) 

      ! determine spatial bin
      call get_mesh_indices(m, source_bank(i)%xyz, ijk, in_mesh)

      ! determine energy bin
      n_groups = size(cmfd%egrid) - 1
      if (source_bank(i) % E < cmfd%egrid(1)) then
        e_bin = 1
        message = 'source pt below energy grid'
        call warning()
      elseif (source_bank(i) % E > cmfd%egrid(n_groups+1)) then
        e_bin = n_groups
        message = 'source pt above energy grid'
        call warning()
      else
        e_bin = binary_search(cmfd%egrid, n_groups + 1, source_bank(i) % E)
      end if

      ! reverese energy bin (lowest grp is highest energy bin)
      e_bin = n_groups - e_bin + 1

      ! check for outside of mesh
      if (.not. in_mesh) then
        message = 'Source site found outside of CMFD mesh!'
        call fatal_error()
      end if

      ! reweight particle
      source_bank(i)%wgt = source_bank(i)%wgt * &
           cmfd%weightfactors(e_bin,ijk(1),ijk(2),ijk(3))

    end do

    ! deallocate
    if (allocated(egrid)) deallocate(egrid)

  end subroutine cmfd_reweight

!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix
!===============================================================================

  function get_matrix_idx(g, i, j, k, ng, nx, ny) result (matidx)

    use global, only: cmfd, cmfd_coremap

    integer :: matidx ! the index location in matrix
    integer :: i      ! current x index
    integer :: j      ! current y index
    integer :: k      ! current z index
    integer :: g      ! current group index
    integer :: nx     ! maximum number of cells in x direction
    integer :: ny     ! maximum number of cells in y direction
    integer :: ng     ! maximum number of energy groups

    ! check if coremap is used
    if (cmfd_coremap) then

      ! get idx from core map
      matidx = ng*(cmfd % coremap(i,j,k)) - (ng - g)

    else

      ! compute index
      matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if

  end function get_matrix_idx

# endif

!===============================================================================
! CMFD_TALLY_RESET
!===============================================================================

  subroutine cmfd_tally_reset()

    use global,  only: n_cmfd_tallies, cmfd_tallies, message
    use output,  only: write_message
    use tally,   only: reset_result

    integer :: i ! loop counter

    ! print message
    message = "CMFD tallies reset"
    call write_message(7)

    ! begin loop around CMFD tallies
    do i = 1, n_cmfd_tallies

      ! reset that tally
      cmfd_tallies(i) % n_realizations = 0
      call reset_result(cmfd_tallies(i) % results)

    end do

  end subroutine cmfd_tally_reset

end module cmfd_execute
