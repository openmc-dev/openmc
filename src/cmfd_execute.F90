module cmfd_execute 

!==============================================================================
! CMFD_EXECUTE -- This module is the highest level cmfd module that controls the
! cross section generation, diffusion calculation, and source re-weighting
!==============================================================================

  use global

  implicit none
  private
  public :: execute_cmfd, cmfd_init_batch

contains

!==============================================================================
! EXECUTE_CMFD runs the CMFD calculation
!==============================================================================

  subroutine execute_cmfd()

    use cmfd_data,              only: set_up_cmfd
    use cmfd_power_solver,      only: cmfd_power_execute
    use cmfd_jfnk_solver,       only: cmfd_jfnk_execute
    use error,                  only: warning, fatal_error 

    ! CMFD single processor on master
    if (master) then

      ! Start cmfd timer
      call time_cmfd % start()

      ! Create cmfd data from OpenMC tallies 
      call set_up_cmfd()

      ! Process solver options
      call process_cmfd_options()

      ! Call solver
      if (trim(cmfd_solver_type) == 'power') then
        call cmfd_power_execute()
      elseif (trim(cmfd_solver_type) == 'jfnk') then
        call cmfd_jfnk_execute()
      else
        message = 'solver type became invalid after input processing'
        call fatal_error() 
      end if

      ! Save k-effective
      cmfd % k_cmfd(current_batch) = cmfd % keff

      ! check to perform adjoint on last batch
      if (current_batch == n_batches .and. cmfd_run_adjoint) then
        if (trim(cmfd_solver_type) == 'power') then
          call cmfd_power_execute(adjoint = .true.)
        elseif (trim(cmfd_solver_type) == 'jfnk') then
          call cmfd_jfnk_execute(adjoint = .true.)
        end if
      end if

    end if

    ! calculate fission source
    call calc_fission_source()

    ! calculate weight factors
    if (cmfd_feedback) call cmfd_reweight(.true.)

    ! stop cmfd timer
    if (master) call time_cmfd % stop()

  end subroutine execute_cmfd

!==============================================================================
! CMFD_INIT_BATCH handles cmfd options at the start of every batch
!==============================================================================

  subroutine cmfd_init_batch()

    use global,            only: cmfd_begin, cmfd_on, cmfd_tally_on,         &
                                 cmfd_inact_flush, cmfd_act_flush, cmfd_run, &
                                 current_batch, cmfd_hold_weights

    ! Check to activate CMFD diffusion and possible feedback
    ! this guarantees that when cmfd begins at least one batch of tallies are
    ! accumulated
    if (cmfd_run .and. cmfd_begin == current_batch) then
      cmfd_on = .true.
      cmfd_tally_on = .true.
    end if

    ! If this is a restart run and we are just replaying batches leave
    if (restart_run .and. current_batch <= restart_batch) return

    ! Check to flush cmfd tallies for active batches, no more inactive flush
    if (cmfd_run .and. cmfd_act_flush == current_batch) then
      call cmfd_tally_reset()
      cmfd_tally_on = .true.
      cmfd_inact_flush(2) = -1
    end if

    ! Check to flush cmfd tallies during inactive batches (>= on number of
    ! flushes important as the code will flush on the first batch which we
    ! dont want to count)
    if (cmfd_run .and. mod(current_batch,cmfd_inact_flush(1))   &
       == 0 .and. cmfd_inact_flush(2) > 0 .and. cmfd_begin < current_batch) then
        cmfd_hold_weights = .true.
        call cmfd_tally_reset()
        cmfd_inact_flush(2) = cmfd_inact_flush(2) - 1
    end if

  end subroutine cmfd_init_batch

!==============================================================================
! PROCESS_CMFD_OPTIONS processes user options that interface with PETSc 
!==============================================================================

  subroutine process_cmfd_options()

#ifdef PETSC
    use global,       only: cmfd_snes_monitor, cmfd_ksp_monitor, mpi_err

    ! Check for snes monitor
    if (cmfd_snes_monitor) call PetscOptionsSetValue("-snes_monitor", &
         "stdout", mpi_err)

    ! Check for ksp monitor
    if (cmfd_ksp_monitor) call PetscOptionsSetValue("-ksp_monitor", &
         "stdout", mpi_err)
#endif

    end subroutine process_cmfd_options

!===============================================================================
! CALC_FISSION_SOURCE calculates the cmfd fission source
!===============================================================================

  subroutine calc_fission_source()

    use constants,  only: CMFD_NOACCEL, ZERO, TWO
    use global,     only: cmfd, cmfd_coremap, master, entropy_on, current_batch

#ifdef MPI
    use global,     only: mpi_err
    use mpi
#endif

    integer :: nx      ! maximum number of cells in x direction
    integer :: ny      ! maximum number of cells in y direction
    integer :: nz      ! maximum number of cells in z direction
    integer :: ng      ! maximum number of energy groups
    integer :: n       ! total size
    integer :: i       ! iteration counter for x
    integer :: j       ! iteration counter for y
    integer :: k       ! iteration counter for z
    integer :: g       ! iteration counter for groups
    integer :: idx     ! index in vector
    real(8) :: hxyz(3) ! cell dimensions of current ijk cell
    real(8) :: vol     ! volume of cell
    real(8),allocatable :: source(:,:,:,:)  ! tmp source array for entropy

    ! Get maximum of spatial and group indices
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)
    n  = ng*nx*ny*nz

    ! Allocate cmfd source if not already allocated and allocate buffer
    if (.not. allocated(cmfd % cmfd_src)) &
       allocate(cmfd % cmfd_src(ng,nx,ny,nz))

    ! Reset cmfd source to 0
    cmfd % cmfd_src = ZERO

    ! Only perform for master
    if (master) then

      ! Loop around indices to map to cmfd object
      ZLOOP: do k = 1, nz

        YLOOP: do j = 1, ny

          XLOOP: do i = 1, nx

            GROUP: do g = 1, ng

              ! Check for core map
              if (cmfd_coremap) then
                if (cmfd % coremap(i,j,k) == CMFD_NOACCEL) then
                  cycle
                end if
              end if

              ! Get dimensions of cell
              hxyz = cmfd % hxyz(:,i,j,k)

              ! Calculate volume
              vol = hxyz(1)*hxyz(2)*hxyz(3)

              ! Get first index
              idx = get_matrix_idx(1,i,j,k,ng,nx,ny)

              ! Compute fission source
              cmfd % cmfd_src(g,i,j,k) = sum(cmfd % nfissxs(:,g,i,j,k) * &
                     cmfd % phi(idx:idx + (ng - 1)))*vol

            end do GROUP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

      ! Normalize source such that it sums to 1.0
      cmfd % cmfd_src = cmfd % cmfd_src/sum(cmfd % cmfd_src)

      ! Compute entropy
      if (entropy_on) then

        ! Allocate tmp array
        if (.not.allocated(source)) allocate(source(ng,nx,ny,nz))

        ! Initialize the source
        source = ZERO

        ! Compute log
        where (cmfd % cmfd_src > ZERO)
          source = cmfd % cmfd_src*log(cmfd % cmfd_src)/log(TWO)
        end where

        ! Sum that source
        cmfd % entropy(current_batch) = -sum(source)

        ! Deallocate tmp array
        if (allocated(source)) deallocate(source)

      end if

      ! Normalize source so average is 1.0
      cmfd % cmfd_src = cmfd % cmfd_src/sum(cmfd % cmfd_src)*cmfd % norm

      ! Calculate differences between normalized sources
      cmfd % src_cmp(current_batch) = sqrt(ONE/cmfd % norm * &
             sum((cmfd % cmfd_src - cmfd % openmc_src)**2))

    end if

#ifdef MPI
    ! Broadcast full source to all procs
    call MPI_BCAST(cmfd % cmfd_src, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
#endif

  end subroutine calc_fission_source

!===============================================================================
! CMFD_REWEIGHT performs weighting of particles in the source bank
!===============================================================================

  subroutine cmfd_reweight(new_weights)

    use constants,   only: ZERO, ONE
    use error,       only: warning, fatal_error
    use global,      only: meshes, source_bank, work, n_user_meshes, message, &
                           cmfd, master
    use mesh_header, only: StructuredMesh
    use mesh,        only: count_bank_sites, get_mesh_indices
    use search,      only: binary_search

#ifdef MPI
    use global,      only: mpi_err
    use mpi
#endif

    logical, intent(in) :: new_weights ! calcualte new weights

    integer :: nx       ! maximum number of cells in x direction
    integer :: ny       ! maximum number of cells in y direction
    integer :: nz       ! maximum number of cells in z direction
    integer :: ng       ! maximum number of energy groups
    integer :: i        ! iteration counter
    integer :: ijk(3)   ! spatial bin location
    integer :: e_bin    ! energy bin of source particle
    integer :: n_groups ! number of energy groups
    logical :: outside  ! any source sites outside mesh
    logical :: in_mesh  ! source site is inside mesh

    type(StructuredMesh), pointer :: m ! point to mesh
    real(8), allocatable :: egrid(:)   ! energy grid

    ! Associate pointer
    m => meshes(n_user_meshes + 1)

    ! Get maximum of spatial and group indices
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! allocate arrays in cmfd object (can take out later extend to multigroup)
    if (.not.allocated(cmfd%sourcecounts)) then 
      allocate(cmfd%sourcecounts(ng,nx,ny,nz))
      cmfd % sourcecounts = 0
    end if
    if (.not.allocated(cmfd % weightfactors)) then 
      allocate(cmfd % weightfactors(ng,nx,ny,nz))
      cmfd % weightfactors = ONE
    end if

    ! Allocate energy grid and reverse cmfd energy grid
    if (.not. allocated(egrid)) allocate(egrid(ng + 1))
    egrid = (/(cmfd % egrid(ng - i + 2), i = 1, ng + 1)/)

    ! Compute new weight factors
    if (new_weights) then

      ! Zero out weights
      cmfd%weightfactors = ZERO

      ! Count bank sites in mesh
      call count_bank_sites(m, source_bank, cmfd%sourcecounts, egrid, &
           sites_outside=outside, size_bank=work)

      ! Check for sites outside of the mesh
      if (master .and. outside) then
        message = "Source sites outside of the CMFD mesh!"
        call fatal_error()
      end if

      ! Have master compute weight factors (watch for 0s)
      if (master) then
        where(cmfd % cmfd_src > ZERO .and. cmfd % sourcecounts > ZERO)
          cmfd % weightfactors = cmfd % cmfd_src/sum(cmfd % cmfd_src)* &
                               sum(cmfd % sourcecounts) / cmfd % sourcecounts
        end where
      end if

      ! Broadcast weight factors to all procs
#ifdef MPI
      call MPI_BCAST(cmfd % weightfactors, ng*nx*ny*nz, MPI_REAL8, 0, &
           MPI_COMM_WORLD, mpi_err)
#endif
   end if

    ! begin loop over source bank
    do i = 1, int(work,4)

      ! Determine spatial bin
      call get_mesh_indices(m, source_bank(i) % xyz, ijk, in_mesh)

      ! Determine energy bin
      n_groups = size(cmfd % egrid) - 1
      if (source_bank(i) % E < cmfd % egrid(1)) then
        e_bin = 1
        message = 'Source pt below energy grid'
        call warning()
      elseif (source_bank(i) % E > cmfd % egrid(n_groups + 1)) then
        e_bin = n_groups
        message = 'Source pt above energy grid'
        call warning()
      else
        e_bin = binary_search(cmfd % egrid, n_groups + 1, source_bank(i) % E)
      end if

      ! Reverese energy bin (lowest grp is highest energy bin)
      e_bin = n_groups - e_bin + 1

      ! Check for outside of mesh
      if (.not. in_mesh) then
        message = 'Source site found outside of CMFD mesh'
        call fatal_error()
      end if

      ! Reweight particle
      source_bank(i) % wgt = source_bank(i) % wgt * &
             cmfd % weightfactors(e_bin, ijk(1), ijk(2), ijk(3))

    end do

    ! Deallocate all
    if (allocated(egrid)) deallocate(egrid)

  end subroutine cmfd_reweight

!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix
!===============================================================================

  function get_matrix_idx(g, i, j, k, ng, nx, ny) result (matidx)

    use global, only: cmfd, cmfd_coremap

    integer :: matidx ! the index location in matrix
    integer, intent(in) :: i  ! current x index
    integer, intent(in) :: j  ! current y index
    integer, intent(in) :: k  ! current z index
    integer, intent(in) :: g  ! current group index
    integer, intent(in) :: nx ! maximum number of cells in x direction
    integer, intent(in) :: ny ! maximum number of cells in y direction
    integer, intent(in) :: ng ! maximum number of energy groups

    ! Check if coremap is used
    if (cmfd_coremap) then

      ! Get idx from core map
      matidx = ng*(cmfd % coremap(i,j,k)) - (ng - g)

    else

      ! Compute index
      matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if

  end function get_matrix_idx

!===============================================================================
! CMFD_TALLY_RESET resets all cmfd tallies
!===============================================================================

  subroutine cmfd_tally_reset()

    use global,  only: n_cmfd_tallies, cmfd_tallies, message
    use output,  only: write_message
    use tally,   only: reset_result

    integer :: i ! loop counter

    ! Print message
    message = "CMFD tallies reset"
    call write_message(7)

    ! Begin loop around CMFD tallies
    do i = 1, n_cmfd_tallies

      ! Reset that tally
      cmfd_tallies(i) % n_realizations = 0
      call reset_result(cmfd_tallies(i) % results)

    end do

  end subroutine cmfd_tally_reset

end module cmfd_execute
