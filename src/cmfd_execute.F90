module cmfd_execute

!==============================================================================
! CMFD_EXECUTE -- This module is the highest level cmfd module that controls the
! cross section generation, diffusion calculation, and source re-weighting
!==============================================================================

  use cmfd_header
  use settings
  use simulation_header

  implicit none
  private
  public :: execute_cmfd, cmfd_init_batch

contains

!==============================================================================
! EXECUTE_CMFD runs the CMFD calculation
!==============================================================================

  subroutine execute_cmfd()

    use cmfd_data,              only: set_up_cmfd
    use cmfd_solver,            only: cmfd_solver_execute
    use error,                  only: warning, fatal_error
    use message_passing,        only: master

    ! CMFD single processor on master
    if (master) then

      ! Start cmfd timer
      call time_cmfd % start()

      ! Create cmfd data from OpenMC tallies
      call set_up_cmfd()

      ! Call solver
      call cmfd_solver_execute()

      ! Save k-effective
      cmfd % k_cmfd(current_batch) = cmfd % keff

      ! check to perform adjoint on last batch
      if (current_batch == n_batches .and. cmfd_run_adjoint) then
        call cmfd_solver_execute(adjoint=.true.)
      end if

    end if

    ! calculate fission source
    call calc_fission_source()

    ! calculate weight factors
    call cmfd_reweight(.true.)

    ! stop cmfd timer
    if (master) call time_cmfd % stop()

  end subroutine execute_cmfd

!==============================================================================
! CMFD_INIT_BATCH handles cmfd options at the start of every batch
!==============================================================================

  subroutine cmfd_init_batch()

    ! Check to activate CMFD diffusion and possible feedback
    ! this guarantees that when cmfd begins at least one batch of tallies are
    ! accumulated
    if (cmfd_run .and. cmfd_begin == current_batch) then
      cmfd_on = .true.
    end if

    ! If this is a restart run and we are just replaying batches leave
    if (restart_run .and. current_batch <= restart_batch) return

    ! Check to reset tallies
    if (cmfd_run .and. cmfd_reset % contains(current_batch)) then
      call cmfd_tally_reset()
    end if

  end subroutine cmfd_init_batch

!===============================================================================
! CALC_FISSION_SOURCE calculates the cmfd fission source
!===============================================================================

  subroutine calc_fission_source()

    use constants, only: CMFD_NOACCEL, ZERO, TWO
    use message_passing
    use string,    only: to_str

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
#ifdef OPENMC_MPI
    integer :: mpi_err ! MPI error code
#endif

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

#ifdef OPENMC_MPI
    ! Broadcast full source to all procs
    call MPI_BCAST(cmfd % cmfd_src, n, MPI_REAL8, 0, mpi_intracomm, mpi_err)
#endif

  end subroutine calc_fission_source

!===============================================================================
! CMFD_REWEIGHT performs weighting of particles in the source bank
!===============================================================================

  subroutine cmfd_reweight(new_weights)

    use algorithm,   only: binary_search
    use bank_header, only: source_bank
    use constants,   only: ZERO, ONE
    use error,       only: warning, fatal_error
    use mesh_header, only: RegularMesh
    use mesh,        only: count_bank_sites
    use message_passing
    use string,      only: to_str

    logical, intent(in) :: new_weights ! calcualte new weights

    integer :: nx       ! maximum number of cells in x direction
    integer :: ny       ! maximum number of cells in y direction
    integer :: nz       ! maximum number of cells in z direction
    integer :: ng       ! maximum number of energy groups
    integer :: i        ! iteration counter
    integer :: g        ! index for group
    integer :: ijk(3)   ! spatial bin location
    integer :: e_bin    ! energy bin of source particle
    integer :: mesh_bin ! mesh bin of soruce particle
    integer :: n_groups ! number of energy groups
    real(8) :: norm     ! normalization factor
    logical :: outside  ! any source sites outside mesh
    logical :: in_mesh  ! source site is inside mesh
#ifdef OPENMC_MPI
    integer :: mpi_err
#endif

    ! Get maximum of spatial and group indices
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! allocate arrays in cmfd object (can take out later extend to multigroup)
    if (.not.allocated(cmfd%sourcecounts)) then
      allocate(cmfd%sourcecounts(ng, nx*ny*nz))
      cmfd % sourcecounts = 0
    end if
    if (.not.allocated(cmfd % weightfactors)) then
      allocate(cmfd % weightfactors(ng,nx,ny,nz))
      cmfd % weightfactors = ONE
    end if

    ! Compute new weight factors
    if (new_weights) then

      ! Set weight factors to a default 1.0
      cmfd%weightfactors = ONE

      ! Count bank sites in mesh and reverse due to egrid structure
      call count_bank_sites(cmfd_mesh, source_bank, cmfd%sourcecounts, &
           cmfd % egrid, sites_outside=outside, size_bank=work)

      ! Check for sites outside of the mesh
      if (master .and. outside) then
        call fatal_error("Source sites outside of the CMFD mesh!")
      end if

      ! Have master compute weight factors (watch for 0s)
      if (master) then
        ! Calculate normalization factor
        norm = sum(cmfd % sourcecounts) / sum(cmfd % cmfd_src)

        do mesh_bin = 1, nx*ny*nz
          call cmfd_mesh % get_indices_from_bin(mesh_bin, ijk)
          do g = 1, ng
            if (cmfd % sourcecounts(ng - g + 1, mesh_bin) > ZERO) then
              if (cmfd % cmfd_src(g,ijk(1),ijk(2),ijk(3)) > ZERO) then
                cmfd % weightfactors(g,ijk(1),ijk(2),ijk(3)) = &
                     cmfd % cmfd_src(g,ijk(1),ijk(2),ijk(3)) * norm &
                     / cmfd % sourcecounts(ng - g + 1, mesh_bin)
              end if
            end if
          end do
        end do
      end if

      if (.not. cmfd_feedback) return

      ! Broadcast weight factors to all procs
#ifdef OPENMC_MPI
      call MPI_BCAST(cmfd % weightfactors, ng*nx*ny*nz, MPI_REAL8, 0, &
           mpi_intracomm, mpi_err)
#endif
    end if

    ! begin loop over source bank
    do i = 1, int(work,4)

      ! Determine spatial bin
      call cmfd_mesh % get_indices(source_bank(i) % xyz, ijk, in_mesh)

      ! Determine energy bin
      n_groups = size(cmfd % egrid) - 1
      if (source_bank(i) % E < cmfd % egrid(1)) then
        e_bin = 1
        if (master) call warning('Source pt below energy grid')
      elseif (source_bank(i) % E > cmfd % egrid(n_groups + 1)) then
        e_bin = n_groups
        if (master) call warning('Source pt above energy grid')
      else
        e_bin = binary_search(cmfd % egrid, n_groups + 1, source_bank(i) % E)
      end if

      ! Reverese energy bin (lowest grp is highest energy bin)
      e_bin = n_groups - e_bin + 1

      ! Check for outside of mesh
      if (.not. in_mesh) then
        call fatal_error('Source site found outside of CMFD mesh')
      end if

      ! Reweight particle
      source_bank(i) % wgt = source_bank(i) % wgt * &
           cmfd % weightfactors(e_bin, ijk(1), ijk(2), ijk(3))
    end do

  end subroutine cmfd_reweight

!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix
!===============================================================================

  function get_matrix_idx(g, i, j, k, ng, nx, ny) result (matidx)

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

    use error,  only: write_message

    integer :: i ! loop counter

    ! Print message
    call write_message("CMFD tallies reset", 6)

    ! Reset CMFD tallies
    do i = 1, size(cmfd_tallies)
      cmfd_tallies(i) % obj % n_realizations = 0
      cmfd_tallies(i) % obj % results(:,:,:) = ZERO
    end do

  end subroutine cmfd_tally_reset

end module cmfd_execute
