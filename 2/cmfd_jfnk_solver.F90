module cmfd_jfnk_solver

  use cmfd_loss_operator,  only: init_loss_matrix, build_loss_matrix
  use cmfd_power_solver,   only: cmfd_power_execute
  use cmfd_prod_operator,  only: init_prod_matrix, build_prod_matrix
  use matrix_header,       only: Matrix
  use solver_interface,    only: JFNKSolver, Jfnk_ctx
  use vector_header,       only: Vector

  implicit none
  private
  public :: cmfd_jfnk_execute

  logical          :: adjoint_calc ! true if an adjoint is to be calculated
  type(Jfnk_ctx)   :: jfnk_data    ! object that holds pointers to user routines
  type(Matrix)     :: jac_prec     ! Jacobian preconditioner object
  type(Matrix)     :: loss         ! CMFD loss matrix
  type(Matrix)     :: prod         ! CMFD production matrix
#ifdef PETSC
  type(JFNKSolver) :: jfnk         ! JFNK solver object
#endif
  type(Vector)     :: resvec       ! JFNK residual vector
  type(Vector)     :: xvec         ! JFNK solution vector

contains

!===============================================================================
! CMFD_JFNK_EXECUTE main routine for JFNK solver
!===============================================================================

  subroutine cmfd_jfnk_execute(adjoint)

    use global,  only: time_cmfdbuild, time_cmfdsolve

    logical, intent(in), optional :: adjoint ! adjoint calculation

    ! Check for adjoint
    adjoint_calc = .false.
    if (present(adjoint)) adjoint_calc = adjoint

    ! Seed with power iteration to help converge to fundamental mode
    call cmfd_power_execute(k_tol=1.E-3_8, s_tol=1.E-3_8, adjoint=adjoint_calc)

    ! Start timer for build
    call time_cmfdbuild % start()

    ! Initialize data
    call init_data()

    ! Initialize solver
#ifdef PETSC
    call jfnk % create()
#endif

    ! Set up residual and jacobian routines
#ifdef PETSC
    jfnk_data % res_proc_ptr => compute_nonlinear_residual
    jfnk_data % jac_proc_ptr => build_jacobian_matrix
    call jfnk % set_functions(jfnk_data, resvec, jac_prec)
#endif

    ! Stop timer for build
    call time_cmfdbuild % stop()

    ! Solve the system
    call time_cmfdsolve % start()
#ifdef PETSC
    call jfnk % solve(xvec)
#endif
    call time_cmfdsolve % stop()

    ! Extracts results to cmfd object
    call extract_results()

    ! Deallocate all slepc data
    call finalize()

  end subroutine cmfd_jfnk_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

    use constants, only: ZERO, ONE
    use global,    only: cmfd, cmfd_adjoint_type, current_batch

    logical :: physical_adjoint ! physical adjoing calculation logical
    integer :: n ! size of matrices

    ! Set up all matrices
    call init_loss_matrix(loss)
    call init_prod_matrix(prod)
    call init_jacobian_matrix()

    ! Check for physical adjoint
    physical_adjoint = .false.
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'physical') &
         physical_adjoint = .true.

    ! Create matrix operators
    call build_loss_matrix(loss, adjoint = physical_adjoint)
    call build_prod_matrix(prod, adjoint = physical_adjoint)

    ! Assemble matrices and use PETSc
    call loss % assemble()
    call prod % assemble()
#ifdef PETSC
    call loss % setup_petsc()
    call prod % setup_petsc()
#endif

    ! Check for mathematical adjoint
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'math') &
         call compute_adjoint()

    ! Get problem size
    n = loss % n

    ! Create problem vectors
    call resvec % create(n + 1)
    call xvec % create(n + 1)

    ! Set flux in guess from rough power iteration
    if (adjoint_calc) then
      xvec % val(1:n) = cmfd % adj_phi
    else
      xvec % val(1:n) = cmfd % phi
    end if

    ! Set keff in guess from rough power iteration
    if (adjoint_calc) then
      xvec % val(n + 1) = ONE/cmfd % adj_keff
    else
      xvec % val(n + 1) = ONE/cmfd % keff
    end if

    ! Set up vectors for PETSc 
#ifdef PETSC
    call resvec   % setup_petsc()
    call xvec     % setup_petsc()
#endif

    ! Build jacobian from initial guess
    call build_jacobian_matrix(xvec)

    ! Set up Jacobian for PETSc
#ifdef PETSC
    call jac_prec % setup_petsc()
#endif

    ! Set dominance ratio to 0
    cmfd % dom(current_batch) = ZERO

  end subroutine init_data

!==============================================================================
! INIT_JACOBIAN_MATRIX preallocates jacobian matrix and initializes it
!==============================================================================

  subroutine init_jacobian_matrix()

    integer :: nnz ! number of nonzeros in matrix
    integer :: n   ! dimension of matrix

    ! Get length of matrix and number of nonzeros total in loss matrix
    nnz = loss % nnz
    n   = loss % n

    ! There is one more nonzero for each row and and additional row of nonzeros
    nnz = nnz + 2*n + 1

    ! We need 1 more row for the jacobian
    n = n + 1

    ! Configure Jacobian matrix
    call jac_prec % create(n, nnz)

  end subroutine init_jacobian_matrix

!===============================================================================
! BUILD_JACOBIAN_MATRIX creates the Jacobian of nonlinear eigenvalue problem
!===============================================================================

  subroutine build_jacobian_matrix(x)

    use constants,  only: ONE

    type(Vector), intent(in) :: x ! solution vector

    integer      :: i      ! loop counter for jacobian rows
    integer      :: jjac   ! loop counter for jacobian cols
    integer      :: jloss  ! loop counter for loss matrix cols
    integer      :: jprod  ! loop counter for prod matrix cols
    integer      :: n      ! problem size
    real(8)      :: lambda ! eigenvalue
    real(8)      :: val    ! temporary real scalar
    type(Vector) :: flux   ! flux vector
    type(Vector) :: fsrc   ! fission source vector

    ! Get the problem size
    n = loss % n

    ! Get flux and eigenvalue
    flux % val => x % val(1:n)
    lambda = x % val(n + 1)

    ! Create fission source vector
    call fsrc % create(n)
    call prod % vector_multiply(flux, fsrc)

    ! Reset counters in Jacobian matrix
    jac_prec % n_count  = 1
    jac_prec % nz_count = 1

    ! Begin loop around rows of Jacobian matrix
    ROWS: do i = 1, n

      ! Add a new row in Jacobian
      call jac_prec % new_row()

      ! Begin loop around columns of loss matrix
      COLS_LOSS: do jloss = loss % get_row(i), loss % get_row(i + 1) - 1

        ! Start with the value in the loss matrix
        val = loss % val(jloss)

        ! Loop around columns of prod matrix
        COLS_PROD: do jprod = prod % get_row(i), prod % get_row(i + 1) - 1

          ! See if columns agree with loss matrix
          if (prod % get_col(jprod) == loss % get_col(jloss)) then
            val = val - lambda*prod % val(jprod)
            exit
          end if

        end do COLS_PROD 

        ! Record value in jacobian
        call jac_prec % add_value(loss % get_col(jloss), val)

      end do COLS_LOSS

      ! Add fission source value
      val = -fsrc % val(i)
      call jac_prec % add_value(n + 1, val)

    end do ROWS

    ! Need to add negative transpose of flux vector on last row
    call jac_prec % new_row()
    do jjac = 1, n
      val = -flux % val(jjac)
      call jac_prec % add_value(jjac, val)
    end do

    ! Add unity on bottom right corner of matrix
    call jac_prec % add_value(n + 1, ONE)

    ! CRS requires a final value in row
    call jac_prec % new_row()

    ! Free all allocated memory
    flux % val => null()
    call fsrc % destroy()

  end subroutine build_jacobian_matrix

!===============================================================================
! COMPUTE_NONLINEAR_RESIDUAL computes the residual of the nonlinear equations
!===============================================================================
#ifdef PETSC
  subroutine compute_nonlinear_residual(x, res)

    use global,       only: cmfd_write_matrices

    type(Vector), intent(in)    :: x   ! solution vector
    type(Vector), intent(inout) :: res ! residual vector

    character(len=25) :: filename
    integer           :: n
    real(8)           :: lambda
    type(Vector)      :: res_loss
    type(Vector)      :: res_prod
    type(Vector)      :: flux

    ! Get problem size
    n = loss % n

    ! Set up temporary vectors
    call res_loss % create(n)
    call res_prod % create(n)

    ! Extract flux
    flux % n = n
    flux % val => x % val(1:n)

    ! Extract eigenvalue
    lambda = x % val(n + 1)

    ! Calculate M*flux then F*flux
    call loss % vector_multiply(flux, res_loss)
    call prod % vector_multiply(flux, res_prod)

    ! Put flux component values in residual vector
    res % val(1:n) = res_loss % val - lambda*res_prod % val

    ! Put eigenvalue component in residual vector
    res % val(n+1) = 0.5_8 - 0.5_8*sum(flux % val * flux % val)

    ! Write out data in binary files (debugging)
    if (cmfd_write_matrices) then

      ! Write out residual vector 
      if (adjoint_calc) then
        filename = 'adj_res.bin'
      else
        filename = 'res.bin'
      end if
#ifdef PETSC
      call res % write_petsc_binary(filename)
#endif

      ! Write out solution vector
      if (adjoint_calc) then
        filename = 'adj_x.bin'
      else
        filename = 'x.bin'
      end if
#ifdef PETSC
      call x % write_petsc_binary(filename)
#endif

    end if

  end subroutine compute_nonlinear_residual
#endif

!===============================================================================
! COMPUTE_ADJOINT calculates mathematical adjoint by taking transposes
!===============================================================================

  subroutine compute_adjoint()

    use global,  only: cmfd_write_matrices

    ! Transpose matrices
#ifdef PETSC
    call loss % transpose()
    call prod % transpose()
#endif

    ! Write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
#ifdef PETSC
      call loss % write_petsc_binary('adj_lossmat.bin')
      call loss % write_petsc_binary('adj_prodmat.bin')
#endif
    end if

  end subroutine compute_adjoint

!===============================================================================
! EXTRACT_RESULTS gets the results and puts them in global CMFD object
!===============================================================================

  subroutine extract_results()

    use constants, only: ZERO, ONE
    use global,    only: cmfd

    integer :: n ! problem size

    ! Get problem size
    n = loss % n

    ! Also allocate in cmfd object
    if (adjoint_calc) then
      if (.not. allocated(cmfd % adj_phi)) allocate(cmfd % adj_phi(n))
    else
      if (.not. allocated(cmfd % phi)) allocate(cmfd % phi(n))
    end if

    ! Get flux and eigenvalue
    if (adjoint_calc) then
      cmfd % adj_phi  = xvec % val(1:n)
      cmfd % adj_keff = ONE / xvec % val(n+1)
    else
      cmfd % phi  = xvec % val(1:n)
      cmfd % keff = ONE / xvec % val(n+1)
    end if

  end subroutine extract_results

!===============================================================================
! FINALIZE frees all memory from a JFNK calculation
!===============================================================================

  subroutine finalize()

    call loss % destroy()
    call prod % destroy()
    call jac_prec % destroy()
    call xvec % destroy()
    call resvec % destroy()
#ifdef PETSC
    call jfnk % destroy()
#endif
    nullify(jfnk_data % res_proc_ptr)
    nullify(jfnk_data % jac_proc_ptr)

  end subroutine finalize

end module cmfd_jfnk_solver
