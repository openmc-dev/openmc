module cmfd_jfnk_solver

  use cmfd_loss_operator,  only: init_loss_matrix, build_loss_matrix
  use cmfd_power_solver,   only: cmfd_power_execute
  use cmfd_prod_operator,  only: init_prod_matrix, build_prod_matrix
  use matrix_header,       only: Matrix
  use petsc_solver,        only: Petsc_jfnk, Jfnk_ctx
  use vector_header,       only: Vector

  implicit none
  private
  public :: cmfd_jfnk_execute

  logical          :: adjoint_calc = .false. ! adjoint calculation
  type(Jfnk_ctx)   :: jfnk_data
  type(Matrix)     :: jac_prec
  type(Matrix)     :: loss
  type(Matrix)     :: prod
  type(Petsc_jfnk) :: jfnk 
  type(Vector)     :: resvec
  type(Vector)     :: xvec

contains

!===============================================================================
! CMFD_JFNK_EXECUTE
!===============================================================================

  subroutine cmfd_jfnk_execute(adjoint)

    logical, optional :: adjoint ! adjoint calculation

    ! check for adjoint
    if (present(adjoint)) adjoint_calc = adjoint

    ! seed with power iteration
    call cmfd_power_execute(k_tol=1.E-3_8, s_tol=1.E-3_8, adjoint=adjoint_calc)

    ! initialize data
    call init_data()

    ! initialize solver
    call jfnk % create()

    ! solve the system
!   call SNESSolve(snes, PETSC_NULL_DOUBLE, xvec, ierr)

    ! extracts results to cmfd object
    call extract_results()

    ! deallocate all slepc data
    call finalize()

  end subroutine cmfd_jfnk_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

    use constants, only: ONE
    use global,    only: cmfd

    integer              :: n         ! problem size

    ! Set up all matrices
    call init_loss_matrix(loss)
    call init_prod_matrix(prod)
    call init_jacobian_matrix()

    ! Set up for use with petsc
    call jac_prec % setup_petsc()

    ! Get problem size
    n = jac_prec % n

    ! Create problem vectors
    call resvec % create(n)
    call xvec % create(n)

    ! set flux in guess
    if (adjoint_calc) then
      xvec % val(1:n) = cmfd % adj_phi
    else
      xvec % val(1:n) = cmfd % phi
    end if

    ! set keff in guess
    if (adjoint_calc) then
      xvec % val(n + 1) = ONE/cmfd % adj_keff
    else
      xvec % val(n + 1) = ONE/cmfd % keff
    end if

    ! set solver names

  end subroutine init_data

!==============================================================================
! INIT_JACOBIAN_MATRIX preallocates jacobian matrix and initializes it
!==============================================================================

  subroutine init_jacobian_matrix()

    integer :: nnz
    integer :: n

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

    type(Vector) :: x

    integer      :: i      ! loop counter for jacobian rows
    integer      :: jjac   ! loop counter for jacobian cols
    integer      :: jloss  ! loop counter for loss matrix cols
    integer      :: jprod  ! loop counter for prod matrix cols
    integer      :: n      ! problem size
    real(8)      :: lambda ! eigenvalue
    real(8)      :: val    ! temporary real scalar
    type(Vector) :: flux   ! flux vector
    type(Vector) :: fsrc   ! fission source vector

    ! get the problem size
    n = loss % n

    ! get flux and eigenvalue
    flux % val => x % val(1:n)
    lambda = x % val(n + 1)

    ! create fission source vector
    call fsrc % create(n)
    call prod % vector_multiply(flux, fsrc)

    ! reset counters in Jacobian matrix
    jac_prec % n_kount  = 1
    jac_prec % nz_kount = 1

    ! begin loop around rows of Jacobian matrix
    ROWS: do i = 1, n

      ! add a new row in Jacobian
      call jac_prec % new_row()

      ! begin loop around columns of loss matrix
      COLS_LOSS: do jloss = loss % row(i), loss % row(i + 1) - 1

        ! start with the value in the loss matrix
        val = loss % val(jloss)

        ! loop around columns of prod matrix
        COLS_PROD: do jprod = prod % row(i), prod % row(i + 1) - 1

          ! see if columns agree with loss matrix
          if (prod % col(jprod) == loss % col(jloss)) then
            val = val - lambda*prod % val(jprod)
            exit
          end if

        end do COLS_PROD 

        ! record value in jacobian
        call jac_prec % add_value(jloss, val)

      end do COLS_LOSS

      ! add fission source value
      val = fsrc % val(n)
      call jac_prec % add_value(n + 1, val)

    end do ROWS

    ! need to add negative transpose of flux vector on last row
    call jac_prec % new_row()
    do jjac = 1, n
      val = -flux % val(jjac)
      call jac_prec % add_value(jjac, val)
    end do

    ! add unity on bottom right corner of matrix
    call jac_prec % add_value(n + 1, ONE)

    ! free all allocated memory
    flux % val => null()
    call fsrc % destroy()

  end subroutine build_jacobian_matrix

!===============================================================================
! COMPUTE_NONLINEAR_RESIDUAL
!===============================================================================

  subroutine compute_nonlinear_residual(x)

    use global,       only: cmfd_write_matrices, cmfd_adjoint_type

    type(Vector)      :: x 

    character(len=25) :: filename
    integer           :: n
    logical           :: physical_adjoint = .false.
    real(8)           :: lambda
    type(Vector)      :: res_loss
    type(Vector)      :: res_prod
    type(Vector)      :: flux

    ! Check for physical adjoint
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'physical') &
         physical_adjoint = .true.

    ! Create operators
    call build_loss_matrix(loss, adjoint = physical_adjoint)
    call build_prod_matrix(prod, adjoint = physical_adjoint)

    ! Check for adjoint
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'math') &
         call compute_adjoint()

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
    resvec % val(1:n) = res_loss % val - lambda*res_prod % val

    ! Put eigenvalue component in residual vector
    resvec % val(n+1) = 0.5_8 - 0.5_8*sum(flux % val * flux % val)

    ! write out data in binary files (debugging)
    if (cmfd_write_matrices) then

      ! write out residual vector 
      if (adjoint_calc) then
        filename = 'adj_res.bin'
      else
        filename = 'res.bin'
      end if
      call resvec % write_petsc_binary(filename)

      ! write out solution vector
      if (adjoint_calc) then
        filename = 'adj_x.bin'
      else
        filename = 'x.bin'
      end if
      call x % write_petsc_binary(filename)

    end if

  end subroutine compute_nonlinear_residual

!===============================================================================
! COMPUTE_ADJOINT
!===============================================================================

  subroutine compute_adjoint()

    use global,  only: cmfd_write_matrices

    ! transpose matrices
    call loss % transpose()
    call prod % transpose()

    ! write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
      call loss % write_petsc_binary('adj_lossmat.bin')
      call loss % write_petsc_binary('adj_prodmat.bin')
    end if

  end subroutine compute_adjoint

!===============================================================================
! EXTRACT_RESULTS 
!===============================================================================

  subroutine extract_results()

    use constants, only: ZERO, ONE
    use global,    only: cmfd

    integer :: n ! problem size

    ! get problem size
    n = jfnk_data % loss % n

    ! also allocate in cmfd object
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
      cmfd % keff = xvec % val(n+1)
    end if

  end subroutine extract_results

!===============================================================================
! FINALIZE frees all memory from a JFNK calculation
!===============================================================================

  subroutine finalize()

    call jfnk_data % loss % destroy()
    call jfnk_data % prod % destroy()
    call jac_prec % destroy()
    call xvec % destroy()
    call resvec % destroy()
    call jfnk % destroy() 

  end subroutine finalize

end module cmfd_jfnk_solver
