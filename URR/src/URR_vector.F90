module URR_vector

  implicit none
  private
  public :: VectorReal1D,&
            VectorReal2D,&
            VectorInt1D


!> Vector of reals
  type :: VectorReal1D

    real(8), allocatable :: dim1(:)

  end type VectorReal1D


!> Vector of Vectors of reals
  type :: VectorReal2D

    type(VectorReal1D), allocatable :: dim2(:)

  end type VectorReal2D


!> Vector of integers
  type :: VectorInt1D

    integer, allocatable :: dim1(:)

  end type VectorInt1D


end module URR_vector
