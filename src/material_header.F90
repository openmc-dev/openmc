module material_header

  use, intrinsic :: ISO_C_BINDING

  implicit none

  private
  public :: material_id
  public :: material_nuclide
  public :: material_nuclide_size
  public :: material_atom_density
  public :: material_density_gpcc

  interface
    function material_id(i_mat) bind(C) result(id)
      import C_INT32_T
      integer(C_INT32_T), value :: i_mat
      integer(C_INT32_T)        :: id
    end function

    function material_nuclide(i_mat, idx) bind(C) result(nuc)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value :: i_mat
      integer(C_INT), value :: idx
      integer(C_INT) :: nuc
    end function

    function material_nuclide_size(i_mat) bind(C) result(n)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value :: i_mat
      integer(C_INT) :: n
    end function

    function material_atom_density(i_mat, idx) bind(C) result(density)
      import C_INT32_T, C_INT, C_DOUBLE
      integer(C_INT32_T), value :: i_mat
      integer(C_INT), value :: idx
      real(C_DOUBLE) :: density
    end function

    function material_density_gpcc(i_mat) bind(C) result(density)
      import C_INT32_T, C_DOUBLE
      integer(C_INT32_T), value :: i_mat
      real(C_DOUBLE) :: density
    end function
  end interface

end module material_header
