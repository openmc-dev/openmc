module mgxs_data

  use, intrinsic :: ISO_C_BINDING

  use constants
  use algorithm,       only: find
  use dict_header,     only: DictCharInt
  use error,           only: fatal_error, write_message
  use geometry_header, only: cells
  use hdf5_interface
  use material_header, only: Material, materials, n_materials
  use mgxs_interface
  use nuclide_header,  only: n_nuclides
  use set_header,      only: SetChar
  use settings
  use stl_vector,      only: VectorReal, VectorChar
  use string,          only: to_lower
  implicit none

contains

!===============================================================================
! CREATE_MACRO_XS generates the macroscopic xs from the microscopic input data
!===============================================================================

  subroutine create_macro_xs() bind(C)
    integer                       :: i_mat ! index in materials array
    type(Material), pointer       :: mat   ! current material
    type(VectorReal), allocatable :: kTs(:)
    character(MAX_WORD_LEN)       :: name  ! name of material

    ! Get temperatures to read for each material
    call get_mat_kTs(kTs)

    ! Force all nuclides in a material to be the same representation.
    ! Therefore type(nuclides(mat % nuclide(1)) % obj) dictates type(macroxs).
    ! At the same time, we will find the scattering type, as that will dictate
    ! how we allocate the scatter object within macroxs.allocate(macro_xs(n_materials))
    do i_mat = 1, n_materials

      ! Get the material
      mat => materials(i_mat)

      name = trim(mat % name) // C_NULL_CHAR

      call create_macro_xs_c(name, mat % n_nuclides, mat % nuclide, &
           kTs(i_mat) % size(), kTs(i_mat) % data, mat % atom_density, &
           temperature_tolerance, temperature_method)
    end do

  end subroutine create_macro_xs

!===============================================================================
! GET_MAT_kTs returns a list of temperatures (in eV) that each
! material appears at in the model.
!===============================================================================

  subroutine get_mat_kTs(kTs)

    type(VectorReal), allocatable, intent(out) :: kTs(:)
    integer :: i, j        ! Cell and material index
    integer :: i_material  ! Index in materials array
    real(8) :: kT          ! temperature in eV

    allocate(kTs(size(materials)))

    do i = 1, size(cells)
      ! Skip non-material cells
      if (cells(i) % fill() /= C_NONE) cycle

      do j = 1, cells(i) % material_size()

        ! Skip void materials
        if (cells(i) % material(j) == MATERIAL_VOID) cycle

        ! Get temperature of cell (rounding to nearest integer)
        if (cells(i) % sqrtkT_size() > 1) then
          kT = cells(i) % sqrtkT(j-1)**2
        else
          kT = cells(i) % sqrtkT(0)**2
        end if

        i_material = cells(i) % material(j)

        ! Add temperature if it hasn't already been added
        if (find(kTs(i_material), kT) == -1) then
          call kTs(i_material) % push_back(kT)
        end if

      end do
    end do

  end subroutine get_mat_kTs

end module mgxs_data
