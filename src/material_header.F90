module material_header

  use, intrinsic :: ISO_C_BINDING

  use constants
  use dict_header, only: DictIntInt
  use error
  use nuclide_header
  use sab_header
  use stl_vector, only: VectorReal, VectorInt
  use string, only: to_str

  implicit none

  private
  public :: free_memory_material
  public :: openmc_extend_materials
  public :: openmc_get_material_index
  public :: openmc_material_add_nuclide
  public :: openmc_material_get_id
  public :: openmc_material_get_densities
  public :: openmc_material_set_density
  public :: openmc_material_set_densities
  public :: openmc_material_set_id

!===============================================================================
! MATERIAL describes a material by its constituent nuclides
!===============================================================================

  type, public :: Material
    integer              :: id              ! unique identifier
    character(len=104)   :: name = ""       ! User-defined name
    integer              :: n_nuclides = 0  ! number of nuclides
    integer, allocatable :: nuclide(:)      ! index in nuclides array
    real(8)              :: density         ! total atom density in atom/b-cm
    real(8), allocatable :: atom_density(:) ! nuclide atom density in atom/b-cm
    real(8)              :: density_gpcc    ! total density in g/cm^3

    ! To improve performance of tallying, we store an array (direct address
    ! table) that indicates for each nuclide in the global nuclides(:) array the
    ! index of the corresponding nuclide in the Material % nuclide(:) array. If
    ! it is not present in the material, the entry is set to zero.
    integer, allocatable :: mat_nuclide_index(:)

    ! S(a,b) data
    integer              :: n_sab = 0         ! number of S(a,b) tables
    integer, allocatable :: i_sab_nuclides(:) ! index of corresponding nuclide
    integer, allocatable :: i_sab_tables(:)   ! index in sab_tables
    real(8), allocatable :: sab_fracs(:)      ! how often to use S(a,b)

    ! Temporary names read during initialization
    character(20), allocatable :: names(:)     ! isotope names
    character(20), allocatable :: sab_names(:) ! name of S(a,b) table

    ! Does this material contain fissionable nuclides? Is it depletable?
    logical :: fissionable = .false.
    logical :: depletable = .false.

    ! enforce isotropic scattering in lab for specific nuclides
    logical :: has_isotropic_nuclides = .false.
    logical, allocatable :: p0(:)

  contains
    procedure :: set_density => material_set_density
    procedure :: init_nuclide_index => material_init_nuclide_index
    procedure :: assign_sab_tables => material_assign_sab_tables
  end type Material

  integer(C_INT32_T), public, bind(C) :: n_materials ! # of materials

  type(Material), public, allocatable, target :: materials(:)

  ! Dictionary that maps user IDs to indices in 'materials'
  type(DictIntInt), public :: material_dict

contains

!===============================================================================
! MATERIAL_SET_DENSITY sets the total density of a material in atom/b-cm.
!===============================================================================

  function material_set_density(this, density) result(err)
    class(Material), intent(inout) :: this
    real(8), intent(in) :: density
    integer :: err

    integer :: i
    real(8) :: sum_percent
    real(8) :: awr

    if (allocated(this % atom_density)) then
      ! Set total density based on value provided
      this % density = density

      ! Determine normalized atom percents
      sum_percent = sum(this % atom_density)
      this % atom_density(:) = this % atom_density / sum_percent

      ! Recalculate nuclide atom densities based on given density
      this % atom_density(:) = density * this % atom_density

      ! Calculate density in g/cm^3.
      this % density_gpcc = ZERO
      do i = 1, this % n_nuclides
        awr = nuclides(this % nuclide(i)) % awr
        this % density_gpcc = this % density_gpcc &
             + this % atom_density(i) * awr * MASS_NEUTRON / N_AVOGADRO
      end do
      err = 0
    else
      err = E_ALLOCATE
      call set_errmsg("Material atom density array hasn't been allocated.")
    end if
  end function material_set_density

!===============================================================================
! INIT_NUCLIDE_INDEX creates a mapping from indices in the global nuclides(:)
! array to the Material % nuclides array
!===============================================================================

  subroutine material_init_nuclide_index(this)
    class(Material), intent(inout) :: this

    integer :: i

    ! Allocate nuclide index array and set to zeros
    if (allocated(this % mat_nuclide_index)) &
         deallocate(this % mat_nuclide_index)
    allocate(this % mat_nuclide_index(n_nuclides))
    this % mat_nuclide_index(:) = 0

    ! Assign entries in the index array
    do i = 1, this % n_nuclides
      this % mat_nuclide_index(this % nuclide(i)) = i
    end do
  end subroutine material_init_nuclide_index

!===============================================================================
! ASSIGN_SAB_TABLES assigns S(alpha,beta) tables to specific nuclides within
! materials so the code knows when to apply bound thermal scattering data
!===============================================================================

  subroutine material_assign_sab_tables(this)
    class(Material), intent(inout) :: this

    integer :: j            ! index over nuclides in material
    integer :: k            ! index over S(a,b) tables in material
    integer :: m            ! position for sorting
    integer :: temp_nuclide ! temporary value for sorting
    integer :: temp_table   ! temporary value for sorting
    real(8) :: temp_frac    ! temporary value for sorting
    logical :: found
    type(VectorInt)  :: i_sab_tables
    type(VectorInt)  :: i_sab_nuclides
    type(VectorReal) :: sab_fracs

    if (.not. allocated(this % i_sab_tables)) return

    ASSIGN_SAB: do k = 1, size(this % i_sab_tables)
      ! In order to know which nuclide the S(a,b) table applies to, we need
      ! to search through the list of nuclides for one which has a matching
      ! name
      found = .false.
      associate (sab => sab_tables(this % i_sab_tables(k)))
        FIND_NUCLIDE: do j = 1, size(this % nuclide)
          if (any(sab % nuclides == nuclides(this % nuclide(j)) % name)) then
            call i_sab_tables % push_back(this % i_sab_tables(k))
            call i_sab_nuclides % push_back(j)
            call sab_fracs % push_back(this % sab_fracs(k))
            found = .true.
          end if
        end do FIND_NUCLIDE
      end associate

      ! Check to make sure S(a,b) table matched a nuclide
      if (.not. found) then
        call fatal_error("S(a,b) table " // trim(this % &
             sab_names(k)) // " did not match any nuclide on material " &
             // trim(to_str(this % id)))
      end if
    end do ASSIGN_SAB

    ! Make sure each nuclide only appears in one table.
    do j = 1, i_sab_nuclides % size()
      do k = j+1, i_sab_nuclides % size()
        if (i_sab_nuclides % data(j) == i_sab_nuclides % data(k)) then
          call fatal_error(trim( &
               nuclides(this % nuclide(i_sab_nuclides % data(j))) % name) &
               // " in material " // trim(to_str(this % id)) // " was found &
               &in multiple S(a,b) tables. Each nuclide can only appear in &
               &one S(a,b) table per material.")
        end if
      end do
    end do

    ! Update i_sab_tables and i_sab_nuclides
    deallocate(this % i_sab_tables)
    deallocate(this % sab_fracs)
    if (allocated(this % i_sab_nuclides)) deallocate(this % i_sab_nuclides)
    m = i_sab_tables % size()
    allocate(this % i_sab_tables(m))
    allocate(this % i_sab_nuclides(m))
    allocate(this % sab_fracs(m))
    this % i_sab_tables(:) = i_sab_tables % data(1:m)
    this % i_sab_nuclides(:) = i_sab_nuclides % data(1:m)
    this % sab_fracs(:) = sab_fracs % data(1:m)

    ! Clear entries in vectors for next material
    call i_sab_tables % clear()
    call i_sab_nuclides % clear()
    call sab_fracs % clear()

    ! If there are multiple S(a,b) tables, we need to make sure that the
    ! entries in i_sab_nuclides are sorted or else they won't be applied
    ! correctly in the cross_section module. The algorithm here is a simple
    ! insertion sort -- don't need anything fancy!

    if (size(this % i_sab_tables) > 1) then
      SORT_SAB: do k = 2, size(this % i_sab_tables)
        ! Save value to move
        m = k
        temp_nuclide = this % i_sab_nuclides(k)
        temp_table   = this % i_sab_tables(k)
        temp_frac    = this % i_sab_tables(k)

        MOVE_OVER: do
          ! Check if insertion value is greater than (m-1)th value
          if (temp_nuclide >= this % i_sab_nuclides(m-1)) exit

          ! Move values over until hitting one that's not larger
          this % i_sab_nuclides(m) = this % i_sab_nuclides(m-1)
          this % i_sab_tables(m)   = this % i_sab_tables(m-1)
          this % sab_fracs(m)      = this % sab_fracs(m-1)
          m = m - 1

          ! Exit if we've reached the beginning of the list
          if (m == 1) exit
        end do MOVE_OVER

        ! Put the original value into its new position
        this % i_sab_nuclides(m) = temp_nuclide
        this % i_sab_tables(m)   = temp_table
        this % sab_fracs(m)      = temp_frac
      end do SORT_SAB
    end if

    ! Deallocate temporary arrays for names of nuclides and S(a,b) tables
    if (allocated(this % names)) deallocate(this % names)
  end subroutine material_assign_sab_tables

!===============================================================================
! FREE_MEMORY_MATERIAL deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_material()
    n_materials = 0
    if (allocated(materials)) deallocate(materials)
    call material_dict % clear()
  end subroutine free_memory_material

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_extend_materials(n, index_start, index_end) result(err) bind(C)
    ! Extend the materials array by n elements
    integer(C_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), optional, intent(out) :: index_start
    integer(C_INT32_T), optional, intent(out) :: index_end
    integer(C_INT) :: err

    type(Material), allocatable :: temp(:) ! temporary materials array

    if (n_materials == 0) then
      ! Allocate materials array
      allocate(materials(n))
    else
      ! Allocate materials array with increased size
      allocate(temp(n_materials + n))

      ! Move original materials to temporary array
      temp(1:n_materials) = materials(:)

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=materials)
    end if

    ! Return indices in materials array
    if (present(index_start)) index_start = n_materials + 1
    if (present(index_end)) index_end = n_materials + n
    n_materials = n_materials + n

    err = 0
  end function openmc_extend_materials

  function openmc_get_material_index(id, index) result(err) bind(C)
    ! Returns the index in the materials array of a material with a given ID
    integer(C_INT32_T), value :: id
    integer(C_INT32_T), intent(out) :: index
    integer(C_INT) :: err

    if (allocated(materials)) then
      if (material_dict % has(id)) then
        index = material_dict % get(id)
        err = 0
      else
        err = E_INVALID_ID
        call set_errmsg("No material exists with ID=" // trim(to_str(id)) // ".")
      end if
    else
      err = E_ALLOCATE
      call set_errmsg("Memory has not been allocated for materials.")
    end if
  end function openmc_get_material_index


  function openmc_material_add_nuclide(index, name, density) result(err) bind(C)
    ! Add a nuclide at a specified density in atom/b-cm to a material
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR) :: name(*)
    real(C_DOUBLE), value, intent(in) :: density
    integer(C_INT) :: err

    integer :: j, k, n
    real(8) :: awr
    integer, allocatable :: new_nuclide(:)
    real(8), allocatable :: new_density(:)
    character(:), allocatable :: name_

    name_ = to_f_string(name)

    err = E_UNASSIGNED
    if (index >= 1 .and. index <= size(materials)) then
      associate (m => materials(index))
        ! Check if nuclide is already in material
        do j = 1, m % n_nuclides
          k = m % nuclide(j)
          if (nuclides(k) % name == name_) then
            awr = nuclides(k) % awr
            m % density = m % density + density - m % atom_density(j)
            m % density_gpcc = m % density_gpcc + (density - &
                 m % atom_density(j)) * awr * MASS_NEUTRON / N_AVOGADRO
            m % atom_density(j) = density
            err = 0
          end if
        end do

        ! If nuclide wasn't found, extend nuclide/density arrays
        if (err /= 0) then
          ! If nuclide hasn't been loaded, load it now
          err = openmc_load_nuclide(name)

          if (err == 0) then
            ! Extend arrays
            n = m % n_nuclides
            allocate(new_nuclide(n + 1))
            if (n > 0) new_nuclide(1:n) = m % nuclide
            call move_alloc(FROM=new_nuclide, TO=m % nuclide)

            allocate(new_density(n + 1))
            if (n > 0) new_density(1:n) = m % atom_density
            call move_alloc(FROM=new_density, TO=m % atom_density)

            ! Append new nuclide/density
            k = nuclide_dict % get(to_lower(name_))
            m % nuclide(n + 1) = k
            m % atom_density(n + 1) = density
            m % density = m % density + density
            m % density_gpcc = m % density_gpcc + &
                 density * nuclides(k) % awr * MASS_NEUTRON / N_AVOGADRO
            m % n_nuclides = n + 1
          end if
        end if
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in materials array is out of bounds.")
    end if

  end function openmc_material_add_nuclide


  function openmc_material_get_densities(index, nuclides, densities, n) &
       result(err) bind(C)
    ! returns an array of nuclide densities in a material
    integer(C_INT32_T), value :: index
    type(C_PTR),        intent(out) :: nuclides
    type(C_PTR),        intent(out) :: densities
    integer(C_INT),     intent(out) :: n
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(materials)) then
      associate (m => materials(index))
        if (allocated(m % atom_density)) then
          nuclides = C_LOC(m % nuclide(1))
          densities = C_LOC(m % atom_density(1))
          n = size(m % atom_density)
          err = 0
        else
          err = E_ALLOCATE
          call set_errmsg("Material atom density array has not been allocated.")
        end if
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in materials array is out of bounds.")
    end if
  end function openmc_material_get_densities


  function openmc_material_get_id(index, id) result(err) bind(C)
    ! returns the ID of a material
    integer(C_INT32_T), value       :: index
    integer(C_INT32_T), intent(out) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(materials)) then
      id = materials(index) % id
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in materials array is out of bounds.")
    end if
  end function openmc_material_get_id


  function openmc_material_set_id(index, id) result(err) bind(C)
    ! Set the ID of a material
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_materials) then
      materials(index) % id = id
      call material_dict % set(id, index)
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in materials array is out of bounds.")
    end if
  end function openmc_material_set_id


  function openmc_material_set_density(index, density) result(err) bind(C)
    ! Set the total density of a material in atom/b-cm
    integer(C_INT32_T), value, intent(in) :: index
    real(C_DOUBLE), value, intent(in) :: density
    integer(C_INT) :: err

    err = E_UNASSIGNED
    if (index >= 1 .and. index <= size(materials)) then
      associate (m => materials(index))
        err = m % set_density(density)
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in materials array is out of bounds.")
    end if
  end function openmc_material_set_density


  function openmc_material_set_densities(index, n, name, density) result(err) bind(C)
    ! Sets the densities for a list of nuclides in a material. If the nuclides
    ! don't already exist in the material, they will be added
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT), value, intent(in) :: n
    type(C_PTR),    intent(in) :: name(n)
    real(C_DOUBLE), intent(in) :: density(n)
    integer(C_INT) :: err

    integer :: i
    integer :: stat
    character(C_CHAR), pointer :: string(:)
    character(len=:, kind=C_CHAR), allocatable :: name_

    if (index >= 1 .and. index <= size(materials)) then
      associate (m => materials(index))
        ! If nuclide/density arrays are not correct size, reallocate
        if (n /= m % n_nuclides) then
          deallocate(m % nuclide, m % atom_density, STAT=stat)
          allocate(m % nuclide(n), m % atom_density(n))
        end if

        do i = 1, n
          ! Convert C string to Fortran string
          call c_f_pointer(name(i), string, [10])
          name_ = to_lower(to_f_string(string))

          if (.not. nuclide_dict % has(name_)) then
            err = openmc_load_nuclide(string)
            if (err < 0) return
          end if

          m % nuclide(i) = nuclide_dict % get(name_)
          m % atom_density(i) = density(i)
        end do
        m % n_nuclides = n

        ! Set total density to the sum of the vector
        err = m % set_density(sum(density))

        ! Assign S(a,b) tables
        call m % assign_sab_tables()
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in materials array is out of bounds.")
    end if

  end function openmc_material_set_densities

end module material_header
