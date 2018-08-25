module material_header

  use, intrinsic :: ISO_C_BINDING

  use constants
  use dict_header, only: DictIntInt
  use error
  use math, only: spline, spline_integrate
  use nuclide_header
  use particle_header, only: Particle
  use photon_header
  use sab_header
  use simulation_header, only: log_spacing
  use stl_vector, only: VectorReal, VectorInt
  use string, only: to_str

  implicit none

  private
  public :: bremsstrahlung_init
  public :: free_memory_material
  public :: openmc_extend_materials
  public :: openmc_get_material_index
  public :: openmc_material_add_nuclide
  public :: openmc_material_get_id
  public :: openmc_material_get_densities
  public :: openmc_material_get_volume
  public :: openmc_material_set_density
  public :: openmc_material_set_densities
  public :: openmc_material_set_id
  public :: material_pointer

  interface
    function material_pointer(mat_ind) bind(C) result(ptr)
      import C_PTR, C_INT32_T
      integer(C_INT32_T), intent(in), value :: mat_ind
      type(C_PTR)                           :: ptr
    end function material_pointer

    function material_id_c(mat_ptr) bind(C, name='material_id') result(id)
      import C_PTR, C_INT32_T
      type(C_PTR), intent(in), value :: mat_ptr
      integer(C_INT32_T)             :: id
    end function material_id_c

    subroutine material_set_id_c(mat_ptr, id, index) &
         bind(C, name='material_set_id')
      import C_PTR, C_INT32_T
      type(C_PTR),        intent(in), value :: mat_ptr
      integer(C_INT32_T), intent(in), value :: id
      integer(C_INT32_T), intent(in), value :: index
    end subroutine material_set_id_c

    subroutine extend_materials_c(n) bind(C)
      import C_INT32_T
      integer(C_INT32_T), intent(in), value :: n
    end subroutine extend_materials_c

    function openmc_material_get_volume(index, volume) result(err) bind(C)
      import C_INT32_T, C_DOUBLE, C_INT
      integer(C_INT32_T), value :: index
      real(C_DOUBLE), intent(out) :: volume
      integer(C_INT) :: err
    end function openmc_material_get_volume
  end interface

!===============================================================================
! MATERIAL describes a material by its constituent nuclides
!===============================================================================

  type, public :: Material
    type(C_PTR) :: ptr
    character(len=104)   :: name = ""       ! User-defined name
    integer              :: n_nuclides = 0  ! number of nuclides
    integer, allocatable :: nuclide(:)      ! index in nuclides array
    integer, allocatable :: element(:)      ! index in elements array
    real(8)              :: density         ! total atom density in atom/b-cm
    real(C_DOUBLE), allocatable :: atom_density(:) ! nuclide atom density in atom/b-cm
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
    procedure :: id => material_id
    procedure :: set_id => material_set_id
    procedure :: set_density => material_set_density
    procedure :: init_nuclide_index => material_init_nuclide_index
    procedure :: assign_sab_tables => material_assign_sab_tables
    procedure :: calculate_xs => material_calculate_xs
    procedure, private :: calculate_neutron_xs
    procedure, private :: calculate_photon_xs
  end type Material

  integer(C_INT32_T), public, bind(C) :: n_materials ! # of materials

  type(Material), public, allocatable, target :: materials(:)

  ! Dictionary that maps user IDs to indices in 'materials'
  type(DictIntInt), public :: material_dict

contains

!===============================================================================
! MATERIAL_SET_DENSITY sets the total density of a material in atom/b-cm.
!===============================================================================

  function material_id(this) result(id)
    class(Material), intent(in) :: this
    integer(C_INT32_T)          :: id
    id = material_id_c(this % ptr)
  end function material_id

  subroutine material_set_id(this, id, index)
    class(Material),    intent(in) :: this
    integer(C_INT32_T), intent(in) :: id
    integer(C_INT32_T), intent(in) :: index
    call material_set_id_c(this % ptr, id, index)
  end subroutine material_set_id

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
          if (sab % has_nuclide(nuclides(this % nuclide(j)) % name)) then
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
             // trim(to_str(this % id())))
      end if
    end do ASSIGN_SAB

    ! Make sure each nuclide only appears in one table.
    do j = 1, i_sab_nuclides % size()
      do k = j+1, i_sab_nuclides % size()
        if (i_sab_nuclides % data(j) == i_sab_nuclides % data(k)) then
          call fatal_error(trim( &
               nuclides(this % nuclide(i_sab_nuclides % data(j))) % name) &
               // " in material " // trim(to_str(this % id())) // " was found &
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
! MATERIAL_CALCULATE_XS determines the macroscopic cross sections for the
! material the particle is currently traveling through.
!===============================================================================

  subroutine material_calculate_xs(this, p)
    class(Material), intent(in) :: this
    type(Particle),  intent(in) :: p

    ! Set all material macroscopic cross sections to zero
    material_xs % total           = ZERO
    material_xs % absorption      = ZERO
    material_xs % fission         = ZERO
    material_xs % nu_fission      = ZERO

    if (p % type == NEUTRON) then
      call this % calculate_neutron_xs(p)
    elseif (p % type == PHOTON) then
      call this % calculate_photon_xs(p)
    end if

  end subroutine material_calculate_xs

!===============================================================================
! CALCULATE_NEUTRON_XS determines the neutron cross section for the material the
! particle is traveling through
!===============================================================================

  subroutine calculate_neutron_xs(this, p)
    class(Material), intent(in) :: this
    type(Particle),  intent(in) :: p

    integer :: i             ! loop index over nuclides
    integer :: i_nuclide     ! index into nuclides array
    integer :: i_sab         ! index into sab_tables array
    integer :: j             ! index in this % i_sab_nuclides
    integer :: i_grid        ! index into logarithmic mapping array or material
                             ! union grid
    real(8) :: atom_density  ! atom density of a nuclide
    real(8) :: sab_frac      ! fraction of atoms affected by S(a,b)
    logical :: check_sab     ! should we check for S(a,b) table?

    ! Find energy index on energy grid
    i_grid = int(log(p % E/energy_min(NEUTRON))/log_spacing)

    ! Determine if this material has S(a,b) tables
    check_sab = (this % n_sab > 0)

    ! Initialize position in i_sab_nuclides
    j = 1

    ! Add contribution from each nuclide in material
    do i = 1, this % n_nuclides
      ! ======================================================================
      ! CHECK FOR S(A,B) TABLE

      i_sab = 0
      sab_frac = ZERO

      ! Check if this nuclide matches one of the S(a,b) tables specified.
      ! This relies on i_sab_nuclides being in sorted order
      if (check_sab) then
        if (i == this % i_sab_nuclides(j)) then
          ! Get index in sab_tables
          i_sab = this % i_sab_tables(j)
          sab_frac = this % sab_fracs(j)

          ! If particle energy is greater than the highest energy for the
          ! S(a,b) table, then don't use the S(a,b) table
          if (p % E > sab_tables(i_sab) % threshold()) then
            i_sab = 0
          end if

          ! Increment position in i_sab_nuclides
          j = j + 1

          ! Don't check for S(a,b) tables if there are no more left
          if (j > size(this % i_sab_tables)) check_sab = .false.
        end if
      end if

      ! ======================================================================
      ! CALCULATE MICROSCOPIC CROSS SECTION

      ! Determine microscopic cross sections for this nuclide
      i_nuclide = this % nuclide(i)

      ! Calculate microscopic cross section for this nuclide
      if (p % E /= micro_xs(i_nuclide) % last_E &
           .or. p % sqrtkT /= micro_xs(i_nuclide) % last_sqrtkT &
           .or. i_sab /= micro_xs(i_nuclide) % index_sab &
           .or. sab_frac /= micro_xs(i_nuclide) % sab_frac) then
        call nuclides(i_nuclide) % calculate_xs(i_sab, p % E, i_grid, &
             p % sqrtkT, sab_frac, micro_xs(i_nuclide))
      end if

      ! ======================================================================
      ! ADD TO MACROSCOPIC CROSS SECTION

      ! Copy atom density of nuclide in material
      atom_density = this % atom_density(i)

      ! Add contributions to material macroscopic total cross section
      material_xs % total = material_xs % total + &
           atom_density * micro_xs(i_nuclide) % total

      ! Add contributions to material macroscopic absorption cross section
      material_xs % absorption = material_xs % absorption + &
           atom_density * micro_xs(i_nuclide) % absorption

      ! Add contributions to material macroscopic fission cross section
      material_xs % fission = material_xs % fission + &
           atom_density * micro_xs(i_nuclide) % fission

      ! Add contributions to material macroscopic nu-fission cross section
      material_xs % nu_fission = material_xs % nu_fission + &
           atom_density * micro_xs(i_nuclide) % nu_fission
    end do

  end subroutine calculate_neutron_xs

!===============================================================================
! CALCULATE_PHOTON_XS determines the macroscopic photon cross sections for the
! material the particle is currently traveling through.
!===============================================================================

  subroutine calculate_photon_xs(this, p)
    class(Material), intent(in) :: this
    type(Particle),  intent(in) :: p

    integer :: i             ! loop index over nuclides
    integer :: i_element     ! index into elements array
    real(8) :: atom_density  ! atom density of a nuclide

    material_xs % coherent        = ZERO
    material_xs % incoherent      = ZERO
    material_xs % photoelectric   = ZERO
    material_xs % pair_production = ZERO

    ! Add contribution from each nuclide in material
    do i = 1, this % n_nuclides
      ! ========================================================================
      ! CALCULATE MICROSCOPIC CROSS SECTION

      ! Determine microscopic cross sections for this nuclide
      i_element = this % element(i)

      ! Calculate microscopic cross section for this nuclide
      if (p % E /= micro_photon_xs(i_element) % last_E) then
        call elements(i_element) % calculate_xs(&
             p % E, micro_photon_xs(i_element))
      end if

      ! ========================================================================
      ! ADD TO MACROSCOPIC CROSS SECTION

      ! Copy atom density of nuclide in material
      atom_density = this % atom_density(i)

      ! Add contributions to material macroscopic total cross section
      material_xs % total = material_xs % total + &
           atom_density * micro_photon_xs(i_element) % total

      ! Add contributions to material macroscopic coherent cross section
      material_xs % coherent = material_xs % coherent + &
           atom_density * micro_photon_xs(i_element) % coherent

      ! Add contributions to material macroscopic incoherent cross section
      material_xs % incoherent = material_xs % incoherent + &
           atom_density * micro_photon_xs(i_element) % incoherent

      ! Add contributions to material macroscopic photoelectric cross section
      material_xs % photoelectric = material_xs % photoelectric + &
           atom_density * micro_photon_xs(i_element) % photoelectric

      ! Add contributions to material macroscopic pair production cross section
      material_xs % pair_production = material_xs % pair_production + &
           atom_density * micro_photon_xs(i_element) % pair_production
    end do

  end subroutine calculate_photon_xs

!===============================================================================
! FREE_MEMORY_MATERIAL deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_material()
    interface
      subroutine free_memory_material_c() bind(C)
      end subroutine free_memory_material_c
    end interface
    call free_memory_material_c()
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

    integer :: i
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

    ! Extend the C++ materials array and get pointers to the C++ objects
    call extend_materials_c(n)
    do i = n_materials - n, n_materials
      materials(i) % ptr = material_pointer(i - 1)
    end do

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
      id = materials(index) % id()
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
      call materials(index) % set_id(id, index)
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

  subroutine bremsstrahlung_init(this, i_material, particle)
    class(BremsstrahlungData), intent(inout) :: this
    integer, intent(in) :: i_material
    integer, intent(in) :: particle

    integer                 :: i, j
    integer                 :: i_k
    integer                 :: n, n_e, n_k
    real(8)                 :: c
    real(8)                 :: k, k_l, k_r
    real(8)                 :: e, e_l, e_r
    real(8)                 :: w, w_l, w_r
    real(8)                 :: x, x_l, x_r
    real(8)                 :: t
    real(8)                 :: r
    real(8)                 :: awr
    real(8)                 :: beta
    real(8)                 :: Z_eq_sq
    real(8)                 :: atom_density
    real(8)                 :: mass_density
    real(8)                 :: sum_density
    real(8), allocatable    :: stopping_power_collision(:)
    real(8), allocatable    :: stopping_power_radiative(:)
    real(8), allocatable    :: stopping_power(:)
    real(8), allocatable    :: dcs(:,:)
    real(8), allocatable    :: f(:)
    real(8), allocatable    :: z(:)
    logical                 :: positron_
    type(Material), pointer :: mat
    type(PhotonInteraction), pointer :: elm

    ! Get pointer to this material
    mat => materials(i_material)

    ! Determine whether we are generating electron or positron data
    positron_ = (particle == POSITRON)

    ! Get the size of the energy grids
    n_k = size(ttb_k_grid)
    n_e = size(ttb_e_grid)

    ! Allocate arrays for TTB data
    allocate(this % pdf(n_e, n_e), source=ZERO)
    allocate(this % cdf(n_e, n_e), source=ZERO)
    allocate(this % yield(n_e))

    ! Allocate temporary arrays
    allocate(stopping_power_collision(n_e), source=ZERO)
    allocate(stopping_power_radiative(n_e), source=ZERO)
    allocate(stopping_power(n_e))
    allocate(dcs(n_k, n_e), source=ZERO)
    allocate(f(n_e))
    allocate(z(n_e))

    Z_eq_sq = ZERO
    sum_density = ZERO

    ! Calculate the molecular DCS and the molecular total stopping power using
    ! Bragg's additivity rule.
    ! TODO: The collision stopping power cannot be accurately calculated using
    ! Bragg's additivity rule since the mean excitation energies and the
    ! density effect corrections cannot simply be summed together. Bragg's
    ! additivity rule fails especially when a higher-density compound is
    ! composed of elements that are in lower-density form at normal temperature
    ! and pressure (at which the NIST stopping powers are given). It will be
    ! used to approximate the collision stopping powers for now, but should be
    ! fixed in the future.
    do i = 1, mat % n_nuclides
      ! Get pointer to current element
      elm => elements(mat % element(i))

      awr = nuclides(mat % nuclide(i)) % awr

      ! Get atomic density and mass density of nuclide given atom percent
      if (mat % atom_density(1) > ZERO) then
        atom_density = mat % atom_density(i)
        mass_density = mat % atom_density(i) * awr
      ! Given weight percent
      else
        atom_density = -mat % atom_density(i) / awr
        mass_density = -mat % atom_density(i)
      end if

      ! Calculate the "equivalent" atomic number Zeq of the material
      Z_eq_sq = Z_eq_sq + atom_density * elm % Z**2
      sum_density = sum_density + atom_density

      ! Accumulate material DCS
      dcs = dcs + atom_density * elm % Z**2 * elm % dcs

      ! Accumulate material collision stopping power
      stopping_power_collision = stopping_power_collision + mass_density &
           * MASS_NEUTRON / N_AVOGADRO * elm % stopping_power_collision

      ! Accumulate material radiative stopping power
      stopping_power_radiative = stopping_power_radiative + mass_density &
           * MASS_NEUTRON / N_AVOGADRO * elm % stopping_power_radiative
    end do
    Z_eq_sq = Z_eq_sq / sum_density

    ! Calculate the positron DCS and radiative stopping power. These are
    ! obtained by multiplying the electron DCS and radiative stopping powers by
    ! a factor r, which is a numerical approximation of the ratio of the
    ! radiative stopping powers for positrons and electrons. Source: F. Salvat,
    ! J. M. Fernández-Varea, and J. Sempau, "PENELOPE-2011: A Code System for
    ! Monte Carlo Simulation of Electron and Photon Transport," OECD-NEA,
    ! Issy-les-Moulineaux, France (2011).
    if (positron_) then
      do i = 1, n_e
        t = log(ONE + 1.0e6_8*ttb_e_grid(i)/(Z_eq_sq*MASS_ELECTRON_EV))
        r = ONE - exp(-1.2359e-1_8*t + 6.1274e-2_8*t**2 - 3.1516e-2_8*t**3 + &
             7.7446e-3_8*t**4 - 1.0595e-3_8*t**5 + 7.0568e-5_8*t**6 - &
             1.808e-6_8*t**7)
        stopping_power_radiative(i) = r*stopping_power_radiative(i)
        dcs(:,i) = r*dcs(:,i)
      end do
    end if

    ! Total material stopping power
    stopping_power = stopping_power_collision + stopping_power_radiative

    ! Loop over photon energies
    do i = 1, n_e - 1
      w = ttb_e_grid(i)

      ! Loop over incident particle energies
      do j = i, n_e
        e = ttb_e_grid(j)

        ! Reduced photon energy
        k = w / e

        ! Find the lower bounding index of the reduced photon energy
        i_k = binary_search(ttb_k_grid, n_k, k)

        ! Get the interpolation bounds
        k_l = ttb_k_grid(i_k)
        k_r = ttb_k_grid(i_k+1)
        x_l = dcs(i_k, j)
        x_r = dcs(i_k+1, j)

        ! Find the value of the DCS using linear interpolation in reduced
        ! photon energy k
        x = x_l + (k - k_l) * (x_r - x_l) / (k_r - k_l)

        ! Ratio of the velocity of the charged particle to the speed of light
        beta = sqrt(e*(e + TWO*MASS_ELECTRON_EV)) / (e + MASS_ELECTRON_EV)

        ! Compute the integrand of the PDF
        f(j) = x / (beta**2 * stopping_power(j) * w)
      end do

      ! Number of points to integrate
      n = n_e - i + 1

      ! Integrate the PDF using cubic spline integration over the incident
      ! particle energy
      if (n > 2) then
        call spline(n, ttb_e_grid(i:), f(i:), z(i:))

        c = ZERO
        do j = i, n_e - 1
          c = c + spline_integrate(n, ttb_e_grid(i:), f(i:), z(i:), &
               ttb_e_grid(j), ttb_e_grid(j+1))
          this % pdf(i,j+1) = c
        end do

      ! Integrate the last two points using trapezoidal rule in log-log space
      else
        e_l = log(ttb_e_grid(i))
        e_r = log(ttb_e_grid(i+1))
        x_l = log(f(i))
        x_r = log(f(i+1))

        this % pdf(i,i+1) = HALF * (e_r - e_l) * (exp(e_l + x_l) + exp(e_r + x_r))
      end if
    end do

    ! Loop over incident particle energies
    do j = 2, n_e
      ! Set last element of PDF to small non-zero value to enable log-log
      ! interpolation
      this % pdf(j,j) = exp(-500.0_8)

      ! Loop over photon energies
      c = ZERO
      do i = 1, j - 1
        ! Integrate the CDF from the PDF using the trapezoidal rule in log-log
        ! space
        w_l = log(ttb_e_grid(i))
        w_r = log(ttb_e_grid(i+1))
        x_l = log(this % pdf(i,j))
        x_r = log(this % pdf(i+1,j))

        c = c + HALF * (w_r - w_l) * (exp(w_l + x_l) + exp(w_r + x_r))
        this % cdf(i+1,j) = c
      end do

      ! Set photon number yield
      this % yield(j) = c
    end do

    ! Use logarithm of number yield since it is log-log interpolated
    where (this % yield > ZERO)
      this % yield = log(this % yield)
    elsewhere
      this % yield = -500.0_8
    end where

  end subroutine bremsstrahlung_init

end module material_header
