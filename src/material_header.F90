module material_header

  use, intrinsic :: ISO_C_BINDING

  use constants
  use dict_header, only: DictIntInt
  use error
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
    integer, allocatable :: element(:)      ! index in elements array
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
          if (p % E > sab_tables(i_sab) % data(1) % threshold_inelastic) then
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

  subroutine bremsstrahlung_init(this, i_material)
    class(Bremsstrahlung), intent(inout) :: this
    integer, intent(in) :: i_material

    integer                 :: i, j
    integer                 :: i_k
    integer                 :: n_e, n_k
    real(8)                 :: e
    real(8)                 :: c
    real(8)                 :: k, k_l, k_r, k_c
    real(8)                 :: x_l, x_r, x_c
    real(8)                 :: awr
    real(8)                 :: density
    real(8)                 :: density_gpcc
    real(8)                 :: Z_eq_sq
    real(8)                 :: beta
    real(8)                 :: atom_sum
    real(8)                 :: mass_sum
    real(8), allocatable    :: atom_fraction(:)
    real(8), allocatable    :: mass_fraction(:)
    real(8), allocatable    :: stopping_power(:)
    real(8), allocatable    :: mfp_inv(:)
    real(8), allocatable    :: z(:)
    type(Material), pointer :: mat
    type(PhotonInteraction), pointer :: elm

    ! Get pointer to this material
    mat => materials(i_material)
    this % i_material = i_material

    ! Allocate and initialize arrays
    n_k = size(ttb_k_grid)
    n_e = size(ttb_e_grid)
    allocate(atom_fraction(mat % n_nuclides))
    allocate(mass_fraction(mat % n_nuclides))
    allocate(stopping_power(n_e))
    allocate(mfp_inv(n_e))
    allocate(this % yield(n_e))
    allocate(this % dcs(n_k, n_e))
    allocate(this % cdf(n_k, n_e))
    allocate(z(n_e))
    stopping_power(:) = ZERO
    mfp_inv(:) = ZERO
    this % dcs(:,:) = ZERO
    this % cdf(:,:) = ZERO

    ! Calculate the "equivalent" atomic number Zeq, the atomic fraction and the
    ! mass fraction of each element, and the material density in atom/b-cm and
    ! in g/cm^3
    Z_eq_sq = ZERO
    do i = 1, mat % n_nuclides
      awr = nuclides(mat % nuclide(i)) % awr

      ! Given atom percent
      if (mat % atom_density(1) > ZERO) then
        atom_fraction(i) = mat % atom_density(i)
        mass_fraction(i) = mat % atom_density(i) * awr

      ! Given weight percent
      else
        atom_fraction(i) = -mat % atom_density(i) / awr
        mass_fraction(i) = -mat % atom_density(i)
      end if

      Z_eq_sq = Z_eq_sq + atom_fraction(i) * nuclides(mat % nuclide(i)) % Z**2
    end do
    atom_sum = sum(atom_fraction)
    mass_sum = sum(mass_fraction)

    ! Given material density in g/cm^3
    if (mat % density < ZERO) then
      density = -mat % density * (atom_sum / mass_sum) * N_AVOGADRO / MASS_NEUTRON
      density_gpcc = -mat % density

    ! Given material density in atom/b-cm
    else
      density = mat % density
      density_gpcc = mat % density * (mass_sum / atom_sum) * MASS_NEUTRON / &
           N_AVOGADRO
    end if

    Z_eq_sq = Z_eq_sq / atom_sum
    atom_fraction = atom_fraction / atom_sum
    mass_fraction = mass_fraction / mass_sum

    ! Calculate the molecular DCS and the molecular total stopping power using
    ! Bragg's additivity rule. Note: the collision stopping power cannot be
    ! accurately calculated using Bragg's additivity rule since the mean
    ! excitation energies and the density effect corrections cannot simply be
    ! summed together. Bragg's additivity rule fails especially when a
    ! higher-density compound is composed of elements that are in lower-density
    ! form at normal temperature and pressure (at which the NIST stopping
    ! powers are given). It will be used to approximate the collision stopping
    ! powers for now, but should be fixed in the future.
    do i = 1, mat % n_nuclides
      ! Get pointer to current element
      elm => elements(mat % element(i))

      ! TODO: for molecular DCS, atom_fraction should actually be the number of
      ! atoms in the molecule.
      ! Accumulate material DCS
      this % dcs = this % dcs + atom_fraction(i) * elm % Z**2 / Z_eq_sq * elm % dcs

      ! Accumulate material total stopping power
      stopping_power = stopping_power + mass_fraction(i) * density_gpcc * &
           (elm % stopping_power_collision + elm % stopping_power_radiative)
    end do

    ! Calculate inverse bremsstrahlung mean free path
    do i = 1, n_e
      e = ttb_e_grid(i)
      if (e <= energy_cutoff(PHOTON)) cycle

      ! Ratio of the velocity of the charged particle to the speed of light
      beta = sqrt(e*(e + TWO*MASS_ELECTRON)) / (e + MASS_ELECTRON)

      ! Integration lower bound
      k_c = energy_cutoff(PHOTON) / e

      ! Find the upper bounding index of the reduced photon cutoff energy
      i_k = binary_search(ttb_k_grid, n_k, k_c) + 1

      ! Get the interpolation bounds
      k_l = ttb_k_grid(i_k-1)
      k_r = ttb_k_grid(i_k)
      x_l = this % dcs(i_k-1, i)
      x_r = this % dcs(i_k, i)

      ! Use linear interpolation in reduced photon energy k to find value of
      ! the DCS at the cutoff energy
      x_c = (x_l * (k_r - k_c) + x_r * (k_c - k_l)) / (k_r - k_l)

      ! Calculate the CDF using the trapezoidal rule in log-log space
      c = HALF * (log(k_r) - log(k_c)) * (x_c + x_r)
      this % cdf(i_k,i) = c
      do j = i_k, n_k - 1
        c = c + HALF * (log(ttb_k_grid(j+1)) - log(ttb_k_grid(j))) * &
             (this % dcs(j,i) + this % dcs(j+1,i))
        this % cdf(j+1,i) = c
      end do

      ! Calculate the inverse bremsstrahlung mean free path
      mfp_inv(i) = c * density * Z_eq_sq / beta**2 * 1.0e-3_8
    end do

    ! Calculate photon number yield
    mfp_inv(:) = mfp_inv(:) / stopping_power(:)
    call spline(ttb_e_grid, mfp_inv, z, n_e)
    do i = 1, n_e
      this % yield(i) = spline_integrate(ttb_e_grid, mfp_inv, z, n_e, &
           energy_cutoff(PHOTON), ttb_e_grid(i))
    end do

    ! Use logarithm of number yield since it is log-log interpolated
    where (this % yield > ZERO)
      this % yield = log(this % yield)
    end where

    deallocate(atom_fraction, mass_fraction, stopping_power, mfp_inv, z)

  end subroutine bremsstrahlung_init

end module material_header
