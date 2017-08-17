module tally_filter

  use hdf5, only: HID_T

  use algorithm,           only: binary_search
  use constants,           only: ONE, NO_BIN_FOUND, FP_PRECISION, ERROR_REAL
  use dict_header,         only: DictIntInt
  use error
  use geometry_header
  use hdf5_interface
  use particle_header,     only: Particle
  use surface_header
  use string,              only: to_str, to_f_string
  use tally_filter_header

  ! Inherit other filters
  use tally_filter_distribcell
  use tally_filter_energy
  use tally_filter_material
  use tally_filter_mesh

  implicit none

!===============================================================================
! UNIVERSEFILTER specifies which geometric universes tally events reside in.
!===============================================================================
  type, extends(TallyFilter) :: UniverseFilter
    integer, allocatable :: universes(:)
    type(DictIntInt)     :: map
  contains
    procedure :: get_all_bins => get_all_bins_universe
    procedure :: to_statepoint => to_statepoint_universe
    procedure :: text_label => text_label_universe
    procedure :: initialize => initialize_universe
  end type UniverseFilter

!===============================================================================
! CELLFILTER specifies which geometric cells tally events reside in.
!===============================================================================
  type, extends(TallyFilter) :: CellFilter
    integer, allocatable :: cells(:)
    type(DictIntInt)     :: map
  contains
    procedure :: get_all_bins => get_all_bins_cell
    procedure :: to_statepoint => to_statepoint_cell
    procedure :: text_label => text_label_cell
    procedure :: initialize => initialize_cell
  end type CellFilter

!===============================================================================
! CELLFROMFILTER specifies which geometric cells particles exit when crossing a
! surface.
!===============================================================================
  type, extends(CellFilter) :: CellFromFilter
  contains
    procedure :: get_all_bins => get_all_bins_cell_from
    procedure :: to_statepoint => to_statepoint_cell_from
    procedure :: text_label => text_label_cell_from
  end type CellFromFilter

!===============================================================================
! CELLBORNFILTER specifies which cell the particle was born in.
!===============================================================================
  type, extends(TallyFilter) :: CellbornFilter
    integer, allocatable :: cells(:)
    type(DictIntInt)     :: map
  contains
    procedure :: get_all_bins => get_all_bins_cellborn
    procedure :: to_statepoint => to_statepoint_cellborn
    procedure :: text_label => text_label_cellborn
    procedure :: initialize => initialize_cellborn
  end type CellbornFilter

!===============================================================================
! SURFACEFILTER specifies which surface particles are crossing
!===============================================================================
  type, extends(TallyFilter) :: SurfaceFilter
    integer, allocatable :: surfaces(:)

    ! True if this filter is used for surface currents
    logical              :: current = .false.
  contains
    procedure :: get_all_bins => get_all_bins_surface
    procedure :: to_statepoint => to_statepoint_surface
    procedure :: text_label => text_label_surface
    procedure :: initialize => initialize_surface
  end type SurfaceFilter

!===============================================================================
! DELAYEDGROUPFILTER bins outgoing fission neutrons in their delayed groups.
! The get_all_bins functionality is not actually used.  The bins are manually
! iterated over in the scoring subroutines.
!===============================================================================
  type, extends(TallyFilter) :: DelayedGroupFilter
    integer, allocatable :: groups(:)
  contains
    procedure :: get_all_bins => get_all_bins_dg
    procedure :: to_statepoint => to_statepoint_dg
    procedure :: text_label => text_label_dg
  end type DelayedGroupFilter

!===============================================================================
! MUFILTER bins the incoming-outgoing direction cosine.  This is only used for
! scatter reactions.
!===============================================================================
  type, extends(TallyFilter) :: MuFilter
    real(8), allocatable :: bins(:)
  contains
    procedure :: get_all_bins => get_all_bins_mu
    procedure :: to_statepoint => to_statepoint_mu
    procedure :: text_label => text_label_mu
  end type MuFilter

!===============================================================================
! POLARFILTER bins the incident neutron polar angle (relative to the global
! z-axis).
!===============================================================================
  type, extends(TallyFilter) :: PolarFilter
    real(8), allocatable :: bins(:)
  contains
    procedure :: get_all_bins => get_all_bins_polar
    procedure :: to_statepoint => to_statepoint_polar
    procedure :: text_label => text_label_polar
  end type PolarFilter

!===============================================================================
! AZIMUTHALFILTER bins the incident neutron azimuthal angle (relative to the
! global xy-plane).
!===============================================================================
  type, extends(TallyFilter) :: AzimuthalFilter
    real(8), allocatable :: bins(:)
  contains
    procedure :: get_all_bins => get_all_bins_azimuthal
    procedure :: to_statepoint => to_statepoint_azimuthal
    procedure :: text_label => text_label_azimuthal
  end type AzimuthalFilter

!===============================================================================
! EnergyFunctionFilter multiplies tally scores by an arbitrary function of
! incident energy described by a piecewise linear-linear interpolation.
!===============================================================================
  type, extends(TallyFilter) :: EnergyFunctionFilter
    real(8), allocatable :: energy(:)
    real(8), allocatable :: y(:)

  contains
    procedure :: get_all_bins => get_all_bins_energyfunction
    procedure :: to_statepoint => to_statepoint_energyfunction
    procedure :: text_label => text_label_energyfunction
  end type EnergyFunctionFilter

contains

!===============================================================================
! METHODS: for a description of these methods, see their counterparts bound to
! the abstract TallyFilter class.
!===============================================================================

!===============================================================================
! UniverseFilter methods
!===============================================================================
  subroutine get_all_bins_universe(this, p, estimator, match)
    class(UniverseFilter), intent(in)  :: this
    type(Particle),        intent(in)  :: p
    integer,               intent(in)  :: estimator
    type(TallyFilterMatch),     intent(inout) :: match

    integer :: i

    ! Iterate over coordinate levels to see which universes match
    do i = 1, p % n_coord
      if (this % map % has_key(p % coord(i) % universe)) then
        call match % bins % push_back(this % map % get_key(p % coord(i) &
             % universe))
        call match % weights % push_back(ONE)
      end if
    end do

  end subroutine get_all_bins_universe

  subroutine to_statepoint_universe(this, filter_group)
    class(UniverseFilter), intent(in) :: this
    integer(HID_T),        intent(in) :: filter_group

    integer :: i
    integer, allocatable :: universe_ids(:)

    call write_dataset(filter_group, "type", "universe")
    call write_dataset(filter_group, "n_bins", this % n_bins)

    allocate(universe_ids(size(this % universes)))
    do i = 1, size(this % universes)
      universe_ids(i) = universes(this % universes(i)) % id
    end do
    call write_dataset(filter_group, "bins", universe_ids)
  end subroutine to_statepoint_universe

  subroutine initialize_universe(this)
    class(UniverseFilter), intent(inout) :: this

    integer :: i, id

    ! Convert ids to indices.
    do i = 1, this % n_bins
      id = this % universes(i)
      if (universe_dict % has_key(id)) then
        this % universes(i) = universe_dict % get_key(id)
      else
        call fatal_error("Could not find universe " // trim(to_str(id)) &
             &// " specified on a tally filter.")
      end if
    end do

    ! Generate mapping from universe indices to filter bins.
    do i = 1, this % n_bins
      call this % map % add_key(this % universes(i), i)
    end do
  end subroutine initialize_universe

  function text_label_universe(this, bin) result(label)
    class(UniverseFilter), intent(in) :: this
    integer,               intent(in) :: bin
    character(MAX_LINE_LEN)           :: label

    label = "Universe " // to_str(universes(this % universes(bin)) % id)
  end function text_label_universe

!===============================================================================
! CellFilter methods
!===============================================================================
  subroutine get_all_bins_cell(this, p, estimator, match)
    class(CellFilter), intent(in)  :: this
    type(Particle),    intent(in)  :: p
    integer,           intent(in)  :: estimator
    type(TallyFilterMatch), intent(inout) :: match

    integer :: i

    ! Iterate over coordinate levels to see with cells match
    do i = 1, p % n_coord
      if (this % map % has_key(p % coord(i) % cell)) then
        call match % bins % push_back(this % map % get_key(p % coord(i) % cell))
        call match % weights % push_back(ONE)
      end if
    end do

  end subroutine get_all_bins_cell

  subroutine to_statepoint_cell(this, filter_group)
    class(CellFilter), intent(in) :: this
    integer(HID_T),    intent(in) :: filter_group

    integer :: i
    integer, allocatable :: cell_ids(:)

    call write_dataset(filter_group, "type", "cell")
    call write_dataset(filter_group, "n_bins", this % n_bins)

    allocate(cell_ids(size(this % cells)))
    do i = 1, size(this % cells)
      cell_ids(i) = cells(this % cells(i)) % id
    end do
    call write_dataset(filter_group, "bins", cell_ids)
  end subroutine to_statepoint_cell

  subroutine initialize_cell(this)
    class(CellFilter), intent(inout) :: this

    integer :: i, id

    ! Convert ids to indices.
    do i = 1, this % n_bins
      id = this % cells(i)
      if (cell_dict % has_key(id)) then
        this % cells(i) = cell_dict % get_key(id)
      else
        call fatal_error("Could not find cell " // trim(to_str(id)) &
             &// " specified on tally filter.")
      end if
    end do

    ! Generate mapping from cell indices to filter bins.
    do i = 1, this % n_bins
      call this % map % add_key(this % cells(i), i)
    end do
  end subroutine initialize_cell

  function text_label_cell(this, bin) result(label)
    class(CellFilter), intent(in) :: this
    integer,           intent(in) :: bin
    character(MAX_LINE_LEN)       :: label

    label = "Cell " // to_str(cells(this % cells(bin)) % id)
  end function text_label_cell

!===============================================================================
! CellFromFilter methods
!===============================================================================
  subroutine get_all_bins_cell_from(this, p, estimator, match)
    class(CellFromFilter), intent(in)  :: this
    type(Particle),    intent(in)  :: p
    integer,           intent(in)  :: estimator
    type(TallyFilterMatch), intent(inout) :: match

    integer :: i

    ! Starting one coordinate level deeper, find the next bin.
    do i = 1, p % last_n_coord
      if (this % map % has_key(p % last_cell(i))) then
        call match % bins % push_back(this % map % get_key(p % last_cell(i)))
        call match % weights % push_back(ONE)
        exit
      end if
    end do

  end subroutine get_all_bins_cell_from

  subroutine to_statepoint_cell_from(this, filter_group)
    class(CellFromFilter), intent(in) :: this
    integer(HID_T),        intent(in) :: filter_group

    integer :: i
    integer, allocatable :: cell_ids(:)

    call write_dataset(filter_group, "type", "cellfrom")
    call write_dataset(filter_group, "n_bins", this % n_bins)

    allocate(cell_ids(size(this % cells)))
    do i = 1, size(this % cells)
      cell_ids(i) = cells(this % cells(i)) % id
    end do
    call write_dataset(filter_group, "bins", cell_ids)
  end subroutine to_statepoint_cell_from

  function text_label_cell_from(this, bin) result(label)
    class(CellFromFilter), intent(in) :: this
    integer,               intent(in) :: bin
    character(MAX_LINE_LEN)           :: label

    label = "Cell from " // to_str(cells(this % cells(bin)) % id)
  end function text_label_cell_from

!===============================================================================
! CellbornFilter methods
!===============================================================================
  subroutine get_all_bins_cellborn(this, p, estimator, match)
    class(CellbornFilter), intent(in)  :: this
    type(Particle),        intent(in)  :: p
    integer,               intent(in)  :: estimator
    type(TallyFilterMatch),     intent(inout) :: match

      if (this % map % has_key(p % cell_born)) then
        call match % bins % push_back(this % map % get_key(p % cell_born))
        call match % weights % push_back(ONE)
      end if

  end subroutine get_all_bins_cellborn

  subroutine to_statepoint_cellborn(this, filter_group)
    class(CellbornFilter), intent(in) :: this
    integer(HID_T),        intent(in) :: filter_group

    integer :: i
    integer, allocatable :: cell_ids(:)

    call write_dataset(filter_group, "type", "cellborn")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    allocate(cell_ids(size(this % cells)))
    do i = 1, size(this % cells)
      cell_ids(i) = cells(this % cells(i)) % id
    end do
    call write_dataset(filter_group, "bins", cell_ids)
  end subroutine to_statepoint_cellborn

  subroutine initialize_cellborn(this)
    class(CellbornFilter), intent(inout) :: this

    integer :: i, id

    ! Convert ids to indices.
    do i = 1, this % n_bins
      id = this % cells(i)
      if (cell_dict % has_key(id)) then
        this % cells(i) = cell_dict % get_key(id)
      else
        call fatal_error("Could not find cell " // trim(to_str(id)) &
             &// " specified on tally filter.")
      end if
    end do

    ! Generate mapping from cell indices to filter bins.
    do i = 1, this % n_bins
      call this % map % add_key(this % cells(i), i)
    end do
  end subroutine initialize_cellborn

  function text_label_cellborn(this, bin) result(label)
    class(CellbornFilter), intent(in) :: this
    integer,               intent(in) :: bin
    character(MAX_LINE_LEN)           :: label

    label = "Birth Cell " // to_str(cells(this % cells(bin)) % id)
  end function text_label_cellborn

!===============================================================================
! SurfaceFilter methods
!===============================================================================
  subroutine get_all_bins_surface(this, p, estimator, match)
    class(SurfaceFilter), intent(in)  :: this
    type(Particle),       intent(in)  :: p
    integer,              intent(in)  :: estimator
    type(TallyFilterMatch),    intent(inout) :: match

    integer :: i

      do i = 1, this % n_bins
        if (abs(p % surface) == this % surfaces(i)) then
          call match % bins % push_back(i)
          if (p % surface < 0) then
            call match % weights % push_back(-ONE)
          else
            call match % weights % push_back(ONE)
          end if
          exit
        end if
      end do

  end subroutine get_all_bins_surface

  subroutine to_statepoint_surface(this, filter_group)
    class(SurfaceFilter), intent(in) :: this
    integer(HID_T),       intent(in) :: filter_group

    call write_dataset(filter_group, "type", "surface")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % surfaces)
  end subroutine to_statepoint_surface

  subroutine initialize_surface(this)
    class(SurfaceFilter), intent(inout) :: this

    integer :: i, id

    ! Convert ids to indices.
    do i = 1, this % n_bins
      id = this % surfaces(i)
      if (surface_dict % has_key(id)) then
        this % surfaces(i) = surface_dict % get_key(id)
      else
        call fatal_error("Could not find surface " // trim(to_str(id)) &
             &// " specified on tally filter.")
      end if
    end do
  end subroutine initialize_surface

  function text_label_surface(this, bin) result(label)
    class(SurfaceFilter), intent(in) :: this
    integer,              intent(in) :: bin
    character(MAX_LINE_LEN)          :: label

    label = "Surface " // to_str(surfaces(this % surfaces(bin)) % obj % id)
  end function text_label_surface

!===============================================================================
! DelayedGroupFilter methods
!===============================================================================
  subroutine get_all_bins_dg(this, p, estimator, match)
    class(DelayedGroupFilter), intent(in)  :: this
    type(Particle),            intent(in)  :: p
    integer,                   intent(in)  :: estimator
    type(TallyFilterMatch),         intent(inout) :: match

    call match % bins % push_back(1)
    call match % weights % push_back(ONE)
  end subroutine get_all_bins_dg

  subroutine to_statepoint_dg(this, filter_group)
    class(DelayedGroupFilter), intent(in) :: this
    integer(HID_T),            intent(in) :: filter_group

    call write_dataset(filter_group, "type", "delayedgroup")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % groups)
  end subroutine to_statepoint_dg

  function text_label_dg(this, bin) result(label)
    class(DelayedGroupFilter), intent(in) :: this
    integer,                   intent(in) :: bin
    character(MAX_LINE_LEN)               :: label

    label = "Delayed Group " // to_str(this % groups(bin))
  end function text_label_dg

!===============================================================================
! MuFilter methods
!===============================================================================
  subroutine get_all_bins_mu(this, p, estimator, match)
    class(MuFilter), intent(in)  :: this
    type(Particle),  intent(in)  :: p
    integer,         intent(in)  :: estimator
    type(TallyFilterMatch), intent(inout) :: match

    integer :: n
    integer :: bin

    n = this % n_bins

    ! Search to find incoming energy bin.
    bin = binary_search(this % bins, n + 1, p % mu)
    if (bin /= NO_BIN_FOUND) then
      call match % bins % push_back(bin)
      call match % weights % push_back(ONE)
    end if
  end subroutine get_all_bins_mu

  subroutine to_statepoint_mu(this, filter_group)
    class(MuFilter), intent(in) :: this
    integer(HID_T),  intent(in) :: filter_group

    call write_dataset(filter_group, "type", "mu")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % bins)
  end subroutine to_statepoint_mu

  function text_label_mu(this, bin) result(label)
    class(MuFilter), intent(in) :: this
    integer,         intent(in) :: bin
    character(MAX_LINE_LEN)     :: label

    real(8) :: E0, E1

    E0 = this % bins(bin)
    E1 = this % bins(bin + 1)
    label = "Change-in-Angle [" // trim(to_str(E0)) // ", " &
         // trim(to_str(E1)) // ")"
  end function text_label_mu

!===============================================================================
! PolarFilter methods
!===============================================================================
  subroutine get_all_bins_polar(this, p, estimator, match)
    class(PolarFilter), intent(in)  :: this
    type(Particle),     intent(in)  :: p
    integer,            intent(in)  :: estimator
    type(TallyFilterMatch),  intent(inout) :: match

    integer :: n
    integer :: bin
    real(8) :: theta

    n = this % n_bins

    ! Make sure the correct direction vector is used.
    if (estimator == ESTIMATOR_TRACKLENGTH) then
      theta = acos(p % coord(1) % uvw(3))
    else
      theta = acos(p % last_uvw(3))
    end if

    ! Search to find polar angle bin.
    bin = binary_search(this % bins, n + 1, theta)
    if (bin /= NO_BIN_FOUND) then
      call match % bins % push_back(bin)
      call match % weights % push_back(ONE)
    end if
  end subroutine get_all_bins_polar

  subroutine to_statepoint_polar(this, filter_group)
    class(PolarFilter), intent(in) :: this
    integer(HID_T),     intent(in) :: filter_group

    call write_dataset(filter_group, "type", "polar")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % bins)
  end subroutine to_statepoint_polar

  function text_label_polar(this, bin) result(label)
    class(PolarFilter), intent(in) :: this
    integer,            intent(in) :: bin
    character(MAX_LINE_LEN)        :: label

    real(8) :: E0, E1

    E0 = this % bins(bin)
    E1 = this % bins(bin + 1)
    label = "Polar Angle [" // trim(to_str(E0)) // ", " // trim(to_str(E1)) &
         // ")"
  end function text_label_polar

!===============================================================================
! AzimuthalFilter methods
!===============================================================================
  subroutine get_all_bins_azimuthal(this, p, estimator, match)
    class(AzimuthalFilter), intent(in)  :: this
    type(Particle),         intent(in)  :: p
    integer,                intent(in)  :: estimator
    type(TallyFilterMatch),      intent(inout) :: match

    integer :: n
    integer :: bin
    real(8) :: phi

    n = this % n_bins

    ! Make sure the correct direction vector is used.
    if (estimator == ESTIMATOR_TRACKLENGTH) then
      phi = atan2(p % coord(1) % uvw(2), p % coord(1) % uvw(1))
    else
      phi = atan2(p % last_uvw(2), p % last_uvw(1))
    end if

    ! Search to find azimuthal angle bin.
    bin = binary_search(this % bins, n + 1, phi)
    if (bin /= NO_BIN_FOUND) then
      call match % bins % push_back(bin)
      call match % weights % push_back(ONE)
    end if

  end subroutine get_all_bins_azimuthal

  subroutine to_statepoint_azimuthal(this, filter_group)
    class(AzimuthalFilter), intent(in) :: this
    integer(HID_T),         intent(in) :: filter_group

    call write_dataset(filter_group, "type", "azimuthal")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % bins)
  end subroutine to_statepoint_azimuthal

  function text_label_azimuthal(this, bin) result(label)
    class(AzimuthalFilter), intent(in) :: this
    integer,                intent(in) :: bin
    character(MAX_LINE_LEN)            :: label

    real(8) :: E0, E1

    E0 = this % bins(bin)
    E1 = this % bins(bin + 1)
    label = "Azimuthal Angle [" // trim(to_str(E0)) // ", " &
         // trim(to_str(E1)) // ")"
  end function text_label_azimuthal

!===============================================================================
! EnergyFunctionFilter methods
!===============================================================================
  subroutine get_all_bins_energyfunction(this, p, estimator, match)
    class(EnergyFunctionFilter), intent(in)  :: this
    type(Particle),              intent(in)  :: p
    integer,                     intent(in)  :: estimator
    type(TallyFilterMatch),      intent(inout) :: match

    integer :: n, indx
    real(8) :: E, f, weight

    select type(this)
    type is (EnergyFunctionFilter)
      n = size(this % energy)

      ! Get pre-collision energy of particle
      E = p % last_E

      ! Search to find incoming energy bin.
      indx = binary_search(this % energy, n, E)

      ! Compute an interpolation factor between nearest bins.
      f = (E - this % energy(indx)) &
           / (this % energy(indx+1) - this % energy(indx))

      ! Interpolate on the lin-lin grid.
      call match % bins % push_back(1)
      weight = (ONE - f) * this % y(indx) + f * this % y(indx+1)
      call match % weights % push_back(weight)
    end select
  end subroutine get_all_bins_energyfunction

  subroutine to_statepoint_energyfunction(this, filter_group)
    class(EnergyFunctionFilter), intent(in) :: this
    integer(HID_T),              intent(in) :: filter_group

    select type(this)
    type is (EnergyFunctionFilter)
      call write_dataset(filter_group, "type", "energyfunction")
      call write_dataset(filter_group, "energy", this % energy)
      call write_dataset(filter_group, "y", this % y)
    end select
  end subroutine to_statepoint_energyfunction

  function text_label_energyfunction(this, bin) result(label)
    class(EnergyFunctionFilter), intent(in) :: this
    integer,                     intent(in) :: bin
    character(MAX_LINE_LEN)                 :: label

    select type(this)
    type is (EnergyFunctionFilter)
      write(label, FMT="(A, ES8.1, A, ES8.1, A, ES8.1, A, ES8.1, A)") &
           "Energy Function f([", this % energy(1), ", ..., ", &
           this % energy(size(this % energy)), "]) = [", this % y(1), &
           ", ..., ", this % y(size(this % y)), "]"
    end select
  end function text_label_energyfunction

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_filter_get_type(index, type) result(err) bind(C)
    ! Get the type of a filter
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR), intent(out) :: type(*)
    integer(C_INT) :: err

    integer :: i
    character(20) :: type_

    if (index >= 1 .and. index <= n_filters) then
      if (allocated(filters(index) % obj)) then
        ! Get type as a Fortran string
        select type (f => filters(index) % obj)
        type is (AzimuthalFilter)
          type_ = 'azimuthal'
        type is (CellFilter)
          type_ = 'cell'
        type is (CellbornFilter)
          type_ = 'cellborn'
        type is (CellfromFilter)
          type_ = 'cellfrom'
        type is (DelayedGroupFilter)
          type_ = 'delayedgroup'
        type is (DistribcellFilter)
          type_ = 'distribcell'
        type is (EnergyFilter)
          type_ = 'energy'
        type is (EnergyoutFilter)
          type_ = 'energyout'
        type is (EnergyFunctionFilter)
          type_ = 'energyfunction'
        type is (MaterialFilter)
          type_ = 'material'
        type is (MeshFilter)
          type_ = 'mesh'
        type is (MuFilter)
          type_ = 'mu'
        type is (PolarFilter)
          type_ = 'polar'
        type is (SurfaceFilter)
          type_ = 'surface'
        type is (UniverseFilter)
          type_ = 'universe'
        end select

        ! Convert Fortran string to null-terminated C string. We assume the
        ! caller has allocated a char array buffer
        do i = 1, len_trim(type_)
          type(i) = type_(i:i)
        end do
        type(len_trim(type_) + 1) = C_NULL_CHAR

        err = 0
      else
        err = E_FILTER_NOT_ALLOCATED
      end if
    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_filter_get_type


  function openmc_filter_set_type(index, type) result(err) bind(C)
    ! Set the type of a filter
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR), intent(in) :: type(*)
    integer(C_INT) :: err

    character(:), allocatable :: type_

    ! Convert C string to Fortran string
    type_ = to_f_string(type)

    err = 0
    if (index >= 1 .and. index <= n_filters) then
      if (allocated(filters(index) % obj)) then
        err = E_ALREADY_ALLOCATED
      else
        select case (type_)
        case ('azimuthal')
          allocate(AzimuthalFilter :: filters(index) % obj)
        case ('cell')
          allocate(CellFilter :: filters(index) % obj)
        case ('cellborn')
          allocate(CellbornFilter :: filters(index) % obj)
        case ('cellfrom')
          allocate(CellfromFilter :: filters(index) % obj)
        case ('delayedgroup')
          allocate(DelayedGroupFilter :: filters(index) % obj)
        case ('distribcell')
          allocate(DistribcellFilter :: filters(index) % obj)
        case ('energy')
          allocate(EnergyFilter :: filters(index) % obj)
        case ('energyout')
          allocate(EnergyoutFilter :: filters(index) % obj)
        case ('energyfunction')
          allocate(EnergyFunctionFilter :: filters(index) % obj)
        case ('material')
          allocate(MaterialFilter :: filters(index) % obj)
        case ('mesh')
          allocate(MeshFilter :: filters(index) % obj)
        case ('mu')
          allocate(MuFilter :: filters(index) % obj)
        case ('polar')
          allocate(PolarFilter :: filters(index) % obj)
        case ('surface')
          allocate(SurfaceFilter :: filters(index) % obj)
        case ('universe')
          allocate(UniverseFilter :: filters(index) % obj)
        case default
          err = E_UNASSIGNED
        end select
      end if
    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_filter_set_type

end module tally_filter
