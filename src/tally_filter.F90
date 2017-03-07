module tally_filter

  use algorithm,           only: binary_search
  use constants,           only: ONE, NO_BIN_FOUND, FP_PRECISION, ERROR_REAL
  use dict_header,         only: DictIntInt
  use geometry_header,     only: root_universe, RectLattice, HexLattice
  use global
  use hdf5_interface
  use mesh_header,         only: RegularMesh
  use mesh,                only: get_mesh_bin, bin_to_mesh_indices, &
                                 get_mesh_indices, mesh_indices_to_bin, &
                                 mesh_intersects_1d, mesh_intersects_2d, &
                                 mesh_intersects_3d
  use particle_header,     only: Particle
  use string,              only: to_str
  use tally_filter_header, only: TallyFilter, TallyFilterContainer

  use hdf5, only: HID_T

  implicit none

!===============================================================================
! MESHFILTER indexes the location of particle events to a regular mesh.  For
! tracklength tallies, it will produce multiple valid bins and the bin weight
! will correspond to the fraction of the track length that lies in that bin.
!===============================================================================
  type, extends(TallyFilter) :: MeshFilter
    integer :: mesh
  contains
    procedure :: get_next_bin => get_next_bin_mesh
    procedure :: to_statepoint => to_statepoint_mesh
    procedure :: text_label => text_label_mesh
  end type MeshFilter

!===============================================================================
! UNIVERSEFILTER specifies which geometric universes tally events reside in.
!===============================================================================
  type, extends(TallyFilter) :: UniverseFilter
    integer, allocatable :: universes(:)
    type(DictIntInt)     :: map
  contains
    procedure :: get_next_bin => get_next_bin_universe
    procedure :: to_statepoint => to_statepoint_universe
    procedure :: text_label => text_label_universe
    procedure :: initialize => initialize_universe
  end type UniverseFilter

!===============================================================================
! MATERIAL specifies which material tally events reside in.
!===============================================================================
  type, extends(TallyFilter) :: MaterialFilter
    integer, allocatable :: materials(:)
    type(DictIntInt)     :: map
  contains
    procedure :: get_next_bin => get_next_bin_material
    procedure :: to_statepoint => to_statepoint_material
    procedure :: text_label => text_label_material
    procedure :: initialize => initialize_material
  end type MaterialFilter

!===============================================================================
! CELLFILTER specifies which geometric cells tally events reside in.
!===============================================================================
  type, extends(TallyFilter) :: CellFilter
    integer, allocatable :: cells(:)
    type(DictIntInt)     :: map
  contains
    procedure :: get_next_bin => get_next_bin_cell
    procedure :: to_statepoint => to_statepoint_cell
    procedure :: text_label => text_label_cell
    procedure :: initialize => initialize_cell
  end type CellFilter

!===============================================================================
! DISTRIBCELLFILTER specifies which distributed geometric cells tally events
! reside in.
!===============================================================================
  type, extends(TallyFilter) :: DistribcellFilter
    integer :: cell
  contains
    procedure :: get_next_bin => get_next_bin_distribcell
    procedure :: to_statepoint => to_statepoint_distribcell
    procedure :: text_label => text_label_distribcell
    procedure :: initialize => initialize_distribcell
  end type DistribcellFilter

!===============================================================================
! CELLBORNFILTER specifies which cell the particle was born in.
!===============================================================================
  type, extends(TallyFilter) :: CellbornFilter
    integer, allocatable :: cells(:)
    type(DictIntInt)     :: map
  contains
    procedure :: get_next_bin => get_next_bin_cellborn
    procedure :: to_statepoint => to_statepoint_cellborn
    procedure :: text_label => text_label_cellborn
    procedure :: initialize => initialize_cellborn
  end type CellbornFilter

!===============================================================================
! SURFACEFILTER is currently not implemented for usual geometric surfaces, but
! it is used as a placeholder for mesh surfaces used in current tallies.
!===============================================================================
  type, extends(TallyFilter) :: SurfaceFilter
    integer, allocatable :: surfaces(:)
  contains
    procedure :: get_next_bin => get_next_bin_surface
    procedure :: to_statepoint => to_statepoint_surface
    procedure :: text_label => text_label_surface
    procedure :: initialize => initialize_surface
  end type SurfaceFilter

!===============================================================================
! ENERGYFILTER bins the incident neutron energy.
!===============================================================================
  type, extends(TallyFilter) :: EnergyFilter
    real(8), allocatable :: bins(:)

    ! True if transport group number can be used directly to get bin number
    logical              :: matches_transport_groups = .false.

  contains
    procedure :: get_next_bin => get_next_bin_energy
    procedure :: to_statepoint => to_statepoint_energy
    procedure :: text_label => text_label_energy
  end type EnergyFilter

!===============================================================================
! ENERGYOUTFILTER bins the outgoing neutron energy.  Only scattering events use
! the get_next_bin functionality.  Nu-fission tallies manually iterate over the
! filter bins.
!===============================================================================
  type, extends(TallyFilter) :: EnergyoutFilter
    real(8), allocatable :: bins(:)

    ! True if transport group number can be used directly to get bin number
    logical              :: matches_transport_groups = .false.

  contains
    procedure :: get_next_bin => get_next_bin_energyout
    procedure :: to_statepoint => to_statepoint_energyout
    procedure :: text_label => text_label_energyout
  end type EnergyoutFilter

!===============================================================================
! DELAYEDGROUPFILTER bins outgoing fission neutrons in their delayed groups.
! The get_next_bin functionality is not actually used.  The bins are manually
! iterated over in the scoring subroutines.
!===============================================================================
  type, extends(TallyFilter) :: DelayedGroupFilter
    integer, allocatable :: groups(:)
  contains
    procedure :: get_next_bin => get_next_bin_dg
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
    procedure :: get_next_bin => get_next_bin_mu
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
    procedure :: get_next_bin => get_next_bin_polar
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
    procedure :: get_next_bin => get_next_bin_azimuthal
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
    procedure :: get_next_bin => get_next_bin_energyfunction
    procedure :: to_statepoint => to_statepoint_energyfunction
    procedure :: text_label => text_label_energyfunction
  end type EnergyFunctionFilter

contains

!===============================================================================
! METHODS: for a description of these methods, see their counterparts bound to
! the abstract TallyFilter class.
!===============================================================================

!===============================================================================
! MeshFilter methods
!===============================================================================
  subroutine get_next_bin_mesh(this, p, estimator, current_bin, next_bin, &
       weight)
    class(MeshFilter), intent(in)  :: this
    type(Particle),    intent(in)  :: p
    integer,           intent(in)  :: estimator
    integer, value,    intent(in)  :: current_bin
    integer,           intent(out) :: next_bin
    real(8),           intent(out) :: weight

    integer, parameter :: MAX_SEARCH_ITER = 100 ! Maximum number of times we can
                                                !  can loop while trying to find
                                                !  the first intersection.

    integer :: j                    ! loop index for direction
    integer :: ijk0(3)              ! indices of starting coordinates
    integer :: ijk1(3)              ! indices of ending coordinates
    integer :: search_iter          ! loop count for intersection search
    real(8) :: uvw(3)               ! cosine of angle of particle
    real(8) :: xyz0(3)              ! starting/intermediate coordinates
    real(8) :: xyz1(3)              ! ending coordinates of particle
    real(8) :: xyz_cross            ! coordinates of next boundary
    real(8) :: d(3)                 ! distance to each bounding surface
    real(8) :: total_distance       ! distance of entire particle track
    real(8) :: distance             ! distance traveled in mesh cell
    logical :: start_in_mesh        ! starting coordinates inside mesh?
    logical :: end_in_mesh          ! ending coordinates inside mesh?
    type(RegularMesh), pointer :: m

    weight = ERROR_REAL

    ! Get a pointer to the mesh.
    m => meshes(this % mesh)

    if (estimator /= ESTIMATOR_TRACKLENGTH) then
      ! If this is an analog or collision tally, then there can only be one
      ! valid mesh bin.
      if (current_bin == NO_BIN_FOUND) then
        call get_mesh_bin(m, p % coord(1) % xyz, next_bin)
        weight = ONE
      else
        next_bin = NO_BIN_FOUND
      end if

    else
      ! A track can span multiple mesh bins so we need to handle a lot of
      ! intersection logic for tracklength tallies.

      ! ========================================================================
      ! Determine if the track intersects the tally mesh.

      ! Copy the starting and ending coordinates of the particle.  Offset these
      ! just a bit for the purposes of determining if there was an intersection
      ! in case the mesh surfaces coincide with lattice/geometric surfaces which
      ! might produce finite-precision errors.
      xyz0 = p % last_xyz + TINY_BIT * p % coord(1) % uvw
      xyz1 = p % coord(1) % xyz - TINY_BIT * p % coord(1) % uvw

      ! Determine indices for starting and ending location.
      call get_mesh_indices(m, xyz0, ijk0(:m % n_dimension), start_in_mesh)
      call get_mesh_indices(m, xyz1, ijk1(:m % n_dimension), end_in_mesh)

      ! If this is the first iteration of the filter loop, check if the track
      ! intersects any part of the mesh.
      if (current_bin == NO_BIN_FOUND) then
        if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) then
          if (m % n_dimension == 1) then
            if (.not. mesh_intersects_1d(m, xyz0, xyz1)) then
              next_bin = NO_BIN_FOUND
              return
            end if
          else if (m % n_dimension == 2) then
            if (.not. mesh_intersects_2d(m, xyz0, xyz1)) then
              next_bin = NO_BIN_FOUND
              return
            end if
          else
            if (.not. mesh_intersects_3d(m, xyz0, xyz1)) then
              next_bin = NO_BIN_FOUND
              return
            end if
          end if
        end if
      end if

      ! ========================================================================
      ! Figure out which mesh cell to tally.

      ! Copy the un-modified coordinates the particle direction.
      xyz0 = p % last_xyz
      xyz1 = p % coord(1) % xyz
      uvw = p % coord(1) % uvw

      ! Compute the length of the entire track.
      total_distance = sqrt(sum((xyz1 - xyz0)**2))

      if (current_bin == NO_BIN_FOUND) then
        ! We are looking for the first valid mesh bin.  Check to see if the
        ! particle starts inside the mesh.
        if (any(ijk0(:m % n_dimension) < 1) &
             .or. any(ijk0(:m % n_dimension) > m % dimension)) then
          ! The particle does not start in the mesh.  Note that we nudged the
          ! start and end coordinates by a TINY_BIT each so we will have
          ! difficulty resolving tracks that are less than 2*TINY_BIT in length.
          ! If the track is that short, it is also insignificant so we can
          ! safely ignore it in the tallies.
          if (total_distance < 2*TINY_BIT) then
            next_bin = NO_BIN_FOUND
            return
          end if

          ! The particle does not start in the mesh so keep iterating the ijk0
          ! indices to cross the nearest mesh surface until we've found a valid
          ! bin.  MAX_SEARCH_ITER prevents an infinite loop.
          search_iter = 0
          do while (any(ijk0(:m % n_dimension) < 1) &
               .or. any(ijk0(:m % n_dimension) > m % dimension))
            if (search_iter == MAX_SEARCH_ITER) then
              call warning("Failed to find a mesh intersection on a tally mesh &
                   &filter.")
              next_bin = NO_BIN_FOUND
              return
            end if

            do j = 1, m % n_dimension
              if (abs(uvw(j)) < FP_PRECISION) then
                d(j) = INFINITY
              else if (uvw(j) > 0) then
                xyz_cross = m % lower_left(j) + ijk0(j) * m % width(j)
                d(j) = (xyz_cross - xyz0(j)) / uvw(j)
              else
                xyz_cross = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
                d(j) = (xyz_cross - xyz0(j)) / uvw(j)
              end if
            end do
            j = minloc(d(:m % n_dimension), 1)
            if (uvw(j) > ZERO) then
              ijk0(j) = ijk0(j) + 1
            else
              ijk0(j) = ijk0(j) - 1
            end if

            search_iter = search_iter + 1
          end do
          distance = d(j)
          xyz0 = xyz0 + distance * uvw

        end if

      else
        ! We have already scored some mesh bins for this track.  Pick up where
        ! we left off and find the next mesh cell that the particle enters.

        ! Get the indices to the last bin we scored.
        call bin_to_mesh_indices(m, current_bin, ijk0(:m % n_dimension))

        ! If the particle track ends in that bin, then we are done.
        if (all(ijk0(:m % n_dimension) == ijk1(:m % n_dimension))) then
          next_bin = NO_BIN_FOUND
          return
        end if

        ! Figure out which face of the previous mesh cell our track exits, i.e.
        ! the closest surface of that cell for which
        ! dot(p % uvw, face_normal) > 0.
        do j = 1, m % n_dimension
          if (abs(uvw(j)) < FP_PRECISION) then
            d(j) = INFINITY
          else if (uvw(j) > 0) then
            xyz_cross = m % lower_left(j) + ijk0(j) * m % width(j)
            d(j) = (xyz_cross - xyz0(j)) / uvw(j)
          else
            xyz_cross = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
            d(j) = (xyz_cross - xyz0(j)) / uvw(j)
          end if
        end do
        j = minloc(d(:m % n_dimension), 1)

        ! Translate the starting coordintes by the distance to that face. This
        ! should be the xyz that we computed the distance to in the last
        ! iteration of the filter loop.
        distance = d(j)
        xyz0 = xyz0 + distance * uvw

        ! Increment the indices into the next mesh cell.
        if (uvw(j) > ZERO) then
          ijk0(j) = ijk0(j) + 1
        else
          ijk0(j) = ijk0(j) - 1
        end if

        ! If the next indices are invalid, then the track has left the mesh and
        ! we are done.
        if (any(ijk0(:m % n_dimension) < 1) &
             .or. any(ijk0(:m % n_dimension) > m % dimension)) then
          next_bin = NO_BIN_FOUND
          return
        end if
      end if

      ! ========================================================================
      ! Compute the length of the track segment in the appropiate mesh cell and
      ! return.

      if (all(ijk0(:m % n_dimension) == ijk1(:m % n_dimension))) then
        ! The track ends in this cell.  Use the particle end location rather
        ! than the mesh surface.
        distance = sqrt(sum((xyz1 - xyz0)**2))
      else
        ! The track exits this cell.  Use the distance to the mesh surface.
        do j = 1, m % n_dimension
          if (abs(uvw(j)) < FP_PRECISION) then
            d(j) = INFINITY
          else if (uvw(j) > 0) then
            xyz_cross = m % lower_left(j) + ijk0(j) * m % width(j)
            d(j) = (xyz_cross - xyz0(j)) / uvw(j)
          else
            xyz_cross = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
            d(j) = (xyz_cross - xyz0(j)) / uvw(j)
          end if
        end do
        distance = minval(d(:m % n_dimension))
      end if

      ! Assign the next tally bin and the score.
      next_bin = mesh_indices_to_bin(m, ijk0(:m % n_dimension))
      weight = distance / total_distance
    endif
  end subroutine get_next_bin_mesh

  subroutine to_statepoint_mesh(this, filter_group)
    class(MeshFilter), intent(in) :: this
    integer(HID_T),    intent(in) :: filter_group

    call write_dataset(filter_group, "type", "mesh")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", meshes(this % mesh) % id)
  end subroutine to_statepoint_mesh

  function text_label_mesh(this, bin) result(label)
    class(MeshFilter), intent(in) :: this
    integer,           intent(in) :: bin
    character(MAX_LINE_LEN)       :: label

    integer, allocatable       :: ijk(:)
    type(RegularMesh), pointer :: m

    m => meshes(this % mesh)
    allocate(ijk(m % n_dimension))
    call bin_to_mesh_indices(m, bin, ijk)
    if (m % n_dimension == 1) then
      label = "Mesh Index (" // trim(to_str(ijk(1))) // ")"
    elseif (m % n_dimension == 2) then
      label = "Mesh Index (" // trim(to_str(ijk(1))) // ", " // &
           trim(to_str(ijk(2))) // ")"
    elseif (m % n_dimension == 3) then
      label = "Mesh Index (" // trim(to_str(ijk(1))) // ", " // &
           trim(to_str(ijk(2))) // ", " // trim(to_str(ijk(3))) // ")"
    end if
  end function text_label_mesh

!===============================================================================
! UniverseFilter methods
!===============================================================================
  subroutine get_next_bin_universe(this, p, estimator, current_bin, next_bin, &
       weight)
    class(UniverseFilter), intent(in)  :: this
    type(Particle),        intent(in)  :: p
    integer,               intent(in)  :: estimator
    integer, value,        intent(in)  :: current_bin
    integer,               intent(out) :: next_bin
    real(8),               intent(out) :: weight

    integer :: i, start

    ! Find the coordinate level of the last bin we found.
    if (current_bin == NO_BIN_FOUND) then
      start = 1
    else
      do i = 1, p % n_coord
        if (p % coord(i) % universe == this % universes(current_bin)) then
          start = i + 1
          exit
        end if
      end do
    end if

    ! Starting one coordinate level deeper, find the next bin.
    next_bin = NO_BIN_FOUND
    weight = ERROR_REAL
    do i = start, p % n_coord
      if (this % map % has_key(p % coord(i) % universe)) then
        next_bin = this % map % get_key(p % coord(i) % universe)
        weight = ONE
        exit
      end if
    end do
  end subroutine get_next_bin_universe

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
! MaterialFilter methods
!===============================================================================
  subroutine get_next_bin_material(this, p, estimator, current_bin, next_bin, &
       weight)
    class(MaterialFilter), intent(in)  :: this
    type(Particle),        intent(in)  :: p
    integer,               intent(in)  :: estimator
    integer, value,        intent(in)  :: current_bin
    integer,               intent(out) :: next_bin
    real(8),               intent(out) :: weight

    next_bin = NO_BIN_FOUND
    weight = ERROR_REAL
    if (current_bin == NO_BIN_FOUND) then
      if (this % map % has_key(p % material)) then
        next_bin = this % map % get_key(p % material)
        weight = ONE
      end if
    end if
  end subroutine get_next_bin_material

  subroutine to_statepoint_material(this, filter_group)
    class(MaterialFilter), intent(in) :: this
    integer(HID_T),        intent(in) :: filter_group

    integer :: i
    integer, allocatable :: material_ids(:)

    call write_dataset(filter_group, "type", "material")
    call write_dataset(filter_group, "n_bins", this % n_bins)

    allocate(material_ids(size(this % materials)))
    do i = 1, size(this % materials)
      material_ids(i) = materials(this % materials(i)) % id
    end do
    call write_dataset(filter_group, "bins", material_ids)
  end subroutine to_statepoint_material

  subroutine initialize_material(this)
    class(MaterialFilter), intent(inout) :: this

    integer :: i, id

    ! Convert ids to indices.
    do i = 1, this % n_bins
      id = this % materials(i)
      if (material_dict % has_key(id)) then
        this % materials(i) = material_dict % get_key(id)
      else
        call fatal_error("Could not find material " // trim(to_str(id)) &
             &// " specified on a tally filter.")
      end if
    end do

    ! Generate mapping from material indices to filter bins.
    do i = 1, this % n_bins
      call this % map % add_key(this % materials(i), i)
    end do
  end subroutine initialize_material

  function text_label_material(this, bin) result(label)
    class(MaterialFilter), intent(in) :: this
    integer,               intent(in) :: bin
    character(MAX_LINE_LEN)           :: label

    label = "Material " // to_str(materials(this % materials(bin)) % id)
  end function text_label_material

!===============================================================================
! CellFilter methods
!===============================================================================
  subroutine get_next_bin_cell(this, p, estimator, current_bin, next_bin, &
       weight)
    class(CellFilter), intent(in)  :: this
    type(Particle),    intent(in)  :: p
    integer,           intent(in)  :: estimator
    integer, value,    intent(in)  :: current_bin
    integer,           intent(out) :: next_bin
    real(8),           intent(out) :: weight

    integer :: i, start

    ! Find the coordinate level of the last bin we found.
    if (current_bin == NO_BIN_FOUND) then
      start = 1
    else
      do i = 1, p % n_coord
        if (p % coord(i) % cell == this % cells(current_bin)) then
          start = i + 1
          exit
        end if
      end do
    end if

    ! Starting one coordinate level deeper, find the next bin.
    next_bin = NO_BIN_FOUND
    weight = ERROR_REAL
    do i = start, p % n_coord
      if (this % map % has_key(p % coord(i) % cell)) then
        next_bin = this % map % get_key(p % coord(i) % cell)
        weight = ONE
        exit
      end if
    end do
  end subroutine get_next_bin_cell

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
! DistribcellFilter methods
!===============================================================================
  subroutine get_next_bin_distribcell(this, p, estimator, current_bin, &
       next_bin, weight)
    class(DistribcellFilter), intent(in)  :: this
    type(Particle),           intent(in)  :: p
    integer,                  intent(in)  :: estimator
    integer, value,           intent(in)  :: current_bin
    integer,                  intent(out) :: next_bin
    real(8),                  intent(out) :: weight

    integer :: distribcell_index, offset, i

    if (current_bin == NO_BIN_FOUND) then
      distribcell_index = cells(this % cell) % distribcell_index
      offset = 0
      do i = 1, p % n_coord
        if (cells(p % coord(i) % cell) % type == FILL_UNIVERSE) then
          offset = offset + cells(p % coord(i) % cell) % &
               offset(distribcell_index)
        elseif (cells(p % coord(i) % cell) % type == FILL_LATTICE) then
          if (lattices(p % coord(i + 1) % lattice) % obj &
               % are_valid_indices([&
               p % coord(i + 1) % lattice_x, &
               p % coord(i + 1) % lattice_y, &
               p % coord(i + 1) % lattice_z])) then
            offset = offset + lattices(p % coord(i + 1) % lattice) % obj % &
                 offset(distribcell_index, &
                 p % coord(i + 1) % lattice_x, &
                 p % coord(i + 1) % lattice_y, &
                 p % coord(i + 1) % lattice_z)
          end if
        end if
        if (this % cell == p % coord(i) % cell) then
          next_bin = offset + 1
          weight = ONE
          return
        end if
      end do
    end if
    next_bin = NO_BIN_FOUND
    weight = ERROR_REAL
  end subroutine get_next_bin_distribcell

  subroutine to_statepoint_distribcell(this, filter_group)
    class(DistribcellFilter), intent(in) :: this
    integer(HID_T),           intent(in) :: filter_group

    call write_dataset(filter_group, "type", "distribcell")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", cells(this % cell) % id)
  end subroutine to_statepoint_distribcell

  subroutine initialize_distribcell(this)
    class(DistribcellFilter), intent(inout) :: this

    integer :: id

    ! Convert id to index.
    id = this % cell
    if (cell_dict % has_key(id)) then
      this % cell = cell_dict % get_key(id)
    else
      call fatal_error("Could not find cell " // trim(to_str(id)) &
           &// " specified on tally filter.")
    end if
  end subroutine initialize_distribcell

  function text_label_distribcell(this, bin) result(label)
    class(DistribcellFilter), intent(in) :: this
    integer,                  intent(in) :: bin
    character(MAX_LINE_LEN)              :: label

    integer                 :: offset
    type(Universe), pointer :: univ

    univ => universes(root_universe)
    offset = 0
    label = ''
    call find_offset(this % cell, univ, bin-1, offset, label)
    label = "Distributed Cell " // label
  end function text_label_distribcell

!===============================================================================
! CellbornFilter methods
!===============================================================================
  subroutine get_next_bin_cellborn(this, p, estimator, current_bin, next_bin, &
       weight)
    class(CellbornFilter), intent(in)  :: this
    type(Particle),        intent(in)  :: p
    integer,               intent(in)  :: estimator
    integer, value,        intent(in)  :: current_bin
    integer,               intent(out) :: next_bin
    real(8),               intent(out) :: weight

    next_bin = NO_BIN_FOUND
    weight = ERROR_REAL
    if (current_bin == NO_BIN_FOUND) then
      if (this % map % has_key(p % cell_born)) then
        next_bin = this % map % get_key(p % cell_born)
        weight = ONE
      end if
    end if
  end subroutine get_next_bin_cellborn

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
  subroutine get_next_bin_surface(this, p, estimator, current_bin, next_bin, &
       weight)
    class(SurfaceFilter), intent(in)  :: this
    type(Particle),       intent(in)  :: p
    integer,              intent(in)  :: estimator
    integer, value,       intent(in)  :: current_bin
    integer,              intent(out) :: next_bin
    real(8),              intent(out) :: weight

    integer :: i

    next_bin = NO_BIN_FOUND
    weight = ERROR_REAL
    if (current_bin == NO_BIN_FOUND) then
      do i = 1, this % n_bins
        if (p % surface == this % surfaces(i)) then
          next_bin = i
          weight = ONE
          exit
        end if
      end do
    end if
  end subroutine get_next_bin_surface

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
! EnergyFilter methods
!===============================================================================
  subroutine get_next_bin_energy(this, p, estimator, current_bin, next_bin, &
       weight)
    class(EnergyFilter), intent(in)  :: this
    type(Particle),      intent(in)  :: p
    integer,             intent(in)  :: estimator
    integer, value,      intent(in)  :: current_bin
    integer,             intent(out) :: next_bin
    real(8),             intent(out) :: weight

    integer :: n
    real(8) :: E

    if (current_bin == NO_BIN_FOUND) then
      n = this % n_bins

      if ((.not. run_CE) .and. this % matches_transport_groups) then
        if (estimator == ESTIMATOR_TRACKLENGTH) then
          next_bin = p % g
        else
          next_bin = p % last_g
        end if

        ! Tallies are ordered in increasing groups, group indices
        ! however are the opposite, so switch
        next_bin = num_energy_groups - next_bin + 1
        weight = ONE

      else
        ! Make sure the correct energy is used.
        if (estimator == ESTIMATOR_TRACKLENGTH) then
          E = p % E
        else
          E = p % last_E
        end if

        ! Check if energy of the particle is within energy bins.
        if (E < this % bins(1) .or. E > this % bins(n + 1)) then
          next_bin = NO_BIN_FOUND
          weight = ERROR_REAL
        else
          ! Search to find incoming energy bin.
          next_bin = binary_search(this % bins, n + 1, E)
          weight = ONE
        end if
      end if

    else
      next_bin = NO_BIN_FOUND
      weight = ERROR_REAL
    end if
  end subroutine get_next_bin_energy

  subroutine to_statepoint_energy(this, filter_group)
    class(EnergyFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "energy")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % bins)
  end subroutine to_statepoint_energy

  function text_label_energy(this, bin) result(label)
    class(EnergyFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    real(8) :: E0, E1

    E0 = this % bins(bin)
    E1 = this % bins(bin + 1)
    label = "Incoming Energy [" // trim(to_str(E0)) // ", " &
         // trim(to_str(E1)) // ")"
  end function text_label_energy

!===============================================================================
! EnergyoutFilter methods
!===============================================================================
  subroutine get_next_bin_energyout(this, p, estimator, current_bin, next_bin, &
       weight)
    class(EnergyoutFilter), intent(in)  :: this
    type(Particle),         intent(in)  :: p
    integer,                intent(in)  :: estimator
    integer, value,         intent(in)  :: current_bin
    integer,                intent(out) :: next_bin
    real(8),                intent(out) :: weight

    integer :: n

    if (current_bin == NO_BIN_FOUND) then
      n = this % n_bins

      if ((.not. run_CE) .and. this % matches_transport_groups) then
        next_bin = p % g

        ! Tallies are ordered in increasing groups, group indices
        ! however are the opposite, so switch
        next_bin = num_energy_groups - next_bin + 1
        weight = ONE

      else
        ! Check if energy of the particle is within energy bins.
        if (p % E < this % bins(1) .or. p % E > this % bins(n + 1)) then
          next_bin = NO_BIN_FOUND
          weight = ERROR_REAL
        else
          ! Search to find incoming energy bin.
          next_bin = binary_search(this % bins, n + 1, p % E)
          weight = ONE
        end if
      end if

    else
      next_bin = NO_BIN_FOUND
      weight = ERROR_REAL
    end if
  end subroutine get_next_bin_energyout

  subroutine to_statepoint_energyout(this, filter_group)
    class(EnergyoutFilter), intent(in) :: this
    integer(HID_T),         intent(in) :: filter_group

    call write_dataset(filter_group, "type", "energyout")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % bins)
  end subroutine to_statepoint_energyout

  function text_label_energyout(this, bin) result(label)
    class(EnergyoutFilter), intent(in) :: this
    integer,                intent(in) :: bin
    character(MAX_LINE_LEN)            :: label

    real(8) :: E0, E1

    E0 = this % bins(bin)
    E1 = this % bins(bin + 1)
    label = "Outgoing Energy [" // trim(to_str(E0)) // ", " &
         // trim(to_str(E1)) // ")"
  end function text_label_energyout

!===============================================================================
! DelayedGroupFilter methods
!===============================================================================
  subroutine get_next_bin_dg(this, p, estimator, current_bin, next_bin, weight)
    class(DelayedGroupFilter), intent(in)  :: this
    type(Particle),            intent(in)  :: p
    integer,                   intent(in)  :: estimator
    integer, value,            intent(in)  :: current_bin
    integer,                   intent(out) :: next_bin
    real(8),                   intent(out) :: weight

    if (current_bin == NO_BIN_FOUND) then
      next_bin = 1
      weight = ONE
    else
      next_bin = NO_BIN_FOUND
      weight = ERROR_REAL
    end if
  end subroutine get_next_bin_dg

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
  subroutine get_next_bin_mu(this, p, estimator, current_bin, next_bin, weight)
    class(MuFilter), intent(in)  :: this
    type(Particle),  intent(in)  :: p
    integer,         intent(in)  :: estimator
    integer, value,  intent(in)  :: current_bin
    integer,         intent(out) :: next_bin
    real(8),         intent(out) :: weight

    integer :: n

    if (current_bin == NO_BIN_FOUND) then
      n = this % n_bins

      ! Check if energy of the particle is within energy bins.
      if (p % mu < this % bins(1) .or. p % mu > this % bins(n + 1)) then
        next_bin = NO_BIN_FOUND
        weight = ERROR_REAL
      else
        ! Search to find incoming energy bin.
        next_bin = binary_search(this % bins, n + 1, p % mu)
        weight = ONE
      end if

    else
      next_bin = NO_BIN_FOUND
      weight = ERROR_REAL
    end if
  end subroutine get_next_bin_mu

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
  subroutine get_next_bin_polar(this, p, estimator, current_bin, next_bin, &
       weight)
    class(PolarFilter), intent(in)  :: this
    type(Particle),     intent(in)  :: p
    integer,            intent(in)  :: estimator
    integer, value,     intent(in)  :: current_bin
    integer,            intent(out) :: next_bin
    real(8),            intent(out) :: weight

    integer :: n
    real(8) :: theta

    if (current_bin == NO_BIN_FOUND) then
      n = this % n_bins

      ! Make sure the correct direction vector is used.
      if (estimator == ESTIMATOR_TRACKLENGTH) then
        theta = acos(p % coord(1) % uvw(3))
      else
        theta = acos(p % last_uvw(3))
      end if

      ! Check if particle is within polar angle bins.
      if (theta < this % bins(1) .or. theta > this % bins(n + 1)) then
        next_bin = NO_BIN_FOUND
        weight = ERROR_REAL
      else
        ! Search to find polar angle bin.
        next_bin = binary_search(this % bins, n + 1, theta)
        weight = ONE
      end if

    else
      next_bin = NO_BIN_FOUND
      weight = ERROR_REAL
    end if
  end subroutine get_next_bin_polar

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
  subroutine get_next_bin_azimuthal(this, p, estimator, current_bin, next_bin, &
       weight)
    class(AzimuthalFilter), intent(in)  :: this
    type(Particle),         intent(in)  :: p
    integer,                intent(in)  :: estimator
    integer, value,         intent(in)  :: current_bin
    integer,                intent(out) :: next_bin
    real(8),                intent(out) :: weight

    integer :: n
    real(8) :: phi

    if (current_bin == NO_BIN_FOUND) then
      n = this % n_bins

      ! Make sure the correct direction vector is used.
      if (estimator == ESTIMATOR_TRACKLENGTH) then
        phi = atan2(p % coord(1) % uvw(2), p % coord(1) % uvw(1))
      else
        phi = atan2(p % last_uvw(2), p % last_uvw(1))
      end if

      ! Check if particle is within azimuthal angle bins.
      if (phi < this % bins(1) .or. phi > this % bins(n + 1)) then
        next_bin = NO_BIN_FOUND
        weight = ERROR_REAL
      else
        ! Search to find azimuthal angle bin.
        next_bin = binary_search(this % bins, n + 1, phi)
        weight = ONE
      end if

    else
      next_bin = NO_BIN_FOUND
      weight = ERROR_REAL
    end if
  end subroutine get_next_bin_azimuthal

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
  subroutine get_next_bin_energyfunction(this, p, estimator, current_bin, &
       next_bin, weight)
    class(EnergyFunctionFilter), intent(in)  :: this
    type(Particle),              intent(in)  :: p
    integer,                     intent(in)  :: estimator
    integer, value,              intent(in)  :: current_bin
    integer,                     intent(out) :: next_bin
    real(8),                     intent(out) :: weight

    integer :: n, indx
    real(8) :: E, f

    select type(this)
    type is (EnergyFunctionFilter)
      if (current_bin == NO_BIN_FOUND) then
        n = size(this % energy)

        ! Make sure the correct energy is used.
        if (estimator == ESTIMATOR_TRACKLENGTH) then
          E = p % E
        else
          E = p % last_E
        end if

        ! Check if energy of the particle is within energy bins.
        if (E < this % energy(1) .or. E > this % energy(n)) then
          next_bin = NO_BIN_FOUND
          weight = ERROR_REAL

        else
          ! Search to find incoming energy bin.
          indx = binary_search(this % energy, n, E)

          ! Compute an interpolation factor between nearest bins.
          f = (E - this % energy(indx)) &
               / (this % energy(indx+1) - this % energy(indx))

          ! Interpolate on the lin-lin grid.
          next_bin = 1
          weight = (ONE - f) * this % y(indx) + f * this % y(indx+1)
        end if

      else
        next_bin = NO_BIN_FOUND
        weight = ERROR_REAL
      end if
    end select
  end subroutine get_next_bin_energyfunction

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
! FIND_OFFSET (for distribcell) uses a given map number, a target cell ID, and
! a target offset to build a string which is the path from the base universe to
! the target cell with the given offset
!===============================================================================

  recursive subroutine find_offset(i_cell, univ, target_offset, offset, path)

    integer, intent(in) :: i_cell         ! The target cell index
    type(Universe), intent(in) :: univ  ! Universe to begin search
    integer, intent(in) :: target_offset        ! Target offset
    integer, intent(inout) :: offset    ! Current offset
    character(*), intent(inout) :: path ! Path to offset

    integer :: map                  ! Index in maps vector
    integer :: i, j                 ! Index over cells
    integer :: k, l, m              ! Indices in lattice
    integer :: old_k, old_l, old_m  ! Previous indices in lattice
    integer :: n_x, n_y, n_z        ! Lattice cell array dimensions
    integer :: n                    ! Number of cells to search
    integer :: cell_index           ! Index in cells array
    integer :: lat_offset           ! Offset from lattice
    integer :: temp_offset          ! Looped sum of offsets
    integer :: i_univ               ! index in universes array
    logical :: this_cell = .false.  ! Advance in this cell?
    logical :: later_cell = .false. ! Fill cells after this one?
    type(Cell), pointer :: c           ! Pointer to current cell
    type(Universe), pointer :: next_univ  ! Next universe to loop through
    class(Lattice), pointer :: lat        ! Pointer to current lattice

    ! Get the distribcell index for this cell
    map = cells(i_cell) % distribcell_index

    n = size(univ % cells)

    ! Write to the geometry stack
    i_univ = universe_dict % get_key(univ % id)
    if (i_univ == root_universe) then
      path = trim(path) // "u" // to_str(univ%id)
    else
      path = trim(path) // "->u" // to_str(univ%id)
    end if

    ! Look through all cells in this universe
    do i = 1, n
      ! If the cell matches the goal and the offset matches final, write to the
      ! geometry stack
      if (univ % cells(i) == i_cell .and. offset == target_offset) then
        c => cells(univ % cells(i))
        path = trim(path) // "->c" // to_str(c % id)
        return
      end if
    end do

    ! Find the fill cell or lattice cell that we need to enter
    do i = 1, n

      later_cell = .false.

      c => cells(univ % cells(i))

      this_cell = .false.

      ! If we got here, we still think the target is in this universe
      ! or further down, but it's not this exact cell.
      ! Compare offset to next cell to see if we should enter this cell
      if (i /= n) then

        do j = i+1, n

          c => cells(univ % cells(j))

          ! Skip normal cells which do not have offsets
          if (c % type == FILL_MATERIAL) cycle

          ! Break loop once we've found the next cell with an offset
          exit
        end do

        ! Ensure we didn't just end the loop by iteration
        if (c % type /= FILL_MATERIAL) then

          ! There are more cells in this universe that it could be in
          later_cell = .true.

          ! Two cases, lattice or fill cell
          if (c % type == FILL_UNIVERSE) then
            temp_offset = c % offset(map)

          ! Get the offset of the first lattice location
          else
            lat => lattices(c % fill) % obj
            temp_offset = lat % offset(map, 1, 1, 1)
          end if

          ! If the final offset is in the range of offset - temp_offset+offset
          ! then the goal is in this cell
          if (target_offset < temp_offset + offset) then
            this_cell = .true.
          end if
        end if
      end if

      if (n == 1 .and. c % type /= FILL_MATERIAL) then
        this_cell = .true.
      end if

      if (.not. later_cell) then
        this_cell = .true.
      end if

      ! Get pointer to THIS cell because target must be in this cell
      if (this_cell) then

        cell_index = univ % cells(i)
        c => cells(cell_index)

        path = trim(path) // "->c" // to_str(c%id)

        ! ====================================================================
        ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
        if (c % type == FILL_UNIVERSE) then

          ! Enter this cell to update the current offset
          offset = c % offset(map) + offset

          next_univ => universes(c % fill)
          call find_offset(i_cell, next_univ, target_offset, offset, path)
          return

        ! ====================================================================
        ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
        elseif (c % type == FILL_LATTICE) then

          ! Set current lattice
          lat => lattices(c % fill) % obj

          select type (lat)

          ! ==================================================================
          ! RECTANGULAR LATTICES
          type is (RectLattice)

            ! Write to the geometry stack
            path = trim(path) // "->l" // to_str(lat%id)

            n_x = lat % n_cells(1)
            n_y = lat % n_cells(2)
            n_z = lat % n_cells(3)
            old_m = 1
            old_l = 1
            old_k = 1

            ! Loop over lattice coordinates
            do k = 1, n_x
              do l = 1, n_y
                do m = 1, n_z

                  if (target_offset >= lat % offset(map, k, l, m) + offset) then
                    if (k == n_x .and. l == n_y .and. m == n_z) then
                      ! This is last lattice cell, so target must be here
                      lat_offset = lat % offset(map, k, l, m)
                      offset = offset + lat_offset
                      next_univ => universes(lat % universes(k, l, m))
                      if (lat % is_3d) then
                        path = trim(path) // "(" // trim(to_str(k-1)) // &
                             "," // trim(to_str(l-1)) // "," // &
                             trim(to_str(m-1)) // ")"
                      else
                        path = trim(path) // "(" // trim(to_str(k-1)) // &
                             "," // trim(to_str(l-1)) // ")"
                      end if
                      call find_offset(i_cell, next_univ, target_offset, offset, path)
                      return
                    else
                      old_m = m
                      old_l = l
                      old_k = k
                      cycle
                    end if
                  else
                    ! Target is at this lattice position
                    lat_offset = lat % offset(map, old_k, old_l, old_m)
                    offset = offset + lat_offset
                    next_univ => universes(lat % universes(old_k, old_l, old_m))
                    if (lat % is_3d) then
                      path = trim(path) // "(" // trim(to_str(old_k-1)) // &
                           "," // trim(to_str(old_l-1)) // "," // &
                           trim(to_str(old_m-1)) // ")"
                    else
                      path = trim(path) // "(" // trim(to_str(old_k-1)) // &
                           "," // trim(to_str(old_l-1)) // ")"
                    end if
                    call find_offset(i_cell, next_univ, target_offset, offset, path)
                    return
                  end if

                end do
              end do
            end do

          ! ==================================================================
          ! HEXAGONAL LATTICES
          type is (HexLattice)

            ! Write to the geometry stack
            path = trim(path) // "->l" // to_str(lat%id)

            n_z = lat % n_axial
            n_y = 2 * lat % n_rings - 1
            n_x = 2 * lat % n_rings - 1
            old_m = 1
            old_l = 1
            old_k = 1

            ! Loop over lattice coordinates
            do m = 1, n_z
              do l = 1, n_y
                do k = 1, n_x

                  ! This array position is never used
                  if (k + l < lat % n_rings + 1) then
                    cycle
                  ! This array position is never used
                  else if (k + l > 3*lat % n_rings - 1) then
                    cycle
                  end if

                  if (target_offset >= lat % offset(map, k, l, m) + offset) then
                    if (k == lat % n_rings .and. l == n_y .and. m == n_z) then
                      ! This is last lattice cell, so target must be here
                      lat_offset = lat % offset(map, k, l, m)
                      offset = offset + lat_offset
                      next_univ => universes(lat % universes(k, l, m))
                      if (lat % is_3d) then
                        path = trim(path) // "(" // &
                             trim(to_str(k - lat % n_rings)) // "," // &
                             trim(to_str(l - lat % n_rings)) // "," // &
                             trim(to_str(m - 1)) // ")"
                      else
                        path = trim(path) // "(" // &
                             trim(to_str(k - lat % n_rings)) // "," // &
                             trim(to_str(l - lat % n_rings)) // ")"
                      end if
                      call find_offset(i_cell, next_univ, target_offset, offset, path)
                      return
                    else
                      old_m = m
                      old_l = l
                      old_k = k
                      cycle
                    end if
                  else
                    ! Target is at this lattice position
                    lat_offset = lat % offset(map, old_k, old_l, old_m)
                    offset = offset + lat_offset
                    next_univ => universes(lat % universes(old_k, old_l, old_m))
                    if (lat % is_3d) then
                      path = trim(path) // "(" // &
                           trim(to_str(old_k - lat % n_rings)) // "," // &
                           trim(to_str(old_l - lat % n_rings)) // "," // &
                           trim(to_str(old_m - 1)) // ")"
                    else
                      path = trim(path) // "(" // &
                           trim(to_str(old_k - lat % n_rings)) // "," // &
                           trim(to_str(old_l - lat % n_rings)) // ")"
                    end if
                    call find_offset(i_cell, next_univ, target_offset, offset, path)
                    return
                  end if

                end do
              end do
            end do

          end select

        end if
      end if
    end do
  end subroutine find_offset

end module tally_filter
