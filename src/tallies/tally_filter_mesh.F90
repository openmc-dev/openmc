module tally_filter_mesh

  use, intrinsic :: ISO_C_BINDING

  use constants
  use dict_header,         only: EMPTY
  use error
  use mesh_header
  use hdf5_interface
  use particle_header,     only: Particle
  use string,              only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private
  public :: openmc_mesh_filter_get_mesh
  public :: openmc_mesh_filter_set_mesh

!===============================================================================
! MESHFILTER indexes the location of particle events to a regular mesh.  For
! tracklength tallies, it will produce multiple valid bins and the bin weight
! will correspond to the fraction of the track length that lies in that bin.
!===============================================================================

  type, public, extends(CppTallyFilter) :: MeshFilter
    integer :: mesh
  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_mesh
    procedure :: to_statepoint => to_statepoint_mesh
    procedure :: text_label => text_label_mesh
  end type MeshFilter

contains

  subroutine from_xml(this, node)
    class(MeshFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: i
    integer :: id
    integer :: n
    integer(C_INT) :: err
    type(RegularMesh) :: m

    call this % from_xml_c(node)

    n = node_word_count(node, "bins")

    if (n /= 1) call fatal_error("Only one mesh can be &
         &specified per mesh filter.")

    ! Determine id of mesh
    call get_node_value(node, "bins", id)

    ! Get pointer to mesh
    err = openmc_get_mesh_index(id, this % mesh)
    if (err /= 0) then
      call fatal_error("Could not find mesh " // trim(to_str(id)) &
           // " specified on filter.")
    end if

    ! Determine number of bins
    m = meshes(this % mesh)
    this % n_bins = 1
    do i = 1, m % n_dimension()
      this % n_bins = this % n_bins * m % dimension(i)
    end do
  end subroutine from_xml

  subroutine get_all_bins_mesh(this, p, estimator, match)
    class(MeshFilter), intent(in)  :: this
    type(Particle),    intent(in)  :: p
    integer,           intent(in)  :: estimator
    type(TallyFilterMatch), intent(inout) :: match

    call this % get_all_bins_c(p, estimator, match)

  end subroutine get_all_bins_mesh

  subroutine to_statepoint_mesh(this, filter_group)
    class(MeshFilter), intent(in) :: this
    integer(HID_T),    intent(in) :: filter_group

    type(RegularMesh) :: m

    m = meshes(this % mesh)
    call write_dataset(filter_group, "type", "mesh")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", m % id())
  end subroutine to_statepoint_mesh

  function text_label_mesh(this, bin) result(label)
    class(MeshFilter), intent(in) :: this
    integer,           intent(in) :: bin
    character(MAX_LINE_LEN)       :: label

    integer, allocatable       :: ijk(:)
    type(RegularMesh) :: m

    m = meshes(this % mesh)
    allocate(ijk(m % n_dimension()))
    call m % get_indices_from_bin(bin, ijk)
    if (m % n_dimension() == 1) then
      label = "Mesh Index (" // trim(to_str(ijk(1))) // ")"
    elseif (m % n_dimension() == 2) then
      label = "Mesh Index (" // trim(to_str(ijk(1))) // ", " // &
            trim(to_str(ijk(2))) // ")"
    elseif (m % n_dimension() == 3) then
      label = "Mesh Index (" // trim(to_str(ijk(1))) // ", " // &
            trim(to_str(ijk(2))) // ", " // trim(to_str(ijk(3))) // ")"
    end if
  end function text_label_mesh

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_mesh_filter_get_mesh(index, index_mesh) result(err) bind(C)
    ! Get the mesh for a mesh filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), intent(out)       :: index_mesh
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MeshFilter)
        index_mesh = f % mesh
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to set mesh on a non-mesh filter.")
      end select
    end if
  end function openmc_mesh_filter_get_mesh


  function openmc_mesh_filter_set_mesh(index, index_mesh) result(err) bind(C)
    ! Set the mesh for a mesh filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: index_mesh
    integer(C_INT) :: err

    type(RegularMesh) :: m
    integer :: i

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MeshFilter)
        if (index_mesh >= 0 .and. index_mesh < n_meshes()) then
          f % mesh = index_mesh
          f % n_bins = 1
          m = meshes(index_mesh)
          do i = 1, m % n_dimension()
            f % n_bins = f % n_bins * m % dimension(i)
          end do
        else
          err = E_OUT_OF_BOUNDS
          call set_errmsg("Index in 'meshes' array is out of bounds.")
        end if
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to set mesh on a non-mesh filter.")
      end select
    end if
  end function openmc_mesh_filter_set_mesh

end module tally_filter_mesh
