module tally_filter_material

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T

  use constants
  use dict_header,     only: DictIntInt, EMPTY
  use error
  use hdf5_interface
  use material_header, only: materials, material_dict
  use particle_header, only: Particle
  use string,          only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private
  public :: openmc_material_filter_get_bins
  public :: openmc_material_filter_set_bins

!===============================================================================
! MATERIAL specifies which material tally events reside in.
!===============================================================================

  type, public, extends(TallyFilter) :: MaterialFilter
    integer, allocatable :: materials(:)
    type(DictIntInt)     :: map
  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_material
    procedure :: to_statepoint => to_statepoint_material
    procedure :: text_label => text_label_material
    procedure :: initialize => initialize_material
  end type MaterialFilter

contains

  subroutine from_xml(this, node)
    class(MaterialFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    n = node_word_count(node, "bins")

    ! Allocate and store bins
    this % n_bins = n
    allocate(this % materials(n))
    call get_node_array(node, "bins", this % materials)
  end subroutine from_xml

  subroutine get_all_bins_material(this, p, estimator, match)
    class(MaterialFilter), intent(in)  :: this
    type(Particle),        intent(in)  :: p
    integer,               intent(in)  :: estimator
    type(TallyFilterMatch),     intent(inout) :: match

    integer :: val

    val = this % map % get(p % material)
    if (val /= EMPTY) then
      call match % bins % push_back(val)
      call match % weights % push_back(ONE)
    end if

  end subroutine get_all_bins_material

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
    integer :: val

    ! Convert ids to indices.
    do i = 1, this % n_bins
      id = this % materials(i)
      val = material_dict % get(id)
      if (val /= EMPTY) then
        this % materials(i) = val
      else
        call fatal_error("Could not find material " // trim(to_str(id)) &
             &// " specified on a tally filter.")
      end if
    end do

    ! Generate mapping from material indices to filter bins.
    do i = 1, this % n_bins
      call this % map % set(this % materials(i), i)
    end do
  end subroutine initialize_material

  function text_label_material(this, bin) result(label)
    class(MaterialFilter), intent(in) :: this
    integer,               intent(in) :: bin
    character(MAX_LINE_LEN)           :: label

    label = "Material " // to_str(materials(this % materials(bin)) % id)
  end function text_label_material

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_material_filter_get_bins(index, bins, n) result(err) bind(C)
    ! Return the bins for a material filter
    integer(C_INT32_T), value :: index
    type(C_PTR), intent(out) :: bins
    integer(C_INT32_T), intent(out) :: n
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_filters) then
      if (allocated(filters(index) % obj)) then
        select type (f => filters(index) % obj)
        type is (MaterialFilter)
          bins = C_LOC(f % materials)
          n = size(f % materials)
          err = 0
        class default
          err = E_INVALID_TYPE
          call set_errmsg("Tried to get material filter bins on a &
               &non-material filter.")
        end select
      else
        err = E_ALLOCATE
        call set_errmsg("Filter type has not been set yet.")
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in filters array out of bounds.")
    end if
  end function openmc_material_filter_get_bins


  function openmc_material_filter_set_bins(index, n, bins) result(err) bind(C)
    ! Set the materials for the filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), intent(in) :: bins(n)
    integer(C_INT) :: err

    integer :: i

    err = 0
    if (index >= 1 .and. index <= n_filters) then
      if (allocated(filters(index) % obj)) then
        select type (f => filters(index) % obj)
        type is (MaterialFilter)
          f % n_bins = n
          if (allocated(f % materials)) deallocate(f % materials)
          allocate(f % materials(n))
          f % materials(:) = bins

          ! Generate mapping from material indices to filter bins.
          call f % map % clear()
          do i = 1, n
            call f % map % set(f % materials(i), i)
          end do

        class default
          err = E_INVALID_TYPE
          call set_errmsg("Tried to set material filter bins on a &
               &non-material filter.")
        end select
      else
        err = E_ALLOCATE
        call set_errmsg("Filter type has not been set yet.")
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in filters array out of bounds.")
    end if
  end function openmc_material_filter_set_bins

end module tally_filter_material
