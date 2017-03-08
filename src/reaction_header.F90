module reaction_header

  use hdf5, only: HID_T, HSIZE_T, SIZE_T, h5gget_info_f, h5lget_name_by_idx_f, &
                  H5_INDEX_NAME_F, H5_ITER_INC_F

  use constants,      only: MAX_WORD_LEN
  use hdf5_interface, only: read_attribute, open_group, close_group, &
       open_dataset, read_dataset, close_dataset, get_shape, get_groups
  use product_header, only: ReactionProduct
  use stl_vector,     only: VectorInt
  use string,         only: to_str, starts_with

  implicit none

!===============================================================================
! REACTION contains the cross-section and secondary energy and angle
! distributions for a single reaction in a continuous-energy ACE-format table
!===============================================================================

  type TemperatureXS
    integer :: threshold             ! Energy grid index of threshold
    real(8), allocatable :: value(:) ! Cross section values
  end type TemperatureXS

  type Reaction
    integer :: MT                      ! ENDF MT value
    real(8) :: Q_value                 ! Reaction Q value
    logical :: scatter_in_cm           ! scattering system in center-of-mass?
    type(TemperatureXS), allocatable :: xs(:)
    type(ReactionProduct), allocatable :: products(:)
  contains
    procedure :: from_hdf5 => reaction_from_hdf5
  end type Reaction

contains

  subroutine reaction_from_hdf5(this, group_id, temperatures)
    class(Reaction), intent(inout) :: this
    integer(HID_T),  intent(in)    :: group_id
    type(VectorInt), intent(in)    :: temperatures

    integer :: i
    integer :: cm
    integer :: n_product
    integer(HID_T) :: pgroup
    integer(HID_T) :: xs, temp_group
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: j
    character(MAX_WORD_LEN) :: temp_str ! temperature dataset name, e.g. '294K'
    character(MAX_WORD_LEN), allocatable :: grp_names(:)

    call read_attribute(this % Q_value, group_id, 'Q_value')
    call read_attribute(this % MT, group_id, 'mt')
    call read_attribute(cm, group_id, 'center_of_mass')
    this % scatter_in_cm = (cm == 1)

    ! Read cross section and threshold_idx data
    allocate(this % xs(temperatures % size()))
    do i = 1, temperatures % size()
      temp_str = trim(to_str(temperatures % data(i))) // "K"
      temp_group = open_group(group_id, temp_str)
      xs = open_dataset(temp_group, 'xs')
      call read_attribute(this % xs(i) % threshold, xs, 'threshold_idx')
      call get_shape(xs, dims)
      allocate(this % xs(i) % value(dims(1)))
      call read_dataset(this % xs(i) % value, xs)
      call close_dataset(xs)
      call close_group(temp_group)
    end do

    ! Determine number of products
    n_product = 0
    call get_groups(group_id, grp_names)
    do j = 1, size(grp_names)
      if (starts_with(grp_names(j), "product_")) n_product = n_product + 1
    end do

    ! Read products
    allocate(this % products(n_product))
    do i = 1, n_product
      pgroup = open_group(group_id, 'product_' // trim(to_str(i - 1)))
      call this % products(i) % from_hdf5(pgroup)
      call close_group(pgroup)
    end do
  end subroutine reaction_from_hdf5

end module reaction_header
