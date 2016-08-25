module reaction_header

  use hdf5, only: HID_T, HSIZE_T, SIZE_T, h5gget_info_f, h5lget_name_by_idx_f, &
                  H5_INDEX_NAME_F, H5_ITER_INC_F

  use constants,      only: MAX_WORD_LEN
  use hdf5_interface, only: read_attribute, open_group, close_group, &
       open_dataset, read_dataset, close_dataset, get_shape
  use product_header, only: ReactionProduct
  use string,         only: to_str, starts_with

  implicit none

!===============================================================================
! REACTION contains the cross-section and secondary energy and angle
! distributions for a single reaction in a continuous-energy ACE-format table
!===============================================================================

  type Reaction
    integer :: MT                      ! ENDF MT value
    real(8) :: Q_value                 ! Reaction Q value
    integer :: threshold               ! Energy grid index of threshold
    logical :: scatter_in_cm           ! scattering system in center-of-mass?
    real(8), allocatable :: sigma(:)   ! Cross section values
    type(ReactionProduct), allocatable :: products(:)
  contains
    procedure :: from_hdf5 => reaction_from_hdf5
  end type Reaction

contains

  subroutine reaction_from_hdf5(this, group_id, temperature)
    class(Reaction), intent(inout) :: this
    integer(HID_T),  intent(in)    :: group_id
    character(6),    intent(in)    :: temperature

    integer :: i
    integer :: cm
    integer :: n_product
    integer :: storage_type
    integer :: max_corder
    integer :: n_links
    integer :: hdf5_err
    integer(HID_T) :: pgroup
    integer(HID_T) :: xs, xs_group
    integer(SIZE_T) :: name_len
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: j
    character(MAX_WORD_LEN) :: name

    call read_attribute(this % Q_value, group_id, 'Q_value')
    call read_attribute(this % MT, group_id, 'mt')
    call read_attribute(cm, group_id, 'center_of_mass')
    this % scatter_in_cm = (cm == 1)

    ! Read cross section and threshold_idx data
    xs_group = open_group(group_id, temperature)
    xs = open_dataset(xs_group, 'xs')
    call read_attribute(this % threshold, xs, 'threshold_idx')
    call get_shape(xs, dims)
    allocate(this % sigma(dims(1)))
    call read_dataset(this % sigma, xs)
    call close_dataset(xs)
    call close_group(xs_group)

    ! Determine number of products
    call h5gget_info_f(group_id, storage_type, n_links, max_corder, hdf5_err)
    n_product = 0
    do j = 0, n_links - 1
      call h5lget_name_by_idx_f(group_id, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, &
           j, name, hdf5_err, name_len)
      if (starts_with(name, "product_")) n_product = n_product + 1
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
