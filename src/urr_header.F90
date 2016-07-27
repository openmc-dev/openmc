module urr_header

  use hdf5,           only: HID_T, HSIZE_T
  use hdf5_interface, only: read_attribute, open_dataset, read_dataset, &
       close_dataset, get_shape

  implicit none

!===============================================================================
! URRDATA contains probability tables for the unresolved resonance range.
!===============================================================================

  type UrrData
    integer :: n_energy        ! # of incident neutron energies
    integer :: n_prob          ! # of probabilities
    integer :: interp          ! inteprolation (2=lin-lin, 5=log-log)
    integer :: inelastic_flag  ! inelastic competition flag
    integer :: absorption_flag ! other absorption flag
    logical :: multiply_smooth ! multiply by smooth cross section?
    real(8), allocatable :: energy(:)   ! incident energies
    real(8), allocatable :: prob(:,:,:) ! actual probabibility tables
  contains
    procedure :: from_hdf5 => urr_from_hdf5
  end type UrrData

contains

  subroutine urr_from_hdf5(this, group_id)
    class(UrrData), intent(inout) :: this
    integer(HID_T), intent(in)    :: group_id

    integer :: i, j, k
    integer(HID_T) :: energy
    integer(HID_T) :: table
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: dims3(3)
    real(8), allocatable :: temp(:,:,:)

    ! Read interpolation and other flags
    call read_attribute(this % interp, group_id, 'interpolation')
    call read_attribute(this % inelastic_flag, group_id, 'inelastic')
    call read_attribute(this % absorption_flag, group_id, 'absorption')
    call read_attribute(i, group_id, 'multiply_smooth')
    this % multiply_smooth = (i == 1)

    ! Read energies at which tables exist
    energy = open_dataset(group_id, 'energy')
    call get_shape(energy, dims)
    this % n_energy = int(dims(1), 4)
    allocate(this % energy(this % n_energy))
    call read_dataset(this % energy, energy)
    call close_dataset(energy)

    ! Read URR tables
    table = open_dataset(group_id, 'table')
    call get_shape(table, dims3)
    this % n_prob = int(dims3(1), 4)
    allocate(temp(this % n_prob, 6, this % n_energy))
    call read_dataset(temp, table)
    call close_dataset(table)

    ! Swap first and last indices
    allocate(this % prob(this % n_energy, 6, this % n_prob))
    do i = 1, this % n_energy
      do j = 1, 6
        do k = 1, this % n_prob
          this % prob(i, j, k) = temp(k, j, i)
        end do
      end do
    end do
  end subroutine urr_from_hdf5

end module urr_header
