module sab_header

  use, intrinsic :: ISO_FORTRAN_ENV

  use constants
  use string, only: to_str

  implicit none

!===============================================================================
! DISTENERGYSAB contains the secondary energy/angle distributions for inelastic
! thermal scattering collisions which utilize a continuous secondary energy
! representation.
!===============================================================================

  type DistEnergySab
    integer              :: n_e_out
    real(8), allocatable :: e_out(:)
    real(8), allocatable :: e_out_pdf(:)
    real(8), allocatable :: e_out_cdf(:)
    real(8), allocatable :: mu(:,:)
  end type DistEnergySab

!===============================================================================
! SALPHABETA contains S(a,b) data for thermal neutron scattering, typically off
! of light isotopes such as water, graphite, Be, etc
!===============================================================================

  type SAlphaBeta
    character(10) :: name     ! name of table, e.g. lwtr.10t
    real(8)       :: awr      ! weight of nucleus in neutron masses
    real(8)       :: kT       ! temperature in MeV (k*T)
    integer       :: n_zaid   ! Number of valid zaids
    integer, allocatable :: zaid(:) ! List of valid Z and A identifiers, e.g. 6012

    ! threshold for S(a,b) treatment (usually ~4 eV)
    real(8) :: threshold_inelastic
    real(8) :: threshold_elastic = ZERO

    ! Inelastic scattering data
    integer :: n_inelastic_e_in  ! # of incoming E for inelastic
    integer :: n_inelastic_e_out ! # of outgoing E for inelastic
    integer :: n_inelastic_mu    ! # of outgoing angles for inelastic
    integer :: secondary_mode    ! secondary mode (equal/skewed/continuous)
    real(8), allocatable :: inelastic_e_in(:)
    real(8), allocatable :: inelastic_sigma(:)
    ! The following are used only if secondary_mode is 0 or 1
    real(8), allocatable :: inelastic_e_out(:,:)
    real(8), allocatable :: inelastic_mu(:,:,:)
    ! The following is used only if secondary_mode is 3
    ! The different implementation is necessary because the continuous
    ! representation has a variable number of outgoing energy points for each
    ! incoming energy
    type(DistEnergySab), allocatable :: inelastic_data(:) ! One for each Ein

    ! Elastic scattering data
    integer :: elastic_mode   ! elastic mode (discrete/exact)
    integer :: n_elastic_e_in ! # of incoming E for elastic
    integer :: n_elastic_mu   ! # of outgoing angles for elastic
    real(8), allocatable :: elastic_e_in(:)
    real(8), allocatable :: elastic_P(:)
    real(8), allocatable :: elastic_mu(:,:)
  contains
      procedure :: print => print_sab_table
  end type SAlphaBeta

  contains

!===============================================================================
! PRINT_SAB_TABLE displays information about a S(a,b) table containing data
! describing thermal scattering from bound materials such as hydrogen in water.
!===============================================================================

    subroutine print_sab_table(this, unit)
      class(SAlphaBeta), intent(in)  :: this
      integer, intent(in), optional :: unit

      integer :: size_sab   ! memory used by S(a,b) table
      integer :: unit_      ! unit to write to
      integer :: i          ! Loop counter for parsing through this % zaid
      integer :: char_count ! Counter for the number of characters on a line

      ! set default unit for writing information
      if (present(unit)) then
        unit_ = unit
      else
        unit_ = OUTPUT_UNIT
      end if

      ! Basic S(a,b) table information
      write(unit_,*) 'S(a,b) Table ' // trim(this % name)
      write(unit_,'(A)',advance="no") '   zaids = '
      ! Initialize the counter based on the above string
      char_count = 11
      do i = 1, this % n_zaid
        ! Deal with a line thats too long
        if (char_count >= 73) then  ! 73 = 80 - (5 ZAID chars + 1 space + 1 comma)
          ! End the line
          write(unit_,*) ""
          ! Add 11 leading blanks
          write(unit_,'(A)', advance="no") "           "
          ! reset the counter to 11
          char_count = 11
        end if
        if (i < this % n_zaid) then
          ! Include a comma
          write(unit_,'(A)',advance="no") trim(to_str(this % zaid(i))) // ", "
          char_count = char_count + len(trim(to_str(this % zaid(i)))) + 2
        else
          ! Don't include a comma, since we are all done
          write(unit_,'(A)',advance="no") trim(to_str(this % zaid(i)))
        end if

      end do
      write(unit_,*) "" ! Move to next line
      write(unit_,*) '  awr = ' // trim(to_str(this % awr))
      write(unit_,*) '  kT = ' // trim(to_str(this % kT))

      ! Inelastic data
      write(unit_,*) '  # of Incoming Energies (Inelastic) = ' // &
           trim(to_str(this % n_inelastic_e_in))
      write(unit_,*) '  # of Outgoing Energies (Inelastic) = ' // &
           trim(to_str(this % n_inelastic_e_out))
      write(unit_,*) '  # of Outgoing Angles (Inelastic) = ' // &
           trim(to_str(this % n_inelastic_mu))
      write(unit_,*) '  Threshold for Inelastic = ' // &
           trim(to_str(this % threshold_inelastic))

      ! Elastic data
      if (this % n_elastic_e_in > 0) then
        write(unit_,*) '  # of Incoming Energies (Elastic) = ' // &
             trim(to_str(this % n_elastic_e_in))
        write(unit_,*) '  # of Outgoing Angles (Elastic) = ' // &
             trim(to_str(this % n_elastic_mu))
        write(unit_,*) '  Threshold for Elastic = ' // &
             trim(to_str(this % threshold_elastic))
      end if

      ! Determine memory used by S(a,b) table and write out
      size_sab = 8 * (this % n_inelastic_e_in * (2 + this % n_inelastic_e_out * &
           (1 + this % n_inelastic_mu)) + this % n_elastic_e_in * &
           (2 + this % n_elastic_mu))
      write(unit_,*) '  Memory Used = ' // trim(to_str(size_sab)) // ' bytes'

      ! Blank line at end
      write(unit_,*)

    end subroutine print_sab_table

end module sab_header
