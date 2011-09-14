module cross_section_header

  implicit none

!===============================================================================
! DISTANGLE contains data for a tabular secondary angle distribution whether it
! be tabular or 32 equiprobable cosine bins
!===============================================================================

  type DistAngle
     integer              :: n_energy    ! # of incoming energies
     real(8), allocatable :: energy(:)   ! incoming energy grid
     integer, allocatable :: type(:)     ! type of distribution
     integer, allocatable :: location(:) ! location of each table
     real(8), allocatable :: data(:)     ! angular distribution data
  end type DistAngle

!===============================================================================
! DISTENERGY contains data for a secondary energy distribution for all
! scattering laws
!===============================================================================

  type DistEnergy
     integer :: law                    ! secondary distribution law
     integer :: n_interp               ! # of interpolation regions
     integer, allocatable :: nbt(:)    ! ENDF interpolation parameters
     integer, allocatable :: int(:)    ! ''
     integer :: n_energy               ! # of energies for law validity
     real(8), allocatable :: energy(:) ! energy grid for law validity
     real(8), allocatable :: pvalid(:) ! probability of law validity
     real(8), allocatable :: data(:)   ! energy distribution data
  end type DistEnergy

!===============================================================================
! REACTION contains the cross-section and secondary energy and angle
! distributions for a single reaction in a continuous-energy ACE-format table
!===============================================================================

  type Reaction
     integer :: MT                     ! ENDF MT value
     real(8) :: Q_value                ! Reaction Q value
     integer :: TY                     ! Number of neutrons released
     integer :: IE                     ! Starting energy grid index
     real(8), allocatable :: sigma(:)  ! Cross section values
     logical :: has_angle_dist         ! Angle distribution present?
     logical :: has_energy_dist        ! Energy distribution present?
     type(DistAngle)  :: adist         ! Secondary angular distribution
     type(DistEnergy) :: edist         ! Secondary energy distribution
  end type Reaction

!===============================================================================
! NUCLIDE contains all the data for an ACE-format continuous-energy cross
! section. The ACE format (A Compact ENDF format) is used in MCNP and several
! other Monte Carlo codes.
!===============================================================================

  type Nuclide
     character(20) :: name
     real(8)       :: awr
     real(8)       :: temp
     logical       :: fissionable

     ! Energy grid information
     integer :: n_grid
     integer, allocatable :: grid_index(:)
     real(8), allocatable :: energy(:)

     ! Cross sections
     real(8), allocatable :: total(:)
     real(8), allocatable :: elastic(:)
     real(8), allocatable :: fission(:)
     real(8), allocatable :: absorption(:)
     real(8), allocatable :: heating(:)

     ! Total fission neutron emission
     integer :: nu_t_type
     real(8), allocatable :: nu_t_data(:)

     ! Prompt fission neutron emission
     integer :: nu_p_type
     real(8), allocatable :: nu_p_data(:)
     
     ! Delayed fission neutron emission
     integer :: nu_d_type
     integer :: n_precursor
     real(8), allocatable :: nu_d_data(:)
     real(8), allocatable :: nu_d_precursor_data(:)
     type(DistEnergy), allocatable :: nu_d_edist(:)

     ! Unresolved resonance data
     logical :: urr_present
     integer, allocatable :: urr_params(:)
     real(8), allocatable :: urr_energy(:)
     real(8), allocatable :: urr_prob(:,:,:)

     ! Reactions
     integer :: n_reaction
     type(Reaction), pointer :: reactions(:) => null()

  end type Nuclide

!===============================================================================
! SAB_TABLE contains S(a,b) data for thermal neutron scattering, typically off
! of light isotopes such as water, graphite, Be, etc
!===============================================================================
     
  type SAB_Table
     character(20) :: name
     real(8) :: awr
     real(8) :: temp
     integer :: n_inelastic_e_in
     integer :: n_inelastic_e_out
     integer :: n_inelastic_mu
     real(8), allocatable :: inelastic_e_in(:)
     real(8), allocatable :: inelastic_sigma(:) 
     real(8), allocatable :: inelastic_e_out(:,:)
     real(8), allocatable :: inelastic_mu(:,:,:)
     integer :: n_elastic_e_in
     integer :: n_elastic_type
     integer :: n_elastic_mu
     real(8), allocatable :: elastic_e_in(:)
     real(8), allocatable :: elastic_P(:)
     real(8), allocatable :: elastic_mu(:,:)
  end type SAB_Table

!===============================================================================
! XSDATA contains data read in from a SERPENT xsdata file
!===============================================================================

  type xsData
     character(10) :: alias
     character(10) :: id
     integer :: type
     integer :: zaid
     integer :: isomeric
     real(8) :: awr
     real(8) :: temp
     integer :: binary
     character(150) :: path
  end type xsData

!===============================================================================
! NUCLIDEMICROXS contains cached microscopic cross sections for a
! particular nuclide at the current energy
!===============================================================================

  type NuclideMicroXS
     integer :: index_grid
     integer :: index_temp
     integer :: last_index_grid
     integer :: last_index_temp
     real(8) :: interp_factor
     real(8) :: total
     real(8) :: elastic
     real(8) :: absorption
     real(8) :: fission
     real(8) :: nu_fission
  end type NuclideMicroXS

!===============================================================================
! MATERIALMACROXS contains cached macroscopic cross sections for the material a
! particle is traveling through
!===============================================================================

  type MaterialMacroXS
     real(8) :: total
     real(8) :: scatter
     real(8) :: elastic
     real(8) :: absorption
     real(8) :: fission
     real(8) :: nu_fission
  end type MaterialMacroXS

end module cross_section_header
