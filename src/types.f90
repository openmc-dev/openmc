module types

  implicit none

!===============================================================================
! UNIVERSE defines a geometry that fills all phase space
!===============================================================================

  type Universe
     integer :: uid                   ! Unique ID
     integer :: type                  ! Type
     integer :: level                 ! Level of universe (0=base)
     integer :: n_cells               ! # of cells within
     integer, allocatable :: cells(:) ! List of cells within
     real(8) :: x0                    ! Translation in x-coordinate
     real(8) :: y0                    ! Translation in y-coordinate
     real(8) :: z0                    ! Translation in z-coordinate
     integer, allocatable :: tallies(:)
  end type Universe

!===============================================================================
! LATTICE is an ordered array of elements (either rectangular, hexagonal, or
! triangular)
!===============================================================================

  type Lattice
     integer :: uid      ! Universe number for lattice
     integer :: type     ! Type of lattice (rectangular, hex, etc)
     integer :: level    ! Level of lattice
     integer :: n_x      ! number of lattice cells in x-direction
     integer :: n_y      ! number of lattice cells in y-direction
     real(8) :: x0       ! x-coordinate of lattice origin
     real(8) :: y0       ! y-coordinate of lattice origin
     real(8) :: width_x  ! width of lattice cell 
     real(8) :: width_y  ! width of lattice cell
     integer, allocatable :: element(:,:) ! specified universes
     integer, allocatable :: tallies(:)
  end type Lattice

!===============================================================================
! SURFACE type defines a first- or second-order surface that can be used to
! construct closed volumes (cells)
!===============================================================================

  type Surface
     integer :: uid                    ! Unique ID
     integer :: type                   ! Type of surface
     real(8), allocatable :: coeffs(:) ! Definition of surface
     integer, allocatable :: & 
          & neighbor_pos(:), &         ! List of cells on positive side
          & neighbor_neg(:)            ! List of cells on negative side
     integer :: bc                     ! Boundary condition
  end type Surface

!===============================================================================
! CELL defines a closed volume by its bounding surfaces
!===============================================================================

  type Cell
     integer :: uid        ! Unique ID
     integer :: type       ! Type of cell (normal, universe, lattice)
     integer :: universe   ! universe # this cell is in
     integer :: fill       ! universe # filling this cell
     integer :: parent     ! cell within which this cell resides
     integer :: material   ! Material within cell (0 for universe)
     integer :: n_surfaces ! Number of surfaces within
     integer, allocatable :: & 
          & surfaces(:)    ! List of surfaces bounding cell -- note that
                           ! parentheses, union, etc operators will be listed
                           ! here too
     integer, allocatable :: tallies(:)
  end type Cell

!===============================================================================
! PARTICLE describes the state of a particle being transported through the
! geometry
!===============================================================================

  type Particle
    integer(8) :: uid          ! Unique ID
    integer    :: type         ! Particle type (n, p, e, etc)
    real(8)    :: xyz(3)       ! location
    real(8)    :: xyz_local(3) ! local location (after transformations)
    real(8)    :: uvw(3)       ! directional cosines
    real(8)    :: wgt          ! particle weight
    real(8)    :: E            ! energy
    integer    :: IE           ! index on energy grid
    real(8)    :: interp       ! interpolation factor for energy grid
    integer    :: cell         ! current cell
    integer    :: universe     ! current universe
    integer    :: lattice      ! current lattice
    integer    :: surface      ! current surface
    integer    :: index_x      ! lattice index for x direction
    integer    :: index_y      ! lattice index for y direction
    logical    :: alive        ! is particle alive?
    integer    :: n_coll       ! # of collisions
  end type Particle

!===============================================================================
! BANK is used for storing fission sites in criticality calculations. Since all
! the state information of a neutron is not needed, this type allows sites to be
! stored with less memory
!===============================================================================

  type Bank
     integer(8) :: uid    ! Unique ID
     real(8)    :: xyz(3) ! location of bank particle
     real(8)    :: uvw(3) ! diretional cosines
     real(8)    :: E      ! energy
  end type Bank

!===============================================================================
! TALLYSCORE
!===============================================================================

  type TallyScore
     integer :: n_events
     real(8) :: val
     real(8) :: val_sq
  end type TallyScore

!===============================================================================
! TALLY
!===============================================================================

  type Tally
     integer :: uid
     integer :: type
     real(8) :: volume
     integer :: cell_type
     integer :: reaction_type
     integer :: material_type
     integer, allocatable :: reactions(:)
     integer, allocatable :: cells(:)
     integer, allocatable :: materials(:)
     integer, allocatable :: universes(:)
     real(8), allocatable :: energies(:)
     real(8) :: xyz_min(3)
     real(8) :: xyz_max(3)
     integer :: n_x
     integer :: n_y
     integer :: n_z
     type(TallyScore), allocatable :: score(:,:,:)
  end type Tally

!===============================================================================
! MATERIAL describes a material by its constituent isotopes
!===============================================================================

  type Material
     integer :: uid
     integer :: n_isotopes
     character(10), allocatable :: names(:)  ! isotope names
     integer, allocatable :: isotopes(:)     ! index in xsdata list
     integer, allocatable :: table(:)        ! index in xs array
     real(8)              :: atom_density    ! total atom density in atom/b-cm
     real(8), allocatable :: atom_percent(:)
     real(8), allocatable :: total_xs(:)     ! macroscopic cross-section
     integer :: sab_table
  end type Material

!===============================================================================
! EXTSOURCE describes an external source of neutrons for a fixed-source problem
! or for the starting source in a criticality problem
!===============================================================================

  type ExtSource
     integer :: type                    ! type of source, e.g. 'box' or 'cell'
     real(8), allocatable :: values(:)  ! values for particular source type
  end type ExtSource

!===============================================================================
! ACEDISTANGLE contains data for a tabular secondary angle distribution whether
! it be tabular or 32 equiprobable cosine bins
!===============================================================================

  type AceDistAngle
     integer              :: n_energy    ! # of incoming energies
     real(8), allocatable :: energy(:)   ! incoming energy grid
     integer, allocatable :: type(:)     ! type of distribution
     integer, allocatable :: location(:) ! location of each table
     real(8), allocatable :: data(:)     ! angular distribution data
  end type AceDistAngle

!===============================================================================
! ACEDISTENERGY contains data for a secondary energy distribution for all
! scattering laws
!===============================================================================

  type AceDistEnergy
     integer :: law                    ! secondary distribution law
     integer :: n_interp               ! # of interpolation regions
     integer, allocatable :: nbt(:)    ! ENDF interpolation parameters
     integer, allocatable :: int(:)    ! ''
     integer :: n_energy               ! # of energies for law validity
     real(8), allocatable :: energy(:) ! energy grid for law validity
     real(8), allocatable :: pvalid(:) ! probability of law validity
     real(8), allocatable :: data(:)   ! energy distribution data
  end type AceDistEnergy

!===============================================================================
! ACEREACTION contains the cross-section and secondary energy and angle
! distributions for a single reaction in a continuous-energy ACE-format table
!===============================================================================

  type AceReaction
     integer :: MT                     ! ENDF MT value
     real(8) :: Q_value                ! Reaction Q value
     integer :: TY                     ! Number of neutrons released
     integer :: IE                     ! Starting energy grid index
     real(8), allocatable :: sigma(:)  ! Cross section values
     logical :: has_angle_dist         ! Angle distribution present?
     logical :: has_energy_dist        ! Energy distribution present?
     type(AceDistAngle)  :: adist      ! Secondary angular distribution
     type(AceDistEnergy) :: edist      ! Secondary energy distribution
  end type AceReaction

!===============================================================================
! ACECONTINUOUS contains all the data for an ACE-format continuous-energy cross
! section. The ACE format (A Compact ENDF format) is used in MCNP and several
! other Monte Carlo codes.
!===============================================================================

  type AceContinuous
     character(20) :: name
     real(8) :: awr
     real(8) :: temp
     integer :: n_grid
     integer, allocatable :: grid_index(:)
     real(8), allocatable :: energy(:)
     real(8), allocatable :: sigma_t(:)
     real(8), allocatable :: sigma_a(:)
     real(8), allocatable :: sigma_el(:)
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
     type(AceDistEnergy), allocatable :: nu_d_edist(:)

     ! Unresolved resonance data
     logical :: urr_present
     integer, allocatable :: urr_params(:)
     real(8), allocatable :: urr_energy(:)
     real(8), allocatable :: urr_prob(:,:,:)

     ! Reactions
     integer :: n_reaction
     type(AceReaction), pointer :: reactions(:) => null()

  end type AceContinuous

!===============================================================================
! ACETHERMAL contains S(a,b) data for thermal neutron scattering, typically off
! of light isotopes such as water, graphite, Be, etc
!===============================================================================
     
  type AceThermal
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
  end type AceThermal

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
! TIMEROBJ represents a timer that can be started and stopped to measure how
! long different routines run. The intrinsic routine system_clock is used to
! measure time rather than cpu_time.
!===============================================================================

  type TimerObj
     logical :: running      = .false. ! is timer running?
     integer :: start_counts = 0       ! counts when started
     real(8) :: elapsed      = 0.      ! total time elapsed in seconds
  end type TimerObj

!===============================================================================
! KEYVALUECI stores the (key,value) pair for a dictionary where the key is a
! string and the value is an integer. Note that we need to store the key in
! addition to the value for collision resolution.
!===============================================================================

  ! Key length for dictionary
  integer, parameter :: DICT_KEY_LENGTH = 20

  type KeyValueCI
     character(len=DICT_KEY_LENGTH) :: key
     integer                        :: value
  end type KeyValueCI

!===============================================================================
! KEYVALUEII stores the (key,value) pair for a dictionary where the key is an
! integer and the value is an integer. Note that we need to store the key in
! addition to the value for collision resolution.
!===============================================================================

  type KeyValueII
     integer :: key
     integer :: value
  end type KeyValueII

!===============================================================================
! LISTKEYVALUECI stores a linked list of (key,value) pairs where the key is a
! character and the value is an integer
!===============================================================================

  type ListKeyValueCI
     type(ListKeyValueCI), pointer :: next => null()
     type(KeyValueCI)              :: data
  end type ListKeyValueCI

!===============================================================================
! LISTKEYVALUEII stores a linked list of (key,value) pairs where the key is a
! character and the value is an integer
!===============================================================================

  type ListKeyValueII
     type(ListKeyValueII), pointer :: next => null()
     type(KeyValueII)              :: data
  end type ListKeyValueII

!===============================================================================
! LISTREAL stores a linked list of real values. This is used for constructing a
! unionized energy grid.
!===============================================================================

  type ListReal
     type(ListReal), pointer :: next => null()
     real(8)                 :: data
  end type ListReal

!===============================================================================
! LISTINT stores a linked list of integer values.
!===============================================================================

  type ListInt
     type(ListInt), pointer :: next => null()
     integer                :: data
  end type ListInt

!===============================================================================
! HASHLISTCI - Since it's not possible to directly do an array of pointers, this
! derived type provides a pointer
!===============================================================================

  type HashListCI
     type(ListKeyValueCI), pointer :: list => null()
  end type HashListCI

!===============================================================================
! HASHLISTII - Since it's not possible to directly do an array of pointers, this
! derived type provides a pointer
!===============================================================================

  type HashListII
     type(ListKeyValueII), pointer :: list => null()
  end type HashListII

!===============================================================================
! DICTIONARYCI provides a dictionary data structure of (key,value) pairs where
! the keys are strings and values are integers.
!===============================================================================

  type DictionaryCI
     type(HashListCI), pointer :: table(:) => null()
  end type DictionaryCI

!===============================================================================
! DICTIONARYII provides a dictionary data structure of (key,value) pairs where
! the keys and values are both integers.
!===============================================================================

  type DictionaryII
     type(HashListII), pointer :: table(:) => null()
  end type DictionaryII

end module types
