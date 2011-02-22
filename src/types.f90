module types

  implicit none

!=====================================================================
! SURFACE type defines a first- or second-order surface that can be
! used to construct closed volumes (cells)
!=====================================================================

  type Surface
     integer :: uid
     integer :: type
     real(8), allocatable :: coeffs(:)
     integer, allocatable :: neighbor_pos(:)
     integer, allocatable :: neighbor_neg(:)
     integer :: bc
  end type Surface

!=====================================================================
! CELL defines a closed volume by its bounding surfaces
!=====================================================================

  type Cell
     integer :: uid
     integer :: n_items
     integer, allocatable :: boundary_list(:)
     integer :: material
  end type Cell

!=====================================================================
! NEUTRON describes the state of a neutron being transported through
! the geometry
!=====================================================================

  type Neutron
    integer :: uid     ! Unique ID
    real(8) :: xyz(3)  ! location
    real(8) :: uvw(3)  ! directional cosines
    real(8) :: E       ! energy
    integer :: IE      ! index on energy grid
    real(8) :: interp  ! interpolation factor for energy grid
    integer :: cell    ! current cell
    integer :: surface ! current surface
    real(8) :: wgt     ! particle weight
    logical :: alive   ! is particle alive?
  end type Neutron

!=====================================================================
! BANK is used for storing fission sites in criticality
! calculations. Since all the state information of a neutron is not
! needed, this type allows sites to be stored with less memory
!=====================================================================

  type Bank
     integer :: uid    ! Unique ID
     real(8) :: xyz(3) ! Location of bank particle
  end type Bank

!=====================================================================
! ISOTOPE describes an isotope, e.g. U-235, within a material. Note
! that two separate variables must be used for the same isotope in two
! different materials since they will generally have different
! densities
!=====================================================================

  type Isotope
     integer :: uid     ! unique identifier
     integer :: zaid    ! ZAID, e.g. 92235
     character(3) :: xs ! cross section identifier, e.g. 70c
     real(8) :: density ! density in atom/b-cm
  end type Isotope

!=====================================================================
! MATERIAL describes a material by its constituent isotopes
!=====================================================================

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

!=====================================================================
! EXTSOURCE describes an external source of neutrons for a
! fixed-source problem or for the starting source in a criticality
! problem
!=====================================================================

  type ExtSource
     integer :: type                    ! type of source, e.g. 'box' or 'cell'
     real(8), allocatable :: values(:)  ! values for particular source type
  end type ExtSource

!=====================================================================
! ACEREACTION contains the cross-section and secondary energy and
! angle distributions for a single reaction in a continuous-energy
! ACE-format table
!=====================================================================

  type AceReaction
     integer :: MT                     ! ENDF MT value
     real(8) :: Q_value                ! Reaction Q value
     real(8) :: TY                     ! Number of neutrons released
     integer :: IE                     ! Starting energy grid index
     real(8), allocatable :: sigma(:)  ! Cross section values
     logical :: has_angle_dist
     logical :: has_energy_dist

     ! Secondary angle distribution
     integer              :: adist_n_energy    ! # of incoming energies
     real(8), allocatable :: adist_energy(:)   ! incoming energy grid
     integer, allocatable :: adist_location(:) ! location of each table
     real(8), allocatable :: adist_data(:)     ! angular distribution data

     ! Secondary energy distribution
     integer :: edist_law                    ! secondary distribution law
     integer :: edist_n_interp               ! # of interpolation regions
     integer, allocatable :: edist_nbt(:)    ! ENDF interpolation parameters
     integer, allocatable :: edist_int(:)    ! ''
     integer :: edist_n_energy               ! # of incoming energies
     real(8), allocatable :: edist_energy(:) ! energy grid for law validity
     real(8), allocatable :: edist_pvalid(:) ! probability of law validity
     real(8), allocatable :: edist_data(:)   ! energy distribution data
  end type AceReaction

!=====================================================================
! ACECONTINUOUS contains all the data for an ACE-format
! continuous-energy cross section. The ACE format (A Compact ENDF
! format) is used in MCNP and several other Monte Carlo codes.
!=====================================================================

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

     integer :: nu_p_type
     real(8), allocatable :: nu_p_energy(:)
     real(8), allocatable :: nu_p_value(:)
     
     integer :: nu_t_type
     real(8), allocatable :: nu_t_energy(:)
     real(8), allocatable :: nu_t_value(:)

     real(8), allocatable :: nu_d_energy(:)
     real(8), allocatable :: nu_d_value(:)
     real(8), allocatable :: nu_d_precursor_const(:,:)
     real(8), allocatable :: nu_d_precursor_energy(:,:)
     real(8), allocatable :: nu_d_precursor_prob(:,:)
     integer :: n_reaction
     type(AceReaction), pointer :: reactions(:) => null()

  end type AceContinuous

!=====================================================================
! ACETHERMAL contains S(a,b) data for thermal neutron scattering,
! typically off of light isotopes such as water, graphite, Be, etc
!=====================================================================
     
  type AceThermal
     character(20) :: name
     real(8) :: awr
     real(8) :: temperature
     real(8), allocatable :: inelastic_e_in(:)
     real(8), allocatable :: inelastic_sigma(:) 
     real(8), allocatable :: inelastic_e_out(:,:)
     real(8), allocatable :: inelastic_mu_out(:,:)
     real(8), allocatable :: elastic_e_in(:)
     real(8), allocatable :: elastic_P(:)
     real(8), allocatable :: elastic_mu_out(:)
  end type AceThermal

!=====================================================================
! ACEDISTANGLE contains data for a tabular secondary angle
! distribution whether it be tabular or 32 equiprobable cosine bins
!=====================================================================

!!$  type AceDistAngle
!!$     integer :: n_energy
!!$     real(8), allocatable :: energy
!!$     integer, allocatable :: location
!!$     real(8), allocatable :: data
!!$  end type AceDistAngle

!=====================================================================
! ACEDISTTAB contains data for a tabular secondary energy distribution
! according to MCNP Law 4 (ENDF Law 1)
!=====================================================================

  type AceDistEnergy
     integer :: law
  end type AceDistEnergy

!=====================================================================
! XSDATA contains data read in from a SERPENT xsdata file
!=====================================================================

  type xsData
     character(10) :: alias
     character(10) :: id
     integer :: type
     integer :: zaid
     integer :: isomeric
     real(8) :: awr
     real(8) :: temp
     integer :: binary
     character(100) :: path
  end type xsData

!=====================================================================
! KEYVALUECI stores the (key,value) pair for a dictionary where the
! key is a string and the value is an integer. Note that we need to
! store the key in addition to the value for collision resolution.
!=====================================================================

  ! Key length for dictionary
  integer, parameter :: DICT_KEY_LENGTH = 20

  type KeyValueCI
     character(len=DICT_KEY_LENGTH) :: key
     integer                        :: value
  end type KeyValueCI

!=====================================================================
! KEYVALUEII stores the (key,value) pair for a dictionary where the
! key is an integer and the value is an integer. Note that we need to
! store the key in addition to the value for collision resolution.
!=====================================================================

  type KeyValueII
     integer :: key
     integer :: value
  end type KeyValueII

!=====================================================================
! LISTKEYVALUECI stores a linked list of (key,value) pairs where the
! key is a character and the value is an integer
!=====================================================================

  type ListKeyValueCI
     type(ListKeyValueCI), pointer :: next => null()
     type(KeyValueCI)              :: data
  end type ListKeyValueCI

!=====================================================================
! LISTKEYVALUEII stores a linked list of (key,value) pairs where the
! key is a character and the value is an integer
!=====================================================================

  type ListKeyValueII
     type(ListKeyValueII), pointer :: next => null()
     type(KeyValueII)              :: data
  end type ListKeyValueII

!=====================================================================
! LISTREAL stores a linked list of real values. This is used for
! constructing a unionized energy grid.
!=====================================================================

  type ListReal
     type(ListReal), pointer :: next => null()
     real(8)                 :: data
  end type ListReal

!=====================================================================
! LISTINT stores a linked list of integer values.
!=====================================================================

  type ListInt
     type(ListInt), pointer :: next => null()
     integer                :: data
  end type ListInt

!=====================================================================
! HASHLISTCI - Since it's not possible to directly do an array of
! pointers, this derived type provides a pointer
!=====================================================================

  type HashListCI
     type(ListKeyValueCI), pointer :: list => null()
  end type HashListCI

!=====================================================================
! HASHLISTII - Since it's not possible to directly do an array of
! pointers, this derived type provides a pointer
!=====================================================================

  type HashListII
     type(ListKeyValueII), pointer :: list => null()
  end type HashListII

!=====================================================================
! DICTIONARYCI provides a dictionary data structure of (key,value)
! pairs where the keys are strings and values are integers.
!=====================================================================

  type DictionaryCI
     type(HashListCI), pointer :: table(:) => null()
  end type DictionaryCI

!=====================================================================
! DICTIONARYII provides a dictionary data structure of (key,value)
! pairs where the keys and values are both integers.
!=====================================================================

  type DictionaryII
     type(HashListII), pointer :: table(:) => null()
  end type DictionaryII

end module types
