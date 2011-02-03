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
    integer  :: uid     ! Unique ID
    real(8)  :: xyz(3)  ! location
    real(8)  :: uvw(3)  ! directional cosines
    integer  :: cell    ! current cell
    integer  :: surface ! current surface
    real(8)  :: wgt     ! particle weight
    logical  :: alive   ! is particle alive?
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
     character(10), allocatable :: names(:)
     integer, allocatable :: isotopes(:)
     integer, allocatable :: table(:)
     real(8), allocatable :: atom_percent(:)
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
     integer :: MT
     real(8) :: Q_value
     real(8) :: TY
     integer :: energy_index
     real(8), allocatable :: sigma(:)
     real(8), allocatable :: ang_cos(:,:)
     real(8), allocatable :: ang_pdf(:,:)
     real(8), allocatable :: ang_cdf(:,:)
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
     type(AceReaction), pointer :: reactions(:)

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
! LISTDATA Data stored in a linked list. In this case, we store the
! (key,value) pair for a dictionary. Note that we need to store the
! key in addition to the value for collision resolution.
!=====================================================================

  ! Key length for dictionary
  integer, parameter :: DICT_KEY_LENGTH = 20

  type ListData
     character(len=DICT_KEY_LENGTH) :: key
     integer                        :: value
  end type ListData

!=====================================================================
! LINKEDLIST stores a simple linked list
!=====================================================================

  type LinkedList
     type(LinkedList), pointer :: next
     type(ListData)            :: data
  end type LinkedList

!=====================================================================
! LINKEDLISTGRID stores a sorted list of energies for the unionized
! energy grid as a linked list
!=====================================================================

  type LinkedListGrid
     type(LinkedListGrid), pointer :: next
     real(8)                   :: energy
  end type LinkedListGrid

!=====================================================================
! HASHLIST - Since it's not possible to directly do an array of
! pointers, this derived type provides a pointer
!=====================================================================

  type HashList
     type(LinkedList), pointer :: list
  end type HashList

!=====================================================================
! DICTIONARY provides a dictionary data structure of (key,value) pairs
!=====================================================================

  type Dictionary
     type(HashList), pointer :: table(:)
  end type Dictionary

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

end module types
