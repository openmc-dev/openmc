module material_header

  implicit none

!===============================================================================
! MATERIAL describes a material by its constituent isotopes
!===============================================================================

  type Material
     integer              :: uid             ! unique identifier
     integer              :: n_nuclides      ! number of nuclides
     character(10), allocatable :: names(:)  ! isotope names
     integer, allocatable :: xsdata(:)       ! index in xsdata list
     integer, allocatable :: nuclide(:)      ! index in nuclides array
     real(8)              :: density         ! total atom density in atom/b-cm
     real(8), allocatable :: atom_density(:) ! nuclide atom density in atom/b-cm
     real(8), allocatable :: atom_percent(:) ! atom/weight percent (negative for weight)
     real(8), allocatable :: total_xs(:)     ! macroscopic cross-section

     ! S(a,b) data references
     integer :: sab_table      ! index in sab_tables
     integer :: sab_nuclide    ! index of nuclide which has S(a,b) table
  end type Material

end module material_header
