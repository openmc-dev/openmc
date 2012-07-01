module material_header

  implicit none

!===============================================================================
! MATERIAL describes a material by its constituent isotopes
!===============================================================================

  type Material
     integer              :: id              ! unique identifier
     integer              :: n_nuclides      ! number of nuclides
     character(12), allocatable :: names(:)  ! isotope names
     integer, allocatable :: nuclide(:)      ! index in nuclides array
     real(8)              :: density         ! total atom density in atom/b-cm
     real(8), allocatable :: atom_density(:) ! nuclide atom density in atom/b-cm

     ! S(a,b) data references
     logical       :: has_sab_table = .false.
     character(12) :: sab_name        ! name of S(a,b) table
     integer       :: sab_table   = 0 ! index in sab_tables
     integer       :: sab_nuclide = 0 ! index of nuclide which has S(a,b) table
  end type Material

end module material_header
