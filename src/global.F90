module global

  use bank_header,      only: Bank
  use cmfd_header
  use constants
  use dict_header,      only: DictCharInt, DictIntInt
  use geometry_header,  only: Cell, Universe, Lattice, LatticeContainer
  use material_header,  only: Material
  use mesh_header,      only: RegularMesh
  use mgxs_header,      only: Mgxs, MgxsContainer
  use nuclide_header
  use plot_header,      only: ObjectPlot
  use sab_header,       only: SAlphaBeta
  use set_header,       only: SetInt
  use surface_header,   only: SurfaceContainer
  use source_header,    only: SourceDistribution
  use tally_header,     only: TallyObject, TallyResult
  use trigger_header,   only: KTrigger
  use timer_header,     only: Timer
  use volume_header,    only: VolumeCalculation

#ifdef MPIF08
  use mpi_f08
#endif

  implicit none

  ! ============================================================================
  ! GEOMETRY-RELATED VARIABLES

  ! Main arrays
  type(Cell),             allocatable, target :: cells(:)
  type(Universe),         allocatable, target :: universes(:)
  type(LatticeContainer), allocatable, target :: lattices(:)
  type(SurfaceContainer), allocatable, target :: surfaces(:)
  type(Material),         allocatable, target :: materials(:)
  type(ObjectPlot),       allocatable, target :: plots(:)

  type(VolumeCalculation), allocatable :: volume_calcs(:)

  ! Size of main arrays
  integer :: n_cells     ! # of cells
  integer :: n_universes ! # of universes
  integer :: n_lattices  ! # of lattices
  integer :: n_surfaces  ! # of surfaces
  integer :: n_materials ! # of materials
  integer :: n_plots     ! # of plots

  ! These dictionaries provide a fast lookup mechanism -- the key is the
  ! user-specified identifier and the value is the index in the corresponding
  ! array
  type(DictIntInt) :: cell_dict
  type(DictIntInt) :: universe_dict
  type(DictIntInt) :: lattice_dict
  type(DictIntInt) :: surface_dict
  type(DictIntInt) :: material_dict
  type(DictIntInt) :: mesh_dict
  type(DictIntInt) :: tally_dict
  type(DictIntInt) :: plot_dict

  ! Number of lost particles
  integer :: n_lost_particles

  ! ============================================================================
  ! ENERGY TREATMENT RELATED VARIABLES
  logical :: run_CE = .true.  ! Run in CE mode?

  ! ============================================================================
  ! CROSS SECTION RELATED VARIABLES NEEDED REGARDLESS OF CE OR MG

  integer :: n_nuclides_total ! Number of nuclide cross section tables

  ! Cross section caches
  type(NuclideMicroXS), allocatable :: micro_xs(:)  ! Cache for each nuclide
  type(MaterialMacroXS)             :: material_xs  ! Cache for current material

  ! Dictionaries to look up cross sections and listings
  type(DictCharInt) :: nuclide_dict

  ! ============================================================================
  ! CONTINUOUS-ENERGY CROSS SECTION RELATED VARIABLES

  ! Cross section arrays
  type(Nuclide), allocatable, target :: nuclides(:)    ! Nuclide cross-sections
  type(SAlphaBeta), allocatable, target :: sab_tables(:)  ! S(a,b) tables

  integer :: n_sab_tables     ! Number of S(a,b) thermal scattering tables

  ! Minimum/maximum energies
  real(8) :: energy_min_neutron = ZERO
  real(8) :: energy_max_neutron = INFINITY

  ! Dictionaries to look up cross sections and listings
  type(DictCharInt) :: sab_dict

  ! Unreoslved resonance probablity tables
  logical :: urr_ptables_on = .true.

  ! What to assume for expanding natural elements
  integer :: default_expand = ENDF_BVII1

  ! Default temperature and method for choosing temperatures
  integer :: temperature_method = TEMPERATURE_NEAREST
  logical :: temperature_multipole = .false.
  real(8) :: temperature_tolerance = 10.0_8
  real(8) :: temperature_default = 293.6_8

  ! ============================================================================
  ! MULTI-GROUP CROSS SECTION RELATED VARIABLES

  ! Cross section arrays
  type(MgxsContainer), allocatable, target :: nuclides_MG(:)

  ! Cross section caches
  type(MgxsContainer), target, allocatable :: macro_xs(:)

  ! Number of energy groups
  integer :: energy_groups

  ! Energy group structure
  real(8), allocatable :: energy_bins(:)

  ! Midpoint of the energy group structure
  real(8), allocatable :: energy_bin_avg(:)

  ! Maximum Data Order
  integer :: max_order

  ! Whether or not to convert Legendres to tabulars
  logical :: legendre_to_tabular = .True.

  ! Number of points to use in the Legendre to tabular conversion
  integer :: legendre_to_tabular_points = 33

  ! ============================================================================
  ! ELEMENT AND NUCLIDE RELATED VARIABLES

  ! List of all possible natural element expansion isotopes
  character (len=5), dimension (NUM_NATURAL_NUCLIDES), parameter :: &
       natural_nuclides = [character(len=5) :: &
       'H1'   , 'H2'   , 'He3'  , 'He4'  , 'Li6'  , 'Li7'  , 'Be9'  , 'B10'  , &
       'B11'  , 'C0'   , 'N14'  , 'N15'  , 'O16'  , 'O17'  , 'O18'  , 'O16'  , &
       'O16'  , 'O17'  , 'F19'  , 'Ne20' , 'Ne21' , 'Ne22' , 'Na23' , 'Mg24' , &
       'Mg25' , 'Mg26' , 'Al27' , 'Si28' , 'Si29' , 'Si30' , 'P31'  , 'S32'  , &
       'S33'  , 'S34'  , 'S35'  , 'Cl35' , 'Cl37' , 'Ar36' , 'Ar38' , 'Ar40' , &
       'K39'  , 'K40'  , 'K41'  , 'Ca40' , 'Ca42' , 'Ca43' , 'Ca44' , 'Ca46' , &
       'Ca48' , 'Sc45' , 'Ti46' , 'Ti47' , 'Ti48' , 'Ti49' , 'Ti50' , 'V0'   , &
       'V50'  , 'V51'  , 'Cr50' , 'Cr52' , 'Cr53' , 'Cr54' , 'Mn55' , 'Fe54' , &
       'Fe56' , 'Fe57' , 'Fe58' , 'Co59' , 'Ni58' , 'Ni60' , 'Ni61' , 'Ni62' , &
       'Ni64' , 'Cu63' , 'Cu65' , 'Zn0'  , 'Zn64' , 'Zn66' , 'Zn67' , 'Zn68' , &
       'Zn70' , 'Ga0'  , 'Ga69' , 'Ga71' , 'Ge70' , 'Ge72' , 'Ge73' , 'Ge74' , &
       'Ge76' , 'As75' , 'Se74' , 'Se76' , 'Se77' , 'Se78' , 'Se80' , 'Se82' , &
       'Br79' , 'Br81' , 'Kr78' , 'Kr80' , 'Kr82' , 'Kr83' , 'Kr84' , 'Kr86' , &
       'Rb85' , 'Rb87' , 'Sr84' , 'Sr86' , 'Sr87' , 'Sr88' , 'Y89'  , 'Zr90' , &
       'Zr91' , 'Zr92' , 'Zr94' , 'Zr96' , 'Nb93' , 'Mo92' , 'Mo94' , 'Mo95' , &
       'Mo96' , 'Mo97' , 'Mo98' , 'Mo100', 'Ru96' , 'Ru98' , 'Ru99' , 'Ru100', &
       'Ru101', 'Ru102', 'Ru104', 'Rh103', 'Pd102', 'Pd104', 'Pd105', 'Pd106', &
       'Pd108', 'Pd110', 'Ag107', 'Ag109', 'Cd106', 'Cd108', 'Cd110', 'Cd111', &
       'Cd112', 'Cd113', 'Cd114', 'Cd116', 'In113', 'In115', 'Sn112', 'Sn114', &
       'Sn115', 'Sn116', 'Sn117', 'Sn118', 'Sn119', 'Sn120', 'Sn122', 'Sn124', &
       'Sb121', 'Sb123', 'Te120', 'Te122', 'Te123', 'Te124', 'Te125', 'Te126', &
       'Te128', 'Te130', 'I127' , 'Xe124', 'Xe126', 'Xe128', 'Xe129', 'Xe130', &
       'Xe131', 'Xe132', 'Xe134', 'Xe136', 'Cs133', 'Ba130', 'Ba132', 'Ba134', &
       'Ba135', 'Ba136', 'Ba137', 'Ba138', 'La138', 'La139', 'Ce136', 'Ce138', &
       'Ce140', 'Ce142', 'Pr141', 'Nd142', 'Nd143', 'Nd144', 'Nd145', 'Nd146', &
       'Nd148', 'Nd150', 'Sm144', 'Sm147', 'Sm148', 'Sm149', 'Sm150', 'Sm152', &
       'Sm154', 'Eu151', 'Eu153', 'Gd152', 'Gd154', 'Gd155', 'Gd156', 'Gd157', &
       'Gd158', 'Gd160', 'Tb159', 'Dy156', 'Dy158', 'Dy160', 'Dy161', 'Dy162', &
       'Dy163', 'Dy164', 'Ho165', 'Er162', 'Er164', 'Er166', 'Er167', 'Er168', &
       'Er170', 'Tm169', 'Yb168', 'Yb170', 'Yb171', 'Yb172', 'Yb173', 'Yb174', &
       'Yb176', 'Lu175', 'Lu176', 'Hf174', 'Hf176', 'Hf177', 'Hf178', 'Hf179', &
       'Hf180', 'Ta181', 'Ta180', 'Ta181', 'W182' , 'W183' , 'W184' , 'W186' , &
       'W180' , 'W182' , 'W183' , 'W184' , 'W186' , 'Re185', 'Re187', 'Os0'  , &
       'Os184', 'Os186', 'Os187', 'Os188', 'Os189', 'Os190', 'Os192', 'Ir191', &
       'Ir193', 'Pt0'  , 'Pt190', 'Pt192', 'Pt194', 'Pt195', 'Pt196', 'Pt198', &
       'Au197', 'Hg196', 'Hg198', 'Hg199', 'Hg200', 'Hg201', 'Hg202', 'Hg204', &
       'Tl0'  , 'Tl203', 'Tl205', 'Pb204', 'Pb206', 'Pb207', 'Pb208', 'Bi209', &
       'Th232', 'Pa231', 'U234' , 'U235' , 'U238']

  character (len=2), dimension (NUM_NATURAL_ELEMENTS), parameter :: &
       natural_elements = [character(len=2) :: &
       'h' , 'he', 'li', 'be', 'b' , 'c' , 'n' , 'o' , 'f' , 'ne', &
       'na', 'mg', 'al', 'si', 'p' , 's' , 'cl', 'ar', 'k' , 'ca', &
       'sc', 'ti', 'v' , 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', &
       'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y' , 'zr', &
       'nb', 'mo', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', &
       'te', 'i' , 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'sm', &
       'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', &
       'ta', 'w' , 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', &
       'bi', 'th', 'pa', 'u']

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_ENDF_BVII0_start = (/ &
       1   , 3   , 5   , 7   , 8   , 10  , 11  , 17  , 19  , 20  , &
       23  , 24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , &
       50  , 51  , 56  , 59  , 63  , 64  , 68  , 69  , 74  , 76  , &
       83  , 85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , &
       117 , 118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , &
       163 , 171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , &
       210 , 212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , &
       250 , 253 , 262 , 265 , 272 , 275 , 281 , 282 , 290 , 292 , &
       296 , 297 , 298 , 299 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_ENDF_BVII0_end = (/ &
       3   , 5   , 7   , 8   , 10  , 11  , 13  , 19  , 20  , 23  , &
       24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , 50  , &
       51  , 56  , 57  , 63  , 64  , 68  , 69  , 74  , 76  , 77  , &
       85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , 117 , &
       118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , 163 , &
       171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , 210 , &
       212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , 250 , &
       251 , 257 , 264 , 272 , 274 , 281 , 282 , 289 , 292 , 296 , &
       297 , 298 , 299 , 302 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_ENDF_BVII1_start = (/ &
       1   , 3   , 5   , 7   , 8   , 10  , 11  , 17  , 19  , 20  , &
       23  , 24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , &
       50  , 51  , 57  , 59  , 63  , 64  , 68  , 69  , 74  , 77  , &
       83  , 85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , &
       117 , 118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , &
       163 , 171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , &
       210 , 212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , &
       251 , 257 , 262 , 265 , 272 , 275 , 281 , 282 , 290 , 292 , &
       296 , 297 , 298 , 299 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_ENDF_BVII1_end = (/ &
       3   , 5   , 7   , 8   , 10  , 11  , 13  , 19  , 20  , 23  , &
       24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , 50  , &
       51  , 56  , 59  , 63  , 64  , 68  , 69  , 74  , 76  , 82  , &
       85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , 117 , &
       118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , 163 , &
       171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , 210 , &
       212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , 250 , &
       253 , 262 , 264 , 272 , 274 , 281 , 282 , 289 , 292 , 296 , &
       297 , 298 , 299 , 302 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JEFF_311_start = (/ &
       1   , 3   , 5   , 7   , 8   , 10  , 11  , 17  , 19  , 20  , &
       23  , 24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , &
       50  , 51  , 56  , 59  , 63  , 64  , 68  , 69  , 74  , 76  , &
       82  , 85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , &
       117 , 118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , &
       163 , 171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , &
       210 , 212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , &
       250 , 253 , 262 , 264 , 272 , 274 , 281 , 282 , 289 , 292 , &
       296 , 297 , 298 , 299 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JEFF_311_end = (/ &
       3   , 5   , 7   , 8   , 10  , 11  , 13  , 19  , 20  , 23  , &
       24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , 50  , &
       51  , 56  , 57  , 63  , 64  , 68  , 69  , 74  , 76  , 77  , &
       83  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , 117 , &
       118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , 163 , &
       171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , 210 , &
       212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , 250 , &
       251 , 257 , 264 , 265 , 274 , 275 , 282 , 289 , 290 , 296 , &
       297 , 298 , 299 , 302 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JEFF_312_start = (/ &
       1   , 3   , 5   , 7   , 8   , 10  , 11  , 17  , 19  , 20  , &
       23  , 24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , &
       50  , 51  , 57  , 59  , 63  , 64  , 68  , 69  , 74  , 76  , &
       82  , 85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , &
       117 , 118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , &
       163 , 171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , &
       210 , 212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , &
       250 , 253 , 262 , 264 , 272 , 274 , 281 , 282 , 289 , 292 , &
       296 , 297 , 298 , 299 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JEFF_312_end = (/ &
       3   , 5   , 7   , 8   , 10  , 11  , 13  , 19  , 20  , 23  , &
       24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , 50  , &
       51  , 56  , 59  , 63  , 64  , 68  , 69  , 74  , 76  , 77  , &
       83  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , 117 , &
       118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , 163 , &
       171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , 210 , &
       212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , 250 , &
       251 , 257 , 264 , 265 , 274 , 275 , 282 , 289 , 290 , 296 , &
       297 , 298 , 299 , 302 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JEFF_32_start = (/ &
       1   , 3   , 5   , 7   , 8   , 10  , 11  , 13  , 19  , 20  , &
       23  , 24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , &
       50  , 51  , 56  , 59  , 63  , 64  , 68  , 69  , 74  , 77  , &
       83  , 85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , &
       117 , 118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , &
       163 , 171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , &
       210 , 212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , &
       251 , 257 , 262 , 265 , 272 , 275 , 281 , 282 , 290 , 292 , &
       296 , 297 , 298 , 299 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JEFF_32_end = (/ &
       3   , 5   , 7   , 8   , 10  , 11  , 13  , 16  , 20  , 23  , &
       24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , 50  , &
       51  , 56  , 57  , 63  , 64  , 68  , 69  , 74  , 76  , 82  , &
       85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , 117 , &
       118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , 163 , &
       171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , 210 , &
       212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , 250 , &
       253 , 262 , 264 , 272 , 274 , 281 , 282 , 289 , 292 , 296 , &
       297 , 298 , 299 , 302 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JENDL_32_start = (/ &
       1   , 3   , 5   , 7   , 8   , 10  , 11  , 16  , 19  , 20  , &
       23  , 24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , &
       50  , 51  , 56  , 59  , 63  , 64  , 68  , 69  , 74  , 77  , &
       83  , 85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , &
       117 , 118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , &
       163 , 171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , &
       210 , 212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , &
       250 , 253 , 262 , 265 , 272 , 275 , 281 , 282 , 290 , 292 , &
       296 , 297 , 298 , 299 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JENDL_32_end = (/ &
       3   , 5   , 7   , 8   , 10  , 11  , 13  , 17  , 20  , 23  , &
       24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , 50  , &
       51  , 56  , 57  , 63  , 64  , 68  , 69  , 74  , 76  , 82  , &
       85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , 117 , &
       118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , 163 , &
       171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , 210 , &
       212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , 250 , &
       251 , 257 , 264 , 272 , 274 , 281 , 282 , 289 , 292 , 296 , &
       297 , 298 , 299 , 302 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JENDL_33_start = (/ &
       1   , 3   , 5   , 7   , 8   , 10  , 11  , 16  , 19  , 20  , &
       23  , 24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , &
       50  , 51  , 56  , 59  , 63  , 64  , 68  , 69  , 74  , 77  , &
       83  , 85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , &
       117 , 118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , &
       163 , 171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , &
       210 , 212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , &
       250 , 253 , 262 , 265 , 272 , 275 , 281 , 282 , 290 , 292 , &
       296 , 297 , 298 , 299 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JENDL_33_end = (/ &
       3   , 5   , 7   , 8   , 10  , 11  , 13  , 17  , 20  , 23  , &
       24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , 50  , &
       51  , 56  , 57  , 63  , 64  , 68  , 69  , 74  , 76  , 82  , &
       85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , 117 , &
       118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , 163 , &
       171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , 210 , &
       212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , 250 , &
       251 , 257 , 264 , 272 , 274 , 281 , 282 , 289 , 292 , 296 , &
       297 , 298 , 299 , 302 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JENDL_40_start = (/ &
       1   , 3   , 5   , 7   , 8   , 10  , 11  , 16  , 19  , 20  , &
       23  , 24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , &
       50  , 51  , 57  , 59  , 63  , 64  , 68  , 69  , 74  , 77  , &
       83  , 85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , &
       117 , 118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , &
       163 , 171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , &
       210 , 212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , &
       250 , 257 , 262 , 265 , 272 , 275 , 281 , 282 , 290 , 292 , &
       296 , 297 , 298 , 299 /)

  integer, dimension(NUM_NATURAL_ELEMENTS) :: &
       natural_nuclides_JENDL_40_end = (/ &
       3   , 5   , 7   , 8   , 10  , 11  , 13  , 17  , 20  , 23  , &
       24  , 27  , 28  , 31  , 32  , 36  , 38  , 41  , 44  , 50  , &
       51  , 56  , 59  , 63  , 64  , 68  , 69  , 74  , 76  , 82  , &
       85  , 90  , 90  , 97  , 99  , 105 , 107 , 111 , 112 , 117 , &
       118 , 125 , 132 , 133 , 139 , 141 , 149 , 151 , 161 , 163 , &
       171 , 172 , 181 , 182 , 189 , 191 , 195 , 196 , 203 , 210 , &
       212 , 219 , 220 , 227 , 228 , 234 , 235 , 242 , 244 , 250 , &
       251 , 262 , 264 , 272 , 274 , 281 , 282 , 289 , 292 , 296 , &
       297 , 298 , 299 , 302 /)

  real(8), dimension (NUM_NATURAL_NUCLIDES) :: natural_nuclides_mf = (/ &
       0.999885_8  , 0.000115_8  , 0.00000134_8, 0.99999866_8, &
       0.0759_8    , 0.9241_8    , ONE         , 0.199_8     , &
       0.801_8     , ONE         , 0.99636_8   , 0.00364_8   , &
       0.99757_8   , 0.00038_8   , 0.00205_8   , ONE         , &
       0.99962_8   , 0.00038_8   , ONE         , 0.9048_8    , &
       0.0027_8    , 0.0925_8    , ONE         , 0.7899_8    , &
       0.1000_8    , 0.1101_8    , ONE         , 0.92223_8   , &
       0.04685_8   , 0.03092_8   , ONE         , 0.9499_8    , &
       0.0075_8    , 0.0425_8    , 0.0001_8    , 0.7576_8    , &
       0.2424_8    , 0.003336_8  , 0.000629_8  , 0.996035_8  , &
       0.932581_8  , 0.000117_8  , 0.067302_8  , 0.96941_8   , &
       0.00647_8   , 0.00135_8   , 0.02086_8   , 0.00004_8   , &
       0.00187_8   , ONE         , 0.0825_8    , 0.0744_8    , &
       0.7372_8    , 0.0541_8    , 0.0518_8    , ONE         , &
       0.0025_8    , 0.9975_8    , 0.04345_8   , 0.83789_8   , &
       0.09501_8   , 0.02365_8   , ONE         , 0.05845_8   , &
       0.91754_8   , 0.02119_8   , 0.00282_8   , ONE         , &
       0.68077_8   , 0.26223_8   , 0.011399_8  , 0.036346_8  , &
       0.009255_8  , 0.6915_8    , 0.3085_8    , ONE         , &
       0.4917_8    , 0.2773_8    , 0.0404_8    , 0.1845_8    , &
       0.0061_8    , ONE         , 0.60108_8   , 0.39892_8   , &
       0.2057_8    , 0.2745_8    , 0.0775_8    , 0.3650_8    , &
       0.0773_8    , ONE         , 0.0089_8    , 0.0937_8    , &
       0.0763_8    , 0.2377_8    , 0.4961_8    , 0.0873_8    , &
       0.5069_8    , 0.4931_8    , 0.00355_8   , 0.02286_8   , &
       0.11593_8   , 0.11500_8   , 0.56987_8   , 0.17279_8   , &
       0.7217_8    , 0.2783_8    , 0.0056_8    , 0.0986_8    , &
       0.0700_8    , 0.8258_8    , ONE         , 0.5145_8    , &
       0.1122_8    , 0.1715_8    , 0.1738_8    , 0.0280_8    , &
       ONE         , 0.1453_8    , 0.0915_8    , 0.1584_8    , &
       0.1667_8    , 0.0960_8    , 0.2439_8    , 0.0982_8    , &
       0.0554_8    , 0.0187_8    , 0.1276_8    , 0.1260_8    , &
       0.1706_8    , 0.3155_8    , 0.1862_8    , ONE         , &
       0.0102_8    , 0.1114_8    , 0.2233_8    , 0.2733_8    , &
       0.2646_8    , 0.1172_8    , 0.51839_8   , 0.48161_8   , &
       0.0125_8    , 0.0089_8    , 0.1249_8    , 0.1280_8    , &
       0.2413_8    , 0.1222_8    , 0.2873_8    , 0.0749_8    , &
       0.0429_8    , 0.9571_8    , 0.0097_8    , 0.0066_8    , &
       0.0034_8    , 0.1454_8    , 0.0768_8    , 0.2422_8    , &
       0.0859_8    , 0.3258_8    , 0.0463_8    , 0.0579_8    , &
       0.5721_8    , 0.4279_8    , 0.0009_8    , 0.0255_8    , &
       0.0089_8    , 0.0474_8    , 0.0707_8    , 0.1884_8    , &
       0.3174_8    , 0.3408_8    , ONE         , 0.000952_8  , &
       0.000890_8  , 0.019102_8  , 0.264006_8  , 0.040710_8  , &
       0.212324_8  , 0.269086_8  , 0.104357_8  , 0.088573_8  , &
       ONE         , 0.00106_8   , 0.00101_8   , 0.02417_8   , &
       0.06592_8   , 0.07854_8   , 0.11232_8   , 0.71698_8   , &
       0.0008881_8 , 0.9991119_8 , 0.00185_8   , 0.00251_8   , &
       0.88450_8   , 0.11114_8   , ONE         , 0.27152_8   , &
       0.12174_8   , 0.23798_8   , 0.08293_8   , 0.17189_8   , &
       0.05756_8   , 0.05638_8   , 0.0307_8    , 0.1499_8    , &
       0.1124_8    , 0.1382_8    , 0.0738_8    , 0.2675_8    , &
       0.2275_8    , 0.4781_8    , 0.5219_8    , 0.0020_8    , &
       0.0218_8    , 0.1480_8    , 0.2047_8    , 0.1565_8    , &
       0.2484_8    , 0.2186_8    , ONE         , 0.00056_8   , &
       0.00095_8   , 0.02329_8   , 0.18889_8   , 0.25475_8   , &
       0.24896_8   , 0.28260_8   , ONE         , 0.00139_8   , &
       0.01601_8   , 0.33503_8   , 0.22869_8   , 0.26978_8   , &
       0.14910_8   , ONE         , 0.00123_8   , 0.02982_8   , &
       0.1409_8    , 0.2168_8    , 0.16103_8   , 0.32026_8   , &
       0.12996_8   , 0.97401_8   , 0.02599_8   , 0.0016_8    , &
       0.0526_8    , 0.1860_8    , 0.2728_8    , 0.1362_8    , &
       0.3508_8    , ONE         , 0.0001201_8 , 0.9998799_8 , &
       0.2662_8    , 0.1431_8    , 0.3064_8    , 0.2843_8    , &
       0.0012_8    , 0.2650_8    , 0.1431_8    , 0.3064_8    , &
       0.2843_8    , 0.3740_8    , 0.6260_8    , ONE         , &
       0.0002_8    , 0.0159_8    , 0.0196_8    , 0.1324_8    , &
       0.1615_8    , 0.2626_8    , 0.4078_8    , 0.373_8     , &
       0.627_8     , ONE         , 0.00012_8   , 0.00782_8   , &
       0.3286_8    , 0.3378_8    , 0.2521_8    , 0.07356_8   , &
       ONE         , 0.0015_8    , 0.0997_8    , 0.1687_8    , &
       0.2310_8    , 0.1318_8    , 0.2986_8    , 0.0687_8    , &
       ONE         , 0.2952_8    , 0.7048_8    , 0.014_8     , &
       0.241_8     , 0.221_8     , 0.524_8     , ONE         , &
       ONE         , ONE         , 0.000054_8  , 0.007204_8  , &
       0.992742_8 /)

  ! ============================================================================
  ! TALLY-RELATED VARIABLES

  type(RegularMesh), allocatable, target :: meshes(:)
  type(TallyObject),    allocatable, target :: tallies(:)
  integer, allocatable :: matching_bins(:)
  real(8), allocatable :: filter_weights(:)

  ! Pointers for different tallies
  type(TallyObject), pointer :: user_tallies(:) => null()
  type(TallyObject), pointer :: cmfd_tallies(:) => null()

  ! Starting index (minus 1) in tallies for each tally group
  integer :: i_user_tallies = -1
  integer :: i_cmfd_tallies = -1

  ! Active tally lists
  type(SetInt) :: active_analog_tallies
  type(SetInt) :: active_tracklength_tallies
  type(SetInt) :: active_current_tallies
  type(SetInt) :: active_collision_tallies
  type(SetInt) :: active_tallies
!$omp threadprivate(active_analog_tallies, active_tracklength_tallies, &
!$omp&              active_current_tallies, active_collision_tallies, &
!$omp&              active_tallies)

  ! Global tallies
  !   1) collision estimate of k-eff
  !   2) absorption estimate of k-eff
  !   3) track-length estimate of k-eff
  !   4) leakage fraction

  type(TallyResult), allocatable, target :: global_tallies(:)

  ! It is possible to protect accumulate operations on global tallies by using
  ! an atomic update. However, when multiple threads accumulate to the same
  ! global tally, it can cause a higher cache miss rate due to
  ! invalidation. Thus, we use threadprivate variables to accumulate global
  ! tallies and then reduce at the end of a generation.
  real(8) :: global_tally_collision   = ZERO
  real(8) :: global_tally_absorption  = ZERO
  real(8) :: global_tally_tracklength = ZERO
  real(8) :: global_tally_leakage     = ZERO
!$omp threadprivate(global_tally_collision, global_tally_absorption, &
!$omp&              global_tally_tracklength, global_tally_leakage)

  integer :: n_meshes       = 0 ! # of structured meshes
  integer :: n_user_meshes  = 0 ! # of structured user meshes
  integer :: n_tallies      = 0 ! # of tallies
  integer :: n_user_tallies = 0 ! # of user tallies

  ! Normalization for statistics
  integer :: n_realizations = 0 ! # of independent realizations
  real(8) :: total_weight       ! total starting particle weight in realization

  ! Flag for turning tallies on
  logical :: tallies_on = .false.
  logical :: active_batches = .false.

  ! Assume all tallies are spatially distinct
  logical :: assume_separate = .false.

  ! Use confidence intervals for results instead of standard deviations
  logical :: confidence_intervals = .false.

  ! ============================================================================
  ! EIGENVALUE SIMULATION VARIABLES

  integer(8) :: n_particles = 0   ! # of particles per generation
  integer    :: n_batches         ! # of batches
  integer    :: n_inactive        ! # of inactive batches
  integer    :: n_active          ! # of active batches
  integer    :: gen_per_batch = 1 ! # of generations per batch
  integer    :: current_batch = 0 ! current batch
  integer    :: current_gen   = 0 ! current generation within a batch
  integer    :: overall_gen   = 0 ! overall generation in the run

  ! ============================================================================
  ! TALLY PRECISION TRIGGER VARIABLES

  integer        :: n_max_batches             ! max # of batches
  integer        :: n_batch_interval = 1      ! batch interval for triggers
  logical        :: pred_batches = .false.    ! predict batches for triggers
  logical        :: trigger_on = .false.      ! flag for turning triggers on/off
  type(KTrigger) :: keff_trigger              ! trigger for k-effective
  logical :: satisfy_triggers = .false.       ! whether triggers are satisfied

  ! External source
  type(SourceDistribution), allocatable :: external_source(:)

  ! Source and fission bank
  type(Bank), allocatable, target :: source_bank(:)
  type(Bank), allocatable, target :: fission_bank(:)
#ifdef _OPENMP
  type(Bank), allocatable, target :: master_fission_bank(:)
#endif
  integer(8) :: n_bank       ! # of sites in fission bank
  integer(8) :: work         ! number of particles per processor
  integer(8), allocatable :: work_index(:) ! starting index in source bank for each process
  integer(8) :: current_work ! index in source bank of current history simulated

  ! Temporary k-effective values
  real(8), allocatable :: k_generation(:) ! single-generation estimates of k
  real(8) :: keff = ONE       ! average k over active batches
  real(8) :: keff_std         ! standard deviation of average k
  real(8) :: k_col_abs = ZERO ! sum over batches of k_collision * k_absorption
  real(8) :: k_col_tra = ZERO ! sum over batches of k_collision * k_tracklength
  real(8) :: k_abs_tra = ZERO ! sum over batches of k_absorption * k_tracklength
  real(8) :: k_combined(2)    ! combined best estimate of k-effective

  ! Shannon entropy
  logical :: entropy_on = .false.
  real(8), allocatable :: entropy(:)         ! shannon entropy at each generation
  real(8), allocatable :: entropy_p(:,:,:,:) ! % of source sites in each cell
  type(RegularMesh), pointer :: entropy_mesh

  ! Uniform fission source weighting
  logical :: ufs = .false.
  type(RegularMesh), pointer :: ufs_mesh => null()
  real(8), allocatable :: source_frac(:,:,:,:)

  ! Write source at end of simulation
  logical :: source_separate = .false.
  logical :: source_write = .true.
  logical :: source_latest = .false.

  ! ============================================================================
  ! PARALLEL PROCESSING VARIABLES

  ! The defaults set here for the number of processors, rank, and master and
  ! mpi_enabled flag are for when MPI is not being used at all, i.e. a serial
  ! run. In this case, these variables are still used at times.

  integer :: n_procs     = 1       ! number of processes
  integer :: rank        = 0       ! rank of process
  logical :: master      = .true.  ! master process?
  logical :: mpi_enabled = .false. ! is MPI in use and initialized?
  integer :: mpi_err               ! MPI error code
#ifdef MPIF08
  type(MPI_Datatype) :: MPI_BANK
  type(MPI_Datatype) :: MPI_TALLYRESULT
#else
  integer :: MPI_BANK              ! MPI datatype for fission bank
  integer :: MPI_TALLYRESULT       ! MPI datatype for TallyResult
#endif

#ifdef _OPENMP
  integer :: n_threads = NONE      ! number of OpenMP threads
  integer :: thread_id             ! ID of a given thread
#endif

  ! No reduction at end of batch
  logical :: reduce_tallies = .true.

  ! ============================================================================
  ! TIMING VARIABLES

  type(Timer) :: time_total         ! timer for total run
  type(Timer) :: time_initialize    ! timer for initialization
  type(Timer) :: time_read_xs       ! timer for reading cross sections
  type(Timer) :: time_unionize      ! timer for material xs-energy grid union
  type(Timer) :: time_bank          ! timer for fission bank synchronization
  type(Timer) :: time_bank_sample   ! timer for fission bank sampling
  type(Timer) :: time_bank_sendrecv ! timer for fission bank SEND/RECV
  type(Timer) :: time_tallies       ! timer for accumulate tallies
  type(Timer) :: time_inactive      ! timer for inactive batches
  type(Timer) :: time_active        ! timer for active batches
  type(Timer) :: time_transport     ! timer for transport only
  type(Timer) :: time_finalize      ! timer for finalization

  ! ===========================================================================
  ! VARIANCE REDUCTION VARIABLES

  logical :: survival_biasing = .false.
  real(8) :: weight_cutoff = 0.25_8
  real(8) :: energy_cutoff = ZERO
  real(8) :: weight_survive = ONE

  ! ============================================================================
  ! MISCELLANEOUS VARIABLES

  ! Mode to run in (fixed source, eigenvalue, plotting, etc)
  integer :: run_mode = NONE

  ! Restart run
  logical :: restart_run = .false.
  integer :: restart_batch

  character(MAX_FILE_LEN) :: path_input            ! Path to input file
  character(MAX_FILE_LEN) :: path_cross_sections   ! Path to cross_sections.xml
  character(MAX_FILE_LEN) :: path_multipole        ! Path to wmp library
  character(MAX_FILE_LEN) :: path_source = ''      ! Path to binary source
  character(MAX_FILE_LEN) :: path_state_point      ! Path to binary state point
  character(MAX_FILE_LEN) :: path_source_point     ! Path to binary source point
  character(MAX_FILE_LEN) :: path_particle_restart ! Path to particle restart
  character(MAX_FILE_LEN) :: path_output = ''      ! Path to output directory

  ! The verbosity controls how much information will be printed to the
  ! screen and in logs
  integer :: verbosity = 7

  ! Flag for enabling cell overlap checking during transport
  logical                  :: check_overlaps = .false.
  integer(8), allocatable  :: overlap_check_cnt(:)

  ! Trace for single particle
  logical    :: trace
  integer    :: trace_batch
  integer    :: trace_gen
  integer(8) :: trace_particle

  ! Particle tracks
  logical :: write_all_tracks = .false.
  integer, allocatable :: track_identifiers(:,:)

  ! Particle restart run
  logical :: particle_restart_run = .false.

  ! Number of distribcell maps
  integer :: n_maps

  ! Write out initial source
  logical :: write_initial_source = .false.

  ! Whether create fission neutrons or not. Only applied for MODE_FIXEDSOURCE
  logical :: create_fission_neutrons = .true.

  ! ============================================================================
  ! CMFD VARIABLES

  ! Main object
  type(cmfd_type) :: cmfd

  ! Is CMFD active
  logical :: cmfd_run = .false.

  ! Timing objects
  type(Timer) :: time_cmfd      ! timer for whole cmfd calculation
  type(Timer) :: time_cmfdbuild ! timer for matrix build
  type(Timer) :: time_cmfdsolve ! timer for solver

  ! Flag for active core map
  logical :: cmfd_coremap = .false.

  ! Flag to reset dhats to zero
  logical :: dhat_reset = .false.

  ! Flag to activate neutronic feedback via source weights
  logical :: cmfd_feedback = .false.

  ! User-defined tally information
  integer :: n_cmfd_meshes  = 1 ! # of structured meshes
  integer :: n_cmfd_tallies = 3 ! # of user-defined tallies

  ! Adjoint method type
  character(len=10) :: cmfd_adjoint_type = 'physical'

  ! Number of incomplete ilu factorization levels
  integer :: cmfd_ilu_levels = 1

  ! Batch to begin cmfd
  integer :: cmfd_begin = 1

  ! Tally reset list
  integer :: n_cmfd_resets
  type(SetInt) :: cmfd_reset

  ! Compute effective downscatter cross section
  logical :: cmfd_downscatter = .false.

  ! Convergence monitoring
  logical :: cmfd_power_monitor = .false.

  ! Cmfd output
  logical :: cmfd_write_matrices = .false.

  ! Run an adjoint calculation (last batch only)
  logical :: cmfd_run_adjoint = .false.

  ! CMFD run logicals
  logical :: cmfd_on             = .false.

  ! CMFD display info
  character(len=25) :: cmfd_display = 'balance'

  ! Estimate of spectral radius of CMFD matrices and tolerances
  real(8) :: cmfd_spectral = ZERO
  real(8) :: cmfd_shift = 1.e6
  real(8) :: cmfd_ktol = 1.e-8_8
  real(8) :: cmfd_stol = 1.e-8_8
  real(8) :: cmfd_atoli = 1.e-10_8
  real(8) :: cmfd_rtoli = 1.e-5_8

  ! Information about state points to be written
  integer :: n_state_points = 0
  type(SetInt) :: statepoint_batch

  ! Information about source points to be written
  integer :: n_source_points = 0
  type(SetInt) :: sourcepoint_batch

  ! Various output options
  logical :: output_summary = .true.
  logical :: output_tallies = .true.

  ! ============================================================================
  ! RESONANCE SCATTERING VARIABLES

  logical :: treat_res_scat = .false. ! is resonance scattering treated?
  integer :: n_res_scatterers_total = 0 ! total number of resonant scatterers
  type(Nuclide0K), allocatable, target :: nuclides_0K(:) ! 0K nuclides info

!$omp threadprivate(micro_xs, material_xs, fission_bank, n_bank, &
!$omp&              trace, thread_id, current_work, matching_bins, &
!$omp&              filter_weights)

contains

!===============================================================================
! FREE_MEMORY deallocates and clears  all global allocatable arrays in the
! program
!===============================================================================

  subroutine free_memory()

    integer :: i ! Loop Index

    ! Deallocate cells, surfaces, materials
    if (allocated(cells)) deallocate(cells)
    if (allocated(universes)) deallocate(universes)
    if (allocated(lattices)) deallocate(lattices)
    if (allocated(surfaces)) deallocate(surfaces)
    if (allocated(materials)) deallocate(materials)
    if (allocated(plots)) deallocate(plots)

    ! Deallocate geometry debugging information
    if (allocated(overlap_check_cnt)) deallocate(overlap_check_cnt)

    ! Deallocate cross section data, listings, and cache
    if (allocated(nuclides)) then
    ! First call the clear routines
      do i = 1, size(nuclides)
        call nuclides(i) % clear()
      end do
      deallocate(nuclides)
    end if

    if (allocated(nuclides_0K)) then
      deallocate(nuclides_0K)
    end if

    if (allocated(nuclides_MG)) then
      deallocate(nuclides_MG)
    end if

    if (allocated(macro_xs)) then
      deallocate(macro_xs)
    end if

    if (allocated(sab_tables)) deallocate(sab_tables)
    if (allocated(micro_xs)) deallocate(micro_xs)

    ! Deallocate external source
    if (allocated(external_source)) deallocate(external_source)

    ! Deallocate k and entropy
    if (allocated(k_generation)) deallocate(k_generation)
    if (allocated(entropy)) deallocate(entropy)
    if (allocated(entropy_p)) deallocate(entropy_p)

    ! Deallocate tally-related arrays
    if (allocated(global_tallies)) deallocate(global_tallies)
    if (allocated(meshes)) deallocate(meshes)
    if (allocated(tallies)) deallocate(tallies)
    if (allocated(matching_bins)) deallocate(matching_bins)
    if (allocated(filter_weights)) deallocate(filter_weights)

    ! Deallocate fission and source bank and entropy
!$omp parallel
    if (allocated(fission_bank)) deallocate(fission_bank)
!$omp end parallel
#ifdef _OPENMP
    if (allocated(master_fission_bank)) deallocate(master_fission_bank)
#endif
    if (allocated(source_bank)) deallocate(source_bank)
    if (allocated(entropy_p)) deallocate(entropy_p)

    ! Deallocate array of work indices
    if (allocated(work_index)) deallocate(work_index)

    ! Deallocate CMFD
    call deallocate_cmfd(cmfd)

    ! Deallocate tally node lists
    call active_analog_tallies % clear()
    call active_tracklength_tallies % clear()
    call active_current_tallies % clear()
    call active_collision_tallies % clear()
    call active_tallies % clear()

    ! Deallocate track_identifiers
    if (allocated(track_identifiers)) deallocate(track_identifiers)

    ! Deallocate dictionaries
    call cell_dict % clear()
    call universe_dict % clear()
    call lattice_dict % clear()
    call surface_dict % clear()
    call material_dict % clear()
    call mesh_dict % clear()
    call tally_dict % clear()
    call plot_dict % clear()
    call nuclide_dict % clear()
    call sab_dict % clear()

    ! Clear statepoint and sourcepoint batch set
    call statepoint_batch % clear()
    call sourcepoint_batch % clear()

    ! Deallocate entropy mesh
    if (associated(entropy_mesh)) then
      if (allocated(entropy_mesh % lower_left)) &
           deallocate(entropy_mesh % lower_left)
      if (allocated(entropy_mesh % upper_right)) &
           deallocate(entropy_mesh % upper_right)
      if (allocated(entropy_mesh % width)) deallocate(entropy_mesh % width)
      deallocate(entropy_mesh)
    end if

    ! Deallocate ufs
    if (allocated(source_frac)) deallocate(source_frac)
    if (associated(ufs_mesh)) then
        if (allocated(ufs_mesh % lower_left)) deallocate(ufs_mesh % lower_left)
        if (allocated(ufs_mesh % upper_right)) &
             deallocate(ufs_mesh % upper_right)
        if (allocated(ufs_mesh % width)) deallocate(ufs_mesh % width)
        deallocate(ufs_mesh)
    end if

  end subroutine free_memory

end module global
