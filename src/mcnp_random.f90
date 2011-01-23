
module mcnp_random
  !=======================================================================
  ! Description:
  !  mcnp_random.F90 -- random number generation routines
  !=======================================================================
  !  This module contains:
  !
  !   * Constants for the RN generator, including initial RN seed for the
  !     problem & the current RN seed
  !
  !   * MCNP interface routines:
  !     - random number function:           rang()
  !     - RN initialization for problem:    RN_init_problem
  !     - RN initialization for particle:   RN_init_particle
  !     - get info on RN parameters:        RN_query
  !     - get RN seed for n-th history:     RN_query_first
  !     - set new RN parameters:            RN_set
  !     - skip-ahead in the RN sequence:    RN_skip_ahead
  !     - Unit tests:        RN_test_basic, RN_test_skip, RN_test_mixed
  !
  !   * For interfacing with the rest of MCNP, arguments to/from these
  !     routines will have types of I8 or I4.
  !     Any args which are to hold random seeds, multipliers,
  !     skip-distance will be type I8, so that 63 bits can be held without
  !     truncation.
  !
  ! Revisions:
  ! * 10-04-2001 - F Brown, initial mcnp version
  ! * 06-06-2002 - F Brown, mods for extended generators
  ! * 12-21-2004 - F Brown, added 3 of LeCuyer's 63-bit mult. RNGs
  ! * 01-29-2005 - J Sweezy, Modify to use mcnp modules prior to automatic
  !                io unit numbers.
  ! * 12-02-2005 - F Brown, mods for consistency with C version
  !=======================================================================

  !-------------------
  ! MCNP output units
  !-------------------
  integer, parameter ::   iuo = 6
  integer, parameter ::  jtty = 6

  PRIVATE
  !---------------------------------------------------
  ! Kinds for LONG INTEGERS (64-bit) & REAL*8 (64-bit)
  !---------------------------------------------------
  integer, parameter :: R8 = selected_real_kind(15,307)
  integer, parameter :: I8 = selected_int_kind(18)

  !-----------------------------------
  ! Public functions and subroutines for this module
  !-----------------------------------
  PUBLIC ::  rang
  PUBLIC ::  RN_init_problem
  PUBLIC ::  RN_init_particle
  PUBLIC ::  RN_set
  PUBLIC ::  RN_query
  PUBLIC ::  RN_query_first
  PUBLIC ::  RN_update_stats
  PUBLIC ::  RN_test_basic
  PUBLIC ::  RN_test_skip
  PUBLIC ::  RN_test_mixed

  !-------------------------------------
  ! Constants for standard RN generators
  !-------------------------------------
  type :: RN_GEN
    integer          :: index
    integer(I8)      :: mult        ! generator (multiplier)
    integer(I8)      :: add         ! additive constant
    integer          :: log2mod     ! log2 of modulus, must be <64
    integer(I8)      :: stride      ! stride for particle skip-ahead
    integer(I8)      :: initseed    ! default seed for problem
    character(len=8) :: name
  end type RN_GEN

  ! parameters for standard generators
  integer,      parameter :: n_RN_GEN = 7
  type(RN_GEN), SAVE      :: standard_generator(n_RN_GEN)
  data standard_generator / &
    & RN_GEN( 1,      19073486328125_I8, 0_I8, 48, 152917_I8, 19073486328125_I8 , 'mcnp std' ), &
    & RN_GEN( 2, 9219741426499971445_I8, 1_I8, 63, 152917_I8, 1_I8,               'LEcuyer1' ), &
    & RN_GEN( 3, 2806196910506780709_I8, 1_I8, 63, 152917_I8, 1_I8,               'LEcuyer2' ), &
    & RN_GEN( 4, 3249286849523012805_I8, 1_I8, 63, 152917_I8, 1_I8,               'LEcuyer3' ), &
    & RN_GEN( 5, 3512401965023503517_I8, 0_I8, 63, 152917_I8, 1_I8,               'LEcuyer4' ), &
    & RN_GEN( 6, 2444805353187672469_I8, 0_I8, 63, 152917_I8, 1_I8,               'LEcuyer5' ), &
    & RN_GEN( 7, 1987591058829310733_I8, 0_I8, 63, 152917_I8, 1_I8,               'LEcuyer6' )  &
    & /

  !-----------------------------------------------------------------
  !   * Linear multiplicative congruential RN algorithm:
  !
  !            RN_SEED = RN_SEED*RN_MULT + RN_ADD  mod RN_MOD
  !
  !   * Default values listed below will be used, unless overridden
  !-----------------------------------------------------------------
  integer,     SAVE :: RN_INDEX  = 1
  integer(I8), SAVE :: RN_MULT   =  19073486328125_I8
  integer(I8), SAVE :: RN_ADD    =               0_I8
  integer,     SAVE :: RN_BITS   = 48
  integer(I8), SAVE :: RN_STRIDE =          152917_I8
  integer(I8), SAVE :: RN_SEED0  =  19073486328125_I8
  integer(I8), SAVE :: RN_MOD    = 281474976710656_I8
  integer(I8), SAVE :: RN_MASK   = 281474976710655_I8
  integer(I8), SAVE :: RN_PERIOD =  70368744177664_I8
  real(R8),    SAVE :: RN_NORM   = 1._R8 / 281474976710656._R8

  !------------------------------------
  ! Private data for a single particle
  !------------------------------------
  integer(I8) :: RN_SEED   = 19073486328125_I8 ! current seed
  integer(I8) :: RN_COUNT  = 0_I8              ! current counter
  integer(I8) :: RN_NPS    = 0_I8              ! current particle number

  common                /RN_THREAD/   RN_SEED, RN_COUNT, RN_NPS
  save                  /RN_THREAD/
  !$OMP THREADprivate ( /RN_THREAD/ )

  !------------------------------------------
  ! Shared data, to collect info on RN usage
  !------------------------------------------
  integer(I8), SAVE :: RN_COUNT_TOTAL   = 0  ! total RN count all particles
  integer(I8), SAVE :: RN_COUNT_STRIDE  = 0  ! count for stride exceeded
  integer(I8), SAVE :: RN_COUNT_MAX     = 0  ! max RN count all particles
  integer(I8), SAVE :: RN_COUNT_MAX_NPS = 0  ! part index for max count

  !---------------------------------------------------------------------
  ! Reference data:  Seeds for case of init.seed = 1,
  !                  Seed numbers for index 1-5, 123456-123460
  !---------------------------------------------------------------------
  integer(I8), dimension(10,n_RN_GEN) ::  RN_CHECK
  data  RN_CHECK /  &
    ! ***** 1 ***** mcnp standard gen *****
    &      19073486328125_I8,      29763723208841_I8,     187205367447973_I8, &
    &     131230026111313_I8,     264374031214925_I8,     260251000190209_I8, &
    &     106001385730621_I8,     232883458246025_I8,      97934850615973_I8, &
    &     163056893025873_I8, &
    ! ***** 2 *****
    & 9219741426499971446_I8,  666764808255707375_I8, 4935109208453540924_I8, &
    & 7076815037777023853_I8, 5594070487082964434_I8, 7069484152921594561_I8, &
    & 8424485724631982902_I8,   19322398608391599_I8, 8639759691969673212_I8, &
    & 8181315819375227437_I8, &
    ! ***** 3 *****
    & 2806196910506780710_I8, 6924308458965941631_I8, 7093833571386932060_I8, &
    & 4133560638274335821_I8,  678653069250352930_I8, 6431942287813238977_I8, &
    & 4489310252323546086_I8, 2001863356968247359_I8,  966581798125502748_I8, &
    & 1984113134431471885_I8, &
    ! ***** 4 *****
    & 3249286849523012806_I8, 4366192626284999775_I8, 4334967208229239068_I8, &
    & 6386614828577350285_I8, 6651454004113087106_I8, 2732760390316414145_I8, &
    & 2067727651689204870_I8, 2707840203503213343_I8, 6009142246302485212_I8, &
    & 6678916955629521741_I8, &
    ! ***** 5 *****
    & 3512401965023503517_I8, 5461769869401032777_I8, 1468184805722937541_I8, &
    & 5160872062372652241_I8, 6637647758174943277_I8,  794206257475890433_I8, &
    & 4662153896835267997_I8, 6075201270501039433_I8,  889694366662031813_I8, &
    & 7299299962545529297_I8, &
    ! ***** 6 *****
    & 2444805353187672469_I8,  316616515307798713_I8, 4805819485453690029_I8, &
    & 7073529708596135345_I8, 3727902566206144773_I8, 1142015043749161729_I8, &
    & 8632479219692570773_I8, 2795453530630165433_I8, 5678973088636679085_I8, &
    & 3491041423396061361_I8, &
    ! ***** 7 *****
    & 1987591058829310733_I8, 5032889449041854121_I8, 4423612208294109589_I8, &
    & 3020985922691845009_I8, 5159892747138367837_I8, 8387642107983542529_I8, &
    & 8488178996095934477_I8,  708540881389133737_I8, 3643160883363532437_I8, &
    & 4752976516470772881_I8  /
  !---------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------

  function rang()
    ! MCNP random number generator
    !
    ! ***************************************
    ! ***** modifies RN_SEED & RN_COUNT *****
    ! ***************************************
    implicit none
    real(R8) ::  rang

    RN_SEED  = iand( iand( RN_MULT*RN_SEED, RN_MASK) + RN_ADD,  RN_MASK)
    rang     = RN_SEED * RN_NORM
    RN_COUNT = RN_COUNT + 1

    return
  end function rang

  !-------------------------------------------------------------------

  function RN_skip_ahead( seed, skip )
    ! advance the seed "skip" RNs:   seed*RN_MULT^n mod RN_MOD
    implicit none
    integer(I8) :: RN_skip_ahead
    integer(I8), intent(in)  :: seed, skip
    integer(I8) :: nskip, gen, g, inc, c, gp, rn, seed_old

    seed_old = seed
    ! add period till nskip>0
    nskip = skip
    do while( nskip<0_I8 )
      if( RN_PERIOD>0_I8 ) then
        nskip = nskip + RN_PERIOD
      else
        nskip = nskip + RN_MASK
        nskip = nskip + 1_I8
      endif
    enddo

    ! get gen=RN_MULT^n,  in log2(n) ops, not n ops !
    nskip = iand( nskip, RN_MASK )
    gen   = 1
    g     = RN_MULT
    inc   = 0
    c     = RN_ADD
    do while( nskip>0_I8 )
      if( btest(nskip,0) )  then
        gen = iand( gen*g, RN_MASK )
        inc = iand( inc*g, RN_MASK )
        inc = iand( inc+c, RN_MASK )
      endif
      gp    = iand( g+1,  RN_MASK )
      g     = iand( g*g,  RN_MASK )
      c     = iand( gp*c, RN_MASK )
      nskip = ishft( nskip, -1 )
    enddo
    rn = iand( gen*seed_old, RN_MASK )
    rn = iand( rn + inc, RN_MASK )
    RN_skip_ahead = rn
    return
  end function RN_skip_ahead

  !-------------------------------------------------------------------

  subroutine RN_init_problem( new_standard_gen, new_seed, &
    &                         new_stride, new_part1,  print_info )
    ! * initialize MCNP random number parameters for problem,
    !   based on user input.  This routine should be called
    !   only from the main thread, if OMP threading is being used.
    !
    ! * for initial & continue runs, these args should be set:
    !     new_standard_gen - index of built-in standard RN generator,
    !                        from RAND gen=   (or dbcn(14)
    !     new_seed   - from RAND seed=        (or dbcn(1))
    !     output     - logical, print RN seed & mult if true
    !
    !     new_stride - from RAND stride=      (or dbcn(13))
    !     new_part1  - from RAND hist=        (or dbcn(8))
    !
    ! * for continue runs only, these should also be set:
    !     new_count_total   - from "rnr"   at end of previous run
    !     new_count_stride  - from nrnh(1) at end of previous run
    !     new_count_max     - from nrnh(2) at end of previous run
    !     new_count_max_nps - from nrnh(3) at end of previous run
    !
    ! * check on size of long-ints & long-int arithmetic
    ! * check the multiplier
    ! * advance the base seed for the problem
    ! * set the initial particle seed
    ! * initialize the counters for RN stats
    implicit none
    integer,     intent(in) :: new_standard_gen
    integer(I8), intent(in) :: new_seed
    integer(I8), intent(in) :: new_stride
    integer(I8), intent(in) :: new_part1
    integer,     intent(in) :: print_info
    character(len=20) :: printseed
    integer(I8)       ::  itemp1, itemp2, itemp3, itemp4

    if( new_standard_gen<1 .or. new_standard_gen>n_RN_GEN ) then
      call expire( 0, 'RN_init_problem', &
        & ' ***** ERROR: illegal index for built-in RN generator')
    endif
      
    ! set defaults, override if input supplied: seed, mult, stride
    RN_INDEX   = new_standard_gen
    RN_MULT    = standard_generator(RN_INDEX)%mult
    RN_ADD     = standard_generator(RN_INDEX)%add
    RN_STRIDE  = standard_generator(RN_INDEX)%stride
    RN_SEED0   = standard_generator(RN_INDEX)%initseed
    RN_BITS    = standard_generator(RN_INDEX)%log2mod
    RN_MOD     = ishft( 1_I8,       RN_BITS )
    RN_MASK    = ishft( not(0_I8),  RN_BITS-64 )
    RN_NORM    = 2._R8**(-RN_BITS)
    if( RN_ADD==0_I8) then
      RN_PERIOD  = ishft( 1_I8, RN_BITS-2 )
    else
      RN_PERIOD  = ishft( 1_I8, RN_BITS )
    endif
    if( new_seed>0_I8 ) then
      RN_SEED0  = new_seed
    endif
    if( new_stride>0_I8 ) then
      RN_STRIDE = new_stride
    endif
    RN_COUNT_TOTAL   = 0
    RN_COUNT_STRIDE  = 0
    RN_COUNT_MAX     = 0
    RN_COUNT_MAX_NPS = 0

    if( print_info /= 0 ) then
      write(printseed,'(i20)') RN_SEED0
      write( iuo,1) RN_INDEX, RN_SEED0, RN_MULT, RN_ADD, RN_BITS, RN_STRIDE
      write(jtty,2) RN_INDEX, adjustl(printseed)
1     format( &
        & /,' ***************************************************', &
        & /,' * Random Number Generator  = ',i20,             ' *', &
        & /,' * Random Number Seed       = ',i20,             ' *', &
        & /,' * Random Number Multiplier = ',i20,             ' *', &
        & /,' * Random Number Adder      = ',i20,             ' *', &
        & /,' * Random Number Bits Used  = ',i20,             ' *', &
        & /,' * Random Number Stride     = ',i20,             ' *', &
        & /,' ***************************************************',/)
2     format(' comment. using random number generator ',i2, &
        &    ', initial seed = ',a20)
    endif

    ! double-check on number of bits in a long int
    if( bit_size(RN_SEED)<64 ) then
      call expire( 0, 'RN_init_problem', &
        & ' ***** ERROR: <64 bits in long-int, can-t generate RN-s')
    endif
    itemp1 = 5_I8**25
    itemp2 = 5_I8**19
    itemp3 = ishft(2_I8**62-1_I8,1) + 1_I8
    itemp4 = itemp1*itemp2
    if( iand(itemp4,itemp3)/=8443747864978395601_I8 ) then
      call expire( 0, 'RN_init_problem', &
        & ' ***** ERROR: can-t do 64-bit integer ops for RN-s')
    endif

    if( new_part1>1_I8 ) then
      ! advance the problem seed to that for part1
      RN_SEED0 = RN_skip_ahead( RN_SEED0, (new_part1-1_I8)*RN_STRIDE )
      itemp1   = RN_skip_ahead( RN_SEED0, RN_STRIDE )
      if( print_info /= 0 ) then
        write(printseed,'(i20)') itemp1
        write( iuo,3) new_part1,  RN_SEED0, itemp1
        write(jtty,4) new_part1,  adjustl(printseed)
3       format( &
          & /,' ***************************************************', &
          & /,' * Random Number Seed will be advanced to that for *', &
          & /,' * previous particle number = ',i20,             ' *', &
          & /,' * New RN Seed for problem  = ',i20,             ' *', &
          & /,' * Next Random Number Seed  = ',i20,             ' *', &
          & /,' ***************************************************',/)
4       format(' comment. advancing random number to particle ',i12, &
          &    ', initial seed = ',a20)
      endif
    endif

    ! set the initial particle seed
    RN_SEED  = RN_SEED0
    RN_COUNT = 0
    RN_NPS   = 0

    return
  end subroutine RN_init_problem

  !-------------------------------------------------------------------

  subroutine RN_init_particle( nps )
    ! initialize MCNP random number parameters for particle "nps"
    !
    !     * generate a new particle seed from the base seed
    !       & particle index
    !     * set the RN count to zero
    implicit none
    integer(I8), intent(in) :: nps

    RN_SEED  = RN_skip_ahead( RN_SEED0, nps*RN_STRIDE )
    RN_COUNT = 0
    RN_NPS   = nps

    return
  end subroutine RN_init_particle

  !-------------------------------------------------------------------

  subroutine RN_set(  key,  value )
    implicit none
    character(len=*), intent(in) :: key
    integer(I8),      intent(in) :: value
    character(len=20) :: printseed
    integer(I8) :: itemp1

    if( key == "stride"        ) then
      if( value>0_I8 ) then
        RN_STRIDE        = value
      endif
    endif
    if( key == "count_total"   )  RN_COUNT_TOTAL   = value
    if( key == "count_stride"  )  RN_COUNT_STRIDE  = value
    if( key == "count_max"     )  RN_COUNT_MAX     = value
    if( key == "count_max_nps" )  RN_COUNT_MAX_NPS = value
    if( key == "seed"          )  then
      if( value>0_I8 ) then
        RN_SEED0 = value
        RN_SEED  = RN_SEED0
        RN_COUNT = 0
        RN_NPS   = 0
      endif
    endif
    if( key == "part1" ) then
      if( value>1_I8 ) then
        ! advance the problem seed to that for part1
        RN_SEED0 = RN_skip_ahead( RN_SEED0, (value-1_I8)*RN_STRIDE )
        itemp1   = RN_skip_ahead( RN_SEED0, RN_STRIDE )
        write(printseed,'(i20)') itemp1
        write( iuo,3) value,  RN_SEED0, itemp1
        write(jtty,4) value,  adjustl(printseed)
3       format( &
          & /,' ***************************************************', &
          & /,' * Random Number Seed will be advanced to that for *', &
          & /,' * previous particle number = ',i20,             ' *', &
          & /,' * New RN Seed for problem  = ',i20,             ' *', &
          & /,' * Next Random Number Seed  = ',i20,             ' *', &
          & /,' ***************************************************',/)
4       format(' comment. advancing random number to particle ',i12, &
          &    ', initial seed = ',a20)
        RN_SEED  = RN_SEED0
        RN_COUNT = 0
        RN_NPS   = 0
      endif
    endif
    return
  end subroutine RN_set

  !-------------------------------------------------------------------

  function RN_query( key )
    implicit none
    integer(I8)                  :: RN_query
    character(len=*), intent(in) :: key
    RN_query = 0_I8
    if( key == "seed"          )  RN_query = RN_SEED
    if( key == "stride"        )  RN_query = RN_STRIDE
    if( key == "mult"          )  RN_query = RN_MULT
    if( key == "add"           )  RN_query = RN_ADD
    if( key == "count"         )  RN_query = RN_COUNT
    if( key == "period"        )  RN_query = RN_PERIOD
    if( key == "count_total"   )  RN_query = RN_COUNT_TOTAL
    if( key == "count_stride"  )  RN_query = RN_COUNT_STRIDE
    if( key == "count_max"     )  RN_query = RN_COUNT_MAX
    if( key == "count_max_nps" )  RN_query = RN_COUNT_MAX_NPS
    if( key == "first"         )  RN_query = RN_SEED0
    return
  end function RN_query
  !-------------------------------------------------------------------

  function RN_query_first( nps )
    implicit none
    integer(I8)                  :: RN_query_first
    integer(I8),      intent(in) :: nps
    RN_query_first = RN_skip_ahead( RN_SEED0, nps*RN_STRIDE )
    return
  end function RN_query_first

  !-------------------------------------------------------------------

  subroutine RN_update_stats()
    ! update overall RN count info
    implicit none

    !$OMP CRITICAL (RN_STATS)

    RN_COUNT_TOTAL = RN_COUNT_TOTAL + RN_COUNT

    if( RN_COUNT>RN_COUNT_MAX ) then
      RN_COUNT_MAX     = RN_COUNT
      RN_COUNT_MAX_NPS = RN_NPS
    endif

    if( RN_COUNT>RN_STRIDE ) then
      RN_COUNT_STRIDE = RN_COUNT_STRIDE + 1
    endif

    !$OMP END CRITICAL (RN_STATS)

    RN_COUNT = 0
    RN_NPS   = 0

    return
  end subroutine RN_update_stats

  !-------------------------------------------------------------------

  subroutine expire( i, c1, c2 )
    integer,          intent(in) :: i
    character(len=*), intent(in) :: c1, c2
    write(*,*) ' ********** error: ',c1
    write(*,*) ' ********** error: ',c2
    stop '**error**'
  end subroutine expire

  !-------------------------------------------------------------------
  !###################################################################
  !#
  !#  Unit tests
  !#
  !###################################################################

  subroutine RN_test_basic( new_gen )
    ! test routine for basic random number generator
    implicit none
    integer, intent(in) :: new_gen
    real(R8)    :: s
    integer(I8) :: seeds(10)
    integer     :: i, j

    write(jtty,"(/,a)")  " ***** random number - basic test *****"

    ! set the seed
    call RN_init_problem( new_gen, 1_I8, 0_I8, 0_I8, 1 )

    ! get the first 5 seeds, then skip a few, get 5 more - directly
    s = 0.0_R8
    do  i = 1,5
      s = s + rang()
      seeds(i) = RN_query( "seed" )
    enddo
    do  i = 6,123455
      s = s + rang()
    enddo
    do  i = 6,10
      s = s + rang()
      seeds(i) = RN_query( "seed" )
    enddo

    ! compare
    do  i = 1,10
      j = i
      if( i>5  ) j = i + 123450
      write(jtty,"(1x,i6,a,i20,a,i20)") &
        &  j, "  reference: ", RN_CHECK(i,new_gen), "  computed: ", seeds(i)
      if( seeds(i)/=RN_CHECK(i,new_gen) ) then
        write(jtty,"(a)")  " ***** basic_test of RN generator failed:"
      endif
    enddo
    return
  end subroutine RN_test_basic

  !-------------------------------------------------------------------

  subroutine RN_test_skip( new_gen )
    ! test routine for basic random number generation & skip-ahead
    implicit none
    integer, intent(in) :: new_gen
    integer(I8) :: seeds(10)
    integer     :: i, j

    ! set the seed
    call RN_init_problem( new_gen, 1_I8, 0_I8, 0_I8, 0 )

    ! use the skip-ahead function to get first 5 seeds, then 5 more
    do i = 1,10
      j = i
      if( i>5 )  j = i + 123450
      seeds(i) = RN_skip_ahead( 1_I8, int(j,I8) )
    enddo

    ! compare
    write(jtty,"(/,a)")  " ***** random number - skip test *****"
    do i = 1,10
      j = i
      if( i>5  ) j = i + 123450
      write(jtty,"(1x,i6,a,i20,a,i20)") &
        &  j, "  reference: ", RN_CHECK(i,new_gen),  "  computed: ", seeds(i)
      if( seeds(i)/=RN_CHECK(i,new_gen) ) then
        write(jtty,"(a)")  " ***** skip_test of RN generator failed:"
      endif
    enddo
    return
  end subroutine RN_test_skip

  !-------------------------------------------------------------------

  subroutine RN_test_mixed( new_gen )
    ! test routine -- print RN's 1-5 & 123456-123460,
    !                 with reference vals
    implicit none
    integer, intent(in) :: new_gen
    integer(I8) :: r
    integer     :: i, j

    write(jtty,"(/,a)")  " ***** random number - mixed test *****"
    ! set the seed & set the stride to 1
    call RN_init_problem( new_gen, 1_I8, 1_I8, 0_I8, 0 )

    write(jtty,"(a,i20,z20)") " RN_MULT   = ", RN_MULT, RN_MULT
    write(jtty,"(a,i20,z20)") " RN_ADD    = ", RN_ADD,  RN_ADD
    write(jtty,"(a,i20,z20)") " RN_MOD    = ", RN_MOD,  RN_MOD
    write(jtty,"(a,i20,z20)") " RN_MASK   = ", RN_MASK, RN_MASK
    write(jtty,"(a,i20)")     " RN_BITS   = ", RN_BITS
    write(jtty,"(a,i20)")     " RN_PERIOD = ", RN_PERIOD
    write(jtty,"(a,es20.14)") " RN_NORM   = ", RN_NORM
    write(jtty,"(a)")  " "
    do i = 1,10
      j = i
      if( i>5  ) j = i + 123450
      call RN_init_particle( int(j,I8) )
      r = RN_query( "seed" )
      write(jtty,"(1x,i6,a,i20,a,i20)") &
        &  j, "  reference: ", RN_CHECK(i,new_gen),"  computed: ", r
      if( r/=RN_CHECK(i,new_gen) ) then
        write(jtty,"(a)")  " ***** mixed test of RN generator failed:"
      endif
    enddo
    return
  end subroutine RN_test_mixed

  !-------------------------------------------------------------------
end module mcnp_random
