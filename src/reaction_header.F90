module reaction_header

  use, intrinsic :: ISO_C_BINDING

  use constants,      only: MAX_WORD_LEN
  use hdf5_interface
  use stl_vector,     only: VectorInt
  use string,         only: to_str, starts_with

  implicit none
  private

!===============================================================================
! REACTION contains the cross-section and secondary energy and angle
! distributions for a single reaction in a continuous-energy ACE-format table
!===============================================================================

  type, public :: Reaction
    type(C_PTR) :: ptr
    integer(C_INT)  :: MT                      ! ENDF MT value
    real(C_DOUBLE)  :: Q_value                 ! Reaction Q value
    logical(C_BOOL) :: scatter_in_cm           ! scattering system in center-of-mass?
    logical(C_BOOL) :: redundant               ! redundant reaction?
  contains
    procedure :: from_hdf5
    procedure :: mt_
    procedure :: q_value_
    procedure :: scatter_in_cm_
    procedure :: redundant_
    procedure :: product_decay_rate
    procedure :: product_emission_mode
    procedure :: product_particle
    procedure :: product_sample
    procedure :: product_yield
    procedure :: products_size
    procedure :: sample_elastic_mu
    procedure :: xs
    procedure :: xs_size
    procedure :: xs_threshold
  end type Reaction

  interface
    function reaction_from_hdf5(group, temperatures, n) result(ptr) bind(C)
      import C_PTR, HID_T, C_INT
      integer(HID_T), value :: group
      integer(C_INT), intent(in) :: temperatures
      integer(C_INT), value :: n
      type(C_PTR) :: ptr
    end function

    function reaction_mt(ptr) result(mt) bind(C)
      import C_PTR, C_INT
      type(C_PTR), value :: ptr
      integer(C_INT) :: mt
    end function

    function reaction_q_value(ptr) result(q_value) bind(C)
      import C_PTR, C_DOUBLE
      type(C_PTR), value :: ptr
      real(C_DOUBLE) :: q_value
    end function

    function reaction_scatter_in_cm(ptr) result(b) bind(C)
      import C_PTR, C_BOOL
      type(C_PTR), value :: ptr
      logical(C_BOOL) :: b
    end function

    function reaction_redundant(ptr) result(b) bind(C)
      import C_PTR, C_BOOL
      type(C_PTR), value :: ptr
      logical(C_BOOL) :: b
    end function

    pure function reaction_product_decay_rate(ptr, product) result(rate) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), intent(in), value :: ptr
      integer(C_INT), intent(in), value :: product
      real(C_DOUBLE) :: rate
    end function

    pure function reaction_product_emission_mode(ptr, product) result(m) bind(C)
      import C_PTR, C_INT
      type(C_PTR), intent(in), value :: ptr
      integer(C_INT), intent(in), value :: product
      integer(C_INT) :: m
    end function

    pure function reaction_product_particle(ptr, product) result(particle) bind(C)
      import C_PTR, C_INT
      type(C_PTR), intent(in), value :: ptr
      integer(C_INT), intent(in), value :: product
      integer(C_INT) :: particle
    end function

    subroutine reaction_product_sample(ptr, product, E_in, E_out, mu) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value :: ptr
      integer(C_INT), value :: product
      real(C_DOUBLE), value :: E_in
      real(C_DOUBLE), intent(out) :: E_out
      real(C_DOUBLE), intent(out) :: mu
    end subroutine

    pure function reaction_product_yield(ptr, product, E) result(val) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), intent(in), value :: ptr
      integer(C_INT), intent(in), value :: product
      real(C_DOUBLE), intent(in), value :: E
      real(C_DOUBLE) :: val
    end function

    pure function reaction_products_size(ptr) result(sz) bind(C)
      import C_PTR, C_INT
      type(C_PTR), intent(in), value :: ptr
      integer(C_INT) :: sz
    end function

    function reaction_sample_elastic_mu(ptr, E) result(mu) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value :: ptr
      real(C_DOUBLE), value :: E
      real(C_DOUBLE) :: mu
    end function

    function reaction_xs(ptr, temperature, energy) result(xs) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value :: ptr
      integer(C_INT), value :: temperature
      integer(C_INT), value :: energy
      real(C_DOUBLE) :: xs
    end function

    function reaction_xs_size(ptr, temperature) result(sz) bind(C)
      import C_PTR, C_INT
      type(C_PTR), value :: ptr
      integer(C_INT), value :: temperature
      integer(C_INT) :: sz
    end function

    function reaction_xs_threshold(ptr, temperature) result(threshold) bind(C)
      import C_PTR, C_INT
      type(C_PTR), value :: ptr
      integer(C_INT), value :: temperature
      integer(C_INT) :: threshold
    end function
  end interface

contains

  subroutine from_hdf5(this, group_id, temperatures)
    class(Reaction), intent(inout) :: this
    integer(HID_T),  intent(in)    :: group_id
    type(VectorInt), intent(in)    :: temperatures

    integer(C_INT) :: dummy
    integer(C_INT) :: n

    n = temperatures % size()
    if (n > 0) then
      this % ptr = reaction_from_hdf5(group_id, temperatures % data(1), n)
    else
      ! In this case, temperatures % data(1) doesn't exist, so we just pass a
      ! dummy value
      this % ptr = reaction_from_hdf5(group_id, dummy, n)
    end if
    this % MT = reaction_mt(this % ptr)
    this % Q_value = reaction_q_value(this % ptr)
    this % scatter_in_cm = reaction_scatter_in_cm(this % ptr)
    this % redundant = reaction_redundant(this % ptr)
  end subroutine from_hdf5

  function mt_(this) result(mt)
    class(Reaction), intent(in) :: this
    integer(C_INT) :: MT

    mt = reaction_mt(this % ptr)
  end function

  function q_value_(this) result(q_value)
    class(Reaction), intent(in) :: this
    real(C_DOUBLE) :: q_value

    q_value = reaction_q_value(this % ptr)
  end function

  function scatter_in_cm_(this) result(cm)
    class (Reaction), intent(in) :: this
    logical(C_BOOL) :: cm

    cm = reaction_scatter_in_cm(this % ptr)
  end function

  function redundant_(this) result(redundant)
    class (Reaction), intent(in) :: this
    logical(C_BOOL) :: redundant

    redundant = reaction_redundant(this % ptr)
  end function

  pure function product_decay_rate(this, product) result(rate)
    class(Reaction), intent(in) :: this
    integer(C_INT), intent(in) :: product
    real(C_DOUBLE) :: rate

    rate = reaction_product_decay_rate(this % ptr, product)
  end function

  pure function product_emission_mode(this, product) result(m)
    class(Reaction), intent(in) :: this
    integer(C_INT), intent(in) :: product
    integer(C_INT) :: m

    m = reaction_product_emission_mode(this % ptr, product)
  end function

  pure function product_particle(this, product) result(p)
    class(Reaction), intent(in) :: this
    integer(C_INT), intent(in) :: product
    integer(C_INT) :: p

    p = reaction_product_particle(this % ptr, product)
  end function

  subroutine product_sample(this, product, E_in, E_out, mu)
    class(Reaction), intent(in) :: this
    integer(C_INT), intent(in) :: product
    real(C_DOUBLE), intent(in) :: E_in
    real(C_DOUBLE), intent(out) :: E_out
    real(C_DOUBLE), intent(out) :: mu

    call reaction_product_sample(this % ptr, product, E_in, E_out, mu)
  end subroutine

  pure function product_yield(this, product, E) result(val)
    class(Reaction), intent(in) :: this
    integer(C_INT), intent(in) :: product
    real(C_DOUBLE), intent(in) :: E
    real(C_DOUBLE) :: val

    val = reaction_product_yield(this % ptr, product, E)
  end function

  pure function products_size(this) result(sz)
    class(Reaction), intent(in) :: this
    integer(C_INT) :: sz

    sz = reaction_products_size(this % ptr)
  end function

  function sample_elastic_mu(this, E) result(mu)
    class(Reaction), intent(in) :: this
    real(C_DOUBLE), intent(in) :: E
    real(C_DOUBLE) :: mu

    mu = reaction_sample_elastic_mu(this % ptr, E)
  end function

  function xs(this, temperature, energy) result(val)
    class(Reaction), intent(in) :: this
    integer(C_INT), intent(in) :: temperature
    integer(C_INT), intent(in) :: energy
    real(C_DOUBLE) :: val

    val = reaction_xs(this % ptr, temperature, energy)
  end function

  function xs_size(this, temperature) result(sz)
    class(Reaction), intent(in) :: this
    integer(C_INT) :: temperature
    integer(C_INT) :: sz

    sz = reaction_xs_size(this % ptr, temperature)
  end function

  function xs_threshold(this, temperature) result(val)
    class(Reaction), intent(in) :: this
    integer(C_INT), intent(in) :: temperature
    integer(C_INT) :: val

    val = reaction_xs_threshold(this % ptr, temperature)
  end function

end module reaction_header
