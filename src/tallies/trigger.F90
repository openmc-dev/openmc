module trigger

  use, intrinsic :: ISO_C_BINDING

  use trigger_header
  use tally_header, only: tallies

  implicit none

contains

  function n_tally_triggers(i_tally) bind(C) result(n)
    integer(C_INT), value :: i_tally
    integer(C_INT) :: n
    n = tallies(i_tally) % obj % n_triggers
  end function

  function get_tally_trigger(i_tally, i_trig) bind(C) result(trigger)
    integer(C_INT), value :: i_tally
    integer(C_INT), value :: i_trig
    type(C_PTR) :: trigger
    trigger = C_LOC(tallies(i_tally) % obj % triggers(i_trig))
  end function

end module trigger
