module URR_probability_table

  use URR_tally, only: Tally

  implicit none
  private
  public :: ProbabilityTable


!> Type containing data for a single probability table
!! (i.e., one isotope, one temperature, one energy, multiple xs magnitudes)
  type ProbabilityTable

    type(Tally), allocatable :: t(:) ! total xs object
    type(Tally), allocatable :: n(:) ! elastic scattering xs object
    type(Tally), allocatable :: g(:) ! radiative capture xs object
    type(Tally), allocatable :: f(:) ! fission xs object
    type(Tally), allocatable :: x(:) ! competitive xs object
    type(Tally) :: avg_n   ! infinite-dilute elastic xs object
    type(Tally) :: avg_g   ! infinite-dilute capture xs object
    type(Tally) :: avg_f   ! infinite-dilute fission xs object
    type(Tally) :: avg_x   ! infinite-dilute competitive xs object
    type(Tally) :: avg_t   ! infinite-dilute total xs object

  end type ProbabilityTable


end module URR_probability_table
