program main

  use global, only: cells, surfaces, materials, inputfile
  use fileio, only: read_input, read_command_line
  use output, only: title, message, warning, error
  use geometry, only: sense, cell_contains

  implicit none

  character(16) :: filename
  character(250) :: msg

  real(8) :: point(3)
  integer :: s(4)

  call title()
  
  call read_command_line()
  call read_input(inputfile)

  point = (/ 4.0, 0.0, 3.1 /)

  call cell_contains( cells(1), point, .true. )


end program main
