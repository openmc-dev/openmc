module endf_reader

  use ace_header, only: Nuclide
  use error,      only: fatal_error
  use global
  use output,     only: write_message

  implicit none

  integer :: in = 11     ! input unit

contains

!===============================================================================
! READ_ENDF reads an ENDF file
!===============================================================================

  subroutine read_endf(filename, i_nuc)

    character(80) :: filename    ! ENDF filename
    integer       :: i_nuc       ! index in global nuclides array
    logical       :: file_exists ! does ENDF file exist?
    character(7)  :: readable    ! is ENDF file readable?
    character(80) :: rec         ! ENDF file record
    integer       :: MAT         ! MAT material #
    integer       :: MF          ! MF file number
    integer       :: MT          ! MT type number
    integer       :: NS          ! line #

    inquire(file = trim(path_endf)//trim(filename), &
      & exist = file_exists, read = readable)
    if (.not. file_exists) then
      message = 'ENDF file '//trim(filename)//' does not exist.'
      call fatal_error()
    else if (readable(1:3) == 'NO') then
      message = 'ENDF file '//trim(filename)// &
        & ' is not readable.  Change file permissions with chmod command.'
      call fatal_error()
    end if

    ! display message
    message = "Loading ENDF file: "//trim(filename)
    call write_message(6)

    open(unit = in, &
      file = trim(path_endf)//trim(filename))

    do
      read(in, 10) rec
10    format(A80)
      read(rec(67:70), '(I4)') MAT
      read(rec(71:72), '(I2)') MF
      read(rec(73:75), '(I3)') MT
      read(rec(76:80), '(I5)') NS
      if (MF == 2 .and. MT == 151) then
        call read_MF2(i_nuc)
      else if (MAT == -1) then
        exit
      end if
    end do

    close(in)

  end subroutine read_endf

  subroutine read_MF2(i_n)

    integer       :: i_n    ! index in global nuclides array
    character(80) :: rec    ! ENDF file record
    character(10) :: blank  ! 10 blank spaces
    integer       :: MAT    ! MAT material #
    integer       :: MF     ! MF file number
    integer       :: MT     ! MT type number
    integer       :: NS
    integer       :: i_L
    integer       :: i_J
    integer       :: i_E
    integer       :: NE_tmp
    type(Nuclide), pointer :: nuc => null()

    nuc => nuclides(i_n)

! TODO: LRF=2 is only supported LRF
    do
      read(in, 10) rec
10    format(A80)

      read(rec(67:70), '(I4)')  MAT
      read(rec(71:72), '(I2)')  MF
      read(rec(73:75), '(I3)')  MT
      read(rec(76:80), '(I5)')  NS
      read(rec(23:32), '(A10)') blank

      if (blank == '          ') then

        read(rec(24:33), '(I10)') nuc % LRU
        read(rec(35:44), '(I10)') nuc % LRF

        if (nuc % LRU == 2) then
          read(rec(1:11),  '(E11.0)') nuc % EL
          read(rec(12:22), '(E11.0)') nuc % EH
          read(rec(45:55),   '(I11)') nuc % NRO
          read(rec(56:66),   '(I11)') nuc % NAPS

          read(in, 10) rec
          read(rec(1:11),  '(E11.0)') nuc % SPI
          read(rec(12:22), '(E11.0)') nuc % AP
          call nuc % pprocess()
          read(rec(23:33),   '(I11)') nuc % LSSF
          read(rec(45:55),   '(I11)') nuc % NLS

          ! allocate total angular momentum quantum #'s
          allocate(nuc % NJS(nuc % NLS))
          allocate(nuc % J_grid(nuc % NLS))

          ! allocate degress of freedom for partial widths
          allocate(nuc % AMUX_grid(nuc % NLS))
          allocate(nuc % AMUN_grid(nuc % NLS))
          allocate(nuc % AMUG_grid(nuc % NLS))
          allocate(nuc % AMUF_grid(nuc % NLS))

          do i_L = 1, nuc % NLS

            read(in, 10) rec
            read(rec(45:55),   '(I11)') nuc % NJS(i_L)
          
            allocate(nuc % J_grid(i_L)    % vals(nuc % NJS(i_L)))
            allocate(nuc % AMUX_grid(i_L) % vals(nuc % NJS(i_L)))
            allocate(nuc % AMUN_grid(i_L) % vals(nuc % NJS(i_L)))
            allocate(nuc % AMUG_grid(i_L) % vals(nuc % NJS(i_L)))
            allocate(nuc % AMUF_grid(i_L) % vals(nuc % NJS(i_L)))

            do i_J = 1, nuc % NJS(i_L)

              read(in, 10) rec
              read(rec( 1:11), '(E11.0)') nuc % J_grid(i_L) % vals(i_J)
              read(rec(23:33), '(I11)') nuc % INT
              read(rec(56:66), '(I11)') NE_tmp

              nuc % NE = NE_tmp

              if (.not. allocated(nuc % ES)) then
                ! allocate mean widths and spacings
                allocate(nuc % ES(nuc % NE))
                allocate(nuc % D_means(nuc % NE, nuc % NLS))
                allocate(nuc % Gam_n_means(nuc % NE, nuc % NLS))
                allocate(nuc % Gam_gam_means(nuc % NE, nuc % NLS))
                allocate(nuc % Gam_f_means(nuc % NE, nuc % NLS))
                allocate(nuc % Gam_x_means(nuc % NE, nuc % NLS))
              end if

              if (.not. allocated(nuc % D_means(1, i_L) % vals)) then
                do i_E = 1, nuc % NE
                  allocate(nuc % D_means(i_E, i_L)       % vals(nuc % NJS(i_L)))
                  allocate(nuc % Gam_n_means(i_E, i_L)   % vals(nuc % NJS(i_L)))
                  allocate(nuc % Gam_gam_means(i_E, i_L) % vals(nuc % NJS(i_L)))
                  allocate(nuc % Gam_f_means(i_E, i_L)   % vals(nuc % NJS(i_L)))
                  allocate(nuc % Gam_x_means(i_E, i_L)   % vals(nuc % NJS(i_L)))
                end do
              end if

              read(in, 10) rec
              read(rec(23:33), '(E11.0)') nuc % AMUX_grid(i_L) % vals(i_J)
              read(rec(34:44), '(E11.0)') nuc % AMUN_grid(i_L) % vals(i_J)
              read(rec(45:55), '(E11.0)') nuc % AMUG_grid(i_L) % vals(i_J)
              read(rec(56:66), '(E11.0)') nuc % AMUF_grid(i_L) % vals(i_J)

              do i_E = 1, nuc % NE
                read(in, 10) rec
                read(rec( 1:11), '(E11.0)') nuc % ES(i_E)
                read(rec(12:22), '(E11.0)') nuc % D_means(i_E, i_L) &
                                              & % vals(i_J)
                read(rec(23:33), '(E11.0)') nuc % Gam_x_means(i_E, i_L) &
                                              & % vals(i_J)
                read(rec(34:44), '(E11.0)') nuc % Gam_n_means(i_E, i_L) &
                                              & % vals(i_J)
                read(rec(45:55), '(E11.0)') nuc % Gam_gam_means(i_E, i_L) &
                                              & % vals(i_J)
                read(rec(56:66), '(E11.0)') nuc % Gam_f_means(i_E, i_L) &
                                              & % vals(i_J)
              end do
            end do
          end do
        end if
      end if
      if (MT == 0) then
        read(in, 10) rec
        exit
      end if
    end do

  end subroutine read_MF2

  subroutine print_shit(i_nuc)

    integer :: i_nuc
    integer :: iL
    integer :: iJ
    integer :: iE
    type(Nuclide), pointer :: nuc => null()

    nuc => nuclides(i_nuc)

    write(*,'(I6,ES10.3,L1,I5)') nuc%zaid,nuc%kT,nuc%otf_urr,nuc%n_resonances
    write(*,'(ES10.3,ES10.3)') nuc%EL,nuc%EH
    write(*,'(I2,I2)') nuc%NRO,nuc%NAPS
    write(*,'(ES10.3,ES10.3,ES10.3)') nuc%SPI,nuc%AP,nuc%ac
    write(*,'(I2)') nuc%LSSF
    write(*,'(I2)') nuc%NLS

    do iL = 1, nuc%NLS
      do iJ = 1, nuc%NJS(iL)
        write(*,'(I2,ES10.3)') iL,nuc%J_grid(iL)%vals(iJ)
        do iE = 1, nuc%NE
          write(*,'(ES10.3,ES10.3,ES10.3,ES10.3,ES10.3,ES10.3)') nuc%ES(iE),nuc%D_means(iE,iL)%vals(iJ),&
            & nuc%Gam_x_means(iE,iL)%vals(iJ), &
            & nuc%Gam_n_means(iE,iL)%vals(iJ), &
            & nuc%Gam_gam_means(iE,iL)%vals(iJ), &
            & nuc%Gam_f_means(iE,iL)%vals(iJ)
        end do
      end do
    end do
  end subroutine print_shit

end module endf_reader
