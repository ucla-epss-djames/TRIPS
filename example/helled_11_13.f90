program helled_11_13

  use iso_fortran_env

  use mod_consts
  use mod_math
  use mod_structure
  use mod_tidal, only : propagator_method, normalize, use_GD_core

  implicit none

  ! parameters for benchmark
  real(real64), parameter :: l = 2
  real(real64), parameter :: OMEGA_U = 1./(17.24*3600.)
  integer(int64), parameter :: index = 5
  integer(int64), parameter :: NMAX = 10002
  character(*), parameter, dimension(index) :: files = ["L34   ", "R10001", "R1001 ", "R34   ", "core  "]
  character(*), parameter :: dpath = "../data/uranus_helled_"
  character(*), parameter :: rpath = "../post/processed/helled_11_13/"

  ! locals
  complex(real64), allocatable :: sd(:,:), tidal(:,:)
  real(real64), dimension(NMAX) :: read_r, read_rho, read_mu, read_eta
  real(real64), allocatable  :: data(:,:)
  real(real64) :: vals(5), track, zero, mass
  integer(int64) :: i, j, k
  logical :: flag
  character(len=:), allocatable :: path
  character(8)  :: date
  character(10)  :: time

  zero = 0

  do i = 1, index

     path = trim(dpath//files(i))
     open(10, file=path, status="old")
     track = 1
     j = 1
     mass = 0

     do while(track .ne. zero)

        read(10,*) vals(1), vals(2), vals(3), vals(4), vals(5)

        read_r(j) = vals(1)
        read_rho(j) = vals(2)
        read_mu(j) = vals(4)
        read_eta(j) = vals(5)

        if(read_eta(j) == 1.) then
           read_eta(j) = 0
        end if

        track = read_r(j)
        j = j + 1

     end do

     j = j - 1
     close(10)

     allocate(sd(j,4), tidal(j,3), data(j,2))

     ! struct files are backwards of how TRIPS reads in data
     sd(:,1) = read_r(j:1:-1)
     sd(1,1) = 0.001
     sd(:,4) = read_rho(j:1:-1)

     data(:,1) = read_mu(j:1:-1)

     ! UNCOMMENT ONE OF THE LINES FOR THE CASE YOU WANT
     ! VISCOELASTIC
     data(:,2) = read_eta(j:1:-1)

     ! FLUID
     ! data(:,2) = 0

     ! check which core to use
     call use_GD_core(sd(1,1)%re, data(1,1), flag)

     call planet_structure(j, OMEGA_U, data, mass, sd)

     call normalize(mass, sd(j,1)%re, sd)

     call date_and_time(DATE=date, TIME=time)
     open(11, file=rpath//date//"-"//time//"-uranus_struct_"//trim(files(i)), status="new")
     write(11, '(5(e16.5))') (sd(k,1)%re, sd(k,2)%re, sd(k,2)%im, sd(k,3)%re, sd(k,4)%re, k = 1, j)
     close(11)

     tidal = propagator_method(l, j, sd, flag)
     open(11, file=rpath//date//"-"//time//"-uranus_"//trim(files(i)), status="new")
     write(11, '(6(e16.5))') (tidal(k,1), tidal(k,2), tidal(k,3), k = 1, j)
     close(11)

     deallocate(sd, tidal, data)

  end do

end program helled_11_13
