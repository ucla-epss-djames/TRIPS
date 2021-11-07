program helled_13

  use iso_fortran_env

  use mod_consts
  use mod_math
  use mod_structure
  use mod_tidal, only : propagator_method, normalize, use_GD_core

  implicit none

  ! parameters for benchmark
  real(real64), parameter :: l = 2
  real(real64), parameter :: ME = 5.9736e27
  real(real64), parameter :: RE = 6380e3
  real(real64), parameter :: OMEGA = 1.012e-4
  real(real64), parameter :: gccm_to_kgm = 1000
  integer(int64), parameter, dimension(2) :: NMAX = [1061, 1285]
  character(*), parameter, dimension(2) :: files = ["1a.dat", "2b.dat"]
  character(*), parameter :: dpath = "../data/table_N"
  character(*), parameter :: rpath = "../post/processed/helled_13/"

  ! locals
  complex(real64), allocatable :: sd(:,:), tidal(:,:)
  real(real64), allocatable  :: r(:), rho(:), data(:,:)
  real(real64) :: vals(5), mass
  integer(int64) :: i, j
  character(len=:), allocatable :: path
  logical :: flag
  character(8)  :: date
  character(10)  :: time

  do i = 1, 2

     path = trim(dpath//files(i))
     open(10, file=path, status="old")
     mass = 0

     allocate(r(NMAX(i)), rho(NMAX(i)), sd(NMAX(i),4), tidal(NMAX(i),3), data(NMAX(i),2))

     do j = 1, NMAX(i)


        read(10,*) vals(1), vals(2), vals(3), vals(4), vals(5)

        r(j) = vals(3) * RE
        rho(j) = vals(5) * gccm_to_kgm

     end do

     close(10)

     sd(:,1) = r(NMAX(i):1:-1)
     sd(:,4) = rho(NMAX(i):1:-1)

     data(:,:) = 0

     ! check which core to use
     call use_GD_core(sd(1,1)%re, data(1,1), flag)
     mass = 0

     call planet_structure(NMAX(i), OMEGA, data, mass, sd)

     call date_and_time(DATE=date, TIME=time)
     open(11, file=rpath//date//"-"//time//"-neptune_struct_"//trim(files(i)), status="new")
     write(11, '(4(e16.5))') (sd(j,1)%re, sd(j,2)%re, sd(j,3)%re, sd(j,4)%re, j = 1,NMAX(i))
     close(11)

     call normalize(mass, sd(NMAX(i),1)%re, sd)

     tidal = propagator_method(l, NMAX(i), sd, flag)
     open(11, file=rpath//date//"-"//time//"-neptune_"//trim(files(i)), status="new")
     write(11, '(3(e16.5))') (tidal(j,1)%re, tidal(j,2)%re, tidal(j,3)%re, j = 1,NMAX(i))
     close(11)

     deallocate(r, rho, sd, tidal, data)

  end do

end program helled_13
