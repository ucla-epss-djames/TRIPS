program prem

  use iso_fortran_env

  use mod_consts
  use mod_math
  use mod_structure
  use mod_tidal, only : propagator_method, normalize, use_GD_core

  implicit none

  ! parameters for benchmark
  real(real64), parameter :: l = 2
  integer(int64), parameter :: index = 6
  integer(int64), parameter, dimension(index) :: NMAX = [6372, 2, 6, 11, 51, 2892]
  character(*), parameter, dimension(index) :: files = ["      ", "_L1   ", "_L5   ", "_L10  ", "_L50  ", "_L2891"]
  character(*), parameter :: dpath = "../data/prem"
  character(*), parameter :: rpath = "../post/processed/prem/"

  ! locals
  procedure(dy), pointer :: df => null()
  complex(real64), allocatable :: sd(:,:), tidal(:,:)
  real(real64), allocatable :: read_r(:), read_rho(:), read_mu(:), read_eta(:)
  real(real64), allocatable :: mu(:), eta(:)
  real(real64) :: vals(5), zero, mass, omega, check
  integer(int64) :: i, j, k, n
  logical :: flag
  character(len=:), allocatable :: path
  character(8)  :: date
  character(10)  :: time

  zero = 0
  check = 1
  omega = 1./(24.*3600.)        !  Rotation rate Earth

  df => planet_m

  do i = 1, 1

     path = trim(dpath//files(i))
     open(10, file=path, status="old")
     mass = 0

     allocate(read_r(NMAX(i)), read_rho(NMAX(i)), read_mu(NMAX(i)), read_eta(NMAX(i)))

     do j = 1, NMAX(i)

        read(10,*) vals(1), vals(2), vals(3), vals(4), vals(5)

        read_r(j) = vals(1)
        read_rho(j) = vals(2)
        read_mu(j) = vals(4)
        read_eta(j) = vals(5)

     end do

     close(10)

     allocate(sd(NMAX(i),4), tidal(NMAX(i),3), mu(NMAX(i)), eta(NMAX(i)))

     ! struct files are backwards of how TRIPS reads in data
     sd(:,1) = read_r(NMAX(i):1:-1)
     sd(:,4) = read_rho(NMAX(i):1:-1)
     mu = read_mu(NMAX(i):1:-1)
     eta = read_eta(NMAX(i):1:-1)

     if(i .eq. 1) then
        sd(1,1) = 0.001
     end if

     ! check which core to use
     call use_GD_core(sd(1,1)%re, zero, flag)
     print *, flag

     do j = 1, NMAX(i)

        if(j == 1) then
           mass = simpson(df, zero, sd(j,1)%re, sd(j,4)%re, NMAX(i))
        else
           mass = mass + simpson(df, sd(j-1,1)%re, sd(j,1)%re, sd(j,4)%re, NMAX(i))
        end if

        sd(j,3) = planet_g(sd(j,1)%re, mass)
        sd(j,2) = planet_cmu(mu(j), omega, eta(j), sd(j,1)%re, sd(j,3)%re, sd(j,4)%re)
     end do

     call normalize(mass, sd(NMAX(i),1)%re, sd)

     call date_and_time(DATE=date, TIME=time)
     open(11, file=rpath//date//"-"//time//"-prem_struct"//trim(files(i)), status="new")
     write(11, '(5(e16.5))') (sd(j,1)%re, sd(j,2)%re, sd(j,2)%im, sd(j,3)%re, sd(j,4)%re, j = 1, NMAX(i))
     close(11)

     tidal = propagator_method(l, NMAX(i), sd, flag)
     open(11, file=rpath//date//"-"//time//"-prem"//trim(files(i)), status="new")
     write(11, '(6(e16.5))') (tidal(j,1), tidal(j,2), tidal(j,3), j = 1, NMAX(i))
     close(11)

     deallocate(read_r, read_rho, read_mu, read_eta, sd, tidal, mu, eta)
  end do

end program prem
