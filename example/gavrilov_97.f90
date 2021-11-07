program gavrilov_97

  use iso_fortran_env

  use mod_consts
  use mod_math
  use mod_structure
  use mod_tidal, only : propagator_method, normalize, use_GD_core

  use aux_gavrilov

  implicit none

  ! parameters for benchmark
  integer(int64), parameter :: n = 10
  integer(int64), parameter :: ln = 5
  real(real64), parameter :: ri = 0.1
  real(real64), parameter :: rf = 69.911e6
  real(real64), parameter :: rhoj = 1326
  character(*), parameter :: path = "../post/processed/gavrilov_97/"

  ! locals
  procedure(dy), pointer :: df => null()
  complex(real64) :: sd(n,4), tidal(n,3)
  real(real64) :: mu, mass, zero, eta, omega
  integer(int64) :: i, j
  real(real64), allocatable :: l(:)
  logical :: flag
  character(8)  :: date
  character(10)  :: time

  allocate(l(ln))
  l = [2, 3, 4, 5, 6]

  zero = 0

  df => planet_m
  mass = 0
  mu = 1
  eta = 0
  omega = 0

  ! LINEAR CALCULATION
  sd(:,1) = linspace(ri, rf, n)
  sd(:,2) = mu

  ! check which core to use
  call use_GD_core(ri, mu, flag)

  print *, "LINEAR STRUCT FILE"
  call date_and_time(DATE=date, TIME=time)
  open(10, file=path//date//"-"//time//"-gavrilov_97_lin.txt", status="new")
  do i = 1, n

     sd(i,4) = lin_rho(rhoj, sd(i,1)%re, rf)

     if(i == 1) then
        mass = simpson(df, zero, sd(i,1)%re, sd(i,4)%re, n)
     else
        mass = mass + simpson(df, sd(i-1,1)%re, sd(i,1)%re, sd(i,4)%re, n)
     end if

     sd(i,3) = lin_g(rhoj, sd(i,1)%re, rf)
     sd(i,2) = planet_cmu(mu, omega, eta, sd(i,1)%re, sd(i,3)%re, sd(i,4)%re)
     write(10, '(4(e16.5))') sd(i,1)%re, sd(i,2)%re, sd(i,3)%re, sd(i,4)%re

  end do
  close(10)

  call normalize(mass, sd(n,1)%re, sd)

  print *, "LINEAR TIDAL CALC"
  do i = 1, ln

     call date_and_time(DATE=date, TIME=time)
     open(10, file=path//date//"-"//time//"-tidal_lin_radius.txt", status="new")

     tidal = propagator_method(l(i), n, sd, flag)
     write(10, '(3(e16.5))') (tidal(j,1)%re, tidal(j,2)%re, tidal(j,3)%re, j = 1, n)
     close(10)

  end do

  ! QUADRATIC CALCULATION
  sd(:,1) = linspace(ri, rf, n)
  sd(:,2) = mu

  print *, "QUAD STRUCT FILE"
  call date_and_time(DATE=date, TIME=time)
  open(10, file=path//date//"-"//time//"-gavrilov_97_quad.txt", status="new")
  do i = 1, n

     sd(i,4) = quad_rho(rhoj, sd(i,1)%re, rf)

     if(i == 1) then
        mass = simpson(df, zero, sd(i,1)%re, sd(i,4)%re, n)
     else
        mass = mass + simpson(df, sd(i-1,1)%re, sd(i,1)%re, sd(i,4)%re, n)
     end if

     sd(i,3) = quad_g(rhoj, sd(i,1)%re, rf)
     sd(i,2) = planet_cmu(mu, omega, eta, sd(i,1)%re, sd(i,3)%re, sd(i,4)%re)
     write(10, '(4(e16.5))') sd(i,1)%re, sd(i,2)%re, sd(i,3)%re, sd(i,4)%re

  end do
  close(10)

  call normalize(mass, sd(n,1)%re, sd)

  print *, "QUAD TIDAL CALC"
  do i = 1, ln

     call date_and_time(DATE=date, TIME=time)
     open(10, file=path//date//"-"//time//"-tidal_quad_radius.txt", status="new")

     tidal = propagator_method(l(i), n, sd, flag)
     write(10, '(3(e16.5))') (tidal(j,1)%re, tidal(j,2)%re, tidal(j,3)%re, j = 1, n)
     close(10)

  end do

end program gavrilov_97
