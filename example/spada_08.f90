program spada_08

  use iso_fortran_env

  use mod_consts
  use mod_math
  use mod_structure
  use mod_tidal, only : propagator_method, normalize, use_GD_core

  implicit none

  ! parameters for benchmark
  integer(int64), parameter :: n = 1000
  integer(int64), parameter :: ln = 18
  integer(int64), parameter :: layers = 4
  real(real64), parameter, dimension(2) :: t = [1e-3, 1e9]
  character(*), parameter :: pathret = "../post/processed/spada_08/"
  character(*), parameter :: pathdat = "../data/billsjames_97"

  ! locals
  procedure(dy), pointer :: df => null()
  complex(real64) :: sd(layers,4), tidal(layers,3)
  real(real64), dimension(layers) :: r, rho, mu, eta
  real(real64) :: mass, zero, omega(2), l(ln)
  integer(int64) :: i
  logical :: flag
  character(8)  :: date
  character(10)  :: time

  l = [2, 3, 4, 5, 6, 8, 10, 20, 30, 40, 50, 60, &
       80, 100, 200, 300, 500, 700]

  zero = 0

  df => planet_m
  mass = 0

  omega = 2.*PI / (t * yrToSec)

  ! opening data
  open(10, file=pathdat, status="old")
  read(10, *) (r(i), rho(i), mu(i), eta(i), i = 1, layers)
  close(10)

  ! check which core to use
  call use_GD_core(r(1), mu(1), flag)

  print *, "ELASTIC STRUCT FILE"
  open(10, file=pathret//"spada_08_elastic.txt", status="new")
  do i = 1, layers

     if(i == 1) then
        mass = simpson(df, zero, r(i), rho(i), n)
     else
        mass = mass + simpson(df, r(i-1), r(i), rho(i), n)
     end if

     sd(i,1) = r(i)
     sd(i,3) = planet_g(r(i), mass)
     sd(i,4) = rho(i)
     sd(i,2) = planet_cmu(mu(i), omega(1), eta(i), sd(i,1)%re, sd(i,3)%re, sd(i,4)%re)
     write(10, '(4(e16.5))') sd(i,1)%re, sd(i,2)%re, sd(i,3)%re, sd(i,4)%re

  end do
  close(10)

  call normalize(mass, sd(layers,1)%re, sd)

  print *, "ELASTIC CALC"
  call date_and_time(DATE=date, TIME=time)
  open(10, file=pathret//date//"-"//time//"-tidal_elastic.txt", status="new")
  do i = 1, ln

     tidal = propagator_method(l(i), layers, sd, flag)
     write(10, '(4(e16.5))') l(i), tidal(layers,1)%re, tidal(layers,2)%re, tidal(layers,3)%re
  end do
  close(10)

  print *, "FLUID STRUCT FILE"
  open(10, file=pathret//"spada_08_fluid.txt", status="new")
  do i = 1, layers

     if(i == 1) then
        mass = simpson(df, zero, r(i), rho(i), n)
     else
        mass = mass + simpson(df, r(i-1), r(i), rho(i), n)
     end if

     sd(i,1) = r(i)
     sd(i,3) = planet_g(r(i), mass)
     sd(i,4) = rho(i)
     sd(i,2) = planet_cmu(mu(i), omega(2), eta(i), sd(i,1)%re, sd(i,3)%re, sd(i,4)%re)
     write(10, '(4(e16.5))') sd(i,1)%re, sd(i,2)%re, sd(i,3)%re, sd(i,4)%re

  end do
  close(10)

  call normalize(mass, sd(layers,1)%re, sd)

  print *, "FLUID CALC"
  call date_and_time(DATE=date, TIME=time)
  open(10, file=pathret//date//"-"//time//"-tidal_fluid.txt", status="new")
  do i = 1, ln

     tidal = propagator_method(l(i), layers, sd, flag)
     write(10, '(4(e16.5))') l(i), tidal(layers,1)%re, tidal(layers,2)%re, tidal(layers,3)%re
  end do
  close(10)

end program spada_08
