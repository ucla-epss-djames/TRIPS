program vermeersen_96

  use iso_fortran_env

  use mod_consts
  use mod_math
  use mod_structure
  use mod_tidal, only : propagator_method, normalize, use_GD_core

  implicit none

  ! parameters for benchmark
  integer(int64), parameter :: n = 1000
  integer(int64), parameter :: ln = 3
  integer(int64), parameter :: layers = 5
  real(real64), parameter :: ti = 1e-5
  real(real64), parameter :: tf = 1e5
  character(*), parameter :: pathret = "../post/processed/vermeersen_96/"
  character(*), parameter :: pathdat = "../data/vermeersen_96.txt"

  ! locals
  procedure(dy), pointer :: df => null()
  complex(real64) :: sd(layers,4), sdsave(layers,4)
  complex(real64) :: tidal(layers,3)
  real(real64), dimension(layers) :: r, rho, mu, eta
  real(real64) :: t(n), omega(n), mass, zero, l(ln)
  integer(int64) :: i, j, k
  logical :: flag
  character(8)  :: date
  character(10)  :: time

  zero = 0

  l = [2, 10, 100]

  t = logspace(ti, tf, n)
  omega = 2.*PI / (t * yrToSec)

  df => planet_m
  mass = 0

  ! opening data
  open(10, file=pathdat, status="old")
  read(10, *) (r(i), rho(i), mu(i), eta(i), i = 1, layers)
  close(10)

  ! check which core to use
  call use_GD_core(r(1), mu(1), flag)

  ! calculating planetary values
  do i = 1, layers

     if(i == 1) then
        mass = simpson(df, zero, r(i), rho(i), n)
     else
        mass = mass + simpson(df, r(i-1), r(i), rho(i), n)
     end if

     sd(i,1) = r(i)
     sd(i,3) = planet_g(r(i), mass)
     sd(i,4) = rho(i)

  end do

  sdsave = sd

  ! love number index
  do i = 1, ln

     call date_and_time(DATE=date, TIME=time)
     open(10, file=pathret//date//"-"//time//"-tidal_omega.txt", status="new")

     ! omega index
     do j = 1, n

        ! producing structure for each omega value
        sd = sdsave
        do k = 1, layers
           sd(k,2) = planet_cmu(mu(k), omega(j), eta(k), sd(k,1)%re, sd(k,3)%re, sd(k,4)%re)
        end do
        call normalize(mass, sd(layers,1)%re, sd)

        ! feeding structure into prop method
        tidal = propagator_method(l(i), layers, sd, flag)
        write(10, '(4(e16.5))') omega(j), tidal(layers,1)%re, tidal(layers,2)%re, tidal(layers,3)%re

     end do
     close(10)

  end do

end program vermeersen_96
