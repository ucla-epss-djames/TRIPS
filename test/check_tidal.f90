program check_tidal

  use iso_fortran_env

  use mod_consts
  use mod_math
  use mod_structure
  use mod_tidal, only : propagator_method, normalize, use_GD_core

  implicit none

  ! parameters for benchmark
  integer(int64), parameter :: layers = 4
  real(real64), parameter :: l = 2
  real(real64), parameter :: OMEGA = 1. / (24. * 3600)
  character(*), parameter :: pathret = "../post/processed/billsjames_97/"
  character(*), parameter :: pathdat = "../data/billsjames_97"

  ! locals
  complex(real64) :: sd(layers,4), tidal(layers,3)
  real(real64), dimension(layers) :: r, rho
  real(real64) :: mass, data(layers,2)
  integer(int64) :: i
  logical :: flag
  character(8)  :: date
  character(10)  :: time

  sd(:,:) = 0

  ! opening data
  open(10, file=pathdat, status="old")
  read(10, *) (sd(i,1)%re, sd(i,4)%re, data(i,1), data(i,2), i = 1, layers)
  close(10)

  ! check which core to use
  call use_GD_core(r(1), data(1,1), flag)

  call date_and_time(DATE=date, TIME=time)

  open(10, file=pathret//date//"-"//time//"-billsjames_97_struct", status="new")
  call planet_structure(layers, OMEGA, data, mass, sd)
  write(10, '(5(e16.5))') (sd(i,1)%re, sd(i,2)%re, sd(i,2)%im, sd(i,3)%re, sd(i,4)%re, i = 1, layers)
  close(10)

  call normalize(mass, sd(layers,1)%re, sd)

  open(11, file=pathret//date//"-"//time//"-billsjames_97_norm", status="new")
  write(11, '(5(e16.5))') (sd(i,1)%re, sd(i,2)%re, sd(i,2)%im, sd(i,3)%re, sd(i,4)%re, i = 1, layers)
  close(11)

  open(10, file=pathret//date//"-"//time//"-billsjames_97_tidal", status="new")
  tidal = propagator_method(l, layers, sd, flag)
  write(10, '(6(e16.5))') (tidal(i,1)%re, tidal(i,1)%im, tidal(i,2)%re, tidal(i,2)%im, tidal(i,3)%re, tidal(i,3)%im, i = 1, layers)
  close(10)

end program check_tidal
