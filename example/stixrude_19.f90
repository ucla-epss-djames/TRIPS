program stixrude_19

  use iso_fortran_env

  use mod_consts
  use mod_math
  use mod_structure
  use mod_tidal, only : propagator_method, normalize, use_GD_core

  use aux_stixrude

  implicit none

  ! parameters for benchmark
  integer(int64), parameter :: n = 100
  integer(int64), parameter :: nc = .1*n
  integer(int64), parameter :: nm = n - nc
  real(real64), parameter :: l = 2
  real(real64), parameter :: Qc = 1e4
  real(real64), parameter :: Re = 6371e3
  real(real64), parameter :: Me = 5.972e24
  real(real64), parameter :: M = 6.55*Me
  real(real64), parameter :: ri = 0.001
  real(real64), parameter :: rf = 2.68*Re
  real(real64), parameter :: OMEGA_U = 2.*pi/(17.24*3600.)
  real(real64), parameter, dimension(2) :: ratio = [0.73, 0.01]
  character(*), parameter :: path = "../post/processed/stixrude_19/"

  ! locals
  complex(real64) :: sd(n,4), tidal(n,3), cmu
  real(real64) :: data(n,2), zero, omega, tau_m, eta, mu, mass, a, rho, rho_m, rho_c, rc, dr
  integer(int64) :: i, j
  logical :: flag
  character(8)  :: date
  character(10)  :: time

  rho = M / vol(rf)
  zero = 0
  omega = 0
  mu = 0
  eta = 0

  data(:,:) = 0

  ! check which core to use
  call use_GD_core(zero, mu, flag)

  print *, "FLUID STRUCTURE"
  do i = 1, 2

     rho_m = ratio(i)*rho

     call date_and_time(DATE=date, TIME=time)
     open(10, file=path//date//"-"//time//"-tidal_fluid_core.txt", status="new")
     do j = 1, n

        rc = rf * (j / (1. * n))

        sd(1:nc,1) = linspace(ri, rc, nc)
        dr = sd(2,1)%re - sd(1,1)%re
        sd(nc+1:n,1) = linspace(rc + dr, rf, nm)

        mass = rho_m * (vol(rf) - vol(rc))
        rho_c = (M - mass) / vol(rc)

        sd(1:nc,4) = rho_c
        sd(nc+1:n,4) = rho_m

        mass = 0

        call planet_structure(n, omega, data, mass, sd)

        ! saving values for output
        cmu = sd(nc,2)
        a = sd(nc,3)%re

        call normalize(mass, sd(n,1)%re, sd)

        tidal = propagator_method(l, n, sd, flag)
        write(10, '(8(e16.5))') sd(nc,1)%re, tidal(n,1)%re, tidal(n,1)%im, rho_c, rho_m, cmu%re, cmu%im, a

     end do
     close(10)

  end do

  omega = OMEGA_U

  print *, "VISCOELASTIC STRUCTURE"
  do i = 1, 2

     rho_m = ratio(i)*rho

     call date_and_time(DATE=date, TIME=time)
     open(10, file=path//date//"-"//time//"-tidal_viscoelastic_core.txt", status="new")
     do j = 1, n

        rc = rf * (j / (1. * n))

        sd(1:nc,1) = linspace(ri, rc, nc)
        dr = sd(2,1)%re - sd(1,1)%re
        sd(nc+1:n,1) = linspace(rc + dr , rf, nm)

        mass = rho_m * (vol(rf) - vol(rc))
        rho_c = (M - mass) / vol(rc)

        sd(1:nc,4) = rho_c
        sd(nc+1:n,4) = rho_m

        a = planet_g(rc, M - mass)
        mu = planet_mu(rc, a, rho_c)
        tau_m = Qc / omega
        eta = mu * tau_m

        data(:,1) = mu
        data(1:nc,2) = eta
        data(nc+1,2) = 0.

        ! check which core to use
        call use_GD_core(ri, mu, flag)
        mass = 0

        call planet_structure(n, omega, data, mass, sd)

        ! saving values for output
        cmu = sd(nc,2)
        a = sd(nc,3)%re

        call normalize(mass, sd(n,1)%re, sd)

        tidal = propagator_method(l, n, sd, flag)
        write(10, '(8(e16.5))') sd(nc,1)%re, tidal(n,1)%re, tidal(n,1)%im, rho_c, rho_m, cmu%re, cmu%im, a

     end do
     close(10)

  end do

end program stixrude_19
