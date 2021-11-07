program check_cmu

  use iso_fortran_env

  use mod_consts
  use mod_math
  use mod_structure, only : cmu_maxwell, cmu_SLS, cmu_andrade, planet_g, planet_mu

  implicit none

  ! parameters for benchmark
  integer(int64), parameter :: layers = 1000
  real(real64), parameter :: mass = 5.97e24
  real(real64), parameter :: a = 152.1e9
  real(real64), parameter :: rho = 5515
  real(real64), parameter :: R = 6371e3
  real(real64), parameter :: xi = 1e-20
  real(real64), parameter :: xf = 1e20
  character(*), parameter :: pathret = "../post/processed/cmu_runs/"

  ! locals
  complex(real64) :: cmu(layers,3)
  real(real64) :: mu, omega, eta, n, gl, p
  real(real64) :: mu_factor, xspace(layers)
  integer(int64) :: i
  character(8)  :: date
  character(10)  :: time

  mu_factor = 60
  eta = 2e21
  gl = planet_g(R, mass)
  mu = planet_mu(R, gl, rho)

  p = G * mass**2 / r**4

  xspace = logspace(xi, xf, layers)

  do i = 1, layers

     cmu(i,1) = cmu_maxwell(mu, xspace(i), eta)
     cmu(i,2) = cmu_SLS(mu, xspace(i), eta, mu_factor)
     cmu(i,3) = cmu_andrade(mu, xspace(i), eta)

  end do

  cmu(:,:) = cmu(:,:) / p

  call date_and_time(DATE=date, TIME=time)
  open(10, file=pathret//date//"-"//time//"-cmu_omega", status="new")
  write(10, '(7(e16.5))') (xspace(i), cmu(i,1)%re, cmu(i,1)%im, cmu(i,2)%re, cmu(i,2)%im, cmu(i,3)%re, cmu(i,3)%im, i = 1, layers)
  close(10)

end program check_cmu
