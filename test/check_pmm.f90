program check_pmm

  use iso_fortran_env

  use mod_consts
  use mod_math
  use mod_structure
  use mod_tidal, only : M, N, normalize, aggregate_matrix, use_GD_core

  implicit none

  ! parameters for benchmark
  integer(int64), parameter :: h = 100
  integer(int64), parameter :: layers = 4
  real(real64), parameter :: l = 2
  real(real64), parameter :: omega = 1. / (24. * 3600)
  character(*), parameter :: pathret = "../post/processed/billsjames_97/"
  character(*), parameter :: pathdat = "../data/billsjames_97"

  ! locals
  procedure(dy), pointer :: df => null()
  complex(real64) :: sd(layers,4), tidal(layers,3)
  real(real64), dimension(layers) :: r, rho, mu, eta
  real(real64) :: mass, zero
  integer(int64) :: i, j, k
  logical :: flag
  character(8)  :: date
  character(10)  :: time

  ! locals for PMM
  complex(real64), dimension(M, N) :: B
  complex(real64), dimension(M, N, layers) :: Bi
  complex(real64), dimension(N) :: bound
  complex(real64), dimension(M) :: yvec

  zero = 0

  df => planet_m
  mass = 0

  ! opening data
  open(10, file=pathdat, status="old")
  read(10, *) (r(i), rho(i), mu(i), eta(i), i = 1, layers)
  close(10)

  ! check which core to use
  call use_GD_core(r(1), mu(1), flag)

  call date_and_time(DATE=date, TIME=time)
  do i = 1, layers

     if(i == 1) then
        mass = simpson(df, zero, r(i), rho(i), h)
     else
        mass = mass + simpson(df, r(i-1), r(i), rho(i), h)
     end if

     sd(i,1) = r(i)
     sd(i,3) = planet_g(r(i), mass)
     sd(i,4) = rho(i)
     sd(i,2) = planet_cmu(mu(i), omega, eta(i), sd(i,1)%re, sd(i,3)%re, sd(i,4)%re)

  end do

  call normalize(mass, sd(layers,1)%re, sd)

  call aggregate_matrix(l, layers, sd, flag, B, Bi)

  ! print *, "PRINTING B"
  ! do i = 1, layers
  !    print *, "i:", i
  !    do j = 1, M
  !       print '(99e21.5)', (Bi(j,k,i), k = 1, N)
  !    end do
  ! end do

end program check_pmm
