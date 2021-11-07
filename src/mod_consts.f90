!==============================================================================
! Author: David A. James, davidabraham@ucla.edu
!
! Background: The following is a module of constants and precision values for
! the rest of the project.
!
!==============================================================================
module mod_consts

  use iso_fortran_env

  implicit none

  ! general constants
  real(real64), parameter :: PI = 4.0*atan(1.0)      ! pi calculated out
  real(real64), parameter :: G = 6.67408e-11         ! gravitational constant [N*m^2/kg^2]

  ! conversions
  real(real64), parameter :: PaToMBar = 1e-11        ! conversion from Pa to MBar
  real(real64), parameter :: MBarToPa = 1e11         ! conversion from MBar to Pa
  real(real64), parameter :: yrToSec = 365*24*60*60

end module mod_consts
