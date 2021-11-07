module aux_stixrude

  use iso_fortran_env

  use mod_consts

  implicit none

contains

  real(real64) function vol(r) &
       result(V)
    ! arguments
    real(real64), intent(in) :: r

    V = 4./3.*PI * r**3

  end function vol

end module aux_stixrude
