module aux_gavrilov

  use iso_fortran_env

  use mod_consts, only : PI, G

  implicit none

contains

  real(real64) pure function lin_rho(rho, s, s1)
    ! arguments
    real(real64), intent(in) :: rho, s, s1

    ! locals
    real(real64) :: x

    x = s / s1

    lin_rho = 4*rho*(1 - x)

  end function lin_rho

  real(real64) pure function lin_g(rho, s, s1)
    ! arguments
    real(real64), intent(in) :: rho, s, s1

    ! locals
    real(real64) :: x

    x = s / s1

    lin_g = 16./3.*PI*G * rho*s * (1 - 3*x / 4)

  end function lin_g

  real(real64) pure function quad_rho(rho, s, s1)
    ! arguments
    real(real64), intent(in) :: rho, s, s1

    ! locals
    real(real64) :: x

    x = s / s1

    quad_rho = 2.5*rho * (1 - x**2)

  end function quad_rho

  real(real64) pure function quad_g(rho, s, s1)
    ! arguments
    real(real64), intent(in) :: rho, s, s1

    ! locals
    real(real64) :: x

    x = s / s1

    quad_g = 10*PI*G * rho*s * (1./3 - x**2/5)

  end function quad_g

end module aux_gavrilov
