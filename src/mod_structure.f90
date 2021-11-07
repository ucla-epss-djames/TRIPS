!==============================================================================
! Author: David A. James, davidabraham@ucla.edu
!
! Background: The following is the planetary structure code containing the ODE
! system to generate a planet by mass and radius. Along with functions to
! calculate other planetary values. Maxwell rheology follow comments on
! Storch & Lai '13.
!
! Module: mod_structure
! * Governing Equations
!   - function dmdr : conservation of mass
!   - function dPdr : hydrostatic equilibrium
!   - function eos  : equation of state
! * Auxiliary Functions
!   - function planet_form        : ODE system of GEs
!   - subroutine planet_structure : structure data for a planet
!   - function planet_m           : mass in a spherical shell
!   - function planet_g           : gravity at a distance r
!   - function planet_mu          : shear modulus
!   - function planet_cmu         : complex shear modulus
!   - function planet_mmotion     : mean motion
!   - function planet_thermo      : collects thermo values
!
!==============================================================================
module mod_structure

  use iso_fortran_env

  use mod_consts, only : PI, G
  use mod_math, only : hunt_loc, poly_int, simpson, dy

  implicit none

  private
  public :: planet_form, planet_g, planet_m, planet_mu, planet_cmu, planet_mmotion, planet_structure, &
       cmu_maxwell, cmu_SLS, cmu_andrade

  complex(real64), parameter :: z = (0, 1)

contains

  ! ***** GOVERNING EQUATIONS *****

  real(real64) pure function dmdr(r, rho)
    !-------------------------------------------------------------------------
    ! dmdr : conservation of mass
    !
    ! Outputs the change in mass per change in distance.
    !
    ! Arguments:
    !   r   - real - radius [m]
    !   rho - real - density [kg / m^3]
    ! Result:
    !   dmdr - real - change in mass [kg / m]
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: r, rho

    dmdr = 4*PI * rho * r**2

  end function dmdr

  real(real64) pure function dPdr(r, rho, m)
    !-------------------------------------------------------------------------
    ! dPdr : hydrostatic equilibrium
    !
    ! Outputs the change in pressure per change in distance.
    !
    ! Arguments:
    !   r   - real - radius [m]
    !   rho - real - density [kg / m^3]
    !   m   - real - mass [kg]
    ! Result:
    !   dPdr - real - change in pressure [Pa / m]
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: r, rho, m

    dPdr = -1*G * m * rho / r**2

  end function dPdr

  real(real64) pure function eos(P, td_P, td_d) &
       result(rho)
    !-------------------------------------------------------------------------
    ! eos : equation of state
    !
    ! Outputs the density at a given pressure. Based on fort.56
    !
    ! Arguments:
    !   P    - real    - pressure [Pa]
    !   td_P - real(:) - pressure values from thermo model [Pa]
    !   td_d - real(:) - density values from thermo model [kg / m^3]
    ! Result:
    !   rho - real - density [kg / m^3]
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: P
    real(real64), intent(in) :: td_P(:), td_d(:)

    ! locals
    integer(int64) :: i, len

    len = size(td_P)
    i = hunt_loc(td_P, P)

    ! need exit flag for when 0 happens
    if(i == 0 .or. i == 1) then
       i = 1
       rho = poly_int(P, td_P(i:i+2), td_d(i:i+2))
    else if(i == len) then
       rho = poly_int(P, td_P(i-2:i), td_d(i-2:i))
    else
       rho = poly_int(P, td_P(i-1:i+1), td_d(i-1:i+1))
    end if

  end function eos

  ! ***** AUXILARY FUNCTIONS *****

  pure function planet_form(r, u, td) &
       result(du)
    !-------------------------------------------------------------------------
    ! planet_form : ODE system of GEs
    !
    ! The ODE system of the GEs where it returns the mass and radius of the
    ! given radius.
    !
    ! Arguments:
    !   rn - real      - radius [m]
    !   u  - real(2)   - IVs for ODE [Pa, kg]
    !   td - real(:,:) - thermo data
    ! Result:
    !   du - real(2) - outputted mass and pressure [Pa, kg]
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: r
    real(real64), intent(in) :: u(:)
    real(real64), intent(in) :: td(:,:)
    real(real64) :: du(size(u))

    ! locals
    real(real64) :: rho

    rho = eos(u(1), td(:,1), td(:,2))

    du(1) = dPdr(r, rho, u(2))
    du(2) = dmdr(r, rho)

  end function planet_form

  subroutine planet_structure(layers, omega, data, mass, sd)
    !-------------------------------------------------------------------------
    ! planet_structure : structure data for a planet
    !
    ! Generates structure data of a planet for tidal code.
    !
    ! Arguments:
    !   layers - real           - layers of planet
    !   omega  - real           - frequency
    !   data   - real(layers,2) - shear and viscosity values
    ! Results:
    !   mass - real              - total mass of planet
    !   sd   - complex(layers,4) - structure data
    !-------------------------------------------------------------------------
    ! parameters
    integer(int64), parameter :: n = 500
    real(real64), parameter :: zero = 0

    ! arguments
    integer(int64), intent(in) :: layers
    real(real64), intent(in) :: omega, data(layers,2)
    real(real64), intent(in out) :: mass
    complex(real64), intent(in out) :: sd(layers,4)

    ! locals
    real(real64) :: rho, mu, eta
    procedure(dy), pointer :: df => null()
    integer(int64) :: i

    mass = 0
    df => planet_m

    do i = 1, layers

       rho = sd(i,4)%re
       mu = data(i,1)
       eta = data(i,2)

       if(i == 1) then
          mass = simpson(df, zero, sd(i,1)%re, rho, n)
       else
          mass = mass + simpson(df, sd(i-1,1)%re, sd(i,1)%re, rho, n)
       end if

       sd(i,3) = planet_g(sd(i,1)%re, mass)
       sd(i,2) = planet_cmu(mu, omega, eta, sd(i,1)%re, sd(i,3)%re, rho)

    end do

  end subroutine planet_structure

  real(real64) pure function planet_m(r, rho) &
       result(m)
    !-------------------------------------------------------------------------
    ! planet_m : mass in a spherical shell
    !
    ! The mass of planet inside a sphere of radius r. This is used in
    ! conjunction with an integration method to capture the mass correctly.
    !
    ! Arguments:
    !   r   - real - radius [m]
    !   rho - real - density [kg / m^3]
    ! Result:
    !   m - real - mass [kg]
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: r, rho

    m = 4*PI * rho * r**2

  end function planet_m

  real(real64) pure function planet_g(r, m) &
       result(a)
    !-------------------------------------------------------------------------
    ! planet_g : gravity at a distance r
    !
    ! The gravity calculated at a distance r containing a mass m.
    !
    ! Arguments:
    !   r - real - radius [m]
    !   m - real - mass [kg]
    ! Result:
    !   a - real - acceleration [m / s^2]
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: r, m

    ! locals
    real(real64) :: r_center

    r_center = 0

    if(r == r_center) then
       a = 0
    else
       a = G * m / r**2
    end if

  end function planet_g

  real(real64) pure function planet_mu(r, a, rho) &
       result(mu)
    !-------------------------------------------------------------------------
    ! planet_mu : shear modulus
    !
    ! The shear modulus calculated using the density.
    !
    ! Arguments:
    !   r   - real - radius [m]
    !   a   - real - acceleration [m / s^2]
    !   rho - real - density [kg / m^3]
    ! Result:
    !   mu - real - shear modulus [Pa]
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: r, a, rho

    mu = rho * a * r

  end function planet_mu

  complex(real64) function planet_cmu(mu, omega, eta, r, a, rho) &
       result(cmu)
    !-------------------------------------------------------------------------
    ! planet_cmu : complex shear modulus
    !
    ! The shear modulus calculated using the density.
    !
    ! Arguments:
    !   mu    - real - shear modulus [Pa]
    !   omega - real - frequency [s^-1]
    !   eta   - real - viscosity [kg / m^3]
    !   r     - real - radius [Pa s]
    !   a     - real - acceleration [m / s^2]
    !   rho   - real - density [kg / m^3]
    ! Result:
    !   mu - complex - complex shear modulus [Pa]
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: mu, omega, eta, r, a, rho

    ! locals
    real(real64) :: omega_m, smu, alpha, dem, mu_factor, mass, orbit_axis, n
    real(real64) :: eta_zero, eta_inf, mu_zero

    mu_factor = 60

    eta_zero = 0
    eta_inf = -1
    mu_zero = 0

    if(eta == eta_zero) then
       ! if eta is small
       ! produce a small shear modulus so tidal code isn't unstable

       smu = 1e-4 * planet_mu(r, a, rho)
       cmu = cmplx(smu, 0, real64)

    else if(eta == eta_inf) then
       ! if eta is infinite
       ! produce the real shear modulus

       cmu = cmplx(mu, 0, real64)

    else
       ! else use the model rheaology
       ! models: Maxwell, SLS, andrade(TBD)

       cmu = cmu_maxwell(mu, omega, eta)
       ! cmu = cmu_SLS(mu, omega, eta, mu_factor)

       mass = 5.97e24
       orbit_axis = 152.1e9

       ! mass = 86.8e24
       ! orbit_axis = 2741.3e9

       ! n = planet_mmotion(mass, orbit_axis)
       ! cmu = cmu_andrade(mu, omega, eta, n)

    end if

    if(cmu%re == mu_zero) then
       ! final check if shear is 0

       cmu = cmplx(1e-9, 0, real64)
    end if

  end function planet_cmu

  complex(real64) function cmu_maxwell(mu, omega, eta) &
       result(cmu)
    ! parameters
    real(real64), intent(in) :: mu, omega, eta

    ! locals
    complex(real64) :: omega_m

    omega_m = mu / eta

    cmu = mu / (1. + z*omega_m/omega)

  end function cmu_maxwell

  complex(real64) function cmu_SLS(mu0, omega, eta, mu_factor) &
       result(cmu)
    ! arguments
    real(real64), intent(in) :: mu0, omega, eta, mu_factor

    ! locals
    real(real64) :: mu1, dmu, tau

    mu1 = mu0 * mu_factor

    dmu = mu0 * (1. - mu1 / (mu0 + mu1))

    tau = eta / (mu0 + mu1)

    cmu = mu0 - dmu / (1. + z*tau*omega)

  end function cmu_SLS

  complex(real64) function cmu_andrade(mu, omega, eta) &
       result(cmu)
    ! parameters
    real(real64), parameter :: alpha = 0.3

    ! arguments
    real(real64), intent(in) :: mu, omega, eta

    ! locals
    complex(real64) :: J
    real(real64) :: beta

    beta = mu**(alpha - 1.) * eta**(-alpha)

    J = 1 / mu - z / (eta * omega) + beta*(z*omega)**(-alpha) * gamma(1. + alpha)

    cmu = 1 / J

  end function cmu_andrade

  function planet_mmotion(m, a) &
       result(n)
    !-------------------------------------------------------------------------
    ! planet_mmotion : mean motion
    !
    ! Calculates the mean motion of a planet.
    !
    ! Arguments:
    !   m - real - mass
    !   a - real - semi-major axis of orbit
    ! Result:
    !   n - real - mean motion
    !-------------------------------------------------------------------------
    ! arguments
    real(real64) :: m, a
    real(real64) :: n

    n = sqrt(G * m / a**3)

  end function planet_mmotion

  function planet_thermo(len, pathdat) &
       result(thermo)
    !-------------------------------------------------------------------------
    ! planet_thermo : collects thermo values
    !
    ! Collects thermo values from fort.56 to be used as a planet formation.
    ! Values produced are used with a planet_form call.
    !
    ! Arguments:
    !   len     - integer   - size of fort.56 file
    !   pathdat - character - path to file
    ! Result:
    !   thermo - real(len,3) - returns columns of interest where the first
    !                          through third are pressure, density, and shear
    !                          modulus, respectively.
    !-------------------------------------------------------------------------
    ! arguments
    integer(int64), intent(in) :: len
    character(*), intent(in) :: pathdat
    real(real64) :: thermo(len,3)

    ! locals
    real(real64) :: val(5)
    integer(int64) :: i

    open(12, file=pathdat)

    do i = 1, len

       read(12, *) val(1), val(2), val(3), val(4), val(5)

       thermo(i,1) = val(1)*1e9
       thermo(i,2) = val(4)*1e5
       thermo(i,3) = val(5)*1e3

    end do
    close(12, status='keep')

  end function planet_thermo

end module mod_structure
