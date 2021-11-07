!==============================================================================
! Author: David A. James, davidabraham@ucla.edu
!
! Background: The following is the tidal love number calculations. Refer to
! Henning & Hurford '14 (HH14) for a full description of the method. Core
! matrix comes from either HH14 or "Global Dynamics" (GD) by Sabadini
!
! Module: mod_tidal
! - function propagator_method    : solve out tidal numbers
! - subroutine aggregate_matrix   : solve out aggregate matrix
! - function core_sabadini_matrix : initial matrix at the core
! - function core_HH_matrix       : initial matrix at the core
! - function propagator_matrix    : matrix to propogate up to the surface
! - function inverse_matrix       : inverse of the propagator matrix
! - function solve_vector         : solves out bound vector
! - function solve_layer          : solves out tidal numbers at given layer
! - subroutine use_GD_core        : determines core matrix
! - subroutine normalize          : calculates normalization units
!
!==============================================================================
module mod_tidal

  use iso_fortran_env

  use mod_consts, only : PI, G

  implicit none

  integer(int64), parameter :: M = 6, N = 3    ! CONTAINER SIZES
  integer(int64), parameter :: I_ONE = 1       ! LAPACK ARG

  interface core_matrix
     module procedure core_sabadini_matrix, core_HH_matrix
  end interface core_matrix

contains

  function propagator_method(l, layers, data, flag) &
       result(tidal)
    !--------------------------------------------------------------------------
    ! propagator_method: solve out tidal numbers
    !
    ! Performs the propogator method from HH'14 to produce the tidal love
    ! numbers for the planetary body.
    !
    ! Arguments:
    !   l      - real         - harmonic degree
    !   layers - integer      - number of layers
    !   data   - complex(:,:) - data for method
    !   flag   - logical      - core choice
    ! Result:
    !   tidal - complex(layers, 3) - tidal love numbers for each layer of the
    !                                body arranged s.t. (kl, hl, ll)
    !--------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: l
    integer(int64), intent(in) :: layers
    complex(real64), intent(in) :: data(:,:)
    logical, intent(in) :: flag
    complex(real64) :: tidal(layers, 3)

    ! locals
    complex(real64), dimension(M, N) :: B
    complex(real64), dimension(M, N, layers) :: Bi
    complex(real64), dimension(N) :: bound

    call aggregate_matrix(l, layers, data, flag, B, Bi)

    bound = solve_vector(l, data(layers,1), B)

    tidal = solve_layer(l, layers, data, bound, Bi)

  end function propagator_method

  subroutine aggregate_matrix(l, layers, data, flag, B, Bi)
    !-------------------------------------------------------------------------
    ! aggregate_matrix : solve out aggregate matrix
    !
    ! Solves out the aggregate matrix s.t. B is the surface aggregate and Bi
    ! is the aggregate for each layer of the planetary body.
    ! Refer to A5 in HH14.
    !
    ! Arguments:
    !   l      - real         - harmonic degree
    !   layers - integer      - number of layers
    !   data   - complex(:,:) - data for method
    !   flag   - logical      - core choice
    ! Result:
    !   B  - complex(M,N)        - surface aggregate matrix
    !   Bi - complex(M,N,layers) - layer aggregate matrix
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: l
    integer(int64), intent(in) :: layers
    complex(real64), intent(in) :: data(:,:)
    logical, intent(in) :: flag
    complex(real64), intent(out) :: B(:,:), Bi(:,:,:)

    ! locals
    complex(real64), allocatable :: r(:), mu(:), a(:), rho(:)
    complex(real64), dimension(M, M) :: Yi, Yi1, Y
    complex(real64), dimension(M, N) :: Bi1
    integer(int64) :: i, j, k

    r = data(:,1)
    mu = data(:,2)
    a = data(:,3)
    rho = data(:,4)

    if(flag) then
       ! GD call
       B = core_matrix(l, r(1), a(1), rho(1))
    else
       ! HH14 call
       B = core_matrix(l, r(1), mu(1), a(1), rho(1))
    end if

    Bi(:,:,1) = B

    do i = 2, layers

       Bi1 = B
       Yi1 = inverse_matrix(l, r(i-1), mu(i), a(i-1), rho(i))
       Yi = propagator_matrix(l, r(i), mu(i), a(i), rho(i))

       Y = matmul(Yi, Yi1)
       B = matmul(Y, Bi1)

       Bi(:,:,i) = B

    end do

  end subroutine aggregate_matrix

  function core_sabadini_matrix(l, r, a, rho) &
       result(B)
    !-------------------------------------------------------------------------
    ! core_sabadini_matrix: initial matrix at the core
    !
    ! Constructs the initial core matrix using eq 2.6 from "Global Dynamics"
    ! by Sabadini
    !
    ! Arguments:
    !   l   - real    - harmonic degree
    !   r   - complex - radius
    !   a   - complex - acceleration due to gravity
    !   rho - complex - density
    ! Result:
    !   B - complex(:,:) - core matrix
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: l
    complex(real64), intent(in) :: r, a, rho
    complex(real64)  :: B(M,N)

    ! locals
    real(real64) :: Ac

    Ac = a / r

    B(:,:) = 0.

    B(1,1) = -r**l / a
    B(1,3) = 1.

    B(2,2) = 1.

    B(3,3) = rho * r * Ac

    B(5,1) = r**l

    B(6,1) = 2.*(l - 1.)*r**(l - 1.)
    B(6,3) = 3. * Ac

  end function core_sabadini_matrix

  function core_HH_matrix(l, r, mu, a, rho) &
       result(B)
    !-------------------------------------------------------------------------
    ! core_HH_matrix: initial matrix at the core
    !
    ! Constructs the initial core matrix. Refer to HH'14 A3's first three
    ! columns
    !
    ! Arguments:
    !   l   - real    - harmonic degree
    !   r   - complex - radius
    !   mu  - complex - complex shear modulus
    !   a   - complex - acceleration due to gravity
    !   rho - complex - density
    ! Result:
    !   B - complex(:,:) - core matrix
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: l
    complex(real64), intent(in) :: r, mu, a, rho
    complex(real64)  :: B(M,N)

    ! locals
    real(real64) :: lp1, lp2, lp3, lm1, lm2
    integer(int64) :: i

    lp1 = l + 1
    lp2 = l + 2
    lp3 = l + 3
    lm1 = l - 1
    lm2 = l - 2

    B(:,:) = 0

    B(1,1) = l * r**lp1 / (2.*(2.*l + 3.))
    B(2,1) = lp3 * r**lp1 / (2.*(2.*l + 3.)*lp1)
    B(3,1) = (l*rho*a*r + 2.*mu*(l**2. - lp3))*r**l / (2.*(2.*l + 3.))
    B(4,1) = l*lp2*mu * r**l / ((2.*l + 3.)*lp1)
    B(6,1) = 2.*PI * rho*l * r**lp1 / (2.*l + 3.)

    B(1,2) = r**lm1
    B(2,2) = r**lm1 / l
    B(3,2) = (r*rho*a + 2.*mu*lm1)*r**lm2
    B(4,2) = 2.*lm1*mu * r**lm2 / l
    B(6,2) = 4.*PI * rho * r**lm1

    B(3,3) = -rho * r**l
    B(5,3) = -r**l
    B(6,3) = -(2.*l + 1.) * r**lm1


  end function core_HH_matrix

  function propagator_matrix(l, r, mu, a, rho) &
       result(Y)
    !-------------------------------------------------------------------------
    ! propagator_matrix: matrix to propogate up to the surface
    !
    ! Constructs the propagator matrix.
    !
    ! Arguments:
    !   l   - real    - harmonic degree
    !   r   - complex - radius
    !   mu  - complex - shear modulus
    !   a   - complex - acceleration due to gravity
    !   rho - complex - density
    ! Result:
    !   Y - complex(:,:) - propagator matrix
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: l
    complex(real64), intent(in) :: r, mu, a, rho
    complex(real64) :: Y(M,M)

    ! locals
    real(real64) :: lp1, lp2, lp3, lm1, lm2

    lp1 = l + 1
    lp2 = l + 2
    lp3 = l + 3
    lm1 = l - 1
    lm2 = l - 2

    Y(:,:) = 0

    Y(1:M,1:N) = core_matrix(l, r, mu, a, rho)

    Y(1,4) = lp1 * r**(-l) / (2.*(2.*l - 1.))
    Y(2,4) = (2. - l)*r**(-l) / (2.*l*(2.*l - 1.))
    Y(3,4) = (r*rho*a*lp1 - 2.*mu*(l**2. + 3.*l - 1.)) / (2.*(2.*l - 1.) * r**lp1)
    Y(4,4) = (l**2. - 1.)*mu / (l*(2.*l - 1.)*r**lp1)
    Y(6,4) = 2.*PI * rho*lp1 / ((2.*l - 1.)*r**l)

    Y(1,5) = r**(-lp2)
    Y(2,5) = -r**(-lp2) / lp1
    Y(3,5) = (r*rho*a - 2.*mu*lp2) / r**lp3
    Y(4,5) = (2.*lp2 * mu) / (lp1 * r**lp3)
    Y(6,5) = 4.*PI * rho / r**lp2

    Y(3,6) = -rho / r**lp1
    Y(5,6) = -1. / r**lp1

  end function propagator_matrix

  function inverse_matrix(l, r, mu, a, rho) &
       result(Y)
    !-------------------------------------------------------------------------
    ! inverse_matrix: inverse of the propagator matrix
    !
    ! Constructs the inverse matrix.
    !
    ! Arguments:
    !   l   - real    - harmonic degree
    !   r   - complex - radius
    !   mu  - complex - shear modulus
    !   a   - complex - acceleration due to gravity
    !   rho - complex - density
    ! Result:
    !   Y - complex(:,:) - propagator matrix
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: l
    complex(real64), intent(in) :: r, mu, a, rho
    complex(real64) :: Y(M,M)

    ! locals
    complex(real64), dimension(M,M) :: B, D
    real(real64) :: lp1, lp2, lp3

    lp1 = l + 1
    lp2 = l + 2
    lp3 = l + 3

    B(:,:) = 0

    B(1,1) = rho*a*r / mu - 2.*lp2
    B(2,1) = -rho*a*r / mu + 2.*(l**2. + 3.*l - 1.) / lp1
    B(3,1) = 4.*PI * rho
    B(4,1) = rho*a*r / mu + 2.*(l - 1.)
    B(5,1) = -rho*a*r / mu - 2.*(l**2. - lp3) / l
    B(6,1) = 4.*PI * rho*r

    B(1,2) = 2.*l * lp2
    B(2,2) = -2.*(l**2. - 1.)
    B(4,2) = 2.*(l**2. - 1.)
    B(5,2) = -2.*l * lp2

    B(1,3) = -r / mu
    B(2,3) = r / mu
    B(4,3) = -r / mu
    B(5,3) = r / mu

    B(1,4) = l*r / mu
    B(2,4) = (2. - l)*r / mu
    B(4,4) = -lp1*r / mu
    B(5,4) = lp3*r / mu

    B(1,5) = rho*r / mu
    B(2,5) = -rho*r / mu
    B(4,5) = rho*r / mu
    B(5,5) = -rho*r / mu
    B(6,5) = 2.*l + 1

    B(3,6) = -1
    B(6,6) = -r

    D(:,:) = 0

    D(1,1) = lp1 * r**-lp1
    D(2,2) = l*lp1 * r**(-l + 1.) / (2.*(2.*l - 1.))
    D(3,3) = r**(-l + 1.)
    D(4,4) = l * r**l
    D(5,5) = l*lp1 * r**lp2 / (2.*(2.*l + 3.))
    D(6,6) = -r**lp1

    D(:,:) = D(:,:) / (2.*l + 1.)

    Y = matmul(D,B)

  end function inverse_matrix

  function solve_vector(l, r, B) &
       result(bound)
    !-------------------------------------------------------------------------
    ! solve_vector : solves out bound vector
    !
    ! Solves out the bound vector for the planetary body. Refer to A6-A7 of
    ! HH'14.
    ! NOTE: this function does not have `pure` because the calls to lapack.
    !
    ! Arguments:
    !   l - real         - harmonic degree
    !   r - complex      - radius
    !   B - complex(:,:) - surface aggregate propagator
    ! Result:
    !   bound - complex(:) - solved vector c from A7
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: l
    complex(real64), intent(in) :: r, B(:,:)
    complex(real64) :: bound(N)
    external :: ZGETRF, ZGETRS

    ! locals
    complex(real64) :: Mat(N,N)
    integer(int64) :: ipiv(N), info

    bound(:) = 0
    bound(3) = -(2.*l + 1.) / r

    Mat(1,:) = B(3,:)
    Mat(2,:) = B(4,:)
    Mat(3,:) = B(6,:)

    call ZGETRF(N, N, Mat, N, ipiv, info)
    call ZGETRS('N', N, I_ONE, Mat, N, ipiv, bound, N, info)

  end function solve_vector

  function solve_layer(l, layers, data, c, Bi) &
       result(tidal)
    !-------------------------------------------------------------------------
    ! solve_layer: solves out tidal numbers for body
    !
    ! Solves out tidal values. Refer to HH'14 A8-A11
    !
    ! Arguments:
    !   l      - real                - harmonic degree
    !   layers - integer             - number of layers
    !   data   - complex(:,:)        - data for method
    !   c      - complex(M)          - solved out bound vector
    !   Bi     - complex(M,N,layers) - layer aggregate matrix
    ! Result:
    !   tidal - complex(layers, 3) - tidal love numbers for each layer of the
    !                                body arranged s.t. (kl, hl, ll)
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: l
    integer(int64), intent(in) :: layers
    complex(real64), intent(in) :: data(:,:), c(:), Bi(:,:,:)
    complex(real64) :: tidal(layers, 3)

    ! locals
    complex(real64), allocatable :: r(:), a(:)
    complex(real64) :: kl, hl, ll, yvec(M)
    integer(int64) :: i

    r = data(:,1)
    a = data(:,3)

    do i = 1, layers

       yvec = matmul(Bi(:,:,i), c)

       kl = -r(i)**l - yvec(5)
       hl = a(i) * yvec(1)
       ll = a(i) * yvec(2)

       tidal(i,:) = [kl, hl, ll]

    end do

  end function solve_layer

  subroutine use_GD_core(r, mu, flag)
    !-------------------------------------------------------------------------
    ! use_GD_core : determines core matrix
    !
    ! If structure is fluid then GD core will be used, otherwise the HH core
    ! will be used for the propagation.
    !
    ! Arguments:
    !   r  - real - first radius value of structure data
    !   mu - real - first shear value of structure data
    ! Results:
    !   flag - logical - true will use the GD core
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: r, mu
    logical, intent(out) :: flag

    if(r > 0.0 .and. mu == 0.0) then
       flag = .true.
    else
       flag = .false.
    end if

  end subroutine use_GD_core

  subroutine normalize(mass, r, data)
    !-------------------------------------------------------------------------
    ! normalize : calculates normalization units
    !
    ! Calculates the normalization units for the data
    !
    ! Arguments:
    !   mass - real         - calculated mass of body
    !   data - complex(:,:) - data to normalize
    ! Result:
    !   data - complex(:,:) - data normalized
    !-------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: mass, r
    complex(real64), intent(in out) :: data(:,:)

    ! locals
    complex(real64) :: a, rho, p

    a = mass / r**2 * G
    rho = mass / r**3
    p = rho * a * r

    data(:,1) = data(:,1) / r
    data(:,2) = data(:,2) / p
    data(:,3) = data(:,3) / a
    data(:,4) = data(:,4) / rho

  end subroutine normalize

end module mod_tidal
