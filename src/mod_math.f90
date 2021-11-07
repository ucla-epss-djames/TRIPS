!===============================================================================
! Author: David A. James, davidabraham@ucla.edu
!
! Background: The following is a module contains mathematical operations
!
! Module: mod_math
! - subroutine rk4: 4th-Order Runge-Kutta Method
! - function hunt_loc: locate index of given value
! - function poly_int: polynomial interpolation
! - function simpson: numerical integration
! - function linspace : equally spaced domain
!
!===============================================================================
module mod_math

  use iso_fortran_env

  use mod_consts

  implicit none

  ! procedure declaration for rk4 function pass
  ! DE System
  abstract interface
     pure function fnd(t, y, d) &
          result(dy)
       use iso_fortran_env
       implicit none
       real(real64), intent(in) :: t
       real(real64), intent(in) :: y(:)
       real(real64), intent(in) :: d(:,:)
       real(real64) :: dy(size(y))
     end function fnd

     ! procedure declaration for simpsons function pass
     ! single function
     real(real64) pure function dy(x, d) &
          result(y)
       use iso_fortran_env
       implicit none
       real(real64), intent(in) :: x
       real(real64), intent(in) :: d
     end function dy
  end interface

contains

  pure function rk4(df, y, xi, xf, n, data) &
       result(yy)
    !-----------------------------------------------------------------------------
    ! rk4 : 4th-Order Runge-Kutta Method
    !
    ! Solves out the ODE system with the IVs given using the RK4 method.
    ! NOTE: Procedure(fnd) is used as a pointer for the df function. df is defined
    ! s.t. f(t,y,d) where t is a scalar, y is a vector, and d is a 2D matrix.
    !
    ! Arguments:
    !   df   - function  - ODE system
    !   y    - real(:)   - IVs for system
    !   xi   - real      - initial boundary
    !   xf   - real      - final boundary
    !   n    - integer   - number of steps
    !   data - real(:,:) - user data for ODE system
    ! Result:
    !   yy - real(:,:) - output of solved system
    !-----------------------------------------------------------------------------
    ! arguments
    procedure(fnd), intent(in), pointer :: df
    real(real64), intent(in) :: y(:)
    real(real64), intent(in) :: xi, xf
    integer(int64), intent(in) :: n
    real(real64), intent(in) :: data(:,:)
    real(real64), allocatable  :: yy(:,:)

    ! locals
    real(real64), allocatable :: yi(:), k1(:), k2(:), k3(:), k4(:)
    real(real64) :: h, x
    integer(int64) :: i, len

    len = size(y)
    allocate(yi(len), k1(len), k2(len), k3(len), k4(len), yy(n+1,len+1))

    h = (xf - xi) / n
    x = xi

    yy(1,1) = x
    yy(1,2:len+1) = y

    yi = y

    do i = 2, n+1

       k1 = h * df(x, yi, data)
       k2 = h * df(x + h / 2, yi + k1 / 2, data)
       k3 = h * df(x + h / 2, yi + k2 / 2, data)
       k4 = h * df(x + h, yi + k3, data)

       yi = yi + (k1 + 2*k2 + 2*k3 + k4) / 6
       x = min(x + h, xf)

       yy(i,1) = x
       yy(i,2:len+1) = yi

    end do

  end function rk4

  integer(int64) pure function hunt_loc(data, mark) &
       result(i)
    !-----------------------------------------------------------------------------
    ! hunt_loc : find index of given value
    !
    ! This will hunt for the value in the given vector. If the value is not there,
    ! it will return the index of the closet value.
    !
    ! Arguments:
    !   data - real(:) - vector of data to search
    !   mark - real    - value of interest
    ! Result:
    !   i - integer - index of closest value in data
    !-----------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: data(:)
    real(real64), intent(in) :: mark

    ! locals
    real(real64) :: diff(size(data))

    diff = abs(data - mark)

    i = minloc(abs(data(:) - mark), dim=1, kind=4)

  end function hunt_loc

  real(real64) pure function poly_int(x0, x, f) &
       result(q0)
    !-----------------------------------------------------------------------------
    ! poly_int : polynomial interpolation
    !
    ! Uses Neville's method to interpolate a value from the function array
    !
    ! Arguments:
    !   x0 - real    - value of interest
    !   x  - real(:) - x-values from f(x)
    !   f  - real(:) - y-values from f(x)
    ! Result:
    !   q0 - real - interpolated value
    !-----------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: x0
    real(real64), intent(in) :: x(:), f(:)

    ! locals
    real(real64) :: num, dem
    real(real64), allocatable :: q(:,:)
    integer(int64) :: n, i, j

    n = size(x)
    allocate(q(n,n))
    q(:,1) = f

    do i = 2, n
       do j = i, n
          num = q(j,i-1)*(x0 - x(j-i+1)) - q(j-1,i-1)*(x0 - x(j))
          dem = x(j) - x(j-i+1)
          q(j,i) = num / dem
       end do
    end do

    q0 = q(n,n)

  end function poly_int

  real(real64) pure function simpson(df, xi, xf, data, n) &
       result(y)
    !-----------------------------------------------------------------------------
    ! simpson : numerical integration
    !
    ! Use Simpson's method to integrate a function over a bound (a,b). Function
    ! should be of the form f(x, d) where x is the independent variable and d is
    ! any data that needs to be passed to the function.
    !
    ! Arguments:
    !   df - function - function integrating
    !   xi - real     - start of bound
    !   xf - real     - end of bound
    ! data - real     - data for the function
    !    n - int      - step count (should be even)
    ! Result:
    !    y - real - integrand
    !-----------------------------------------------------------------------------
    ! arguments
    procedure(dy), intent(in), pointer :: df
    real(real64), intent(in) :: xi, xf, data
    integer(int64), intent(in) :: n

    ! locals
    real(real64) :: h, x, xi0, xi1, xi2
    integer(int64) :: i

    h = (xf - xi) / n

    xi0 = df(xi, data) + df(xf, data)
    xi1 = 0
    xi2 = 0

    do i = 1, n - 1

       x = xi + i*h

       if(mod(i,2) == 0) then
          xi2 = xi2 + df(x, data)
       else
          xi1 = xi1 + df(x, data)
       end if

    end do

    y = h*(xi0 + 2*xi2 + 4*xi1) / 3

  end function simpson

  pure function linspace(a, b, n) &
       result(arr)
    !-----------------------------------------------------------------------------
    ! linspace : equally spaced domain
    !
    ! Generate a domain over a bound (a,b) s.t. there are n spaces equally spaced
    ! out. The bound should be defined s.t. a < b.
    !
    ! Arguments:
    !   a - real - start of the bound
    !   b - real - end of the bound
    !   n - real - number of spaces
    ! Result:
    !   arr - real(n) - domain
    !-----------------------------------------------------------------------------
    ! arguments
    real(real64), intent(in) :: a, b
    integer(int64), intent(in) :: n
    real(real64) :: arr(n)

    ! locals
    real(real64) :: range
    integer(int64) :: i

    range = b - a

    if(n == 0) return

    if(n == 1) then
       arr(1) = a
       return
    end if

    do i = 1, n
       arr(i) = a + range*(i - 1) / (n - 1)
    end do

  end function linspace

  pure function logspace(a, b, n) &
       result(arr)
    ! arguments
    real(real64), intent(in) :: a, b
    integer(int64), intent(in) :: n
    real(real64) :: arr(n)

    ! locals
    real(real64) :: logarr(n)
    integer(int64) :: i

    logarr = linspace(log10(a), log10(b), n)

    do i = 1, n

       arr(i) = 10**logarr(i)

    end do

  end function logspace


end module mod_math
