program check_identity

  use iso_fortran_env

  use mod_tidal

  implicit none

  ! locals
  real(real64) :: l, randre, randim
  complex(real64) :: vals(4)
  complex(real64), dimension(M,M) :: Y, YINV, RES
  integer(int64) :: i, j

  l = 2.

  do i = 1, 4
     call random_number(randre)
     call random_number(randim)

     vals(i) = cmplx(randre, randim, real64)
  end do

  Y = propagator_matrix(l, vals(1), vals(2), vals(3), vals(4))
  YINV = inverse_matrix(l, vals(1), vals(2), vals(3), vals(4))

  print *, "Displaying Matrix Re and Im pairs"

  print *, "Y"
  do i = 1, M
     print '(99f12.5)', (Y(i,j), j = 1, M)
  end do

  print *, "YINV"
  do i = 1, M
     print '(99f12.5)', (YINV(i,j), j = 1, M)
  end do

  print *, "Identity Check: Y*YINV"
  RES = matmul(Y, YINV)
  do i = 1, M
     print '(99f12.5)', (RES(i,j), j = 1, M)
  end do

  print *, "Identity Check: YINV*Y"
  RES = matmul(YINV, Y)
  do i = 1, M
     print '(99f12.5)', (RES(i,j), j = 1, M)
  end do

end program check_identity
