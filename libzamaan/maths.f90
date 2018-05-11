module maths

use datatypes
use constants

implicit none

real(double), parameter  :: six = 6.0_double
real(double), parameter  :: twenty = 20.0_double
real(double), parameter  :: fortytwo = 42.0_double
real(double), parameter  :: seventytwo = 72.0_double
real(double), parameter  :: sh2 = one/six
real(double), parameter  :: sh4 = sh2/twenty
real(double), parameter  :: sh6 = sh4/fortytwo
real(double), parameter  :: sh8 = sh6/seventytwo

contains

! Taylor expansion of sinh(x)/x
function poly_sinhx_x(x) result(f)
  real(double)  :: x, f
  x = x**2
  f = (((sh8*x + sh6)*x + sh4)*x + sh2)*x + one
end function poly_sinhx_x

! 3x3 matrix multiplication using DGEMM
function mmult3x3(a, b) result (c)
  real(double), dimension(3,3)  :: a, b, c

  call dgemm('N', 'N', 3, 3, 3, one, a, 3, b, 3, one, c, 3)
end function mmult3x3

! 3x3 matrix multiplication using do loop
function mmult3x3_loop(a, b) result (c)
  real(double), dimension(3,3)  :: a, b, c
  integer                       :: i, j, k

  c = zero
  do i=1,3
    do j=1,3
      do k=1,3
        c(j,i) = c(j,i) + a(k,j)*b(k,i)
      end do
    end do
  end do
end function mmult3x3_loop

end module maths
