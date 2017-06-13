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
function exp_sinhx_x(x) result(f)
  real(double)  :: x, f
  x = x**2
  f = (((sh8*x + sh6)*x + sh4)*x + sh2)*x + one
end function exp_sinhx_x

end module maths
