module constants

  use datatypes

  implicit none
  save

  real(double), parameter     :: zero = 0.0_double
  real(double), parameter     :: one = 1.0_double
  real(double), parameter     :: two = 2.0_double
  real(double), parameter     :: three = 3.0_double
  real(double), parameter     :: four = 4.0_double
  real(double), parameter     :: eight = 8.0_double
  real(double), parameter     :: half = one/two
  real(double), parameter     :: third = one/three
  real(double), parameter     :: quarter = one/four
  real(double), parameter     :: eighth = one/eight
  real(double), parameter     :: pi = 3.14159265359_double
  real(double), parameter     :: twopi = two*pi
  real(double), parameter     :: rad2deg = 360.0_double/twopi
  real(double), parameter     :: deg2rad = one/rad2deg
  real(double), parameter     :: small = 1.0E-8_double
  real(double), parameter     :: k_B = 1.3806485279E-23_double
  real(double), parameter     :: q_e = 1.602176620898E-19_double

end module constants
