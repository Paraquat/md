module constants

  use datatypes

  implicit none
  save

  real(double), parameter     :: zero = 0.0_double
  real(double), parameter     :: one = 1.0_double
  real(double), parameter     :: two = 2.0_double
  real(double), parameter     :: three = 3.0_double
  real(double), parameter     :: half = 0.5_double
  real(double), parameter     :: pi = 3.14159265359_double
  real(double), parameter     :: rad2deg = 57.29577951_double
  real(double), parameter     :: deg2rad = one/rad2deg
  real(double), parameter     :: small = 1.0E-8_double
  real(double), parameter     :: k_B = 1.3806485279E-23_double

end module constants
