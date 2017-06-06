module barostat_module

use datatypes
use constants
use rng

implicit none

type type_barostat

  real(double)                  :: W_g    ! box mass
  real(double), dimension(3,3)  :: h      ! box matrix (a, b, c)
  real(double), dimension(3,3)  :: sigma  ! inverse lattice matrix
  real(double), dimension(3,3)  :: p_g    ! box momentum
  real(double), dimension(3,3)  :: stress ! stress tensor

contains

  procedure :: weak_coupling_p

end type type_barostat

contains

  subroutine weak_coupling_p
  end subroutine weak_coupling_p

end module barostat_module
