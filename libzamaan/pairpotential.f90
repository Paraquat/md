module pairpotential

use datatypes
use vector
use constants

implicit none

private
public :: type_pp

type type_pp
  integer       :: nspec
  real(double), allocatable, dimension(:,:)  :: sigma
  real(double), allocatable, dimension(:,:)  :: epsilon
  real(double), allocatable, dimension(:,:)  :: s6,s12,c_f,c_e
  contains
    procedure :: init_pp
    procedure :: lj_energy
    procedure :: lj_force
end type type_pp

contains

  subroutine init_pp(pp, nspec, sigma, epsilon)

    class(type_pp), intent(inout) :: pp

    integer, intent(in)           :: nspec
    real(double), dimension(:,:), intent(in)  :: sigma
    real(double), dimension(:,:), intent(in)  :: epsilon

    allocate(pp%sigma(nspec,nspec),pp%s6(nspec,nspec), &
             pp%s12(nspec,nspec), pp%c_f(nspec,nspec), &
             pp%c_e(nspec,nspec))
    pp%sigma = sigma
    pp%epsilon = epsilon
    pp%s6 = sigma**6
    pp%s12 = sigma**12
    pp%c_f = 24.0_double*epsilon
    pp%c_e = 4.0_double*epsilon

  end subroutine init_pp

  function lj_energy(pp, r_ij, s_i, s_j) result(e)

    class(type_pp), intent(in)  :: pp

    ! passed variables
    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j

    real(double)                :: e
    ! local variables
    real(double)                :: inv_r6, inv_r12

    inv_r6 = one/(r_ij**6)
    inv_r12 = one/(r_ij**12)

    e = pp%c_e(s_i,s_j)*(pp%s12(s_i,s_j)*inv_r12-pp%s6(s_i,s_j)*inv_r6)

  end function lj_energy

  function lj_force(pp, r_ij, s_i, s_j) result(f)

    class(type_pp), intent(in)  :: pp

    ! passed variables
    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j

    real(double)                :: f
    ! local variables
    real(double)                :: inv_r7, inv_r13

    inv_r7 = one/(r_ij**7)
    inv_r13 = one/(r_ij**13)

    f = pp%c_f(s_i,s_j)*(2*pp%s12(s_i,s_j)*inv_r13-pp%s6(s_i,s_j)*inv_r7)

  end function lj_force

end module pairpotential
