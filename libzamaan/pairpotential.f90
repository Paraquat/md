module pairpotential

use datatypes
use vector
use constants

implicit none

private
public :: type_pp

type type_pp
  integer       :: ns
  real(double), allocatable, dimension(:,:)  :: r_cut
  character(2), allocatable, dimension(:)    :: species_label
  real(double), allocatable, dimension(:,:)  :: sigma
  real(double), allocatable, dimension(:,:)  :: epsilon
  real(double), allocatable, dimension(:,:)  :: f_shift, e_shift
  real(double), allocatable, dimension(:,:)  :: s6,s12,c_f,c_e
  contains
    procedure :: init_pp_from_file
    procedure :: init_pp
    procedure :: lj_energy
    procedure :: lj_force
end type type_pp

contains

  subroutine init_pp_from_file(pp, filename, shift)

    class(type_pp), intent(inout) :: pp

    character(40), intent(in)     :: filename
    logical, intent(in)           :: shift

    integer                       :: funit, ios, i, ns
    character(80)                 :: junk
    real(double), allocatable, dimension(:,:) :: sigma, epsilon, r_cut

    funit = 11

    open(funit, file=filename, iostat=ios)
    read(funit,*) ns
    read(funit,*) pp%species_label

    allocate(sigma(ns,ns), epsilon(ns,ns), r_cut(ns,ns))

    read(funit,*) junk
    do i=1,ns
      read(funit,*) r_cut(i,:)
    end do

    read(funit,*) junk
    do i=1,ns
      read(funit,*) epsilon(i,:)
    end do

    read(funit,*) junk
    do i=1,ns
      read(funit,*) sigma(i,:)
    end do

    close(funit)

    call pp%init_pp(ns, sigma, epsilon, r_cut, shift)

  end subroutine init_pp_from_file

  subroutine init_pp(pp, ns, sigma, epsilon, r_cut, shift)

    class(type_pp), intent(inout) :: pp

    integer, intent(in)           :: ns
    real(double), dimension(:,:), intent(in)  :: r_cut
    real(double), dimension(:,:), intent(in)  :: sigma
    real(double), dimension(:,:), intent(in)  :: epsilon
    logical, intent(in)                       :: shift

    integer                       :: i,j

    allocate(pp%sigma(ns,ns),pp%s6(ns,ns), pp%s12(ns,ns), pp%c_f(ns,ns), &
             pp%c_e(ns,ns), pp%f_shift(ns,ns), pp%e_shift(ns,ns), &
             pp%species_label(ns), pp%r_cut(ns,ns))
    pp%sigma = sigma
    pp%epsilon = epsilon
    pp%s6 = sigma**6
    pp%s12 = sigma**12
    pp%c_f = 24.0_double*epsilon
    pp%c_e = 4.0_double*epsilon
    pp%r_cut = r_cut

    pp%e_shift = zero
    pp%f_shift = zero

    ! Compute the force and energy shift
    if (shift .eqv. .true.) then
      do i=1,3
        do j=1,3
          pp%e_shift(i,j) = pp%lj_energy(pp%r_cut(i,j), i, j, .false.)
          pp%f_shift(i,j) = pp%lj_force(pp%r_cut(i,j), i, j, .false.)
        end do
      end do
    end if
  end subroutine init_pp

  function lj_energy(pp, r_ij, s_i, s_j, shift) result(e)

    ! passed variables
    class(type_pp), intent(in)  :: pp

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift

    real(double)                :: e
    ! local variables
    real(double)                :: inv_r6, inv_r12

    inv_r6 = one/(r_ij**6)
    inv_r12 = one/(r_ij**12)

    e = pp%c_e(s_i,s_j)*(pp%s12(s_i,s_j)*inv_r12-pp%s6(s_i,s_j)*inv_r6)
    if (shift .eqv. .true.) then
        e = e + pp%e_shift(s_i,s_j)
    end if

  end function lj_energy

  function lj_force(pp, r_ij, s_i, s_j, shift) result(f)

    ! passed variables
    class(type_pp), intent(in)  :: pp

    real(double), intent(in)    :: r_ij
    integer, intent(in)         :: s_i, s_j
    logical, intent(in)         :: shift

    real(double)                :: f
    ! local variables
    real(double)                :: inv_r7, inv_r13

    inv_r7 = one/(r_ij**7)
    inv_r13 = one/(r_ij**13)

    f = pp%c_f(s_i,s_j)*(2*pp%s12(s_i,s_j)*inv_r13-pp%s6(s_i,s_j)*inv_r7)
    if (shift .eqv. .true.) then
      f = f + pp%f_shift(s_i,s_j)
    end if

  end function lj_force

end module pairpotential
