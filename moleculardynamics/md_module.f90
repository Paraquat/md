module md_module

use datatypes
use constants
use vector
use cell
use pairpotential

implicit none

type type_md
  type(type_cell) :: p_t, p_t_dt
  type(type_pp)   :: pp
  integer         :: nstep, nspec, nat
  logical         :: shift
  integer, allocatable, dimension(:)          :: species
  real(double), allocatable, dimension(:,:)   :: v_t, v_t_dt
  real(double), allocatable, dimension(:,:)   :: f_t, f_t_dt
  real(double), dimension(3,3)                :: p_g_t, p_g_t_dt
  real(double), dimension(3,3)                :: P_int_t, P_int_t_dt, P_ext
  real(double)   :: pe_t, pe_t_dt, ke_t, ke_t_dt, T_int_t, T_int_t_dt, T_ext

  contains
    procedure :: init_md
    procedure :: get_force_and_energy
end type type_md

contains

  subroutine init_md(mdr, init_cell, pp, nstep, shift)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    type(type_cell), intent(in)     :: init_cell
    type(type_pp), intent(in)       :: pp
    logical, intent(in)             :: shift
    integer, intent(in)             :: nstep

    mdr%p_t = init_cell
    mdr%nspec = mdr%p_t%nspec
    call mdr%p_t_dt%init(mdr%p_t%nat, mdr%p_t%nspec, mdr%p_t%h, &
                         mdr%p_t%r, mdr%p_t%species)
    mdr%pp = pp
    mdr%nat = init_cell%nat
    mdr%nstep = nstep
    mdr%shift = shift

    allocate(mdr%species(mdr%nat))
    allocate(mdr%v_t(mdr%nat,3), mdr%v_t_dt(mdr%nat,3))
    allocate(mdr%f_t(mdr%nat,3), mdr%f_t_dt(mdr%nat,3))
    mdr%species = init_cell%spec_int
  end subroutine init_md

  subroutine get_force_and_energy(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer                         :: iat, jat, s_i, s_j
    real(double), dimension(3)      :: r_ij, r_ij_cart
    real(double)                    :: mod_r_ij, mod_f, pe

    mdr%f_t = zero
    mdr%pe_t = zero

    do iat=1,mdr%nat
      do jat=1,mdr%nat
        if (iat == jat) cycle
        s_i = mdr%species(iat)
        s_j = mdr%species(jat)
        r_ij = mdr%p_t%r(jat,:) - mdr%p_t%r(iat,:)
        r_ij_cart = mdr%p_t%disp_frac2cart_noshift(r_ij)
        mod_r_ij = modulus(r_ij_cart)
        if (mod_r_ij < mdr%pp%r_cut(s_i,s_j)) then
          r_ij_cart = norm(r_ij_cart)
          mod_f = mdr%pp%lj_force(mod_r_ij, s_i, s_j, mdr%shift)
          mdr%f_t(iat,:) = mdr%f_t(iat,:) + mod_f*r_ij_cart
          pe = mdr%pp%lj_energy(mod_r_ij, s_i, s_j, mdr%shift)
          mdr%pe_t = mdr%pe_t + pe
        end if
      end do
    end do
    pe = pe*half
  end subroutine get_force_and_energy

end module md_module
