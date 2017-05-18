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

  subroutine init_md(mdr, init_cell, pp, ns, shift)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    type(type_cell), intent(in)     :: init_cell
    type(type_pp), intent(in)       :: pp
    logical, intent(in)             :: shift

    integer, intent(in)             :: ns

    mdr%p_t = init_cell
    call mdr%p_t_dt%init(mdr%p_t, mdr%p_t%nat, mdr%p_t%nspec, mdr%p_t%h, &
                         mdr%p_t%r, mdr%p_t%species)
    mdr%pp = pp
    mdr%h_t = init_cell%h
    mdr%nat = init_cell%nat
    mdr%nstep = ns
    mdr%shift = shift

    allocate(mdr%species(nat))
    allocate(mdr%r_t(mdr%nat,3), mdr%r_t_dt(mdr%nat,3))
    allocate(mdr%v_t(mdr%nat,3), mdr%v_t_dt(mdr%nat,3))
    allocate(mdr%f_t(mdr%nat,3), mdr%f_t_dt(mdr%nat,3))
  end subroutine init_md

  subroutine get_force_and_energy(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer                         :: iat, jat
    real(double), dimension(3)      :: r_ij, r_ij_cart
    real(double)                    :: mod_r_ij, mod_f, pe

    mdr%f_t = zero

    do iat=1,mdr%nat
      do jat=1,mdr%nat
        if (iat == jat) cycle
        r_ij = mdr%p_dt(jat,:) - mdr%p_dt(iat,:)
        r_ij_cart = mdr%p_dt%disp_frac2cart_noshift(r_ij)
        mod_r_ij = modulus(r_ij_cart)
        r_ij_cart = norm(r_ij_cart)
        mod_f = mdr%pp%lj_force(mod_r_ij, mdr%species(iat), mdr%species(jat), &
                                shift)
        mdr%f_t(iat) = mdr%f_t(iat) + mod_f*r_ij_cart
        pe = mdr%pp%lf_energy(mod_r_ij, mdr%species(iat), mdr%species(jat), &
                              shift)
        mdr%pe_t = mdr%pe_t +
      end do
    end do
  end subroutine get_force_and_energy

end module md_module
