module md_module

use datatypes
use constants
use vector
use cell
use pairpotential
use rng

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
  real(double)   :: pe_t, pe_t_dt, ke_t, ke_t_dt, T_int_t, T_int_t_dt, T_ext, dt

  contains
    procedure :: init_md
    procedure :: get_force_and_energy
    procedure :: get_kinetic_energy
    procedure :: vVerlet_v_half
    procedure :: vVerlet_r
    procedure :: md_run
end type type_md

contains

  ! Initialise variables/velocities/allocate matrices for MD run
  subroutine init_md(mdr, init_cell, pp, nstep, dt, T_ext, shift)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    type(type_cell), intent(in)     :: init_cell
    type(type_pp), intent(in)       :: pp
    integer, intent(in)             :: nstep
    real(double)                    :: dt
    logical, intent(in)             :: shift

    ! local variables
    integer                         :: i, j
    real(double), dimension(3)      :: sumv
    real(double)                    :: sumv2, sfac

    mdr%p_t = init_cell
    mdr%nspec = mdr%p_t%nspec
    call mdr%p_t_dt%init(mdr%p_t%nat, mdr%p_t%nspec, mdr%p_t%h, &
                         mdr%p_t%r, mdr%p_t%species)
    mdr%pp = pp
    mdr%nat = init_cell%nat
    mdr%nstep = nstep
    mdr%dt = dt
    mdr%T_ext = T_ext
    mdr%shift = shift

    allocate(mdr%species(mdr%nat))
    allocate(mdr%v_t(mdr%nat,3), mdr%v_t_dt(mdr%nat,3))
    allocate(mdr%f_t(mdr%nat,3), mdr%f_t_dt(mdr%nat,3))
    mdr%species = init_cell%spec_int

    ! initialise rng
    call init_rand

    ! initialise velocities (uniform random distribution), scale according
    ! to T_ext
    sumv2 = zero
    do i=1,mdr%nat
      do j=1,3
        call rand(mdr%v_t(i,j))
        mdr%v_t(i,j) = mdr%v_t(i,j) - half
      end do
      sumv2 = sumv2 + sum(mdr%v_t(i)**2)
    end do

    do i=1,3
      sumv(i) = sum(mdr%v_t(:,i))
    end do
    sumv = sumv/mdr%nat
    sfac = sqrt(three*(mdr%nat-1)*mdr%T_ext/sumv2) ! nat-1 because COM is fixed
    do i=1,mdr%nat
      mdr%v_t(i,:) = sfac*(mdr%v_t(i,:) - sumv) ! remove COM velocity
    end do

  end subroutine init_md

  ! Compute the potential energy and force on each atom
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

  subroutine get_kinetic_energy(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer   :: i

    mdr%ke_t = 0

    do i=1,mdr%nat
      mdr%ke_t = mdr%ke_t + sum((mdr%v_t(i,:))**2)
    end do
    mdr%ke_t = mdr%ke_t/two

  end subroutine get_kinetic_energy

  subroutine vVerlet_v_half(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    mdr%v_t_dt = mdr%v_t + mdr%dt*half*(mdr%f_t + mdr%f_t_dt)

  end subroutine vVerlet_v_half

  subroutine vVerlet_r(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    mdr%p_t_dt%r = mdr%p_t%r + mdr%dt*mdr%v_t_dt

  end subroutine vVerlet_r

  subroutine md_run(mdr, s_start, s_end)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    integer, intent(in)             :: s_start, s_end

    ! local variables
    integer   :: s
    subroutine init_md(mdr, init_cell, pp, nstep, dt, T_ext, shift)


    call mdr%get_force_and_energy
    do s=s_start,s_end
      call mdr%vVerlet_v_half
      call mdr%vVerlet_r
      call mdr%vVerlet_v_half
      call mdr%get_kinetic_energy
      call mdr%get_force_and_energy
      mdr%v_t = mdr%v_t_dt
      mdr%p_t%r = mdr%p_t_dt%r
    end do

  end subroutine md_run

end module md_module
