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
  character(3)    :: ensemble
  integer         :: nstep, nspec, nat, ndof, iprint
  logical         :: shift, remove_com_v, dump
  character(40)   :: position_file, dump_file
  integer, allocatable, dimension(:)          :: species
  real(double), allocatable, dimension(:,:)   :: v_t, v_t_dt
  real(double), allocatable, dimension(:,:)   :: f_t, f_t_dt
  real(double), dimension(3,3)                :: p_g_t, p_g_t_dt
  real(double), dimension(3,3)                :: P_int_t, P_int_t_dt, P_ext
  real(double), dimension(3)                  :: sumv
  real(double)   :: pe_t, pe_t_dt, ke_t, ke_t_dt, T_int_t, T_int_t_dt, T_ext, &
                    dt, sumv2

  contains
    procedure :: init_md
    procedure :: get_force_and_energy
    procedure :: get_kinetic_energy
    procedure :: vVerlet_v_half
    procedure :: vVerlet_r
    procedure :: md_run
    procedure :: update_v
    procedure :: dump_atom_arr
    procedure :: md_dump
end type type_md

contains

  ! Initialise variables/velocities/allocate matrices for MD run
  subroutine init_md(mdr, init_cell, pp, init_cell_cart, ensemble, nstep, dt, &
                     T_ext, shift, remove_com_v)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    type(type_cell), intent(in)     :: init_cell
    type(type_pp), intent(in)       :: pp
    character(3), intent(in)        :: ensemble
    integer, intent(in)             :: nstep
    real(double)                    :: dt
    real(double)                    :: T_ext
    logical, intent(in)             :: init_cell_cart
    logical, intent(in)             :: shift
    logical, intent(in)             :: remove_com_v

    ! local variables
    integer                         :: i, j
    real(double), dimension(3)      :: sumv
    real(double)                    :: sumv2, sfac

    mdr%position_file = "trajectory.xsf"
    mdr%dump_file = "dump.out"
    write(*,'(a)') "Starting MD run"
    write(*,*)
    mdr%p_t = init_cell
    if (init_cell_cart .eqv. .false.) call mdr%p_t%cell_frac2cart
    mdr%nspec = mdr%p_t%nspec
    call mdr%p_t_dt%init(mdr%p_t%nat, mdr%p_t%nspec, mdr%p_t%h, &
                         mdr%p_t%r, mdr%p_t%species)
    mdr%pp = pp
    mdr%iprint = 0
    mdr%ensemble = ensemble
    mdr%nat = init_cell%nat
    mdr%nstep = nstep
    mdr%dt = dt
    mdr%T_ext = T_ext
    mdr%shift = shift
    mdr%remove_com_v = remove_com_v
    mdr%dump = .true.

    allocate(mdr%species(mdr%nat))
    allocate(mdr%v_t(mdr%nat,3), mdr%v_t_dt(mdr%nat,3))
    allocate(mdr%f_t(mdr%nat,3), mdr%f_t_dt(mdr%nat,3))
    mdr%species = init_cell%spec_int

    write(*,'(a)') "Simulation parameters:"
    write(*,'("Ensemble               ",a8)') mdr%ensemble
    write(*,'("Number of steps        ",i8)') mdr%nstep
    write(*,'("Time step              ",f8.4)') mdr%dt
    write(*,'("Initial/external T     ",f8.4)') mdr%T_ext
    write(*,'("Remove COM velocity    ",l8)') mdr%remove_com_v
    write(*,*)
    write(*,'(a)') "Species and pair potential details:"
    write(*,'("Pair potential cutoff  ",f8.4)') mdr%pp%r_cut
    write(*,'("Pair potential shift   ",l8)') mdr%shift
    write(*,'(2x,a)') "Sigma"
    do i=1,mdr%nspec
      write(*,'(4x,3f10.6)') mdr%pp%sigma(i,:)
    end do
    write(*,'(a)') "Epsilon"
    do i=1,mdr%nspec
      write(*,'(4x,3f10.6)') mdr%pp%epsilon(i,:)
    end do
    write(*,*)
    write(*,'(2x,a)') "label, species, count"
    do i=1,mdr%nspec
      write(*,'(4x,i4,a4,i8)') mdr%p_t%spec_int(i), mdr%p_t%spec(i), &
                               mdr%p_t%spec_count(i)
    end do
    write(*,*)
    write(*,'(2x,a)') "Initial unit cell:"
    do i=1,3
      write(*,'(4x,3f12.6)') mdr%p_t%h(i,:)
    end do
    write(*,*)
    write(*,'(2x,a)') "Initial atomic positions:"
    do i=1,mdr%nat
      write(*,'(4x,2i6,3f14.8)') i, mdr%species(i), mdr%p_t%rcart(i,:)
    end do

    ! ensemble specifics
    select case (ensemble)
    case ('nve')
      mdr%ndof = 3*mdr%nat
    end select
    if (mdr%remove_com_v .eqv. .true.) mdr%ndof = mdr%ndof-1

    ! initialise rng
    call init_rand

    ! initialise velocities w/ uniform random distribution
    sumv2 = zero
    do i=1,mdr%nat
      do j=1,3
        call rand(mdr%v_t(i,j))
        mdr%v_t(i,j) = mdr%v_t(i,j) - half
      end do
    end do

    call mdr%update_v
    sfac = sqrt(mdr%ndof*mdr%T_ext/mdr%sumv2)
    mdr%v_t = sfac*mdr%v_t ! Scale velocities according to T_ext

    write(*,*)
    write(*,'(2x,a)') "Initial velocities:"
    do i=1,mdr%nat
      write(*,'(4x,2i6,3f14.8)') i, mdr%species(i), mdr%v_t(i,:)
    end do

    write(*,*)
    call mdr%get_kinetic_energy
    call mdr%get_force_and_energy

    write(*,*)
    write(*,'(2x,a)') "Initial forces:"
    do i=1,mdr%nat
      write(*,'(4x,2i6,3f14.8)') i, mdr%species(i), mdr%f_t(i,:)
    end do
    write(*,*)

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
        r_ij = mdr%p_t%rcart(jat,:) - mdr%p_t%rcart(iat,:)
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
    mdr%pe_t = mdr%pe_t*half
    write(*,'("  Potential energy = ",e16.8)') mdr%pe_t
  end subroutine get_force_and_energy

  ! Compute the kinetic energy
  subroutine get_kinetic_energy(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer   :: i

    mdr%ke_t = mdr%sumv2/two
    write(*,'("  Kinetic energy   = ",e16.8)') mdr%ke_t

  end subroutine get_kinetic_energy

  ! Velocity Verlet dt/2 step for velocities
  subroutine vVerlet_v_half(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    if (mdr%iprint == 0) write(*,'(a)') "Velocity Verlet velocity dt/2 update"
    mdr%v_t_dt = mdr%v_t + mdr%dt*half*(mdr%f_t + mdr%f_t_dt)

  end subroutine vVerlet_v_half

  ! Velocity Verlet dt step for positions
  subroutine vVerlet_r(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    if (mdr%iprint == 0) write(*,'(a)') "Velocity Verlet position dt update"
    mdr%p_t_dt%rcart = mdr%p_t%rcart + mdr%dt*mdr%v_t_dt
    ! wrap the updated positions back into the unit cell
    if (mdr%iprint == 0) write(*,'(a)') "Wrapping positions into unit cell"
    call mdr%p_t_dt%wrap_positions_cart

  end subroutine vVerlet_r

  ! The main MD loop
  subroutine md_run(mdr, s_start, s_end)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    integer, intent(in)             :: s_start, s_end

    ! local variables
    integer   :: s, traj_unit, dump_unit

    traj_unit = 101
    dump_unit = 102
    s = 0

    open(unit=traj_unit, file=mdr%position_file, status='replace')
    open(unit=dump_unit, file=mdr%dump_file, status='replace')
    call mdr%p_t%write_xsf(traj_unit, .true., 1, s_end+1)
    if (mdr%dump .eqv. .true.) call mdr%md_dump(dump_unit, s)
    do s=s_start,s_end
      write(*,*)
      write(*,'("MD step ",i6," of ",i6)')  s, s_end
      call mdr%vVerlet_v_half
      call mdr%vVerlet_r
      call mdr%p_t_dt%write_xsf(traj_unit, .true., s+1, s_end+1)
      call mdr%vVerlet_v_half
      call mdr%update_v
      call mdr%get_kinetic_energy
      call mdr%get_force_and_energy
      if (mdr%dump .eqv. .true.) call mdr%md_dump(dump_unit, s)
      mdr%v_t = mdr%v_t_dt
    end do
    close(traj_unit)
    close(dump_unit)

  end subroutine md_run

  ! Reomve centre of mass velocity, compute velocity and sum of velocities
  ! squared
  subroutine update_v(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer :: i

    ! update the COM velocity
    do i=1,3
      mdr%sumv(i) = sum(mdr%v_t(:,i))
    end do
    mdr%sumv = mdr%sumv/mdr%nat

    do i=1,mdr%nat
      mdr%sumv2 = mdr%sumv2 + sum(mdr%v_t(i,:)**2)
    end do

    ! remove COM velocity
    if (mdr%remove_com_v .eqv. .true.) then
      do i=1,mdr%nat
        mdr%v_t(i,:) = mdr%v_t(i,:) - mdr%sumv
      end do
    end if

    ! update the square of the velocity (for kinetic energy etc)
    mdr%sumv2 = zero
    do i=1,mdr%nat
      mdr%sumv2 = mdr%sumv2 + sum(mdr%v_t(i,:)**2)
    end do
  end subroutine update_v

  ! Dump the an atom array (position, force, velocity)
  subroutine dump_atom_arr(mdr, iou, arr)

    ! passed variables
    class(type_md), intent(inout)             :: mdr
    integer, intent(in)                       :: iou
    real(double), dimension(:,:), intent(in)  :: arr

    ! local variables
    integer                         :: i
    character(80)                   :: fmt

    fmt = "(2i5,3f16.10)"
    do i=1,mdr%nat
      write(iou,fmt) i, mdr%species(i), arr(i,:)
    end do

  end subroutine dump_atom_arr

  subroutine md_dump(mdr, iunit, step)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    integer, intent(in)             :: iunit
    integer, intent(in)             :: step

    ! local variables
    integer                         :: i

    write(iunit,'("step ",i8)') step
    write(iunit,'(a)') "cell_vectors"
    do i=1,3
      write(iunit,'(3f12.6)') mdr%p_t%h(i,:)
    end do
    write(iunit,'(a)') "end cell_vectors"
    write(iunit,'(a)') "positions"
    call mdr%dump_atom_arr(iunit, mdr%p_t%rcart)
    write(iunit,'(a)') "end positions"
    write(iunit,'(a)') "velocities"
    call mdr%dump_atom_arr(iunit, mdr%v_t)
    write(iunit,'(a)') "end velocities"
    write(iunit,'(a)') "forces"
    call mdr%dump_atom_arr(iunit, mdr%f_t)
    write(iunit,'(a)') "end forces"

  end subroutine md_dump

end module md_module
