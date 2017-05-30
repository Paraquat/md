module md_module

use datatypes
use constants
use vector
use cell
use pairpotential
use rng
use thermostat_module

implicit none

type type_md
  type(type_cell)           :: p_t, p_t_dt
  type(type_pairpotential)  :: pp
  type(type_thermostat)     :: th
  character(3)    :: ensemble
  integer         :: nstep, nspec, nat, ndof, iprint, dump_freq, tau_T
  logical         :: shift, remove_com_v, dump
  character(40)   :: position_file, dump_file, stat_file, therm_type, &
                     init_distr, units
  integer, allocatable, dimension(:)          :: species
  real(double), allocatable, dimension(:)     :: mass
  real(double), allocatable, dimension(:,:)   :: v_t, v_t_dt
  real(double), allocatable, dimension(:,:)   :: f
  real(double), dimension(3,3)                :: p_g_t, p_g_t_dt
  real(double), dimension(3,3)                :: P_int_t, P_int_t_dt, P_ext
  real(double), dimension(3)                  :: sumv
  real(double)   :: pe, ke, T_int, T_ext, dt, sumv2, k_B_md

  contains
    procedure :: init_md
    procedure :: init_velocities
    procedure :: update_v
    procedure :: get_force_and_energy
    procedure :: get_kinetic_energy
    procedure :: get_temperature
    procedure :: get_pressure
    procedure :: vVerlet_v_half
    procedure :: vVerlet_r
    procedure :: md_run
    procedure :: fire
    procedure :: dump_atom_arr
    procedure :: md_dump
    procedure :: stat_dump
end type type_md

contains

  ! Initialise variables/velocities/allocate matrices for MD run
  subroutine init_md(mdr, init_cell, pp, init_cell_cart, ensemble, nstep, dt, &
                     T_ext, vdistr, shift, remove_com_v, thermo_type)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    type(type_cell), intent(in)     :: init_cell
    type(type_pairpotential), intent(in)  :: pp
    character(3), intent(in)        :: ensemble
    integer, intent(in)             :: nstep
    real(double)                    :: dt
    real(double)                    :: T_ext
    logical, intent(in)             :: init_cell_cart
    character(40), intent(in)       :: vdistr
    character(40), intent(in)       :: thermo_type
    logical, intent(in)             :: shift
    logical, intent(in)             :: remove_com_v

    ! local variables
    integer                         :: i, j

    mdr%position_file = "trajectory.xsf"
    mdr%dump_file = "dump.out"
    mdr%stat_file = "stat.out"
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
    mdr%init_distr = vdistr
    mdr%T_ext = T_ext
    mdr%shift = shift
    mdr%remove_com_v = remove_com_v
    mdr%dump = .true.
    mdr%dump_freq = 1
    mdr%therm_type = 'none'
    mdr%units = 'reduced-lj'

    select case (mdr%units)
    case ('reduced-lj')
      mdr%k_B_md = one
    case default
      mdr%k_B_md = one
    end select

    allocate(mdr%species(mdr%nat))
    allocate(mdr%mass(mdr%nat))
    allocate(mdr%v_t(mdr%nat,3), mdr%v_t_dt(mdr%nat,3))
    allocate(mdr%f(mdr%nat,3))
    mdr%species = init_cell%spec_int
    do i=1,mdr%nat
      mdr%mass(i) = mdr%p_t%mass(mdr%p_t%spec_int(i))
    end do

    ! ensemble specifics
    mdr%ndof = 3*mdr%nat
    select case (ensemble)
    case ('nvt')
      ! set ndof for extended lagrangian systems here
      call mdr%th%init_thermostat(thermo_type, mdr%nat, mdr%nat, mdr%T_ext, &
                                  mdr%iprint)
    end select
    if (mdr%remove_com_v .eqv. .true.) mdr%ndof = mdr%ndof-3

    write(*,'(a)') "Simulation parameters:"
    write(*,'("Ensemble               ",a8)') mdr%ensemble
    write(*,'("Number of steps        ",i8)') mdr%nstep
    write(*,'("Time step              ",f8.4)') mdr%dt
    write(*,'("Remove COM velocity    ",l8)') mdr%remove_com_v
    if (mdr%ensemble == 'nvt' .or. mdr%ensemble == 'npt') then
      write(*,'("Thermostat             ",a16)') mdr%therm_type
      write(*,'("Thermostat period      ",i8)') mdr%tau_T
    end if
    write(*,*)
    write(*,'(a)') "Initialisation:"
    write(*,'("Velocity distribution  ",a20)') mdr%init_distr
    write(*,'("Initial/external T     ",f8.4)') mdr%T_ext
    write(*,*)
    write(*,'(a)') "Species and pair potential details:"
    write(*,'("Pair potential type    ",a)') mdr%pp%potential_type
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
    if (mdr%pp%potential_type == "morse") then
      write(*,'(2x,a)') "r_e"
      do i=1,mdr%nspec
        write(*,'(4x,3f10.6)') mdr%pp%r_e(i,:)
      end do
    end if
    write(*,*)
    write(*,'(2x,3a8,a10)') "label", "species", "count", "mass"
    do i=1,mdr%nspec
      write(*,'(2x,i8,a8,i8,f10.4)') mdr%p_t%spec_int(i), mdr%p_t%spec(i), &
                               mdr%p_t%spec_count(i), mdr%p_t%mass(i)
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

    ! initialise velocity
    call mdr%init_velocities(mdr%init_distr)

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
      write(*,'(4x,2i6,3f14.8)') i, mdr%species(i), mdr%f(i,:)
    end do
    write(*,*)

  end subroutine init_md

  subroutine init_velocities(mdr, distr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    character(40), intent(in)       :: distr

    ! local variables
    real(double)                    :: sfac, sigma, mu, z1, z2
    integer                         :: i, j, k
    real(double), allocatable, dimension(:) :: rn

    ! initialise rng
    call init_rand

    ! initialise velocities w/ uniform random distribution
    select case (distr)
    case ('uniform')
      do i=1,mdr%nat
        do j=1,3
          call rand(mdr%v_t(i,j))
          mdr%v_t(i,j) = mdr%v_t(i,j) - half
        end do
      end do
    case ('maxwell-boltzmann')
      allocate(rn(3*mdr%nat))
      sigma = sqrt(mdr%T_ext/two)
      mu = zero
      i = 1
      do
        call boxmuller_polar(sigma, mu, z1, z2)
        if (i <= 3*mdr%nat) rn(i) = z1
        if (i <= 3*mdr%nat) rn(i) = z2
        i = i+2
        if (i > mdr%nat) exit
      end do
      k = 1
      do i=1,mdr%nat
        do j=1,3
          mdr%v_t(i,j) = rn(k)
          k = k+1
        end do
      end do
      deallocate(rn)
    end select

    call mdr%update_v(mdr%remove_com_v)
    write(*,'(a)') "Before rescaling:"
    call mdr%get_kinetic_energy
    call mdr%get_temperature
    sfac = sqrt(mdr%T_ext/mdr%T_int)
    write(*,*)
    ! sfac = sqrt(mdr%ndof*mdr%T_ext/mdr%sumv2)
    mdr%v_t = sfac*mdr%v_t ! Scale velocities according to T_ext
    write(*,'(a)') "After rescaling:"
    call mdr%get_kinetic_energy
    call mdr%get_temperature
  end subroutine init_velocities

  ! Reomve centre of mass velocity, compute velocity and sum of velocities
  ! squared
  subroutine update_v(mdr, remove_com_v)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    logical, intent(in)             :: remove_com_v

    ! local variables
    integer                         :: i

    ! update the COM velocity
    do i=1,3
      mdr%sumv(i) = sum(mdr%v_t(:,i))
    end do
    mdr%sumv = mdr%sumv/mdr%nat

    ! remove COM velocity
    if (remove_com_v .eqv. .true.) then
      do i=1,mdr%nat
        mdr%v_t(i,:) = mdr%v_t(i,:) - mdr%sumv
      end do
    end if

    ! update the sum of the squared velocity (for kinetic energy etc)
    mdr%sumv2 = zero
    do i=1,mdr%nat
      mdr%sumv2 = mdr%sumv2 + sum(mdr%v_t(i,:)**2)
    end do
  end subroutine update_v

  ! Compute the potential energy and force on each atom
  subroutine get_force_and_energy(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer                         :: iat, jat, s_i, s_j
    real(double), dimension(3)      :: r_ij_cart
    real(double)                    :: mod_r_ij, mod_f, pe

    mdr%pe = zero
    mdr%f = zero

    do iat=1,mdr%nat
      do jat=1,mdr%nat
        if (iat == jat) cycle
        s_i = mdr%species(iat)
        s_j = mdr%species(jat)
        r_ij_cart = mdr%p_t_dt%mic(mdr%p_t_dt%rcart(iat,:), &
                                   mdr%p_t_dt%rcart(jat,:))
        mod_r_ij = modulus(r_ij_cart)
        if (mod_r_ij < mdr%pp%r_cut(s_i,s_j)) then
          r_ij_cart = norm(r_ij_cart)
          call mdr%pp%pp_force_and_energy(mod_r_ij, s_i, s_j, mdr%shift, &
                                          mod_f, pe)
          mdr%f(iat,:) = mdr%f(iat,:) + mod_f*r_ij_cart
          mdr%pe = mdr%pe + pe
        end if
      end do
    end do
    mdr%pe = mdr%pe*half
    write(*,'("  Potential energy = ",e16.8)') mdr%pe
  end subroutine get_force_and_energy

  ! Compute the kinetic energy
  subroutine get_kinetic_energy(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer   :: i

    mdr%ke = zero
    do i=1,mdr%nat
      mdr%ke = mdr%ke + mdr%mass(i)*sum(mdr%v_t(i,:)**2)
    end do
    mdr%ke = mdr%ke/two
    write(*,'("  Kinetic energy   = ",e16.8)') mdr%ke

  end subroutine get_kinetic_energy

  ! Compute the temperature
  subroutine get_temperature(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    mdr%T_int = two*mdr%ke/mdr%k_B_md/real(mdr%ndof, double)
    write(*,'("  Temperature      = ",f16.8)') mdr%T_int
  end subroutine get_temperature

  ! Compute the pressure using the virial
  subroutine get_pressure(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer                         :: i

    do i=1,mdr%nat
    end do

  end subroutine get_pressure

  ! Velocity Verlet dt/2 step for velocities
  subroutine vVerlet_v_half(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer                         :: i

    if (mdr%iprint == 0) write(*,'(a)') "Velocity Verlet velocity dt/2 update"
    do i=1,mdr%nat
      mdr%v_t_dt(i,:) = mdr%v_t(i,:) + mdr%dt*half*mdr%f(i,:)/mdr%mass(i)
    end do

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
    integer       :: s, traj_unit, dump_unit, stat_unit
    real(double)  :: total_energy
    logical       :: d

    traj_unit = 101
    dump_unit = 102
    stat_unit = 103
    s = 0

    open(unit=traj_unit, file=mdr%position_file, status='replace')
    open(unit=dump_unit, file=mdr%dump_file, status='replace')
    open(unit=stat_unit, file=mdr%stat_file, status='replace')
    call mdr%p_t%write_xsf(traj_unit, .true., 1, (s_end/mdr%dump_freq)+1)
    call mdr%stat_dump(stat_unit, s)
    if (mdr%dump .eqv. .true.) then
      if (mdr%dump_freq == 0) call mdr%md_dump(dump_unit, s)
    end if

    ! Main MD loop
    do s=s_start,s_end
      d = .false.
      if (mod(s,mdr%dump_freq) == 0) d = .true.
      write(*,*)
      write(*,'("MD step ",i10," of ",i10)')  s, s_end
      call mdr%vVerlet_v_half
      mdr%v_t = mdr%v_t_dt

      ! velocity Verlet algorithm
      call mdr%vVerlet_r
      if (d .eqv. .true.) then
        call mdr%p_t_dt%write_xsf(traj_unit, .true., s+1, s_end+1)
      end if
      call mdr%get_force_and_energy
      call mdr%vVerlet_v_half

      ! Update arrays and thermodyanmics quantities
      mdr%v_t = mdr%v_t_dt

      call mdr%update_v(mdr%remove_com_v)

      ! Thermostat: velocity update
      if (mdr%ensemble == 'nvt' .or. mdr%ensemble == 'npt') then
        if (mod(s,mdr%tau_T) == 0) then
          call mdr%th%propagate_thermostat(mdr%T_int, mdr%v_t)
        end if
      end if

      call mdr%get_kinetic_energy
      write(*,'("  Total energy     = ",e16.8)') total_energy
      call mdr%get_temperature
      total_energy = mdr%ke + mdr%pe
      mdr%p_t%rcart = mdr%p_t_dt%rcart
      if (mdr%ensemble == 'npt' .or. mdr%ensemble == 'nph') then
        mdr%p_t%h = mdr%p_t_dt%h
      end if

      call mdr%stat_dump(stat_unit, s)
      if (mdr%dump .eqv. .true.) then
        if (d .eqv. .true.) then
          call mdr%md_dump(dump_unit, s)
          call flush(traj_unit)
          call flush(dump_unit)
          call flush(stat_unit)
        end if
      end if
    end do
    close(traj_unit)
    close(dump_unit)
    close(stat_unit)

  end subroutine md_run

  ! Perform geometry optimisation using FIRE
  subroutine fire(mdr, alpha)

    ! passed variables
    class(type_md), intent(inout)             :: mdr
    real(double) , intent(in)                 :: alpha ! F/v mixing paramter

    ! local variables
    real(double)                              :: p, norm_v, norm_f
    integer                                   :: i

    ! dot product of force and velocity, norms
    p = sum(mdr%v_t*mdr%f)
    norm_v = sqrt(sum(mdr%v_t**2))
    norm_f = sqrt(sum(mdr%f**2))

    if (p > zero) then ! we're still going downhill
      ! update velocities
      mdr%v_t = (one - alpha)*mdr%v_t + alpha*norm_v*mdr%f
    else  ! decrease time step, reset alpha
    end if

  end subroutine fire

  ! Dump the an atom array (position, force, velocity)
  subroutine dump_atom_arr(mdr, iou, arr)

    ! passed variables
    class(type_md), intent(inout)             :: mdr
    integer, intent(in)                       :: iou
    real(double), dimension(:,:), intent(in)  :: arr

    ! local variables
    integer                         :: i
    character(80)                   :: fmt

    fmt = "(2i5,3e20.10)"
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
    call mdr%dump_atom_arr(iunit, mdr%f)
    write(iunit,'(a)') "end forces"

  end subroutine md_dump

  ! Dump thermodyanmic statistics
  subroutine stat_dump(mdr, iunit, step)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    integer, intent(in)             :: iunit
    integer, intent(in)             :: step

    if (step == 0) then
      select case (mdr%ensemble)
      case ('nve')
        write(iunit,'(a10,4a16)') "step", "potential", "kinetic", "total", "T"
      end select
    end if
    select case (mdr%ensemble)
    case ('nve')
      write(iunit,'(i10,4e16.6)') step, mdr%pe, mdr%ke, mdr%pe+mdr%ke, mdr%T_int
    end select
  end subroutine stat_dump

end module md_module
