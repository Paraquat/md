module md_module

use datatypes
use constants
use vector
use cell
use pairpotential
use rng
use thermostat_module
use barostat_module
use isotropic_barostat_module

implicit none

type type_md
  type(type_cell)           :: p_t, p_t_dt
  type(type_pairpotential)  :: pp
  type(type_thermostat)     :: th
  type(type_barostat)       :: baro
  type(type_iso_barostat)   :: iso_baro
  character(3)    :: ensemble
  integer         :: nstep, nspec, nat, ndof, iprint, dump_freq, tau_T, tau_P
  integer         :: n_nhc
  logical         :: shift, remove_com_v, dump
  character(40)   :: position_file, dump_file, stat_file, thermo_type, &
                     baro_type, init_distr, units
  integer, allocatable, dimension(:)          :: species
  real(double), allocatable, dimension(:)     :: mass
  real(double), allocatable, dimension(:,:)   :: v_t, v_t_dt
  real(double), allocatable, dimension(:,:)   :: f
  real(double), dimension(3,3)                :: stress, virial_ke, &
                                                 virial_tensor
  real(double), dimension(3)                  :: sumv
  real(double)   :: pe, ke, T_int, T_ext, dt, sumv2, k_B_md, virial, P_int, &
                    P_ext, V, e_nhc, ke_box, pv, e_npt, e_nvt, e_nve

  contains
    procedure :: init_md
    procedure :: init_velocities
    procedure :: update_v
    procedure :: get_force_and_energy
    procedure :: get_kinetic_energy
    procedure :: get_pv
    procedure :: get_temperature
    procedure :: get_pressure
    procedure :: get_stress
    procedure :: vVerlet_v
    procedure :: vVerlet_r
    procedure :: md_run
    procedure :: propagate_npt_mttk
    procedure :: propagate_iso_npt_mttk
    procedure :: fire
    procedure :: dump_atom_arr
    procedure :: md_dump
    procedure :: stat_dump
end type type_md

contains

  ! Initialise variables/velocities/allocate matrices for MD run
  ! TODO: Altogether too many arguments here, fix this.
  subroutine init_md(mdr, init_cell, pp, init_cell_cart, ensemble, nstep, dt, &
                     T_ext, vdistr, shift, remove_com_v, thermo_type, tau_T, &
                     n_nhc, nhc_mass, baro_type, P_ext, box_mass)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    type(type_cell), intent(in)     :: init_cell
    type(type_pairpotential), intent(in)  :: pp
    character(3), intent(in)        :: ensemble
    integer, intent(in)             :: nstep
    real(double), intent(in)        :: dt
    real(double), intent(in)        :: T_ext
    integer, intent(in)             :: tau_T
    logical, intent(in)             :: init_cell_cart
    character(40), intent(in)       :: vdistr
    character(40), intent(in)       :: thermo_type
    integer, intent(in)             :: n_nhc
    real(double), dimension(:), intent(in) :: nhc_mass
    logical, intent(in)             :: shift
    logical, intent(in)             :: remove_com_v
    character(40), intent(in)       :: baro_type
    real(double), intent(in)        :: P_ext
    real(double), intent(in)        :: box_mass

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
    mdr%p_t_dt%V = mdr%p_t%V
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
    mdr%units = 'reduced-lj'
    mdr%V = mdr%p_t%volume()
    mdr%thermo_type = thermo_type
    mdr%tau_T = tau_T
    mdr%n_nhc = n_nhc
    mdr%baro_type = baro_type
    mdr%P_ext = P_ext
    mdr%tau_P = 1

    select case (mdr%units)
    case ('reduced-lj')
      mdr%k_B_md = one
    case default
      mdr%k_B_md = k_B
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
    if (mdr%remove_com_v .eqv. .true.) mdr%ndof = mdr%ndof-3
    ! constant temperature
    if (mdr%ensemble(3:3) == 't') then
      ! set ndof for extended lagrangian systems
      if (mdr%thermo_type == 'nhc') mdr%ndof = mdr%ndof + mdr%n_nhc
      call mdr%th%init_thermostat(thermo_type, mdr%nat, mdr%ndof, mdr%T_ext, &
                                  mdr%tau_T, mdr%n_nhc, mdr%iprint, mdr%ke)
      mdr%th%Q = nhc_mass
    end if
    ! constant pressure
    if (mdr%ensemble(2:2) == 'p') then
      if (mdr%baro_type == 'mttk') then ! fully flexible cell
        mdr%ndof = mdr%ndof + 9
        call mdr%baro%init_barostat(mdr%baro_type, mdr%P_ext, mdr%nat, &
                                    mdr%ndof, box_mass, mdr%iprint)
      else if (mdr%baro_type == 'iso-mttk') then ! isotropic cell variation
        mdr%ndof = mdr%ndof + 1
        call mdr%iso_baro%init_barostat_iso(mdr%baro_type, mdr%p_t%h, &
                                            mdr%P_ext, mdr%nat, mdr%ndof, &
                                            box_mass, mdr%V, mdr%iprint)
      end if
    end if

    write(*,'(a)') "Simulation parameters:"
    write(*,'("Ensemble               ",a8)') mdr%ensemble
    write(*,'("Number of steps        ",i8)') mdr%nstep
    write(*,'("Time step              ",f8.4)') mdr%dt
    write(*,'("Remove COM velocity    ",l8)') mdr%remove_com_v
    if (mdr%ensemble(3:3) == 't') then
      write(*,*)
      write(*,'(a)') "Thermostat details:"
      write(*,'("Thermostat type        ",a16)') mdr%thermo_type
      write(*,'("Thermostat period      ",i8)') mdr%tau_T
      if (mdr%thermo_type == 'nhc') then
        write(*,'("NH chain length        ",i8)') mdr%n_nhc
        if (mdr%iprint == 0) then
        write(*,'(a)') "NHC heat bath parameters:"
          write(*,'(3a12)') 'eta', 'v_eta', 'Q'
          do i=1,mdr%n_nhc
            write(*,'(3f12.4)') mdr%th%eta(i), mdr%th%v_eta(i), mdr%th%Q(i)
          end do
        end if
      end if
    end if
    if (mdr%ensemble(2:2) == 'p') then
      write(*,*)
      write(*,'(a)') "Barostat details:"
      write(*,'("Barostat type          ",a16)') mdr%baro_type
      write(*,'("Barostat period        ",i8)') mdr%tau_P
      if (mdr%baro_type == 'mttk') then
        write(*,'("Barostat mass          ",f8.4)') box_mass
      end if
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
    write(*,'(2x,"Volume = ",f16.6)') mdr%V
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
    call mdr%get_force_and_energy
    call mdr%get_pressure
    call mdr%get_stress

    write(*,*)
    write(*,'(2x,a)') "Initial forces:"
    do i=1,mdr%nat
      write(*,'(4x,2i6,3f14.8)') i, mdr%species(i), mdr%f(i,:)
    end do
    write(*,*)
    call mdr%get_pressure
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
        ! Box-Muller transform generates 2 random numbers in a Gaussian
        ! distribution, so I'm initialising two components at a time.
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

  ! Compute the potential energy and force on each atom, update the virial
  subroutine get_force_and_energy(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer                         :: iat, jat, s_i, s_j, mu, nu
    real(double), dimension(3)      :: r_ij_cart, f_ij
    real(double)                    :: mod_r_ij, mod_f, pe

    mdr%pe = zero
    mdr%f = zero
    mdr%virial = zero
    mdr%virial_tensor = zero

    do iat=1,mdr%nat
      do jat=iat+1,mdr%nat
        s_i = mdr%species(iat)
        s_j = mdr%species(jat)
        r_ij_cart = mdr%p_t_dt%mic(mdr%p_t_dt%rcart(iat,:), &
                                   mdr%p_t_dt%rcart(jat,:))
        mod_r_ij = modulus(r_ij_cart)
        if (mod_r_ij < mdr%pp%r_cut(s_i,s_j)) then
          call mdr%pp%pp_force_and_energy(mod_r_ij, s_i, s_j, mdr%shift, &
                                          mod_f, pe)
          f_ij = mod_f*norm(r_ij_cart)
          mdr%f(iat,:) = mdr%f(iat,:) + f_ij
          mdr%f(jat,:) = mdr%f(jat,:) - f_ij
          mdr%pe = mdr%pe + pe
          ! compute virial
          mdr%virial = mdr%virial + dot_product(f_ij, r_ij_cart)
          mdr%virial_tensor = mdr%virial_tensor + &
                              tensor_product(f_ij, r_ij_cart)
        end if
      end do
    end do
    ! mdr%pe = mdr%pe*half
    mdr%virial = mdr%virial/mdr%V/three ! The configurational part of the &
                                      ! virial, used in iso-mttk barostat
    write(*,'("  Potential energy = ",e16.8)') mdr%pe
  end subroutine get_force_and_energy

  ! Compute the kinetic energy
  subroutine get_kinetic_energy(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer   :: i, mu, nu

    mdr%ke = zero
    mdr%virial_ke = zero
    do i=1,mdr%nat
      mdr%ke = mdr%ke + mdr%mass(i)*sum(mdr%v_t(i,:)**2)
      mdr%virial_ke = mdr%virial_ke + &
                      mdr%mass(i)*tensor_product(mdr%v_t(i,:), mdr%v_t(i,:))
    end do
    mdr%ke = mdr%ke/two
    write(*,'("  Kinetic energy   = ",e16.8)') mdr%ke

  end subroutine get_kinetic_energy

  ! Compute the pV term for the enthalpy
  subroutine get_pv(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    mdr%V = mdr%p_t%volume()
    mdr%pv = mdr%v*mdr%P_ext ! Should use the target pressure according to MTTK

  end subroutine get_pv

  ! Compute the temperature
  subroutine get_temperature(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    mdr%T_int = two*mdr%ke/mdr%k_B_md/real(mdr%ndof, double)
    write(*,'("  Temperature      = ",f16.8)') mdr%T_int
  end subroutine get_temperature

  ! Compute the pressure using the virial
  ! NB only correct in the NVT ensemble. FIX THIS!
  subroutine get_pressure(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! use f_ij*r_ij computed from pp_get_force_and_energy
    mdr%P_int = mdr%nat*mdr%k_B_md*mdr%T_int/mdr%V + mdr%virial
    write(*,'("  Virial           = ",f16.8)') mdr%virial
    write(*,'("  Pressure         = ",f16.8)') mdr%P_int

  end subroutine get_pressure

  ! Compute the stress tensor from the virial tensor from get_force_and_energy
  ! and the kinetic energy from get_kinetic_energy
  subroutine get_stress(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer                         :: i

    mdr%stress = (one/mdr%V)*(-mdr%virial_ke + mdr%virial_tensor)
    if (mdr%iprint == 0) then
      write(*,*)
      write(*,'(4x,a36,4x,a36)') "Virial KE", "Virial r.F"
      do i=1,3
        write(*,'(4x,3f12.4,4x,3f12.4)') mdr%virial_ke(i,:), &
                                         mdr%virial_tensor(i,:)
      end do
      write(*,*)
    end if

    write(*,'(2x,a)') "Stress tensor:"
    do i=1,3
      write(*,'(4x,3f16.8)') mdr%stress(i,:)
    end do
  end subroutine get_stress

  ! Velocity Verlet dt/2 step for velocities
  subroutine vVerlet_v(mdr, dt, dtfac)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    real(double), intent(in)        :: dt, dtfac

    ! local variables
    integer                         :: i

    if (mdr%iprint == 0) write(*,'(4x,a)') "Velocity Verlet velocity update"
    do i=1,mdr%nat
      mdr%v_t_dt(i,:) = mdr%v_t(i,:) + dt*dtfac*mdr%f(i,:)/mdr%mass(i)
    end do

  end subroutine vVerlet_v

  ! Velocity Verlet dt step for positions
  subroutine vVerlet_r(mdr, dt, dtfac)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    real(double), intent(in)        :: dt, dtfac

    if (mdr%iprint == 0) write(*,'(4x,a)') "Velocity Verlet position update"
    mdr%p_t_dt%rcart = mdr%p_t%rcart + dt*dtfac*mdr%v_t_dt
    ! wrap the updated positions back into the unit cell
    if (mdr%iprint == 0) write(*,'(4x,a)') "Wrapping positions into unit cell"
    call mdr%p_t_dt%wrap_positions_cart

  end subroutine vVerlet_r

  ! The main MD loop
  subroutine md_run(mdr, s_start, s_end)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    integer, intent(in)             :: s_start, s_end

    ! local variables
    integer       :: s, traj_unit, dump_unit, stat_unit, debug_unit
    logical       :: d


    traj_unit = 101
    dump_unit = 102
    stat_unit = 103
    debug_unit = 104
    s = 0

    open(unit=traj_unit, file=mdr%position_file, status='replace')
    open(unit=dump_unit, file=mdr%dump_file, status='replace')
    open(unit=stat_unit, file=mdr%stat_file, status='replace')
    if (mdr%iprint == 0) open(debug_unit, file='baro_state.out', &
                              status='replace')
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
      if (mdr%ensemble == 'nvt') then
        if (mdr%thermo_type == 'nhc') then
          call mdr%th%propagate_nvt_nhc(mdr%dt, mdr%v_t_dt)
        end if
      else if (mdr%ensemble == 'npt') then
        if (mdr%thermo_type == 'nhc') then
          if (mdr%baro_type == 'mttk') then
            call mdr%propagate_npt_mttk(mdr%dt)
          else if (mdr%baro_type == 'iso-mttk') then
            call mdr%propagate_iso_npt_mttk(mdr%dt)
          end if
        end if
      end if
      call mdr%vVerlet_v(mdr%dt, half)   ! v_t -> v_t_dt
      mdr%v_t = mdr%v_t_dt

      if (mdr%ensemble(2:2) == 'p') then ! constant pressure
        if (mdr%baro_type == 'mttk') then
          call mdr%baro%vVerlet_r_h_npt(mdr%dt, mdr%p_t_dt%h, mdr%v_t, &
                                        mdr%p_t%rcart, mdr%th%v_eta(1))
          call mdr%p_t_dt%wrap_positions_cart
        else if (mdr%baro_type == 'iso-mttk') then
          call mdr%iso_baro%propagate_r_sys(mdr%dt, half, mdr%v_t_dt, &
                                            mdr%p_t%rcart, mdr%p_t_dt%rcart)
          call mdr%iso_baro%propagate_eps_1(mdr%dt, one)
          call mdr%iso_baro%propagate_box(mdr%p_t_dt%h, mdr%p_t_dt%V)
          call mdr%p_t_dt%wrap_positions_cart
        end if
      else
        call mdr%vVerlet_r(mdr%dt, one) ! velocity Verlet r update
      end if
      mdr%p_t%rcart = mdr%p_t_dt%rcart

      ! mdr%p_t%rcart = mdr%p_t_dt%rcart

      ! dump configuration to .xsf file
      if (d .eqv. .true.) then
        call mdr%p_t_dt%write_xsf(traj_unit, .true., s+1, s_end+1)
      end if
      call mdr%get_force_and_energy
      call mdr%vVerlet_v(mdr%dt, half)   ! v_t -> v_t_dt

      ! Thermostat: velocity update
      if (mdr%ensemble == 'nvt') then
        if (mdr%thermo_type == 'nhc') then
          call mdr%th%propagate_nvt_nhc(mdr%dt, mdr%v_t_dt)
          call mdr%th%get_nhc_ke
        else
          if (mod(s,mdr%tau_T) == 0) then
            call mdr%th%propagate_vr_thermostat(mdr%T_int, mdr%v_t_dt)
          end if
        end if
      else if (mdr%ensemble == 'npt') then
        if (mdr%thermo_type == 'nhc') then
          if (mdr%baro_type == 'mttk') then
            call mdr%propagate_npt_mttk(mdr%dt)
          else if (mdr%baro_type == 'iso-mttk') then
            call mdr%propagate_iso_npt_mttk(mdr%dt)
          end if
          call mdr%th%get_nhc_ke
        end if
      end if

      ! Update arrays and thermodyanmics quantities
      mdr%v_t = mdr%v_t_dt
      call mdr%update_v(mdr%remove_com_v)
      call mdr%get_kinetic_energy
      mdr%e_nve = mdr%ke + mdr%pe
      if (mdr%ensemble(3:3) == 't') mdr%th%ke_system = mdr%ke

      write(*,'("  Total energy     = ",e16.8)') mdr%e_nve
      if (mdr%thermo_type == 'nhc') then
        if (mdr%thermo_type == 'nhc') then
          mdr%e_nhc = mdr%th%e_nhc
          mdr%e_nvt = mdr%e_nve + mdr%e_nhc
          write(*,'("  NHC energy       = ",e16.8)') mdr%e_nhc
          write(*,'("  NVT conserved E  = ",e16.8)') mdr%e_nvt
        end if
      end if
      if (mdr%ensemble == 'npt') then
        if (mdr%baro_type == 'mttk') then
          mdr%ke_box = mdr%baro%ke_box
        else if (mdr%baro_type == 'iso-mttk') then
          mdr%ke_box = mdr%iso_baro%ke_box
        end if
        call mdr%get_pv
        mdr%e_npt = mdr%e_nve + mdr%e_nhc + mdr%pv + mdr%ke_box 
        write(*,'("  NHC energy       = ",e16.8)') mdr%e_nhc
        write(*,'("  Box energy       = ",e16.8)') mdr%ke_box
        write(*,'("  PV               = ",e16.8)') mdr%pv
        write(*,'("  NPT conserved E  = ",e16.8)') mdr%e_npt
      end if
      call mdr%get_temperature
      call mdr%get_pressure
      call mdr%get_stress
      mdr%p_t%rcart = mdr%p_t_dt%rcart
      mdr%p_t%h = mdr%p_t_dt%h
      mdr%p_t%V = mdr%p_t_dt%V
      mdr%V = mdr%p_t_dt%V

      call mdr%stat_dump(stat_unit, s)
      if (mdr%dump .eqv. .true.) then
        if (d .eqv. .true.) then
          call mdr%md_dump(dump_unit, s)
          call flush(traj_unit)
          call flush(dump_unit)
          call flush(stat_unit)
        end if
      end if
      if (mdr%iprint == 0) then
        if (mdr%ensemble(2:2) == 'p') then
          if (mdr%baro_type == 'mttk') then
            call mdr%baro%dump_baro_state(s, debug_unit)
          else if (mdr%baro_type == 'iso-mttk') then
            call mdr%iso_baro%dump_baro_state_iso(s, debug_unit)
          end if
        end if
      end if
    end do
    close(traj_unit)
    close(dump_unit)
    close(stat_unit)
    if (mdr%iprint == 0) close(debug_unit)

  end subroutine md_run

  ! The velocity part of the update for the MTTK integrator. I'm putting this
  ! here because it requires thermostat/barostat coupling
  subroutine propagate_npt_mttk(mdr, dt)

    ! passed variables
    class(type_md), intent(inout)             :: mdr
    real(double), intent(in)                  :: dt

    ! local variables
    integer         :: i,j, k
    real(double)    :: dtys ! Yoshida-Suzuki time step
    real(double)    :: v_sfac ! velocity scaling factor
    real(double)    :: G_nhc_1 ! force on first NHC thermostat
    real(double), dimension(3,3) :: P_int

    ! Get the kinetic energy
    ! Update forces on thermostat and barostat
    call mdr%baro%update_G_h(mdr%ke, mdr%virial_ke, mdr%stress, mdr%V, &
                             mdr%T_ext, mdr%th%Q(1), G_nhc_1)
    mdr%th%G_nhc(1) = G_nhc_1
    do i=1,mdr%th%mts_nhc ! MTS loop
      do j=1,mdr%th%n_ys ! Yoshida-Suzuki loop
        ! Update thermostat velocities
        do k=mdr%th%n_nhc,1,-1
          if (k==mdr%th%n_nhc) then
            call mdr%th%update_G_k(k, mdr%baro%ke_box) ! should this be here?
            call mdr%th%propagate_v_eta_k_1(k, dt, quarter)
          else
            call mdr%th%propagate_v_eta_k_2(k, dt, eighth)
            call mdr%th%propagate_v_eta_k_1(k, dt, quarter)
            call mdr%th%update_G_k(k, mdr%baro%ke_box) ! should this be here?
            call mdr%th%propagate_v_eta_k_2(k, dt, eighth)
          end if
        end do
        ! Update box velocities
        call mdr%baro%propagate_v_h_2(dt, eighth, mdr%th%v_eta(1))
        call mdr%baro%propagate_v_h_1(dt, quarter)
        call mdr%baro%propagate_v_h_2(dt, eighth, mdr%th%v_eta(1))
        ! Update thermostat positions
        do k=1,mdr%th%n_nhc
          call mdr%th%propagate_eta_k(k, dt, half)
        end do
        ! Update Particle velocities
        call mdr%baro%propagate_v_sys(dt, half, mdr%th%v_eta(1), mdr%v_t)
        mdr%v_t = mdr%v_t_dt
        ! Get total ke
        call mdr%get_kinetic_energy
        ! Update the box forces
        P_int = mdr%virial_tensor/mdr%V
        call mdr%baro%update_G_h(mdr%ke, mdr%virial_ke, P_int, mdr%V, &
                                 mdr%T_ext, mdr%th%Q(1), G_nhc_1)
        ! Update the box velocities
        call mdr%baro%propagate_v_h_2(dt, eighth, mdr%th%v_eta(1))
        call mdr%baro%propagate_v_h_1(dt, quarter)
        call mdr%baro%propagate_v_h_2(dt, eighth, mdr%th%v_eta(1))
        ! Update the box ke
        call mdr%baro%get_box_ke
        ! Update the thermostat forces
        do k=1,mdr%th%n_nhc
          call mdr%th%update_G_k(k, mdr%baro%ke_box)
        end do
        ! Update the thermostat velocities
        do k=1,mdr%th%n_nhc
          if (k<mdr%th%n_nhc) then
            call mdr%th%propagate_v_eta_k_2(k, dt, eighth)
            call mdr%th%update_G_k(k, mdr%baro%ke_box)
            call mdr%th%propagate_v_eta_k_1(k, dt, quarter)
            call mdr%th%propagate_v_eta_k_2(k, dt, eighth)
          else
            call mdr%th%update_G_k(k, mdr%baro%ke_box)
            call mdr%th%propagate_v_eta_k_1(k, dt, quarter)
          end if
        end do
      end do ! MTS loop
    end do ! Yoshida-Suzuki loop

  end subroutine propagate_npt_mttk

  ! The velocity part of the update for the MTTK integrator. I'm putting this
  ! here because it requires thermostat/barostat coupling
  subroutine propagate_iso_npt_mttk(mdr, dt)

    ! passed variables
    class(type_md), intent(inout)             :: mdr
    real(double), intent(in)                  :: dt

    ! local variables
    integer         :: i,j, k
    real(double)    :: dtys ! Yoshida-Suzuki time step
    real(double)    :: G_nhc_1 ! force on first NHC thermostat
    real(double)    :: v_sfac, vscale

    v_sfac = one
    ! Get the kinetic energy
    call mdr%get_kinetic_energy
    ! Update forces on thermostat and barostat
    call mdr%iso_baro%update_G_eps(mdr%ke, mdr%virial, mdr%V, mdr%T_ext, &
                                   mdr%th%Q(1), G_nhc_1)
    mdr%th%G_nhc(1) = G_nhc_1
    do i=1,mdr%th%mts_nhc ! MTS loop
      do j=1,mdr%th%n_ys ! Yoshida-Suzuki loop
        ! Update thermostat velocities
        do k=mdr%th%n_nhc,1,-1
          if (k==mdr%th%n_nhc) then
            ! call mdr%th%update_G_k(k, mdr%baro%ke_box) ! should this be here?
            call mdr%th%propagate_v_eta_k_1(k, dt, quarter)
          else
            call mdr%th%propagate_v_eta_k_2(k, dt, eighth)
            call mdr%th%propagate_v_eta_k_1(k, dt, quarter)
            call mdr%th%update_G_k(k, mdr%baro%ke_box) ! should this be here?
            call mdr%th%propagate_v_eta_k_2(k, dt, eighth)
          end if
        end do
        ! Update box velocities
        call mdr%iso_baro%propagate_v_eps_2(dt, eighth, mdr%th%v_eta(1))
        call mdr%iso_baro%propagate_v_eps_1(dt, quarter)
        call mdr%iso_baro%propagate_v_eps_2(dt, eighth, mdr%th%v_eta(1))
        ! Update Particle velocities
        vscale = exp(-half*dt*(mdr%th%v_eta(1) + &
                 (one + three/mdr%ndof)*mdr%iso_baro%v_eps))
        v_sfac = v_sfac*vscale
        mdr%ke = mdr%th%ke_system*vscale**2
        mdr%th%ke_system = mdr%ke
        ! Update the box forces
        call mdr%iso_baro%update_G_eps(mdr%ke, mdr%virial, mdr%V, mdr%T_ext, &
                                       mdr%th%Q(1), G_nhc_1)
        mdr%th%G_nhc(1) = G_nhc_1
        ! Update thermostat positions
        do k=1,mdr%th%n_nhc
          call mdr%th%propagate_eta_k(k, dt, half)
        end do
        ! Get total ke
        call mdr%get_kinetic_energy
        ! Update the box velocities
        call mdr%iso_baro%propagate_v_eps_2(dt, eighth, mdr%th%v_eta(1))
        call mdr%iso_baro%propagate_v_eps_1(dt, quarter)
        call mdr%iso_baro%propagate_v_eps_2(dt, eighth, mdr%th%v_eta(1))
        ! Update the box ke
        call mdr%iso_baro%get_box_ke_iso
        ! Update the thermostat velocities
        do k=1,mdr%th%n_nhc
          if (k<mdr%th%n_nhc) then
            call mdr%th%propagate_v_eta_k_2(k, dt, eighth)
            call mdr%th%update_G_k(k, mdr%baro%ke_box)
            call mdr%th%propagate_v_eta_k_1(k, dt, quarter)
            call mdr%th%propagate_v_eta_k_2(k, dt, eighth)
          else
            call mdr%th%update_G_k(k, mdr%baro%ke_box)
            call mdr%th%propagate_v_eta_k_1(k, dt, quarter)
          end if
        end do
      end do ! MTS loop
    end do ! Yoshida-Suzuki loop

    ! scale the particle velocities
    mdr%v_t_dt = mdr%v_t*v_sfac

  end subroutine propagate_iso_npt_mttk

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

    write(iunit,'("frame ",i8)') step
    write(iunit,'(a)') "cell_vectors"
    do i=1,3
      write(iunit,'(3f12.6)') mdr%p_t%h(i,:)
    end do
    write(iunit,'(a)') "end cell_vectors"
    write(iunit,'(a)') "stress_tensor"
    do i=1,3
      write(iunit,'(3f12.6)') mdr%stress(i,:)
    end do
    write(iunit,'(a)') "end stress_tensor"
    write(iunit,'(a)') "positions"
    call mdr%dump_atom_arr(iunit, mdr%p_t%rcart)
    write(iunit,'(a)') "end positions"
    write(iunit,'(a)') "velocities"
    call mdr%dump_atom_arr(iunit, mdr%v_t)
    write(iunit,'(a)') "end velocities"
    write(iunit,'(a)') "forces"
    call mdr%dump_atom_arr(iunit, mdr%f)
    write(iunit,'(a)') "end forces"
    write(iunit,'(a)') "end frame"

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
        write(iunit,'(a10,5a16)') "step", "pe", "ke", "total", "T", "P"
      case ('nvt')
        if (mdr%thermo_type == 'nhc') then
          write(iunit,'(a10,6a16)') "step", "pe", "ke", "nhc", "total", "T", "P"
        else
          write(iunit,'(a10,5a16)') "step", "pe", "ke", "total", "T", "P"
        end if
      case ('npt')
        if (mdr%thermo_type == 'nhc') then
          write(iunit,'(a10,8a16)') "step", "pe", "ke", "nhc", "box", "pV", "total", "T", "P"
        end if
      end select
    end if
    select case (mdr%ensemble)
    case ('nve')
      write(iunit,'(i10,5e16.6)') step, mdr%pe, mdr%ke, mdr%e_nve, &
                                  mdr%T_int, mdr%P_int
    case ('nvt')
      if (mdr%thermo_type == 'nhc') then
        write(iunit,'(i10,6e16.6)') step, mdr%pe, mdr%ke, mdr%e_nhc, &
                                    mdr%e_nvt, mdr%T_int, mdr%P_int
      else
        write(iunit,'(i10,5e16.6)') step, mdr%pe, mdr%ke, mdr%pe+mdr%ke, &
                                    mdr%T_int, mdr%P_int
      end if
    case ('npt')
      if (mdr%thermo_type == 'nhc') then
        write(iunit,'(i10,8e16.6)') step, mdr%pe, mdr%ke, mdr%e_nhc, &
                                    mdr%e_nvt, mdr%ke_box, mdr%pv, &
                                    mdr%T_int, mdr%P_int
      end if
    end select
  end subroutine stat_dump

end module md_module
