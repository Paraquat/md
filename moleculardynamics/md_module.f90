module md_module

use datatypes
use constants
use vector
use cell
use pairpotential
use rng
use md_model
use md_control
use pairdist

implicit none

type type_md
  type(type_cell)           :: p_t
  type(type_model)          :: mdl
  type(type_pairpotential)  :: pp
  type(type_thermostat)     :: th
  type(type_barostat)       :: baro
  type(type_pairdist)       :: pd
  character(3)    :: ensemble
  integer         :: nsteps, step, nspec, nat, ndof
  integer         :: iprint, dump_freq, nrestart, cp_freq
  integer         :: n_nhc, n_mts, n_ys
  logical         :: shift, remove_com_v, dump, restart, append, rdf
  character(40)   :: position_file, dump_file, stat_file, cp_file, gr_file, &
                     thermo_type, baro_type, init_distr, units, pbc_method
  integer, allocatable, dimension(:)          :: species
  real(double), allocatable, dimension(:)     :: mass
  real(double), allocatable, dimension(:,:)   :: v_t
  real(double), allocatable, dimension(:,:)   :: f
  real(double), dimension(3,3)                :: stress, kinetic_stress, &
                                                 static_stress
  real(double)   :: T_ext, dt, k_B_md, P_ext, &
                    V, PV, h_prime, tau_T, tau_P, bulkmod_est

  contains
    procedure, public  :: init_md
    procedure, private :: init_velocities
    procedure, private :: remove_com_velocity
    procedure, private :: get_force_and_energy
    procedure, private :: get_kinetic_energy_and_stress
    procedure, private :: get_static_stress
    procedure, private :: get_pv
    procedure, private :: vVerlet_v
    procedure, private :: vVerlet_r
    procedure, public  :: md_run
    procedure, private :: get_cons_qty
    procedure, public  :: fire
    procedure, private :: dump_atom_arr
    procedure, private :: md_dump
    procedure, private :: write_checkpoint
    procedure, private :: read_checkpoint

end type type_md

contains

  subroutine init_md(mdr, inp, init_cell, pp)

    use input_module,   only: md_input

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    type(md_input), intent(in)      :: inp
    type(type_cell), intent(in)     :: init_cell
    type(type_pairpotential), intent(in)  :: pp

    ! local variables
    integer                         :: i, j
    real(double)                    :: la, lb, lc, cut_new
    integer                         :: cp_unit
    character(40)                   :: dumpfile_prefix, statfile_prefix, &
                                       posfile_prefix, suffix

    write(*,'(a)') "Initialising MD"
    mdr%dump_file = "Frames"
    mdr%stat_file = "Stats"
    mdr%position_file = "trajectory.xsf"
    mdr%gr_file = "gr.dat"

    mdr%nrestart = 0
    mdr%restart = inp%restart
    mdr%ensemble = inp%ensemble
    mdr%nsteps = inp%nsteps
    mdr%iprint = inp%iprint
    mdr%p_t = init_cell
    mdr%nspec = mdr%p_t%nspec
    mdr%nat = mdr%p_t%nat
    mdr%dt = inp%dt
    mdr%T_ext = inp%T_ext
    mdr%shift = inp%shift
    mdr%init_distr = inp%v_distr
    mdr%remove_com_v = inp%comv
    mdr%thermo_type = inp%thermo_type
    mdr%baro_type = inp%baro_type
    mdr%n_nhc = inp%n_nhc
    mdr%P_ext = inp%P_ext
    mdr%tau_T = inp%tau_T
    mdr%tau_P = inp%tau_P
    mdr%n_mts = inp%n_mts
    mdr%n_ys = inp%n_ys
    mdr%rdf = inp%rdf

    if (inp%cart .eqv. .false.) call mdr%p_t%cell_frac2cart
    mdr%pp = pp

    mdr%pbc_method = init_cell%pbc_method
    mdr%p_t%pbc_method = mdr%pbc_method
    mdr%V = mdr%p_t%volume()

    ! degrees of freedom
    mdr%ndof = 3*mdr%nat
    if (mdr%remove_com_v .eqv. .true.) mdr%ndof = mdr%ndof-3

    call mdr%mdl%init_model(inp, mdr%p_t, mdr%v_t, mdr%f, mdr%th, &
                            mdr%baro, mdr%ndof)

    ! Dumping information
    mdr%dump = .true.
    mdr%dump_freq = inp%dump_freq
    mdr%units = 'reduced-lj'

    select case (mdr%units)
    case ('reduced-lj')
      mdr%k_B_md = one
    case default
      mdr%k_B_md = k_B
    end select

    allocate(mdr%species(mdr%nat))
    allocate(mdr%mass(mdr%nat))
    allocate(mdr%v_t(mdr%nat,3))
    allocate(mdr%f(mdr%nat,3))
    mdr%species = init_cell%spec_int
    do i=1,mdr%nat
      mdr%mass(i) = mdr%p_t%mass(mdr%p_t%spec_int(i))
    end do

    write(*,'(a)') "Simulation parameters:"
    write(*,'("Ensemble               ",a8)') mdr%ensemble
    write(*,'("Number of steps        ",i12)') mdr%nsteps
    write(*,'("Time step              ",f8.4)') mdr%dt
    write(*,'("Remove COM velocity    ",l8)') mdr%remove_com_v
    write(*,*) 

    ! thermostat and barostat (just for temperature and pressure if type==none)
    call mdr%th%init_thermostat(inp, mdr%ndof, mdr%mdl%E_k)
    call mdr%baro%init_barostat(inp, mdr%ndof, mdr%mdl%E_k, init_cell)

    ! Adjust the pair potential cutoff if necessary
    cut_new = maxval(mdr%pp%r_cut) + one
    if (mdr%pbc_method == 'mic') then
      call mdr%p_t%cut_ortho(cut_new)
    else if (mdr%pbc_method == 'frac') then
      call mdr%p_t%cut_frac(cut_new)
    end if
    if (cut_new < maxval(mdr%pp%r_cut)) then
        write(*,'(2x,"Cutoff in pp.in is larger than shortest box side!")')
        write(*,'(2x,"Adjusting r_cut to ",f8.4)') cut_new
        write(*,*)
        do i=1,mdr%nspec
          do j=1,mdr%nspec
            if (mdr%pp%r_cut(i,j) > cut_new) mdr%pp%r_cut(i,j) = cut_new
          end do
        end do
    end if
    write(*,'(a)') "Species information:"
    write(*,'(2x,3a8,a10)') "label", "species", "count", "mass"
    do i=1,mdr%nspec
      write(*,'(2x,i8,a8,i8,f10.4)') mdr%p_t%spec_int(i), mdr%p_t%spec(i), &
                               mdr%p_t%spec_count(i), mdr%p_t%mass(i)
    end do
    write(*,*)

    ! radial distribution function
    if (mdr%rdf) then
      call mdr%pd%init_pd(mdr%p_t, inp%rdfcut, inp%gwidth, inp%dr, inp%rmin)
    end if

    if (mdr%restart .eqv. .false.) then
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
      write(*,*)

      ! initialise velocity
      call mdr%init_velocities(mdr%th, mdr%init_distr)
      write(*,*)
      write(*,'(2x,a)') "Initial velocities:"
      do i=1,mdr%nat
        write(*,'(4x,2i6,3f14.8)') i, mdr%species(i), mdr%v_t(i,:)
      end do
      write(*,*)

      call mdr%get_force_and_energy
      write(*,*)
      write(*,'(2x,a)') "Initial forces:"
      do i=1,mdr%nat
        write(*,'(4x,2i6,3f14.8)') i, mdr%species(i), mdr%f(i,:)
      end do
      write(*,*)
      call mdr%baro%update_stress(mdr%static_stress, mdr%kinetic_stress)
      call mdr%baro%get_stress_and_pressure
      call mdr%get_pv
      write(*,*)
    end if

    mdr%step = 1
    ! Restart from checkpoint
    mdr%cp_file = "checkpoint"
    if (mdr%restart .eqv. .true.) then
      write(*,'("Reading checkpoint file ", a)') mdr%cp_file
      cp_unit   = 104
      call mdr%read_checkpoint(cp_unit)
      write(*,'("Restarting from step    ", i8)') mdr%step
      write(*, '("This is restart number ", i8)') mdr%nrestart
      call mdr%get_force_and_energy
      call mdr%get_kinetic_energy_and_stress
      mdr%V = mdr%p_t%volume()
      if (mdr%thermo_type == 'nhc')  then
        mdr%th%ke_atoms = mdr%mdl%E_k
        call mdr%th%get_nhc_energy
      end if
      if (mdr%baro_type == 'mttk' .or. mdr%baro_type == 'iso-mttk' .or. &
          mdr%baro_type == 'ortho-mttk') then
        call mdr%baro%get_box_ke
        mdr%mdl%E_box = mdr%baro%ke_box
      end if
      call mdr%get_cons_qty
      call mdr%baro%update_stress(mdr%static_stress, mdr%kinetic_stress)
      call mdr%baro%get_stress_and_pressure
      call mdr%th%get_temperature(mdr%mdl%E_k)
      call mdr%get_pv
    end if
    write(*,*)

  end subroutine init_md

  subroutine init_velocities(mdr, th, distr)

    ! passed variables
    class(type_md), intent(inout)         :: mdr
    type(type_thermostat), intent(inout)  :: th
    character(40), intent(in)             :: distr

    ! local variables
    real(double)                    :: sfac, sigma, mu, z1, z2
    integer                         :: i, j, k
    real(double), allocatable, dimension(:) :: rn

    if (mdr%iprint > 1) write(*,'(a)') "init_velocities"
    write(*,'(a)') "Initialising velocities:"
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

    write(*,'("Velocity distribution  ",a20)') mdr%init_distr
    write(*,'("Initial/external T     ",f8.4)') mdr%T_ext
    call mdr%remove_com_velocity(mdr%remove_com_v)
    write(*,'(a)') "Before rescaling:"
    call mdr%get_kinetic_energy_and_stress
    call th%get_temperature(mdr%mdl%E_k)
    sfac = sqrt(mdr%T_ext/mdr%th%T_int)
    write(*,*)
    mdr%v_t = sfac*mdr%v_t ! Scale velocities according to T_ext
    write(*,'(a)') "After rescaling:"
    call mdr%get_kinetic_energy_and_stress
    call th%get_temperature(mdr%mdl%E_k)
    write(*,*)

  end subroutine init_velocities

  ! Reomve centre of mass velocity, compute velocity and sum of velocities
  ! squared
  subroutine remove_com_velocity(mdr, remove_com_v)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    logical, intent(in)             :: remove_com_v

    ! local variables
    integer                         :: i
    real(double), dimension(3)      :: sumv

    if (mdr%iprint > 1) write(*,'(2x,a)') "remove_com_velocity"
    ! update the COM velocity
    sumv = zero
    do i=1,3
      sumv(i) = sum(mdr%v_t(:,i))
    end do
    sumv = sumv/mdr%nat

    ! remove COM velocity
    if (remove_com_v .eqv. .true.) then
      do i=1,mdr%nat
        mdr%v_t(i,:) = mdr%v_t(i,:) - sumv
      end do
    end if

  end subroutine remove_com_velocity

  ! wrapper to get force and energy
  subroutine get_force_and_energy(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    if (mdr%iprint > 1) write(*,'(2x,a)') "get_force_and_energy"
    call mdr%pp%get_force_and_energy(mdr%p_t, mdr%mdl%E_p, mdr%f)
    write(*,'("  Potential energy = ",e16.8)') mdr%mdl%E_p
    call mdr%get_static_stress

  end subroutine get_force_and_energy

  ! Compute the kinetic energy
  subroutine get_kinetic_energy_and_stress(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer   :: i

    mdr%mdl%E_k = zero
    mdr%kinetic_stress = zero
    do i=1,mdr%nat
      mdr%mdl%E_k = mdr%mdl%E_k + mdr%mass(i)*sum(mdr%v_t(i,:)**2)
      mdr%kinetic_stress = mdr%kinetic_stress + &
                      mdr%mass(i)*tensor_product(mdr%v_t(i,:), mdr%v_t(i,:))
    end do
    mdr%mdl%E_k = mdr%mdl%E_k/two
    write(*,'("  Kinetic energy   = ",e16.8)') mdr%mdl%E_k

  end subroutine get_kinetic_energy_and_stress

  ! Compute the static stress from the forces
  subroutine get_static_stress(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    ! local variables
    integer   :: i

    if (mdr%iprint > 1) write(*,'(2x,a)') "get_static_stress"
    mdr%static_stress = zero

    do i=1,mdr%nat
      mdr%static_stress = mdr%static_stress + &
                          tensor_product(mdr%f(i,:), mdr%p_t%rcart(i,:))
    end do

  end subroutine get_static_stress

  ! Compute the pV term for the enthalpy
  subroutine get_pv(mdr)

    ! passed variables
    class(type_md), intent(inout)   :: mdr

    mdr%V = mdr%p_t%volume()
    mdr%mdl%PV = mdr%V*mdr%P_ext ! Should use the target pressure according to MTTK
    write(*,'("  volume           = ",e16.8)') mdr%V
    write(*,'("  PV               = ",e16.8)') mdr%PV

  end subroutine get_pv

  ! Velocity Verlet dt/2 step for velocities
  subroutine vVerlet_v(mdr, dt, dtfac)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    real(double), intent(in)        :: dt, dtfac

    ! local variables
    integer                         :: i

    if (mdr%iprint > 1) write(*,'(4x,a)') "Velocity Verlet velocity update"
    do i=1,mdr%nat
      mdr%v_t(i,:) = mdr%v_t(i,:) + dt*dtfac*mdr%f(i,:)/mdr%mass(i)
    end do

  end subroutine vVerlet_v

  ! Velocity Verlet dt step for positions
  subroutine vVerlet_r(mdr, dt, dtfac)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    real(double), intent(in)        :: dt, dtfac

    if (mdr%iprint > 1) write(*,'(4x,a)') "Velocity Verlet position update"
    mdr%p_t%rcart = mdr%p_t%rcart + dt*dtfac*mdr%v_t
    ! wrap the updated positions back into the unit cell
    if (mdr%iprint > 1) write(*,'(4x,a)') "Wrapping positions into unit cell"
    call mdr%p_t%wrap_positions_cart

  end subroutine vVerlet_r

  ! The main MD loop
  subroutine md_run(mdr, s_start, s_end)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    integer, intent(in)             :: s_start, s_end

    ! local variables
    integer       :: s, traj_unit, dump_unit, stat_unit, cp_unit, &
                     debug_unit_1, debug_unit_2, i
    logical       :: d

    write(*,'(a)') "Starting main MD loop"
    traj_unit = 101
    dump_unit = 102
    stat_unit = 103
    cp_unit   = 104
    debug_unit_1 = 109
    debug_unit_1 = 110

    if (mdr%restart) then
      open(unit=traj_unit,file=mdr%position_file,status='old',position='append')
      open(unit=dump_unit,file=mdr%dump_file, status='old',position='append')
      open(unit=stat_unit,file=mdr%stat_file, status='old',position='append')
      if (mdr%iprint > 2) then
        open(debug_unit_1,file='baro_state', status='old',position='append')
        open(debug_unit_2,file='thermo_state', status='old',position='append')
      end if
    else
      open(unit=traj_unit,file=mdr%position_file,status='replace')
      open(unit=dump_unit,file=mdr%dump_file,status='replace')
      open(unit=stat_unit,file=mdr%stat_file,status='replace')
      if (mdr%iprint > 2) then
        open(debug_unit_1,file='baro_state',status='replace')
        open(debug_unit_2,file='thermo_state',status='replace')
      end if
    end if

    if (mdr%restart .eqv. .false.) then
      ! Only dump the initial configuration before the first step
      s = 0
      call mdr%p_t%write_xsf(traj_unit, .true., 1, (s_end/mdr%dump_freq)+1)
      call mdr%mdl%stat_dump(stat_unit, s)
      if (mdr%dump .eqv. .true.) then
        call mdr%md_dump(dump_unit, 0)
        if (mdr%rdf) then
          call mdr%pd%update_rdist(mdr%p_t)
        end if 
      end if
    else
      mdr%step = mdr%step + 1
    end if

    if (mdr%iprint > 2) then
      if (mdr%ensemble(3:3) == 't') then
        if (mdr%thermo_type == 'nhc') then
          call mdr%th%dump_thermo_state(s, debug_unit_2)
        end if
      end if
      if (mdr%ensemble(2:2) == 'p') then
        call mdr%baro%dump_baro_state(s, debug_unit_1)
      end if
    end if

    do s=s_start,s_end ! Main MD loop
      mdr%step = s
      d = .false.
      if (mod(s,mdr%dump_freq) == 0) d = .true.
      write(*,*)
      write(*,'("MD step ",i10," of ",i10)')  s, s_end

      ! Propagate thermostat/barostat
      select case (mdr%ensemble)
      case('nvt')
        if (mdr%thermo_type == 'nhc') then
          call mdr%th%propagate_nvt_nhc(mdr%dt, mdr%v_t)
        else if (mdr%thermo_type == 'berendsen') then
          call mdr%th%get_berendsen_thermo_sf
        end if
      case('npt')
        if (mdr%baro_type == 'mttk' .or. mdr%baro_type == 'iso-mttk' .or. &
            mdr%baro_type == 'ortho-mttk') then
          call mdr%baro%propagate_npt_mttk(mdr%th, mdr%dt, mdr%mass, mdr%v_t, mdr%mdl%E_k)
        end if
        if (mdr%thermo_type == 'berendsen') then
          call mdr%th%get_berendsen_thermo_sf
        end if
      end select
      ! Velocity Verlet dt/2 step
      call mdr%vVerlet_v(mdr%dt, half)

      ! Propagate atoms (and box for NPT)
      if (mdr%ensemble(2:2) == 'p') then ! constant pressure
        select case(mdr%baro_type)
        case('mttk')
          ! call mdr%baro%diag_vbox(mdr%th%v_eta(1))
          call mdr%baro%get_Ie(mdr%dt, one)
          call mdr%baro%get_Is(mdr%dt, half)
          call mdr%baro%propagate_r_ions(mdr%dt, one, mdr%p_t%rcart, mdr%v_t)
          call mdr%baro%propagate_h(mdr%dt, one, mdr%p_t%h)
          mdr%p_t%V = mdr%p_t%volume()
          mdr%V = mdr%p_t%V
        case('ortho-mttk')
          call mdr%baro%propagate_r_ions(mdr%dt, half, mdr%v_t, &
                                            mdr%p_t%rcart)
          call mdr%baro%propagate_box(mdr%dt, one, mdr%p_t%h, mdr%p_t%V)
          mdr%p_t%V = mdr%p_t%volume()
          mdr%V = mdr%p_t%V
        case('iso-mttk')
          call mdr%baro%propagate_r_ions(mdr%dt, half, mdr%v_t, &
                                            mdr%p_t%rcart)
          call mdr%baro%propagate_eps_lin(mdr%dt, one)
          call mdr%baro%propagate_box(mdr%dt, one, mdr%p_t%h, mdr%p_t%V)
          mdr%p_t%V = mdr%p_t%volume()
          mdr%V = mdr%p_t%V
        case('berendsen')
          call mdr%baro%get_berendsen_baro_sf
          call mdr%vVerlet_r(mdr%dt, one) ! velocity Verlet r update
        end select
      else
        call mdr%vVerlet_r(mdr%dt, one) ! velocity Verlet r update
      end if

      ! dump configuration to .xsf file
      if (d .eqv. .true.) then
        call mdr%p_t%write_xsf(traj_unit, .true., s+1, s_end+1)
      end if
      call mdr%get_force_and_energy
      if (mdr%iprint > 3) then
        write(*,'(2x,a)') "Forces on atoms:"
        do i=1,mdr%nat
          write(*,'(4x,2i6,3e16.8)') i, mdr%species(i), mdr%f(i,:)
        end do
      end if
      call mdr%baro%update_stress(mdr%static_stress, mdr%kinetic_stress)
      call mdr%baro%get_stress_and_pressure
      call mdr%th%get_temperature(mdr%mdl%E_k)

      ! Velocity Verlet dt/2 step
      call mdr%vVerlet_v(mdr%dt, half)
      ! Propagate thermostat/barostat
      select case (mdr%ensemble)
      case('nvt')
        if (mdr%thermo_type == 'nhc') then
          call mdr%th%propagate_nvt_nhc(mdr%dt, mdr%v_t)
          call mdr%th%get_nhc_energy
        else if (mdr%thermo_type == 'berendsen') then
          call mdr%th%berendsen_thermo_propagate(mdr%v_t)
        end if
      case('npt')
        if (mdr%baro_type=='mttk' .or. mdr%baro_type=='iso-mttk' .or. &
            mdr%baro_type=='ortho-mttk') then
          call mdr%baro%propagate_npt_mttk(mdr%th, mdr%dt, mdr%mass, mdr%v_t, mdr%mdl%E_k)
        end if
        if (mdr%baro_type == 'berendsen') then
          ! Propagate the Berendsen barostat after the velocity update so that
          ! rescaling does not affect velocities
          call mdr%baro%berendsen_baro_propagate(mdr%p_t%rcart, mdr%p_t%h)
          call mdr%th%berendsen_thermo_propagate(mdr%v_t)
        end if
        mdr%baro%h = mdr%p_t%h
        mdr%baro%V = mdr%V
        if (mdr%thermo_type == 'nhc') call mdr%th%get_nhc_energy
      end select

      ! Update arrays and thermodyanmics quantities
      call mdr%get_kinetic_energy_and_stress
      call mdr%remove_com_velocity(mdr%remove_com_v)
      ! call mdr%p_t%wrap_positions_cart

      write(*,'(2x,a,i10)') "End of MD step", s

      call mdr%get_cons_qty ! Compute the conserved quantity
      call mdr%baro%update_stress(mdr%static_stress, mdr%kinetic_stress)
      call mdr%baro%get_stress_and_pressure
      call mdr%th%get_temperature(mdr%mdl%E_k)
      if (mdr%ensemble(3:3) == 't') mdr%th%ke_atoms = mdr%mdl%E_k
      if (mdr%ensemble(3:3) == 'p') mdr%baro%ke_atoms = mdr%mdl%E_k
      if (mdr%ensemble(2:2) == 'p') call mdr%get_pv

      ! Write the checkpoint
      if (mod(s,mdr%cp_freq) == 0) call mdr%write_checkpoint(cp_unit, s)

      ! Dump the statistics
      call mdr%mdl%stat_dump(stat_unit, s)
      if (mdr%dump .eqv. .true.) then
        if (d .eqv. .true.) then
          call mdr%md_dump(dump_unit, s)
          call flush(traj_unit)
          call flush(dump_unit)
          call flush(stat_unit)
        end if
      end if
      if (mdr%iprint > 2) then
        if (mdr%ensemble(3:3) == 't') then
          if (mdr%thermo_type == 'nhc') then
            call mdr%th%dump_thermo_state(s, debug_unit_2)
          end if
        end if
        if (mdr%ensemble(2:2) == 'p') then
          call mdr%baro%dump_baro_state(s, debug_unit_1)
        end if
      end if
    end do ! Main MD loop
    write(*,'(a)')  "Completed MD run"
    call mdr%p_t%cell_cart2frac
    call mdr%p_t%write_cell("cell.out")
    if (mdr%rdf) then
      call mdr%pd%norm_rdist
      call mdr%pd%write_gr(mdr%gr_file)
    end if
    close(traj_unit)
    close(dump_unit)
    close(stat_unit)
    if (mdr%iprint > 2) then
      close(debug_unit_1)
      close(debug_unit_2)
    end if

  end subroutine md_run

  ! Get the conserved quantity for the dynamics
  subroutine get_cons_qty(mdr)

    ! passed variables
    class(type_md), intent(inout)             :: mdr

    write(*,*)
    write(*,'("  Components of conserved quantity:")')
    ! For the NVE ensemble, the conserved quantity is just the total energy
    mdr%mdl%H_prime = mdr%mdl%E_k + mdr%mdl%E_p
    write(*,'("  Kinetic energy   = ",e16.8)') mdr%mdl%E_k
    write(*,'("  Potential energy = ",e16.8)') mdr%mdl%E_p
    write(*,'("  Total energy     = ",e16.8)') mdr%mdl%H_prime
    ! For NVT (NHC), add the NHC energy
    if (mdr%thermo_type == 'nhc') then
      mdr%mdl%H_prime = mdr%mdl%H_prime + mdr%mdl%E_nhc
      write(*,'("  NHC energy       = ",e16.8)') mdr%mdl%E_nhc
    end if
    ! For NPT, add the NHC energy, the box kinetic energy and the PV term
    if (mdr%ensemble(2:2) == 'p') then
      mdr%mdl%H_prime = mdr%mdl%H_prime + mdr%mdl%E_box
      mdr%mdl%H_prime = mdr%mdl%H_prime + mdr%mdl%PV
      write(*,'("  Box energy       = ",e16.8)') mdr%mdl%E_box
      write(*,'("  PV               = ",e16.8)') mdr%mdl%PV
    end if
    write(*,'(2x,a3,a7,a6,a3,e16.8)') mdr%ensemble, ' energy', '       ', ' = ', mdr%mdl%H_prime
    write(*,*)

  end subroutine get_cons_qty

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

    fmt = "(2i5,3f20.12)"
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

  ! create a checkpoint for restarting
  subroutine write_checkpoint(mdr, iunit, step)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    integer, intent(in)             :: iunit
    integer, intent(in)             :: step

    ! local variables
    integer                         :: i

    open(iunit, file=mdr%cp_file, status='replace')
    write(iunit,*) mdr%nrestart + 1
    write(iunit,*) step
    do i=1,3
      write(iunit,*) mdr%p_t%h(i,:)
    end do
    do i=1,mdr%nat
      write(iunit,*) mdr%species(i), mdr%p_t%rcart(i,:)
    end do
    do i=1,mdr%nat
      write(iunit,*) mdr%v_t(i,:)
    end do
    do i=1,mdr%nat
      write(iunit,*) mdr%f(i,:)
    end do
    if (mdr%thermo_type == 'nhc') then
      write(iunit,*) mdr%n_nhc      
      write(iunit,*) mdr%th%eta
      write(iunit,*) mdr%th%v_eta
      write(iunit,*) mdr%th%G_nhc
    end if
    if (mdr%baro_type == 'iso-mttk') then
      write(iunit,*) mdr%baro%V_ref
      write(iunit,*) mdr%baro%V
      write(iunit,*) mdr%baro%eps
      write(iunit,*) mdr%baro%v_eps
      write(iunit,*) mdr%baro%G_eps
    else if (mdr%baro_type == 'ortho-mtk') then
      write(iunit,*) mdr%baro%Q_ref
      write(iunit,*) mdr%baro%Q
      write(iunit,*) mdr%baro%v_Q
      write(iunit,*) mdr%baro%G_Q
    else if (mdr%baro_type == 'mttk') then
      do i=1,3
        write(iunit,*) mdr%baro%h_ref(i,:)
      end do
      do i=1,3
        write(iunit,*) mdr%baro%v_h(i,:)
      end do
      do i=1,3
        write(iunit,*) mdr%baro%G_h(i,:)
      end do
    end if
    close(iunit)

  end subroutine write_checkpoint

  ! read a checkpoint file
  subroutine read_checkpoint(mdr, iunit)

    ! passed variables
    class(type_md), intent(inout)   :: mdr
    integer, intent(in)             :: iunit

    ! local variables
    integer                         :: i

    open(iunit, file=mdr%cp_file, status='old')
    read(iunit,*) mdr%nrestart
    read(iunit,*) mdr%step
    do i=1,3
      read(iunit,*) mdr%p_t%h(i,:)
    end do
    do i=1,mdr%nat
      read(iunit,*) mdr%species(i), mdr%p_t%rcart(i,:)
    end do
    do i=1,mdr%nat
      read(iunit,*) mdr%v_t(i,:)
    end do
    do i=1,mdr%nat
      read(iunit,*) mdr%f(i,:)
    end do
    if (mdr%thermo_type == 'nhc') then
      read(iunit,*) mdr%n_nhc      
      read(iunit,*) mdr%th%eta
      read(iunit,*) mdr%th%v_eta
      read(iunit,*) mdr%th%G_nhc
    end if
    if (mdr%baro_type == 'iso-mttk') then
      read(iunit,*) mdr%baro%V_ref
      read(iunit,*) mdr%baro%V
      read(iunit,*) mdr%baro%eps
      read(iunit,*) mdr%baro%v_eps
      read(iunit,*) mdr%baro%G_eps
    else if (mdr%baro_type == 'ortho-mtk') then
      read(iunit,*) mdr%baro%Q_ref
      read(iunit,*) mdr%baro%Q
      read(iunit,*) mdr%baro%v_Q
      read(iunit,*) mdr%baro%G_Q
    else if (mdr%baro_type == 'mttk') then
      do i=1,3
        read(iunit,*) mdr%baro%h_ref(i,:)
      end do
      do i=1,3
        read(iunit,*) mdr%baro%v_h(i,:)
      end do
      do i=1,3
        read(iunit,*) mdr%baro%G_h(i,:)
      end do
    end if
    close(iunit)
    
  end subroutine read_checkpoint

end module md_module
