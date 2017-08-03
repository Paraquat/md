module thermostat_module

use datatypes
use constants
use rng

implicit none

type type_thermostat

  real(double)        :: T_ext, k_B_md, ke_system
  integer             :: nat, ndof, iprint
  integer             :: n_ys    ! Yoshida-Suzuki order
  integer             :: n_nhc   ! number of nh chains
  integer             :: mts_nhc ! number of NHC loops per time step
  real(double)        :: e_nhc   ! NHC contribution to the conserved term
  real(double)        :: tau_T   ! temperature coupling time period
  real(double)        :: dt      ! time step
  real(double)        :: lambda  ! Berendsen scaling factor
  character(40)       :: th_type
  real(double), dimension(:), allocatable   :: eta   ! thermostat k position
  real(double), dimension(:), allocatable   :: v_eta ! thermostat k velocity
  real(double), dimension(:), allocatable   :: G_nhc ! force on thermostat k
  real(double), dimension(:), allocatable   :: Q     ! thermostat k mass

contains

  procedure :: init_thermostat
  procedure :: velocity_rescale
  procedure :: get_berendsen_thermo_sf
  procedure :: berendsen_thermo_propagate
  procedure :: init_nhc
  procedure :: update_G_k
  procedure :: propagate_eta_k
  procedure :: propagate_v_eta_k_1
  procedure :: propagate_v_eta_k_2
  procedure :: propagate_nvt_nhc
  procedure :: get_nhc_ke
  procedure :: dump_thermo_state

end type type_thermostat

contains

  ! Initialise the thermostat
  subroutine init_thermostat(th, th_type, nat, ndof, dt, T_ext, tau_T, n_nhc, &
                             iprint, ke_system)

    ! passed variables
    class(type_thermostat), intent(inout)     :: th
    character(40), intent(in)                 :: th_type
    integer, intent(in)                       :: nat, ndof, iprint, n_nhc
    real(double), intent(in)                  :: T_ext, ke_system, dt, tau_T

    th%th_type = th_type
    th%nat = nat
    th%ndof = ndof
    th%dt = dt
    th%T_ext = T_ext
    th%iprint = iprint
    th%tau_T = tau_T
    th%k_B_md = one
    th%n_nhc = n_nhc
    th%mts_nhc = 1
    th%n_ys = 1
    th%ke_system = ke_system
    if (th_type == 'nhc') call th%init_nhc(th%n_nhc)

  end subroutine init_thermostat

  ! Isokinetic thermostat: maintain temperature by rescaling velocities by
  ! sqrt(T_int/T_ext) every tau_T steps
  subroutine velocity_rescale(th, T_int, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), intent(in)                    :: T_int
    real(double), dimension(:,:), intent(inout) :: v

    ! local variables
    real(double)  :: lambda

    lambda = sqrt(th%T_ext/T_int)
    v = v*lambda
    if (th%iprint == 0) write(*,'(4x,"Velocity rescale, lambda = ",f8.4)') lambda
  end subroutine velocity_rescale

  subroutine get_berendsen_thermo_sf(th, T_int)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), intent(in)                    :: T_int

    th%lambda = sqrt(one + (th%dt/th%tau_T)*(th%T_ext/T_int - one))
    if (th%iprint == 0) then
      write(*,'(4x,"Weak coupling, thermostat lambda = ",f8.4)') th%lambda
    end if

  end subroutine get_berendsen_thermo_sf

  ! Berendsen weak coupling thermostat
  subroutine berendsen_thermo_propagate(th, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v

    v = v*th%lambda

  end subroutine berendsen_thermo_propagate

  ! Initialise the NHC with n_nhc heat baths
  subroutine init_nhc(th, n_nhc)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: n_nhc

    th%n_nhc = n_nhc
    allocate(th%eta(n_nhc))
    allocate(th%v_eta(n_nhc))
    allocate(th%G_nhc(n_nhc))
    allocate(th%Q(n_nhc))

    th%eta = zero
    if (th%n_nhc > 0) then
      th%Q = one
      th%v_eta = sqrt(two*th%T_ext/th%Q(1))
    else
      th%v_eta = zero
    end if
    th%G_nhc = zero
    call th%get_nhc_ke

  end subroutine init_nhc

  ! get the force on thermostat k
  subroutine update_G_k(th, k, box_ke)

    ! passed variables
    class(type_thermostat), intent(inout)   :: th
    integer, intent(in)                     :: k
    real(double), intent(in)                :: box_ke ! for P-T coupling

    if (k == 1) then
      th%G_nhc(k) = 2*th%ke_system - th%ndof*th%k_B_md*th%T_ext + box_ke
    else
      th%G_nhc(k) = (th%Q(k-1)*th%v_eta(k-1)**2 - th%k_B_md*th%T_ext)
    end if
    th%G_nhc(k) = th%G_nhc(k)/th%Q(k)

  end subroutine update_G_k

  ! propagate eta (this is the same for all k)
  subroutine propagate_eta_k(th, k, dt, dtfac)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer                               :: k
    real(double)                          :: dt, dtfac

    if (th%iprint == 0) write(*,'(6x,"NHC: propagating eta,                k = ",i2)') k
    th%eta(k) = th%eta(k) + dtfac*dt*th%v_eta(k)

  end subroutine propagate_eta_k

  ! propagate v_eta_k, simple shift step
  subroutine propagate_v_eta_k_1(th, k, dt, dtfac)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: k
    real(double), intent(in)              :: dt, dtfac

    if (th%iprint == 0) write(*,'(6x,"NHC: propagating v_eta linear shift, k = ",i2)') k
    ! th%v_eta(k) = th%v_eta(k) + half*dt*((th%v_eta(k-1)**2)*th%Q(k-1) - &
    !              th%k_B_md*th%T_ext)/th%Q(k)
    th%v_eta(k) = th%v_eta(k) + dtfac*dt*th%G_nhc(k)

  end subroutine propagate_v_eta_k_1

  ! propagate v_eta_k, exponential shift step
  subroutine propagate_v_eta_k_2(th, k, dt, dtfac)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: k
    real(double), intent(in)              :: dt, dtfac

    if (th%iprint == 0) write(*,'(6x,"NHC: propagating v_eta exp factor,   k = ",i2)') k
    th%v_eta(k) = th%v_eta(k)*exp(-dtfac*dt*th%v_eta(k+1))

  end subroutine propagate_v_eta_k_2

  ! Propagate the NHC
  subroutine propagate_nvt_nhc(th, dt, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), intent(in)                    :: dt
    real(double), dimension(:,:), intent(inout) :: v

    ! local variables
    integer       :: i, j, k
    real(double)  :: F_sys
    real(double)  :: dtys ! Yoshida-Suzuki time step
    real(double)  :: v_sfac ! velocity scaling factor

    if (th%iprint == 0) write(*,'(4x,a)') "Thermostat: propagating NHC dt/2 update"

    v_sfac = one
    do i=1,th%mts_nhc ! MTS loop
      do j=1,th%n_ys ! Yoshida-Suzuki loop
        ! Reverse part of expansion: update forces and thermostat velocities
        do k=th%n_nhc,1,-1
          if (k==th%n_nhc) then
            call th%update_G_k(k, zero)
            call th%propagate_v_eta_k_1(k, dt, quarter)
          else
            ! Trotter expansion to avoid sinh singularity (see MTTK paper)
            call th%propagate_v_eta_k_2(k, dt, eighth)
            call th%propagate_v_eta_k_1(k, dt, quarter)
            call th%update_G_k(k, zero)
            call th%propagate_v_eta_k_2(k, dt, eighth)
          end if
        end do
        v_sfac = v_sfac * exp(-half*dt*th%v_eta(1))
        th%ke_system = th%ke_system*exp(-dt*th%v_eta(1))
        ! update the nhc "positions" eta
        do k=1,th%n_nhc
          call th%propagate_eta_k(k, dt, half)
        end do
        ! Forward part of expansion: update forces and thermostat velocities
        do k=1,th%n_nhc
          if (k<th%n_nhc) then
            ! Trotter expansion to avoid sinh singularity (see MTTK paper)
            call th%propagate_v_eta_k_2(k, dt, eighth)
            call th%update_G_k(k, zero)
            call th%propagate_v_eta_k_1(k, dt, quarter)
            call th%propagate_v_eta_k_2(k, dt, eighth)
          else
            call th%update_G_k(k, zero)
            call th%propagate_v_eta_k_1(k, dt, quarter)
          end if
        end do
      end do ! Yoshida-Suzuki loop
    end do ! MTS loop

    v = v_sfac*v

  end subroutine propagate_nvt_nhc

  ! Compute the "NHC energy" for monitoring the conserved quantity (KE + PE + NHC energy)
  subroutine get_nhc_ke(th)

    ! passed variables
    class(type_thermostat), intent(inout) :: th

    if (th%n_nhc > 0) then
      th%e_nhc = half*sum(th%Q*th%v_eta**2)
      th%e_nhc = th%e_nhc + th%ndof*th%k_B_md*th%T_ext*th%eta(1)
      th%e_nhc = th%e_nhc + th%k_B_md*th%T_ext*sum(th%eta(2:))
    else
      th%e_nhc = zero
    end if
    
  end subroutine get_nhc_ke

  ! Debugging routine: dump the state of the thermostat
  subroutine dump_thermo_state(th, step, funit)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: step, funit

    ! local variables
    integer                               :: i
    character(40)                         :: fmt


    write(fmt,'("(a8,",i4,"e16.4)")') th%n_nhc
    write(funit,'("step   ",i12)') step
    if (th%th_type == 'nhc') then
      write(funit,fmt) "eta:    ", th%eta
      write(funit,fmt) "v_eta:  ", th%v_eta
      write(funit,fmt) "G_nhc:  ", th%G_nhc
      write(funit,'("e_nhc:  ",e16.4)') th%e_nhc
    else if (th%th_type == 'berendsen') then
      write(funit,'("lambda: ",e16.4)') th%lambda
    end if
    write(funit,*)

  end subroutine dump_thermo_state

end module thermostat_module
