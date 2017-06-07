module thermostat_module

use datatypes
use constants
use rng

implicit none

type type_thermostat

  real(double)          :: T_ext, k_B_md
  integer               :: nat, ndof, iprint
  integer               :: tau_T ! coupling time period of v-rescale thermostat
  character(40)         :: th_type
  integer                                   :: n_nhc  ! number of nh chains
  real(double), dimension(:), allocatable   :: eta   ! thermostat k position
  real(double), dimension(:), allocatable   :: v_eta ! thermostat k momentum
  real(double), dimension(:), allocatable   :: Q    ! thermostat k mass

contains

  procedure :: init_thermostat
  procedure :: velocity_rescale
  procedure :: weak_coupling_t
  procedure :: propagate_vr_thermostat
  procedure :: init_nhc
  procedure :: propagate_v
  procedure :: propagate_eta_k
  procedure :: propagate_v_eta_1
  procedure :: propagate_v_eta_k_1
  procedure :: propagate_v_eta_k_2
  procedure :: propagate_nhc
  procedure :: get_nhc_energy

end type type_thermostat

contains

  ! Initialise the thermostat
  subroutine init_thermostat(th, th_type, nat, ndof, T_ext, tau_T, n_nhc, &
                             iprint)

    ! passed variables
    class(type_thermostat), intent(inout)     :: th
    character(40), intent(in)                 :: th_type
    integer, intent(in)                       :: nat, ndof, iprint, tau_T, n_nhc
    real(double), intent(in)                  :: T_ext

    th%th_type = th_type
    th%nat = nat
    th%ndof = ndof
    th%T_ext = T_ext
    th%iprint = iprint
    th%tau_T = tau_T
    th%k_B_md = one
    th%n_nhc = n_nhc
    if (th_type == 'nhc') call th%init_nhc(th%n_nhc)

  end subroutine init_thermostat

  ! Isokinetic thermostat: maintain temperature by rescaling velocities by
  ! sqrt(T_int/T_ext) every tau_T steps
  subroutine velocity_rescale(th, T_int, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v
    real(double), intent(in)                    :: T_int

    ! local variables
    real(double)  :: lambda

    lambda = sqrt(th%T_ext/T_int)
    v = v*lambda
    if (th%iprint == 0) write(*,'(6x,"Velocity rescale, lambda = ",f8.4)') lambda
  end subroutine velocity_rescale

  ! Berendsen weak coupling thermostat - a different type of velocity rescaling
  subroutine weak_coupling_t(th, T_int, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v
    real(double), intent(in)                    :: T_int

    ! local variables
    real(double)  :: lambda

    lambda = sqrt(one + (one/th%tau_T)*(T_int/th%T_ext - one))
    v = v*lambda
    if (th%iprint == 0) write(*,'(6x,"Weak coupling, lambda = ",f8.4)') lambda
  end subroutine weak_coupling_t

  ! Update the thermostat
  subroutine propagate_vr_thermostat(th, T_int, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v
    real(double), intent(in)                    :: T_int

    if (th%iprint == 0) write(*,'(4x,a)') "Thermostat: rescaling velocities"
    select case(th%th_type)
    case ('velocity_rescale')
      call th%velocity_rescale(T_int, v)
    case ('weak_coupling')
      call th%weak_coupling_t(T_int, v)
    end select

  end subroutine propagate_vr_thermostat

  ! Initialise the NHC with n_nhc heat baths
  subroutine init_nhc(th, n_nhc)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: n_nhc

    th%n_nhc = n_nhc
    allocate(th%eta(n_nhc))
    allocate(th%v_eta(n_nhc))
    allocate(th%Q(n_nhc))

    th%eta = zero
    th%v_eta = zero
    th%Q = one

  end subroutine init_nhc

  ! Propagate the velocity (coupling with NHC heat bath k=1)
  subroutine propagate_v(th, dt, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v
    real(double)                                :: dt

    ! local variables
    integer                                   :: i

    v = v - half*dt*th%v_eta(1)

  end subroutine propagate_v

  ! propagate eta (this is the same for all k)
  subroutine propagate_eta_k(th, k, dt)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer                               :: k
    real(double)                          :: dt

    if (th%iprint == 0) write(*,'(6x,"NHC: propagating eta,                k = ",i2)') k
    th%eta(k) = th%eta(k) + half*dt*th%v_eta(k)

  end subroutine propagate_eta_k

  ! propagate v_eta_1 (k=1 only, couples NHC 1 to the system)
  ! I'm using the kinetic energy from the previous step in both updates
  subroutine propagate_v_eta_1(th, dt, ke_t)

    ! passed variables
    class(type_thermostat), intent(inout)     :: th
    real(double), intent(in)                  :: dt
    real(double), intent(in)                  :: ke_t

    ! local variables
    integer                                   :: i

    if (th%iprint == 0) write(*,'(6x,"NHC: propagating v_eta linear shift, k = ",i2)') 1
    th%v_eta(1) = th%v_eta(1) + half*dt*(ke_t - &
                  3*th%nat*th%k_B_md*th%T_ext)/th%Q(1)

  end subroutine propagate_v_eta_1

  ! propagate v_eta_k, simple shift step
  subroutine propagate_v_eta_k_1(th, k, dt)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: k
    real(double), intent(in)              :: dt

    if (th%iprint == 0) write(*,'(6x,"NHC: propagating v_eta linear shift, k = ",i2)') k
    th%v_eta(k) = th%v_eta(k) + half*dt*((th%v_eta(k-1)**2)*th%Q(k-1) - &
                 th%k_B_md*th%T_ext)/th%Q(k)

  end subroutine propagate_v_eta_k_1

  ! propagate v_eta_k, exponential shift step
  subroutine propagate_v_eta_k_2(th, k, dt)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: k
    real(double), intent(in)              :: dt

    if (th%iprint == 0) write(*,'(6x,"NHC: propagating v_eta exp factor,   k = ",i2)') k
    th%v_eta(k) = th%v_eta(k)*exp(-half*dt*th%v_eta(k+1))

  end subroutine propagate_v_eta_k_2

  ! Propagate the NHC
  subroutine propagate_nhc(th, dt, ke_t, reverse, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), intent(in)                    :: dt
    real(double), intent(in)                    :: ke_t
    logical, intent(in)                         :: reverse
    real(double), dimension(:,:), intent(inout) :: v

    ! local variables
    integer                                     :: k

    if (th%iprint == 0) write(*,'(4x,a)') "Thermostat: propagating NHC dt/2 update"
    if (reverse .eqv. .false.) then
      ! propagate chain starting from NHC heat bath M (end of chain)
      do k=th%n_nhc,1,-1
        if (k < th%n_nhc) call th%propagate_v_eta_k_2(k,dt)
        if (k==1) then
          call th%propagate_v_eta_1(dt, ke_t)
        else
          call th%propagate_v_eta_k_1(k, dt)
        end if
        call th%propagate_eta_k(k, dt)
      end do
      call th%propagate_v(dt, v)
    else
      ! propagate chain starting from NHC heat bath 1 (coupled to system)
      call th%propagate_v(dt, v)
      do k=1,th%n_nhc
        call th%propagate_eta_k(k, dt)
        if (k==1) then
          call th%propagate_v_eta_1(dt, ke_t)
        else
          call th%propagate_v_eta_k_1(k, dt)
        end if
        if (k < th%n_nhc) call th%propagate_v_eta_k_2(k,dt)
      end do
    end if

  end subroutine propagate_nhc

  ! Compute the "NHC energy" for monitoring the conserved quantity (KE + PE + NHC energy)
  subroutine get_nhc_energy(th, nhc_energy)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    real(double), intent(out)             :: nhc_energy

    nhc_energy = half*sum((th%v_eta**2)/th%Q)
    nhc_energy = nhc_energy + 3*th%nat*th%k_B_md*th%T_ext*th%eta(1)
    nhc_energy = nhc_energy + th%k_B_md*th%T_ext*sum(th%eta(2:))
  end subroutine get_nhc_energy

end module thermostat_module
