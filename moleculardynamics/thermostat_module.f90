module thermostat_module

use datatypes
use constants
use rng

implicit none

type type_thermostat

  real(double)          :: T_ext
  integer               :: nat, ndof, iprint
  integer               :: tau_T ! coupling time period of v-rescale thermostat
  character(40)         :: th_type

contains

  procedure :: init_thermostat
  procedure :: velocity_rescale
  procedure :: weak_coupling_t
  procedure :: propagate_thermostat

end type type_thermostat

contains

  ! Initialise the thermostat
  subroutine init_thermostat(th, th_type, nat, ndof, T_ext, tau_T, iprint)

    ! passed variables
    class(type_thermostat), intent(inout)     :: th
    character(40), intent(in)                 :: th_type
    integer, intent(in)                       :: nat, ndof, iprint, tau_T
    real(double), intent(in)                  :: T_ext

    th%th_type = th_type
    th%nat = nat
    th%ndof = ndof
    th%T_ext = T_ext
    th%iprint = iprint
    th%tau_T = tau_T

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
    if (th%iprint == 0) write(*,'("Velocity rescale, lambda = ",f8.4)') lambda
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
    if (th%iprint == 0) write(*,'("Weak coupling, lambda = ",f8.4)') lambda
  end subroutine weak_coupling_t

  ! Update the thermostat
  subroutine propagate_thermostat(th, T_int, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v
    real(double), intent(in)                    :: T_int

    if (th%iprint == 0) write(*,'(a)') "Propagating thermostat"
    select case(th%th_type)
    case ('velocity_rescale')
      call th%velocity_rescale(T_int, v)
    case ('weak_coupling')
      call th%weak_coupling_t(T_int, v)
    end select

  end subroutine propagate_thermostat

end module thermostat_module
