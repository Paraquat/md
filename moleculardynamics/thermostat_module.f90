module thermostat_module

use datatypes
use constants
use rng

implicit none

type type_thermostat

  real(double)          :: sumv2, T_int, T_ext
  integer               :: nat, ndof
  character(40)         :: th_type

contains

  procedure :: init_thermostat
  procedure :: get_T
  procedure :: velocity_rescale

end type type_thermostat

contains

  ! Initialise the thermostat
  subroutine init_thermostat(th, th_type, nat, ndof, T_ext)

    ! passed variables
    class(type_thermostat), intent(inout)     :: th
    character(40), intent(in)                 :: th_type
    integer, intent(in)                       :: nat, ndof
    real(double), intent(in)                  :: T_ext ! target temperature

    th%th_type = th_type
    th%nat = nat
    th%ndof = ndof
    th%T_ext = T_ext

  end subroutine init_thermostat

  ! Compute the temperature
  subroutine get_T(th, sumv2)

    ! passed variables
    class(type_thermostat), intent(inout)     :: th
    real(double), intent(in)                  :: sumv2

    th%T_int = sumv2/th%ndof
  end subroutine get_T

  ! Isokinetic thermostat: maintain temperature by rescaling velocities by
  ! sqrt(T_int/T_ext) every tau_T steps
  subroutine velocity_rescale(th, v, sumv2)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v
    real(double), intent(in)                    :: sumv2

    ! local variables
    integer       :: i
    real(double)  :: lambda

    call th%get_T(sumv2)
    lambda = sqrt(th%T_ext/th%T_int)
    v = v*lambda
  end subroutine velocity_rescale

  ! Update the thermostat
  subroutine propagate_thermostat(th, v, sumv2)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v
    real(double), intent(in)                    :: sumv2

  select case(th%th_type)
  case ('velocity_rescale')
    call th%velocity_rescale(v, sumv2)
  end select

  end subroutine propagate_thermostat

end module thermostat_module
