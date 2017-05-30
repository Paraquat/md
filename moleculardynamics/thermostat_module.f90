module thermostat_module

use datatypes
use constants
use rng

implicit none

type type_thermostat

  real(double)          :: T_ext
  integer               :: nat, ndof, iprint
  character(40)         :: th_type

contains

  procedure :: init_thermostat
  procedure :: velocity_rescale
  procedure :: propagate_thermostat

end type type_thermostat

contains

  ! Initialise the thermostat
  subroutine init_thermostat(th, th_type, nat, ndof, T_ext, iprint)

    ! passed variables
    class(type_thermostat), intent(inout)     :: th
    character(40), intent(in)                 :: th_type
    integer, intent(in)                       :: nat, ndof, iprint
    real(double), intent(in)                  :: T_ext ! target temperature

    th%th_type = th_type
    th%nat = nat
    th%ndof = ndof
    th%T_ext = T_ext
    th%iprint = iprint

  end subroutine init_thermostat

  ! Isokinetic thermostat: maintain temperature by rescaling velocities by
  ! sqrt(T_int/T_ext) every tau_T steps
  subroutine velocity_rescale(th, T_int, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v
    real(double), intent(in)                    :: T_int

    ! local variables
    integer       :: i
    real(double)  :: lambda

    lambda = sqrt(th%T_ext/T_int)
    v = v*lambda
    if (th%iprint == 0) write(*,'("Velocity rescale, lambda = ",f8.4)') lambda
  end subroutine velocity_rescale

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
    end select

  end subroutine propagate_thermostat

end module thermostat_module
