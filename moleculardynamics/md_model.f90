module md_model

use datatypes

implicit none

type model
  
  ! Positions, velocities, forces
  real(double), dimension(:,:), allocatable :: r, r0
  real(double), dimension(:,:), allocatable :: r_cart, r_cart0
  real(double), dimension(:,:), allocatable :: v, v0
  real(double), dimension(:,:), allocatable :: f, f0

  ! box
  real(double), dimension(3,3)              :: h, h0

  !  MD variables
  integer                                   :: nstep
  real(double)                              :: dt

  ! Microscopic variables
  integer                                   :: nat
  integer                                   :: ndof

  ! Thermodynamic variables
  real(double)                              :: T, T_ext
  real(double)                              :: P, P_ext
  real(double)                              :: vol, vol0
  real(double)                              :: Ek, Ek0
  real(double)                              :: Ep, Ep0
  real(double)                              :: enthalpy, enthalpy0
  real(double)                              :: pV, pV0
  real(double)                              :: Hprime
  real(double), dimension(3,3)              :: stress

  ! Thermostat
  character(40)                             :: thermo_type
  integer                                   :: n_nhc
  real(double), dimension(:), allocatable   :: Q
  real(double), dimension(:), allocatable   :: eta
  real(double), dimension(:), allocatable   :: v_eta
  real(double)                              :: G_eta
  real(double)                              :: E_nhc

  ! Barostat
  character(40)                             :: baro_type
  real(double)                              :: W_eps
  real(double)                              :: eps
  real(double)                              :: v_eps
  real(double)                              :: G_eps
  real(double)                              :: E_baro
  real(double), dimension(3,3)              :: c_g

contains
  procedure :: init_model

end type model

contains

  subroutine init_model(mdl, nat)

    ! passed variables
    class(model), intent(inout)     :: mdl
    integer, intent(in)             :: nat

    mdl%nat = nat
    mdl%ndof = 3*nat

    allocate(mdl%r(nat,3), mdl%r0(nat,3))
    allocate(mdl%r_cart(nat,3), mdl%r_cart0(nat,3))
    allocate(mdl%v(nat,3), mdl%v0(nat,3))
    allocate(mdl%f(nat,3), mdl%f0(nat,3))

  end subroutine init_model

end module md_model