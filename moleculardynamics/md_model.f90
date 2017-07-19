module md_model

use datatypes

implicit none

type type_model
  
  ! Positions, velocities, forces
  real(double), dimension(:,:), allocatable :: r, r0
  real(double), dimension(:,:), allocatable :: rcart, rcart0
  real(double), dimension(:,:), allocatable :: v, v0
  real(double), dimension(:,:), allocatable :: f, f0

  ! Species
  integer, allocatable, dimension(:)        :: species
  character(2), allocatable, dimension(:)   :: spec_label
  real(double), allocatable, dimension(:)   :: mass

  ! box
  real(double), dimension(3,3)              :: h, h0

  !  MD variables
  integer                                   :: nstep
  real(double)                              :: dt
  character(3)                              :: ensemble

  ! Microscopic variables
  integer                                   :: nat
  integer                                   :: nspec
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

end type type_model

contains

  subroutine init_model(mdl, nat, nspec)

    ! passed variables
    class(type_model), intent(inout)  :: mdl
    integer, intent(in)               :: nat, nspec

    mdl%nat = nat
    mdl%nspec = nspec
    mdl%ndof = 3*nat

    allocate(mdl%r(nat,3), mdl%r0(nat,3))
    allocate(mdl%rcart(nat,3), mdl%rcart0(nat,3))
    allocate(mdl%v(nat,3), mdl%v0(nat,3))
    allocate(mdl%f(nat,3), mdl%f0(nat,3))
    allocate(mdl%species(nat))
    allocate(mdl%spec_label(nspec))
    allocate(mdl%mass(nspec))

  end subroutine init_model

end module md_model