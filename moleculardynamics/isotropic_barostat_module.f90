module isotropic_barostat_module

use datatypes
use constants
use maths, only: exp_sinhx_x

implicit none

type type_iso_barostat

  character(40)                 :: baro_type
  character(3)                  :: ensemble
  real(double)                  :: dt
  integer                       :: iprint
  integer                       :: nat
  integer                       :: ndof   ! degrees of freedom
  real(double)                  :: k_B_md
  real(double), dimension(3,3)  :: h      ! current cell
  real(double), dimension(3,3)  :: h_0    ! reference cell
  real(double), dimension(3,3)  :: stress
  real(double)                  :: V      ! box volume
  real(double)                  :: V_ref  ! reference box volume
  real(double)                  :: P_int  ! internal pressure excluding kinetic
  real(double)                  :: P_ext  ! applied pressure

  ! MTTK variables
  real(double)                  :: eps    ! third*log(V/V_ref)
  real(double)                  :: eps_ref ! reference cell epsilon
  real(double)                  :: W_eps  ! box mass
  real(double)                  :: v_eps  ! box velocity
  real(double)                  :: G_eps  ! box force
  real(double)                  :: G_nhc_1 ! force from first NHC heat bath
  real(double)                  :: ke_box ! box kinetic energy
  real(double)                  :: odnf   ! 1 + d/N_f

  ! Berendsen variables
  real(double)                  :: tau_P  ! pressure coupling time constant
  real(double)                  :: beta   ! Isothermal compressibility
  real(double)                  :: mu     ! Berensen scaling factor

contains

  procedure :: init_barostat_iso
  procedure :: get_berendsen_baro_sf
  procedure :: berendsen_baro_propagate
  procedure :: get_box_ke_iso
  procedure :: update_G_eps
  procedure :: propagate_eps_1
  procedure :: propagate_eps_2
  procedure :: propagate_v_eps_1
  procedure :: propagate_v_eps_2
  procedure :: propagate_v_sys_iso
  procedure :: propagate_r_sys
  procedure :: propagate_box
  procedure :: dump_baro_state_iso

end type type_iso_barostat

contains

  subroutine init_barostat_iso(baro, baro_type, ensemble, dt, h_ref, P_ext, &
                               nat, ndof, box_mass, volume, tau_P, iprint)

    ! Passed variables
    class(type_iso_barostat), intent(inout)   :: baro
    character(40), intent(in)                 :: baro_type
    character(3), intent(in)                  :: ensemble
    real(double), dimension(3,3), intent(in)  :: h_ref
    real(double), intent(in)                  :: P_ext, box_mass, volume, &
                                                 tau_P, dt
    integer, intent(in)                       :: nat, ndof, iprint

    ! Local variables
    integer       :: i

    baro%iprint = iprint
    baro%baro_type = baro_type
    baro%ensemble = ensemble
    baro%dt = dt
    baro%ndof = ndof
    baro%odnf = one + three/baro%ndof
    baro%k_B_md = one
    baro%h_0 = h_ref
    baro%h = h_ref
    baro%P_ext = P_ext

    ! initialise MTTK
    baro%V_ref = volume
    baro%V = baro%V_ref
    baro%eps_ref = third*log(baro%V/baro%V_ref)
    baro%eps = baro%eps_ref
    baro%v_eps = zero
    baro%W_eps = box_mass
    baro%ke_box = zero

    ! Intialise Berendsen
    baro%beta = one
    baro%tau_P = tau_P

  end subroutine init_barostat_iso

  ! Compute the pressure scaling factor. This is done separately because the
  ! rescaling should not affect the velocity
  subroutine get_berendsen_baro_sf(baro, P_int)

    ! Passed variables
    class(type_iso_barostat), intent(inout)     :: baro
    real(double), intent(in)                    :: P_int

    baro%mu = (one - (baro%dt/baro%tau_P)*baro%beta*(baro%P_ext - P_int))**third
    if (baro%iprint == 0) then
      write(*,'(4x,"Weak coupling barostat, mu = ",f8.4)') baro%mu
    end if

  end subroutine get_berendsen_baro_sf

  ! Scale the particle coordinates and box according to the Berendsen weak
  ! coupling method
  subroutine berendsen_baro_propagate(baro, r, h)

    ! Passed variables
    class(type_iso_barostat), intent(inout)     :: baro
    real(double), dimension(:,:), intent(inout) :: r
    real(double), dimension(3,3), intent(inout) :: h

    h = h*baro%mu
    r = r*baro%mu

  end subroutine berendsen_baro_propagate

  ! Get the box kinetic energy
  subroutine get_box_ke_iso(baro)

    ! Passed variables
    class(type_iso_barostat), intent(inout)  :: baro

    baro%ke_box = baro%W_eps*baro%v_eps**2

  end subroutine get_box_ke_iso

  ! Update the box forces
  subroutine update_G_eps(baro, ke, P_int, volume)

    ! Passed variables
    class(type_iso_barostat), intent(inout)   :: baro
    real(double), intent(in)                  :: ke
    real(double), intent(in)                  :: P_int
    real(double), intent(in)                  :: volume

    ! local variables
    real(double)                              :: akin

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: updating box force G_eps"
    akin = ke*two
    baro%P_int = P_int
    baro%G_eps = (baro%odnf*akin + three*(P_int - &
                  baro%P_ext)*volume)/baro%W_eps

    write(14,*) baro%P_int, baro%odnf*akin, three*(P_int - baro%P_ext)*volume, baro%G_eps, baro%v_eps

  end subroutine update_G_eps

  ! linear shift in epsilon
  subroutine propagate_eps_1(baro, dt, dtfac)

    ! Passed variables
    class(type_iso_barostat), intent(inout)   :: baro
    real(double), intent(in)                  :: dt
    real(double), intent(in)                  :: dtfac

    if (baro%iprint == 0) write(*,'(4x,a)') "MTTK: propagating eps linear shift"
    baro%eps = baro%eps + dtfac*dt*baro%v_eps

  end subroutine propagate_eps_1

  ! exponential shift in epsilon
  subroutine propagate_eps_2(baro, dt, dtfac, v_eta_1)

    ! Passed variables
    class(type_iso_barostat), intent(inout)   :: baro
    real(double), intent(in)                  :: dt
    real(double), intent(in)                  :: dtfac
    real(double), intent(in)                  :: v_eta_1

    if (baro%iprint == 0) write(*,'(4x,a)') "MTTK: propagating eps exp factor"
    baro%eps = baro%eps*exp(-dtfac*dt*v_eta_1)

  end subroutine propagate_eps_2

  ! linear shift in box velocity
  subroutine propagate_v_eps_1(baro, dt, dtfac)

    ! Passed variables
    class(type_iso_barostat), intent(inout)   :: baro
    real(double), intent(in)                  :: dt, dtfac

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: propagating v_eps linear shift"
    baro%v_eps = baro%v_eps + dtfac*dt*baro%G_eps

  end subroutine propagate_v_eps_1

  ! exponential shift in box velocity from thermostat coupling
  subroutine propagate_v_eps_2(baro, dt, dtfac, v_eta)

    ! Passed variables
    class(type_iso_barostat), intent(inout)   :: baro
    real(double), intent(in)                  :: dt, dtfac
    real(double), intent(in)                  :: v_eta ! velocity of
                                                       ! thermostat 1

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: propagating v_eps exp factor"
    baro%v_eps = baro%v_eps*exp(-dtfac*dt*v_eta)

  end subroutine propagate_v_eps_2

  ! coupling box and particle velocities
  subroutine propagate_v_sys_iso(baro, dt, dtfac, v_eta_1, v)

    ! Passed variables
    class(type_iso_barostat), intent(inout)     :: baro
    real(double), intent(in)                    :: dt, dtfac
    real(double), intent(in)                    :: v_eta_1
    real(double), dimension(:,:), intent(inout) :: v

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: propagating particle velocities"

    v = v*exp(-dtfac*dt*(v_eta_1 + baro%odnf*baro%v_eps))

  end subroutine propagate_v_sys_iso

  ! coupling box and particle positions. This is a little complicated and
  ! contains a polynomial expansion of sinh(x)/x, so I will do it in one
  ! subroutine for now, as in the MTTK paper.
  subroutine propagate_r_sys(baro, dt, dtfac, v, r)

    ! Passed variables
    class(type_iso_barostat), intent(inout)     :: baro
    real(double), intent(in)                    :: dt, dtfac
    real(double), dimension(:,:), intent(in)    :: v
    real(double), dimension(:,:), intent(inout) :: r

    ! local variables
    integer                                     :: i, j
    real(double), dimension(3)                  :: u, uv
    real(double)                                :: sinhx_x, rscale, vscale, &
                                                   exp_v_eps, eps

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: propagating particle positions"
    exp_v_eps = exp(dt*dtfac*baro%v_eps)
    sinhx_x = exp_sinhx_x(dtfac*dt*baro%v_eps)
    rscale = exp_v_eps**2
    vscale = exp_v_eps*sinhx_x*dt
    r = r*rscale + v*vscale

  end subroutine propagate_r_sys

  ! Update the lattice vectors
  subroutine propagate_box(baro, h, V)

    ! Passed variables
    class(type_iso_barostat), intent(inout)     :: baro
    real(double), dimension(3,3), intent(inout) :: h
    real(double), intent(inout)                 :: V

    ! Local variables
    integer                                     :: i
    real(double)                                :: V_new, V_old, hscale

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: propagating box"
    V_old = V
    v_new = exp(three*baro%eps)*baro%V_ref
    hscale = (V_new/V_old)**third
    h = hscale*h
    baro%h = h
    baro%V = V_new

  end subroutine propagate_box

  ! Debugging routine: dump the state of the barostat
  subroutine dump_baro_state_iso(baro, step, funit)

    ! Passed variables
    class(type_iso_barostat), intent(inout)     :: baro
    integer, intent(in)                         :: step, funit

    ! local variables
    integer                                     :: i

    write(funit,'("step   ",i12)') step
    write(funit,'(8x,a)') "Lattice vectors"
    do i=1,3
      write(funit,'(4x,3f12.4)') baro%h(i,:)
    end do
    write(funit,'(8x,a)') "Stress tensor"
    do i=1,3
      write(funit,'(4x,3f12.4)') baro%stress(i,:)
    end do
    if (baro%baro_type == 'mttk' .or. baro%baro_type == 'iso-mttk') then
      write(funit,'("eps:   ",e16.4)') baro%eps
      write(funit,'("v_eps: ",e16.4)') baro%v_eps
      write(funit,'("ke_box:",e16.4)') baro%ke_box
      write(funit,'("G_eps: ",e16.4)') baro%G_eps
      write(funit,'("G_nhc: ",e16.4)') baro%G_nhc_1
    else if (baro%baro_type == 'berendsen') then
      write(funit,'("mu:    ",e16.4)') baro%mu
    end if
    write(funit,'("P_int: ",e16.4)') baro%P_int
    write(funit,'("V:     ",e16.4)') baro%V
    write(funit,*)

  end subroutine dump_baro_state_iso

end module isotropic_barostat_module
