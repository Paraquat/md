module isotropic_barostat_module

use datatypes
use constants
use maths, only: exp_sinhx_x

implicit none

type type_iso_barostat

  character(40)                 :: baro_type
  real(double), dimension(3,3)  :: h      ! current cell
  real(double), dimension(3,3)  :: h_0    ! reference cell
  real(double)                  :: V      ! box volume
  real(double)                  :: V_ref  ! reference box volume
  real(double)                  :: eps    ! third*log(V/V_ref)
  real(double)                  :: P_int  ! internal pressure excluding kinetic
  real(double)                  :: P_ext  ! applied pressure
  real(double)                  :: W_eps  ! box mass
  real(double)                  :: v_eps  ! box velocity
  real(double)                  :: G_eps  ! box force
  real(double)                  :: G_nhc_1 ! force from first NHC heat bath
  integer                       :: iprint
  integer                       :: nat
  integer                       :: ndof   ! degrees of freedom
  real(double)                  :: ke_box ! box kinetic energy
  real(double)                  :: k_B_md
  real(double)                  :: odnf   ! 1 + d/N_f

contains

  procedure :: init_barostat_iso
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

  subroutine init_barostat_iso(baro, baro_type, h_ref, P_ext, nat, ndof, &
                               box_mass, volume, iprint)

    ! Passed variables
    class(type_iso_barostat), intent(inout)   :: baro
    character(40), intent(in)                 :: baro_type
    real(double), dimension(3,3), intent(in)  :: h_ref
    real(double), intent(in)                  :: P_ext, box_mass, volume
    integer, intent(in)                       :: nat, ndof, iprint

    ! Local variables
    integer       :: i
    real(double)  :: W_eps

    baro%baro_type = baro_type
    baro%h_0 = h_ref
    baro%h = h_ref

    baro%V = volume
    baro%V_ref = volume
    baro%nat = nat
    baro%ndof = ndof
    baro%odnf = (one + three/baro%ndof)
    baro%eps = zero
    baro%iprint = iprint

    baro%v_eps = zero
    baro%k_B_md = one

    baro%W_eps = box_mass
    baro%ke_box = zero

  end subroutine init_barostat_iso


  ! Get the box kinetic energy
  subroutine get_box_ke_iso(baro)

    ! Passed variables
    class(type_iso_barostat), intent(inout)  :: baro

    baro%ke_box = baro%W_eps*baro%v_eps**2

  end subroutine get_box_ke_iso

  ! Update the box forces
  subroutine update_G_eps(baro, ke, P_int, volume, T_ext, Q_1, G_nhc_1)

    ! Passed variables
    class(type_iso_barostat), intent(inout)   :: baro
    real(double), intent(in)                  :: ke
    real(double), intent(in)                  :: P_int
    real(double), intent(in)                  :: volume
    real(double), intent(in)                  :: T_ext
    real(double), intent(in)                  :: Q_1 ! mass of NHC k=1
    real(double), intent(out)                 :: G_nhc_1 ! force on NHC k=1

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: updating box force G_eps"
    baro%P_int = P_int
    baro%G_eps = (baro%odnf*ke + three*(P_int - &
                  baro%P_ext)*volume)/baro%W_eps
    ! The force on the first NHC heat bath is required for thermostat coupling
    G_nhc_1 = (ke + baro%ke_box - baro%ndof*baro%k_B_md*T_ext)/Q_1
    baro%G_nhc_1 = G_nhc_1

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
    real(double), intent(in)                  :: dt, dtfac, v_eta

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: propagating v_eps exp factor"
    baro%v_eps = baro%v_eps*exp(-dtfac*dt*v_eta)

  end subroutine propagate_v_eps_2

  ! coupling box and particle velocities
  subroutine propagate_v_sys_iso(baro, dt, dtfac, v_eta_1, vin, vout)

    ! Passed variables
    class(type_iso_barostat), intent(inout)     :: baro
    real(double), intent(in)                    :: dt, dtfac
    real(double), intent(in)                    :: v_eta_1
    real(double), dimension(:,:), intent(in)    :: vin
    real(double), dimension(:,:), intent(out)   :: vout

    ! local variables
    real(double)                                :: vscale

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: propagating particle velocities"

    vscale = exp(-dtfac*dt*v_eta_1)
    vout = vin*vscale

  end subroutine propagate_v_sys_iso

  ! coupling box and particle positions. This is a little complicated and
  ! contains a polynomial expansion of sinh(x)/x, so I will do it in one
  ! subroutine for now, as in the MTTK paper.
  subroutine propagate_r_sys(baro, dt, dtfac, v, rin, rout)

    ! Passed variables
    class(type_iso_barostat), intent(inout)     :: baro
    real(double), intent(in)                    :: dt, dtfac
    real(double), dimension(:,:), intent(in)    :: v
    real(double), dimension(:,:), intent(in)    :: rin
    real(double), dimension(:,:), intent(out)   :: rout

    ! local variables
    integer                                     :: i, j
    real(double), dimension(3)                  :: u, uv
    real(double)                                :: sinhx_x, rscale, vscale, &
                                                   exp_v_eps, eps

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: propagating particle positions"
    exp_v_eps = exp(dt*dtfac*baro%v_eps)
    sinhx_x = exp_sinhx_x((dtfac*dt*baro%v_eps)**2)
    rscale = exp_v_eps**2
    vscale = exp_v_eps*sinhx_x*dt
    do i=1,baro%nat
      do j=1,3
        rout(i,j) = rin(i,j)*rscale + v(i,j)*vscale
      end do 
    end do

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
    write(funit,'("eps:   ",e16.4)') baro%eps
    write(funit,'("v_eps: ",e16.4)') baro%v_eps
    write(funit,'("G_eps: ",e16.4)') baro%G_eps
    write(funit,'("G_nhc: ",e16.4)') baro%G_nhc_1
    write(funit,'("ke_box:",e16.4)') baro%ke_box
    write(funit,'("P_int: ",e16.4)') baro%P_int
    write(funit,'("V:     ",e16.4)') baro%V
    write(funit,*)

  end subroutine dump_baro_state_iso

end module isotropic_barostat_module
