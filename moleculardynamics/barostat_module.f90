module barostat_module

use datatypes
use constants
use maths, only: poly_sinhx_x

implicit none

type type_barostat

  character(40)                 :: baro_type
  character(3)                  :: ensemble
  real(double)                  :: dt
  integer                       :: iprint
  integer                       :: nat
  integer                       :: ndof
  real(double)                  :: k_B_md

  real(double)                  :: P_int    ! internal pressure
  real(double)                  :: P_ext    ! applied pressure
  real(double)                  :: V
  real(double), dimension(3,3)  :: h        ! current cell
  real(double), dimension(3,3)  :: h_0      ! refernce cell
  real(double), dimension(3,3)  :: stress   ! stress tensor
  real(double), dimension(3,3)  :: virial_ke ! ke contribution to stress tensor

  ! MTTK variables
  real(double)                  :: W_h      ! box mass
  real(double), dimension(3,3)  :: v_h      ! box velocity
  real(double), dimension(3,3)  :: G_h      ! box force
  real(double), dimension(3)    :: G_nhc_1  !
  real(double), dimension(3,3)  :: ident    ! 3x3 identity matrix
  real(double), dimension(3,3)  :: onfm     ! ident*(1/ndof)
  real(double), dimension(3,3)  :: c_g      ! eigenvectors 
  real(double), dimension(3)    :: lambda   ! eigenvalues
  real(double)                  :: ke_box

  ! Berendsen variables
  real(double)                  :: tau_P    ! pressure coupling time constatn
  real(double)                  :: beta     ! isothermal compressibility
  real(double)                  :: mu       ! Berendsen scaling factor

contains

  procedure :: init_barostat
  procedure :: get_lambda_cg
  procedure :: get_lambda_sfac
  procedure :: get_box_ke
  procedure :: update_G_h
  procedure :: propagate_h
  procedure :: propagate_v_h_1
  procedure :: propagate_v_h_2
  procedure :: propagate_v_sys
  procedure :: propagate_r_sys
  procedure :: vVerlet_r_h_npt
  procedure :: dump_baro_state

end type type_barostat

contains

  subroutine init_barostat(baro, baro_type, ensemble, dt, h_ref, P_ext, nat, &
                           ndof, box_mass, volume, tau_P, iprint)

    ! Passed variables
    class(type_barostat), intent(inout)       :: baro
    character(40), intent(in)                 :: baro_type
    character(3), intent(in)                  :: ensemble
    real(double), dimension(3,3), intent(in)  :: h_ref
    real(double), intent(in)                  :: dt, P_ext, box_mass, tau_P, &
                                                 volume
    integer, intent(in)                       :: nat, ndof, iprint

    ! Local variables
    integer       :: i
    real(double)  :: W_h

    baro%iprint = iprint
    baro%baro_type = baro_type
    baro%ensemble = ensemble
    baro%dt = dt
    baro%nat = nat
    baro%ndof = ndof
    baro%k_B_md = one
    baro%h_0 = h_ref
    baro%h = h_ref
    baro%P_ext = P_ext
    baro%V = volume

    ! initialise MTTK
    baro%v_h = zero
    baro%W_h = box_mass
    baro%ke_box = zero
    baro%ident = zero
    do i=1,3
      baro%ident(i,i) = one
    end do
    baro%onfm = (one/baro%ndof)*baro%ident

    ! initialise Berendsen
    baro%beta = one
    baro%tau_P = tau_P

  end subroutine init_barostat

  ! Get the eigenvalues/vectors of v_h(0) + Tr[v_h(0)]/N_f + v_eta(1)
  subroutine get_lambda_cg(baro, v_eta)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro
    real(double), intent(in)                  :: v_eta

    ! Local variables
    real(double), dimension(:), allocatable   :: work
    real(double)                              :: tr_vh
    integer                                   :: i, lwork, info

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: Computing eigenvalues of v_h"
    lwork = 27
    allocate(work(lwork))

    baro%c_g = zero
    tr_vh = zero
    do i=1,3
      tr_vh = tr_vh + baro%v_h(i,i)
    end do

    ! baro%c_g = baro%v_h + (tr_vh/baro%ndof + v_eta)*baro%ident
    baro%c_g = baro%v_h + tr_vh/baro%ndof + v_eta
    ! diagonalise
    call dsyev('V', 'U', 3, baro%c_g, 3, baro%lambda, work, lwork, info)

  end subroutine get_lambda_cg

  subroutine get_lambda_sfac(baro, dt, dtfac, rscale, vscale)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro
    real(double), intent(in)                  :: dt, dtfac
    real(double), dimension(3), intent(out)   :: rscale, vscale

    ! local variables
    integer                                   :: i
    real(double)                              :: exp_eig, sinhx_x

    do i=1,3
      exp_eig = exp(dtfac*dt*baro%lambda(i))
      sinhx_x = poly_sinhx_x(dtfac*dt*baro%lambda(i))
      rscale(i) = exp_eig**2
      vscale(i) = exp_eig*dt*sinhx_x
    end do

  end subroutine get_lambda_sfac

  ! Get the box kinetic energy
  subroutine get_box_ke(baro)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro

    ! local variables
    real(double), dimension(3,3)              :: temp
    integer                                   :: i

    temp = matmul(transpose(baro%v_h), baro%v_h)
    baro%ke_box = zero
    do i=1,3
      baro%ke_box = baro%ke_box + temp(i,i)
    end do
    baro%ke_box = baro%W_h*baro%ke_box

  end subroutine get_box_ke

  ! Update the box forces
  subroutine update_G_h(baro, total_ke, virial_ke, stress, volume, &
                        T_ext, Q_1, G_nhc_1)

    ! Passed variables
    class(type_barostat), intent(inout)       :: baro
    real(double), dimension(3,3), intent(in)  :: virial_ke
    real(double), dimension(3,3), intent(in)  :: stress
    real(double), intent(in)                  :: total_ke
    real(double), intent(in)                  :: volume
    real(double), intent(in)                  :: T_ext
    real(double), intent(in)                  :: Q_1 ! mass of NHC k=1
    real(double), intent(out)                 :: G_nhc_1 ! force on NHC k=1

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: updating box force G_h"
    baro%stress = stress
    baro%virial_ke = virial_ke
    baro%G_h = (baro%ident*total_ke/baro%ndof + virial_ke + &
                volume*(stress - baro%ident*baro%P_ext))/baro%W_h
    ! The force on the first NHC heat bath is required for thermostat coupling
    G_nhc_1 = (total_ke + baro%ke_box - baro%ndof*baro%k_B_md*T_ext)/Q_1

  end subroutine update_G_h

  ! Update the box
  subroutine propagate_h(baro, rscale, h)

    ! Passed variables
    class(type_barostat), intent(inout)        :: baro
    real(double), dimension(3), intent(in)      :: rscale
    real(double), dimension(3,3), intent(inout) :: h

    ! local variables
    integer                                   :: i, j
    real(double), dimension(3,3)              :: htemp

    if (baro%iprint == 0) write(*,'(4x,a)') "MTTK: propagating h"
    htemp = matmul(transpose(baro%c_g), h)
    do i=1,3
      do j=1,3
        htemp(j,i) = htemp(j,i)*rscale(j)
      end do 
    end do
    h = matmul(baro%c_g, htemp)
    baro%h = h

  end subroutine propagate_h

  ! linear shift in box velocity
  subroutine propagate_v_h_1(baro, dt, dtfac)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro
    real(double), intent(in)                  :: dt, dtfac

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: propagating v_h linear shift"
    baro%v_h = baro%v_h + dtfac*dt*baro%G_h

  end subroutine propagate_v_h_1

  ! exponential shift in box velocity from thermostat coupling
  subroutine propagate_v_h_2(baro, dt, dtfac, v_eta)

    ! Passed variables
    class(type_barostat), intent(inout)      :: baro
    real(double), intent(in)                  :: dt, dtfac, v_eta

    if (baro%iprint == 0) write(*,'(6x,a)') "MTTK: propagating v_h exp factor"
    baro%v_h = baro%v_h*exp(-dtfac*dt*v_eta)

  end subroutine propagate_v_h_2

  ! coupling box and particle velocities
  subroutine propagate_v_sys(baro, dt, dtfac, v_eta_1, v)

    ! Passed variables
    class(type_barostat), intent(inout)         :: baro
    real(double), intent(in)                    :: dt, dtfac
    real(double), intent(in)                    :: v_eta_1
    real(double), dimension(:,:), intent(inout) :: v

    ! local variables
    real(double), dimension(3)                :: vscale, vtemp
    integer                                   :: i, j

    if (baro%iprint == 0) write(*,'(4x,a)') "MTTK: propagating particle velocities"
    call baro%get_lambda_cg(v_eta_1)

    do i=1,3
      vscale(i) = exp(-baro%lambda(i)*dt*dtfac)
    end do

    do i=1,baro%nat
      do j=1,3
        vtemp(j) = dot_product(v(i,:), baro%c_g(:,j))
      end do
      vtemp = vtemp*vscale
      do j=1,3
        v(i,j) = dot_product(vtemp, baro%c_g(:,j))
      end do
    end do

  end subroutine propagate_v_sys

  ! coupling box and particle positions. This is a little complicated and
  ! contains a polynomial expansion of sinh(x)/x, so I will do it in one
  ! subroutine for now, as in the MTTK paper.
  subroutine propagate_r_sys(baro, dt, rscale, vscale, r,  v)

    ! Passed variables
    class(type_barostat), intent(inout)        :: baro
    real(double), intent(in)                    :: dt
    real(double), dimension(:,:), intent(inout) :: r
    real(double), dimension(:,:), intent(in)    :: v
    real(double), dimension(3), intent(in)      :: rscale, vscale

    ! local variables
    integer                                     :: i, j
    real(double), dimension(3)                  :: u, uv
    real(double)                                :: exp_eig, sinhx_x

    if (baro%iprint == 0) write(*,'(4x,a)') "MTTK: propagating particle positions"
    do i=1,baro%nat
      do j=1,3
        u(j) = dot_product(r(i,:), baro%c_g(:,j))
        uv(j) = dot_product(v(i,:), baro%c_g(:,j))
        u(j) = u(j)*rscale(j) + uv(j)*vscale(j)
      end do
      do j=1,3
        r(i,j) = dot_product(u, baro%c_g(j,:))
      end do
    end do

  end subroutine propagate_r_sys

  ! Propagate the ionic positions and box of the extended system. Temperature
  ! coupling determined by v_eta
  subroutine vVerlet_r_h_npt(baro, dt, h, v, r, v_eta)

    ! Passed variables
    class(type_barostat), intent(inout)        :: baro
    real(double), intent(in)                    :: dt
    real(double), intent(in)                    :: v_eta ! thermostat velocity
    real(double), dimension(3,3), intent(inout) :: h
    real(double), dimension(:,:), intent(inout) :: r
    real(double), dimension(:,:), intent(inout) :: v

    ! local variables
    real(double), dimension(3)                  :: rscale, vscale

    ! Preceded by velocity Verlet dt/2 step
    call baro%get_lambda_cg(v_eta)
    call baro%get_lambda_sfac(dt, half, rscale, vscale)
    call baro%propagate_r_sys(dt, rscale, vscale, r, v)
    call baro%propagate_h(rscale, h)
    ! Proceded by velocity Verlet dt/2 step

  end subroutine vVerlet_r_h_npt

  ! Debugging routine: dump the state of the barostat
  subroutine dump_baro_state(baro, step, funit)

    ! Passed variables
    class(type_barostat), intent(inout)        :: baro
    integer, intent(in)                        :: step, funit

    ! local variables
    integer                                    :: i

    write(funit,'("step   ",i12)') step
    write(funit,'("box ke ",f12.4)') baro%ke_box
    write(funit,'(8x,a)') "Lattice vectors"
    do i=1,3
      write(funit,'(4x,3f12.4)') baro%h(i,:)
    end do
    write(funit,'(8x,a,29x,a)') "Virial KE", "Stress"
    do i=1,3
      write(funit,'(4x,3f12.4,4x,3f12.4)') baro%virial_ke(i,:), &
                                           baro%stress(i,:)
    end do
    write(funit,'(8x,a,12x,a)') 'lambda', 'c_g'
    do i=1,3
      write(funit,'(4x,f12.4,6x,3f12.4)') baro%lambda(i), baro%c_g(i,:)
    end do
    write(funit,'(8x,a)') "G_h"
    do i=1,3
      write(funit,'(4x,3f12.4)') baro%G_h(i,:)
    end do
    write(funit,*)

  end subroutine dump_baro_state

end module barostat_module
